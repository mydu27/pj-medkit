import os
import re
import traceback
from loguru import logger
from datetime import date
from pydicom import dcmread
from pydicom.errors import InvalidDicomError

# 期相规则版本：判断规则或阈值有变化时加一，随结果一起返回，便于区分存下来的结果是哪版规则判的
PHASE_RULE_VERSION = 1

# 期相编码
PHASE_PLAIN = 0
PHASE_ARTERIAL = 1
PHASE_VENOUS = 2
PHASE_DELAY = 3
PHASE_NAMES = {
    PHASE_PLAIN: "平扫",
    PHASE_ARTERIAL: "动脉期",
    PHASE_VENOUS: "静脉期",
    PHASE_DELAY: "延迟期",
}
# 注药后秒数分界：< 动脉期上限为动脉期，< 静脉期上限为静脉期，其余为延迟期。
# 实测动脉期 22~40s、静脉期 43~72s、延迟期 91s 以后
ARTERIAL_MAX_SECONDS = 42
VENOUS_MAX_SECONDS = 85
# 无注药时间时，首个增强期肝脏比平扫高不到该值视为动脉期(动脉期肝脏强化不明显)
ARTERIAL_MAX_RISE_HU = 15
# 整例都无造影剂字段时，肝脏比最早序列高不到该值且无注药时间的视为重复扫的平扫
# (实测两次平扫相差 0.8)
REPEAT_PLAIN_MAX_RISE_HU = 2
# ContrastBolusAgent 取这些值时视为未用造影剂
_NO_CONTRAST_VALUES = frozenset({"", "NONE", "NO", "N", "0", "无"})
# 整例都没填造影剂字段时，序列描述里表示增强期的词(英文整词、中文子串)。
# 只看序列描述不看协议名：协议名整次检查共用，带 CE 时平扫序列也会被误判
_ENHANCED_DESC_TOKENS = frozenset({
    "A", "AP", "ART", "ARTERIAL", "P", "PV", "PORTAL", "V", "VENOUS", "D",
    "DL", "DELAY", "DELAYED", "C+", "CE"
})
_ENHANCED_DESC_WORDS = ("动脉", "静脉", "门脉", "延迟", "增强")
# 增强序列其他规则都判断不了时，序列描述或协议名里写明的期相(英文整词、中文子串)。
# 此时已确定是增强，协议名整次检查共用也不影响；单字母容易误判，不认
_PHASE_NAME_TOKENS = {
    PHASE_ARTERIAL: frozenset({"ART", "ARTERIAL"}),
    PHASE_VENOUS: frozenset({"PV", "PORTAL", "VENOUS"}),
    PHASE_DELAY: frozenset({"DELAY", "DELAYED"}),
}
_PHASE_NAME_WORDS = {
    PHASE_ARTERIAL: ("动脉", ),
    PHASE_VENOUS: ("静脉", "门脉"),
    PHASE_DELAY: ("延迟", ),
}
# 期相判断要读的文件头字段，只读这些比读整个文件头快
_HEADER_TAGS = [
    "Modality", "SeriesDate", "SeriesNumber", "ProtocolName", "InstanceNumber",
    "AcquisitionDate", "AcquisitionTime", "SeriesDescription",
    "SeriesInstanceUID", "ContrastBolusAgent", "ContrastBolusStartTime"
]


def _dicom_seconds(dicom_date, dicom_time) -> float | None:
    """
    DICOM DA(YYYYMMDD) + TM(HHMMSS.ffffff，兼容带冒号的旧写法)转为可比较的秒数；
    没有日期时只算当天秒数，时间缺失或格式不对返回 None
    """
    text = str(dicom_time or "").replace(":", "").strip()
    try:
        seconds = (int(text[0:2]) * 3600 + int(text[2:4]) * 60 +
                   float(text[4:] or 0))
    except ValueError:
        return None
    day = str(dicom_date or "").strip()
    if len(day) == 8 and day.isdigit():
        seconds += date(int(day[:4]), int(day[4:6]), int(
            day[6:])).toordinal() * 86400
    return seconds


def _series_header(ds) -> dict:
    """
    从一层图像的文件头取期相判断要用的字段：序列描述、协议名、是否用造影剂、
    采集时间、注药开始时间。
    采集时间用 AcquisitionTime 而非 SeriesTime：后者常是重建时间，会打乱先后。
    """
    acq_date = ds.get("AcquisitionDate") or ds.get("SeriesDate")
    acq = _dicom_seconds(acq_date, ds.get("AcquisitionTime"))
    bolus = _dicom_seconds(None, ds.get("ContrastBolusStartTime"))
    if bolus is not None and acq is not None:
        # 注药时间只有时分秒，按采集当天补齐，跨零点的往前退一天
        bolus += acq - acq % 86400
        if bolus > acq:
            bolus -= 86400
    agent = str(ds.get("ContrastBolusAgent") or "").strip().upper()
    return {
        "desc": ds.get("SeriesDescription"),
        "protocol": ds.get("ProtocolName"),
        "contrast": agent not in _NO_CONTRAST_VALUES,
        "acq": acq,
        "bolus": bolus,
    }


def read_series_header(path) -> dict:
    """
    读一个 DICOM 文件的文件头，返回 classify_phases 要用的
    {desc, protocol, contrast, acq, bolus}，调用方补上 key 后即可传入。
    一个序列读一层即可，建议取中间一层。文件读不了时抛出 pydicom 的异常
    """
    ds = dcmread(path, stop_before_pixels=True, specific_tags=_HEADER_TAGS)
    return _series_header(ds)


def _phase_by_bolus(series: dict) -> int | None:
    """
    按注药后秒数判断期相，没有注药时间返回 None
    """
    if series["bolus"] is None or series["acq"] is None:
        return None
    elapsed = series["acq"] - series["bolus"]
    if elapsed < ARTERIAL_MAX_SECONDS:
        return PHASE_ARTERIAL
    if elapsed < VENOUS_MAX_SECONDS:
        return PHASE_VENOUS
    return PHASE_DELAY


def _desc_marks_enhanced(desc: str | None) -> bool:
    """
    序列描述是否写明是增强期。英文按整词匹配，避免把 B30f、D30f 这类重建核名误当成期相
    """
    text = str(desc or "").upper()
    if any(word in text for word in _ENHANCED_DESC_WORDS):
        return True
    return any(token in _ENHANCED_DESC_TOKENS
               for token in re.split(r"[\s\-_/,()]+", text))


def _phase_by_name(series: dict) -> int | None:
    """
    序列描述、协议名里写明的增强期相；两处合起来只指向一个期相时返回它，
    没写或写了不止一个期相返回 None。只用于已确定是增强、其他规则都判断不了的序列
    """
    found = set()
    for text in (series.get("desc"), series.get("protocol")):
        text = str(text or "").upper()
        tokens = set(re.split(r"[\s\-_/,()]+", text))
        for phase, names in _PHASE_NAME_TOKENS.items():
            if tokens & names or any(word in text
                                     for word in _PHASE_NAME_WORDS[phase]):
                found.add(phase)
    return found.pop() if len(found) == 1 else None


def classify_phases(series: list[dict], liver_hu: dict | None = None) -> dict:
    """
    按规则判断各序列期相，series 每项为 {key, desc, protocol, contrast, acq, bolus}，
    key 是调用方给序列起的标识(序列号或 SeriesInstanceUID 均可)；liver_hu 为
    {key: 肝脏平均 CT 值}，没有可不传。返回 {key: 期相}，判断不了的序列不出现：
    1. 无造影剂字段的为平扫；整例都无该字段却有多个序列时，可能是医院不填该字段：
       序列描述里有增强期用词的，带词的为增强、其余为平扫(检查可能根本没有平扫)；
       都没有时，采集最早的为平扫，其余有注药时间或肝脏比它高出
       REPEAT_PLAIN_MAX_RISE_HU 的按增强处理，差不多的视为重复扫的平扫，
       两样都没有的判断不了
    2. 增强序列按采集时间排序：3 个及以上依次为动脉、静脉，其余为延迟；
       1~2 个时有注药时间按注药后秒数判断，没有则首个看肝脏比平扫高多少，
       第二个顺延(动脉之后为静脉，静脉之后为延迟)
    3. 增强序列按以上规则判断不了的(含缺采集时间的)，最后看序列描述、协议名里
       写明的期相(_phase_by_name)
    采集时间相同的序列按传入顺序排先后
    """
    liver_hu = liver_hu or {}
    if any(s["contrast"] for s in series):
        plain = [s for s in series if not s["contrast"]]
        enhanced = [s for s in series if s["contrast"]]
    elif len(series) > 1:
        if any(s["acq"] is None for s in series):
            return {}
        ordered = sorted(series, key=lambda s: s["acq"])
        marked = [_desc_marks_enhanced(s.get("desc")) for s in ordered]
        if any(marked):
            plain = [s for s, m in zip(ordered, marked) if not m]
            enhanced = [s for s, m in zip(ordered, marked) if m]
        else:
            plain, enhanced = ordered[:1], []
            first_hu = liver_hu.get(ordered[0]["key"])
            for s in ordered[1:]:
                if s["bolus"] is not None:
                    enhanced.append(s)
                elif first_hu is not None and s["key"] in liver_hu:
                    if (liver_hu[s["key"]] - first_hu
                            >= REPEAT_PLAIN_MAX_RISE_HU):
                        enhanced.append(s)
                    else:
                        plain.append(s)
    else:
        plain, enhanced = series, []
    phases = {s["key"]: PHASE_PLAIN for s in plain}
    if any(s["acq"] is None for s in enhanced):
        for s in enhanced:
            phase = _phase_by_bolus(s)
            if phase is None:
                phase = _phase_by_name(s)
            if phase is not None:
                phases[s["key"]] = phase
        return phases
    enhanced = sorted(enhanced, key=lambda s: s["acq"])
    if len(enhanced) >= 3:
        for i, s in enumerate(enhanced):
            phases[s["key"]] = ((PHASE_ARTERIAL,
                                 PHASE_VENOUS)[i] if i < 2 else PHASE_DELAY)
        return phases
    if not enhanced:
        return phases
    first = _phase_by_bolus(enhanced[0])
    base = min((liver_hu[s["key"]] for s in plain if s["key"] in liver_hu),
               default=None)
    if first is None and base is not None and enhanced[0]["key"] in liver_hu:
        rise = liver_hu[enhanced[0]["key"]] - base
        first = (PHASE_ARTERIAL
                 if rise < ARTERIAL_MAX_RISE_HU else PHASE_VENOUS)
    if first is None:
        first = _phase_by_name(enhanced[0])
    if first is not None:
        phases[enhanced[0]["key"]] = first
    if len(enhanced) == 2:
        second = _phase_by_bolus(enhanced[1])
        if second is None and first is not None:
            second = (PHASE_VENOUS if first == PHASE_ARTERIAL else PHASE_DELAY)
        if second is None:
            second = _phase_by_name(enhanced[1])
        if second is not None:
            phases[enhanced[1]["key"]] = second
    return phases


def guess_phases(series: list[dict], phases: dict) -> dict:
    """
    给 classify_phases 判断不了的序列按扫描先后强行补一个期相，series 同
    classify_phases，phases 是它的返回值。返回 {key: 期相}，只含 phases 里没有的序列。
    补的期相没有依据，只用于必须每个序列都有期相才能继续的场合：
    1. 先后顺序：序列都有采集时间时按采集时间，有缺的整体按传入顺序
    2. 判断不了的序列按增强处理。整例都没填造影剂字段时例外：序列描述里有增强期
       用词的，不带词的补为平扫；都没有时最早的补为平扫
    3. 增强序列按先后补位：已判出的保留，没判出的接在前一个增强期之后(动脉、静脉、
       延迟)，且不晚于后面最近一个已判出的期相
    4. 整例只有一个增强序列且判断不了时补为静脉期(单期增强多为静脉期)
    """
    if all(s["key"] in phases for s in series):
        return {}
    if all(s["acq"] is not None for s in series):
        ordered = sorted(series, key=lambda s: s["acq"])
    else:
        ordered = list(series)
    if any(s["contrast"] for s in series):
        plain = []
    else:
        plain = [s for s in ordered if not _desc_marks_enhanced(s.get("desc"))]
        if len(plain) == len(ordered):
            plain = ordered[:1]
    guessed = {s["key"]: PHASE_PLAIN for s in plain if s["key"] not in phases}
    enhanced = [
        s for s in ordered
        if phases.get(s["key"], guessed.get(s["key"])) != PHASE_PLAIN
    ]
    if len(enhanced) == 1:
        guessed.setdefault(enhanced[0]["key"], PHASE_VENOUS)
        return guessed
    # 期相编码按动脉、静脉、延迟递增，前一个加一即下一个期相
    previous = PHASE_PLAIN
    for i, s in enumerate(enhanced):
        phase = phases.get(s["key"])
        if phase is None:
            later = next((phases[x["key"]]
                          for x in enhanced[i + 1:] if x["key"] in phases),
                         PHASE_DELAY)
            phase = min(previous + 1, later)
            guessed[s["key"]] = phase
        previous = phase
    return guessed


def _scan_ct_series(dicom_dir, series_uids: list[str] | None) -> dict:
    """
    扫目录(含子目录)下的 CT 图像文件头，按 SeriesInstanceUID 分组返回 {uid: [ds]}。
    不是 DICOM 的文件、非 CT 的文件直接跳过；series_uids 不为空时只留这些序列
    """
    wanted = set(series_uids) if series_uids else None
    groups: dict[str, list] = {}
    stack_logged = False
    for root, _, names in os.walk(dicom_dir):
        for name in names:
            path = os.path.join(root, name)
            try:
                ds = dcmread(path,
                             stop_before_pixels=True,
                             specific_tags=_HEADER_TAGS)
            except InvalidDicomError:
                continue
            except Exception as e:
                # 文件读不了时往往一批都读不了，一次扫描只有第一个带栈
                if stack_logged:
                    logger.warning(f"期相识别读取文件头失败,跳过: {path} ({e!r})")
                else:
                    stack_logged = True
                    logger.warning(f"期相识别读取文件头失败,跳过: {path}\n"
                                   f"{traceback.format_exc()}")
                continue
            uid = ds.get("SeriesInstanceUID")
            if ds.get("Modality") != "CT" or not uid:
                continue
            if wanted is not None and uid not in wanted:
                continue
            groups.setdefault(str(uid), []).append(ds)
    return groups


def classify_dir(dicom_dir,
                 liver_hu: dict | None = None,
                 series_uids: list[str] | None = None,
                 fallback: bool = False) -> dict:
    """
    判断一个检查目录(含子目录)里各 CT 序列的期相。目录里的 CT 序列都参与判断，
    序列个数会影响结果，所以要么传筛选后的目录，要么用 series_uids 指定参与的序列。
    liver_hu 为 {SeriesInstanceUID: 肝脏平均 CT 值}，没有可不传。
    每个序列取中间一层的文件头；采集时间相同的序列按序列号排先后。
    返回 {"rule_version": 规则版本, "series": [{series_uid, series_number,
    description, phase, guessed}]}，series 按序列号排序。判断不了的序列 phase 为
    None；fallback 为 True 时改为按扫描先后补一个(guess_phases，缺采集时间时按
    序列号排先后)，补的 guessed 为 True
    """
    groups = _scan_ct_series(dicom_dir, series_uids)
    series = []
    for uid, datasets in groups.items():
        datasets.sort(key=lambda ds: int(ds.get("InstanceNumber") or 0))
        middle = datasets[len(datasets) // 2]
        number = middle.get("SeriesNumber")
        series.append({
            "key": uid,
            "number": None if number is None else int(number),
            **_series_header(middle)
        })
    series.sort(
        key=lambda s: (s["number"] is None, s["number"] or 0, s["key"]))
    phases = classify_phases(series, liver_hu)
    guessed = guess_phases(series, phases) if fallback else {}
    result = {
        "rule_version":
        PHASE_RULE_VERSION,
        "series": [{
            "series_uid": s["key"],
            "series_number": s["number"],
            "description": s["desc"],
            "phase": phases.get(s["key"], guessed.get(s["key"])),
            "guessed": s["key"] in guessed,
        } for s in series],
    }
    summary = {s["series_number"]: s["phase"] for s in result["series"]}
    guessed_numbers = [
        s["series_number"] for s in result["series"] if s["guessed"]
    ]
    logger.info(f"期相识别完成: {dicom_dir} {summary} 兜底补的序列: {guessed_numbers}")
    return result
