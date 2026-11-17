import pytest
from pydicom.uid import CTImageStorage
from pydicom.uid import ExplicitVRLittleEndian
from pydicom.dataset import Dataset, FileMetaDataset
from pj_medkit.phase import (PHASE_PLAIN, PHASE_DELAY, PHASE_VENOUS,
                             PHASE_ARTERIAL, PHASE_RULE_VERSION,
                             _dicom_seconds, classify_dir, classify_phases,
                             read_series_header)

# 08:00:00 的当天秒数，用例里的采集时间都在它之后
T0 = 8 * 3600


def _series(key,
            contrast=False,
            acq=None,
            bolus=None,
            desc=None,
            protocol=None) -> dict:
    return {
        "key": key,
        "desc": desc,
        "protocol": protocol,
        "contrast": contrast,
        "acq": acq,
        "bolus": bolus,
    }


def _write_dicom(path,
                 *,
                 series_uid,
                 series_number,
                 instance_number,
                 modality="CT",
                 **tags) -> None:
    """
    写一个只有文件头的 DICOM 文件，tags 为要额外写入的字段
    """
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = CTImageStorage
    meta.MediaStorageSOPInstanceUID = f"{series_uid}.{instance_number}"
    meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds = Dataset()
    ds.file_meta = meta
    ds.SpecificCharacterSet = "ISO_IR 192"
    ds.SOPClassUID = CTImageStorage
    ds.SOPInstanceUID = f"{series_uid}.{instance_number}"
    ds.Modality = modality
    ds.SeriesInstanceUID = series_uid
    ds.SeriesNumber = series_number
    ds.InstanceNumber = instance_number
    for name, value in tags.items():
        setattr(ds, name, value)
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.save_as(path, enforce_file_format=True)


def _write_series(root, uid, number, count=3, **tags) -> None:
    for i in range(1, count + 1):
        _write_dicom(root / f"s{number}" / f"{i}.dcm",
                     series_uid=uid,
                     series_number=number,
                     instance_number=i,
                     **tags)


def test_dicom_seconds():
    assert _dicom_seconds(None, "080000") == T0
    assert _dicom_seconds(None, "08:00:00") == T0
    assert _dicom_seconds(None, "080000.500000") == T0 + 0.5
    assert _dicom_seconds(None, "0800") == T0
    assert _dicom_seconds(None, None) is None
    assert _dicom_seconds(None, "abc") is None
    one_day = _dicom_seconds("20260102", "080000") - _dicom_seconds(
        "20260101", "080000")
    assert one_day == 86400


def test_single_series_is_plain():
    assert classify_phases([_series("1", acq=T0)]) == {"1": PHASE_PLAIN}
    assert classify_phases([_series("1")]) == {"1": PHASE_PLAIN}


def test_three_or_more_enhanced_by_acquisition_order():
    series = [
        _series("5", True, T0 + 300),
        _series("2", False, T0),
        _series("4", True, T0 + 70),
        _series("3", True, T0 + 30),
        _series("6", True, T0 + 600),
    ]
    assert classify_phases(series) == {
        "2": PHASE_PLAIN,
        "3": PHASE_ARTERIAL,
        "4": PHASE_VENOUS,
        "5": PHASE_DELAY,
        "6": PHASE_DELAY,
    }


@pytest.mark.parametrize("elapsed, phase", [
    (22, PHASE_ARTERIAL),
    (41.9, PHASE_ARTERIAL),
    (42, PHASE_VENOUS),
    (84.9, PHASE_VENOUS),
    (85, PHASE_DELAY),
    (180, PHASE_DELAY),
])
def test_one_enhanced_by_seconds_after_bolus(elapsed, phase):
    series = [
        _series("1", False, T0),
        _series("2", True, T0 + 100 + elapsed, bolus=T0 + 100),
    ]
    assert classify_phases(series) == {"1": PHASE_PLAIN, "2": phase}


def test_two_enhanced_by_seconds_after_bolus():
    series = [
        _series("1", False, T0),
        _series("2", True, T0 + 130, bolus=T0 + 100),
        _series("3", True, T0 + 160, bolus=T0 + 100),
    ]
    assert classify_phases(series) == {
        "1": PHASE_PLAIN,
        "2": PHASE_ARTERIAL,
        "3": PHASE_VENOUS,
    }


@pytest.mark.parametrize("rise, first, second", [
    (5, PHASE_ARTERIAL, PHASE_VENOUS),
    (14.9, PHASE_ARTERIAL, PHASE_VENOUS),
    (15, PHASE_VENOUS, PHASE_DELAY),
    (50, PHASE_VENOUS, PHASE_DELAY),
])
def test_two_enhanced_without_bolus_by_liver_rise(rise, first, second):
    series = [
        _series("1", False, T0),
        _series("2", True, T0 + 100),
        _series("3", True, T0 + 200),
    ]
    liver_hu = {"1": 55.0, "2": 55.0 + rise, "3": 110.0}
    assert classify_phases(series, liver_hu) == {
        "1": PHASE_PLAIN,
        "2": first,
        "3": second,
    }


def test_enhanced_without_bolus_or_liver_falls_back_to_name():
    series = [
        _series("1", False, T0),
        _series("2", True, T0 + 100, desc="Abdomen Arterial 5.0"),
        _series("3", True, T0 + 200, protocol="肝脏增强 门脉期"),
    ]
    assert classify_phases(series) == {
        "1": PHASE_PLAIN,
        "2": PHASE_ARTERIAL,
        "3": PHASE_VENOUS,
    }


def test_enhanced_undetermined_is_left_out():
    series = [
        _series("1", False, T0),
        _series("2", True, T0 + 100, desc="Abdomen 5.0"),
    ]
    assert classify_phases(series) == {"1": PHASE_PLAIN}
    # 同时写了两个期相的名字不采用
    series[1]["desc"] = "Arterial Venous"
    assert classify_phases(series) == {"1": PHASE_PLAIN}


def test_enhanced_missing_acquisition_time_uses_name():
    series = [
        _series("1", False, T0),
        _series("2", True, None, desc="Delay"),
        _series("3", True, T0 + 100),
    ]
    assert classify_phases(series) == {"1": PHASE_PLAIN, "2": PHASE_DELAY}


def test_no_contrast_field_uses_description_words():
    series = [
        _series("1", acq=T0, desc="Abdomen B30f"),
        _series("2", acq=T0 + 100, desc="Abdomen A"),
        _series("3", acq=T0 + 160, desc="Abdomen P"),
        _series("4", acq=T0 + 300, desc="Abdomen D"),
    ]
    assert classify_phases(series) == {
        "1": PHASE_PLAIN,
        "2": PHASE_ARTERIAL,
        "3": PHASE_VENOUS,
        "4": PHASE_DELAY,
    }


def test_no_contrast_field_reconstruction_kernel_is_not_a_phase():
    series = [
        _series("1", acq=T0, desc="Abdomen B30f"),
        _series("2", acq=T0 + 100, desc="Abdomen D30f"),
    ]
    # 两个序列都没写期相，也没有注药时间和肝脏 CT 值，第二个判断不了
    assert classify_phases(series) == {"1": PHASE_PLAIN}


def test_no_contrast_field_no_description_uses_bolus_and_liver():
    series = [
        _series("1", acq=T0),
        _series("2", acq=T0 + 60),
        _series("3", acq=T0 + 400, bolus=T0 + 370),
        _series("4", acq=T0 + 430, bolus=T0 + 370),
    ]
    # 序列 2 肝脏和最早的序列差不多，是重复扫的平扫
    liver_hu = {"1": 55.0, "2": 55.8}
    assert classify_phases(series, liver_hu) == {
        "1": PHASE_PLAIN,
        "2": PHASE_PLAIN,
        "3": PHASE_ARTERIAL,
        "4": PHASE_VENOUS,
    }
    # 肝脏明显升高则按增强处理，三个增强序列按先后排
    liver_hu = {"1": 55.0, "2": 70.0}
    assert classify_phases(series, liver_hu) == {
        "1": PHASE_PLAIN,
        "2": PHASE_ARTERIAL,
        "3": PHASE_VENOUS,
        "4": PHASE_DELAY,
    }


def test_no_contrast_field_missing_acquisition_time_gives_nothing():
    series = [_series("1", acq=T0), _series("2")]
    assert classify_phases(series) == {}


def test_read_series_header(tmp_path):
    path = tmp_path / "a.dcm"
    _write_dicom(path,
                 series_uid="1.2.3",
                 series_number=4,
                 instance_number=1,
                 SeriesDescription="Abdomen A",
                 ProtocolName="Liver",
                 ContrastBolusAgent="Iohexol",
                 AcquisitionDate="20260102",
                 AcquisitionTime="000010",
                 ContrastBolusStartTime="235950")
    header = read_series_header(path)
    assert header["desc"] == "Abdomen A"
    assert header["protocol"] == "Liver"
    assert header["contrast"] is True
    # 注药在零点前、采集在零点后，间隔仍是 20 秒
    assert header["acq"] - header["bolus"] == 20


@pytest.mark.parametrize("agent", ["", "NONE", "no", "0", "无"])
def test_read_series_header_no_contrast_values(tmp_path, agent):
    path = tmp_path / "a.dcm"
    _write_dicom(path,
                 series_uid="1.2.3",
                 series_number=4,
                 instance_number=1,
                 ContrastBolusAgent=agent)
    assert read_series_header(path)["contrast"] is False


def test_classify_dir(tmp_path):
    common = {"AcquisitionDate": "20260102"}
    _write_series(tmp_path, "1.2.2", 2, AcquisitionTime="080000", **common)
    _write_series(tmp_path,
                  "1.2.3",
                  3,
                  AcquisitionTime="080530",
                  ContrastBolusAgent="Iohexol",
                  ContrastBolusStartTime="080500",
                  SeriesDescription="动脉期",
                  **common)
    _write_series(tmp_path,
                  "1.2.4",
                  4,
                  AcquisitionTime="080600",
                  ContrastBolusAgent="Iohexol",
                  ContrastBolusStartTime="080500",
                  **common)
    # 不是 DICOM 的文件、非 CT 的文件不参与
    (tmp_path / "readme.txt").write_text("not dicom")
    _write_dicom(tmp_path / "seg" / "1.dcm",
                 series_uid="1.2.9",
                 series_number=9,
                 instance_number=1,
                 modality="SEG")
    result = classify_dir(tmp_path)
    assert result == {
        "rule_version":
        PHASE_RULE_VERSION,
        "series": [{
            "series_uid": "1.2.2",
            "series_number": 2,
            "description": None,
            "phase": PHASE_PLAIN,
        }, {
            "series_uid": "1.2.3",
            "series_number": 3,
            "description": "动脉期",
            "phase": PHASE_ARTERIAL,
        }, {
            "series_uid": "1.2.4",
            "series_number": 4,
            "description": None,
            "phase": PHASE_VENOUS,
        }],
    }


def test_classify_dir_series_uids_and_undetermined(tmp_path):
    common = {"AcquisitionDate": "20260102"}
    _write_series(tmp_path, "1.2.2", 2, AcquisitionTime="080000", **common)
    _write_series(tmp_path,
                  "1.2.3",
                  3,
                  AcquisitionTime="080530",
                  ContrastBolusAgent="Iohexol",
                  **common)
    _write_series(tmp_path,
                  "1.2.4",
                  4,
                  AcquisitionTime="080600",
                  ContrastBolusAgent="Iohexol",
                  **common)
    _write_series(tmp_path,
                  "1.2.5",
                  5,
                  AcquisitionTime="080900",
                  ContrastBolusAgent="Iohexol",
                  **common)
    full = classify_dir(tmp_path)
    assert [s["phase"] for s in full["series"]
            ] == [PHASE_PLAIN, PHASE_ARTERIAL, PHASE_VENOUS, PHASE_DELAY]
    # 只留一个增强序列：没有注药时间、肝脏 CT 值和期相名，判断不了
    part = classify_dir(tmp_path, series_uids=["1.2.2", "1.2.3"])
    assert [(s["series_uid"], s["phase"])
            for s in part["series"]] == [("1.2.2", PHASE_PLAIN),
                                         ("1.2.3", None)]
    # 补上肝脏 CT 值后能判断
    part = classify_dir(tmp_path,
                        liver_hu={
                            "1.2.2": 55.0,
                            "1.2.3": 60.0
                        },
                        series_uids=["1.2.2", "1.2.3"])
    assert [s["phase"]
            for s in part["series"]] == [PHASE_PLAIN, PHASE_ARTERIAL]


def test_classify_dir_empty(tmp_path):
    assert classify_dir(tmp_path) == {
        "rule_version": PHASE_RULE_VERSION,
        "series": []
    }
