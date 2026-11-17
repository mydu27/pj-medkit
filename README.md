# pj-medkit

医学影像公共工具包：把各项目都要用的方法收在一处，统一安装、统一版本，不用每个项目各写一份。

## 文档目录

按要做的事找对应的模块和文档：

|我想做什么|用哪个模块|提供的函数|文档|
|:---|:---|:---|:---|
|判断 CT 检查里各序列是平扫、动脉期、静脉期还是延迟期|`pj_medkit.phase`|`classify_dir()`、`classify_phases()`、`read_series_header()`|[期相识别](docs/期相识别.md)|

每个模块一篇文档，都在 [docs/](docs/) 目录下。各版本改了什么见 [CHANGELOG.md](CHANGELOG.md)。

## 安装

Python 3.10 及以上，依赖 pydicom、loguru。

每个版本打成一个 `.whl` 文件，挂在仓库的发布页。使用方下载后放进自己项目的 `wheels/` 目录，在 `requirements.txt` 里按文件路径安装，构建时不需要访问本仓库：

```
./wheels/pj_medkit-1.0.0.0-py3-none-any.whl
```

用 Docker 构建的项目，在 `pip install` 之前把 `wheels/` 目录拷进镜像：

```dockerfile
COPY requirements.txt .
COPY wheels ./wheels
RUN pip install -r requirements.txt
```

## 使用

```python
from pj_medkit.phase import classify_dir

result = classify_dir("/data/检查目录")
# {
#     "rule_version": 1,
#     "series": [
#         {"series_uid": "1.2.840...", "series_number": 2, "description": "Abdomen", "phase": 0},
#         {"series_uid": "1.2.840...", "series_number": 3, "description": "Abdomen A", "phase": 1},
#     ],
# }
```

期相编码：0 平扫，1 动脉期，2 静脉期，3 延迟期；判断不了为 `None`。详细用法和判断规则见 [docs/期相识别.md](docs/期相识别.md)。

## 目录结构

```
pj_medkit/      包代码，一个功能一个模块
tests/          测试，文件名与模块对应（test_模块名.py）
docs/           各模块的参考文档，一个模块一篇
CHANGELOG.md    更新日志
```

参考文档的写法仿照 Python 官方文档：模块简介和快速示例 → 每个函数一节（签名、说明、参数、返回值、示例）→ 常量 → 补充说明 → 示例。函数是哪个版本加的、哪个版本改过，在该函数末尾标"在 x 版本加入""在 x 版本发生变更"。

## 开发

跑测试：

```
python -m pytest
```

新增一个方法时：在 `pj_medkit/` 下加模块，在 `tests/` 下加对应测试，在 `docs/` 下加参考文档，并在上面的"文档目录"表里登记一行。

代码约定：

- 核心方法写成不依赖具体项目的纯函数：入参是普通的值，不绑定某个项目的目录结构、数据库或接口格式。
- 代码按 yapf 默认风格格式化；import 按整行长度从短到长排列。
- 改动会影响结果的规则或阈值时，升对应模块的规则版本号（如 `PHASE_RULE_VERSION`），并补测试用例。

## 发版

1. 改 `pj_medkit/__init__.py` 里的 `__version__`，在 `CHANGELOG.md` 写明更新内容。
2. 测试全部通过后打包，产物在 `dist/` 下：

   ```
   python -m pip wheel . --no-deps -w dist
   ```

3. 提交、打标签（如 `v1.0.0.0`），把 `dist/` 下的 `.whl` 文件挂到仓库发布页。
4. 通知使用方换新文件并改 `requirements.txt` 里的文件名。
