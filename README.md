# PS1 小行星历史数据回溯

本项目提供一个 Python 脚本，用于检索和整理 **Pan-STARRS1（PS1）** 巡天中的小行星历史观测。脚本以目标编号为输入，通过加拿大天文数据中心（CADC）的 SSOS 服务查询观测记录；可选地利用 JPL Horizons 更新目标星历；随后取得 PS1 图像切片、匹配 PS1 目录检测结果，并导出数据与可视化图件。

> 当前实现以小行星 `2016 OR17` 为示例目标。运行前请在脚本中调整 `source` 变量，并确认对所使用的数据服务具有适当的访问权限。

## 功能概览

| 功能 | 说明 |
|---|---|
| 历史观测检索 | 通过 CADC SSOS 查询 PS1 的移动天体观测记录。 |
| 高精度星历 | 可通过 JPL Horizons 获取观测时刻的高精度赤经和赤纬，并在服务不可用时回退至 CADC 提供的位置。 |
| 图像检索 | 生成 PS1 cutout 链接，下载 JPG、PNG 或 FITS 格式的图像切片。 |
| 目录匹配 | 可选地使用 MAST CasJobs 查询 Pan-STARRS DR2 Detection 表，匹配邻近检测源。 |
| 文件输出 | 输出观测表、检测表、轨迹图、缩略图拼图、单张缩略图和原始 warp FITS 文件。 |

## 前置条件

建议使用 Python 3.10 或更高版本，并先创建独立虚拟环境。脚本需要访问互联网，且依赖以下外部服务：

| 服务 | 用途 | 是否必需 |
|---|---|---|
| [CADC SSOS](https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/ssos/) | 检索移动天体观测记录。 | 是 |
| [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/) | 提供高精度星历。 | 否；可将 `use_jpl` 设为 `False`。 |
| [PS1 图像服务](https://ps1images.stsci.edu/) | 获取 cutout 与原始 warp FITS 图像。 | 是（如需图像输出）。 |
| [MAST CasJobs](https://mastweb.stsci.edu/ps1casjobs/) | 查询 PS1 目录检测源。 | 否；不可用时脚本会继续运行，但不生成检测匹配结果。 |

## 安装

克隆仓库后，在项目根目录执行以下命令创建并启用虚拟环境：

```bash
git clone https://github.com/yizhizhangxiu/PS1-data-retroanalysis.git
cd PS1-data-retroanalysis
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell：.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

`tkinter` 属于 Python 标准库，但部分 Linux 发行版需通过系统包管理器单独安装。例如 Ubuntu/Debian 可执行 `sudo apt install python3-tk`。

## 配置 CasJobs（可选）

若希望查询 PS1 Detection 表，请先创建 [MAST CasJobs](https://mastweb.stsci.edu/ps1casjobs/) 账户，并以环境变量提供凭据。请勿将用户名、密码或令牌写入脚本、`README` 或提交到仓库。

```bash
export CASJOBS_USERID="你的 CasJobs 用户名"
export CASJOBS_PW="你的 CasJobs 密码"
```

在 Windows PowerShell 中，可在当前会话设置为：

```powershell
$env:CASJOBS_USERID = "你的 CasJobs 用户名"
$env:CASJOBS_PW = "你的 CasJobs 密码"
```

如果未设置上述变量，脚本会尝试通过图形对话框请求输入。无图形界面的服务器环境建议始终使用环境变量。

## 使用方法

1. 打开 `Pan-STARRS1小行星数据回溯.py`，找到并修改目标编号，例如：

   ```python
   source = "2016 OR17"
   ```

2. 根据需要调整顶部配置项：

   ```python
   saveplots = True  # 是否保存图表
   save_data = True  # 是否保存数据、缩略图和原始 FITS
   use_jpl = True    # 是否通过 JPL Horizons 更新星历
   ```

3. 运行脚本：

   ```bash
   python "Pan-STARRS1小行星数据回溯.py"
   ```

脚本中的默认 MJD 查询范围为 `54985` 至 `57079`。如需其他时间段，请在 `cadc_ssos_query` 调用及绘图部分同步调整时间范围。若只希望快速验证查询流程，可将 `save_data` 设为 `False`，避免下载大量原始 FITS 文件。

## 输出说明

默认情况下，程序会创建类似 `PS1_data_2016 OR17/` 的目录。主要输出如下：

| 路径或文件 | 内容 |
|---|---|
| `sky_path.png` | 目标轨道路径与 PS1 观测位置图。 |
| `position-diffs.png` | 星历位置与目录检测位置的差异散点图；仅在 CasJobs 查询成功时生成。 |
| `images.png` | 观测图像的缩略图网格，可能包含检测位置、质量标记和估计星等。 |
| `observations_<目标>.csv` / `.fits` | CADC 返回并经位置、sky cell 等信息补充的观测表。 |
| `detections_<目标>.csv` / `.fits` | PS1 Detection 表的匹配结果；仅在 CasJobs 查询成功时生成。 |
| `thumbnails/` | 每个有效观测的 PNG 缩略图。 |
| `fits_data/` | 下载的原始 PS1 warp FITS 文件；文件数量和体积可能较大。 |

## 注意事项

脚本会访问多个第三方数据服务，运行时间、结果数量和下载体积取决于目标、查询时间范围与服务状态。部分查询或下载失败时，脚本会输出警告并尽可能继续执行；使用结果前应检查控制台日志与各输出表的完整性。

该项目面向科研数据检索与初步分析。使用、发布或解释结果前，请核验天体标识、时间尺度、星历精度、目录匹配条件及各数据服务的最新使用政策。图中星等由脚本中的 PSF flux 换算公式得到，应在科研使用前进行独立校验。

## 致谢与引用

使用 CADC 设施或 SSOS 服务开展研究时，请遵循其引用与致谢要求。脚本保留了以下建议致谢：

> This research used the facilities of the Canadian Astronomy Data Centre operated by the National Research Council of Canada with the support of the Canadian Space Agency.
>
> Please cite: Gwyn, Hill & Kavelaars (2012).

此外，请根据实际使用的数据产品、JPL Horizons、Pan-STARRS、MAST/CasJobs 以及所在机构或期刊的要求补充相应引用。 

## 依赖

完整的 Python 依赖列表见 [`requirements.txt`](requirements.txt)。

## 贡献

欢迎提交问题和改进建议。提交变更时，请避免上传 CasJobs 凭据、下载的原始 FITS、批量图像或其他大体积生成文件；这些内容已通过 `.gitignore` 排除。
