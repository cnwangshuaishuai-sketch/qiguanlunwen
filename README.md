# organ-crosstalk-eclinm

从常规体检数据重建时序性器官串扰网络 — 论文代码

## 论文信息

- **标题**: A Temporal Organ Crosstalk Network Reconstructed from Routine Health Screening Data
- **期刊**: eClinicalMedicine (The Lancet)
- **作者**: 王帅
- **单位**: 大庆铭德医院

## 文件说明

| 文件 | 用途 |
|------|------|
| `crp_sensitivity.py` | hs-CRP缺失敏感性分析三合一（MAR检验 + IPW逆概率加权 + Tipping Point分析） |
| `run_missing_dock.py` | 分子对接辅助脚本 |

## Methods 概览

论文使用28项方法，覆盖：

- **数据工程**: DuckDB + OMOP CDM v5.4 标准化入库，正则提取hs-CRP
- **预处理**: 个体内z-score标准化，log转换，±5SD异常值剔除
- **核心模型**: 交叉滞后面板模型(CLPM)，30条有向边，Bonferroni校正
- **模型诊断**: HC2稳健标准误，VIF共线性，Durbin-Watson自相关
- **敏感性分析**: 替代指标验证，完整CRP亚组，IPW，Tipping Point
- **亚组**: 性别分层，年龄中位数

## 数据可用性

去标识化的个体级数据可根据合理请求从通讯作者处获取。
