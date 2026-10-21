# Figure_4H — Tumor到Macrophage信号 / Tumor-to-macrophage signaling

## 中文

新编重绘代码读取保存的汇总通讯网络。历史通讯推断过程未完整恢复；4D/4I 的历史通路筛选规则及 4G 的原始样式参数仍有缺口。

本公开目录仅包含代码和说明。绘图 Source data、结果图及核验材料暂存[私有资料库](https://github.com/chymlltrn/SLC7A11_JICHANG_SourceData)，稍后公开，目前需由通讯作者安排授权访问。源数据应保持 `Figure_4H_Tumor到Macrophage信号/Source_data/` 结构。

## English

New redraw code reads saved aggregate communication networks. Historical inference was not fully recovered; historical pathway-selection rules for 4D/4I and original circle styling for 4G remain unresolved.

This public directory contains code and documentation only. Plot-level source data, figures, and audit materials are retained in the [private archive](https://github.com/chymlltrn/SLC7A11_JICHANG_SourceData) for later public release. Authorized access is currently arranged through the corresponding author. Retain the `Figure_4H_Tumor到Macrophage信号/Source_data/` layout.

```bash
python3 scripts/reproduce_figures.py --panel Figure_4H --data-root "/absolute/path/to/SLC7A11_JICHANG_SourceData" --rscript /absolute/path/to/Rscript
```

`Original_*`：历史来源脚本 / historical provenance scripts. `Recovery_export*`, `export_*`, `extract_*`：提取脚本，依赖原始对象或历史 PDF / extraction scripts requiring original objects or historical PDFs. 默认入口不会执行这些脚本 / the launcher does not execute these scripts.
