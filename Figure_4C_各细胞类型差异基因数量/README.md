# Figure_4C — 各细胞类型差异基因数量 / Differentially expressed gene counts

## 中文

新编 R 重绘代码读取已导出的图级数据，不重新运行上游单细胞分析。 只恢复差异基因数量汇总，未恢复基因级结果和统计参数。

本公开目录仅包含代码和说明。绘图 Source data、结果图及核验材料暂存[私有资料库](https://github.com/chymlltrn/SLC7A11_JICHANG_SourceData)，稍后公开，目前需由通讯作者安排授权访问。源数据应保持 `Figure_4C_各细胞类型差异基因数量/Source_data/` 结构。

## English

New R redraw code reads exported plot-level data without rerunning upstream single-cell analysis. Only DEG-count summaries were recovered, not gene-level results or historical statistical parameters.

This public directory contains code and documentation only. Plot-level source data, figures, and audit materials are retained in the [private archive](https://github.com/chymlltrn/SLC7A11_JICHANG_SourceData) for later public release. Authorized access is currently arranged through the corresponding author. Retain the `Figure_4C_各细胞类型差异基因数量/Source_data/` layout.

```bash
python3 scripts/reproduce_figures.py --panel Figure_4C --data-root "/absolute/path/to/SLC7A11_JICHANG_SourceData" --rscript /absolute/path/to/Rscript
```

`Original_*`：历史来源脚本 / historical provenance scripts. `Recovery_export*`, `export_*`, `extract_*`：提取脚本，依赖原始对象或历史 PDF / extraction scripts requiring original objects or historical PDFs. 默认入口不会执行这些脚本 / the launcher does not execute these scripts.
