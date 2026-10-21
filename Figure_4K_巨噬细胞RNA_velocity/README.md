# Figure_4K — 巨噬细胞RNA_velocity / Macrophage RNA velocity

## 中文

已恢复图级矢量绘图来源，使用 Python 标准库重绘历史图形；不重新拟合 RNA velocity，不代表原生细胞速度或转移矩阵。历史 velocyto.R 脚本仅供来源核查。

本公开目录仅包含代码和说明。绘图 Source data、结果图及核验材料暂存[私有资料库](https://github.com/chymlltrn/SLC7A11_JICHANG_SourceData)，稍后公开，目前需由通讯作者安排授权访问。源数据应保持 `Figure_4K_巨噬细胞RNA_velocity/Source_data/` 结构。

## English

Recovered plot-level vectors support Python-standard-library replay of the historical graphic. This does not refit RNA velocity or provide native cell velocity/transition matrices. Historical velocyto.R scripts are provided for provenance.

This public directory contains code and documentation only. Plot-level source data, figures, and audit materials are retained in the [private archive](https://github.com/chymlltrn/SLC7A11_JICHANG_SourceData) for later public release. Authorized access is currently arranged through the corresponding author. Retain the `Figure_4K_巨噬细胞RNA_velocity/Source_data/` layout.

```bash
python3 scripts/reproduce_figures.py --panel Figure_4K --data-root "/absolute/path/to/SLC7A11_JICHANG_SourceData" --rscript /absolute/path/to/Rscript
```

`Original_*`：历史来源脚本 / historical provenance scripts. `Recovery_export*`, `export_*`, `extract_*`：提取脚本，依赖原始对象或历史 PDF / extraction scripts requiring original objects or historical PDFs. 默认入口不会执行这些脚本 / the launcher does not execute these scripts.
