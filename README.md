# SLC7A11_JICHANG

[中文](#中文) · [English](#english)

## 中文

SLC7A11 单细胞图的公开代码库，按子图整理历史脚本、提取代码和新编重绘代码。代码已公开；绘图 Source data、结果图及完整核验材料暂保留于 [SLC7A11_JICHANG_SourceData 私有资料库](https://github.com/chymlltrn/SLC7A11_JICHANG_SourceData)，计划稍后公开，目前可由通讯作者安排编辑和审稿人授权访问。私有链接本身不授予访问权限。

共登记 18 个子图，其中 13 个具有匹配的图级来源及重绘入口；Figure 4B、4E、4F、5G、5H 的最终来源仍有缺项。本库没有包含 Source data、结果图、单细胞原始数据或含数据的旧提交历史。

先核对代码完整性；获得授权数据后，将私有资料库作为 `--data-root`：

```bash
python3 scripts/verify_repository.py
python3 scripts/reproduce_figures.py --set all --dry-run
python3 scripts/reproduce_figures.py --set core --data-root "/absolute/path/to/SLC7A11_JICHANG_SourceData" --rscript /absolute/path/to/Rscript
python3 scripts/reproduce_figures.py --set cellchat --data-root "/absolute/path/to/SLC7A11_JICHANG_SourceData" --rscript /absolute/path/to/Rscript
```

`core` 需要 Python 3、R 和 ggplot2；CellChat 子图另需 CellChat 2.1.2。数据目录名称与子图目录应保持一致。也可设置环境变量 `SLC7A11_SOURCE_DATA_ROOT`。尚无数据时可检查代码和生成命令，但不能完成重绘。输出默认写入已忽略的 `validation_outputs/`。

历史 `Original_*` 脚本保留原路径和依赖，供来源核查；`Replot_*`、CellChat recovery 和 velocity replay 是新编重绘代码。默认入口不执行上游分析或导出脚本。Figure 5B 对应 C1/T2 单一样本，并非组均值，原题目 Average 需作者核对。Velocity replay 恢复历史图形，不重新推断速度；CellChat 重绘不能补齐缺失的历史推断和筛选规则。

## English

Public code for the SLC7A11 single-cell figure panels, organized by panel into historical scripts, extraction code, and newly written redraw code. Code is public. Plot-level source data, figures, and full audit materials remain in the [private SLC7A11_JICHANG_SourceData archive](https://github.com/chymlltrn/SLC7A11_JICHANG_SourceData) for later public release. Editors and reviewers may currently obtain authorized access through the corresponding author. A private URL alone does not grant access.

The inventory covers 18 panels: 13 have matching plot-level sources and redraw entry points; final source materials for Figure 4B, 4E, 4F, 5G, and 5H remain incomplete. This public repository excludes source data, result figures, raw single-cell data, and the old data-bearing commit history.

Use the commands above to verify code integrity and inspect planned commands. Once authorized data are available locally, supply the private archive root with `--data-root`, or set `SLC7A11_SOURCE_DATA_ROOT`. Retain the panel directory names. The core set requires Python 3, R, and ggplot2; CellChat panels additionally require CellChat 2.1.2. Redraws cannot run without their source data. Outputs default to ignored `validation_outputs/`.

`Original_*` scripts retain historical paths and dependencies for provenance. `Replot_*`, CellChat recovery, and velocity replay are newly written visualization code. The launcher does not execute upstream analysis or extraction scripts. Figure 5B matches individual C1/T2 samples, not group means; its original Average title needs author review. Velocity replay reconstructs the historical graphic without refitting velocity. CellChat redraws do not recover missing historical inference or selection rules.

## 子图目录 / Panel inventory

| 子图 / Panel | 内容 / Description | 状态 / Status |
|---|---|---|
| [Figure_4A](Figure_4A_全细胞UMAP/README.md) | 全细胞UMAP / Whole-cell UMAP | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Figure_4B](Figure_4B_细胞类型富集log2OR/README.md) | 细胞类型富集log2OR / Cell-type enrichment (log2 odds ratio) | 缺失 / Unavailable |
| [Figure_4C](Figure_4C_各细胞类型差异基因数量/README.md) | 各细胞类型差异基因数量 / Differentially expressed gene counts | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Figure_4D](Figure_4D_CellChat差异信号通路/README.md) | CellChat差异信号通路 / CellChat pathway comparison | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Figure_4E](Figure_4E_巨噬细胞M1_M2评分/README.md) | 巨噬细胞M1_M2评分 / Macrophage M1/M2 scores | 缺失 / Unavailable |
| [Figure_4F](Figure_4F_巨噬细胞KEGG/README.md) | 巨噬细胞KEGG / Macrophage KEGG enrichment | 缺失 / Unavailable |
| [Figure_4G](Figure_4G_CellChat通讯强度圈图/README.md) | CellChat通讯强度圈图 / CellChat interaction-strength circles | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Figure_4H](Figure_4H_Tumor到Macrophage信号/README.md) | Tumor到Macrophage信号 / Tumor-to-macrophage signaling | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Figure_4I](Figure_4I_Macrophage到Tumor信号/README.md) | Macrophage到Tumor信号 / Macrophage-to-tumor signaling | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Figure_4J](Figure_4J_巨噬细胞亚群UMAP/README.md) | 巨噬细胞亚群UMAP / Macrophage-subtype UMAP | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Figure_4K](Figure_4K_巨噬细胞RNA_velocity/README.md) | 巨噬细胞RNA_velocity / Macrophage RNA velocity | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Figure_5A](Figure_5A_T细胞亚群UMAP/README.md) | T细胞亚群UMAP / T-cell-subtype UMAP | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Figure_5B](Figure_5B_T细胞亚群比例/README.md) | T细胞亚群比例 / T-cell-subtype proportions | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Figure_5G](Figure_5G_GO富集/README.md) | GO富集 / GO enrichment | 缺失 / Unavailable |
| [Figure_5H](Figure_5H_KEGG富集/README.md) | KEGG富集 / KEGG enrichment | 缺失 / Unavailable |
| [Figure_5I](Figure_5I_T细胞RNA_velocity/README.md) | T细胞RNA_velocity / T-cell RNA velocity | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Supplementary_Figure_5A](Supplementary_Figure_5A_细胞标记基因气泡图/README.md) | 细胞标记基因气泡图 / Cell-marker dot plot | 重绘代码可用；数据暂私有 / Redraw code available; data private |
| [Supplementary_Figure_5B](Supplementary_Figure_5B_各样本细胞类型比例/README.md) | 各样本细胞类型比例 / Per-sample cell-type proportions | 重绘代码可用；数据暂私有 / Redraw code available; data private |

## 进一步说明 / Further documentation

- [重绘方法、环境和限制 / Reproducibility, environment, and limits](REPRODUCIBILITY.md)
- [数据和代码可用性 / Data and code availability](DATA_AVAILABILITY.md)
- [代码来源与修改记录 / Code provenance and modifications](CODE_PROVENANCE.tsv)

`FILE_MANIFEST.tsv` 仅记录本公开库文件的 SHA256；不能证明数据完整或上游分析已复现。 / `FILE_MANIFEST.tsv` records SHA256 hashes of public code files; it does not establish data completeness or upstream-analysis reproducibility.

历史脚本及方法参数保留原样，注释保留原文；导航、用法、子图说明与限制均提供中英对照。 / Historical scripts and method parameters retain their original wording; navigation, usage, panel notes, and limitations are bilingual.
