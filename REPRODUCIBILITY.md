# 重绘与可复现性 / Reproducibility

## 中文

公开库的独立检查：`python3 scripts/verify_repository.py`。该检查验证代码清单、哈希、Python 语法、子图配置和禁入数据文件，不读取私有数据，也不声称完成科学结果复现。

获得授权的私有资料库后，在公开代码库根目录运行以下命令（路径有空格时加引号）。从其他工作目录启动时，使用 `scripts/reproduce_figures.py` 的绝对路径：

```bash
python3 scripts/reproduce_figures.py --set core --data-root "/absolute/path/to/SLC7A11_JICHANG_SourceData" --rscript /absolute/path/to/Rscript --out "/absolute/path/to/output"
python3 scripts/reproduce_figures.py --set cellchat --data-root "/absolute/path/to/SLC7A11_JICHANG_SourceData" --rscript /absolute/path/to/Rscript --out "/absolute/path/to/output"
python3 scripts/reproduce_figures.py --panel Figure_5B --data-root "/absolute/path/to/SLC7A11_JICHANG_SourceData" --rscript /absolute/path/to/Rscript
python3 scripts/reproduce_figures.py --set all --dry-run
```

也可用 `SLC7A11_SOURCE_DATA_ROOT` 环境变量指定私有资料库根目录；各子图仍使用 `<子图目录>/Source_data/`。直接运行 `Replot_*.R` 时也应设置该环境变量。`--dry-run` 不读取数据、不要求 R 已安装。实际重绘会先检查输入是否齐全，再执行并写入 stdout、stderr 和 `run_status.json`；失败返回非零状态。

已测试环境：核心 R 重绘为 R 4.6.0、ggplot2 4.0.3（macOS）；CellChat 为 R 4.5.1、CellChat 2.1.2、ggplot2（Linux）；velocity 矢量重放仅需 Python 3 标准库。可选提取脚本另依赖 Seurat 或 PyMuPDF/NumPy/pandas；默认重绘不执行提取。未提供完整历史环境锁定文件，不自动安装依赖。

7 个 R 入口和 2 个 velocity 子图构成 core（两组 velocity 分别执行），4 个子图构成 CellChat。5 个缺失子图被明确拒绝，不生成推测结果。Figure 4C 只有 DEG 数量汇总；Figure 5B 为 C1/T2 单一样本组成，不是组均值。UMAP 使用已保存坐标，历史筛选及最终注释流程未完全恢复；velocity 重放不拟合模型；CellChat 的历史筛选/样式缺口仍在。重绘的字体、颜色、点大小及布局可与稿图不同。

## English

Run `python3 scripts/verify_repository.py` for a data-independent check of the public file inventory, hashes, Python syntax, panel configuration, and exclusion of data files. This does not read private data or establish scientific reproducibility.

After obtaining an authorized copy of the private archive, run the commands above from the public code repository root. From another working directory, use the absolute path to `scripts/reproduce_figures.py`. Quote paths containing spaces. `--data-root` or `SLC7A11_SOURCE_DATA_ROOT` points to the private archive root, preserving `<panel-directory>/Source_data/`. Set that environment variable when invoking a `Replot_*.R` script directly. `--dry-run` needs neither data nor R. Actual redraws preflight input files, then record stdout, stderr, and `run_status.json`; failures return a nonzero status.

Tested environments: core R redraws use R 4.6.0 and ggplot2 4.0.3 on macOS; CellChat uses R 4.5.1, CellChat 2.1.2, and ggplot2 on Linux; velocity vector replay uses only the Python 3 standard library. Optional extraction scripts require Seurat or PyMuPDF/NumPy/pandas and are not run by the launcher. No complete historical lockfile is claimed; dependencies are not installed automatically.

The core set comprises seven R entry points and two velocity panels (two group commands per velocity panel); four panels form the CellChat set. Five incomplete panels are explicitly rejected. Figure 4C contains DEG-count summaries only. Figure 5B uses individual C1/T2 samples rather than group means. UMAP uses saved coordinates; historical filtering and final annotation were not fully recovered. Velocity replay does not refit a model, and CellChat historical selection/style gaps remain. Redrawn fonts, colors, point sizes, and layout may differ from the manuscript.
