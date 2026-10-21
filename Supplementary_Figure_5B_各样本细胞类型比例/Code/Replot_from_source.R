# Created 2026-10-20: independent redraw from exported plotting tables.
# This is new recovery code, not the historical analysis script.
# Usage: Rscript /absolute/path/to/Code/Replot_from_source.R [output_directory]
suppressPackageStartupMessages(library(ggplot2))
file_arg <- grep("^--file=", commandArgs(), value=TRUE)
code_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
panel_dir <- dirname(code_dir)
# 外部 Source data 根目录；保留私有资料库的子图目录结构。
# External source-data root; retain the private archive's panel directory names.
data_root <- Sys.getenv("SLC7A11_SOURCE_DATA_ROOT", unset="")
source_dir <- if(nzchar(data_root)) file.path(data_root, basename(panel_dir), "Source_data") else file.path(panel_dir, "Source_data")
if(!dir.exists(source_dir)) stop("Source data are not included in the public code repository. 请设置 SLC7A11_SOURCE_DATA_ROOT 指向已授权下载的私有资料库。")
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else file.path(panel_dir,"Reconstructed_figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
src <- source_dir
d <- read.csv(file.path(src,"SupplementaryFigure5B_manuscript_selected_C1_T2_counts_proportions.csv"))
d$sample <- factor(d$sample,levels=c("C1","T2"))
p <- ggplot(d,aes(sample,proportion_percent,fill=sample)) + geom_col(width=.7) +
  facet_wrap(~celltype,scales="free_y",ncol=4) +
  scale_fill_manual(values=c(C1="#6D8EA9",T2="#AFD8CB")) +
  theme_classic(base_size=11) + guides(fill="none") +
  labs(x="Manuscript-visible sample",y="Percent of all cells in that sample",title="C1 / T2: individual samples")
w <- 12; h <- 9

ggsave(file.path(out,"RECONSTRUCTED_from_source.pdf"),p,width=w,height=h,limitsize=FALSE)
ggsave(file.path(out,"RECONSTRUCTED_from_source.png"),p,width=w,height=h,dpi=150,limitsize=FALSE)
cat("REDRAW_OK",normalizePath(out),"\n")
