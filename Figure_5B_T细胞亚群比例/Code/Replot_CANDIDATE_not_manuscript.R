# Created 2026-10-20: independent redraw from exported plotting tables.
# This is new recovery code, not the historical analysis script.
# Usage: Rscript /absolute/path/to/Code/Replot_CANDIDATE_not_manuscript.R [output_directory]
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
out <- if(length(args)) args[1] else file.path(panel_dir,"Reference","候选_3vs3均值_与稿图不一致")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
src <- source_dir
d <- read.csv(file.path(src,"候选_与稿图不一致","T_cells_group_mean_of_sample_proportions.csv"))
d$group <- factor(d$group,levels=c("Control","1A4_Antibody"))
p <- ggplot(d,aes(group,mean_sample_proportion_percent,fill=celltype)) + geom_col(width=.65) +
  theme_classic(base_size=11) +
  labs(x="Group",y="Mean of 3 within-sample proportions (%)",fill="Saved subtype",
       title="CANDIDATE ONLY: differs from manuscript Figure 5B")
w <- 8; h <- 6

ggsave(file.path(out,"RECONSTRUCTED_from_source.pdf"),p,width=w,height=h,limitsize=FALSE)
ggsave(file.path(out,"RECONSTRUCTED_from_source.png"),p,width=w,height=h,dpi=150,limitsize=FALSE)
cat("REDRAW_OK",normalizePath(out),"\n")
