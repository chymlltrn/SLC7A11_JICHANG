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
d <- read.csv(file.path(src,"Macrophages_UMAP_plot_source.csv"),check.names=FALSE)
stopifnot(!anyDuplicated(d$cell_id),all(is.finite(d$saved_umap_1)),all(is.finite(d$saved_umap_2)))
d$label <- d[["celltype"]]
centers <- aggregate(d[,c("saved_umap_1","saved_umap_2")],list(label=d$label),median)
p <- ggplot(d,aes(saved_umap_1,saved_umap_2,color=label)) +
  geom_point(size=.25,alpha=.75,shape=16) +
  geom_text(data=centers,aes(label=label),color="black",size=3,check_overlap=TRUE) +
  coord_equal() + theme_classic(base_size=11) +
  labs(x="UMAP1",y="UMAP2",color="Saved label",title="Saved UMAP: independent redraw") +
  guides(color=guide_legend(override.aes=list(size=3,alpha=1)))
w <- 10; h <- 7

ggsave(file.path(out,"RECONSTRUCTED_from_source.pdf"),p,width=w,height=h,limitsize=FALSE)
ggsave(file.path(out,"RECONSTRUCTED_from_source.png"),p,width=w,height=h,dpi=150,limitsize=FALSE)
cat("REDRAW_OK",normalizePath(out),"\n")
