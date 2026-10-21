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
d <- read.csv(file.path(src,"SupplementaryFigure5A_original_DotPlot_layer_data.csv"),check.names=FALSE)
markers <- read.csv(file.path(src,"SupplementaryFigure5A_original_marker_list.csv"))
d$features.plot <- factor(d$features.plot,levels=markers$gene)
d$id <- factor(d$id,levels=sort(unique(as.integer(as.character(d$id)))))
p <- ggplot(d,aes(features.plot,id)) + geom_point(aes(size=pct.exp,color=avg.exp.scaled)) +
  scale_size(range=c(0,8),limits=c(0,100),name="Percent expressed") +
  scale_color_gradient(low="lightgrey",high="red",name="Average expression") +
  theme_classic(base_size=10) + theme(axis.text.x=element_text(size=6,angle=90,hjust=1,vjust=.5)) +
  labs(x="Marker",y="Seurat cluster",title="Marker dotplot: saved data / original DotPlot call")
w <- 14; h <- 9

ggsave(file.path(out,"RECONSTRUCTED_from_source.pdf"),p,width=w,height=h,limitsize=FALSE)
ggsave(file.path(out,"RECONSTRUCTED_from_source.png"),p,width=w,height=h,dpi=150,limitsize=FALSE)
cat("REDRAW_OK",normalizePath(out),"\n")
