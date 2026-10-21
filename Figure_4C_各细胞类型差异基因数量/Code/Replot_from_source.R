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
d <- read.csv(file.path(src,"DEG_statistics.csv"))
z <- rbind(data.frame(Cell_Type=d$Cell_Type,direction="Up_DEGs",count=d$Up_DEGs),
           data.frame(Cell_Type=d$Cell_Type,direction="Down_DEGs",count=d$Down_DEGs))
z$Cell_Type <- factor(z$Cell_Type,levels=rev(d$Cell_Type))
p <- ggplot(z,aes(direction,Cell_Type,fill=count)) + geom_tile() + geom_text(aes(label=count),size=3.8) +
  scale_fill_gradient(low="white",high="#007F81") + theme_minimal(base_size=12) +
  labs(x=NULL,y=NULL,fill="DEG count",title="Recovered DEG counts")
w <- 6; h <- 7

ggsave(file.path(out,"RECONSTRUCTED_from_source.pdf"),p,width=w,height=h,limitsize=FALSE)
ggsave(file.path(out,"RECONSTRUCTED_from_source.png"),p,width=w,height=h,dpi=150,limitsize=FALSE)
cat("REDRAW_OK",normalizePath(out),"\n")
