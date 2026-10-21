# Created 2026-10-20: redraw from saved annotations, not recovered historical code.
# The manuscript bars match C1 and T2, rather than three-sample group averages.
# Usage: Rscript /absolute/path/to/Code/Replot_from_source.R [output_directory]
suppressPackageStartupMessages(library(ggplot2))
s <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
panel_dir <- dirname(dirname(normalizePath(s)))
# 外部 Source data 根目录；保留私有资料库的子图目录结构。
# External source-data root; retain the private archive's panel directory names.
data_root <- Sys.getenv("SLC7A11_SOURCE_DATA_ROOT", unset="")
source_dir <- if(nzchar(data_root)) file.path(data_root, basename(panel_dir), "Source_data") else file.path(panel_dir, "Source_data")
if(!dir.exists(source_dir)) stop("Source data are not included in the public code repository. 请设置 SLC7A11_SOURCE_DATA_ROOT 指向已授权下载的私有资料库。")
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else file.path(panel_dir,'Reconstructed_figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(file.path(source_dir,'Figure5B_manuscript_C1_T2_source_data.csv'))
stopifnot(all(tapply(d$count,d$sample,sum)==tapply(d$denominator,d$sample,unique)))
d$sample <- factor(d$sample,levels=c('C1','T2'))
pal <- c(CD8_Tex='#ec2637',CD8_Tem='#52dacb',CD4_Tn='#31cb72',CD8_Proliferation='#2b4baf',
         CD4_Treg='#ef7b8b',CD8_Tcm='#7796c8',CD4_Th='#84efb8','mt-high_T'='#65573c',GammaDelta_T='#aaa797')
d$celltype <- factor(d$celltype,levels=names(pal))
p <- ggplot(d,aes(sample,proportion,fill=celltype)) + geom_col(width=.65) +
  scale_fill_manual(values=pal) + scale_y_continuous(limits=c(0,100),expand=c(0,0)) +
  theme_classic(base_size=11) + labs(x='C1 / T2 individual samples',y='T cell proportion (%)',
   fill='Saved T subtype',title='C1 / T2: manuscript-matching composition')
ggsave(file.path(out,'RECONSTRUCTED_C1_T2_from_source.pdf'),p,width=9,height=6)
ggsave(file.path(out,'RECONSTRUCTED_C1_T2_from_source.png'),p,width=9,height=6,dpi=150)
cat('REDRAW_OK: manuscript-matching C1/T2 composition\n')
