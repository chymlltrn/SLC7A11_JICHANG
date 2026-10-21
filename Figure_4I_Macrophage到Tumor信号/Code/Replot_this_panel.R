# New recovery wrapper; historical plotting code was not recovered.
s <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
panel_dir <- dirname(dirname(normalizePath(s)))
# 外部 Source data 根目录；保留私有资料库的子图目录结构。
# External source-data root; retain the private archive's panel directory names.
data_root <- Sys.getenv("SLC7A11_SOURCE_DATA_ROOT", unset="")
source_dir <- if(nzchar(data_root)) file.path(data_root, basename(panel_dir), "Source_data") else file.path(panel_dir, "Source_data")
if(!dir.exists(source_dir)) stop("Source data are not included in the public code repository. 请设置 SLC7A11_SOURCE_DATA_ROOT 指向已授权下载的私有资料库。")
rscript <- file.path(R.home('bin'),'Rscript')
args <- c(shQuote(file.path(panel_dir,'Code','Recovery_replot_from_aggregate_networks.R')),
          shQuote(file.path(source_dir,'aggregate_networks_NO_SINGLE_CELL_EXPRESSION.rds')),
          shQuote(file.path(panel_dir,'Reconstructed_figures')), 'I')
status <- system2(rscript,args)
if(status != 0) stop('CellChat redraw failed; use the recorded server CellChat 2.1.2 environment.')
