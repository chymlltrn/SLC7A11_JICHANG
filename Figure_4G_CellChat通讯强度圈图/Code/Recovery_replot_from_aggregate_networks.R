# Newly written diagnostic reconstruction code (2026-10-20), not recovered historical plotting code.
# Usage: Rscript replot_from_aggregate_networks.R aggregate_networks_NO_SINGLE_CELL_EXPRESSION.rds output_dir [D|G|H|I|all]
# Input contains only aggregate networks and cell-type counts; no expression matrices or single-cell identifiers.
suppressPackageStartupMessages(library(CellChat))
suppressPackageStartupMessages(library(ggplot2))
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<2) stop('Provide aggregate RDS path and output directory.')
s<-readRDS(a[1]);out<-a[2];panel<-if(length(a)>=3)a[3] else 'all';dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(identical(s$contains_single_cell_expression,FALSE))
m<-new('CellChat');m@netP<-setNames(lapply(s$groups,function(z)z$netP),s$group_label_mapping$figure_group);m@options<-list(mode='merged')
redraw_rank<-function(tag,paths=NULL,src=NULL,tgt=NULL,colors,title) {
 z<-rankNet(m,mode='comparison',slot.name='netP',measure='weight',comparison=c(1,2),stacked=TRUE,signaling=paths,sources.use=src,targets.use=tgt,return.data=TRUE,do.stat=FALSE,color.use=colors,title=title,font.size=10)
 write.csv(z$signaling.contribution,file.path(out,paste0(tag,'_rankNet_return.csv')),row.names=FALSE)
 ggsave(file.path(out,paste0(tag,'_RECONSTRUCTED.png')),z$gg.obj,width=7,height=8,dpi=180)
}
if(panel %in% c('D','all')) {
 # These 20 labels are transcribed from the manuscript. They also equal Control's top 20 raw pathway flows.
 # This agreement is evidence for the subset; the historical selection code was not recovered.
 d<-c('TENASCIN','CD34','SELE','LAMININ','SEMA3','SEMA6','EPHA','FGF','FN1','SPP1','COLLAGEN','THBS','CSF','GALECTIN','MIF','CD45','APP','CADM','MHC-I','CCL')
 redraw_rank('F4D_target20',d,colors=c('#2d81db','#e42b36'),title='Top 20 pathways: reconstruction')
}
if(panel %in% c('H','all')) redraw_rank('F4H',src='Tumor cell',tgt='Macrophage',colors=c('#fb5483','#4afbae'),title='Tumor - Macrophage: reconstruction')
if(panel %in% c('I','all')) {
 # The raster shows 14 labels. Default rankNet also returns SEMA3, CD86 and LAMININ.
 # The historical reason for excluding those three was not recovered.
 i<-c('MIF','MHC-I','GALECTIN','VISTA','CSF','CD45','VISFATIN','VCAM','CADM','IL2','ICAM','SN','CD34','APP')
 redraw_rank('F4I_target14',i,src='Macrophage',tgt='Tumor cell',colors=c('#fb5483','#4afbae'),title='Macrophage - Tumor: reconstruction')
}
if(panel %in% c('G','all')) {
 # Diagnostic parameters, not asserted historical parameters. Edge transparency, max width,
 # node sizes, color palette and scaling may differ from the manuscript raster.
 png(file.path(out,'F4G_RECONSTRUCTED.png'),width=2400,height=1200,res=160)
 par(mfrow=c(1,2),xpd=TRUE)
 for(j in seq_along(s$groups)) {
  z<-s$groups[[j]];netVisual_circle(z$net$weight,vertex.weight=z$nodes$n_cells,weight.scale=TRUE,label.edge=FALSE,title.name=paste('Interaction strength -',s$group_label_mapping$figure_group[j]))
 }
 dev.off()
}
writeLines(capture.output(sessionInfo()),file.path(out,'replot_sessionInfo.txt'))
