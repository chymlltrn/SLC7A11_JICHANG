# Read existing UMAP/metadata only; no fitting, expression export or new embedding.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(Seurat))
root<-"/data/public_database/Home_data/scRNAseq/1A4_anti/1A4_scRNA"
out<-"/tmp/slc7a11_velocity_export_20261020"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
for(kind in c("Macrophages","T_cells")) {
 p<-file.path(root,paste0(kind,"_final_annotated.rds"))
 s<-readRDS(p)
 stopifnot("umap" %in% names(s@reductions))
 e<-Embeddings(s,"umap")[Cells(s),,drop=FALSE]
 m<-s@meta.data[Cells(s),,drop=FALSE]
 stopifnot(identical(rownames(e),rownames(m)))
 g<-ifelse(m$old.ident %in% paste0("C",1:3),"Control",ifelse(m$old.ident %in% paste0("T",1:3),"1A4_Antibody",NA_character_))
 labels<-as.character(m$celltype)
 if(kind=="T_cells")labels[labels=="γδ_T"]<-"gammadelta_T"
 z<-data.frame(cell_id=rownames(e),UMAP1=e[,1],UMAP2=e[,2],celltype=labels,sample=as.character(m$old.ident),group=g,row.names=NULL)
 write.csv(z,file.path(out,paste0(kind,"_existing_umap_metadata.csv")),row.names=FALSE)
 writeLines(p,file.path(out,paste0(kind,"_existing_umap_source_path.txt")))
 cat(kind,nrow(z),"\n");print(table(z$group,useNA="always"))
 rm(s,m,e,z);gc()
}
writeLines(capture.output(sessionInfo()),file.path(out,"UMAP_export_sessionInfo.txt"))
