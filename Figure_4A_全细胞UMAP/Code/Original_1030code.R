#读数据
library(Seurat)
library(harmony)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(patchwork)
library(patchwork)  # 用于图形组合
setwd("/data/public_database/Home_data/scRNAseq/1A4_anti/originaldata")
file_names <- list.files()
file_path<-paste("./",file_names,sep = "")
file_names<-gsub("_","-",file_names)
names(file_path)=file_names
#
file_path_name<-file_path
counts <- Read10X(data.dir = file_path_name)
scRNAsub = CreateSeuratObject(counts,
                              min.features = 200,
                              min.cells = 3)
unique(scRNAsub$orig.ident)
ncol(scRNAsub)
nrow(scRNAsub)
scRNAsub[["RNA3"]] <- as(object = scRNAsub[["RNA"]], Class = "Assay")
scRNAsub[["RNA"]]<-scRNAsub[["RNA3"]]
scRNAsub[["RNA3"]]<-NULL

##细胞周期
s_genes_mouse <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2",
                   "Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Mlf1ip","Hells","Rfc2",
                   "Rpa2","Nasps","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7",
                   "Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1",
                   "Tipin","Dscc1","Blm","Casp8ap2","Uhrf1bp1","H2afx","Cdc7","Ddx11")

g2m_genes_mouse <- c("Hmgb2","Cdc20","Gtse1","Top2a","Nusap1","Birc5","Tpx2",
                     "Cdc25c","Kif11","Anp32e","Tubb4b","Cks2","Nuf2","Ckap2l",
                     "Ckap2","Kif20b","Hjurp","Cdca3","Kif2c","Cks1b","Cdca8",
                     "Kif23","Ttk","Bub1","Cenpe","Kif18b","Cenpf","Kif4","Smc4",
                     "Ndc80","Ckap5","Cdca2","Cenpa")
#大群注释
scRNAsub[["percent.mt"]] <- PercentageFeatureSet(object = scRNAsub, pattern = "^mt-")#小鼠换成mt
VlnPlot(object = scRNAsub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,raster=FALSE)
scRNAsub<- subset(x = scRNAsub, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt <15)
scRNAsub <- NormalizeData(scRNAsub)
scRNAsub <- CellCycleScoring(scRNAsub,
                             s.features = s_genes_mouse,
                             g2m.features = g2m_genes_mouse,
                             set.ident = TRUE)
scRNAsub <- FindVariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub)
scRNAsub <- RunPCA(scRNAsub)
scRNAsub_harmony<-RunHarmony(scRNAsub,group.by.vars = "orig.ident",lambda=2)


scRNAsub_harmony <- scRNAsub_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:25) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:25) %>% 
  FindClusters(resolution = 1)  
gc()
##细胞注释
markers<-c(
  "Mki67","Top2a","Birc5",#Tumor cells
  "Col1a2", "Col3a1", "Dcn",#Fibroblasts
  "Ereg", "Krt18","Btc",#Epithelial
  "Cdh5", "Pecam1", "Esm1",#Endothelial
  "Ctss", "C1qa","C1qb",#Macrophage
  "Cd68","Adgre4","Treml4",#Monocyte
  "Cd3d", "Cd3e", "Cd3g",#T
  "Nkg7", "Gzma", "Gzmb",#NK
  "Csf3r","S100a9","S100a8",#Neutrophil
  "Clec9a", "Xcr1", "Cd24a",#DC
  "Cd79a", "Cd79b", "Ms4a1",#B cell
  "Cpa3", "Mcpt4", "Tpsb2",#Mast cells
  "Dclk1","Col1a1","Col5a1"#Tuft cells
)
#
Idents(scRNAsub_harmony)<-scRNAsub_harmony$seurat_clusters



DimPlot(
  object = scRNAsub_harmony, 
  reduction = "umap",  # 指定降维结果为UMAP（你的数据已运行RunUMAP）
  group.by = "seurat_clusters",  # 按聚类结果分组（未注释的原始簇）
  label = TRUE,  # 显示簇编号
  label.size = 4,  # 簇编号字体大小，可根据图调整
  repel = TRUE  # 避免簇编号重叠
) + 
  ggtitle("未注释的细胞分簇（UMAP）") +  # 添加图标题
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中
# 1. 先绘制UMAP分簇图，并赋值给变量p
p <- DimPlot(
  object = scRNAsub_harmony, 
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 4,
  repel = TRUE
) + 
  ggtitle("未注释的细胞分簇（UMAP）") +
  theme(plot.title = element_text(hjust = 0.5))

# 2. 保存图片到指定路径（以PNG格式为例，最常用）
ggsave(
  filename = "/data/public_database/Home_data/scRNAseq/1A4_anti/originaldata/umap_unannotated_clusters.png",  # 保存路径+文件名，建议和数据放同目录
  plot = p,  # 要保存的图（即上面赋值的p）
  width = 10,  # 图片宽度（单位：厘米，可根据需求调整）
  height = 8,  # 图片高度（单位：厘米）
  dpi = 300  # 分辨率（300dpi适合论文/报告，150dpi适合快速预览）
)


# 绘制单张气泡图（包含所有marker基因，不分面）
p_bubble <- DotPlot(
  object = scRNAsub_harmony,
  features = markers,              # 所有marker基因
  group.by = "seurat_clusters",    # 按未注释的簇分组
  cols = c("lightgrey", "red"),    # 表达量颜色梯度（浅灰→红）
  dot.scale = 8                    # 气泡大小缩放
) +
  RotatedAxis() +                  # 基因名垂直显示（避免重叠）
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题居中加粗
    axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5),     # 簇编号水平居中
    axis.text.y = element_text(size = 6),                              # 基因名大小
    legend.title = element_text(size = 12),                            # 图例标题
    legend.text = element_text(size = 10)                              # 图例文字
  ) +
  ggtitle("marker dotplot") +  # 图标题
  labs(
    x = "Marker",
    y = "Seurat Clusters"
  )

# 保存气泡图（调整尺寸适配所有基因）
ggsave(
  filename = "/data/public_database/Home_data/scRNAseq/1A4_anti/originaldata/marker_bubbleplot_all.pdf",
  plot = p_bubble,
  width = 14,  # 宽度足够展示28个簇
  height = 10, # 高度适配36个marker基因
  dpi = 300
)

##注释
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(0,2,4,6,11))) <- "Tumor cell"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(1,7,14,8,10))) <- "Fibroblast"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(3,12,9))) <- "Epithelial"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(27))) <- "Endothelial"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(5,16,18,23))) <- "Macrophage"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(15,19,21))) <- "T cell"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(17))) <- "NK"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(20))) <- "Neutrophil"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(24))) <- "DC"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(25))) <- "B cell"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(26))) <- "Mast cell"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(22))) <- "Tuft cell"
Idents(scRNAsub_harmony, cells = WhichCells(scRNAsub_harmony, idents = c(13))) <- "Monocyte"
scRNAsub_harmony$ALLTYPE<-as.character(scRNAsub_harmony@active.ident)
scRNAsub_harmony$ALLTYPE<-factor(scRNAsub_harmony$ALLTYPE,levels = c("Tumor cell","Fibroblast","Epithelial","Endothelial",
                                                                     "Macrophage","T cell","NK", "Neutrophil",
                                                                     "DC","B cell","Mast cell","Tuft cell","Monocyte"))


# 1. 绘制注释后的UMAP图（按细胞类型着色，显示标签）
p_annotated_umap <- DimPlot(
  object = scRNAsub_harmony,
  reduction = "umap",        # 基于UMAP降维结果
  group.by = "ALLTYPE",      # 按注释后的细胞类型分组
  label = TRUE,              # 显示细胞类型标签
  label.size = 3.5,          # 标签字体大小（避免过大遮挡）
  repel = TRUE,              # 避免标签重叠
  cols = NULL                # 自动分配颜色（也可自定义，如cols = scales::hue_pal()(13)）
) +
  ggtitle("注释后的细胞分簇UMAP图") +  # 图标题
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题居中加粗
    legend.title = element_text(size = 12, face = "bold"),             # 图例标题加粗
    legend.text = element_text(size = 10)                              # 图例文字大小
  )

# 2. 保存UMAP图（PNG+PDF双格式，适配不同场景）
# PNG格式（适合快速查看、插入PPT）
ggsave(
  filename = "/data/public_database/Home_data/scRNAseq/1A4_anti/annotated_umap.png",
  plot = p_annotated_umap,
  width = 12,
  height = 10,
  dpi = 300
)

# PDF格式（矢量图，适合论文排版、放大不失真）
ggsave(
  filename = "/data/public_database/Home_data/scRNAseq/1A4_anti/originaldata/annotated_umap.pdf",
  plot = p_annotated_umap,
  width = 12,
  height = 10
)

# 3. 显示UMAP图
print(p_annotated_umap)


# 保存注释后的Seurat对象（后缀为.rds，R专用格式）
saveRDS(
  object = scRNAsub_harmony,
  file = "/data/public_database/Home_data/scRNAseq/1A4_anti/originaldata/annotated_scRNAsub_harmony.rds"
)

