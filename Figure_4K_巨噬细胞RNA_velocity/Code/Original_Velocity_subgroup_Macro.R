## =========================================================
## RNA velocity（仅基于已有巨噬细胞Seurat对象 + 现有loom）
## - 读取已做完聚类与亚群注释的RDS（仅巨噬细胞）
## - 与loom对齐：样本前缀 + 条形码（C1:BARCODEx ↔ C1_BARCODE-1）
## - 估计velocity；输出：按亚群（cluster）着色、（可选）按分组着色/分组内独立拟合
## - 不做重聚类；UMAP若缺失则临时计算
## =========================================================
options(stringsAsFactors = FALSE)
set.seed(2025)

suppressPackageStartupMessages({
  library(Seurat); library(velocyto.R); library(Matrix)
  library(RColorBrewer); library(scales); library(dplyr); library(glue); library(ggplot2)
})

## ---------- 0) 路径 & 参数 ----------
loom_dir   <- "/data/public_database/Home_data/scRNAseq/1A4_anti/1A4_scRNA/loom_unzip"
seurat_rds <- "/data/public_database/Home_data/scRNAseq/1A4_anti/1A4_scRNA/Macrophages_final_annotated.rds"
outdir     <- "/data/public_database/Home_data/scRNAseq/1A4_anti/1A4_scRNA/velocity_macrophages_251104"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## （占位）你的“亚群列名”（若自动侦测失败请改成正确列名）
idents_col <- "celltype"

## 你数据中的“分组列名”（存在则启用相应绘图/计算；不存在会自动跳过）
group_col  <- "group"   # 期望可能是 "1A4_Antibody" / "Control"

## velocity 拟合参数（沿用你之前风格）
deltaT <- 1; kCells <- 30; fit.quantile <- 0.02; n.cores <- 16

## 是否写回loom（本次默认关闭）
write_back_loom <- FALSE

## ---------- 1) 收集loom ----------
loom_files <- list.files(loom_dir, pattern = "\\.loom$", full.names = TRUE, recursive = TRUE)
stopifnot(length(loom_files) > 0)

## ---------- 2) 读Seurat对象（已是巨噬细胞子集 & 已注释） ----------
message("[*] Loading Seurat object (Macrophages only) ...")
seu <- readRDS(seurat_rds)
## 定义映射
ctrl_levels <- paste0("C", 1:3)  # C1–C3
anti_levels <- paste0("T", 1:3)  # T1–T3
## 生成新的分组变量 group
seu$group <- dplyr::case_when(
  seu$old.ident %in% ctrl_levels ~ "Control",
  seu$old.ident %in% anti_levels ~ "1A4_Antibody",
  TRUE ~ NA_character_
)
seu$group <- factor(seu$group, levels = c("Control", "1A4_Antibody"))
## 快速核对
table(seu$group, useNA = "ifany")
addmargins(table(seu$old.ident, seu$group, useNA = "ifany"))

## ---------- 2a) 设定Idents（自动侦测 + 占位回退） ----------
subcluster_candidates <- c(
  idents_col,                    # 你手动指定/占位的优先
  "macro_subcluster", "subcluster", "Macrophage_subtype",
  "annot", "annotation", "Cluster", "clusters",
  "celltype_sub", "celltype_detail",
  "seurat_clusters"              # 万能兜底
)
valid_cols <- intersect(subcluster_candidates, colnames(seu@meta.data))
if (length(valid_cols) > 0) {
  idents_col <- valid_cols[1]
  message(glue("[*] Using '{idents_col}' as Idents."))
  Idents(seu) <- idents_col
} else {
  message("[*] No suitable subcluster column found; keep current Idents(seu).")
}

## ---------- 2b) 确保有UMAP（若无，则基于现有Assay/PCA临时计算一份） ----------
has_umap <- "umap" %in% names(seu@reductions)
if (!has_umap) {
  message("[*] UMAP not found in object; computing a temporary UMAP...")
  DefaultAssay(seu) <- DefaultAssay(seu) %||% "RNA"
  if (!"pca" %in% names(seu@reductions)) {
    seu <- NormalizeData(seu, verbose = FALSE)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    seu <- ScaleData(seu, verbose = FALSE)
    seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  }
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
}
emb_umap <- Embeddings(seu, "umap")[Cells(seu), , drop = FALSE]

## ---------- 2c) 分组列存在性检测 ----------
has_group <- group_col %in% colnames(seu@meta.data)
if (!has_group) {
  message(glue("[*] NOTE: column '{group_col}' not found; group-based plotting will be skipped."))
}

## ---------- 3) 样本/条形码解析与映射工具 ----------
extract_sample <- function(v) {
  has_sep <- grepl("[:_]", v)
  out <- ifelse(has_sep, sub("[:_].*$", "", v), NA_character_)
  out
}
extract_barcode <- function(v) {
  core <- sub("^.*[:_]", "", v)  # 去前缀
  core <- sub("x$", "", core)    # 10x loom尾部 x
  core <- sub("-1$", "", core)   # Seurat尾部 -1
  core <- sub("\\.1$", "", core) # 某些流程 .1
  core
}
make_keys <- function(vec) {
  data.frame(sample = extract_sample(vec),
             bc     = extract_barcode(vec),
             stringsAsFactors = FALSE)
}

## ---------- 4) 读入各loom并与Seurat映射 ----------
emat_list <- list(); nmat_list <- list()
for (f in loom_files) {
  sid <- tools::file_path_sans_ext(basename(f))
  message(sprintf("[*] Reading loom: %s", sid))
  ld <- velocyto.R::read.loom.matrices(f)  # list(spliced, unspliced, ambiguous)
  em <- ld$spliced; nu <- ld$unspliced; rm(ld)
  
  cand <- colnames(seu)
  key_loom <- make_keys(colnames(em))
  key_cand <- make_keys(cand)
  
  s1 <- paste(key_loom$sample, key_loom$bc, sep = "|")
  s2 <- paste(key_cand$sample, key_cand$bc, sep = "|")
  m  <- match(s1, s2)
  
  need_fallback <- which(is.na(m))
  if (length(need_fallback) > 0) {
    bc_cand <- key_cand$bc
    dup_bc  <- duplicated(bc_cand) | duplicated(bc_cand, fromLast = TRUE)
    m2 <- match(key_loom$bc[need_fallback], bc_cand)
    valid <- which(!is.na(m2) & !dup_bc[m2])
    if (length(valid) > 0) m[need_fallback[valid]] <- m2[valid]
  }
  
  keep <- which(!is.na(m))
  message(sprintf("    mapped cells: %d / %d", length(keep), ncol(em)))
  if (length(keep) < 50) warning(sprintf("    WARNING: matched cells only %d for %s", length(keep), sid))
  
  colnames(em)[keep] <- cand[m[keep]]
  colnames(nu)[keep] <- cand[m[keep]]
  em <- em[, keep, drop = FALSE]
  nu <- nu[, keep, drop = FALSE]
  
  ## 去重 & 对齐
  em <- em[, unique(colnames(em)), drop = FALSE]
  nu <- nu[, colnames(em),         drop = FALSE]
  
  emat_list[[sid]] <- if (!inherits(em, "dgCMatrix")) as(em, "dgCMatrix") else em
  nmat_list[[sid]] <- if (!inherits(nu, "dgCMatrix")) as(nu, "dgCMatrix") else nu
  invisible(gc())
}

## ---------- 5) 基因对齐合并 ----------
cbind_sparse_by_rownames <- function(mats) {
  stopifnot(length(mats) >= 1)
  mats <- lapply(mats, function(m) if (!inherits(m, "dgCMatrix")) as(m, "dgCMatrix") else m)
  genes <- unique(unlist(lapply(mats, rownames), use.names = FALSE))
  mats2 <- lapply(mats, function(m) {
    miss <- setdiff(genes, rownames(m))
    if (length(miss) > 0) {
      add <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                  dims = c(length(miss), ncol(m)),
                                  dimnames = list(miss, colnames(m)))
      m <- rbind(m, add)
    }
    m[genes, , drop = FALSE]
  })
  Reduce(function(a, b) Matrix::cbind2(a, b), mats2)
}
emat_all <- cbind_sparse_by_rownames(emat_list)
nmat_all <- cbind_sparse_by_rownames(nmat_list)

## 仅保留当前Seurat细胞，并按 Cells(seu) 排序
cells_common <- intersect(Cells(seu), colnames(emat_all))
if (length(cells_common) <= 200) stop("交集细胞数过少：", length(cells_common))
seu      <- subset(seu, cells = cells_common)
emat_all <- emat_all[, Cells(seu), drop = FALSE]
nmat_all <- nmat_all[, Cells(seu), drop = FALSE]
emb_umap <- emb_umap[Cells(seu), , drop = FALSE]  # 重新对齐

## ---------- 6) 基因过滤（按当前Idents稳健过滤） ----------
clusters <- Idents(seu)
emat_f <- filter.genes.by.cluster.expression(emat_all, clusters, min.max.cluster.average = 0.05)
nmat_f <- filter.genes.by.cluster.expression(nmat_all,  clusters, min.max.cluster.average = 0.02)

## ---------- 7) 挂到Seurat（留痕） ----------
seu[["spliced"]]   <- CreateAssayObject(counts = emat_f)
seu[["unspliced"]] <- CreateAssayObject(counts = nmat_f)

## ---------- 8) 估计velocity ----------
message("[*] Estimating gene-relative velocity ...")
rvel <- gene.relative.velocity.estimates(
  emat = emat_f,
  nmat = nmat_f,
  deltaT = deltaT, kCells = kCells, fit.quantile = fit.quantile,
  n.cores = n.cores,
  cell.emb = emb_umap
)

## ---------- 9) 调色（JAMA风格 + 亚群色） ----------
ac <- function(col, alpha = 1) scales::alpha(col, alpha)
jama <- c(
  ink   = "#374E55",  gold  = "#DFBF44",
  cyan  = "#00A1D5",  red   = "#B24745",
  teal  = "#79AF97",  purple= "#6A5599",
  gray  = "#80796B"
)

clusters <- Idents(seu)  # 确保与你前面一致
pal_cluster <- setNames(scales::hue_pal()(length(levels(clusters))), levels(clusters))
cols_cluster <- pal_cluster[as.character(clusters)]; names(cols_cluster) <- Cells(seu)

## —— 工具1：在 velocity 基底图上标注分群名称（用簇中心点，鲁棒=中位数）——
label_clusters_on_embedding <- function(emb, clusters, cex = 0.95, col = "#2B2B2B") {
  stopifnot(nrow(emb) == length(clusters))
  cl <- factor(clusters)
  levs <- levels(cl)
  centers <- sapply(levs, function(lev){
    idx <- which(cl == lev)
    c(median(emb[idx,1]), median(emb[idx,2]))
  })
  centers <- t(centers); rownames(centers) <- levs
  text(centers[,1], centers[,2], labels = rownames(centers),
       cex = cex, col = col, font = 2)
}

## —— 工具2：画每个簇的“凸包”作为分群轮廓（极浅填充 + 稍深边界）——
draw_cluster_hulls <- function(emb, clusters, col_map, alpha_fill = 0.05, alpha_border = 0.7, lwd = 0.8) {
  cl <- factor(clusters)
  for (lev in levels(cl)) {
    pts <- emb[cl == lev, , drop = FALSE]
    if (nrow(pts) < 3) next
    h <- try(chull(pts), silent = TRUE)
    if (inherits(h, "try-error") || length(h) < 3) next
    polygon(pts[h,1], pts[h,2],
            col = ac(col_map[lev], alpha_fill),
            border = ac(col_map[lev], alpha_border),
            lwd = lwd)
  }
}

## 若有group（两组）则给定配色；否则略过
pal_group <- c("1A4_Antibody" = jama["red"], "Control" = jama["teal"])

## ---------- 10) 输出图：按“亚群/cluster”着色（弱化点，突出箭头 + 标出分群） ----------
message("[*] Plotting velocity on UMAP (by subcluster) ...")

## 可调参数：弱化点 & 强化箭头 & 更密网格
cell_alpha  <- 0.08   # 细胞点透明度（越小越淡）
cell_cex    <- 0.28   # 细胞点大小（越小越不显眼）
arrow_scale <- 5.0    # 箭头更醒目（长度缩放）
arrow_lwd   <- 1.2    # 箭头线更粗
grid_n      <- 60     # 网格更密，箭头分布更“满”
min_mass    <- 0.35   # 网格最小细胞质量阈值（略降可更均匀）

pdf(file.path(outdir, "velocity_umap_by_subcluster.pdf"), width = 8, height = 6)
velocyto.R::show.velocity.on.embedding.cor(
  emb = emb_umap,
  vel = rvel,
  n = 200, scale = "sqrt",
  cell.colors = ac(cols_cluster, alpha = cell_alpha),
  cex = cell_cex,
  ## 关键：函数内部用 pch=16，是无描边实心点；这里靠透明度+尺寸来“隐身”
  arrow.scale = arrow_scale,
  show.grid.flow = TRUE, min.grid.cell.mass = min_mass,
  grid.n = grid_n, arrow.lwd = arrow_lwd,
  do.par = FALSE,
  main = "Velocity on UMAP (colored by subcluster)"
)
## 在同一底图上标注“分群”：先画簇凸包，再打标签（让箭头保持视觉主角）
draw_cluster_hulls(emb_umap, clusters, pal_cluster, alpha_fill = 0.06, alpha_border = 0.8, lwd = 0.9)
label_clusters_on_embedding(emb_umap, clusters, cex = 0.95, col = "#2B2B2B")
dev.off()

## UMAP互证（亚群标签）——保持简洁，但显式标注簇名
pdf(file.path(outdir, "UMAP_subcluster_labeled.pdf"), width = 7.2, height = 6)
print(DimPlot(seu, reduction = "umap", group.by = idents_col, label = TRUE, repel = TRUE))
dev.off()

## ---------- 10b) （可选）按group着色/分组内独立拟合（沿用同样风格） ----------
if (has_group) {
  grp_vals <- as.character(seu@meta.data[[group_col]])
  if (!all(sort(unique(grp_vals)) %in% names(pal_group))) {
    message("[*] NOTE: group levels are: ", paste(sort(unique(grp_vals)), collapse = ", "))
    levs <- sort(unique(grp_vals))
    pal_group <- setNames(scales::hue_pal()(length(levs)), levs)
  }
  cols_group <- setNames(unname(pal_group[grp_vals]), Cells(seu))
  
  ## 按group着色的velocity整图（弱化点、凸包 + 标签仍按亚群）
  pdf(file.path(outdir, "velocity_umap_by_group.pdf"), width = 8, height = 6)
  velocyto.R::show.velocity.on.embedding.cor(
    emb = emb_umap,
    vel = rvel,
    n = 200, scale = "sqrt",
    cell.colors = ac(cols_group, alpha = cell_alpha),
    cex = cell_cex, arrow.scale = arrow_scale,
    show.grid.flow = TRUE, min.grid.cell.mass = min_mass,
    grid.n = grid_n, arrow.lwd = arrow_lwd, do.par = FALSE,
    main = "Velocity on UMAP (colored by group)"
  )
  ## 仍然叠加分群轮廓+标签（让读者在group配色下也能看到簇结构）
  draw_cluster_hulls(emb_umap, clusters, pal_cluster, alpha_fill = 0.06, alpha_border = 0.8, lwd = 0.9)
  label_clusters_on_embedding(emb_umap, clusters, cex = 0.95, col = "#2B2B2B")
  dev.off()
  
  ## UMAP互证（group）
  pdf(file.path(outdir, "UMAP_by_group.pdf"), width = 7.2, height = 6)
  print(DimPlot(seu, reduction = "umap", group.by = group_col) +
          scale_color_manual(values = pal_group))
  dev.off()
  
  ## —— 分组内“独立拟合”保持相同视觉规范 —— ##
  compute_and_plot_velocity_subset <- function(
    seu, emat_all, nmat_all, cells_subset, emb_umap,
    outfile, title,
    min_spliced = 0.05, min_unspliced = 0.02,
    kCells = 30, fit.quantile = 0.02, n.cores = 16,
    arrow.scale = 5.0, grid.n = 60,
    cell_alpha = 0.08, cell_cex = 0.28
  ) {
    stopifnot(length(cells_subset) > 100)
    emat <- emat_all[, cells_subset, drop = FALSE]
    nmat <- nmat_all[, cells_subset, drop = FALSE]
    clust_sub <- Idents(seu)[cells_subset]
    emat_f <- filter.genes.by.cluster.expression(emat, clust_sub, min.max.cluster.average = min_spliced)
    nmat_f <- filter.genes.by.cluster.expression(nmat, clust_sub, min.max.cluster.average = min_unspliced)
    emb_sub <- emb_umap[cells_subset, , drop = FALSE]
    rvel_sub <- gene.relative.velocity.estimates(
      emat = emat_f, nmat = nmat_f,
      deltaT = 1, kCells = kCells, fit.quantile = fit.quantile,
      n.cores = n.cores, cell.emb = emb_sub
    )
    ## 颜色与凸包以“全集的簇调色”为准，保证视觉一致
    cols_sub <- pal_cluster[as.character(clust_sub)]; names(cols_sub) <- cells_subset
    
    pdf(outfile, width = 8, height = 6)
    velocyto.R::show.velocity.on.embedding.cor(
      emb = emb_sub,
      vel = rvel_sub,
      n = 200, scale = "sqrt",
      cell.colors = ac(cols_sub, alpha = cell_alpha),
      cex = cell_cex, arrow.scale = arrow.scale,
      show.grid.flow = TRUE, min.grid.cell.mass = min_mass,
      grid.n = grid.n, arrow.lwd = arrow_lwd, do.par = FALSE,
      main = title
    )
    draw_cluster_hulls(emb_sub, clust_sub, pal_cluster, alpha_fill = 0.06, alpha_border = 0.8, lwd = 0.9)
    label_clusters_on_embedding(emb_sub, clust_sub, cex = 0.95, col = "#2B2B2B")
    dev.off()
    invisible(rvel_sub)
  }
  
  grp_vec <- seu@meta.data[[group_col]]
  levs_g  <- sort(unique(grp_vals))
  cells_A <- rownames(seu@meta.data)[grp_vec == levs_g[1]]
  cells_B <- rownames(seu@meta.data)[grp_vec == levs_g[2]]
  
  if (length(cells_A) >= 100) {
    rvel_A <- compute_and_plot_velocity_subset(
      seu, emat_all, nmat_all, cells_A, emb_umap,
      outfile = file.path(outdir, glue("velocity_umap_group_{levs_g[1]}_CALC.pdf")),
      title   = glue("Velocity (独立拟合) — {levs_g[1]}")
    )
    saveRDS(rvel_A, file.path(outdir, glue("velocity_rvel_group_{levs_g[1]}.rds")))
  } else {
    warning(glue("组 {levs_g[1]} 细胞数过少（{length(cells_A)}），跳过独立拟合。"))
  }
  
  if (length(cells_B) >= 100) {
    rvel_B <- compute_and_plot_velocity_subset(
      seu, emat_all, nmat_all, cells_B, emb_umap,
      outfile = file.path(outdir, glue("velocity_umap_group_{levs_g[2]}_CALC.pdf")),
      title   = glue("Velocity (独立拟合) — {levs_g[2]}")
    )
    saveRDS(rvel_B, file.path(outdir, glue("velocity_rvel_group_{levs_g[2]}.rds")))
  } else {
    warning(glue("组 {levs_g[2]} 细胞数过少（{length(cells_B)}），跳过独立拟合。"))
  }
}
## ========================= 仅绘图：两个 group 的 velocity 可视化（顺序配色版） =========================
## 依赖：seu, emb_umap, outdir, group_col, idents_col 已在上文定义
suppressPackageStartupMessages({
  library(glue); library(scales)
})

## ---------- 0) 固定的颜色面板（按顺序循环分配；不与名称匹配） ----------
celltype_l4_levels <- c("Tumor", "Macrophage","Monocyte",
                        "TA Cell", "T cell","B cell",
                        "Neutrophil", "Epithelial","NK", "Fibroblast","Endothelial", "Mast",
                        "DC")
celltype_l4_colors <- c("#E64B35FF", "#4DBBD5FF","#2F70AF", "#F1BAA1", "#F39B7FFF",
                        "#3C5488FF", "#00A087FF", "#8491B4FF", "#91D1C2FF", "#AA1803",
                        "#7E6148FF", "#B09C85FF", "black")

## ---------- 1) 辅助函数（标签/凸包按给定 col_map 着色；不做名称匹配） ----------
if (!exists("label_clusters_on_embedding")) {
  label_clusters_on_embedding <- function(emb, clusters, cex = 0.95,
                                          col_map = NULL, default_col = "#2B2B2B") {
    stopifnot(nrow(emb) == length(clusters))
    cl   <- factor(clusters)
    levs <- levels(cl)
    centers <- sapply(levs, function(lev){
      idx <- which(cl == lev)
      c(median(emb[idx,1]), median(emb[idx,2]))
    })
    centers <- t(centers); rownames(centers) <- levs
    
    col_vec <- if (!is.null(col_map)) {
      unname(col_map[levs])
    } else {
      rep(default_col, length(levs))
    }
    text(centers[,1], centers[,2], labels = rownames(centers),
         cex = cex, col = col_vec, font = 2)
  }
}

if (!exists("draw_cluster_hulls")) {
  draw_cluster_hulls <- function(emb, clusters, col_map,
                                 alpha_fill = 0.06, alpha_border = 0.8, lwd = 0.9,
                                 default_col = "#2B2B2B") {
    cl <- factor(clusters)
    for (lev in levels(cl)) {
      pts <- emb[cl == lev, , drop = FALSE]
      if (nrow(pts) < 3) next
      h <- try(chull(pts), silent = TRUE)
      if (inherits(h, "try-error") || length(h) < 3) next
      col_this <- if (!is.null(col_map) && !is.na(col_map[lev])) col_map[lev] else default_col
      polygon(pts[h,1], pts[h,2],
              col    = scales::alpha(col_this, alpha_fill),
              border = scales::alpha(col_this, alpha_border),
              lwd    = lwd)
    }
  }
}

## ---------- 2) 选择“亚群注释”来源 ----------
stopifnot(group_col %in% colnames(seu@meta.data))
if (exists("idents_col") && is.character(idents_col) && idents_col %in% colnames(seu@meta.data)) {
  cluster_labels <- as.character(seu@meta.data[[idents_col]])   # 你的 Macrophage 亚群注释列
} else {
  cluster_labels <- as.character(Idents(seu))                    # 或者用 active.ident
}
names(cluster_labels) <- Cells(seu)

## ---------- 3) 基于“亚群的 factor 顺序”构建顺序配色映射 ----------
## 说明：不和名称匹配，直接按 levels 顺序把 celltype_l4_colors 轮流分给各亚群
levs_all <- levels(factor(cluster_labels))  # 保留原有 factor 顺序（建议你在上游就把 Idents 的 levels 排好）
pal_for_clusters <- setNames(
  rep(celltype_l4_colors, length.out = length(levs_all)),
  levs_all
)

## 是否让“文字标签颜色”使用 panel（TRUE）或统一色（FALSE）
label_use_palette <- TRUE
label_default_col <- "#2B2B2B"
## 文字描边（无依赖版）：先画一圈黑边，再画彩色正字
text_stroke <- function(x, y, labels,
                        col,                  # 正文字颜色
                        cex = 0.95, font = 2,
                        stroke_col = "black", # 描边颜色
                        stroke = 0.9,         # 描边粗细（相对，调这一个就行）
                        ...) {
  # 根据字符尺寸自适应偏移，确保不同设备上粗细一致
  dx <- strwidth("W",  cex = cex) * stroke / 25
  dy <- strheight("W", cex = cex) * stroke / 25
  offs <- rbind(
    c(-1, 0), c(1, 0), c(0, -1), c(0, 1),
    c(-1,-1), c(-1, 1), c(1,-1), c(1, 1)
  )
  # 允许越界描边（比如刚好在边上）
  old_xpd <- par("xpd"); on.exit(par(xpd = old_xpd), add = TRUE); par(xpd = NA)
  apply(offs, 1, function(o)
    text(x + o[1]*dx, y + o[2]*dy, labels, cex = cex, font = font,
         col = stroke_col, ...))
  # 正字覆盖
  text(x, y, labels, cex = cex, font = font, col = col, ...)
}

## 每个细胞的点颜色（弱化显示，但仍与亚群对应）
cols_per_cell <- setNames(unname(pal_for_clusters[cluster_labels]), Cells(seu))

## ---------- 4) 可调视觉参数 ----------
cell_alpha  <- 0.3   # 点透明度（0~1）
cell_cex    <- 0.25    # 点大小
arrow_scale <- 3.0    # 箭头长度缩放
arrow_lwd   <- 1    # 箭头线宽
grid_n      <- 30     # 箭头网格密度
min_mass    <- 0.35   # 网格最小细胞质量阈值

## —— 工具：透明度包装（避免依赖 scales）——
alpha_col <- function(col, alpha = 1) grDevices::adjustcolor(col, alpha.f = alpha)

## —— 工具：对过浅的颜色自动加深（提高清晰度）——
safe_darken_col <- function(cols, factor = 0.8, min_lum = 0.65) {
  cols <- as.character(cols)
  vapply(cols, function(cl) {
    rgbv <- grDevices::col2rgb(cl) / 255
    lum  <- 0.299 * rgbv[1] + 0.587 * rgbv[2] + 0.114 * rgbv[3]
    if (lum > min_lum) {
      grDevices::adjustcolor(cl, red.f = factor, green.f = factor, blue.f = factor, alpha.f = 1)
    } else {
      grDevices::adjustcolor(cl, alpha.f = 1)
    }
  }, FUN.VALUE = character(1))
}

## —— 工具：简单“排斥”算法，避免标签重叠（基于文本包围盒）——
repel_labels <- function(x, y, labels, cex = 0.95, max_iter = 200, step = 0.02, margin = 0.01) {
  stopifnot(length(x) == length(y), length(labels) == length(x))
  xy <- cbind(x, y)
  wx <- strwidth(labels,  cex = cex)
  wy <- strheight(labels, cex = cex)
  
  usr <- par("usr"); dx <- usr[2] - usr[1]; dy <- usr[4] - usr[3]
  for (it in seq_len(max_iter)) {
    moved <- FALSE
    for (i in seq_len(nrow(xy) - 1)) {
      for (j in (i + 1):nrow(xy)) {
        if (abs(xy[i,1] - xy[j,1]) < (wx[i] + wx[j]) / 2 &&
            abs(xy[i,2] - xy[j,2]) < (wy[i] + wy[j]) / 2) {
          v <- xy[i, ] - xy[j, ]
          if (!any(v != 0)) v <- c(runif(1, -1, 1), runif(1, -1, 1))
          v <- v / sqrt(sum(v^2))
          xy[i, ] <- xy[i, ] + step * v
          xy[j, ] <- xy[j, ] - step * v
          moved <- TRUE
        }
      }
    }
    if (!moved) break
    ## 保持在绘图区内，留一点边距
    xy[,1] <- pmin(pmax(xy[,1], usr[1] + margin * dx), usr[2] - margin * dx)
    xy[,2] <- pmin(pmax(xy[,2], usr[3] + margin * dy), usr[4] - margin * dy)
  }
  xy
}

## ---------- 5) 遍历 group 绘图（修复轴范围/裁切/空level/箭头位置） ----------
grp_vals <- as.character(seu@meta.data[[group_col]])
levs_g   <- sort(unique(grp_vals))

for (g in levs_g) {
  cells_g <- rownames(seu@meta.data)[grp_vals == g]
  if (length(cells_g) < 50) {
    warning(glue("分组 {g} 细胞数过少（{length(cells_g)}），跳过绘图。"))
    next
  }
  emb_g <- emb_umap[cells_g, , drop = FALSE]
  
  ## 统一坐标范围：extendrange + 适度留白（避免裁切）
  xr  <- range(emb_g[,1]); yr <- range(emb_g[,2])
  pad <- 0.1
  xlim <- extendrange(xr, f = pad)
  ylim <- extendrange(yr, f = pad)
  
  ## 设备尺寸根据纵横比自适应
  asp <- diff(ylim) / diff(xlim)
  W <- 7
  H <- max(6, W * asp)
  
  ## 读入/选择 rvel
  safe_g  <- gsub("[^A-Za-z0-9._-]", "_", g)
  rvel_fp <- file.path(outdir, glue("velocity_rvel_group_{safe_g}.rds"))
  if (!file.exists(rvel_fp)) {
    if (!exists("rvel")) {
      warning(glue("对象 rvel 不存在，且未找到 {basename(rvel_fp)}，跳过 {g}。"))
      next
    }
    rvel_use <- rvel
  } else {
    rvel_use <- readRDS(rvel_fp)
  }
  
  ## 该组内实际出现的亚群 + 颜色映射（去掉空 level 防止 NA 质心）
  cl_g0      <- factor(cluster_labels[cells_g], levels = levs_all)
  cl_g       <- droplevels(cl_g0)
  levs_here  <- levels(cl_g)
  pal_g_map  <- pal_for_clusters[levs_here]
  
  ## 出图：先开空画布（固定坐标系），再叠 velocity/点/凸包/标签
  pdf(file.path(outdir, glue("velocity_umap_group_{safe_g}_PLOT_only.pdf")),
      width = W, height = H, pointsize = 12)
  
  ## 固定坐标，不要额外扩展；留出适中的边距
  par(xaxt = "n", yaxt = "n", bty = "n",
      xaxs = "i", yaxs = "i",
      mar  = c(4.5, 4.5, 2.0, 2.0))
  
  ## 空面板，锁定坐标范围
  plot(NA, NA, type = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE)
  
  ## 背景 velocity（与面板使用完全一致的 xlim/ylim）
  velocyto.R::show.velocity.on.embedding.cor(
    emb = emb_g,
    vel = rvel_use,
    n = 200,
    scale = "sqrt",
    cell.colors = alpha_col(pal_for_clusters[as.character(cl_g)], alpha = 0.0), # 我们自己叠点
    cex = cell_cex,
    arrow.scale = arrow_scale,
    show.grid.flow = TRUE,
    min.grid.cell.mass = min_mass,
    grid.n = grid_n,
    arrow.lwd = arrow_lwd,
    do.par = FALSE,
    xlim = xlim, ylim = ylim,
    main = glue("Velocity on UMAP — {g}（Macrophage 亚群，顺序配色）")
  )
  
  ## 叠点（圈+填充同色系）
  cols_core   <- pal_for_clusters[as.character(cl_g)]
  cols_border <- alpha_col(safe_darken_col(cols_core, factor = 0.8, min_lum = 0.70), alpha = 0.65)
  cols_fill   <- alpha_col(cols_core, alpha = max(0.22, cell_alpha * 1.5))
  
  points(emb_g[,1], emb_g[,2],
         pch = 21, cex = cell_cex,
         col = cols_border, bg = cols_fill, lwd = 0.2)
  
  ## 凸包（只遍历真实出现的 level）
  for (lev in levs_here) {
    idx <- which(cl_g == lev)
    if (length(idx) >= 3) {
      pts <- emb_g[idx, , drop = FALSE]
      ch  <- grDevices::chull(pts)
      polygon(pts[ch, , drop = FALSE],
              border = alpha_col(pal_g_map[lev], 0.95),
              lwd = 1.0, lty = 2,
              col = alpha_col(pal_g_map[lev], 0.05))
    }
  }
  
  ## 标签（质心 -> 排斥 -> 描边）
  centers <- t(sapply(levs_here, function(lev) colMeans(emb_g[cl_g == lev, , drop = FALSE])))
  rownames(centers) <- levs_here
  
  label_cols <- alpha_col(safe_darken_col(unname(pal_g_map), factor = 0.55, min_lum = 0.8), alpha = 1)
  names(label_cols) <- levs_here
  
  usr <- par("usr"); dx <- usr[2] - usr[1]; dy <- usr[4] - usr[3]
  lab_pos <- repel_labels(centers[,1], centers[,2],
                          labels = rownames(centers),
                          cex = 0.95, max_iter = 300,
                          step = 0.02 * min(dx, dy), margin = 0.02)
  
  text_stroke(x = lab_pos[,1], y = lab_pos[,2],
              labels = rownames(centers),
              cex = 0.95, font = 2,
              col = label_cols[rownames(centers)],
              stroke_col = "black", stroke = 1.0)
  
  ## 轴向小箭头：移到绘图区“内侧”3%，避免被 PDF 裁掉
  x0 <- xlim[1] + 0.03 * dx
  y0 <- ylim[1] + 0.03 * dy
  arrows(x0, y0, x0 + 0.09 * dx, y0,        length = 0.08, lwd = 1.1)
  arrows(x0, y0, x0, y0 + 0.09 * dy,        length = 0.08, lwd = 1.1)
  text(x0 + 0.10 * dx, y0,        "UMAP1", adj = c(0, 0.5), cex = 0.8, col = "#666666")
  text(x0,        y0 + 0.10 * dy, "UMAP2", adj = c(0.5, 0), cex = 0.8, col = "#666666")
  
  dev.off()
}
