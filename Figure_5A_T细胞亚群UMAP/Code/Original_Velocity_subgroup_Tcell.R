## =========================================================
## RNA velocity（仅基于已有 T 细胞 Seurat 对象 + 现有 loom）
## =========================================================
options(stringsAsFactors = FALSE)
set.seed(2025)

suppressPackageStartupMessages({
  library(Seurat); library(velocyto.R); library(Matrix)
  library(RColorBrewer); library(scales); library(dplyr); library(glue); library(ggplot2)
})

## ---------- 0) 路径 & 参数 ----------
loom_dir   <- "/data/public_database/Home_data/scRNAseq/1A4_anti/1A4_scRNA/loom_unzip"
seurat_rds <- "/data/public_database/Home_data/scRNAseq/1A4_anti/1A4_scRNA/T_cells_final_annotated.rds"
outdir     <- "/data/public_database/Home_data/scRNAseq/1A4_anti/1A4_scRNA/velocity_Tcells_251104"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

idents_col <- "celltype"  # 你的亚群列
group_col  <- "group"

deltaT <- 1; kCells <- 30; fit.quantile <- 0.02; n.cores <- 16
write_back_loom <- FALSE

## ---------- 1) 收集 loom ----------
loom_files <- list.files(loom_dir, pattern = "\\.loom$", full.names = TRUE, recursive = TRUE)
stopifnot(length(loom_files) > 0)

## ---------- 2) 读 Seurat 对象（已是 T 细胞子集 & 已注释） ----------
message("[*] Loading Seurat object (T cells only) ...")
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


## ---------- 2a) 设定 Idents（自动侦测 + 占位回退） ----------
subcluster_candidates <- c(
  idents_col, "T_cell_subtype", "Tcell_subtype", "T_subtype", "T_subcluster",
  "T_state", "T_cluster", "T.cluster", "t_subcluster",
  "annot", "annotation", "Cluster", "clusters",
  "celltype_sub", "celltype_detail",
  "seurat_clusters"
)
valid_cols <- intersect(subcluster_candidates, colnames(seu@meta.data))
if (length(valid_cols) > 0) {
  idents_col <- valid_cols[1]
  message(glue("[*] Using '{idents_col}' as Idents."))
  Idents(seu) <- idents_col
} else {
  message("[*] No suitable subcluster column found; keep current Idents(seu).")
}

## ---------- 2a-b) 仅做一处重命名：γδ_T -> gammadelta_T（其余不改） ----------
if (idents_col %in% colnames(seu@meta.data)) {
  vec <- as.character(seu@meta.data[[idents_col]])
  if (any(vec == "γδ_T", na.rm = TRUE)) {
    vec[vec == "γδ_T"] <- "gammadelta_T"   # 修改
    seu@meta.data[[idents_col]] <- vec
    Idents(seu) <- factor(vec)
  }
} else {
  ident_vec <- as.character(Idents(seu))
  if (any(ident_vec == "γδ_T", na.rm = TRUE)) {
    ident_vec[ident_vec == "γδ_T"] <- "gammadelta_T"   # 修改
    Idents(seu) <- factor(ident_vec)
  }
}

## ---------- 2b) 确保有 UMAP ----------
has_umap <- "umap" %in% names(seu@reductions)
if (!has_umap) {
  message("[*] UMAP not found in object; computing a temporary UMAP...")
  if (is.null(DefaultAssay(seu))) DefaultAssay(seu) <- "RNA"
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
  ifelse(has_sep, sub("[:_].*$", "", v), NA_character_)
}
extract_barcode <- function(v) {
  core <- sub("^.*[:_]", "", v)
  core <- sub("x$", "", core)
  core <- sub("-1$", "", core)
  core <- sub("\\.1$", "", core)
  core
}
make_keys <- function(vec) data.frame(sample = extract_sample(vec), bc = extract_barcode(vec), stringsAsFactors = FALSE)

## ---------- 4) 读入各 loom 并与 Seurat 映射 ----------
emat_list <- list(); nmat_list <- list()
for (f in loom_files) {
  sid <- tools::file_path_sans_ext(basename(f))
  message(sprintf("[*] Reading loom: %s", sid))
  ld <- velocyto.R::read.loom.matrices(f)
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

## 仅保留当前 Seurat 细胞，并按 Cells(seu) 排序
cells_common <- intersect(Cells(seu), colnames(emat_all))
if (length(cells_common) <= 200) stop("交集细胞数过少：", length(cells_common))
seu      <- subset(seu, cells = cells_common)
emat_all <- emat_all[, Cells(seu), drop = FALSE]
nmat_all <- nmat_all[, Cells(seu), drop = FALSE]
emb_umap <- emb_umap[Cells(seu), , drop = FALSE]

## ---------- 6) 基因过滤 ----------
clusters <- Idents(seu)
emat_f <- filter.genes.by.cluster.expression(emat_all, clusters, min.max.cluster.average = 0.05)
nmat_f <- filter.genes.by.cluster.expression(nmat_all,  clusters, min.max.cluster.average = 0.02)

## ---------- 7) 挂到 Seurat ----------
seu[["spliced"]]   <- CreateAssayObject(counts = emat_f)
seu[["unspliced"]] <- CreateAssayObject(counts = nmat_f)

## ---------- 8) 估计 velocity ----------
message("[*] Estimating gene-relative velocity ...")
rvel <- gene.relative.velocity.estimates(
  emat = emat_f,
  nmat = nmat_f,
  deltaT = deltaT, kCells = kCells, fit.quantile = fit.quantile,
  n.cores = n.cores,
  cell.emb = emb_umap
)

## ---------- 9) 调色 ----------
ac <- function(col, alpha = 1) scales::alpha(col, alpha)
jama <- c(ink="#374E55", gold="#DFBF44", cyan="#00A1D5", red="#B24745", teal="#79AF97", purple="#6A5599", gray="#80796B")

clusters <- Idents(seu)
pal_cluster <- setNames(scales::hue_pal()(length(levels(clusters))), levels(clusters))
cols_cluster <- pal_cluster[as.character(clusters)]; names(cols_cluster) <- Cells(seu)

label_clusters_on_embedding <- function(emb, clusters, cex = 0.95, col = "#2B2B2B") {
  stopifnot(nrow(emb) == length(clusters))
  cl <- factor(clusters)
  levs <- levels(cl)
  centers <- sapply(levs, function(lev){
    idx <- which(cl == lev)
    c(median(emb[idx,1]), median(emb[idx,2]))
  })
  centers <- t(centers); rownames(centers) <- levs
  text(centers[,1], centers[,2], labels = rownames(centers), cex = cex, col = col, font = 2)
}
draw_cluster_hulls <- function(emb, clusters, col_map, alpha_fill = 0.05, alpha_border = 0.3, lwd = 0.3) {
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

pal_group <- c("1A4_Antibody" = jama["red"], "Control" = jama["teal"])

## ---------- 10) 输出图：按“亚群/cluster”着色 ----------
message("[*] Plotting velocity on UMAP (by subcluster) ...")

## —— 这里仅微调三处参数以减少空白区箭头 ——  #
cell_alpha  <- 0.08
cell_cex    <- 0.28
arrow_scale <- 3.0   # 修改：缩短箭头
arrow_lwd   <- 1.2
grid_n      <- 30    # 修改：降低网格密度
min_mass    <- 1.0   # 修改：提高网格细胞质量阈值，抑制空白区箭头

pdf(file.path(outdir, "velocity_umap_by_subcluster.pdf"), width = 8, height = 6)
velocyto.R::show.velocity.on.embedding.cor(
  emb = emb_umap,
  vel = rvel,
  n = 200, scale = "sqrt",
  cell.colors = ac(cols_cluster, alpha = cell_alpha),
  cex = cell_cex,
  arrow.scale = arrow_scale,
  show.grid.flow = TRUE, min.grid.cell.mass = min_mass,
  grid.n = grid_n, arrow.lwd = arrow_lwd,
  do.par = FALSE,
  main = "Velocity on UMAP (colored by T-cell subclusters)"
)
draw_cluster_hulls(emb_umap, clusters, pal_cluster, alpha_fill = 0.06, alpha_border = 0.8, lwd = 0.9)
label_clusters_on_embedding(emb_umap, clusters, cex = 0.95, col = "#2B2B2B")
dev.off()

## UMAP 互证（亚群标签）
pdf(file.path(outdir, "UMAP_subcluster_labeled.pdf"), width = 7.2, height = 6)
print(DimPlot(seu, reduction = "umap", group.by = idents_col, label = TRUE, repel = TRUE))
dev.off()

## ---------- 10b) （可选）按 group 着色/分组内独立拟合 ----------
if (has_group) {
  grp_vals <- as.character(seu@meta.data[[group_col]])
  if (!all(sort(unique(grp_vals)) %in% names(pal_group))) {
    levs <- sort(unique(grp_vals))
    pal_group <- setNames(scales::hue_pal()(length(levs)), levs)
  }
  cols_group <- setNames(unname(pal_group[grp_vals]), Cells(seu))
  
  pdf(file.path(outdir, "velocity_umap_by_group.pdf"), width = 8, height = 6)
  velocyto.R::show.velocity.on.embedding.cor(
    emb = emb_umap,
    vel = rvel,
    n = 200, scale = "sqrt",
    cell.colors = ac(cols_group, alpha = cell_alpha),
    cex = cell_cex, arrow.scale = arrow_scale,
    show.grid.flow = TRUE, min.grid.cell.mass = min_mass,  # 使用同一抑制策略
    grid.n = grid_n, arrow.lwd = arrow_lwd, do.par = FALSE,
    main = "Velocity on UMAP (colored by group, T cells)"
  )
  draw_cluster_hulls(emb_umap, clusters, pal_cluster, alpha_fill = 0.06, alpha_border = 0.8, lwd = 0.9)
  label_clusters_on_embedding(emb_umap, clusters, cex = 0.95, col = "#2B2B2B")
  dev.off()
  
  pdf(file.path(outdir, "UMAP_by_group.pdf"), width = 7.2, height = 6)
  print(DimPlot(seu, reduction = "umap", group.by = group_col) +
          scale_color_manual(values = pal_group))
  dev.off()
  
  compute_and_plot_velocity_subset <- function(
    seu, emat_all, nmat_all, cells_subset, emb_umap,
    outfile, title,
    min_spliced = 0.05, min_unspliced = 0.02,
    kCells = 30, fit.quantile = 0.02, n.cores = 16,
    arrow.scale = 3.0, grid.n = 30,   # 与上保持一致（修改）
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
    cols_sub <- pal_cluster[as.character(clust_sub)]; names(cols_sub) <- cells_subset
    
    pdf(outfile, width = 8, height = 6)
    velocyto.R::show.velocity.on.embedding.cor(
      emb = emb_sub,
      vel = rvel_sub,
      n = 200, scale = "sqrt",
      cell.colors = ac(cols_sub, alpha = cell_alpha),
      cex = cell_cex, arrow.scale = arrow.scale,
      show.grid.flow = TRUE, min.grid.cell.mass = min_mass,  # 使用同一阈值（修改）
      grid.n = grid.n, arrow.lwd = 1.2, do.par = FALSE,
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
      title   = glue("Velocity (独立拟合) — {levs_g[1]} — T cells")
    )
    saveRDS(rvel_A, file.path(outdir, glue("velocity_rvel_group_{levs_g[1]}.rds")))
  } else {
    warning(glue("组 {levs_g[1]} 细胞数过少（{length(cells_A)}），跳过独立拟合。"))
  }
  
  if (length(cells_B) >= 100) {
    rvel_B <- compute_and_plot_velocity_subset(
      seu, emat_all, nmat_all, cells_B, emb_umap,
      outfile = file.path(outdir, glue("velocity_umap_group_{levs_g[2]}_CALC.pdf")),
      title   = glue("Velocity (独立拟合) — {levs_g[2]} — T cells")
    )
    saveRDS(rvel_B, file.path(outdir, glue("velocity_rvel_group_{levs_g[2]}.rds")))
  } else {
    warning(glue("组 {levs_g[2]} 细胞数过少（{length(cells_B)}），跳过独立拟合。"))
  }
}

## ========================= 仅绘图：两个 group 的 velocity 可视化（顺序配色版，T cells） =========================
suppressPackageStartupMessages({ library(glue); library(scales) })

celltype_l4_levels <- c("Tumor","Macrophage","Monocyte","TA Cell","T cell","B cell",
                        "Neutrophil","Epithelial","NK","Fibroblast","Endothelial","Mast","DC")
celltype_l4_colors <- c("#E64B35FF","#4DBBD5FF","#2F70AF","#F1BAA1","#F39B7FFF",
                        "#3C5488FF","#00A087FF","#8491B4FF","#91D1C2FF","#AA1803",
                        "#7E6148FF","#B09C85FF","black")

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
    col_vec <- if (!is.null(col_map)) unname(col_map[levs]) else rep(default_col, length(levs))
    text(centers[,1], centers[,2], labels = rownames(centers), cex = cex, col = col_vec, font = 2)
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
              col = scales::alpha(col_this, alpha_fill),
              border = scales::alpha(col_this, alpha_border),
              lwd = lwd)
    }
  }
}

stopifnot(group_col %in% colnames(seu@meta.data))
cluster_labels <- if (idents_col %in% colnames(seu@meta.data)) {
  as.character(seu@meta.data[[idents_col]])
} else {
  as.character(Idents(seu))
}
names(cluster_labels) <- Cells(seu)

levs_all <- levels(factor(cluster_labels))
pal_for_clusters <- setNames(rep(celltype_l4_colors, length.out = length(levs_all)), levs_all)

label_default_col <- "#2B2B2B"
text_stroke <- function(x, y, labels, col, cex = 0.95, font = 2, stroke_col = "black", stroke = 0.9, ...) {
  dx <- strwidth("W",  cex = cex) * stroke / 25
  dy <- strheight("W", cex = cex) * stroke / 25
  offs <- rbind(c(-1,0), c(1,0), c(0,-1), c(0,1), c(-1,-1), c(-1,1), c(1,-1), c(1,1))
  old_xpd <- par("xpd"); on.exit(par(xpd = old_xpd), add = TRUE); par(xpd = NA)
  apply(offs, 1, function(o) text(x + o[1]*dx, y + o[2]*dy, labels, cex = cex, font = font, col = stroke_col, ...))
  text(x, y, labels, cex = cex, font = font, col = col, ...)
}

cols_per_cell <- setNames(unname(pal_for_clusters[cluster_labels]), Cells(seu))

## —— 仅绘图部分也同步采用更“克制”的箭头参数 ——  #
cell_alpha  <- 0.3
cell_cex    <- 0.25
arrow_scale <- 3.0   # 修改
arrow_lwd   <- 1
grid_n      <- 30    # 修改
min_mass    <- 1.0   # 修改

alpha_col <- function(col, alpha = 1) grDevices::adjustcolor(col, alpha.f = alpha)
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
    xy[,1] <- pmin(pmax(xy[,1], usr[1] + margin * dx), usr[2] - margin * dx)
    xy[,2] <- pmin(pmax(xy[,2], usr[3] + margin * dy), usr[4] - margin * dy)
  }
  xy
}

grp_vals <- as.character(seu@meta.data[[group_col]])
levs_g   <- sort(unique(grp_vals))

for (g in levs_g) {
  cells_g <- rownames(seu@meta.data)[grp_vals == g]
  if (length(cells_g) < 50) {
    warning(glue("分组 {g} 细胞数过少（{length(cells_g)}），跳过绘图。"))
    next
  }
  emb_g <- emb_umap[cells_g, , drop = FALSE]
  
  xr  <- range(emb_g[,1]); yr <- range(emb_g[,2])
  pad <- 0.1
  xlim <- extendrange(xr, f = pad)
  ylim <- extendrange(yr, f = pad)
  asp <- diff(ylim) / diff(xlim)
  W <- 7
  H <- max(6, W * asp)
  
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
  
  cl_g0      <- factor(cluster_labels[cells_g], levels = levs_all)
  cl_g       <- droplevels(cl_g0)
  levs_here  <- levels(cl_g)
  pal_g_map  <- pal_for_clusters[levs_here]
  
  pdf(file.path(outdir, glue("velocity_umap_group_{safe_g}_PLOT_only.pdf")),
      width = W, height = H, pointsize = 12)
  par(xaxt = "n", yaxt = "n", bty = "n", xaxs = "i", yaxs = "i", mar = c(4.5, 4.5, 2.0, 2.0))
  plot(NA, NA, type = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE)
  
  velocyto.R::show.velocity.on.embedding.cor(
    emb = emb_g,
    vel = rvel_use,
    n = 200,
    scale = "sqrt",
    cell.colors = alpha_col(pal_for_clusters[as.character(cl_g)], alpha = 0.0),
    cex = cell_cex,
    arrow.scale = arrow_scale,
    show.grid.flow = TRUE,
    min.grid.cell.mass = min_mass,   # 修改
    grid.n = grid_n,                 # 修改
    arrow.lwd = arrow_lwd,
    do.par = FALSE,
    xlim = xlim, ylim = ylim,
    main = glue("Velocity on UMAP — {g}（T-cell subclusters，顺序配色）")
  )
  
  cols_core   <- pal_for_clusters[as.character(cl_g)]
  cols_border <- alpha_col(safe_darken_col(cols_core, factor = 0.8, min_lum = 0.70), alpha = 0.65)
  cols_fill   <- alpha_col(cols_core, alpha = max(0.22, cell_alpha * 1.5))
  points(emb_g[,1], emb_g[,2], pch = 21, cex = cell_cex, col = cols_border, bg = cols_fill, lwd = 0.2)
  
  for (lev in levs_here) {
    idx <- which(cl_g == lev)
    if (length(idx) >= 3) {
      pts <- emb_g[idx, , drop = FALSE]
      ch  <- grDevices::chull(pts)
      polygon(pts[ch, , drop = FALSE],
              border = alpha_col(pal_g_map[lev], 0.95),
              lwd = 0.4, lty = 2,
              col = alpha_col(pal_g_map[lev], 0.05))
    }
  }
  
  centers <- t(sapply(levs_here, function(lev) colMeans(emb_g[cl_g == lev, , drop = FALSE])))
  rownames(centers) <- levs_here
  label_cols <- alpha_col(safe_darken_col(unname(pal_g_map), factor = 0.55, min_lum = 0.8), alpha = 1)
  names(label_cols) <- levs_here
  
  usr <- par("usr"); dx <- usr[2] - usr[1]; dy <- usr[4] - usr[3]
  lab_pos <- repel_labels(centers[,1], centers[,2], labels = rownames(centers),
                          cex = 0.95, max_iter = 300, step = 0.02 * min(dx, dy), margin = 0.02)
  text_stroke(x = lab_pos[,1], y = lab_pos[,2],
              labels = rownames(centers), cex = 0.95, font = 2,
              col = label_cols[rownames(centers)], stroke_col = "black", stroke = 0.4)
  
  x0 <- xlim[1] + 0.03 * dx
  y0 <- ylim[1] + 0.03 * dy
  arrows(x0, y0, x0 + 0.09 * dx, y0, length = 0.08, lwd = 1.1)
  arrows(x0, y0, x0, y0 + 0.09 * dy, length = 0.08, lwd = 1.1)
  text(x0 + 0.10 * dx, y0, "UMAP1", adj = c(0, 0.5), cex = 0.8, col = "#666666")
  text(x0, y0 + 0.10 * dy, "UMAP2", adj = c(0.5, 0), cex = 0.8, col = "#666666")
  
  dev.off()
}
