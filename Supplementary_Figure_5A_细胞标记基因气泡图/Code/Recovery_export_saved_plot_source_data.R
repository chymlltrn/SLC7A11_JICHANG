#!/usr/bin/env Rscript
# Export only saved plotting coordinates/annotations and plot-level summaries.
# No clustering, normalization, cell typing, polarization scoring, or raw/count matrix export.
suppressPackageStartupMessages(library(Seurat))
outdir <- "/tmp/slc7a11_export_tables_20261020"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)
options(warn = 1)
log <- file("export_audit.log", "wt")
sink(log, split = TRUE); sink(log, type = "message")
on.exit <- function() {}
write_table <- function(x, filename) {
  write.csv(x, file.path(outdir, filename), row.names = FALSE, na = "")
  cat("WROTE", filename, "rows", nrow(x), "cols", ncol(x), "\n")
}
source_rows <- list()
write_object_audit <- function(obj, path, prefix) {
  fi <- file.info(path)
  source_rows[[prefix]] <<- data.frame(object_id = prefix, source_path = path,
    size_bytes = fi$size, mtime = format(fi$mtime, "%Y-%m-%d %H:%M:%S %Z"),
    object_class = paste(class(obj), collapse = ";"), n_features = nrow(obj), n_cells = ncol(obj),
    default_assay = DefaultAssay(obj), assays = paste(Assays(obj), collapse = ";"),
    reductions = paste(Reductions(obj), collapse = ";"),
    export_status = "derived_from_saved_object", stringsAsFactors = FALSE)
  md <- obj[[]]
  cols <- lapply(names(md), function(nm) {
    x <- md[[nm]]
    uv <- unique(as.character(x))
    data.frame(field = nm, field_class = paste(class(x), collapse = ";"),
      n_unique_including_na = length(uv), n_na = sum(is.na(x)),
      values_if_le_100 = if (length(uv) <= 100) paste(uv, collapse = " | ") else "[high_cardinality]",
      stringsAsFactors = FALSE)
  })
  write_table(do.call(rbind, cols), paste0(prefix, "_metadata_fields.csv"))
  cat("OBJECT", prefix, "class", class(obj), "dimensions", dim(obj), "\n")
  cat("METADATA", names(md), "\n")
  cat("ASSAYS", Assays(obj), "REDUCTIONS", Reductions(obj), "\n")
}
export_embedding <- function(obj, prefix, columns = NULL) {
  md <- obj[[]]
  if (!is.null(columns)) md <- md[, intersect(columns, names(md)), drop = FALSE]
  md <- data.frame(cell_id = rownames(md), md, check.names = FALSE)
  if ("umap" %in% Reductions(obj)) {
    emb <- Embeddings(obj, reduction = "umap")
    stopifnot(identical(md$cell_id, rownames(emb)))
    e <- as.data.frame(emb)
    names(e) <- paste0("saved_", names(e))
    md <- cbind(md, e)
  } else {
    cat("WARNING: no saved umap reduction in", prefix, "\n")
  }
  md$saved_active_ident <- as.character(Idents(obj))[match(md$cell_id, names(Idents(obj)))]
  write_table(md, paste0(prefix, "_saved_umap_metadata.csv"))
  invisible(md)
}
sample_counts <- function(samples, classes, sample_levels = unique(samples), class_levels = unique(classes)) {
  x <- as.data.frame(table(sample = factor(samples, levels = sample_levels),
                           celltype = factor(classes, levels = class_levels)), stringsAsFactors = FALSE)
  names(x)[3] <- "cell_count"
  den <- as.data.frame(table(sample = factor(samples, levels = sample_levels)), stringsAsFactors = FALSE)
  names(den)[2] <- "denominator_cells"
  x$denominator_cells <- den$denominator_cells[match(x$sample, den$sample)]
  x$proportion_fraction <- x$cell_count / x$denominator_cells
  x$proportion_percent <- 100 * x$proportion_fraction
  x
}
main_path <- "/data/public_database/Home_data/scRNAseq/1A4_anti/originaldata/annotated_scRNAsub_harmony.rds"
cat("READ", main_path, "\n")
main <- readRDS(main_path)
write_object_audit(main, main_path, "main")
export_embedding(main, "Figure4A", c("orig.ident", "sample", "sample_id", "ALLTYPE", "seurat_clusters", "group"))
markers <- c(
  "Mki67","Top2a","Birc5", "Col1a2","Col3a1","Dcn", "Ereg","Krt18","Btc",
  "Cdh5","Pecam1","Esm1", "Ctss","C1qa","C1qb", "Cd68","Adgre4","Treml4",
  "Cd3d","Cd3e","Cd3g", "Nkg7","Gzma","Gzmb", "Csf3r","S100a9","S100a8",
  "Clec9a","Xcr1","Cd24a", "Cd79a","Cd79b","Ms4a1", "Cpa3","Mcpt4","Tpsb2",
  "Dclk1","Col1a1","Col5a1")
write_table(data.frame(marker_order = seq_along(markers), gene = markers,
  present_in_default_assay = markers %in% rownames(main)), "SupplementaryFigure5A_original_marker_list.csv")
# Exact plotting call from originaldata/1030code.R, using saved normalized assay data.
# This recalculates only DotPlot layer statistics. No NormalizeData or upstream analysis is run.
p <- DotPlot(object = main, features = markers, group.by = "seurat_clusters",
             cols = c("lightgrey", "red"), dot.scale = 8)
write_table(p$data, "SupplementaryFigure5A_original_DotPlot_layer_data.csv")
md <- main[[]]
stopifnot(all(c("orig.ident", "ALLTYPE") %in% names(md)), !anyNA(md$orig.ident), !anyNA(md$ALLTYPE))
main_counts <- sample_counts(as.character(md$orig.ident), as.character(md$ALLTYPE),
  sample_levels = sort(unique(as.character(md$orig.ident))),
  class_levels = if (is.factor(md$ALLTYPE)) levels(md$ALLTYPE) else unique(as.character(md$ALLTYPE)))
main_counts$group_from_sample_id <- ifelse(grepl("^C[1-3]$", main_counts$sample), "Control",
  ifelse(grepl("^T[1-3]$", main_counts$sample), "1A4_Antibody", NA_character_))
write_table(main_counts, "SupplementaryFigure5B_all_six_samples_counts_proportions.csv")
write_table(main_counts[main_counts$sample %in% c("C1", "T2"), ],
  "SupplementaryFigure5B_manuscript_selected_C1_T2_counts_proportions.csv")
original_prop_path <- "/data/public_database/Home_data/scRNAseq/1A4_anti/originaldata/cell_proportion_data.csv"
orig <- read.csv(original_prop_path, stringsAsFactors = FALSE, check.names = FALSE)
write_table(orig, "SupplementaryFigure5B_original_cell_proportion_data.csv")
derived_key <- paste(main_counts$sample, main_counts$celltype, sep = "::")
original_key <- paste(orig$orig.ident, orig$ALLTYPE, sep = "::")
idx <- match(original_key, derived_key)
comparison <- data.frame(sample = orig$orig.ident, celltype = orig$ALLTYPE,
  original_cell_count = orig$count, derived_cell_count = main_counts$cell_count[idx],
  original_proportion_percent = orig$proportion,
  derived_proportion_percent = main_counts$proportion_percent[idx],
  count_delta = main_counts$cell_count[idx] - orig$count,
  percent_delta = main_counts$proportion_percent[idx] - orig$proportion,
  matched = !is.na(idx))
write_table(comparison, "SupplementaryFigure5B_original_vs_saved_object_validation.csv")
cat("PROPORTION_VALIDATION all rows matched", all(comparison$matched),
  "all counts match", all(comparison$count_delta == 0),
  "max abs pct difference", max(abs(comparison$percent_delta)), "\n")
stopifnot(all(comparison$matched), all(comparison$count_delta == 0),
  max(abs(comparison$percent_delta)) < 1e-8)
rm(main, md, p); gc()
for (objname in c("T_cells", "Macrophages")) {
  path <- paste0("/data/public_database/Home_data/scRNAseq/1A4_anti/1A4_scRNA/", objname, "_final_annotated.rds")
  cat("READ", path, "\n")
  obj <- readRDS(path)
  write_object_audit(obj, path, objname)
  export_embedding(obj, objname)
  md <- obj[[]]
  sample_candidates <- c("old.ident", "orig.ident", "sample", "sample_id")
  has_sample <- sample_candidates[sample_candidates %in% names(md)]
  if (length(has_sample)) {
    valid_samples <- vapply(has_sample, function(nm) {
      vals <- unique(as.character(md[[nm]])); length(vals) == 6 && setequal(vals, c("C1","C2","C3","T1","T2","T3"))
    }, logical(1))
    sample_field <- if (any(valid_samples)) has_sample[which(valid_samples)[1]] else NA_character_
  } else sample_field <- NA_character_
  class_candidates <- c("celltype", "cell_type", "CellType", "subcelltype", "T_cell_type", "T_cell_subtype", "ALLTYPE")
  has_class <- class_candidates[class_candidates %in% names(md)]
  if (length(has_class) && !is.na(sample_field)) {
    class_field <- has_class[1]
    samples <- as.character(md[[sample_field]])
    labels <- as.character(md[[class_field]])
    x <- sample_counts(samples, labels, sample_levels = c("C1","C2","C3","T1","T2","T3"))
    x$group <- ifelse(grepl("^C", x$sample), "Control", "1A4_Antibody")
    x$source_sample_field <- sample_field; x$source_celltype_field <- class_field
    x$provenance <- "derived_from_saved_object_no_original_plot_script"
    write_table(x, paste0(objname, "_per_sample_celltype_counts_proportions.csv"))
    gs <- aggregate(x[, c("proportion_fraction", "proportion_percent")],
      by = list(group = x$group, celltype = x$celltype), FUN = mean)
    names(gs)[3:4] <- c("mean_sample_proportion_fraction", "mean_sample_proportion_percent")
    gs$n_samples <- 3L; gs$sample_denominator <- paste0(objname, " within each sample")
    gs$provenance <- "derived_from_saved_object_no_original_plot_script"
    write_table(gs, paste0(objname, "_group_mean_of_sample_proportions.csv"))
    cat("GROUP_TABLE", objname, "sample_field", sample_field, "class_field", class_field, "\n")
  } else cat("GROUP_TABLE_NOT_EXPORTED", objname, "sample_field", sample_field, "class_fields", has_class, "\n")
  polarization_fields <- names(md)[vapply(md, function(v) {
    vals <- unique(as.character(v));
    any(grepl("(^|[^A-Za-z])M[12]([^A-Za-z]|$)", vals)) && length(vals) <= 100
  }, logical(1))]
  write_table(data.frame(field_with_existing_M1_M2_label = polarization_fields),
              paste0(objname, "_existing_M1_M2_label_fields.csv"))
  if (length(polarization_fields) && !is.na(sample_field)) for (pf in polarization_fields) {
    x <- sample_counts(as.character(md[[sample_field]]), as.character(md[[pf]]),
      sample_levels = c("C1","C2","C3","T1","T2","T3"))
    x$source_label_field <- pf
    x$provenance <- "existing_annotation_no_new_polarization_scoring"
    write_table(x, paste0(objname, "_existing_", gsub("[^A-Za-z0-9]", "_", pf), "_counts_proportions.csv"))
  }
  rm(obj, md); gc()
}
write_table(do.call(rbind, source_rows), "source_objects_provenance.csv")
writeLines(capture.output(sessionInfo()), "R_sessionInfo.txt")
cat("EXPORT_COMPLETE\n")
sink(type = "message"); sink(); close(log)
