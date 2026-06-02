# Script 1: Bulk microarray integration + limma reproduction script
# ────────────────────────────────────────────────────────────────────────────────
# Purpose
#   Build the 4-cohort integration pipeline.
#
# Caution
#   Archived exploratory/debugging sections are retained as part of the analysis
#   record and are not required for interpreting the reported results.
#────────────────────────────────────────────────────────────────────────────────


###─────────────────────────────────────── A: Integrate 4 cohorts and detect DEGs---------------------------------------

#  0. User parameters------------------------------------------------------------

proj_root_bulk              <- "C:/DMD_project"
raw_dir_bulk                <- file.path(proj_root_bulk, "exprdata")
cache_dir_bulk              <- file.path(proj_root_bulk, "cache_bulk_revision")
res_dir_bulk                <- file.path(proj_root_bulk, "results_bulk_revision")
original_cache_dir_bulk     <- file.path(proj_root_bulk, "cache")

age_cutoff_bulk             <- 5
treat_lfc_bulk              <- 0.15

mapping_mode_bulk           <- "legacy_select"

min_measured_samples_bulk   <- 6
min_present_each_cohort_bulk<- 1
knn_rowmax_bulk             <- 1
knn_colmax_bulk             <- 0.99

combat_mean_only_bulk       <- FALSE
compare_to_original_bulk    <- TRUE
write_stage_rds_bulk        <- TRUE
write_deg_table_bulk        <- TRUE
write_qc_plots_bulk         <- TRUE
auto_install_bulk           <- FALSE

seed_bulk                   <- 123
comparison_tolerance_bulk   <- 1e-10


#  1. Directory setup------------------------------------------------------------

dir.create(cache_dir_bulk, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir_bulk,   showWarnings = FALSE, recursive = TRUE)

if ("package:conflicted" %in% search()) {
  conflicted::conflict_prefer("unname",    "base", quiet = TRUE)
  conflicted::conflict_prefer("union",     "base", quiet = TRUE)
  conflicted::conflict_prefer("intersect", "base", quiet = TRUE)
  conflicted::conflict_prefer("setdiff",   "base", quiet = TRUE)
  conflicted::conflict_prefer("Reduce",    "base", quiet = TRUE)
  conflicted::conflict_prefer("Filter",    "base", quiet = TRUE)
  conflicted::conflict_prefer("setequal",  "base", quiet = TRUE)
}


#  2. Package management---------------------------------------------------------

cran_pkgs_bulk <- c(
  "digest",
  "matrixStats"
)

bioc_pkgs_bulk <- c(
  "Biobase",
  "limma",
  "sva",
  "impute",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "hgu133plus2.db",
  "hgu133a.db",
  "hgu133b.db",
  "hgu95av2.db",
  "hgu95b.db",
  "hgu95c.db",
  "hgu95d.db",
  "hgu95e.db"
)

install_if_missing_bulk <- function(pkgs_bulk, installer_bulk, ...) {
  need_bulk <- base::setdiff(pkgs_bulk, rownames(installed.packages()))
  if (length(need_bulk) > 0) {
    installer_bulk(need_bulk, ...)
  }
}

if (auto_install_bulk) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cran.rstudio.com")
  }
  install_if_missing_bulk(
    pkgs_bulk   = cran_pkgs_bulk,
    installer_bulk = install.packages,
    repos       = "https://cran.rstudio.com",
    dependencies= TRUE
  )
  install_if_missing_bulk(
    pkgs_bulk   = bioc_pkgs_bulk,
    installer_bulk = BiocManager::install,
    ask         = FALSE,
    update      = FALSE
  )
}

missing_cran_bulk <- cran_pkgs_bulk[
  !vapply(cran_pkgs_bulk, requireNamespace, logical(1), quietly = TRUE)
]
missing_bioc_bulk <- bioc_pkgs_bulk[
  !vapply(bioc_pkgs_bulk, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran_bulk) > 0 || length(missing_bioc_bulk) > 0) {
  stop(
    "Missing packages for the revised bulk script.\n",
    "CRAN: ", paste(missing_cran_bulk, collapse = ", "), "\n",
    "Bioconductor: ", paste(missing_bioc_bulk, collapse = ", "), "\n",
    "Either install them first or set auto_install_bulk <- TRUE."
  )
}

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(sva)
  library(impute)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(hgu133plus2.db)
  library(hgu133a.db)
  library(hgu133b.db)
  library(hgu95av2.db)
  library(hgu95b.db)
  library(hgu95c.db)
  library(hgu95d.db)
  library(hgu95e.db)
  library(matrixStats)
  library(digest)
})

options(stringsAsFactors = FALSE, scipen = 999)


#  3. Helper functions-----------------------------------------------------------

assert_file_exists_bulk <- function(path_bulk) {
  if (!file.exists(path_bulk)) {
    stop("Required file not found: ", path_bulk)
  }
}

strip_gsm_bulk <- function(x_bulk) {
  sub("(GSM[0-9]+).*", "\\1", x_bulk)
}

needs_log2_bulk <- function(mat_bulk, cut_bulk = 100) {
  max(mat_bulk, na.rm = TRUE) > cut_bulk
}

match_gsm_indices_bulk <- function(expected_gsm_bulk, observed_sample_names_bulk, label_bulk) {
  idx_bulk <- match(
    strip_gsm_bulk(expected_gsm_bulk),
    strip_gsm_bulk(observed_sample_names_bulk)
  )
  if (anyNA(idx_bulk)) {
    missing_bulk <- expected_gsm_bulk[is.na(idx_bulk)]
    stop(
      "Some expected GSM IDs were not found in ", label_bulk, ": ",
      paste(missing_bulk, collapse = ", ")
    )
  }
  idx_bulk
}

read_sm_bulk <- function(fname_bulk, raw_dir_bulk) {
  path_bulk <- file.path(raw_dir_bulk, fname_bulk)
  assert_file_exists_bulk(path_bulk)
  
  con_bulk <- gzfile(path_bulk, "rt")
  txt_bulk <- suppressWarnings(readLines(con_bulk, warn = FALSE))
  close(con_bulk)
  
  beg_bulk <- grep("!series_matrix_table_begin", txt_bulk, fixed = TRUE)[1] + 1
  end_bulk <- grep("!series_matrix_table_end",   txt_bulk, fixed = TRUE)[1] - 1
  
  if (is.na(beg_bulk) || is.na(end_bulk) || end_bulk < beg_bulk) {
    stop("Failed to detect the series-matrix table in: ", fname_bulk)
  }
  
  mat_bulk <- read.delim(
    textConnection(txt_bulk[beg_bulk:end_bulk]),
    header           = TRUE,
    sep              = "\t",
    quote            = "\"",
    row.names        = 1,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  
  strip_quote_bulk <- function(x_bulk) gsub('^"|"$', "", x_bulk)
  
  rownames(mat_bulk) <- strip_quote_bulk(rownames(mat_bulk))
  colnames(mat_bulk) <- strip_quote_bulk(colnames(mat_bulk))
  
  Biobase::ExpressionSet(as.matrix(mat_bulk))
}

pad_rows_bulk <- function(mat_bulk, genes_bulk) {
  out_bulk <- matrix(
    NA_real_,
    nrow     = length(genes_bulk),
    ncol     = ncol(mat_bulk),
    dimnames = list(genes_bulk, colnames(mat_bulk))
  )
  hit_bulk <- base::intersect(rownames(mat_bulk), genes_bulk)
  out_bulk[hit_bulk, ] <- mat_bulk[hit_bulk, , drop = FALSE]
  out_bulk
}

pad_rows_cols_bulk <- function(mat_bulk, genes_bulk, cols_bulk) {
  out_bulk <- matrix(
    NA_real_,
    nrow     = length(genes_bulk),
    ncol     = length(cols_bulk),
    dimnames = list(genes_bulk, cols_bulk)
  )
  row_hit_bulk <- base::intersect(rownames(mat_bulk), genes_bulk)
  col_hit_bulk <- base::intersect(colnames(mat_bulk), cols_bulk)
  out_bulk[row_hit_bulk, col_hit_bulk] <- mat_bulk[row_hit_bulk, col_hit_bulk, drop = FALSE]
  out_bulk
}

rename_cols_by_pid_bulk <- function(mat_bulk, gsm_vec_bulk, pid_vec_bulk) {
  gsm_map_bulk <- setNames(pid_vec_bulk, strip_gsm_bulk(gsm_vec_bulk))
  cn_bulk      <- strip_gsm_bulk(colnames(mat_bulk))
  
  mapped_bulk <- gsm_map_bulk[cn_bulk]
  
  if (anyNA(mapped_bulk)) {
    stop(
      "Some columns could not be mapped to PID: ",
      paste(colnames(mat_bulk)[is.na(mapped_bulk)], collapse = ", ")
    )
  }
  
  colnames(mat_bulk) <- base::unname(mapped_bulk)
  mat_bulk
}

avg_two_bulk <- function(mat_a_bulk, mat_b_bulk) {
  genes_bulk <- base::union(rownames(mat_a_bulk), rownames(mat_b_bulk))
  ids_bulk   <- base::union(colnames(mat_a_bulk), colnames(mat_b_bulk))
  
  x_bulk <- pad_rows_cols_bulk(mat_a_bulk, genes_bulk, ids_bulk)
  y_bulk <- pad_rows_cols_bulk(mat_b_bulk, genes_bulk, ids_bulk)
  
  cnt_bulk <- (!is.na(x_bulk)) + (!is.na(y_bulk))
  sum_bulk <- replace(x_bulk, is.na(x_bulk), 0) + replace(y_bulk, is.na(y_bulk), 0)
  
  out_bulk <- sum_bulk / pmax(cnt_bulk, 1)
  out_bulk[cnt_bulk == 0] <- NA_real_
  
  out_bulk
}

avg_multi_bulk <- function(mat_list_bulk) {
  stopifnot(length(mat_list_bulk) >= 1)
  
  genes_bulk <- base::Reduce(base::union, lapply(mat_list_bulk, rownames))
  cols_bulk  <- base::Reduce(base::union, lapply(mat_list_bulk, colnames))
  
  padded_bulk <- lapply(mat_list_bulk, function(m_bulk) {
    pad_rows_cols_bulk(m_bulk, genes_bulk, cols_bulk)
  })
  
  sum_bulk <- base::Reduce(
    `+`,
    lapply(padded_bulk, function(m_bulk) replace(m_bulk, is.na(m_bulk), 0))
  )
  
  cnt_bulk <- base::Reduce(
    `+`,
    lapply(padded_bulk, function(m_bulk) !is.na(m_bulk))
  )
  
  out_bulk <- sum_bulk / pmax(cnt_bulk, 1)
  out_bulk[cnt_bulk == 0] <- NA_real_
  
  out_bulk
}

find_dup_cols_bulk <- function(mat_bulk, algo_bulk = "xxhash64") {
  h_bulk <- apply(mat_bulk, 2, digest::digest, algo = algo_bulk)
  which(duplicated(h_bulk) | duplicated(h_bulk, fromLast = TRUE))
}

map_to_entrez_bulk <- function(expr_bulk, chip_pkg_bulk, mapping_mode_bulk) {
  mapping_mode_bulk <- match.arg(
    mapping_mode_bulk,
    choices = c("legacy_select", "mapids_first", "mapids_asNA")
  )
  
  if (all(grepl("^[0-9]+$", rownames(expr_bulk)))) {
    return(expr_bulk)
  }
  
  chip_db_bulk <- get(chip_pkg_bulk, envir = asNamespace(chip_pkg_bulk))
  
  if (mapping_mode_bulk == "legacy_select") {
    mp_bulk <- AnnotationDbi::select(
      x       = chip_db_bulk,
      keys    = rownames(expr_bulk),
      columns = "ENTREZID",
      keytype = "PROBEID"
    )
    mp_bulk <- mp_bulk[!is.na(mp_bulk$ENTREZID), , drop = FALSE]
    
    if (nrow(mp_bulk) == 0) {
      stop("No probe -> Entrez mappings were retained for ", chip_pkg_bulk)
    }
    
    expr_sub_bulk <- expr_bulk[mp_bulk$PROBEID, , drop = FALSE]
    rownames(expr_sub_bulk) <- mp_bulk$ENTREZID
    
    expr_sum_bulk <- rowsum(expr_sub_bulk, group = rownames(expr_sub_bulk), reorder = TRUE)
    expr_n_bulk   <- table(factor(rownames(expr_sub_bulk), levels = rownames(expr_sum_bulk)))
    expr_out_bulk <- expr_sum_bulk / as.numeric(expr_n_bulk)
    
    return(expr_out_bulk)
  }
  
  multi_vals_bulk <- if (mapping_mode_bulk == "mapids_first") "first" else "asNA"
  
  entrez_bulk <- AnnotationDbi::mapIds(
    x         = chip_db_bulk,
    keys      = rownames(expr_bulk),
    column    = "ENTREZID",
    keytype   = "PROBEID",
    multiVals = multi_vals_bulk
  )
  
  keep_bulk <- !is.na(entrez_bulk)
  
  if (!any(keep_bulk)) {
    stop("No probe -> Entrez mappings were retained for ", chip_pkg_bulk)
  }
  
  expr_sub_bulk <- expr_bulk[keep_bulk, , drop = FALSE]
  rownames(expr_sub_bulk) <- base::unname(entrez_bulk[keep_bulk])
  
  expr_sum_bulk <- rowsum(expr_sub_bulk, group = rownames(expr_sub_bulk), reorder = TRUE)
  expr_n_bulk   <- table(factor(rownames(expr_sub_bulk), levels = rownames(expr_sum_bulk)))
  expr_out_bulk <- expr_sum_bulk / as.numeric(expr_n_bulk)
  
  expr_out_bulk
}

plotPCA_micro_bulk <- function(expr_bulk, group_bulk, cohort_bulk, main_bulk = "PCA after ComBat") {
  pc_bulk <- prcomp(t(expr_bulk))
  
  col_vec_bulk <- c(Early = "#2C7BB6", Late = "#D7191C")[group_bulk]
  pch_map_bulk <- c(C1 = 17, C2 = 7, C3 = 3, C4 = 5)
  pch_vec_bulk <- pch_map_bulk[as.character(cohort_bulk)]
  
  plot(
    x    = pc_bulk$x[, 1],
    y    = pc_bulk$x[, 2],
    col  = col_vec_bulk,
    pch  = pch_vec_bulk,
    xlab = sprintf("PC1 (%.1f%%)", 100 * pc_bulk$sdev[1]^2 / sum(pc_bulk$sdev^2)),
    ylab = sprintf("PC2 (%.1f%%)", 100 * pc_bulk$sdev[2]^2 / sum(pc_bulk$sdev^2)),
    main = main_bulk
  )
  
  legend(
    "topleft",
    legend = levels(group_bulk),
    col    = c("#2C7BB6", "#D7191C"),
    pch    = 16,
    title  = "Group",
    bty    = "n"
  )
  
  legend(
    "topright",
    legend = levels(cohort_bulk),
    pch    = pch_map_bulk[levels(cohort_bulk)],
    title  = "Cohort",
    bty    = "n"
  )
}

plotRLE_micro_bulk <- function(expr_bulk, group_bulk, cohort_bulk, main_bulk = "RLE after ComBat") {
  med_bulk <- matrixStats::rowMedians(expr_bulk, na.rm = TRUE)
  rle_bulk <- sweep(expr_bulk, 1, med_bulk, "-")
  
  boxplot(
    as.data.frame(rle_bulk),
    las     = 2,
    outline = FALSE,
    col     = c(C1 = "black", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")[cohort_bulk],
    ylab    = "Relative log expression",
    main    = main_bulk
  )
  
  legend(
    "topright",
    legend = levels(cohort_bulk),
    fill   = c(C1 = "black", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")[levels(cohort_bulk)],
    title  = "Cohort",
    bty    = "n"
  )
}

compare_matrix_bulk <- function(
    mat_new_bulk,
    mat_old_bulk,
    label_new_bulk,
    label_old_bulk,
    tol_bulk = 1e-10
) {
  if (is.null(mat_old_bulk)) {
    return(data.frame(
      label_new_bulk      = label_new_bulk,
      label_old_bulk      = label_old_bulk,
      available_old_bulk  = FALSE,
      same_dim_bulk       = NA,
      same_rownames_bulk  = NA,
      same_colnames_bulk  = NA,
      n_common_rows_bulk  = 0,
      n_common_cols_bulk  = 0,
      one_side_na_n_bulk  = NA,
      max_abs_diff_bulk   = NA_real_,
      mean_abs_diff_bulk  = NA_real_,
      all_within_tol_bulk = NA,
      exact_match_bulk    = NA,
      stringsAsFactors    = FALSE
    ))
  }
  
  same_dim_bulk      <- identical(dim(mat_new_bulk), dim(mat_old_bulk))
  same_rownames_bulk <- identical(rownames(mat_new_bulk), rownames(mat_old_bulk))
  same_colnames_bulk <- identical(colnames(mat_new_bulk), colnames(mat_old_bulk))
  
  common_rows_bulk <- base::intersect(rownames(mat_new_bulk), rownames(mat_old_bulk))
  common_cols_bulk <- base::intersect(colnames(mat_new_bulk), colnames(mat_old_bulk))
  
  one_side_na_n_bulk <- NA_integer_
  max_abs_diff_bulk  <- NA_real_
  mean_abs_diff_bulk <- NA_real_
  all_within_tol_bulk<- NA
  
  if (length(common_rows_bulk) > 0 && length(common_cols_bulk) > 0) {
    x_bulk <- mat_new_bulk[common_rows_bulk, common_cols_bulk, drop = FALSE]
    y_bulk <- mat_old_bulk[common_rows_bulk, common_cols_bulk, drop = FALSE]
    
    one_side_na_bulk <- xor(is.na(x_bulk), is.na(y_bulk))
    both_obs_bulk    <- !(is.na(x_bulk) | is.na(y_bulk))
    
    one_side_na_n_bulk <- sum(one_side_na_bulk)
    
    if (any(both_obs_bulk)) {
      diff_bulk <- abs(x_bulk[both_obs_bulk] - y_bulk[both_obs_bulk])
      max_abs_diff_bulk  <- max(diff_bulk, na.rm = TRUE)
      mean_abs_diff_bulk <- mean(diff_bulk, na.rm = TRUE)
      all_within_tol_bulk<- (one_side_na_n_bulk == 0) && all(diff_bulk <= tol_bulk)
    } else {
      max_abs_diff_bulk  <- 0
      mean_abs_diff_bulk <- 0
      all_within_tol_bulk<- (one_side_na_n_bulk == 0)
    }
  }
  
  exact_match_bulk <- same_dim_bulk &&
    same_rownames_bulk &&
    same_colnames_bulk &&
    isTRUE(all_within_tol_bulk)
  
  data.frame(
    label_new_bulk      = label_new_bulk,
    label_old_bulk      = label_old_bulk,
    available_old_bulk  = TRUE,
    same_dim_bulk       = same_dim_bulk,
    same_rownames_bulk  = same_rownames_bulk,
    same_colnames_bulk  = same_colnames_bulk,
    n_common_rows_bulk  = length(common_rows_bulk),
    n_common_cols_bulk  = length(common_cols_bulk),
    one_side_na_n_bulk  = one_side_na_n_bulk,
    max_abs_diff_bulk   = max_abs_diff_bulk,
    mean_abs_diff_bulk  = mean_abs_diff_bulk,
    all_within_tol_bulk = all_within_tol_bulk,
    exact_match_bulk    = exact_match_bulk,
    stringsAsFactors    = FALSE
  )
}

compare_deg_tables_bulk <- function(
    deg_new_bulk,
    deg_old_bulk,
    tol_bulk = 1e-10
) {
  if (is.null(deg_old_bulk)) {
    return(data.frame(
      available_old_deg_bulk      = FALSE,
      n_new_deg_bulk              = nrow(deg_new_bulk),
      n_old_deg_bulk              = NA_integer_,
      same_n_deg_bulk             = NA,
      same_row_order_deg_bulk     = NA,
      same_row_set_deg_bulk       = NA,
      n_common_deg_bulk           = 0,
      jaccard_deg_bulk            = NA_real_,
      precision_deg_bulk          = NA_real_,
      recall_deg_bulk             = NA_real_,
      direction_cons_deg_bulk     = NA_real_,
      max_abs_logFC_diff_bulk     = NA_real_,
      max_abs_adjP_diff_bulk      = NA_real_,
      exact_table_match_deg_bulk  = NA,
      stringsAsFactors            = FALSE
    ))
  }
  
  ids_new_bulk <- rownames(deg_new_bulk)
  ids_old_bulk <- rownames(deg_old_bulk)
  
  common_ids_bulk <- base::intersect(ids_new_bulk, ids_old_bulk)
  union_ids_bulk  <- base::union(ids_new_bulk, ids_old_bulk)
  
  same_n_deg_bulk         <- nrow(deg_new_bulk) == nrow(deg_old_bulk)
  same_row_order_deg_bulk <- identical(ids_new_bulk, ids_old_bulk)
  same_row_set_deg_bulk   <- base::setequal(ids_new_bulk, ids_old_bulk)
  
  jaccard_deg_bulk   <- if (length(union_ids_bulk) == 0) NA_real_ else length(common_ids_bulk) / length(union_ids_bulk)
  precision_deg_bulk <- if (nrow(deg_new_bulk) == 0) NA_real_ else length(common_ids_bulk) / nrow(deg_new_bulk)
  recall_deg_bulk    <- if (nrow(deg_old_bulk) == 0) NA_real_ else length(common_ids_bulk) / nrow(deg_old_bulk)
  
  direction_cons_deg_bulk <- NA_real_
  if (length(common_ids_bulk) > 0 &&
      "logFC" %in% colnames(deg_new_bulk) &&
      "logFC" %in% colnames(deg_old_bulk)) {
    direction_cons_deg_bulk <- mean(
      sign(deg_new_bulk[common_ids_bulk, "logFC"]) ==
        sign(deg_old_bulk[common_ids_bulk, "logFC"])
    )
  }
  
  max_abs_logFC_diff_bulk <- NA_real_
  max_abs_adjP_diff_bulk  <- NA_real_
  
  if (length(common_ids_bulk) > 0) {
    if ("logFC" %in% colnames(deg_new_bulk) && "logFC" %in% colnames(deg_old_bulk)) {
      max_abs_logFC_diff_bulk <- max(
        abs(deg_new_bulk[common_ids_bulk, "logFC"] - deg_old_bulk[common_ids_bulk, "logFC"]),
        na.rm = TRUE
      )
    }
    if ("adj.P.Val" %in% colnames(deg_new_bulk) && "adj.P.Val" %in% colnames(deg_old_bulk)) {
      max_abs_adjP_diff_bulk <- max(
        abs(deg_new_bulk[common_ids_bulk, "adj.P.Val"] - deg_old_bulk[common_ids_bulk, "adj.P.Val"]),
        na.rm = TRUE
      )
    }
  }
  
  exact_table_match_deg_bulk <- same_n_deg_bulk &&
    same_row_order_deg_bulk &&
    isTRUE(max_abs_logFC_diff_bulk <= tol_bulk || is.na(max_abs_logFC_diff_bulk)) &&
    isTRUE(max_abs_adjP_diff_bulk  <= tol_bulk || is.na(max_abs_adjP_diff_bulk))
  
  data.frame(
    available_old_deg_bulk      = TRUE,
    n_new_deg_bulk              = nrow(deg_new_bulk),
    n_old_deg_bulk              = nrow(deg_old_bulk),
    same_n_deg_bulk             = same_n_deg_bulk,
    same_row_order_deg_bulk     = same_row_order_deg_bulk,
    same_row_set_deg_bulk       = same_row_set_deg_bulk,
    n_common_deg_bulk           = length(common_ids_bulk),
    jaccard_deg_bulk            = jaccard_deg_bulk,
    precision_deg_bulk          = precision_deg_bulk,
    recall_deg_bulk             = recall_deg_bulk,
    direction_cons_deg_bulk     = direction_cons_deg_bulk,
    max_abs_logFC_diff_bulk     = max_abs_logFC_diff_bulk,
    max_abs_adjP_diff_bulk      = max_abs_adjP_diff_bulk,
    exact_table_match_deg_bulk  = exact_table_match_deg_bulk,
    stringsAsFactors            = FALSE
  )
}

read_stage_component_bulk <- function(stage_path_bulk, component_name_bulk) {
  if (!file.exists(stage_path_bulk)) {
    return(NULL)
  }
  stage_obj_bulk <- readRDS(stage_path_bulk)
  if (is.list(stage_obj_bulk) && component_name_bulk %in% names(stage_obj_bulk)) {
    return(stage_obj_bulk[[component_name_bulk]])
  }
  NULL
}


#  4. Sample tables--------------------------------------------------------------

plus_tbl_bulk <- data.frame(
  PID = sprintf("P%02d", 1:17),
  GSM = c(
    "GSM2934819","GSM2934820","GSM2934821","GSM2934822","GSM2934823",
    "GSM2934824","GSM2934825","GSM2934826","GSM2934827","GSM2934828",
    "GSM2934829","GSM2934830","GSM2934831","GSM2934832","GSM2934833",
    "GSM2934834","GSM2934835"
  ),
  Age = c(7, 0.9, 4, 1.6, 4, 8, 5, 6, 1.9, 4, 3, 3, 1.9, 1, 2, 3.5, 7)
)

ab_tbl_bulk <- data.frame(
  PID   = sprintf("P%02d", 18:24),
  GSM_A = c(
    "GSM74377","GSM74378","GSM74379","GSM74380",
    "GSM121357","GSM121361","GSM121363"
  ),
  GSM_B = c(
    "GSM120786","GSM120777","GSM120763","GSM120760",
    "GSM121329","GSM121331","GSM121333"
  ),
  Age   = c(9, 8, 7, 6.5, 5, 9, 7)
)

gse1764_tbl_bulk <- data.frame(
  PID   = sprintf("P%02d", 25:27),
  GSM_A = c("GSM30669","GSM30670","GSM30671"),
  GSM_B = c("GSM30675","GSM30676","GSM30677"),
  Age   = c(8.3, 10.4, 16.7)
)

u95_tbl_bulk <- data.frame(
  PID   = sprintf("P%02d", 28:32),
  GSM_A = c("GSM15833","GSM15834","GSM15836","GSM15838","GSM15839"),
  GSM_B = c("GSM15923","GSM15924","GSM15926","GSM15928","GSM15929"),
  GSM_C = c("GSM16215","GSM16216","GSM16218","GSM16219","GSM16220"),
  GSM_D = c("GSM16237","GSM16244","GSM16245","GSM16249","GSM16250"),
  GSM_E = c("GSM16299","GSM16302","GSM16304","GSM16305","GSM16306"),
  Age   = c(1.0, 1.5, 3.0, 1.0, 0.8)
)

meta_c1_bulk <- data.frame(
  PID    = plus_tbl_bulk$PID,
  GSM    = plus_tbl_bulk$GSM,
  Age    = plus_tbl_bulk$Age,
  Cohort = "C1"
)

meta_c2_bulk <- data.frame(
  PID    = ab_tbl_bulk$PID,
  GSM    = ab_tbl_bulk$GSM_A,
  Age    = ab_tbl_bulk$Age,
  Cohort = "C2"
)

meta_c3_bulk <- data.frame(
  PID    = gse1764_tbl_bulk$PID,
  GSM    = gse1764_tbl_bulk$GSM_A,
  Age    = gse1764_tbl_bulk$Age,
  Cohort = "C3"
)

meta_c4_bulk <- data.frame(
  PID    = u95_tbl_bulk$PID,
  GSM    = u95_tbl_bulk$GSM_A,
  Age    = u95_tbl_bulk$Age,
  Cohort = "C4"
)

sample_meta_bulk <- rbind(
  meta_c1_bulk,
  meta_c2_bulk,
  meta_c3_bulk,
  meta_c4_bulk
)


#  5. Read local series-matrix files---------------------------------------------

message("===== Bulk pipeline: reading series-matrix files =====")

c1_file_bulk <- "GSE109178_series_matrix.txt.gz"
c2a_file_bulk <- "GSE3307-GPL96_series_matrix.txt.gz"
c2b_file_bulk <- "GSE3307-GPL97_series_matrix.txt.gz"
c3a_file_bulk <- "GSE1764-GPL96_series_matrix.txt.gz"
c3b_file_bulk <- "GSE1764-GPL97_series_matrix.txt.gz"

u95_files_bulk <- c(
  A = "GSE1004-GPL8300_series_matrix.txt.gz",
  B = "GSE1007-GPL92_series_matrix.txt.gz",
  C = "GSE1007-GPL93_series_matrix.txt.gz",
  D = "GSE1007-GPL94_series_matrix.txt.gz",
  E = "GSE1007-GPL95_series_matrix.txt.gz"
)

eset_c1_bulk  <- read_sm_bulk(c1_file_bulk,  raw_dir_bulk)
eset_c2a_bulk <- read_sm_bulk(c2a_file_bulk, raw_dir_bulk)
eset_c2b_bulk <- read_sm_bulk(c2b_file_bulk, raw_dir_bulk)
eset_c3a_bulk <- read_sm_bulk(c3a_file_bulk, raw_dir_bulk)
eset_c3b_bulk <- read_sm_bulk(c3b_file_bulk, raw_dir_bulk)

eset_u95_bulk <- lapply(u95_files_bulk, function(fname_bulk) {
  read_sm_bulk(fname_bulk, raw_dir_bulk)
})

idx_c1_bulk <- match_gsm_indices_bulk(
  expected_gsm_bulk        = plus_tbl_bulk$GSM,
  observed_sample_names_bulk = Biobase::sampleNames(eset_c1_bulk),
  label_bulk               = "C1 / GSE109178"
)

idx_c2a_bulk <- match_gsm_indices_bulk(
  expected_gsm_bulk        = ab_tbl_bulk$GSM_A,
  observed_sample_names_bulk = Biobase::sampleNames(eset_c2a_bulk),
  label_bulk               = "C2A / GSE3307 GPL96"
)

idx_c2b_bulk <- match_gsm_indices_bulk(
  expected_gsm_bulk        = ab_tbl_bulk$GSM_B,
  observed_sample_names_bulk = Biobase::sampleNames(eset_c2b_bulk),
  label_bulk               = "C2B / GSE3307 GPL97"
)

idx_c3a_bulk <- match_gsm_indices_bulk(
  expected_gsm_bulk        = gse1764_tbl_bulk$GSM_A,
  observed_sample_names_bulk = Biobase::sampleNames(eset_c3a_bulk),
  label_bulk               = "C3A / GSE1764 GPL96"
)

idx_c3b_bulk <- match_gsm_indices_bulk(
  expected_gsm_bulk        = gse1764_tbl_bulk$GSM_B,
  observed_sample_names_bulk = Biobase::sampleNames(eset_c3b_bulk),
  label_bulk               = "C3B / GSE1764 GPL97"
)

mat_c1_raw_bulk <- Biobase::exprs(eset_c1_bulk)[, idx_c1_bulk, drop = FALSE]
if (needs_log2_bulk(mat_c1_raw_bulk)) {
  mat_c1_raw_bulk <- log2(mat_c1_raw_bulk + 1)
}

mat_c2a_raw_bulk <- Biobase::exprs(eset_c2a_bulk)[, idx_c2a_bulk, drop = FALSE]
if (needs_log2_bulk(mat_c2a_raw_bulk)) {
  mat_c2a_raw_bulk <- log2(mat_c2a_raw_bulk + 1)
}

mat_c2b_raw_bulk <- Biobase::exprs(eset_c2b_bulk)[, idx_c2b_bulk, drop = FALSE]
if (needs_log2_bulk(mat_c2b_raw_bulk)) {
  mat_c2b_raw_bulk <- log2(mat_c2b_raw_bulk + 1)
}

mat_c3a_raw_bulk <- Biobase::exprs(eset_c3a_bulk)[, idx_c3a_bulk, drop = FALSE]
if (needs_log2_bulk(mat_c3a_raw_bulk)) {
  mat_c3a_raw_bulk <- log2(mat_c3a_raw_bulk + 1)
}

mat_c3b_raw_bulk <- Biobase::exprs(eset_c3b_bulk)[, idx_c3b_bulk, drop = FALSE]
if (needs_log2_bulk(mat_c3b_raw_bulk)) {
  mat_c3b_raw_bulk <- log2(mat_c3b_raw_bulk + 1)
}

mat_u95_raw_list_bulk <- list()

for (chip_bulk in names(u95_files_bulk)) {
  gsm_vec_bulk <- u95_tbl_bulk[[paste0("GSM_", chip_bulk)]]
  idx_u95_bulk <- match_gsm_indices_bulk(
    expected_gsm_bulk          = gsm_vec_bulk,
    observed_sample_names_bulk = Biobase::sampleNames(eset_u95_bulk[[chip_bulk]]),
    label_bulk                 = paste0("C4 / U95-", chip_bulk)
  )
  
  mat_tmp_bulk <- Biobase::exprs(eset_u95_bulk[[chip_bulk]])[, idx_u95_bulk, drop = FALSE]
  if (needs_log2_bulk(mat_tmp_bulk)) {
    mat_tmp_bulk <- log2(mat_tmp_bulk + 1)
  }
  mat_u95_raw_list_bulk[[chip_bulk]] <- mat_tmp_bulk
}


#  6. Probe -> Entrez mapping----------------------------------------------------

message("===== Bulk pipeline: probe -> Entrez mapping =====")

expr_c1_entrez_bulk <- map_to_entrez_bulk(
  expr_bulk         = mat_c1_raw_bulk,
  chip_pkg_bulk     = "hgu133plus2.db",
  mapping_mode_bulk = mapping_mode_bulk
)

expr_c2a_entrez_bulk <- map_to_entrez_bulk(
  expr_bulk         = mat_c2a_raw_bulk,
  chip_pkg_bulk     = "hgu133a.db",
  mapping_mode_bulk = mapping_mode_bulk
)

expr_c2b_entrez_bulk <- map_to_entrez_bulk(
  expr_bulk         = mat_c2b_raw_bulk,
  chip_pkg_bulk     = "hgu133b.db",
  mapping_mode_bulk = mapping_mode_bulk
)

expr_c3a_entrez_bulk <- map_to_entrez_bulk(
  expr_bulk         = mat_c3a_raw_bulk,
  chip_pkg_bulk     = "hgu133a.db",
  mapping_mode_bulk = mapping_mode_bulk
)

expr_c3b_entrez_bulk <- map_to_entrez_bulk(
  expr_bulk         = mat_c3b_raw_bulk,
  chip_pkg_bulk     = "hgu133b.db",
  mapping_mode_bulk = mapping_mode_bulk
)

expr_u95a_entrez_bulk <- map_to_entrez_bulk(
  expr_bulk         = mat_u95_raw_list_bulk[["A"]],
  chip_pkg_bulk     = "hgu95av2.db",
  mapping_mode_bulk = mapping_mode_bulk
)

expr_u95b_entrez_bulk <- map_to_entrez_bulk(
  expr_bulk         = mat_u95_raw_list_bulk[["B"]],
  chip_pkg_bulk     = "hgu95b.db",
  mapping_mode_bulk = mapping_mode_bulk
)

expr_u95c_entrez_bulk <- map_to_entrez_bulk(
  expr_bulk         = mat_u95_raw_list_bulk[["C"]],
  chip_pkg_bulk     = "hgu95c.db",
  mapping_mode_bulk = mapping_mode_bulk
)

expr_u95d_entrez_bulk <- map_to_entrez_bulk(
  expr_bulk         = mat_u95_raw_list_bulk[["D"]],
  chip_pkg_bulk     = "hgu95d.db",
  mapping_mode_bulk = mapping_mode_bulk
)

expr_u95e_entrez_bulk <- map_to_entrez_bulk(
  expr_bulk         = mat_u95_raw_list_bulk[["E"]],
  chip_pkg_bulk     = "hgu95e.db",
  mapping_mode_bulk = mapping_mode_bulk
)


#  7. Rename columns by PID and build patient-level cohort matrices--------------

message("===== Bulk pipeline: patient-level averaging =====")

expr_c1_pid_bulk <- rename_cols_by_pid_bulk(
  mat_bulk     = expr_c1_entrez_bulk,
  gsm_vec_bulk = plus_tbl_bulk$GSM,
  pid_vec_bulk = plus_tbl_bulk$PID
)

expr_c2a_pid_bulk <- rename_cols_by_pid_bulk(
  mat_bulk     = expr_c2a_entrez_bulk,
  gsm_vec_bulk = ab_tbl_bulk$GSM_A,
  pid_vec_bulk = ab_tbl_bulk$PID
)

expr_c2b_pid_bulk <- rename_cols_by_pid_bulk(
  mat_bulk     = expr_c2b_entrez_bulk,
  gsm_vec_bulk = ab_tbl_bulk$GSM_B,
  pid_vec_bulk = ab_tbl_bulk$PID
)

expr_c3a_pid_bulk <- rename_cols_by_pid_bulk(
  mat_bulk     = expr_c3a_entrez_bulk,
  gsm_vec_bulk = gse1764_tbl_bulk$GSM_A,
  pid_vec_bulk = gse1764_tbl_bulk$PID
)

expr_c3b_pid_bulk <- rename_cols_by_pid_bulk(
  mat_bulk     = expr_c3b_entrez_bulk,
  gsm_vec_bulk = gse1764_tbl_bulk$GSM_B,
  pid_vec_bulk = gse1764_tbl_bulk$PID
)

expr_u95a_pid_bulk <- rename_cols_by_pid_bulk(
  mat_bulk     = expr_u95a_entrez_bulk,
  gsm_vec_bulk = u95_tbl_bulk$GSM_A,
  pid_vec_bulk = u95_tbl_bulk$PID
)

expr_u95b_pid_bulk <- rename_cols_by_pid_bulk(
  mat_bulk     = expr_u95b_entrez_bulk,
  gsm_vec_bulk = u95_tbl_bulk$GSM_B,
  pid_vec_bulk = u95_tbl_bulk$PID
)

expr_u95c_pid_bulk <- rename_cols_by_pid_bulk(
  mat_bulk     = expr_u95c_entrez_bulk,
  gsm_vec_bulk = u95_tbl_bulk$GSM_C,
  pid_vec_bulk = u95_tbl_bulk$PID
)

expr_u95d_pid_bulk <- rename_cols_by_pid_bulk(
  mat_bulk     = expr_u95d_entrez_bulk,
  gsm_vec_bulk = u95_tbl_bulk$GSM_D,
  pid_vec_bulk = u95_tbl_bulk$PID
)

expr_u95e_pid_bulk <- rename_cols_by_pid_bulk(
  mat_bulk     = expr_u95e_entrez_bulk,
  gsm_vec_bulk = u95_tbl_bulk$GSM_E,
  pid_vec_bulk = u95_tbl_bulk$PID
)

expr_C1_bulk <- expr_c1_pid_bulk
expr_C2_bulk <- avg_two_bulk(expr_c2a_pid_bulk, expr_c2b_pid_bulk)
expr_C3_bulk <- avg_two_bulk(expr_c3a_pid_bulk, expr_c3b_pid_bulk)

expr_C4_bulk <- avg_multi_bulk(
  list(
    expr_u95a_pid_bulk,
    expr_u95b_pid_bulk,
    expr_u95c_pid_bulk,
    expr_u95d_pid_bulk,
    expr_u95e_pid_bulk
  )
)

stopifnot(identical(sort(colnames(expr_C1_bulk)), sort(plus_tbl_bulk$PID)))
stopifnot(identical(sort(colnames(expr_C2_bulk)), sort(ab_tbl_bulk$PID)))
stopifnot(identical(sort(colnames(expr_C3_bulk)), sort(gse1764_tbl_bulk$PID)))
stopifnot(identical(sort(colnames(expr_C4_bulk)), sort(u95_tbl_bulk$PID)))


#  8. Merge the 4 cohorts--------------------------------------------------------

message("===== Bulk pipeline: 4-cohort merge =====")

cohort_list_bulk <- list(
  C1 = expr_C1_bulk,
  C2 = expr_C2_bulk,
  C3 = expr_C3_bulk,
  C4 = expr_C4_bulk
)

genes_all_bulk <- base::Reduce(base::union, lapply(cohort_list_bulk, rownames))
genes_common_bulk <- base::Reduce(base::intersect, lapply(cohort_list_bulk, rownames))

expr_mat_bulk <- do.call(
  cbind,
  lapply(cohort_list_bulk, function(m_bulk) {
    pad_rows_bulk(m_bulk, genes_all_bulk)
  })
)

PID_bulk <- colnames(expr_mat_bulk)

Cohort_bulk <- factor(
  rep(names(cohort_list_bulk), times = vapply(cohort_list_bulk, ncol, integer(1))),
  levels = c("C1", "C2", "C3", "C4")
)

meta_aligned_bulk <- sample_meta_bulk[match(PID_bulk, sample_meta_bulk$PID), , drop = FALSE]
stopifnot(all(meta_aligned_bulk$PID == PID_bulk))

Age_bulk <- meta_aligned_bulk$Age
names(Age_bulk) <- meta_aligned_bulk$PID

Group_bulk <- factor(
  ifelse(Age_bulk[PID_bulk] < age_cutoff_bulk, "Early", "Late"),
  levels = c("Early", "Late")
)

stopifnot(length(find_dup_cols_bulk(expr_mat_bulk)) == 0)

gene_measured_n_bulk <- rowSums(!is.na(expr_mat_bulk))
na_total_bulk        <- sum(is.na(expr_mat_bulk))
na_fraction_bulk     <- na_total_bulk / prod(dim(expr_mat_bulk))

message("Merged matrix bulk dimensions: ", nrow(expr_mat_bulk), " genes x ", ncol(expr_mat_bulk), " samples")
message("Merged matrix NA fraction: ", sprintf("%.4f", na_fraction_bulk))


#  9. Within-cohort quantile normalization---------------------------------------

message("===== Bulk pipeline: within-cohort QN =====")

expr_filt_bulk <- expr_mat_bulk[gene_measured_n_bulk >= min_measured_samples_bulk, , drop = FALSE]

PID_filt_bulk <- colnames(expr_filt_bulk)
Cohort_filt_bulk <- Cohort_bulk[match(PID_filt_bulk, PID_bulk)]
stopifnot(all(!is.na(Cohort_filt_bulk)))
stopifnot(identical(as.character(Cohort_filt_bulk), as.character(Cohort_bulk[match(PID_filt_bulk, PID_bulk)])))

genes_filt_bulk <- rownames(expr_filt_bulk)

expr_QN_list_bulk <- lapply(levels(Cohort_filt_bulk), function(cohort_name_bulk) {
  pid_subset_bulk <- PID_filt_bulk[Cohort_filt_bulk == cohort_name_bulk]
  
  mat_subset_bulk <- expr_filt_bulk[, pid_subset_bulk, drop = FALSE]
  mat_subset_bulk <- mat_subset_bulk[rowSums(!is.na(mat_subset_bulk)) > 0, , drop = FALSE]
  
  mat_qn_bulk <- limma::normalizeBetweenArrays(mat_subset_bulk, method = "quantile")
  pad_rows_bulk(mat_qn_bulk, genes_filt_bulk)
})

names(expr_QN_list_bulk) <- levels(Cohort_filt_bulk)

expr_QN_bulk <- do.call(cbind, expr_QN_list_bulk)
expr_QN_bulk <- expr_QN_bulk[, PID_filt_bulk, drop = FALSE]

stopifnot(identical(rownames(expr_QN_bulk), genes_filt_bulk))
stopifnot(identical(colnames(expr_QN_bulk), PID_filt_bulk))


# 10. Gene selection across cohorts + kNN imputation----------------------------

message("===== Bulk pipeline: gene selection + kNN imputation =====")

keep_present_all_cohorts_bulk <- base::Reduce(
  `&`,
  lapply(levels(Cohort_filt_bulk), function(cohort_name_bulk) {
    rowSums(!is.na(expr_QN_bulk[, Cohort_filt_bulk == cohort_name_bulk, drop = FALSE])) >=
      min_present_each_cohort_bulk
  })
)

expr_sel_bulk <- expr_QN_bulk[keep_present_all_cohorts_bulk, , drop = FALSE]

set.seed(seed_bulk)
expr_imp_knn_bulk <- impute::impute.knn(
  data   = as.matrix(expr_sel_bulk),
  rowmax = knn_rowmax_bulk,
  colmax = knn_colmax_bulk
)$data

na_pos_bulk <- which(is.na(expr_imp_knn_bulk), arr.ind = TRUE)

if (nrow(na_pos_bulk) > 0) {
  row_mean_bulk <- rowMeans(expr_imp_knn_bulk, na.rm = TRUE)
  expr_imp_knn_bulk[na_pos_bulk] <- row_mean_bulk[na_pos_bulk[, "row"]]
}

expr_imp_bulk <- expr_imp_knn_bulk[rowSums(is.na(expr_imp_knn_bulk)) == 0, , drop = FALSE]

stopifnot(length(find_dup_cols_bulk(expr_imp_bulk)) == 0)
stopifnot(sum(is.na(expr_imp_bulk)) == 0)

message("Imputed matrix bulk dimensions: ", nrow(expr_imp_bulk), " genes x ", ncol(expr_imp_bulk), " samples")


# 11. Save stage caches-------------------------------------------------

stamp_bulk <- list(
  built_at_bulk        = Sys.time(),
  project_bulk         = normalizePath(proj_root_bulk, winslash = "/", mustWork = FALSE),
  age_cutoff_bulk      = age_cutoff_bulk,
  treat_lfc_bulk       = treat_lfc_bulk,
  mapping_mode_bulk    = mapping_mode_bulk,
  min_measured_samples_bulk = min_measured_samples_bulk,
  min_present_each_cohort_bulk = min_present_each_cohort_bulk,
  knn_rowmax_bulk      = knn_rowmax_bulk,
  knn_colmax_bulk      = knn_colmax_bulk,
  combat_mean_only_bulk= combat_mean_only_bulk,
  seed_bulk            = seed_bulk,
  session_bulk         = utils::sessionInfo()
)

if (write_stage_rds_bulk) {
  samples_meta_stage_bulk <- list(
    plus_tbl_bulk        = plus_tbl_bulk,
    ab_tbl_bulk          = ab_tbl_bulk,
    gse1764_tbl_bulk     = gse1764_tbl_bulk,
    u95_tbl_bulk         = u95_tbl_bulk,
    sample_meta_bulk     = sample_meta_bulk,
    meta_aligned_bulk    = meta_aligned_bulk,
    PID_bulk             = PID_bulk,
    Cohort_bulk          = Cohort_bulk,
    Age_bulk             = Age_bulk,
    Group_bulk           = Group_bulk,
    stamp_bulk           = stamp_bulk
  )
  saveRDS(
    samples_meta_stage_bulk,
    file = file.path(cache_dir_bulk, "samples_meta_bulk.rds"),
    compress = "xz",
    version = 3
  )
  
  raw_mats_stage_bulk <- list(
    mat_c1_raw_bulk         = mat_c1_raw_bulk,
    mat_c2a_raw_bulk        = mat_c2a_raw_bulk,
    mat_c2b_raw_bulk        = mat_c2b_raw_bulk,
    mat_c3a_raw_bulk        = mat_c3a_raw_bulk,
    mat_c3b_raw_bulk        = mat_c3b_raw_bulk,
    mat_u95_raw_list_bulk   = mat_u95_raw_list_bulk,
    stamp_bulk              = stamp_bulk
  )
  saveRDS(
    raw_mats_stage_bulk,
    file = file.path(cache_dir_bulk, "raw_mats_bulk.rds"),
    compress = "xz",
    version = 3
  )
  
  entrez_stage_bulk <- list(
    expr_C1_bulk           = expr_C1_bulk,
    expr_C2_bulk           = expr_C2_bulk,
    expr_C3_bulk           = expr_C3_bulk,
    expr_C4_bulk           = expr_C4_bulk,
    cohort_list_bulk       = cohort_list_bulk,
    genes_common_bulk      = genes_common_bulk,
    genes_all_bulk         = genes_all_bulk,
    stamp_bulk             = stamp_bulk
  )
  saveRDS(
    entrez_stage_bulk,
    file = file.path(cache_dir_bulk, "entrez_and_cohorts_bulk.rds"),
    compress = "xz",
    version = 3
  )
  
  merged_stage_bulk <- list(
    expr_mat_bulk          = expr_mat_bulk,
    gene_measured_n_bulk   = gene_measured_n_bulk,
    na_total_bulk          = na_total_bulk,
    na_fraction_bulk       = na_fraction_bulk,
    genes_common_bulk      = genes_common_bulk,
    genes_all_bulk         = genes_all_bulk,
    stamp_bulk             = stamp_bulk
  )
  saveRDS(
    merged_stage_bulk,
    file = file.path(cache_dir_bulk, "merged_expr_mat_bulk.rds"),
    compress = "xz",
    version = 3
  )
  
  normalized_stage_bulk <- list(
    expr_filt_bulk        = expr_filt_bulk,
    Cohort_filt_bulk      = Cohort_filt_bulk,
    expr_QN_bulk          = expr_QN_bulk,
    expr_sel_bulk         = expr_sel_bulk,
    expr_imp_bulk         = expr_imp_bulk,
    stamp_bulk            = stamp_bulk
  )
  saveRDS(
    normalized_stage_bulk,
    file = file.path(cache_dir_bulk, "normalized_imputed_bulk.rds"),
    compress = "xz",
    version = 3
  )
}


# 12. ComBat + limma-treat------------------------------------------------------

message("===== Bulk pipeline: ComBat + limma =====")

meta_imp_bulk <- meta_aligned_bulk[match(colnames(expr_imp_bulk), meta_aligned_bulk$PID), , drop = FALSE]
stopifnot(all(meta_imp_bulk$PID == colnames(expr_imp_bulk)))

Group_imp_bulk <- factor(
  ifelse(meta_imp_bulk$Age < age_cutoff_bulk, "Early", "Late"),
  levels = c("Early", "Late")
)

Cohort_imp_bulk <- factor(
  meta_imp_bulk$Cohort,
  levels = c("C1", "C2", "C3", "C4")
)

design_bulk <- model.matrix(~ 0 + Group_imp_bulk)
colnames(design_bulk) <- levels(Group_imp_bulk)

combat_bulk <- sva::ComBat(
  dat       = expr_imp_bulk,
  batch     = Cohort_imp_bulk,
  mod       = model.matrix(~ Group_imp_bulk),
  par.prior = TRUE,
  mean.only = combat_mean_only_bulk
)

contrast_matrix_bulk <- matrix(
  c(-1, 1),
  ncol     = 1,
  dimnames = list(colnames(design_bulk), "Late_vs_Early_bulk")
)

v_bulk <- limma::vooma(combat_bulk, design_bulk, plot = FALSE)

fit_bulk <- limma::lmFit(v_bulk, design_bulk)
fit_bulk <- limma::contrasts.fit(fit_bulk, contrast_matrix_bulk)
fit_bulk <- limma::eBayes(fit_bulk, trend = TRUE)
fit_treat_bulk <- limma::treat(fit_bulk, lfc = treat_lfc_bulk)

deg_bulk <- limma::topTreat(
  fit_treat_bulk,
  coef    = "Late_vs_Early_bulk",
  p.value = 0.05,
  number  = Inf
)

pc_bulk <- prcomp(t(combat_bulk))$x
rle_center_bulk <- matrixStats::rowMedians(combat_bulk, na.rm = TRUE)
rle_iqr_bulk <- apply(sweep(combat_bulk, 1, rle_center_bulk, "-"), 2, IQR, na.rm = TRUE)

message("DEGs in bulk pipeline: ", nrow(deg_bulk))
# DEGs in bulk pipeline: 358


stage6_bulk <- list(
  combat_bulk         = combat_bulk,
  Group_imp_bulk      = Group_imp_bulk,
  Cohort_imp_bulk     = Cohort_imp_bulk,
  design_bulk         = design_bulk,
  v_bulk              = v_bulk,
  fit_bulk            = fit_bulk,
  fit_treat_bulk      = fit_treat_bulk,
  deg_bulk            = deg_bulk,
  pca_scores_bulk     = pc_bulk,
  rle_iqr_bulk        = rle_iqr_bulk,
  stamp_bulk          = stamp_bulk
)

if (write_stage_rds_bulk) {
  saveRDS(
    stage6_bulk,
    file = file.path(cache_dir_bulk, "stage6_combat_only_bulk.rds"),
    compress = "xz",
    version = 3
  )
}


# 13. Write DEG table and QC plots----------------------------------------------

if (write_deg_table_bulk) {
  deg_out_bulk <- deg_bulk
  deg_out_bulk$EntrezID_bulk <- rownames(deg_out_bulk)
  deg_out_bulk <- deg_out_bulk[, c("EntrezID_bulk", base::setdiff(colnames(deg_out_bulk), "EntrezID_bulk")), drop = FALSE]
  
  utils::write.table(
    deg_out_bulk,
    file      = file.path(res_dir_bulk, "deg_bulk.tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
}

if (write_qc_plots_bulk) {
  pdf(file.path(res_dir_bulk, "qc_pca_rle_bulk.pdf"), width = 11, height = 5.5)
  par(mfrow = c(1, 2))
  plotPCA_micro_bulk(
    expr_bulk   = combat_bulk,
    group_bulk  = Group_imp_bulk,
    cohort_bulk = Cohort_imp_bulk,
    main_bulk   = sprintf("PCA after ComBat (cutoff %s yr)", age_cutoff_bulk)
  )
  plotRLE_micro_bulk(
    expr_bulk   = combat_bulk,
    group_bulk  = Group_imp_bulk,
    cohort_bulk = Cohort_imp_bulk,
    main_bulk   = sprintf("RLE after ComBat (cutoff %s yr)", age_cutoff_bulk)
  )
  dev.off()
}

# 14. PCA plot display------------------------------------------------------------------

proj_root_bulk   <- "C:/DMD_project"
cache_dir_bulk   <- file.path(proj_root_bulk, "cache_bulk_revision")
stage6_path_bulk <- file.path(cache_dir_bulk, "stage6_combat_only_bulk.rds")

assert_file_exists_bulk <- function(path_bulk) {
  if (!file.exists(path_bulk)) stop("Required file was not found: ", path_bulk)
}

get_first_existing_name_bulk <- function(x_bulk, candidates_bulk) {
  hit_bulk <- candidates_bulk[candidates_bulk %in% names(x_bulk)]
  if (length(hit_bulk) == 0) return(NULL)
  hit_bulk[1]
}

assert_file_exists_bulk(stage6_path_bulk)
stage6_bulk <- readRDS(stage6_path_bulk)

combat_bulk <- stage6_bulk$combat_bulk
Group_name_bulk  <- get_first_existing_name_bulk(stage6_bulk, c("Group_imp_bulk", "Group_bulk"))
Cohort_name_bulk <- get_first_existing_name_bulk(stage6_bulk, c("Cohort_imp_bulk", "Cohort_bulk"))

if (is.null(Group_name_bulk) || is.null(Cohort_name_bulk)) {
  stop("Group/Cohort vectors were not found in stage6_combat_only_bulk.rds")
}

Group_plot_bulk  <- factor(stage6_bulk[[Group_name_bulk]], levels = c("Early", "Late"))
Cohort_plot_bulk <- factor(stage6_bulk[[Cohort_name_bulk]], levels = c("C1", "C2", "C3", "C4"))

stopifnot(
  !is.null(combat_bulk),
  ncol(combat_bulk) == length(Group_plot_bulk),
  ncol(combat_bulk) == length(Cohort_plot_bulk)
)

{
  pc_bulk <- stats::prcomp(t(combat_bulk))
  
  col_map_bulk <- c(Early = "#2C7BB6", Late = "#D7191C")
  pch_map_bulk <- c(C1 = 17, C2 = 7, C3 = 3, C4 = 5)
  
  col_vec_bulk <- base::unname(col_map_bulk[as.character(Group_plot_bulk)])
  pch_vec_bulk <- base::unname(pch_map_bulk[as.character(Cohort_plot_bulk)])
  
  graphics::plot(
    x    = pc_bulk$x[, 1],
    y    = pc_bulk$x[, 2],
    col  = col_vec_bulk,
    pch  = pch_vec_bulk,
    xlab = sprintf("PC1 (%.1f%%)", 100 * pc_bulk$sdev[1]^2 / sum(pc_bulk$sdev^2)),
    ylab = sprintf("PC2 (%.1f%%)", 100 * pc_bulk$sdev[2]^2 / sum(pc_bulk$sdev^2)),
    main = "PCA after ComBat"
  )
  
  graphics::legend(
    "topleft",
    legend = levels(Group_plot_bulk),
    col    = col_map_bulk[levels(Group_plot_bulk)],
    pch    = 16,
    title  = "Group",
    bty    = "n"
  )
  
  graphics::legend(
    "topright",
    legend = levels(Cohort_plot_bulk),
    pch    = pch_map_bulk[levels(Cohort_plot_bulk)],
    title  = "Cohort",
    bty    = "n"
  )
}

# 15. RLE plot display------------------------------------------------------------------
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  stop("Package 'matrixStats' is required.")
}
{
  proj_root_bulk   <- "C:/DMD_project"
  cache_dir_bulk   <- file.path(proj_root_bulk, "cache_bulk_revision")
  stage6_path_bulk <- file.path(cache_dir_bulk, "stage6_combat_only_bulk.rds")
  
  assert_file_exists_bulk <- function(path_bulk) {
    if (!file.exists(path_bulk)) stop("Required file was not found: ", path_bulk)
  }
  
  get_first_existing_name_bulk <- function(x_bulk, candidates_bulk) {
    hit_bulk <- candidates_bulk[candidates_bulk %in% names(x_bulk)]
    if (length(hit_bulk) == 0) return(NULL)
    hit_bulk[1]
  }
  
  assert_file_exists_bulk(stage6_path_bulk)
  stage6_bulk <- readRDS(stage6_path_bulk)
  
  combat_bulk <- stage6_bulk$combat_bulk
  Cohort_name_bulk <- get_first_existing_name_bulk(stage6_bulk, c("Cohort_imp_bulk", "Cohort_bulk"))
  
  if (is.null(Cohort_name_bulk)) {
    stop("Cohort vector was not found in stage6_combat_only_bulk.rds")
  }
  
  Cohort_plot_bulk <- factor(stage6_bulk[[Cohort_name_bulk]], levels = c("C1", "C2", "C3", "C4"))
  
  stopifnot(
    !is.null(combat_bulk),
    ncol(combat_bulk) == length(Cohort_plot_bulk)
  )
  
  med_bulk <- matrixStats::rowMedians(combat_bulk, na.rm = TRUE)
  rle_bulk <- sweep(combat_bulk, 1, med_bulk, "-")
  fill_map_bulk <- c(C1 = "black", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")
  
  graphics::boxplot(
    as.data.frame(rle_bulk),
    las     = 2,
    outline = FALSE,
    col     = fill_map_bulk[as.character(Cohort_plot_bulk)],
    ylab    = "Relative log expression",
    main    = "RLE after ComBat"
  )
  
  graphics::legend(
    "topright",
    legend = levels(Cohort_plot_bulk),
    fill   = fill_map_bulk[levels(Cohort_plot_bulk)],
    title  = "Cohort",
    bty    = "n"
  )
}

# 16. Heatmap display-------------------------------------------------------------------

required_pkgs_bulk <- c("matrixStats", "AnnotationDbi", "org.Hs.eg.db", "pheatmap")
missing_pkgs_bulk <- required_pkgs_bulk[!vapply(required_pkgs_bulk, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs_bulk) > 0) {
  stop("Missing packages:\n  ", paste(missing_pkgs_bulk, collapse = "\n  "))
}
{
  proj_root_bulk   <- "C:/DMD_project"
  cache_dir_bulk   <- file.path(proj_root_bulk, "cache_bulk_revision")
  stage6_path_bulk <- file.path(cache_dir_bulk, "stage6_combat_only_bulk.rds")
  
  fdr_top_n_bulk      <- 400
  row_sd_cut_bulk     <- 0.15
  heatmap_clip_bulk   <- 2.5
  show_rownames_bulk  <- FALSE
  
  assert_file_exists_bulk <- function(path_bulk) {
    if (!file.exists(path_bulk)) stop("Required file was not found: ", path_bulk)
  }
  
  get_first_existing_name_bulk <- function(x_bulk, candidates_bulk) {
    hit_bulk <- candidates_bulk[candidates_bulk %in% names(x_bulk)]
    if (length(hit_bulk) == 0) return(NULL)
    hit_bulk[1]
  }
  
  assert_file_exists_bulk(stage6_path_bulk)
  stage6_bulk <- readRDS(stage6_path_bulk)
  
  combat_bulk <- stage6_bulk$combat_bulk
  deg_bulk    <- stage6_bulk$deg_bulk
  Group_name_bulk  <- get_first_existing_name_bulk(stage6_bulk, c("Group_imp_bulk", "Group_bulk"))
  Cohort_name_bulk <- get_first_existing_name_bulk(stage6_bulk, c("Cohort_imp_bulk", "Cohort_bulk"))
  
  if (is.null(Group_name_bulk) || is.null(Cohort_name_bulk)) {
    stop("Group/Cohort vectors were not found in stage6_combat_only_bulk.rds")
  }
  
  Group_plot_bulk  <- factor(stage6_bulk[[Group_name_bulk]], levels = c("Early", "Late"))
  Cohort_plot_bulk <- factor(stage6_bulk[[Cohort_name_bulk]], levels = c("C1", "C2", "C3", "C4"))
  
  stopifnot(
    !is.null(combat_bulk),
    !is.null(deg_bulk),
    ncol(combat_bulk) == length(Group_plot_bulk),
    ncol(combat_bulk) == length(Cohort_plot_bulk)
  )
  
  deg_df_bulk <- as.data.frame(deg_bulk)
  fdr_candidates_bulk <- c("FDR", "adj.P.Val", "padj", "qvalue", "p.adj", "FDR.BH")
  fdr_col_bulk <- base::intersect(fdr_candidates_bulk, colnames(deg_df_bulk))
  if (length(fdr_col_bulk) == 0) {
    stop("No FDR-like column was found in deg_bulk.")
  }
  fdr_col_bulk <- fdr_col_bulk[1]
  
  deg_df_sorted_bulk <- deg_df_bulk[order(deg_df_bulk[[fdr_col_bulk]]), , drop = FALSE]
  top100_id_bulk     <- head(rownames(deg_df_sorted_bulk), fdr_top_n_bulk)
  top100_bulk        <- base::intersect(top100_id_bulk, rownames(combat_bulk))
  
  expr_top_bulk <- combat_bulk[top100_bulk, , drop = FALSE]
  expr_top_bulk <- expr_top_bulk[matrixStats::rowSds(expr_top_bulk) > row_sd_cut_bulk, , drop = FALSE]
  n_gene_bulk   <- nrow(expr_top_bulk)
  
  if (n_gene_bulk == 0) {
    stop("No rows remained after the row SD filter.")
  }
  
  expr_z_bulk <- t(scale(t(expr_top_bulk)))
  expr_z_bulk <- pmin(pmax(expr_z_bulk, -heatmap_clip_bulk), heatmap_clip_bulk)
  
  sym_bulk <- AnnotationDbi::mapIds(
    x         = org.Hs.eg.db::org.Hs.eg.db,
    keys      = rownames(expr_z_bulk),
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )
  
  keep_bulk <- !is.na(sym_bulk) & sym_bulk != ""
  expr_z_bulk <- expr_z_bulk[keep_bulk, , drop = FALSE]
  rownames(expr_z_bulk) <- make.unique(base::unname(sym_bulk[keep_bulk]))
  
  ann_col_bulk <- data.frame(
    Group  = Group_plot_bulk,
    Cohort = factor(Cohort_plot_bulk, levels = c("C1", "C2", "C3", "C4")),
    row.names = colnames(expr_z_bulk),
    stringsAsFactors = FALSE
  )
  
  ann_cols_bulk <- list(
    Group  = c(Early = "#2C7BB6", Late = "#D7191C"),
    Cohort = c(C1 = "black", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")
  )
  
  pheatmap::pheatmap(
    expr_z_bulk,
    scale = "none",
    color = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100),
    breaks = seq(-3.5, 3.5, length.out = 101),
    annotation_col    = ann_col_bulk,
    annotation_colors = ann_cols_bulk,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    show_rownames = show_rownames_bulk,
    border_color  = NA,
    fontsize_row  = 5,
    main = sprintf("DEG heatmap | All DEGs -> %d genes x %d samples",
                   n_gene_bulk, ncol(expr_z_bulk))
  )
}

# 17. Volcano plot display--------------------------------------------------------------

required_pkgs_bulk <- c("ggplot2", "ggrepel", "AnnotationDbi", "org.Hs.eg.db", "limma")
missing_pkgs_bulk <- required_pkgs_bulk[!vapply(required_pkgs_bulk, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs_bulk) > 0) {
  stop("Missing packages:\n  ", paste(missing_pkgs_bulk, collapse = "\n  "))
}

{
  proj_root_bulk   <- "C:/DMD_project"
  cache_dir_bulk   <- file.path(proj_root_bulk, "cache_bulk_revision")
  stage6_path_bulk <- file.path(cache_dir_bulk, "stage6_combat_only_bulk.rds")
  
  lfc_cut_bulk      <- 0.15
  fdr_cut_bulk      <- 0.05
  top_label_n_bulk  <- 10
  label_exclude_regex_bulk <- paste0(
    "^LINC",
    "|^LOC",
    "|^MIR|^MIR\\d",
    "|^SNORD|^SNORA",
    "|^AC\\d|^AL\\d",
    "|-AS\\d+$"
  )
  
  assert_file_exists_bulk <- function(path_bulk) {
    if (!file.exists(path_bulk)) stop("Required file was not found: ", path_bulk)
  }
  
  assert_file_exists_bulk(stage6_path_bulk)
  stage6_bulk <- readRDS(stage6_path_bulk)
  fit_treat_bulk <- stage6_bulk$fit_treat_bulk
  
  if (is.null(fit_treat_bulk)) {
    stop("fit_treat_bulk was not found in stage6_combat_only_bulk.rds")
  }
  
  tbl_all_bulk <- limma::topTreat(fit_treat_bulk, p.value = 1, number = Inf)
  
  tbl_plot_bulk <- as.data.frame(tbl_all_bulk)
  tbl_plot_bulk$EntrezID_bulk <- rownames(tbl_plot_bulk)
  
  sym_bulk <- AnnotationDbi::mapIds(
    x         = org.Hs.eg.db::org.Hs.eg.db,
    keys      = tbl_plot_bulk$EntrezID_bulk,
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )
  
  tbl_plot_bulk$GeneSymbol_bulk <- ifelse(
    is.na(sym_bulk) | sym_bulk == "",
    tbl_plot_bulk$EntrezID_bulk,
    base::unname(sym_bulk)
  )
  
  tbl_plot_bulk$signif_bulk <- ifelse(
    tbl_plot_bulk$adj.P.Val < fdr_cut_bulk & tbl_plot_bulk$logFC >=  lfc_cut_bulk,
    "Up",
    ifelse(
      tbl_plot_bulk$adj.P.Val < fdr_cut_bulk & tbl_plot_bulk$logFC <= -lfc_cut_bulk,
      "Down",
      "NS"
    )
  )
  tbl_plot_bulk$signif_bulk <- factor(tbl_plot_bulk$signif_bulk, levels = c("Up", "Down", "NS"))
  
  top_labs_bulk <- tbl_plot_bulk[
    tbl_plot_bulk$signif_bulk != "NS" &
      !grepl(label_exclude_regex_bulk, tbl_plot_bulk$GeneSymbol_bulk),
    , drop = FALSE
  ]
  
  if (nrow(top_labs_bulk) > 0) {
    top_labs_bulk <- top_labs_bulk[order(top_labs_bulk$adj.P.Val), , drop = FALSE]
    top_labs_bulk <- utils::head(top_labs_bulk, top_label_n_bulk)
  }
  
  p_bulk <- ggplot2::ggplot(tbl_plot_bulk, ggplot2::aes(logFC, -log10(adj.P.Val))) +
    ggplot2::geom_point(ggplot2::aes(colour = signif_bulk), size = 1.6, alpha = 0.8) +
    ggplot2::scale_colour_manual(values = c(Up = "#D7191C", Down = "#2C7BB6", NS = "grey70")) +
    ggplot2::geom_vline(xintercept = c(-lfc_cut_bulk, lfc_cut_bulk), linetype = "dashed", colour = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(fdr_cut_bulk), linetype = "dashed", colour = "grey40") +
    ggrepel::geom_text_repel(
      data = top_labs_bulk,
      ggplot2::aes(label = GeneSymbol_bulk),
      size = 5,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2,
      min.segment.length = 0
    ) +
    ggplot2::labs(
      title = sprintf("Volcano plot | %d DEGs (|logFC| >= %.2f, FDR < %.2f)",
                      sum(tbl_plot_bulk$signif_bulk != "NS"), lfc_cut_bulk, fdr_cut_bulk),
      x = "log2 fold-change",
      y = expression(-log[10]~FDR)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  
  print(p_bulk)
}
# 18. All DEGs display on console----------------------------------------------------------

required_pkgs_bulk <- c("AnnotationDbi", "org.Hs.eg.db")
missing_pkgs_bulk <- required_pkgs_bulk[!vapply(required_pkgs_bulk, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs_bulk) > 0) {
  stop("Missing packages:\n  ", paste(missing_pkgs_bulk, collapse = "\n  "))
}

{
  proj_root_bulk   <- "C:/DMD_project"
  cache_dir_bulk   <- file.path(proj_root_bulk, "cache_bulk_revision")
  stage6_path_bulk <- file.path(cache_dir_bulk, "stage6_combat_only_bulk.rds")
  
  sort_by_fdr_bulk <- TRUE
  
  assert_file_exists_bulk <- function(path_bulk) {
    if (!file.exists(path_bulk)) stop("Required file was not found: ", path_bulk)
  }
  
  assert_file_exists_bulk(stage6_path_bulk)
  stage6_bulk <- readRDS(stage6_path_bulk)
  deg_bulk <- stage6_bulk$deg_bulk
  
  if (is.null(deg_bulk)) {
    stop("deg_bulk was not found in stage6_combat_only_bulk.rds")
  }
  
  out_bulk <- as.data.frame(deg_bulk)
  out_bulk$EntrezID_bulk <- rownames(out_bulk)
  
  sym_bulk <- AnnotationDbi::mapIds(
    x         = org.Hs.eg.db::org.Hs.eg.db,
    keys      = out_bulk$EntrezID_bulk,
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )
  out_bulk$SYMBOL_bulk <- base::unname(sym_bulk)
  
  if (sort_by_fdr_bulk && "adj.P.Val" %in% colnames(out_bulk)) {
    out_bulk <- out_bulk[order(out_bulk$adj.P.Val), , drop = FALSE]
  }
  
  keep_cols_bulk <- c("EntrezID_bulk", "SYMBOL_bulk", "logFC", "adj.P.Val", "P.Value")
  keep_cols_bulk <- keep_cols_bulk[keep_cols_bulk %in% colnames(out_bulk)]
  out_bulk <- out_bulk[, keep_cols_bulk, drop = FALSE]
  
  colnames(out_bulk) <- c("EntrezID", "SYMBOL", "logFC", "FDR", "Pvalue")[seq_len(ncol(out_bulk))]
  
  old_max_print_bulk <- getOption("max.print")
  on.exit(options(max.print = old_max_print_bulk), add = TRUE)
  options(max.print = max(old_max_print_bulk, nrow(out_bulk) * ncol(out_bulk) + 1000))
  print(out_bulk, row.names = FALSE, digits = 3)
}


###─────────────────────────────────────── B: Pathway analysis (Reactome, GO-BP, KEGG) and GSEA ---------------------------------
#  0. User settings--------------------------------------------------------------
{
  proj_root_bulk      <- "C:/DMD_project"
  cache_dir_bulk      <- file.path(proj_root_bulk, "cache_bulk_revision")
  stage6_path_bulk    <- file.path(cache_dir_bulk, "stage6_combat_only_bulk.rds")
  
  ## ORA parameters
  ora_p_cut_bulk      <- 0.05
  ora_q_cut_bulk      <- 0.05
  go_simplify_cut_bulk<- 0.50
  ora_top_n_bulk      <- 20
  
  ## GSEA parameters
  lfc_cutoff_bulk     <- 0.15
  gsea_p_cut_bulk     <- 0.05
  gsea_show_top_bulk  <- 20
  gsea_x_breaks_bulk  <- seq(-2, 2, 1)
  ridge_label_sz_bulk <- 10
  
  ## Console output
  print_gsea_hits_bulk <- TRUE
}

#  1. Packages-------------------------------------------------------------------

suppressPackageStartupMessages({
  library(limma)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(BiocParallel)
  library(enrichplot)
})


#  2. Load stage6 objects if needed--------------------------------------

need_stage6_bulk <- c(
  "combat_bulk",
  "deg_bulk",
  "v_bulk",
  "Group_imp_bulk",
  "fit_treat_bulk"
)

missing_stage6_bulk <- need_stage6_bulk[
  !vapply(need_stage6_bulk, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)
]

if (length(missing_stage6_bulk) > 0) {
  if (!file.exists(stage6_path_bulk)) {
    stop("Revised stage6 file was not found: ", stage6_path_bulk)
  }
  stage6_bulk <- readRDS(stage6_path_bulk)
  list2env(stage6_bulk, envir = .GlobalEnv)
}

stopifnot(
  exists("combat_bulk",     envir = .GlobalEnv, inherits = FALSE),
  exists("deg_bulk",        envir = .GlobalEnv, inherits = FALSE),
  exists("v_bulk",          envir = .GlobalEnv, inherits = FALSE),
  exists("Group_imp_bulk",  envir = .GlobalEnv, inherits = FALSE),
  exists("fit_treat_bulk",  envir = .GlobalEnv, inherits = FALSE)
)

if (!exists("design_bulk", envir = .GlobalEnv, inherits = FALSE)) {
  design_bulk <- model.matrix(~ 0 + Group_imp_bulk)
  colnames(design_bulk) <- levels(Group_imp_bulk)
}


#  3. ORA: KEGG / Reactome / GO-BP-----------------------------------------------

message("===== Revised bulk pathway analysis: ORA =====")

bg_ids_bulk <- unique(sub("_.*", "", rownames(combat_bulk)))
bg_ids_bulk <- bg_ids_bulk[!is.na(bg_ids_bulk) & nzchar(bg_ids_bulk)]

deg_ids_bulk <- unique(sub("_.*", "", rownames(deg_bulk)))
deg_ids_bulk <- deg_ids_bulk[!is.na(deg_ids_bulk) & nzchar(deg_ids_bulk)]

if (length(bg_ids_bulk) == 0L) {
  stop("Background Entrez IDs could not be derived from combat_bulk rownames.")
}
if (length(deg_ids_bulk) == 0L) {
  stop("DEG Entrez IDs could not be derived from deg_bulk rownames.")
}

kegg_res_bulk <- clusterProfiler::enrichKEGG(
  gene         = as.character(deg_ids_bulk),
  universe     = as.character(bg_ids_bulk),
  organism     = "hsa",
  keyType      = "ncbi-geneid",
  pvalueCutoff = ora_p_cut_bulk,
  qvalueCutoff = ora_q_cut_bulk
)

kegg_clean_bulk <- kegg_res_bulk
kegg_clean_df_bulk <- as.data.frame(kegg_res_bulk@result)
if (nrow(kegg_clean_df_bulk) > 0) {
  kegg_clean_df_bulk <- dplyr::filter(
    kegg_clean_df_bulk,
    !grepl("^hsa05|disease", ID, ignore.case = TRUE)
  )
}
kegg_clean_bulk@result <- kegg_clean_df_bulk

react_res_bulk <- ReactomePA::enrichPathway(
  gene         = as.character(deg_ids_bulk),
  universe     = as.character(bg_ids_bulk),
  organism     = "human",
  pvalueCutoff = ora_p_cut_bulk,
  qvalueCutoff = ora_q_cut_bulk
)

go_res_raw_bulk <- clusterProfiler::enrichGO(
  gene          = as.character(deg_ids_bulk),
  universe      = as.character(bg_ids_bulk),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = ora_q_cut_bulk,
  readable      = FALSE
)

go_res_bulk <- if (nrow(as.data.frame(go_res_raw_bulk)) > 0) {
  clusterProfiler::simplify(
    go_res_raw_bulk,
    cutoff     = go_simplify_cut_bulk,
    by         = "p.adjust",
    select_fun = min
  )
} else {
  go_res_raw_bulk
}

pub_dotplot_bulk <- function(eres_bulk,
                             n_bulk = 15,
                             fdr_cut_bulk = 0.05,
                             title_bulk = "",
                             wrap_len_bulk = 35) {
  df_bulk <- as.data.frame(eres_bulk@result)
  df_bulk <- df_bulk[df_bulk$p.adjust <= fdr_cut_bulk, , drop = FALSE]
  
  if (!nrow(df_bulk)) {
    message("<", title_bulk, "> no significant pathway (FDR < ", fdr_cut_bulk, ")")
    return(NULL)
  }
  
  df_bulk <- df_bulk[order(df_bulk$p.adjust), , drop = FALSE]
  df_bulk <- df_bulk[seq_len(min(n_bulk, nrow(df_bulk))), , drop = FALSE]
  df_bulk$GeneRatio_num_bulk <- vapply(
    strsplit(df_bulk$GeneRatio, "/", fixed = TRUE),
    function(x) as.numeric(x[1]) / as.numeric(x[2]),
    numeric(1)
  )
  wrapped_bulk <- stringr::str_wrap(df_bulk$Description, width = wrap_len_bulk)
  df_bulk$Description_bulk <- factor(wrapped_bulk, levels = rev(wrapped_bulk))
  
  ggplot(df_bulk, aes(x = GeneRatio_num_bulk, y = Description_bulk)) +
    geom_point(aes(size = Count, colour = p.adjust)) +
    scale_colour_gradient(low = "#b2182b", high = "#2166ac", name = "FDR") +
    scale_size(range = c(3, 8)) +
    labs(title = title_bulk, x = "Gene ratio", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title  = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9)
    )
}

pub_barplot_bulk <- function(eres_bulk,
                             n_bulk = 15,
                             fdr_cut_bulk = 0.05,
                             title_bulk = "",
                             wrap_len_bulk = 40) {
  df_bulk <- as.data.frame(eres_bulk@result)
  df_bulk <- df_bulk[df_bulk$p.adjust <= fdr_cut_bulk, , drop = FALSE]
  
  if (!nrow(df_bulk)) {
    message("<", title_bulk, "> no significant pathway (FDR < ", fdr_cut_bulk, ")")
    return(NULL)
  }
  
  df_bulk <- df_bulk[order(df_bulk$p.adjust), , drop = FALSE]
  df_bulk <- df_bulk[seq_len(min(n_bulk, nrow(df_bulk))), , drop = FALSE]
  df_bulk$GeneRatio_num_bulk <- vapply(
    strsplit(df_bulk$GeneRatio, "/", fixed = TRUE),
    function(x) as.numeric(x[1]) / as.numeric(x[2]),
    numeric(1)
  )
  wrapped_bulk <- stringr::str_wrap(df_bulk$Description, width = wrap_len_bulk)
  df_bulk$Description_bulk <- factor(wrapped_bulk, levels = rev(wrapped_bulk))
  
  ggplot(df_bulk, aes(x = GeneRatio_num_bulk, y = Description_bulk)) +
    geom_col(aes(fill = p.adjust), width = 0.7) +
    scale_fill_gradient(low = "#b2182b", high = "#2166ac", name = "FDR") +
    labs(title = title_bulk, x = "Gene ratio", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", hjust = 0.5),
      axis.text.y        = element_text(size = 12),
      panel.grid.major.y = element_blank()
    )
}

p_reactome_dot_bulk <- pub_dotplot_bulk(
  react_res_bulk,
  n_bulk       = ora_top_n_bulk,
  fdr_cut_bulk = ora_p_cut_bulk,
  title_bulk   = "Reactome enrichment"
)

p_reactome_bar_bulk <- pub_barplot_bulk(
  react_res_bulk,
  n_bulk       = ora_top_n_bulk,
  fdr_cut_bulk = ora_p_cut_bulk,
  title_bulk   = "Reactome enrichment"
)

p_go_dot_bulk <- pub_dotplot_bulk(
  go_res_bulk,
  n_bulk       = ora_top_n_bulk,
  fdr_cut_bulk = ora_p_cut_bulk,
  title_bulk   = "GO-BP enrichment"
)

p_go_bar_bulk <- pub_barplot_bulk(
  go_res_bulk,
  n_bulk       = ora_top_n_bulk,
  fdr_cut_bulk = ora_p_cut_bulk,
  title_bulk   = "GO-BP enrichment"
)

p_kegg_dot_bulk <- pub_dotplot_bulk(
  kegg_clean_bulk,
  n_bulk       = ora_top_n_bulk,
  fdr_cut_bulk = ora_p_cut_bulk,
  title_bulk   = "KEGG enrichment"
)

p_kegg_bar_bulk <- pub_barplot_bulk(
  kegg_clean_bulk,
  n_bulk       = ora_top_n_bulk,
  fdr_cut_bulk = ora_p_cut_bulk,
  title_bulk   = "KEGG enrichment"
)

reactome_sig_table_bulk <- as.data.frame(react_res_bulk@result)
reactome_sig_table_bulk <- reactome_sig_table_bulk[
  reactome_sig_table_bulk$p.adjust <= ora_p_cut_bulk,
  , drop = FALSE
]
reactome_sig_table_bulk <- reactome_sig_table_bulk[order(reactome_sig_table_bulk$p.adjust), , drop = FALSE]

go_sig_table_bulk <- as.data.frame(go_res_bulk@result)
go_sig_table_bulk <- go_sig_table_bulk[
  go_sig_table_bulk$p.adjust <= ora_p_cut_bulk,
  , drop = FALSE
]
go_sig_table_bulk <- go_sig_table_bulk[order(go_sig_table_bulk$p.adjust), , drop = FALSE]

kegg_sig_table_bulk <- as.data.frame(kegg_clean_bulk@result)
kegg_sig_table_bulk <- kegg_sig_table_bulk[
  kegg_sig_table_bulk$p.adjust <= ora_p_cut_bulk,
  , drop = FALSE
]
kegg_sig_table_bulk <- kegg_sig_table_bulk[order(kegg_sig_table_bulk$p.adjust), , drop = FALSE]

message("Reactome ORA significant pathways: ", nrow(reactome_sig_table_bulk))
message("GO-BP ORA significant pathways  : ", nrow(go_sig_table_bulk))
message("KEGG ORA significant pathways   : ", nrow(kegg_sig_table_bulk))



#  4. GSEA-----------------------------------------------------------------------
message("===== Revised bulk pathway analysis: Reactome GSEA =====")

design_gsea_bulk <- model.matrix(~ 0 + Group_imp_bulk)
colnames(design_gsea_bulk) <- levels(Group_imp_bulk)

fit_lfc_bulk <- v_bulk |>
  limma::lmFit(design_gsea_bulk) |>
  limma::contrasts.fit(cbind(Late_vs_Early = c(-1, 1))) |>
  limma::eBayes(trend = TRUE) |>
  limma::treat(lfc = lfc_cutoff_bulk)

coef_mat_bulk <- fit_lfc_bulk$coefficients
if (is.matrix(coef_mat_bulk)) {
  if ("Late_vs_Early" %in% colnames(coef_mat_bulk)) {
    fc_vec_bulk <- drop(coef_mat_bulk[, "Late_vs_Early"])
  } else {
    fc_vec_bulk <- drop(coef_mat_bulk[, 1])
  }
  names(fc_vec_bulk) <- rownames(coef_mat_bulk)
} else {
  fc_vec_bulk <- as.numeric(coef_mat_bulk)
  names(fc_vec_bulk) <- rownames(combat_bulk)[seq_along(fc_vec_bulk)]
}

names(fc_vec_bulk) <- sub("_.*", "", names(fc_vec_bulk))
geneList_bulk <- sort(fc_vec_bulk, decreasing = TRUE)

BiocParallel::register(BiocParallel::SerialParam())
set.seed(123)

gsea_react_bulk <- ReactomePA::gsePathway(
  geneList     = geneList_bulk,
  organism     = "human",
  pvalueCutoff = gsea_p_cut_bulk,
  maxGSSize    = 3000,
  minGSSize    = 10,
  eps          = 0
)

ridgeplot2_bulk <- function(gsea_obj_bulk,
                            showCategory_bulk = 20,
                            pathway_text_size_bulk = 10,
                            x_breaks_bulk = seq(-2, 2, 1)) {
  enrichplot::ridgeplot(
    gsea_obj_bulk,
    showCategory = showCategory_bulk,
    fill         = "p.adjust"
  ) +
    scale_fill_gradient(low = "#b2182b", high = "#2166ac", name = "FDR") +
    scale_x_continuous(
      name   = "Preranked log₂ fold-change metric",
      breaks = x_breaks_bulk,
      labels = sprintf("%.0f", x_breaks_bulk)
    ) +
    theme(
      axis.text.y    = element_text(size = pathway_text_size_bulk),
      axis.title.y   = element_blank(),
      legend.position= "right",
      plot.title     = element_text(hjust = 0.5)
    )
}

p_reactome_gsea_ridge_bulk <- NULL
if (nrow(as.data.frame(gsea_react_bulk)) > 0) {
  p_reactome_gsea_ridge_bulk <- ridgeplot2_bulk(
    gsea_obj_bulk          = gsea_react_bulk,
    showCategory_bulk      = gsea_show_top_bulk,
    pathway_text_size_bulk = ridge_label_sz_bulk,
    x_breaks_bulk          = gsea_x_breaks_bulk
  ) +
    ggtitle(sprintf(
      "FDR Top %d Pathway: Reactome Enrichment (|log2FC| >= %.2f)",
      gsea_show_top_bulk,
      lfc_cutoff_bulk
    )) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 14, margin = margin(b = 6)),
      axis.text.x  = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      legend.text  = element_text(size = 9),
      legend.title = element_text(size = 10)
    )
}

gsea_hits_bulk <- as.data.frame(gsea_react_bulk)
gsea_hits_bulk <- dplyr::filter(gsea_hits_bulk, p.adjust < gsea_p_cut_bulk)
gsea_hits_bulk <- dplyr::select(gsea_hits_bulk, Description, NES, p.adjust)
gsea_hits_bulk <- dplyr::arrange(gsea_hits_bulk, p.adjust)

sig3_bulk <- function(x, sig_bulk = 3) {
  vapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    if (v == 0)  return("0")
    dec_bulk <- max(sig_bulk - 1 - floor(log10(abs(v))), 0)
    formatC(v, digits = dec_bulk, format = "f", drop0trailing = FALSE)
  }, FUN.VALUE = character(1))
}

gsea_hits_print_bulk <- gsea_hits_bulk %>%
  dplyr::mutate(
    NES      = signif(NES, 3),
    p.adjust = sig3_bulk(p.adjust)
  )

if (isTRUE(print_gsea_hits_bulk)) {
  if (nrow(gsea_hits_print_bulk) == 0L) {
    message("Reactome GSEA significant pathways: 0")
  } else {
    message("Reactome GSEA significant pathways: ", nrow(gsea_hits_print_bulk))
    print(gsea_hits_print_bulk, row.names = FALSE)
  }
}

# Reactome GSEA significant pathways: 39
# Description                                                              NES    p.adjust
# Phosphorylation of CD3 and TCR zeta chains                              2.60 0.000000150
# Co-inhibition by PD-1                                                   2.56 0.000000919
# Translocation of ZAP-70 to Immunological synapse                        2.52  0.00000119
# Aerobic respiration and respiratory electron transport                 -2.07  0.00000119
# Generation of second messenger molecules                                2.31    0.000529
# Extracellular matrix organization                                       1.74    0.000603
# Mitochondrial protein degradation                                      -2.05    0.000823
# Degradation of the extracellular matrix                                 1.92     0.00304
# Activation of Matrix Metalloproteinases                                 2.19     0.00366
# Regulation of T cell activation by CD28 family                          2.06     0.00366
# Myogenesis                                                             -2.09     0.00534
# Mitochondrial Fatty Acid Beta-Oxidation                                -2.11     0.00551
# Interferon gamma signaling                                              1.91     0.00691
# Innate Immune System                                                    1.36     0.00754
# Keratan sulfate degradation                                             2.15      0.0112
# Assembly of collagen fibrils and other multimeric structures            1.93      0.0145
# Transcriptional regulation of white adipocyte differentiation           1.87      0.0145
# Protein localization                                                   -1.74      0.0150
# Downstream TCR signaling                                                1.92      0.0150
# Dissolution of Fibrin Clot                                              2.08      0.0178
# Adipogenesis                                                            1.83      0.0178
# Collagen degradation                                                    1.93      0.0220
# Epigenetic regulation of gene expression by MLL3 and MLL4 complexes     1.77      0.0253
# MLL4 and MLL3 complexes regulate expression of PPARG target genes       1.77      0.0253
# Epigenetic regulation of adipogenesis genes by MLL3 and MLL4 complexes  1.77      0.0253
# Neutrophil degranulation                                                1.45      0.0253
# Scavenging of heme from plasma                                          2.02      0.0321
# NR1H2 and NR1H3-mediated signaling                                      1.98      0.0321
# Peroxisomal protein import                                             -1.82      0.0321
# Collagen formation                                                      1.77      0.0321
# Respiratory electron transport                                         -1.71      0.0321
# Immune System                                                           1.22      0.0321
# Cell Cycle, Mitotic                                                    -1.42      0.0390
# mitochondrial fatty acid beta-oxidation of saturated fatty acids       -1.86      0.0436
# MHC class II antigen presentation                                       1.70      0.0436
# Mitochondrial biogenesis                                               -1.69      0.0436
# Pyruvate metabolism                                                    -1.81      0.0459
# Cell Cycle                                                             -1.36      0.0475
# TCR signaling                                                           1.69      0.0488


#  5. Plots----------------------------------------------------------------------
{
  message("===== Bulk pathway analysis: display =====")
  
  show_selected_plot_bulk <- function(plot_name_bulk) {
    valid_bulk <- c(
      "none",
      "reactome_dot", "reactome_bar",
      "go_dot",       "go_bar",
      "kegg_dot",     "kegg_bar",
      "reactome_gsea_ridge"
    )
    
    if (!plot_name_bulk %in% valid_bulk) {
      stop(
        "plot_to_show_bulk must be one of: ",
        paste(valid_bulk, collapse = ", ")
      )
    }
    
    if (identical(plot_name_bulk, "none")) {
      message("No plot displayed. Change plot_to_show_bulk to display one figure.")
      return(invisible(NULL))
    }
    
    plot_obj_bulk <- switch(
      plot_name_bulk,
      reactome_dot        = p_reactome_dot_bulk,
      reactome_bar        = p_reactome_bar_bulk,
      go_dot              = p_go_dot_bulk,
      go_bar              = p_go_bar_bulk,
      kegg_dot            = p_kegg_dot_bulk,
      kegg_bar            = p_kegg_bar_bulk,
      reactome_gsea_ridge = p_reactome_gsea_ridge_bulk
    )
    
    if (is.null(plot_obj_bulk)) {
      message("Selected plot is NULL (likely no significant pathway for that result).")
      return(invisible(NULL))
    }
    
    print(plot_obj_bulk)
    invisible(plot_obj_bulk)
  }
}

## Display plot
## Choices:
##   "none"
##   "reactome_dot", "reactome_bar"
##   "go_dot",       "go_bar"
##   "kegg_dot",     "kegg_bar"
##   "reactome_gsea_ridge"
## In this script, the author choices 3 dotplots for ORA and GSEA ridgeplot.


# Reactome dotplot
{
  reactome_plot_to_show_bulk   <- "reactome_dot"
  show_selected_plot_bulk(reactome_plot_to_show_bulk)
}


# GO-BP dotplot
{
  go_plot_to_show_bulk <- "go_dot"
  show_selected_plot_bulk(go_plot_to_show_bulk)
}


# KEGG dotplot
{
  kegg_plot_to_show_bulk <- "kegg_dot"
  show_selected_plot_bulk(kegg_plot_to_show_bulk)
}


# GSEA ridgeplot
{
  gsea_plot_to_show_bulk   <- "reactome_gsea_ridge"
  show_selected_plot_bulk(gsea_plot_to_show_bulk)
}


###─────────────────────────────────────── C: Sensitivity analysis-----------------------------------------------------

#  0. Parameters / packages / helper functions----------------------------------

proj_root_bulk   <- "C:/DMD_project"
cache_dir_bulk   <- file.path(proj_root_bulk, "cache_bulk_revision")
res_dir_sense_bulk <- file.path(proj_root_bulk, "results_sensitivity_bulk")
dir.create(res_dir_sense_bulk, showWarnings = FALSE, recursive = TRUE)

age_cutoff_sense_bulk    <- 6
lfc_cutoff_sense_bulk    <- 0.15
fdr_cutoff_sense_bulk    <- 0.05
volcano_lfc_sense_bulk   <- 0.30
heatmap_cap_sense_bulk   <- 217
show_top_gsea_sense_bulk <- 20
x_breaks_gsea_sense_bulk <- seq(-2, 2, 1)
seed_sense_bulk          <- 123
panel_cutoffs_sense_bulk <- c(4, 6, 7, 8)

suppressPackageStartupMessages({
  library(limma)
  library(sva)
  library(matrixStats)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(ReactomePA)
  library(BiocParallel)
  library(enrichplot)
  library(stringr)
})

if ("package:conflicted" %in% search()) {
  conflicted::conflict_prefer("select",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("filter",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("mutate",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("arrange",   "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("slice_min", "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("slice_head", "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("union",     "base",  quiet = TRUE)
  conflicted::conflict_prefer("intersect", "base",  quiet = TRUE)
  conflicted::conflict_prefer("setdiff",   "base",  quiet = TRUE)
  conflicted::conflict_prefer("setequal",  "base",  quiet = TRUE)
  conflicted::conflict_prefer("unname",    "base",  quiet = TRUE)
}

plotPCA_micro_sense_bulk <- function(expr_bulk, group_bulk, cohort_bulk,
                                     main_bulk = "PCA after ComBat") {
  pc_bulk <- prcomp(t(expr_bulk))
  col_vec_bulk <- c(Early = "#2C7BB6", Late = "#D7191C")[group_bulk]
  pch_vec_bulk <- c(C1 = 17, C2 = 7, C3 = 3, C4 = 5)[cohort_bulk]
  
  plot(
    pc_bulk$x[, 1], pc_bulk$x[, 2],
    col  = col_vec_bulk,
    pch  = pch_vec_bulk,
    xlab = sprintf("PC1 (%.1f%%)", 100 * pc_bulk$sdev[1]^2 / sum(pc_bulk$sdev^2)),
    ylab = sprintf("PC2 (%.1f%%)", 100 * pc_bulk$sdev[2]^2 / sum(pc_bulk$sdev^2)),
    main = main_bulk
  )
  
  legend("topleft",
         legend = levels(group_bulk),
         col    = c("#2C7BB6", "#D7191C"),
         pch    = 16,
         title  = "Group",
         bty    = "n")
  legend("topright",
         legend = levels(cohort_bulk),
         pch    = c(17, 7, 3, 5),
         title  = "Cohort",
         bty    = "n")
}

plotRLE_micro_sense_bulk <- function(expr_bulk, cohort_bulk,
                                     main_bulk = "RLE after ComBat") {
  med_bulk <- matrixStats::rowMedians(expr_bulk, na.rm = TRUE)
  rle_bulk <- sweep(expr_bulk, 1, med_bulk, "-")
  
  boxplot(
    as.data.frame(rle_bulk),
    las     = 2,
    outline = FALSE,
    col     = c(C1 = "black", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")[cohort_bulk],
    ylab    = "Relative log expression",
    main    = main_bulk
  )
  
  legend("topright",
         legend = levels(cohort_bulk),
         fill   = c(C1 = "black", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")[levels(cohort_bulk)],
         title  = "Cohort",
         bty    = "n")
}

safe_div_sense_bulk <- function(num_bulk, den_bulk) {
  if (is.na(den_bulk) || den_bulk == 0) return(NA_real_)
  num_bulk / den_bulk
}

deg_overlap_metrics_sense_bulk <- function(ref_deg_bulk, ref_label_bulk,
                                           test_deg_bulk, test_label_bulk) {
  ref_set_bulk  <- rownames(ref_deg_bulk)
  test_set_bulk <- rownames(test_deg_bulk)
  
  inter_bulk <- length(base::intersect(ref_set_bulk, test_set_bulk))
  union_bulk <- length(base::union(ref_set_bulk, test_set_bulk))
  common_bulk <- base::intersect(ref_set_bulk, test_set_bulk)
  
  dir_cons_bulk <- if (length(common_bulk) == 0) {
    NA_real_
  } else {
    mean(sign(ref_deg_bulk[common_bulk, "logFC"]) ==
           sign(test_deg_bulk[common_bulk, "logFC"]))
  }
  
  data.frame(
    ref       = ref_label_bulk,
    test      = test_label_bulk,
    n_ref     = length(ref_set_bulk),
    n_test    = length(test_set_bulk),
    n_inter   = inter_bulk,
    jaccard   = round(safe_div_sense_bulk(inter_bulk, union_bulk), 3),
    precision = round(safe_div_sense_bulk(inter_bulk, length(test_set_bulk)), 3),
    recall    = round(safe_div_sense_bulk(inter_bulk, length(ref_set_bulk)), 3),
    dir_cons  = round(dir_cons_bulk, 3),
    stringsAsFactors = FALSE
  )
}

sig3_sense_bulk <- function(x, sig = 3) {
  vapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    if (v == 0) return("0")
    dec_bulk <- max(sig - 1 - floor(log10(abs(v))), 0)
    formatC(v, digits = dec_bulk, format = "f", drop0trailing = FALSE)
  }, FUN.VALUE = character(1))
}

pub_dotplot_sense_bulk <- function(eres_bulk, n_bulk = 15, fdr_cut_bulk = 0.05,
                                   title_bulk = "", wrap_len_bulk = 35) {
  if (is.null(eres_bulk)) {
    message("＜", title_bulk, "＞ analysis object is NULL ")
    return(invisible(NULL))
  }
  df_bulk <- as.data.frame(eres_bulk)
  if (!nrow(df_bulk)) {
    message("＜", title_bulk, "＞ no significant pathway")
    return(invisible(NULL))
  }
  df_bulk <- df_bulk[df_bulk$p.adjust <= fdr_cut_bulk, , drop = FALSE]
  if (!nrow(df_bulk)) {
    message("＜", title_bulk, "＞ no significant pathway (FDR < ", fdr_cut_bulk, ")")
    return(invisible(NULL))
  }
  df_bulk <- df_bulk[order(df_bulk$p.adjust), , drop = FALSE]
  df_bulk <- df_bulk[seq_len(min(n_bulk, nrow(df_bulk))), , drop = FALSE]
  df_bulk$GeneRatio <- sapply(strsplit(df_bulk$GeneRatio, "/"), function(x) {
    as.numeric(x[1]) / as.numeric(x[2])
  })
  df_bulk$Description <- factor(
    stringr::str_wrap(df_bulk$Description, wrap_len_bulk),
    levels = rev(stringr::str_wrap(df_bulk$Description, wrap_len_bulk))
  )
  
  ggplot(df_bulk, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = Count, colour = p.adjust)) +
    scale_colour_gradient(low = "#b2182b", high = "#2166ac", name = "FDR") +
    scale_size(range = c(3, 8)) +
    labs(title = title_bulk, x = "Gene ratio", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9)
    )
}

ridgeplot2_sense_bulk <- function(gsea_obj_bulk,
                                  showCategory_bulk = 20,
                                  pathway_text_size_bulk = 10,
                                  x_breaks_bulk = seq(-2, 2, 1)) {
  if (is.null(gsea_obj_bulk) || !nrow(as.data.frame(gsea_obj_bulk))) {
    message("Reactome GSEA: no significant pathway")
    return(invisible(NULL))
  }
  enrichplot::ridgeplot(
    gsea_obj_bulk,
    showCategory = showCategory_bulk,
    fill         = "p.adjust"
  ) +
    scale_fill_gradient(low = "#b2182b", high = "#2166ac", name = "FDR") +
    scale_x_continuous(
      name   = "Preranked log₂ fold-change metric",
      breaks = x_breaks_bulk,
      labels = sprintf("%.0f", x_breaks_bulk)
    ) +
    theme(
      axis.text.y  = element_text(size = pathway_text_size_bulk),
      axis.title.y = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )
}

safe_enrichPathway_sense_bulk <- function(gene_bulk, universe_bulk) {
  tryCatch(
    ReactomePA::enrichPathway(
      gene         = gene_bulk,
      universe     = universe_bulk,
      organism     = "human",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.10
    ),
    error = function(e) NULL
  )
}

safe_enrichGO_sense_bulk <- function(gene_bulk, universe_bulk) {
  out_bulk <- tryCatch(
    clusterProfiler::enrichGO(
      gene          = gene_bulk,
      universe      = universe_bulk,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05
    ),
    error = function(e) NULL
  )
  if (is.null(out_bulk)) return(NULL)
  if (!nrow(as.data.frame(out_bulk))) return(out_bulk)
  tryCatch(
    clusterProfiler::simplify(out_bulk, cutoff = 0.10, by = "p.adjust", select_fun = min),
    error = function(e) out_bulk
  )
}

safe_enrichKEGG_sense_bulk <- function(gene_bulk, universe_bulk) {
  out_bulk <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene         = gene_bulk,
      universe     = universe_bulk,
      organism     = "hsa",
      keyType      = "ncbi-geneid",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.10
    ),
    error = function(e) NULL
  )
  if (is.null(out_bulk)) return(NULL)
  df_bulk <- as.data.frame(out_bulk)
  if (!nrow(df_bulk)) return(out_bulk)
  df_bulk <- dplyr::filter(df_bulk, !grepl("^hsa05|disease", ID))
  out_bulk@result <- df_bulk
  out_bulk
}

safe_gsePathway_sense_bulk <- function(geneList_bulk, seed_bulk = 123) {
  BiocParallel::register(BiocParallel::SerialParam())
  set.seed(seed_bulk)
  tryCatch(
    ReactomePA::gsePathway(
      geneList     = geneList_bulk,
      organism     = "human",
      pvalueCutoff = 0.05,
      maxGSSize    = 3000,
      minGSSize    = 10,
      eps          = 0
    ),
    error = function(e) NULL
  )
}

count_sig_terms_sense_bulk <- function(eres_bulk, fdr_bulk = 0.05) {
  if (is.null(eres_bulk)) return(0L)
  df_bulk <- tryCatch(as.data.frame(eres_bulk), error = function(e) data.frame())
  if (!nrow(df_bulk) || !("p.adjust" %in% colnames(df_bulk))) return(0L)
  sum(!is.na(df_bulk$p.adjust) & df_bulk$p.adjust < fdr_bulk)
}


#  1.  Rerun ComBat in "6-year-old version"-----------------------------------------------

stage1_bulk <- readRDS(file.path(cache_dir_bulk, "samples_meta_bulk.rds"))
list2env(stage1_bulk, .GlobalEnv)

stage5_bulk <- readRDS(file.path(cache_dir_bulk, "normalized_imputed_bulk.rds"))
list2env(stage5_bulk, .GlobalEnv)

# Note: "stage6" here is the cache stage number carried over from the main bulk pipeline.
#       It does not indicate a 6-year cutoff.
bulk_main_stage_cache <- readRDS(file.path(cache_dir_bulk, "stage6_combat_only_bulk.rds"))
list2env(bulk_main_stage_cache, .GlobalEnv)

meta_sense_bulk <- meta_aligned_bulk[match(colnames(expr_imp_bulk), meta_aligned_bulk$PID), , drop = FALSE]
stopifnot(all(meta_sense_bulk$PID == colnames(expr_imp_bulk)))

Group_sense_bulk <- factor(
  ifelse(meta_sense_bulk$Age < age_cutoff_sense_bulk, "Early", "Late"),
  levels = c("Early", "Late")
)

Cohort_sense_bulk <- factor(
  meta_sense_bulk$Cohort,
  levels = c("C1", "C2", "C3", "C4")
)

design_sense_bulk <- model.matrix(~ 0 + Group_sense_bulk)
colnames(design_sense_bulk) <- levels(Group_sense_bulk)

combat_sense_bulk <- sva::ComBat(
  dat       = expr_imp_bulk,
  batch     = Cohort_sense_bulk,
  mod       = model.matrix(~ Group_sense_bulk),
  par.prior = TRUE,
  mean.only = FALSE
)

pc_sense_bulk  <- prcomp(t(combat_sense_bulk))$x
med_sense_bulk <- apply(combat_sense_bulk, 1, median, na.rm = TRUE)

cat(sprintf("%dyrs cutoff  Group: Early=%d / Late=%d\n",
            age_cutoff_sense_bulk,
            sum(Group_sense_bulk == "Early"),
            sum(Group_sense_bulk == "Late")))
# Group: Early=19 / Late=13

#  2.  Run the limma vooma-treat pipeline with the "_sense_bulk" suffix--------------

designGLM_sense_bulk <- model.matrix(~ 0 + Group_sense_bulk)
colnames(designGLM_sense_bulk) <- levels(Group_sense_bulk)

contrast_sense_bulk <- matrix(
  c(-1, 1),
  ncol = 1,
  dimnames = list(colnames(designGLM_sense_bulk), "Late_vs_Early_sense_bulk")
)

v_sense_bulk <- limma::vooma(combat_sense_bulk, designGLM_sense_bulk, plot = FALSE)

fit_sense_bulk <- v_sense_bulk |>
  limma::lmFit(designGLM_sense_bulk) |>
  limma::contrasts.fit(contrast_sense_bulk) |>
  limma::eBayes(trend = TRUE) |>
  limma::treat(lfc = lfc_cutoff_sense_bulk)

deg_sense_bulk <- limma::topTreat(
  fit_sense_bulk,
  coef    = "Late_vs_Early_sense_bulk",
  p.value = fdr_cutoff_sense_bulk,
  number  = Inf
)

cat(sprintf("%dyrs cutoff  DEGs: %d\n", age_cutoff_sense_bulk, nrow(deg_sense_bulk)))
# DEGs: 217

#  3.  QC plot---------------------------------------------------------------------

plotRLE_micro_sense_bulk(
  expr_bulk   = combat_sense_bulk,
  cohort_bulk = Cohort_sense_bulk,
  main_bulk   = sprintf("RLE (6-y cutoff, after ComBat)")
)

plotPCA_micro_sense_bulk(
  expr_bulk   = combat_sense_bulk,
  group_bulk  = Group_sense_bulk,
  cohort_bulk = Cohort_sense_bulk,
  main_bulk   = sprintf("PCA (6-y cutoff, after ComBat)")
)


#  4.  FDR significant heatmap--------------------------------------------------

{
  deg_df_sense_bulk <- as.data.frame(deg_sense_bulk)
  fdr_candidates_sense_bulk <- c("FDR", "adj.P.Val", "padj", "qvalue", "p.adj", "FDR.BH")
  fdr_col_sense_bulk <- base::intersect(fdr_candidates_sense_bulk, colnames(deg_df_sense_bulk))
  if (length(fdr_col_sense_bulk) == 0) stop("The column corresponding to FDR was not found.")
  fdr_col_sense_bulk <- fdr_col_sense_bulk[1]
  
  deg_df_sorted_sense_bulk <- deg_df_sense_bulk[order(deg_df_sense_bulk[[fdr_col_sense_bulk]]), , drop = FALSE]
  top_n_sense_bulk <- min(heatmap_cap_sense_bulk, nrow(deg_df_sorted_sense_bulk))
  top_id_sense_bulk <- head(rownames(deg_df_sorted_sense_bulk), top_n_sense_bulk)
  
  expr_top_sense_bulk <- combat_sense_bulk[top_id_sense_bulk, , drop = FALSE]
  expr_top_sense_bulk <- expr_top_sense_bulk[matrixStats::rowSds(expr_top_sense_bulk) > 0.15, , drop = FALSE]
  n_gene_sense_bulk <- nrow(expr_top_sense_bulk)
  
  expr_z_sense_bulk <- t(scale(t(expr_top_sense_bulk)))
  expr_z_sense_bulk <- pmin(pmax(expr_z_sense_bulk, -2.5), 2.5)
  
  sym_sense_bulk <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = rownames(expr_z_sense_bulk),
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )
  keep_sense_bulk <- !is.na(sym_sense_bulk) & sym_sense_bulk != ""
  expr_z_sense_bulk <- expr_z_sense_bulk[keep_sense_bulk, , drop = FALSE]
  rownames(expr_z_sense_bulk) <- sym_sense_bulk[keep_sense_bulk]
  
  ann_col_sense_bulk <- data.frame(
    Group  = Group_sense_bulk,
    Cohort = factor(Cohort_sense_bulk, levels = c("C1", "C2", "C3", "C4"))
  )
  rownames(ann_col_sense_bulk) <- colnames(expr_z_sense_bulk)
  
  ann_cols_sense_bulk <- list(
    Group  = c(Early = "#2C7BB6", Late = "#D7191C"),
    Cohort = c(C1 = "black", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")
  )
  
  pheatmap::pheatmap(
    expr_z_sense_bulk,
    scale = "none",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    breaks = seq(-2.5, 2.5, length = 101),
    annotation_col    = ann_col_sense_bulk,
    annotation_colors = ann_cols_sense_bulk,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    show_rownames = FALSE,
    border_color  = FALSE,
    fontsize_row  = 5,
    main = sprintf("DEG heatmap | FDR significant  ->  %d genes x %d samples (6yr)",
                   n_gene_sense_bulk, ncol(expr_z_sense_bulk))
  )
}


#  5.  Comparing 6-year-old cuts based on the primary analysis at age 5----------------------------------------

{
  deg_metrics6_sense_bulk <- deg_overlap_metrics_sense_bulk(
    ref_deg_bulk    = deg_bulk,
    ref_label_bulk  = "5yr",
    test_deg_bulk   = deg_sense_bulk,
    test_label_bulk = "6yr"
  )
  print(deg_metrics6_sense_bulk)
}
#   ref  test  n_ref  n_test  n_inter   jaccard   precision   recall   dir_cons
# 1 5yr   6yr    358     217      156     0.372       0.719    0.436          1

#  6.  Common DEGs heatmap (Main vs 6yr)----------------------------------------

{
  common_ent_sense_bulk <- base::intersect(rownames(deg_bulk), rownames(deg_sense_bulk))
  
  idx_common_sense_bulk <- sub("_.*", "", rownames(combat_bulk)) %in% common_ent_sense_bulk
  mat_common_sense_bulk <- combat_bulk[idx_common_sense_bulk, , drop = FALSE]
  matZ_common_sense_bulk <- t(scale(t(mat_common_sense_bulk)))
  
  entrez_common_sense_bulk <- sub("_.*", "", rownames(matZ_common_sense_bulk))
  sym_common_sense_bulk <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = entrez_common_sense_bulk,
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )
  
  keep_common_sense_bulk <- !is.na(sym_common_sense_bulk) & sym_common_sense_bulk != ""
  matZ_common_sense_bulk <- matZ_common_sense_bulk[keep_common_sense_bulk, , drop = FALSE]
  rownames(matZ_common_sense_bulk) <- make.unique(sym_common_sense_bulk[keep_common_sense_bulk])
  
  ann_col_common_sense_bulk <- data.frame(
    Group  = Group_imp_bulk,
    Cohort = factor(Cohort_imp_bulk, levels = c("C1", "C2", "C3", "C4"))
  )
  rownames(ann_col_common_sense_bulk) <- colnames(matZ_common_sense_bulk)
  
  pheatmap::pheatmap(
    matZ_common_sense_bulk,
    annotation_col    = ann_col_common_sense_bulk,
    annotation_colors = list(
      Group  = c(Early = "#2C7BB6", Late = "#D7191C"),
      Cohort = c(C1 = "#000000", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")
    ),
    color         = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    cluster_rows  = TRUE,
    cluster_cols  = TRUE,
    show_rownames = TRUE,
    fontsize_row  = 5,
    border_color  = "grey60",
    main = sprintf("Common %d DEGs (6yr vs 5yr)", nrow(matZ_common_sense_bulk))
  )
}


#  7.  Volcano plot — 6-year cutoff analysis---------------------------------------

{
  tbl_all_sense_bulk <- limma::topTreat(
    fit_sense_bulk,
    coef    = "Late_vs_Early_sense_bulk",
    p.value = 1,
    number  = Inf
  )
  
  tbl_plot_sense_bulk <- tbl_all_sense_bulk |>
    tibble::rownames_to_column("EntrezID") |>
    dplyr::mutate(
      GeneSymbol = AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys      = EntrezID,
        column    = "SYMBOL",
        keytype   = "ENTREZID",
        multiVals = "first"
      ),
      GeneSymbol = ifelse(is.na(GeneSymbol), EntrezID, GeneSymbol),
      signif_sense_bulk = dplyr::case_when(
        adj.P.Val < fdr_cutoff_sense_bulk & logFC >=  volcano_lfc_sense_bulk ~ "Up",
        adj.P.Val < fdr_cutoff_sense_bulk & logFC <= -volcano_lfc_sense_bulk ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  pattern_exclude_sense_bulk <- "^(LOC|LINC|CH[0-9]|AC[0-9]|AL[0-9]|CT[0-9]|RP[0-9])"
  
  top_labs_sense_bulk <- dplyr::bind_rows(
    tbl_plot_sense_bulk |>
      dplyr::filter(signif_sense_bulk == "Up", !grepl(pattern_exclude_sense_bulk, GeneSymbol)) |>
      dplyr::slice_min(order_by = adj.P.Val, n = 5, with_ties = FALSE),
    tbl_plot_sense_bulk |>
      dplyr::filter(signif_sense_bulk == "Down", !grepl(pattern_exclude_sense_bulk, GeneSymbol)) |>
      dplyr::slice_min(order_by = adj.P.Val, n = 5, with_ties = FALSE)
  )
  
  print(
    ggplot(tbl_plot_sense_bulk, aes(logFC, -log10(adj.P.Val))) +
      geom_point(aes(colour = signif_sense_bulk), size = 1.6, alpha = 0.8) +
      scale_colour_manual(values = c(Up = "#D7191C", Down = "#2C7BB6", NS = "grey70")) +
      geom_vline(xintercept = c(-volcano_lfc_sense_bulk, volcano_lfc_sense_bulk),
                 linetype = "dashed", colour = "grey40") +
      geom_hline(yintercept = -log10(fdr_cutoff_sense_bulk),
                 linetype = "dashed", colour = "grey40") +
      ggrepel::geom_text_repel(data = top_labs_sense_bulk,
                               aes(label = GeneSymbol),
                               size = 3,
                               max.overlaps = Inf) +
      labs(title = sprintf("Volcano plot | %d DEGs (|logFC| >= %.2f, FDR <= %.2f, 6-yr cutoff)",
                           sum(tbl_plot_sense_bulk$signif_sense_bulk != "NS"),
                           volcano_lfc_sense_bulk,
                           fdr_cutoff_sense_bulk),
           x = "log2 fold-change",
           y = expression(-log[10]~FDR)) +
      theme_bw(base_size = 12) +
      theme(legend.title = element_blank())
  )
}


#  8.  List of common DEG gene names (5yr ∩ 6yr)-------------------------------------

{
  common_id_sense_bulk <- base::intersect(rownames(deg_bulk), rownames(deg_sense_bulk))
  cat(sprintf("Common DEGs (5yr ∩ 6yr): %d\n", length(common_id_sense_bulk)))
  
  tbl_common_sense_bulk <- data.frame(
    SYMBOL = AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = common_id_sense_bulk,
      column    = "SYMBOL",
      keytype   = "ENTREZID",
      multiVals = "first"
    ),
    logFC     = deg_bulk[common_id_sense_bulk, "logFC"],
    adj.P.Val = deg_bulk[common_id_sense_bulk, "adj.P.Val"],
    row.names = NULL,
    check.names = FALSE
  )
  
  tbl_common_sense_bulk <- tbl_common_sense_bulk[order(tbl_common_sense_bulk$adj.P.Val), , drop = FALSE]
  print(tbl_common_sense_bulk, digits = 3, row.names = FALSE)
  
  utils::write.table(
    tbl_common_sense_bulk,
    file      = file.path(res_dir_sense_bulk, "common_deg_5yr_vs_6yr_sense_bulk.tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
}
# Common DEGs (5yr ∩ 6yr): 156

#  9.  KEGG, Reactome, GO-BP (6-yr cutoff)--------------------------------------

{
  bg_ids_sense_bulk  <- unique(sub("_.*", "", rownames(combat_sense_bulk)))
  deg_ids_sense_bulk <- rownames(deg_sense_bulk)
  
  react_res_sense_bulk <- safe_enrichPathway_sense_bulk(deg_ids_sense_bulk, bg_ids_sense_bulk)
  go_res_sense_bulk    <- safe_enrichGO_sense_bulk(deg_ids_sense_bulk, bg_ids_sense_bulk)
  kegg_clean_sense_bulk <- safe_enrichKEGG_sense_bulk(deg_ids_sense_bulk, bg_ids_sense_bulk)
}

print(pub_dotplot_sense_bulk(
  react_res_sense_bulk,
  n_bulk       = 20,
  fdr_cut_bulk = 0.05,
  title_bulk   = "Reactome enrichment (6-yr)"
))

print(pub_dotplot_sense_bulk(
  go_res_sense_bulk,
  n_bulk       = 20,
  fdr_cut_bulk = 0.05,
  title_bulk   = "GO-BP enrichment (6-yr)"
))

print(pub_dotplot_sense_bulk(
  kegg_clean_sense_bulk,
  n_bulk       = 20,
  fdr_cut_bulk = 0.05,
  title_bulk   = "KEGG enrichment (6-yr)"
))


# 10.  GSEA (Reactome) (6-yr cutoff)-------------------------------------------

{
  fit_lfc_sense_bulk <- v_sense_bulk |>
    limma::lmFit(designGLM_sense_bulk) |>
    limma::contrasts.fit(contrast_sense_bulk) |>
    limma::eBayes(trend = TRUE) |>
    limma::treat(lfc = lfc_cutoff_sense_bulk)
  
  fc_vec_sense_bulk <- drop(fit_lfc_sense_bulk$coefficients[, "Late_vs_Early_sense_bulk"])
  names(fc_vec_sense_bulk) <- sub("_.*", "", rownames(fit_lfc_sense_bulk$coefficients))
  geneList_sense_bulk <- sort(fc_vec_sense_bulk, decreasing = TRUE)
  
  gsea_react_sense_bulk <- safe_gsePathway_sense_bulk(geneList_sense_bulk, seed_bulk = seed_sense_bulk)
}

{
  p_ridge_sense_bulk <- ridgeplot2_sense_bulk(
    gsea_obj_bulk          = gsea_react_sense_bulk,
    showCategory_bulk      = show_top_gsea_sense_bulk,
    pathway_text_size_bulk = 10,
    x_breaks_bulk          = x_breaks_gsea_sense_bulk
  )
  
  if (!is.null(p_ridge_sense_bulk)) {
    print(
      p_ridge_sense_bulk +
        ggtitle(sprintf("Reactome enrichment (|log2FC| > %.2f, 6-yr)", lfc_cutoff_sense_bulk)) +
        theme(
          plot.title   = element_text(hjust = 0.5, size = 14, margin = margin(b = 6)),
          axis.text.x  = element_text(size = 10),
          axis.title.x = element_text(size = 11),
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 10)
        )
    )
  }
}

{
  tbl_short_sense_bulk <- if (is.null(gsea_react_sense_bulk)) {
    data.frame()
  } else {
    as.data.frame(gsea_react_sense_bulk) |>
      dplyr::filter(p.adjust < 0.05) |>
      dplyr::select(Description, NES, p.adjust) |>
      dplyr::arrange(p.adjust)
  }
  
  tbl_print_sense_bulk <- tbl_short_sense_bulk |>
    dplyr::mutate(p.adjust = sig3_sense_bulk(p.adjust))
  
  print(tbl_print_sense_bulk, row.names = FALSE)
  
  utils::write.table(
    tbl_print_sense_bulk,
    file      = file.path(res_dir_sense_bulk, "reactome_gsea_hits_6yr_sense_bulk.tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
}


# 11.  Main vs 6yr : Reactome GSEA concordance---------------------------------

{
  fc_vec_main_sense_bulk <- drop(fit_treat_bulk$coefficients[, "Late_vs_Early_bulk"])
  names(fc_vec_main_sense_bulk) <- sub("_.*", "", rownames(fit_treat_bulk$coefficients))
  geneList_bulk <- sort(fc_vec_main_sense_bulk, decreasing = TRUE)
  
  gsea_react_bulk <- safe_gsePathway_sense_bulk(geneList_bulk, seed_bulk = seed_sense_bulk)
  
  tbl_main_bulk_gsea <- if (is.null(gsea_react_bulk)) {
    data.frame(ID = character(), NES = numeric(), p.adjust = numeric())
  } else {
    as.data.frame(gsea_react_bulk) |>
      dplyr::select(ID, NES, p.adjust)
  }
  
  tbl_sense_bulk_gsea <- if (is.null(gsea_react_sense_bulk)) {
    data.frame(ID = character(), NES_sense = numeric(), p.sense = numeric())
  } else {
    as.data.frame(gsea_react_sense_bulk) |>
      dplyr::select(ID, NES, p.adjust) |>
      dplyr::rename(NES_sense = NES, p.sense = p.adjust)
  }
  
  sig_main_bulk_gsea  <- dplyr::filter(tbl_main_bulk_gsea, p.adjust < 0.05)
  sig_sense_bulk_gsea <- dplyr::filter(tbl_sense_bulk_gsea, p.sense < 0.05)
  
  overlap_sense_bulk <- dplyr::inner_join(sig_main_bulk_gsea, sig_sense_bulk_gsea, by = "ID")
  denom_sense_bulk <- nrow(sig_main_bulk_gsea) + nrow(sig_sense_bulk_gsea) - nrow(overlap_sense_bulk)
  jaccard_sense_bulk <- safe_div_sense_bulk(nrow(overlap_sense_bulk), denom_sense_bulk)
  cor_NES_sense_bulk <- if (nrow(overlap_sense_bulk) >= 2) {
    cor(overlap_sense_bulk$NES, overlap_sense_bulk$NES_sense, method = "spearman")
  } else {
    NA_real_
  }
  
  cat(sprintf("shared pathway %d  | Jaccard %.2f | NES correlation rho = %.3f\n",
              nrow(overlap_sense_bulk), jaccard_sense_bulk, cor_NES_sense_bulk))
  
  if (nrow(overlap_sense_bulk) > 0) {
    print(
      ggplot(overlap_sense_bulk, aes(NES, NES_sense)) +
        geom_point(size = 2, alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, colour = "red") +
        labs(x = "NES (Main)", y = "NES (6-yr)", title = "Pathway NES concordance") +
        theme_bw()
    )
  }
  
  if (length(geneList_bulk) > 0 && length(geneList_sense_bulk) > 0) {
    comp_sense_bulk <- clusterProfiler::compareCluster(
      geneCluster  = stats::setNames(list(geneList_bulk, geneList_sense_bulk), c("Main", "6yr")),
      fun          = "gsePathway",
      pvalueCutoff = 0.05,
      organism     = "human"
    )
    
    print(
      enrichplot::dotplot(comp_sense_bulk, showCategory = 20, font.size = 6) +
        ggtitle("Reactome GSEA – Main vs 6-yr") +
        theme(
          axis.text.y = element_text(size = 6, margin = margin(r = 2)),
          plot.title  = element_text(hjust = 0.5, size = 13)
        )
    )
  }
}
# shared pathway 27  | Jaccard 0.34 | NES correlation rho = 0.947


# 12.  Exploratory cutoff panel (4 / 6 / 7 / 8 years)--------------------------

## This section is intentionally compact. However, it takes significant computational time.
## Purpose:
##   - retain information that higher cutoffs may produce many DEGs,
##     yet may lose pathway-level coherence with the 5yr/6yr result;
##   - document that pattern without turning 7yr / 8yr into full figure blocks.
## Output:
##   - one summary table printed to console
##   - one TSV written to results_sensitivity_bulk/

{
  if (!exists("gsea_react_bulk", inherits = FALSE) || is.null(gsea_react_bulk)) {
    fc_vec_main_panel_bulk <- drop(fit_treat_bulk$coefficients[, "Late_vs_Early_bulk"])
    names(fc_vec_main_panel_bulk) <- sub("_.*", "", rownames(fit_treat_bulk$coefficients))
    geneList_bulk <- sort(fc_vec_main_panel_bulk, decreasing = TRUE)
    gsea_react_bulk <- safe_gsePathway_sense_bulk(geneList_bulk, seed_bulk = seed_sense_bulk)
  }
  
  sig_main_panel_bulk <- if (is.null(gsea_react_bulk)) {
    data.frame(ID = character(), NES = numeric(), p.adjust = numeric())
  } else {
    as.data.frame(gsea_react_bulk) |>
      dplyr::filter(p.adjust < 0.05) |>
      dplyr::select(ID, NES, p.adjust)
  }
  
  run_cutoff_panel_sense_bulk <- function(cutoff_bulk) {
    Group_panel_bulk <- factor(
      ifelse(meta_sense_bulk$Age < cutoff_bulk, "Early", "Late"),
      levels = c("Early", "Late")
    )
    Cohort_panel_bulk <- factor(meta_sense_bulk$Cohort, levels = c("C1", "C2", "C3", "C4"))
    
    design_panel_bulk <- model.matrix(~ 0 + Group_panel_bulk)
    colnames(design_panel_bulk) <- levels(Group_panel_bulk)
    contrast_panel_bulk <- matrix(
      c(-1, 1),
      ncol = 1,
      dimnames = list(colnames(design_panel_bulk), "Late_vs_Early_panel_bulk")
    )
    
    combat_panel_bulk <- sva::ComBat(
      dat       = expr_imp_bulk,
      batch     = Cohort_panel_bulk,
      mod       = model.matrix(~ Group_panel_bulk),
      par.prior = TRUE,
      mean.only = FALSE
    )
    
    v_panel_bulk <- limma::vooma(combat_panel_bulk, design_panel_bulk, plot = FALSE)
    fit_panel_bulk <- v_panel_bulk |>
      limma::lmFit(design_panel_bulk) |>
      limma::contrasts.fit(contrast_panel_bulk) |>
      limma::eBayes(trend = TRUE) |>
      limma::treat(lfc = lfc_cutoff_sense_bulk)
    
    deg_panel_bulk <- limma::topTreat(
      fit_panel_bulk,
      coef    = "Late_vs_Early_panel_bulk",
      p.value = fdr_cutoff_sense_bulk,
      number  = Inf
    )
    
    deg_metrics_panel_bulk <- deg_overlap_metrics_sense_bulk(
      ref_deg_bulk    = deg_bulk,
      ref_label_bulk  = "5yr",
      test_deg_bulk   = deg_panel_bulk,
      test_label_bulk = sprintf("%dyr", cutoff_bulk)
    )
    
    bg_ids_panel_bulk  <- unique(sub("_.*", "", rownames(combat_panel_bulk)))
    deg_ids_panel_bulk <- rownames(deg_panel_bulk)
    
    react_panel_bulk <- safe_enrichPathway_sense_bulk(deg_ids_panel_bulk, bg_ids_panel_bulk)
    go_panel_bulk    <- safe_enrichGO_sense_bulk(deg_ids_panel_bulk, bg_ids_panel_bulk)
    kegg_panel_bulk  <- safe_enrichKEGG_sense_bulk(deg_ids_panel_bulk, bg_ids_panel_bulk)
    
    fc_vec_panel_bulk <- drop(fit_panel_bulk$coefficients[, "Late_vs_Early_panel_bulk"])
    names(fc_vec_panel_bulk) <- sub("_.*", "", rownames(fit_panel_bulk$coefficients))
    geneList_panel_bulk <- sort(fc_vec_panel_bulk, decreasing = TRUE)
    gsea_panel_bulk <- safe_gsePathway_sense_bulk(geneList_panel_bulk, seed_bulk = seed_sense_bulk)
    
    sig_gsea_panel_bulk <- if (is.null(gsea_panel_bulk)) {
      data.frame(ID = character(), NES_panel = numeric(), p.adjust_panel = numeric())
    } else {
      as.data.frame(gsea_panel_bulk) |>
        dplyr::filter(p.adjust < 0.05) |>
        dplyr::select(ID, NES, p.adjust) |>
        dplyr::rename(NES_panel = NES, p.adjust_panel = p.adjust)
    }
    
    overlap_gsea_panel_bulk <- dplyr::inner_join(sig_main_panel_bulk, sig_gsea_panel_bulk, by = "ID")
    denom_gsea_panel_bulk <- nrow(sig_main_panel_bulk) + nrow(sig_gsea_panel_bulk) - nrow(overlap_gsea_panel_bulk)
    rho_panel_bulk <- if (nrow(overlap_gsea_panel_bulk) >= 2) {
      cor(overlap_gsea_panel_bulk$NES, overlap_gsea_panel_bulk$NES_panel, method = "spearman")
    } else {
      NA_real_
    }
    
    note_panel_bulk <- if (cutoff_bulk == 6) {
      "Detailed sensitivity cutoff used in the main sensitivity blocks above"
    } else if (cutoff_bulk >= 7) {
      "Exploratory high-cutoff reference"
    } else {
      "Exploratory low-cutoff reference"
    }
    
    data.frame(
      cutoff_bulk                = cutoff_bulk,
      Early_bulk                 = sum(Group_panel_bulk == "Early"),
      Late_bulk                  = sum(Group_panel_bulk == "Late"),
      DEG_n_bulk                 = nrow(deg_panel_bulk),
      DEG_common_vs_5yr_bulk     = deg_metrics_panel_bulk$n_inter,
      DEG_jaccard_vs_5yr_bulk    = deg_metrics_panel_bulk$jaccard,
      DEG_precision_vs_5yr_bulk  = deg_metrics_panel_bulk$precision,
      DEG_recall_vs_5yr_bulk     = deg_metrics_panel_bulk$recall,
      DEG_dir_cons_vs_5yr_bulk   = deg_metrics_panel_bulk$dir_cons,
      Reactome_ORA_n_bulk        = count_sig_terms_sense_bulk(react_panel_bulk, fdr_bulk = 0.05),
      GO_BP_ORA_n_bulk           = count_sig_terms_sense_bulk(go_panel_bulk,    fdr_bulk = 0.05),
      KEGG_ORA_n_bulk            = count_sig_terms_sense_bulk(kegg_panel_bulk,  fdr_bulk = 0.05),
      Reactome_GSEA_n_bulk       = nrow(sig_gsea_panel_bulk),
      GSEA_shared_vs_5yr_bulk    = nrow(overlap_gsea_panel_bulk),
      GSEA_jaccard_vs_5yr_bulk   = round(safe_div_sense_bulk(nrow(overlap_gsea_panel_bulk), denom_gsea_panel_bulk), 3),
      GSEA_NES_rho_vs_5yr_bulk   = round(rho_panel_bulk, 3),
      note_bulk                  = note_panel_bulk,
      stringsAsFactors = FALSE
    )
  }
  
  panel_summary_sense_bulk <- do.call(
    rbind,
    lapply(panel_cutoffs_sense_bulk, run_cutoff_panel_sense_bulk)
  )
  
  print(panel_summary_sense_bulk, row.names = FALSE)
  
  utils::write.table(
    panel_summary_sense_bulk,
    file      = file.path(res_dir_sense_bulk, "cutoff_panel_summary_sense_bulk.tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
}

# cutoff_bulk  Early_bulk  Late_bulk  DEG_n_bulk  DEG_common_vs_5yr_bulk  DEG_jaccard_vs_5yr_bulk  DEG_precision_vs_5yr_bulk  DEG_recall_vs_5yr_bulk
#       4 yrs        14         18          16                     15                   0.042                     0.938                  0.042
#       6 yrs        19         13         217                    156                   0.372                     0.719                  0.436
#       7 yrs        21         11         171                     78                   0.173                     0.456                  0.218
#       8 yrs        25          7         431                     37                   0.049                     0.086                  0.103

# cutoff_bulk  DEG_dir_cons_vs_5yr_bulk  Reactome_ORA_n_bulk  GO_BP_ORA_n_bulk  KEGG_ORA_n_bulk  Reactome_GSEA_n_bulk  GSEA_shared_vs_5yr_bulk  GSEA_jaccard_vs_5yr_bulk
#       4 yrs        1                        18                   2                  4                   36                      13                    0.210
#       6 yrs        1                         8                   1                  5                   67                      27                    0.342
#       7 yrs        1                        11                   4                  2                   45                      23                    0.377
#       8 yrs        1                         0                   0                  0                   21                       8                    0.154

# Conclusions
# 4yrs: Exploratory low-cutoff reference
# 6yrs: Detailed sensitivity cutoff used in the main sensitivity blocks above
# 7yrs: Exploratory high-cutoff reference
# 8yrs: Exploratory high-cutoff reference


###─────────────────────────────────────── D: Subgroup analysis--------------------------------------------------------

#  0. Settings and packages-----------------------------------------------------

proj_root_bulk <- "C:/DMD_project"
cache_dir_bulk <- file.path(proj_root_bulk, "cache_bulk_revision")

cutoff_C1_bulk        <- 5
p_cutoff_gsea_C1_bulk <- 0.05
show_top_gsea_C1_bulk <- 30
x_breaks_gsea_C1_bulk <- seq(-2, 2, 1)
seed_C1_bulk          <- 123

suppressPackageStartupMessages({
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(matrixStats)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(ReactomePA)
  library(BiocParallel)
  library(enrichplot)
  library(ggridges)
  library(gridExtra)
  library(grid)
  library(VennDiagram)
})

sig3_C1_bulk <- function(x, sig = 3) {
  vapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    if (v == 0) return("0")
    dec_C1_bulk <- max(sig - 1 - floor(log10(abs(v))), 0)
    formatC(v, digits = dec_C1_bulk, format = "f", drop0trailing = FALSE)
  }, FUN.VALUE = character(1))
}

safe_div_C1_bulk <- function(num_bulk, den_bulk) {
  if (is.na(den_bulk) || den_bulk == 0) return(NA_real_)
  num_bulk / den_bulk
}

ridgeplot_C1_bulk <- function(gsea_obj_bulk,
                              showCategory_bulk = 20,
                              pathway_text_size_bulk = 4,
                              x_breaks_bulk = seq(-2, 2, 1)) {
  enrichplot::ridgeplot(
    gsea_obj_bulk,
    showCategory = showCategory_bulk,
    fill         = "p.adjust"
  ) +
    ggplot2::scale_fill_gradient(
      low  = "#b2182b",
      high = "#2166ac",
      name = "FDR"
    ) +
    ggplot2::scale_x_continuous(
      name   = "Running enrichment score",
      breaks = x_breaks_bulk,
      labels = sprintf("%.0f", x_breaks_bulk)
    ) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = pathway_text_size_bulk),
      axis.title.y = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}


#  1. Load revised bulk stages--------------------------------------------------

stage1_bulk <- readRDS(file.path(cache_dir_bulk, "samples_meta_bulk.rds"))
stage5_bulk <- readRDS(file.path(cache_dir_bulk, "normalized_imputed_bulk.rds"))
stage6_bulk <- readRDS(file.path(cache_dir_bulk, "stage6_combat_only_bulk.rds"))

meta_aligned_bulk <- stage1_bulk$meta_aligned_bulk
expr_QN_bulk      <- stage5_bulk$expr_QN_bulk
Cohort_filt_bulk  <- stage5_bulk$Cohort_filt_bulk
combat_bulk       <- stage6_bulk$combat_bulk
Group_imp_bulk    <- stage6_bulk$Group_imp_bulk
design_bulk       <- stage6_bulk$design_bulk
v_bulk            <- stage6_bulk$v_bulk
fit_treat_bulk    <- stage6_bulk$fit_treat_bulk
deg_bulk          <- stage6_bulk$deg_bulk


#  2. Extract the C1 subgroup from the within-cohort quantile-normalized matrix-----

expr_C1_bulk <- expr_QN_bulk[, Cohort_filt_bulk == "C1", drop = FALSE]

meta_C1_bulk <- meta_aligned_bulk[match(colnames(expr_C1_bulk), meta_aligned_bulk$PID), , drop = FALSE]
stopifnot(all(meta_C1_bulk$PID == colnames(expr_C1_bulk)))

Group_C1_bulk <- factor(
  ifelse(meta_C1_bulk$Age < cutoff_C1_bulk, "Early", "Late"),
  levels = c("Early", "Late")
)

cat("C1 subgroup counts:\n")
print(table(Group_C1_bulk))
cat("Genes in expr_C1_bulk:", nrow(expr_C1_bulk), "\n")
cat("Samples in expr_C1_bulk:", ncol(expr_C1_bulk), "\n")
# Group_C1_bulk
# Early  Late 
# 12     5
# Genes in expr_C1_bulk: 22171
# Samples in expr_C1_bulk: 17

#  3. limma-vooma for the C1 subgroup-------------------------------------------

design_C1_bulk <- model.matrix(~ 0 + Group_C1_bulk)
colnames(design_C1_bulk) <- levels(Group_C1_bulk)

contrast_C1_bulk <- matrix(
  c(-1, 1),
  ncol = 1,
  dimnames = list(colnames(design_C1_bulk), "Late_vs_Early_C1_bulk")
)

v_C1_bulk <- limma::vooma(expr_C1_bulk, design_C1_bulk, plot = FALSE)

fit_C1_bulk <- limma::lmFit(v_C1_bulk, design_C1_bulk)
fit_C1_bulk <- limma::contrasts.fit(fit_C1_bulk, contrast_C1_bulk)
fit_C1_bulk <- limma::eBayes(fit_C1_bulk, trend = TRUE)

deg_C1_bulk <- limma::topTable(
  fit_C1_bulk,
  coef    = "Late_vs_Early_C1_bulk",
  p.value = 0.05,
  number  = Inf
)

cat("DEG (C1 only):", nrow(deg_C1_bulk), "\n")
# DEG (C1 only): 103


ann_col_C1_bulk <- data.frame(
  Group  = Group_C1_bulk,
  Cohort = factor("C1", levels = c("C1", "C2", "C3", "C4")),
  row.names = colnames(expr_C1_bulk)
)

ann_cols_C1_bulk <- list(
  Group  = c(Early = "#2C7BB6", Late = "#D7191C"),
  Cohort = c(C1 = "black", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")
)


#  4. Heatmap of subgroup DEGs--------------------------------------------------

{
  hit_genes_C1_bulk <- base::intersect(rownames(deg_C1_bulk), rownames(expr_C1_bulk))
  if (!length(hit_genes_C1_bulk)) {
    stop("No subgroup DEG rows were found in expr_C1_bulk.")
  }
  
  hit_sub_C1_bulk <- expr_C1_bulk[hit_genes_C1_bulk, , drop = FALSE]
  hit_sub_z_C1_bulk <- t(scale(t(hit_sub_C1_bulk)))
  
  sym_C1_bulk <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = rownames(hit_sub_z_C1_bulk),
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )
  
  rownames(hit_sub_z_C1_bulk) <- make.unique(
    ifelse(is.na(sym_C1_bulk) | sym_C1_bulk == "",
           rownames(hit_sub_z_C1_bulk),
           sym_C1_bulk)
  )
  
  pheatmap::pheatmap(
    hit_sub_z_C1_bulk,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    annotation_col    = ann_col_C1_bulk,
    annotation_colors = ann_cols_C1_bulk,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = TRUE,
    fontsize_row      = 4.5,
    main = sprintf("DEGs (%d genes) - C1 subgroup", nrow(hit_sub_z_C1_bulk))
  )
}


#  5. Heatmap of DEGs shared by the main analysis and the subgroup analysis-----

{
  common_ent_C1_bulk <- base::intersect(rownames(deg_bulk), rownames(deg_C1_bulk))
  if (!length(common_ent_C1_bulk)) {
    stop("No common DEGs were found between the main analysis and the subgroup analysis.")
  }
  
  expr_common_C1_bulk <- expr_C1_bulk[common_ent_C1_bulk, , drop = FALSE]
  expr_common_z_C1_bulk <- t(scale(t(expr_common_C1_bulk)))
  expr_common_z_C1_bulk <- pmin(pmax(expr_common_z_C1_bulk, -2.5), 2.5)
  
  sym_common_C1_bulk <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = rownames(expr_common_z_C1_bulk),
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )
  
  keep_common_C1_bulk <- !is.na(sym_common_C1_bulk) & sym_common_C1_bulk != ""
  expr_common_z_C1_bulk <- expr_common_z_C1_bulk[keep_common_C1_bulk, , drop = FALSE]
  sym_common_C1_bulk <- sym_common_C1_bulk[keep_common_C1_bulk]
  rownames(expr_common_z_C1_bulk) <- make.unique(sym_common_C1_bulk)
  
  pheatmap::pheatmap(
    expr_common_z_C1_bulk,
    annotation_col    = ann_col_C1_bulk[colnames(expr_common_z_C1_bulk), , drop = FALSE],
    annotation_colors = ann_cols_C1_bulk,
    color             = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = TRUE,
    fontsize_row      = 4.5,
    border_color      = "grey60",
    main = sprintf("Common %d DEGs (Main and Subgroup)", nrow(expr_common_z_C1_bulk))
  )
}


#  6. logFC concordance between the main analysis and the subgroup analysis-----

{
  deg_main_ent_bulk <- deg_bulk
  rownames(deg_main_ent_bulk) <- sub("_.*", "", rownames(deg_main_ent_bulk))
  
  deg_C1_ent_bulk <- deg_C1_bulk
  rownames(deg_C1_ent_bulk) <- sub("_.*", "", rownames(deg_C1_ent_bulk))
  
  tab_main_C1_bulk <- data.frame(
    EntrezID = rownames(deg_main_ent_bulk),
    logFC_main_bulk = deg_main_ent_bulk$logFC,
    stringsAsFactors = FALSE
  )
  
  tab_C1_bulk <- data.frame(
    EntrezID = rownames(deg_C1_ent_bulk),
    logFC_C1_bulk = deg_C1_ent_bulk$logFC,
    stringsAsFactors = FALSE
  )
  
  merge_tab_C1_bulk <- merge(
    tab_main_C1_bulk,
    tab_C1_bulk,
    by = "EntrezID",
    sort = FALSE
  )
  
  plot(
    merge_tab_C1_bulk$logFC_main_bulk,
    merge_tab_C1_bulk$logFC_C1_bulk,
    xlab = "logFC (Main)",
    ylab = "logFC (C1)",
    pch  = 16,
    col  = "grey50",
    main = "logFC concordance: Main vs C1"
  )
  abline(0, 1, col = "red")
}


#  7. Overlap metrics and Venn-style summary------------------------------------

{
  common_id_C1_bulk <- base::intersect(rownames(deg_main_ent_bulk), rownames(deg_C1_ent_bulk))
  n_common_C1_bulk <- length(common_id_C1_bulk)
  
  jaccard_C1_bulk <- safe_div_C1_bulk(
    n_common_C1_bulk,
    nrow(deg_main_ent_bulk) + nrow(deg_C1_ent_bulk) - n_common_C1_bulk
  )
  
  same_dir_C1_bulk <- if (n_common_C1_bulk == 0) {
    NA_real_
  } else {
    sum(sign(deg_main_ent_bulk[common_id_C1_bulk, "logFC"]) ==
          sign(deg_C1_ent_bulk[common_id_C1_bulk, "logFC"])) / n_common_C1_bulk
  }
  
  precision_C1_bulk <- safe_div_C1_bulk(n_common_C1_bulk, nrow(deg_C1_ent_bulk))
  recall_C1_bulk    <- safe_div_C1_bulk(n_common_C1_bulk, nrow(deg_main_ent_bulk))
  
  cat(sprintf(
    "Jaccard = %.3f   Concordance of effect direction = %.2f%%   Precision = %.2f%%   Recall = %.2f%%\n",
    jaccard_C1_bulk,
    same_dir_C1_bulk * 100,
    precision_C1_bulk * 100,
    recall_C1_bulk * 100
  ))
  
  cat(sprintf(
    "Common DEGs : %d   (%.1f%% of C1-DEG, %.1f%% of Main-DEG)\n",
    n_common_C1_bulk,
    precision_C1_bulk * 100,
    recall_C1_bulk * 100
  ))
  
  metrics_C1_bulk <- data.frame(
    Metric = c("Jaccard", "Precision", "Recall", "Concordance"),
    Value  = c(jaccard_C1_bulk, precision_C1_bulk, recall_C1_bulk, same_dir_C1_bulk),
    stringsAsFactors = FALSE
  )
  
  metrics_plot_C1_bulk <- ggplot2::ggplot(metrics_C1_bulk, ggplot2::aes(Metric, Value * 100)) +
    ggplot2::geom_col(fill = "#1f78b4") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", Value * 100)), vjust = -0.3, size = 4) +
    ggplot2::scale_y_continuous(limits = c(0, 100), expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      y = "Percentage / Score",
      x = NULL,
      title = "Overlap and consistency between Main-DEG and C1-DEG"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  
  venn_list_C1_bulk <- VennDiagram::draw.pairwise.venn(
    area1      = nrow(deg_main_ent_bulk),
    area2      = nrow(deg_C1_ent_bulk),
    cross.area = n_common_C1_bulk,
    category   = c(
      sprintf("Main-DEG: %d", nrow(deg_main_ent_bulk)),
      sprintf("C1-DEG: %d", nrow(deg_C1_ent_bulk))
    ),
    fill       = c("#4daf4a", "#377eb8"),
    alpha      = 0.5,
    scaled     = TRUE,
    cat.cex    = 1.2,
    cex        = 1.2,
    cat.pos    = c(180, 0),
    cat.dist   = c(0.03, 0.01),
    cat.just   = list(c(0, 0.5), c(1, 0.5)),
    ind        = FALSE
  )
  
  venn_grob_C1_bulk <- grid::gTree(children = do.call(grid::gList, venn_list_C1_bulk))
  gridExtra::grid.arrange(
    metrics_plot_C1_bulk,
    venn_grob_C1_bulk,
    ncol   = 2,
    widths = c(1.2, 1)
  )
}

# Jaccard = 0.210   Concordance of effect direction = 100.00%   Precision = 77.67%   Recall = 22.35%

#  8. Reactome GSEA for the C1 subgroup-----------------------------------------

{
  coef_main_name_bulk <- colnames(fit_treat_bulk$coefficients)[1]
  geneList_main_bulk <- stats::setNames(
    drop(fit_treat_bulk$coefficients[, coef_main_name_bulk]),
    sub("_.*", "", rownames(fit_treat_bulk$coefficients))
  )
  geneList_main_bulk <- sort(geneList_main_bulk, decreasing = TRUE)
  
  fit_lfc_C1_bulk <- limma::lmFit(v_C1_bulk, design_C1_bulk)
  fit_lfc_C1_bulk <- limma::contrasts.fit(fit_lfc_C1_bulk, contrast_C1_bulk)
  fit_lfc_C1_bulk <- limma::eBayes(fit_lfc_C1_bulk, trend = TRUE)
  
  t_name_C1_bulk <- colnames(fit_lfc_C1_bulk$t)[1]
  geneList_C1_bulk <- stats::setNames(
    drop(fit_lfc_C1_bulk$t[, t_name_C1_bulk]),
    sub("_.*", "", rownames(fit_lfc_C1_bulk$t))
  )
  geneList_C1_bulk <- sort(geneList_C1_bulk, decreasing = TRUE)
  
  BiocParallel::register(BiocParallel::SerialParam())
  
  set.seed(123)
  gsea_react_bulk <- ReactomePA::gsePathway(
    geneList     = geneList_main_bulk,
    organism     = "human",
    pvalueCutoff = p_cutoff_gsea_C1_bulk,
    maxGSSize    = 3000,
    minGSSize    = 10,
    eps          = 0
  )
  
  set.seed(seed_C1_bulk)
  gsea_react_C1_bulk <- ReactomePA::gsePathway(
    geneList     = geneList_C1_bulk,
    organism     = "human",
    pvalueCutoff = p_cutoff_gsea_C1_bulk,
    maxGSSize    = 3000,
    minGSSize    = 10,
    eps          = 0
  )
}

if (nrow(as.data.frame(gsea_react_C1_bulk)) > 0) {
  p_ridge_C1_bulk <- ridgeplot_C1_bulk(
    gsea_obj_bulk          = gsea_react_C1_bulk,
    showCategory_bulk      = show_top_gsea_C1_bulk,
    pathway_text_size_bulk = 4,
    x_breaks_bulk          = x_breaks_gsea_C1_bulk
  ) +
    ggplot2::ggtitle("Reactome enrichment (Subgroup C1)") +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, size = 14, margin = ggplot2::margin(b = 6)),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.title.x= ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 9),
      legend.title= ggplot2::element_text(size = 10)
    )
  print(p_ridge_C1_bulk)
} else {
  message("No Reactome GSEA result was returned for the C1 subgroup.")
}


#  9. Print significant C1 Reactome pathways------------------------------------

{
  tbl_short_C1_bulk <- as.data.frame(gsea_react_C1_bulk) |>
    dplyr::filter(p.adjust < 0.05) |>
    dplyr::select(Description, NES, p.adjust) |>
    dplyr::arrange(p.adjust)
  
  if (nrow(tbl_short_C1_bulk) == 0) {
    message("No significant Reactome pathway was detected in the C1 subgroup at FDR < 0.05.")
  } else {
    tbl_print_C1_bulk <- tbl_short_C1_bulk |>
      dplyr::mutate(p.adjust = sig3_C1_bulk(p.adjust))
    
    print(tbl_print_C1_bulk, row.names = FALSE)
  }
}


# 10. Compare Reactome GSEA between the main analysis and the C1 subgroup------

{
  tbl_main_bulk <- as.data.frame(gsea_react_bulk) |>
    dplyr::select(ID, NES, p.adjust)
  
  tbl_C1_bulk <- as.data.frame(gsea_react_C1_bulk) |>
    dplyr::select(ID, NES, p.adjust) |>
    dplyr::rename(NES_C1_bulk = NES, p_C1_bulk = p.adjust)
  
  sig_main_bulk <- tbl_main_bulk |>
    dplyr::filter(p.adjust < 0.05)
  
  sig_C1_bulk <- tbl_C1_bulk |>
    dplyr::filter(p_C1_bulk < 0.05)
  
  overlap_path_bulk <- dplyr::inner_join(sig_main_bulk, sig_C1_bulk, by = "ID")
  
  jaccard_path_C1_bulk <- safe_div_C1_bulk(
    nrow(overlap_path_bulk),
    nrow(sig_main_bulk) + nrow(sig_C1_bulk) - nrow(overlap_path_bulk)
  )
  
  cor_NES_C1_bulk <- if (nrow(overlap_path_bulk) < 2) {
    NA_real_
  } else {
    stats::cor(overlap_path_bulk$NES, overlap_path_bulk$NES_C1_bulk, method = "spearman")
  }
  
  cat(sprintf(
    "Shared pathways %d | Jaccard %.3f | NES correlation rho = %.3f\n",
    nrow(overlap_path_bulk),
    jaccard_path_C1_bulk,
    cor_NES_C1_bulk
  ))
  
  if (nrow(overlap_path_bulk) > 0) {
    print(
      ggplot2::ggplot(overlap_path_bulk, ggplot2::aes(NES, NES_C1_bulk)) +
        ggplot2::geom_point(size = 2, alpha = 0.7) +
        ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red") +
        ggplot2::labs(
          x = "NES (Main)",
          y = "NES (C1)",
          title = "Pathway NES concordance"
        ) +
        ggplot2::theme_bw()
    )
  }
  
  comp_C1_bulk <- clusterProfiler::compareCluster(
    geneCluster = list(Main = geneList_main_bulk, C1 = geneList_C1_bulk),
    fun         = "gsePathway",
    pvalueCutoff= 0.05,
    organism    = "human"
  )
  
  print(
    enrichplot::dotplot(comp_C1_bulk, showCategory = 20, font.size = 8) +
      ggplot2::ggtitle("Reactome GSEA - Main vs Subgroup (C1)") +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 7, margin = ggplot2::margin(r = 2)),
        plot.title  = ggplot2::element_text(hjust = 0.5, size = 13)
      )
  )
}

# Shared pathways 34 | Jaccard 0.264 | NES correlation rho = 0.960

## Save the entire current workspace-------------------------------------------
getwd()         # example: "C:/DMD_project"
# Specify the save directory (e.g., /cache under the project root)
#save_dir <- "C:/DMD_project/cache"
#dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
#save.image(file = file.path(save_dir, "workspace_2026-04-24.RData"))
#savehistory(file = file.path(save_dir, "workspace_2026-04-24.Rhistory"))

#load("cache/cache/workspace_2026.0424.RData")
