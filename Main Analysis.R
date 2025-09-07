rm(list=ls())

setwd("C:/Users/neuro/OneDrive/大学院研究/統計解析用フォルダ")
setwd("C:/Users/PC_User/OneDrive/大学院研究/統計解析用フォルダ")
getwd()
##  DMD multi‑cohort 解析 : パッケージ & フォルダ初期設定（初回のみ）-----------

## ---------- A. プロジェクトルート（PC 別に選択)
# proj_root <- "C:/Users/neuro/OneDrive/大学院研究/DMD_project"       # ← PC A
# proj_root <- "C:/Users/PC_User/OneDrive/大学院研究/DMD_project"    # ← PC B

# raw_dir <- file.path(proj_root, "exprdata")
# res_dir <- file.path(proj_root, "results")
# dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
# dir.create(res_dir , recursive = TRUE, showWarnings = FALSE)

## ---------- B. R パッケージ保存先（ユーザー権限で書込可な場所
# user_lib <- "C:/R/win-library/4.5-user"  # Windows 用例
# dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
# .libPaths(user_lib)                      # 優先的にここに入れる

## ---------- C. パッケージ一覧
# cran_pkgs <- c("tidyverse","pheatmap","stringr","ggplot2",
#                "igraph","ggraph","digest","RColorBrewer")

# bioc_pkgs <- c("limma","sva","impute","GEOquery","Biobase",
#                "AnnotationDbi","org.Hs.eg.db",
#                "clusterProfiler","ReactomePA","STRINGdb")

# chip_pkgs <- c("hgu133plus2.db","hgu133a.db","hgu133b.db")

## ---------- D. インストール関数
# install_if_missing <- function(pkgs, installer, ...) {
#   need <- setdiff(pkgs, rownames(installed.packages()))
#   if (length(need)) installer(need, ...)
# }

## ---------- E. CRAN
# install_if_missing(
#   cran_pkgs,
#   install.packages,
#   dependencies = TRUE,
#   repos        = "https://cran.rstudio.com",
#   type         = "binary",   # Windows でソースビルド回避
#   quiet        = TRUE
# )

## ---------- F. Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager", repos = "https://cran.rstudio.com")

# install_if_missing(
#   bioc_pkgs,
#   BiocManager::install,
#   ask    = FALSE,
#   update = FALSE
# )

# install_if_missing(
#   chip_pkgs,
#   BiocManager::install,
#   ask    = FALSE,
#   update = FALSE
# )

# message("✅ インストール完了 — 次回からは 01_session_init.R だけ実行してください")



##### A: -------------------------------------------------------データセット構築
# A-1. DMD multi‑cohort 解析 : セッション初期化（毎回実行）------------------------

## ---------- A. プロジェクトルート（PC 別に選択)
proj_root <- "C:/DMD_project"      # ← PC A
# proj_root <- "C:/Users/PC_User/OneDrive/大学院研究/DMD_project"    # ← PC B
{
  raw_dir <- file.path(proj_root, "exprdata")
  res_dir <- file.path(proj_root, "results")
  
  ## ---------- B. パッケージ保存先を優先登録
  user_lib <- "C:/R/win-library/4.5-user"
  .libPaths(c(user_lib, .libPaths()))   # 先頭にユーザーライブラリを挿入
  
  ## ---------- C. 必要ライブラリ読み込み
  suppressPackageStartupMessages({
    library(tidyverse)
    library(limma);          library(sva);            library(impute)
    library(AnnotationDbi);  library(org.Hs.eg.db)
    library(clusterProfiler);library(ReactomePA);     library(STRINGdb)
    library(pheatmap);       library(igraph);         library(ggraph)
  })
  
  suppressPackageStartupMessages({
    library(hgu133plus2.db)
    library(hgu133a.db)
    library(hgu133b.db)
  })
  
  ## ---------- D. 名前空間衝突解消（dplyr 優先)
  if (requireNamespace("conflicted", quietly = TRUE)) {
    conflicted::conflict_prefer("select",     "dplyr")
    conflicted::conflict_prefer("filter",     "dplyr")
    conflicted::conflict_prefer("rowMedians", "matrixStats")
    conflicted::conflict_prefer("rowSds",     "matrixStats")   # ★ 追加
    conflicted::conflict_prefer("union",      "base")
    conflicted::conflict_prefer("intersect",  "base")
    conflicted::conflict_prefer("setdiff",    "base")
    conflicted::conflict_prefer("simplify", "clusterProfiler")
  }
  
  
  ## ---------- E. 共通オプション
  options(stringsAsFactors = FALSE,
          scipen = 999)        # 出力で指数表記を抑制
  
  message("✅ セッション初期化完了 — 解析スクリプトを実行できます")
  
}

# A-2. サンプルテーブル ---------------------------------------------------------
{
  library(GEOquery)
  
  ## -- C1 (GSE109178 / U133Plus2.0) 17 例
  plus_tbl <- data.frame(
    PID = sprintf("P%02d", 1:17),
    GSM = c("GSM2934819","GSM2934820","GSM2934821","GSM2934822","GSM2934823",
            "GSM2934824","GSM2934825","GSM2934826","GSM2934827","GSM2934828",
            "GSM2934829","GSM2934830","GSM2934831","GSM2934832","GSM2934833",
            "GSM2934834","GSM2934835"),
    Age = c(7, 0.9, 4, 1.6, 4, 8, 5, 6, 1.9, 4,3, 3, 1.9, 1,2, 3.5, 7)
  )
  
  ## -- C2 (GSE3307 / U133A+B) 7 例
  ab_tbl <- data.frame(
    PID   = sprintf("P%02d", 18:24),
    GSM_A = c("GSM74377","GSM74378","GSM74379","GSM74380",
              "GSM121357","GSM121361","GSM121363"),
    GSM_B = c("GSM120786","GSM120777","GSM120763","GSM120760",
              "GSM121329","GSM121331","GSM121333"),
    Age   = c(9, 8, 7, 6.5, 5, 9 ,7)
  )
  
  ## C3.  GSE1764  (U133A/B) 3 症例
  ## ────────────────────────────────────────
  ##  Note: 年齢は SOFT 上 "8.3‑year‑old" のみ明示。
  ##        他 2 症例は description に具体値が無いため、
  ##        ひとまず 8.3 で仮置きし、後で自動取得。
  gse1764_tbl <- data.frame(
    PID   = sprintf("P%02d", 25:27),      # 既存 P25 の次から
    GSM_A = c("GSM30669","GSM30670","GSM30671"),
    GSM_B = c("GSM30675","GSM30676","GSM30677"),
    Age   = c(8.3, 8.3, 8.3)              # 不明分は暫定値
  )
  
  ## C4.  GSE1004 / 1007  (U95A–E) 5 症例
  ## ────────────────────────────────────────
  ##  採用症例: 1,2,3,4,6,7,8（患者 3 は D 欠損; 5,9,10,11 除外）
  ##  欠損チップは NA を設定
  u95_tbl <- data.frame(
    PID   = sprintf("P%02d", 28:32),       # P25–29 を割当
    GSM_A = c("GSM15833","GSM15834","GSM15836",
              "GSM15838","GSM15839"),
    GSM_B = c("GSM15923","GSM15924","GSM15926",
              "GSM15928","GSM15929"),
    GSM_C = c("GSM16215","GSM16216","GSM16218",
              "GSM16219","GSM16220"),
    GSM_D = c("GSM16237","GSM16244","GSM16245",
              "GSM16249","GSM16250"),
    GSM_E = c("GSM16299","GSM16302","GSM16304",
              "GSM16305","GSM16306"),
    Age   = c(1.0, 1.5, 3.0, 1.0, 0.8)
  )
  
  ## ────────────────────────────────────────
  ## 3.  確認
  ## ────────────────────────────────────────
  list(gse1764_tbl, u95_tbl)
}
# A-3. Series‑matrix 読み込み-----------------------------------------------------
## ──────────────────────────────────────────────
## 1.  これまでの C1・C2 読み込み（変更なし）
## ──────────────────────────────────────────────
## 先にグローバル変数 raw_dir を作っておく

  proj_root <- "C:/DMD_project"
  raw_dir <- file.path(proj_root, "exprdata")     # ← ここを先に定義
  cache_dir <- file.path(proj_root, "cache")      # 中間キャッシュ
  

  ## 循環しない read_sm() 本体
  read_sm <- function(fname, raw_dir = NULL) {
    if (is.null(raw_dir)) raw_dir <- get("raw_dir", envir = .GlobalEnv)
    path <- file.path(raw_dir, fname)
    
    con  <- gzfile(path, "rt")
    txt  <- suppressWarnings(readLines(con, warn = FALSE))
    close(con)
    
    beg <- grep("!series_matrix_table_begin", txt, fixed = TRUE)[1] + 1
    end <- grep("!series_matrix_table_end",   txt, fixed = TRUE)[1] - 1
    if (is.na(beg) || is.na(end) || end < beg)
      stop("series‑matrix テーブルを検出できません: ", fname)
    
    ## ★ quote = "\"",  stringsAsFactors = FALSE を指定
    mat <- read.delim(textConnection(txt[beg:end]),
                      header = TRUE, sep = "\t",
                      quote  = "\"",                # ← クォートを正しく除去
                      row.names = 1, check.names = FALSE,
                      stringsAsFactors = FALSE)
    
    ## ★ 行名・列名の残余クォートを一括除去
    stripQ <- \(v) gsub('^"|"$', "", v)
    rownames(mat) <- stripQ(rownames(mat))
    colnames(mat) <- stripQ(colnames(mat))
    
    Biobase::ExpressionSet(as.matrix(mat))
  }
  
  
  # 従来コード (C1〜C3)  —— オブジェクト名そのまま
  message("…loading series‑matrix files (C1 / C2)")
  eset_1 <- read_sm("GSE109178_series_matrix.txt.gz")              # C1
  eset_A <- read_sm("GSE3307-GPL96_series_matrix.txt.gz")          # C2A
  eset_B <- read_sm("GSE3307-GPL97_series_matrix.txt.gz")          # C2B
  
  strip_gsm <- function(x) sub("(GSM[0-9]+).*", "\\1", x)
  needs_log2 <- function(m, cut = 100) max(m, na.rm = TRUE) > cut
  
  idx_1 <- match(plus_tbl$GSM, strip_gsm(sampleNames(eset_1)))
  mat_1 <- exprs(eset_1)[, idx_1]; if (needs_log2(mat_1)) mat_1 <- log2(mat_1 + 1)
  
  idx_A <- match(ab_tbl$GSM_A, strip_gsm(sampleNames(eset_A)))
  idx_B <- match(ab_tbl$GSM_B, strip_gsm(sampleNames(eset_B)))
  mat_A <- exprs(eset_A)[, idx_A]; if (needs_log2(mat_A)) mat_A <- log2(mat_A + 1)
  mat_B <- exprs(eset_B)[, idx_B]; if (needs_log2(mat_B)) mat_B <- log2(mat_B + 1)
  
  message("…loading series‑matrix files (C3: GSE1764)")
  eset_C3_A <- read_sm("GSE1764-GPL96_series_matrix.txt.gz")       # U133A
  eset_C3_B <- read_sm("GSE1764-GPL97_series_matrix.txt.gz")       # U133B
  
  idx_C3_A <- match(gse1764_tbl$GSM_A, strip_gsm(sampleNames(eset_C3_A)))
  idx_C3_B <- match(gse1764_tbl$GSM_B, strip_gsm(sampleNames(eset_C3_B)))
  
  mat_C3_A <- exprs(eset_C3_A)[, idx_C3_A, drop = FALSE]
  mat_C3_B <- exprs(eset_C3_B)[, idx_C3_B, drop = FALSE]
  if (needs_log2(mat_C3_A)) mat_C3_A <- log2(mat_C3_A + 1)
  if (needs_log2(mat_C3_B)) mat_C3_B <- log2(mat_C3_B + 1)

  
  # C4 (U95 A–E) も同じ read_sm() でオフライン読込
  message("…loading series‑matrix files (C4: U95 A–E)")
  u95_files <- c(
    A = "GSE1004-GPL8300_series_matrix.txt.gz",
    B = "GSE1007-GPL92_series_matrix.txt.gz",
    C = "GSE1007-GPL93_series_matrix.txt.gz",
    D = "GSE1007-GPL94_series_matrix.txt.gz",
    E = "GSE1007-GPL95_series_matrix.txt.gz"
  )
  eset_U95 <- lapply(u95_files, read_sm)
  
  mat_U95 <- mapply(
    FUN = \(eset, gsm_vec) {
      idx <- match(gsm_vec, strip_gsm(sampleNames(eset)))
      m   <- exprs(eset)[, idx, drop = FALSE]
      if (needs_log2(m)) m <- log2(m + 1)
      m
    },
    eset_U95,
    list(u95_tbl$GSM_A, u95_tbl$GSM_B, u95_tbl$GSM_C,
         u95_tbl$GSM_D, u95_tbl$GSM_E),
    SIMPLIFY = FALSE
  )
  names(mat_U95) <- c("U95_A", "U95_B", "U95_C", "U95_D", "U95_E")

  # 自動選択ヘルパーもローカル版に置換
  read_sm_auto <- function(prefix, platforms) {
    for (gpl in platforms) {
      fn <- file.path(raw_dir, sprintf("%s-%s_series_matrix.txt.gz", prefix, gpl))
      if (file.exists(fn)) return(read_sm(basename(fn)))
    }
    stop("No series‑matrix file found for ", prefix)
  }
  
  


raw_dir <- "C:/DMD_project/exprdata"


  eset_GSE1004 <- read_sm_auto("GSE1004", c("GPL8300","GPL91"))
  test_eset <- read_sm("GSE1004-GPL8300_series_matrix.txt.gz")
  head(sampleNames(test_eset))
  
  # 期待する GSM（縦ベクトルにまとめる）
  gsm_expect <- unlist(u95_tbl[, paste0("GSM_", c("A","B","C","D","E"))], use.names = FALSE)
  gsm_expect <- na.omit(gsm_expect)                         # 欠損 NA を除く
  gsm_expect <- strip_gsm(gsm_expect)                       # “.CEL.gz” 等を除去
  
  # 実際にファイルに含まれている GSM
  gsm_in_file <- strip_gsm(sampleNames(test_eset))
  
  # ファイルに無い GSM があるか？
  setdiff(gsm_expect, gsm_in_file)
  
  # ファイルごとに期待 GSM を切り分け
  expect_A <- u95_tbl$GSM_A
  expect_B <- u95_tbl$GSM_B
  expect_C <- u95_tbl$GSM_C
  expect_D <- u95_tbl$GSM_D
  expect_E <- u95_tbl$GSM_E
  
  # 例) チップ A (GSE1004‑GPL8300) のみで確認
  missing_A <- setdiff(strip_gsm(expect_A), gsm_in_file)  # ← 空になれば OK
  
# Entrez 変換 & 行合わせ
  library(AnnotationDbi)  # map_to_entrez 内で使用
  
  ## ── 4‑1  ヘルパー ───────────────────────────────────────────────────────────
  map_to_entrez <- function(expr, chip_pkg){
    if (all(grepl("^[0-9]+$", rownames(expr)))) return(expr)   # 既に Entrez の場合
    if (!requireNamespace(chip_pkg, quietly = TRUE))
      stop(chip_pkg, " が読み込めません（Bioconductor install 必要）")
    chip_db <- get(chip_pkg, envir = asNamespace(chip_pkg))
    mp <- AnnotationDbi::select(chip_db,
                                keys   = rownames(expr),
                                columns= "ENTREZID",
                                keytype= "PROBEID") |>
      na.omit()
    expr <- expr[mp$PROBEID, , drop = FALSE]
    rownames(expr) <- mp$ENTREZID
    rowsum(expr, group = rownames(expr)) / as.numeric(table(mp$ENTREZID))
  }
  
  pad <- \(mat, genes){
    o <- matrix(NA_real_, length(genes), ncol(mat),
                dimnames = list(genes, colnames(mat)))
    hit <- intersect(rownames(mat), genes); o[hit, ] <- mat[hit, ]; o
  }
  
  ## ── 4‑2  コホートごとの Entrez 変換 ──────────────────────────────────────────
  expr_1     <- map_to_entrez(mat_1,      "hgu133plus2.db") # C1
  
  expr_c2_A  <- map_to_entrez(mat_A,      "hgu133a.db")     # C2  (U133A)
  expr_c2_B  <- map_to_entrez(mat_B,      "hgu133b.db")     #     (U133B)
  
  expr_c3_A  <- map_to_entrez(mat_C3_A,   "hgu133a.db")     # C3  (U133A)
  expr_c3_B  <- map_to_entrez(mat_C3_B,   "hgu133b.db")     #     (U133B)
  
  needed_pkgs <- c(
    # 既存 U133 系
    "hgu133plus2.db", "hgu133a.db", "hgu133b.db",
    # 追加 U95 系
    "hgu95av2.db", "hgu95b.db", "hgu95c.db", "hgu95d.db", "hgu95e.db"
  )
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")          # CRAN 版を導入
  
  needed_pkgs <- c("hgu133plus2.db", "hgu133a.db", "hgu133b.db",
                   "hgu95av2.db", "hgu95b.db", "hgu95c.db",
                   "hgu95d.db", "hgu95e.db")
  
  # 未インストールだけを抽出
  to_install <- needed_pkgs[!vapply(needed_pkgs, requireNamespace,
                                    FUN.VALUE = logical(1), quietly = TRUE)]
  
  if (length(to_install)) {
    BiocManager::install(to_install, ask = FALSE, update = FALSE)
  }
  
  expr_u95_A <- map_to_entrez(mat_U95$U95_A, "hgu95av2.db") # C4  (U95A)
  expr_u95_B <- map_to_entrez(mat_U95$U95_B, "hgu95b.db")   #     (U95B)
  expr_u95_C <- map_to_entrez(mat_U95$U95_C, "hgu95c.db")   #     (U95C)
  expr_u95_D <- map_to_entrez(mat_U95$U95_D, "hgu95d.db")   #     (U95D)
  expr_u95_E <- map_to_entrez(mat_U95$U95_E, "hgu95e.db")   #     (U95E)
  
  ## ── 4‑3  U133 系（C2・C3）を A+B 平均で統合 ────────────────────────────────
  genes_AB_c2 <- union(rownames(expr_c2_A), rownames(expr_c2_B))
  expr_c2_AB  <- (pad(expr_c2_A, genes_AB_c2) +
                    pad(expr_c2_B, genes_AB_c2)) / 2
  
  genes_AB_c3 <- union(rownames(expr_c3_A), rownames(expr_c3_B))
  expr_c3_AB  <- (pad(expr_c3_A, genes_AB_c3) +
                    pad(expr_c3_B, genes_AB_c3)) / 2
  
  length(genes_AB_c2) # 19771
  length(genes_AB_c3) # 19771
  
  genes_u133_all <- Reduce(union, list(rownames(expr_1),
                                       genes_AB_c2,
                                       genes_AB_c3))
  length(genes_u133_all)       # 22171
  
  extra_genes <- setdiff(rownames(expr_1), genes_AB_c2)   # Plus2.0 独自分
  length(extra_genes)  #  2400                                  
  
  ## ── 4‑4  U95 系（C4）を A–E 平均で統合 ───────────────────────────────────
  expr_u95_list <- list(expr_u95_A, expr_u95_B, expr_u95_C,
                        expr_u95_D, expr_u95_E)
  genes_u95     <- Reduce(union, lapply(expr_u95_list, rownames))

  ## すべてを同じ行数に pad してから平均
  expr_u95_pad  <- lapply(expr_u95_list, pad, genes = genes_u95)
  
  ## 行列を「足し算＋カウント」で NA を無視して平均
  expr_u95_sum   <- Reduce(`+`, lapply(expr_u95_pad, \(m){ m[is.na(m)] <- 0; m }))
  expr_u95_count <- Reduce(`+`, lapply(expr_u95_pad, \(m) !is.na(m)))
  
  expr_u95_ABCD <- expr_u95_sum / expr_u95_count
  expr_u95_ABCD[expr_u95_count == 0] <- NA  # 全チップ欠損セルを NA
  
  ## ── 4‑5  全コホート共通遺伝子集合 & 結合 ────────────────────────────────
  genes_all <- Reduce(union, list(rownames(expr_1),        # C1
                                  rownames(expr_c2_AB),    # C2
                                  rownames(expr_c3_AB),    # C3
                                  rownames(expr_u95_ABCD)) # C4
  )
  
  expr_mat <- cbind(
    pad(expr_1,        genes_all),     
    pad(expr_c2_AB,    genes_all),     
    pad(expr_u95_ABCD, genes_all),     
    pad(expr_c3_AB,    genes_all)      
  )
  
  ## ── 4‑6  列名を PID に統一 ────────────────────────────────────────────────
  colnames(expr_mat) <- c(plus_tbl$PID,            
                          ab_tbl$PID,              
                          u95_tbl$PID,             
                          gse1764_tbl$PID)         
  
# 5.  メタデータ & グループ定義 (! 年齢カットオフ設定含む-----------------------
  
  ## ── 5‑1  コホート別メタデータ作成 ──────────────────────────────────────────
  ##  representative GSM は各患者の “チップ A” を採用
  ##  （統計解析には PID が主キーになるため GSM はラベル用途のみ）
  
  meta_c1 <- transform(plus_tbl,  GSM = GSM,   Cohort = "C1")[, c("PID","GSM","Age","Cohort")]
  meta_c2 <- transform(ab_tbl,    GSM = GSM_A, Cohort = "C2")[, c("PID","GSM","Age","Cohort")]
  meta_c3 <- transform(gse1764_tbl, GSM = GSM_A, Cohort = "C3")[, c("PID","GSM","Age","Cohort")]
  meta_c4 <- transform(u95_tbl,   GSM = GSM_A, Cohort = "C4")[, c("PID","GSM","Age","Cohort")]
  
  ## 連結
  sample_meta <- rbind(meta_c1, meta_c2, meta_c3, meta_c4)
  
  ## ── 5‑2  expr_mat の列順に合わせて整列 ──────────────────────────────────
  meta_aligned <- sample_meta[ match(colnames(expr_mat), sample_meta$PID), ]
  
  ## ── 5‑3  因子の作成 ──────────────────────────────────────────────────────
  Cohort <- factor(meta_aligned$Cohort,
                   levels = c("C1","C2","C3","C4"))   # 好みで並び替え可


# A-4. 年齢指定①---------------------------------------------------------------------
Group  <- factor(ifelse(meta_aligned$Age < 5, "Early", "Late"),
                 levels = c("Early","Late"))

## ── 5‑4  共通遺伝子 or 測定率フィルタ ────────────────────────────────────
## ① 完全共通 Entrez（4 コホート全てに存在）
  genes_common <- Reduce(intersect, list(rownames(expr_1),
                                        rownames(expr_c2_AB),
                                        rownames(expr_c3_AB),
                                        rownames(expr_u95_ABCD)))
  
  expr_mat_common <- expr_mat[genes_common, ]
  length(expr_mat_common)
  
  ## ② ある程度欠測を許容する場合（少なくとも 6 症例で測定）
  #gene_keep <- rowSums(!is.na(expr_mat)) >= 6
  #expr_mat_filtered <- expr_mat[gene_keep, ]
  
  ## 以降の解析では ①または②、お好みの行列オブジェクトを使用してください
  
  pad <- function(mat, genes){
    out <- matrix(NA_real_, length(genes), ncol(mat),
                  dimnames = list(genes, colnames(mat)))
    hit <- intersect(rownames(mat), genes)
    out[hit, ] <- mat[hit, ]
    out
  }
  
  strip_gsm <- function(x) sub("(GSM[0-9]+).*", "\\1", x)
  
  # 1‑A 5×7 = 35 列行列 u95_long
  genes_u95  <- Reduce(union, lapply(list(expr_u95_A,expr_u95_B,
                                          expr_u95_C,expr_u95_D,expr_u95_E),
                                     rownames))
  u95_long   <- do.call(cbind, lapply(
    list(expr_u95_A,expr_u95_B,expr_u95_C,
         expr_u95_D,expr_u95_E),
    pad, genes = genes_u95))
  colnames(u95_long) <- strip_gsm(colnames(u95_long))
  
  ## --- U95 症例平均 ― 欠測 GSM があっても列を残す
  pid_all <- u95_tbl$PID                                # P25–P30 の6症例
  
  expr_u95_avg <- sapply(pid_all, function(pid){
    gsm_vec <- strip_gsm(unlist(
      u95_tbl[u95_tbl$PID == pid,
              paste0("GSM_", c("A","B","C","D","E"))]))
    gsm_vec <- gsm_vec[gsm_vec != "" & !is.na(gsm_vec)] # 空白除外
    idx     <- match(gsm_vec, colnames(u95_long))
    idx     <- idx[!is.na(idx)]                         # 見つかった列だけ
    if (length(idx) == 0) {
      warning("症例 ", pid,
              " は有効 GSM が 0 枚。全 NA で保持します。")
      return(rep(NA_real_, nrow(u95_long)))
    }
    rowMeans(u95_long[, idx, drop = FALSE], na.rm = TRUE)
  })
  colnames(expr_u95_avg) <- pid_all                     # 列名順固定
  
  # 1‑B GSM→PID
  gsm_pid <- unlist(u95_tbl[,paste0("GSM_",c("A","B","C","D","E"))])
  gsm_pid <- strip_gsm(gsm_pid); pid_vec <- rep(u95_tbl$PID, each = 5)
  map      <- setNames(pid_vec, gsm_pid)[colnames(u95_long)]   # NA は自動除去
  
  # 1‑C 平均
  expr_C4 <- sapply(split(seq_along(map), map), function(idx){
    rowMeans(u95_long[, idx, drop=FALSE], na.rm=TRUE)
  })
  
  ## ── ヘルパー: GSM → PID 変換（U133A/B 共通） ──
  rename_cols_by_pid <- function(mat, gsm_vec, pid_vec){
    stopifnot(length(gsm_vec) == length(pid_vec))          # 同じ長さ
    m <- match(colnames(mat), strip_gsm(gsm_vec))          # GSM →行番号
    newnames <- pid_vec[m]
    if (anyNA(newnames))
      stop("列名に対応する PID が見つかりません: ",
           paste(colnames(mat)[is.na(newnames)], collapse = ", "))
    colnames(mat) <- newnames
    mat
  }
  
  ##―― C2（7 症例）
  expr_c2_A <- rename_cols_by_pid(expr_c2_A,
                                  ab_tbl$GSM_A, ab_tbl$PID)
  expr_c2_B <- rename_cols_by_pid(expr_c2_B,
                                  ab_tbl$GSM_B, ab_tbl$PID)
  
  ##―― C3（3 症例）
  expr_c3_A <- rename_cols_by_pid(expr_c3_A,
                                  gse1764_tbl$GSM_A, gse1764_tbl$PID)
  expr_c3_B <- rename_cols_by_pid(expr_c3_B,
                                  gse1764_tbl$GSM_B, gse1764_tbl$PID)
  avg_AB <- function(A, B){
    ## すでに列名は PID になっている前提
    g   <- union(rownames(A), rownames(B))      # 遺伝子の全集合
    pid <- union(colnames(A), colnames(B))      # 患者 ID の全集合
    
    X <- pad(A, g)[ , pid, drop = FALSE]        # 列順を pid に固定
    Y <- pad(B, g)[ , pid, drop = FALSE]
    
    cnt <- (!is.na(X)) + (!is.na(Y))            # 実測数 0/1/2
    out <- (replace(X, is.na(X), 0) +
              replace(Y, is.na(Y), 0)) / pmax(cnt, 1)
    out[cnt == 0] <- NA                         # 両方 NA は NA
    out
  }
  
  
  expr_C2 <- avg_AB(expr_c2_A, expr_c2_B)
  expr_C3 <- avg_AB(expr_c3_A, expr_c3_B)
  
  stopifnot(all(colSums(!is.na(expr_C2)) > 0),
            all(colSums(!is.na(expr_C3)) > 0))
  
  ## データ統合(! 年齢指定)-----------------------------------------------------
  ##─────────────────────────────────────
  ##  A. ユーティリティ
  ##─────────────────────────────────────
  strip_gsm <- function(x) sub("(GSM[0-9]+).*", "\\1", x)
  
  pad_mat <- function(mat, genes){
    out <- matrix(NA_real_, length(genes), ncol(mat),
                  dimnames = list(genes, colnames(mat)))
    hit <- intersect(rownames(mat), genes)
    out[hit, ] <- mat[hit, ]
    out
  }
  
  ##─────────────────────────────────────
  ##  B. GSM → PID 置換
  ##─────────────────────────────────────
  rename_cols <- function(mat, gsm_vec, pid_vec){
    gsm_map <- setNames(pid_vec, strip_gsm(gsm_vec))
    cn      <- strip_gsm(colnames(mat))
    
    ## 置換
    is_gsm  <- cn %in% names(gsm_map)
    cn[is_gsm] <- gsm_map[cn[is_gsm]]
    
    ## 見つからなかった GSM を警告
    lost <- !(cn %in% pid_vec)
    if (any(lost))
      warning("未マッピング列: ", paste(colnames(mat)[lost], collapse = ", "))
    
    colnames(mat) <- cn
    mat
  }
  
  
  ##─────────────────────────────────────
  ##  C. 2チップ平均（U133A/B など）
  ##─────────────────────────────────────
  avg_two <- function(A, B){
    genes <- union(rownames(A), rownames(B))
    ids   <- union(colnames(A), colnames(B))
    X <- pad_mat(A, genes)[, ids, drop = FALSE]
    Y <- pad_mat(B, genes)[, ids, drop = FALSE]
    cnt <- (!is.na(X)) + (!is.na(Y))
    out <- (replace(X, is.na(X), 0) + replace(Y, is.na(Y), 0)) / pmax(cnt, 1)
    out[cnt == 0] <- NA
    out
  }
  
  head(colnames(expr_c2_B))
  # → もし "P18" など PID が並んでいれば、すでにリネーム済み
  
  
  ##―― C1 (U133 Plus2)
  expr_C1 <- map_to_entrez(mat_1, "hgu133plus2.db") |>
    rename_cols(plus_tbl$GSM, plus_tbl$PID)
  
  ##―― C2 (U133 A + B)
  mat_C2A <- rename_cols(expr_c2_A, ab_tbl$GSM_A, ab_tbl$PID)  # 2 回目ならそのまま返る
  mat_C2B <- rename_cols(expr_c2_B, ab_tbl$GSM_B, ab_tbl$PID)
  expr_C2 <- avg_two(mat_C2A, mat_C2B)                         # 7 列 × 遺伝子数
  
  ##―― C3 (U133 A + B)
  mat_C3A <- rename_cols(expr_c3_A, gse1764_tbl$GSM_A, gse1764_tbl$PID)
  mat_C3B <- rename_cols(expr_c3_B, gse1764_tbl$GSM_B, gse1764_tbl$PID)
  expr_C3 <- avg_two(mat_C3A, mat_C3B)                         # 3 列
  

  ##―― C4 (U95 A–E) ────────────────────
  gsm_lists <- list(
    u95_tbl$GSM_A,
    u95_tbl$GSM_B,
    u95_tbl$GSM_C,
    u95_tbl$GSM_D,
    u95_tbl$GSM_E
  )
  
  u95_raw  <- list(expr_u95_A, expr_u95_B, expr_u95_C,
                   expr_u95_D, expr_u95_E)
  
  u95_list <- Map(
    \(mat, gsm) rename_cols(mat, gsm, u95_tbl$PID),
    u95_raw, gsm_lists
  )
  
  ## ❶ 遺伝子の全集合
  genes_u95 <- Reduce(union, lapply(u95_list, rownames))
  length(genes_u95) # 20856
  
  ## ❷ 列（PID）の全集合 
  pid_u95 <- u95_tbl$PID
  
  ## ❸ 行・列ともに pad する関数
  pad_rc <- function(mat, genes, pid){
    ## ─ 行をそろえる ─
    m <- matrix(NA_real_, length(genes), ncol(mat),
                dimnames = list(genes, colnames(mat)))
    hit <- intersect(rownames(mat), genes)
    m[hit, ] <- mat[hit, ]
    
    ## ─ 列をそろえる（欠け PID を NA 列で補完） ─
    miss <- setdiff(pid, colnames(m))
    if (length(miss))
      m <- cbind(m,
                 matrix(NA_real_, nrow(m), length(miss),
                        dimnames = list(rownames(m), miss)))
    m[ , pid, drop = FALSE]        # 列順を固定
  }
  
  ## ❹ 5チップすべてを行・列 pad
  u95_pad <- lapply(u95_list, pad_rc, genes = genes_u95, pid = pid_u95)
  
  ## ❺ 合計とカウントを並行して計算
  expr_sum   <- Reduce(`+`, lapply(u95_pad, \(m) replace(m, is.na(m), 0)))
  expr_count <- Reduce(`+`, lapply(u95_pad, \(m) !is.na(m)))
  
  expr_C4 <- expr_sum / pmax(expr_count, 1)
  expr_C4[expr_count == 0] <- NA   # 5 枚とも欠測なら NA
  
  
  cohorts <- list(C1 = expr_C1, C2 = expr_C2, C3 = expr_C3, C4 = expr_C4)
  
  genes_all <- Reduce(union, lapply(cohorts, rownames))
  expr_mat  <- do.call(cbind,
                       lapply(cohorts, pad_mat, genes = genes_all))

  ## コホート情報・PID ベクトル
  PID    <- colnames(expr_mat)
  Cohort <- factor(
    rep(names(cohorts), vapply(cohorts, ncol, integer(1))),
    levels = c("C1","C2","C3","C4")
  )


# A-5. 年齢指定②---------------------------------------------------------------------
AgeL  <- setNames(meta_aligned$Age, meta_aligned$PID)
Group  <- factor(ifelse(AgeL[PID] < 5, "Early", "Late"),
                 levels = c("Early","Late"))

stopifnot(!any(colSums(!is.na(expr_mat)) == 0))  # 全 NA 列なし

# ★ ここでは「6検体以上で実測」フィルタを掛けるだけ
# expr <- expr_mat[rowSums(!is.na(expr_mat)) >= 6, ]

# ★ kNN 補完はまだ行わない

## 重複チェック（NA を含んでもハッシュは取れる）

library(digest)
find_dup_cols <- function(mat, algo = "xxhash64"){
  h <- apply(mat, 2, digest, algo = algo)
  which(duplicated(h) | duplicated(h, fromLast = TRUE))
}
stopifnot(length(find_dup_cols(expr_mat)) == 0)  # 重複があればここで停止

length(expr_mat) # 709472

# expr は 18509 × 32 の統合 ExpressionSet
total_cells  <- length(expr_mat)          # 22171 × 32 = 592288
na_cells     <- sum(is.na(expr_mat))  ;  na_cells   # 38365
na_fraction  <- na_cells / total_cells; na_fraction    # 0.05407543 (5.4 %)

  ##  STEP‑1  コホート別 QN ＋ 行合わせ
  
  ## ① 6症例以上で測定された遺伝子を残す
  expr_filt <- expr_mat[rowSums(!is.na(expr_mat)) >= 6, ]
  genes_all <- rownames(expr_filt)
  
  ## ② expr_filt に対応する Cohort ベクトルを作り直す
  pid_in_mat <- colnames(expr_filt)
  Cohort_filt <- Cohort[match(pid_in_mat, PID)]      # 同一長・同順
  
  ## ③ コホートごとに QN
  if (!requireNamespace("preprocessCore", quietly = TRUE))
    BiocManager::install("preprocessCore", ask = FALSE, update = FALSE)
  library(preprocessCore)
  
  pad_mat <- function(mat, genes){
    out <- matrix(NA_real_, length(genes), ncol(mat),
                  dimnames = list(genes, colnames(mat)))
    hit <- intersect(rownames(mat), genes)
    out[hit, ] <- mat[hit, ]
    out
  }
  
  table(Cohort_filt)       # expr_filt に残った列数を確認
  # C1 C2 C3 C4 
  # 17  7  3  5
  
  
  expr_QN_list <- lapply(levels(Cohort_filt), function(b){
    pids <- pid_in_mat[Cohort_filt == b]
    m    <- expr_filt[ , pids, drop = FALSE]
    m    <- m[rowSums(!is.na(m)) > 0, , drop = FALSE]
    limma::normalizeBetweenArrays(m, method = "quantile") |>
      pad_mat(genes_all)
  })
  
  expr_QN <- do.call(cbind, expr_QN_list)[ , pid_in_mat]   # 列順を維持
  
  ## QC
  stopifnot(
    dim(expr_QN) == c(length(genes_all), length(pid_in_mat)),
    identical(rownames(expr_QN), genes_all)
  )
  
  
  ##  STEP‑2  kNN 補完 → 重複列チェック
  
  if (!requireNamespace("qsmooth", quietly = TRUE))
    BiocManager::install("qsmooth", ask = FALSE, update = FALSE)
  library(impute); library(digest)
  
  ## ❶ 各コホートで ≥1 実測遺伝子のみ残す
  keep <- Reduce(`&`, lapply(levels(Cohort_filt), function(b)
    rowSums(!is.na(expr_QN[, Cohort_filt == b, drop = FALSE])) > 0))
  expr_sel <- expr_QN[keep, ]
  
  ## ❷ kNN (colmax = 0.99 で全欠測列を拒否)
  expr_imp <- impute.knn(expr_sel, colmax = 0.99, rowmax = 1)$data
  
  ## ❸ 完全複製列チェック
  find_dup_cols <- function(mat, algo = "xxhash64"){
    h <- apply(mat, 2, digest, algo = algo)
    which(duplicated(h) | duplicated(h, fromLast = TRUE))
  }
  stopifnot(length(find_dup_cols(expr_imp)) == 0)
  
  ## ❹ 残 NA を行平均で埋め、なお残れば行除去
  na_pos <- which(is.na(expr_imp), arr.ind = TRUE)
  if (nrow(na_pos))
    expr_imp[na_pos] <- rowMeans(expr_imp, na.rm = TRUE)[na_pos[, "row"]]
  expr_imp <- expr_imp[rowSums(is.na(expr_imp)) == 0, ]
  
  dup1 <- find_dup_cols(expr_imp)
  stopifnot(length(dup1) == 0)

## 以降: qsmooth → ComBat → limma へ

# A-6. qsmooth へ ------------------------------------------------------------
{
library(qsmooth)
expr_qs <- qsmoothData(qsmooth(expr_imp, group_factor = Cohort))
length(expr_qs)
combat   <- ComBat(expr_qs,
                   batch = Cohort,
                   mod   = model.matrix(~ Group),
                   ref.batch = "C1",
                   par.prior = TRUE,
                   mean.only = FALSE)

cat("NA なし遺伝子行数:", nrow(expr_qs), "\n")  # 例: ~18k
# NA なし遺伝子行数: 18509

library(sva)
combat1 <- ComBat(expr_qs,
                  batch     = Cohort,
                  mod       = model.matrix(~ Group),    # Group を明示
                  ref.batch = "C1",                     # 最大枚数を基準
                  par.prior = TRUE,
                  mean.only = FALSE)                    # 分散も補正


#  limma::removeBatchEffect で残滓を微調整  --------------------
library(limma)
combat <- removeBatchEffect(combat1,
                            batch = Cohort,
                            design = model.matrix(~ Group))

##―― ComBat 後
dup2 <- find_dup_cols(combat)
stopifnot(length(dup2) == 0)
}
# A-7. QC プロット  ------------------------------------------------------
# PCA
{
  plotPCA_micro <- function(expr_qs, group_vec, cohort_vec,
                            col_pal = c(Early="#2C7BB6", Late="#D7191C"),
                            pch_pal = c(C1=17, C2=0, C3=15, C4=1),
                            main = "PCA after ComBat"){
    
    pc <- prcomp(t(expr_qs))
    plot(pc$x[,1], pc$x[,2],
         col = col_pal[as.character(group_vec)],
         pch = pch_pal[as.character(cohort_vec)],
         xlab = sprintf("PC1 (%.1f%%)", 100*pc$sdev[1]^2/sum(pc$sdev^2)),
         ylab = sprintf("PC2 (%.1f%%)", 100*pc$sdev[2]^2/sum(pc$sdev^2)),
         main = main, cex = 1.2)
    legend("topleft",  legend = names(col_pal), col = col_pal, pch = 16, bty = "n")
    legend("bottomright", legend = names(pch_pal), pch = pch_pal, bty = "n")
  }
  plotPCA_micro(combat, Group, Cohort)
}


# RLE
{
  plotRLE_micro <- function(expr_qs, group_vec,
                            col_pal = c(Early="#2C7BB6", Late="#D7191C"),
                            main = "RLE (after normalisation / ComBat)"){
    
    med <- apply(expr_qs, 1, median, na.rm = TRUE)
    rle <- sweep(expr_qs, 1, med, "-")
    boxplot(rle, outline = FALSE,
            col = col_pal[as.character(group_vec)],
            las = 2, ylab = "Relative log2 expression", main = main)
    abline(h = 0, lwd = 2)
  }
  plotRLE_micro(combat, Group)
}


# IQRs
{
  med   <- apply(combat, 1, median, na.rm = TRUE)
  rle   <- sweep(combat, 1, med, "-")
  IQRs  <- apply(rle, 2, IQR)
  summary(IQRs) }  # IQR が 0.25–0.6 程度なら OK
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1677  0.2314  0.3064  0.2950  0.3439  0.4303

# 6yr 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2446  0.3076  0.3334  0.3311  0.3600  0.4258 

##### B: -------------------------------------------------------------limma 解析------------------------------------------------------------------------
# B-1. limma-treat--------------------------------------------------------------
{
  design <- model.matrix(~0+Group); colnames(design) <- levels(Group)
  v   <- vooma(combat, design, plot = FALSE)
  fit <- lmFit(v, design) |>
    contrasts.fit(cbind(Late_vs_Early = c(-1,1))) |>
    eBayes() |>
    treat(lfc = 0.30)
  deg <- topTreat(fit, p.value = 0.05, number = Inf)
  cat("DEG:", nrow(deg), "\n")
  }
# 5yr: DEG 261
# 4yr: DEG: 7
# 6yr: DEG: 108
# 7yr: DEG 86
# 8yr: DEG: 399

table(Group)
#    Early  Late 
# 5yr  17     15
# 4yr  
# 6yr  19     13
# 7yr  21     11
# 8yr  25      7

## ====== 複製チェック
{
  M <- as.data.frame(t(expr_qs))          # sample × gene
  ident <- combn(colnames(expr_qs), 2, \(p){
    all(expr_qs[,p[1]] == expr_qs[,p[2]])
  })
  which(ident)}      # TRUE が出たペアが完全複製

# graphics.off() 



# B-2. heatmap: DEG 300 行 & “ばらつきフィルタ -------------------------------
{ library(pheatmap)
  library(matrixStats)         # rowSds 用
  
  N <- 300                       # 好みで 100–500
  deg.sel <- head(deg, N)
  sel_id  <- rownames(deg.sel)
  
  ## rowSds で極端に平坦な行を除去
  expr_sub <- combat[sel_id, ]
  expr_sub <- expr_sub[ rowSds(expr_sub) > 0.15 , ]   # SD 0.15 未満を除外
  
  ## --- ② 列アノテーション
  ann_col <- data.frame(
    Group  = Group,
    Cohort = factor(Cohort, levels = c("C1","C2","C3","C4"))
  )
  rownames(ann_col) <- colnames(expr_sub)
  
  ann_cols <- list(
    Group  = c(Early = "#2C7BB6", Late = "#D7191C"),
    Cohort = c(C1 = "black", C2 = "grey40",
               C3 = "#4DAF4A", C4 = "purple")
  )
  
  ## --- ③ 色と break 設定
  my_col <- colorRampPalette(c("navy", "white", "firebrick3"))(101)
  my_brk <- seq(-2.5, 2.5, length = 101)   # 行 Z‑score が大体 ±2 に入る
  
  ## --- ④ pheatmap
  col_NA_rate <- colMeans(is.na(expr_mat))     # 1.00 なら全 NA
  col_NA_rate[ col_NA_rate == 1 ]
  
  dup_cols <- duplicated(as.data.frame(t(expr_sub)))
  if(any(dup_cols)){
    message("複製列が ", sum(dup_cols), " 本あり除外します")
    expr_sub <- expr_sub[ , !dup_cols]
    ann_col  <- ann_col [!dup_cols, ]
  }
  
  ids <- rownames(expr_sub)
  sym <- mapIds(org.Hs.eg.db, ids, "SYMBOL", "ENTREZID", multiVals="first")
  rownames(expr_sub) <- ifelse(is.na(sym), ids, sym)
  
  pheatmap(
    expr_sub,
    scale              = "row",          # 行 Z‑score
    color              = my_col,
    breaks             = my_brk,
    annotation_col     = ann_col,
    annotation_colors  = ann_cols,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    show_colnames      = T,
    show_rownames      = F,           
    fontsize_row       = 3,   
    fontsize_col       = 8,
    border_color    = F,
    main = sprintf("DEG heatmap | %d genes × %d samples",
                   nrow(expr_sub), ncol(expr_sub))
  )
}


# B-3. FDR TOP 100 Heatmap------------------------------------------------------------
{
suppressPackageStartupMessages({
  library(matrixStats); library(AnnotationDbi); library(org.Hs.eg.db); library(pheatmap)
})

### 0. 前処理
deg_df <- as.data.frame(deg)          # ← deg を data.frame 化

### 1. FDR 列名を安全に取得
fdr_candidates <- c("FDR", "adj.P.Val", "padj", "qvalue", "p.adj", "FDR.BH")
fdr_col <- intersect(fdr_candidates, colnames(deg_df))

if (length(fdr_col) == 0)
  stop("FDR に相当する列が見つかりません。列名を確認してください。")
fdr_col <- fdr_col[1]                 # 複数あれば最初を採用

### 2. FDR 昇順で並べ替え & 上位 100
deg_df_sorted <- deg_df[order(deg_df[[fdr_col]]), , drop = FALSE]
top100_id     <- head(rownames(deg_df_sorted), 100)

## ここで “top100” を定義しておく
top100 <- intersect(top100_id, rownames(combat))   # 念のため一致確認

## 2. 発現行列を抽出し、行 SD でフィルタ
expr_top <- combat[top100, , drop = FALSE]
expr_top <- expr_top[rowSds(expr_top) > 0.15, ] # 視覚的に平坦な行を除外
n_gene   <- nrow(expr_top)                      # 実際に残った行数

## 3. 行 Z‑score  & クリップ
expr_z <- t(scale(t(expr_top)))
expr_z <- pmin(pmax(expr_z, -2.5), 2.5)

## 4. Entrez → SYMBOL
sym <- mapIds(org.Hs.eg.db, rownames(expr_z), "SYMBOL", "ENTREZID", multiVals="first")
expr_z <- expr_z[!is.na(sym), ]
rownames(expr_z) <- sym[!is.na(sym)]

## 5. アノテーション
ann_col <- data.frame(Group = Group,
                      Cohort = factor(Cohort, levels = c("C1","C2","C3","C4")))
rownames(ann_col) <- colnames(expr_z)
ann_cols <- list(
  Group  = c(Early = "#2C7BB6", Late = "#D7191C"),
  Cohort = c(C1 = "black", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")
)

## 6. ヒートマップ
pheatmap(
  expr_z,
  scale = "none",
  color = colorRampPalette(c("navy","white","firebrick3"))(100),
  breaks = seq(-3, 3, length = 101),
  annotation_col    = ann_col,
  annotation_colors = ann_cols,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_rownames = F,
  border_color    = F,
  fontsize_row  = 5,   # 行名が潰れないよう小さめ
  main = sprintf("DEG heatmap | FDR Top 100  →  %d genes × %d samples",
                 n_gene, ncol(expr_z))
)
}


# 5歳ではバランスが良い。
# 7歳だとEarlyで高発現の遺伝子は減少、Lateで高発現の遺伝子が増加するというアンバランスになる
# 6歳でもバランスは良い。


# B-4. 遺伝子名表示: 上位 N 行を切り出す---------------------------------------------
# 0.  前提：deg10 は treat(lfc = 0.15) で抽出したデータフレーム
#            行名 = EntrezID, 列 = logFC, adj.P.Val, P.Value, ...
{
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  N <- 300                     # 必要なら変更
  deg.top <- head(deg, N)
  deg.top$EntrezID <- rownames(deg.top)
  
  ## 1) SYMBOL 付与
  sym.map <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys     = deg.top$EntrezID,
                                   columns  = "SYMBOL",
                                   keytype  = "ENTREZID")
  
  ## 2) マージ
  deg.top <- merge(deg.top, sym.map,
                   by.x = "EntrezID", by.y = "ENTREZID",
                   all.x = TRUE, sort = FALSE)
  
  ## 3) 必要列を整形
  out <- deg.top[ , c("EntrezID", "SYMBOL",
                      "logFC", "adj.P.Val", "P.Value")]
  colnames(out)[3:5] <- c("logFC", "FDR", "Pvalue")
  
  ## 4) コンソールへ全件表示
  old_opt <- getOption("max.print")          # 現行値を保存
  options(max.print = nrow(out) * ncol(out)) # 十分大きい値に拡張
  print(out, row.names = FALSE, digits = 3)
  options(max.print = old_opt)               # 元に戻す
}



# B-5. volcano plot -----------------------------------------------------------------
{
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  ## 0. limma‑treat の全テーブルを取得
  tbl_all <- topTreat(fit, p.value = 1, number = Inf)   # ← ここを 1 に
  
  ## 1. ラベル列を追加（先ほどと同じ）
  tbl_plot <- tbl_all %>%
    mutate(EntrezID   = rownames(.),
           GeneSymbol = AnnotationDbi::mapIds(org.Hs.eg.db,
                                              keys      = EntrezID,
                                              column    = "SYMBOL",
                                              keytype   = "ENTREZID",
                                              multiVals = "first"),
           GeneSymbol = ifelse(is.na(GeneSymbol), EntrezID, GeneSymbol),
           signif = case_when(
             adj.P.Val < 0.05 & logFC >=  0.30 ~ "Up",
             adj.P.Val < 0.05 & logFC <= -0.30 ~ "Down",
             TRUE                              ~ "NS"))
  
  ## 2. ラベルに使う上位 10 遺伝子だけ抽出（有意点から）
  top_labs <- tbl_plot %>%
    filter(signif != "NS") %>%
    slice_min(adj.P.Val, n = 10)
  
  ## 3. Volcano プロット
  ggplot(tbl_plot, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(colour = signif), size = 1.6, alpha = 0.8) +
    scale_colour_manual(values = c(Up="#D7191C", Down="#2C7BB6", NS="grey70")) +
    geom_vline(xintercept = c(-0.30, 0.30), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05),    linetype = "dashed", colour = "grey40") +
    geom_text_repel(data = top_labs,
                    aes(label = GeneSymbol), size = 3, max.overlaps = Inf) +
    labs(title = sprintf("Volcano plot | %d DEGs (|logFC| ≥ 0.30, FDR ≤ 0.05)",
                         sum(tbl_plot$signif != "NS")),
         x = "log2 fold-change", y = expression(-log[10]~FDR)) +
    theme_bw(base_size = 12) +
    theme(legend.title = element_blank())
}


# B-6. 相関プロット -----------------------------------------------------------------
{
  corr_deg <- cor(expr_sub)
  pheatmap(corr_deg, clustering_distance_rows="euclidean",
           main="Sample–sample correlation")
}

##### C: ----------------------------------------------------Network Analysis-----------------------------------------------------------
# C-1.  STRING PPI ── Off-line (3 files)-------------------------------
{
  suppressPackageStartupMessages({
  library(dplyr);  library(tibble)
  library(AnnotationDbi); library(org.Hs.eg.db)
  library(STRINGdb);      library(igraph); library(ggraph)
})

library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)   # AnnotationDbi の後でもう一度呼ぶと dplyr::select が上書き

cache_dir <- "C:/DMD_project"
# cache_dir <- "C:/Users/PC_User/OneDrive/大学院研究/統計解析用フォルダ"

needed_files <- c(
  "9606.protein.aliases.v11.5.txt.gz",
  "9606.protein.info.v11.5.txt.gz",
  "9606.protein.links.detailed.v11.5.txt.gz",
  "9606.protein.links.v11.5.txt.gz"       # ← 追加
)
missing <- needed_files[!file.exists(file.path(cache_dir, needed_files))]
if (length(missing) > 0)
  stop("❌ 以下のファイルがありません:\n  ",
       paste(missing, collapse = "\n  "),
       "\nブラウザでダウンロードし、cache_dir に置いてください。")
}

{
  str_db <- STRINGdb$new(version = "11.5",
                        species = 9606,
                        score_threshold = 500,
                        input_directory = cache_dir)
  
  ## 4. Entrez → SYMBOL 変換 & ケラチン系除外
  skip_pattern <- paste0(
    "^KRT","|^KRTAP","|^LCE","|^FLG$","|^SPRR","|^IVL$",
    "|^DSG","|^DSP$","|^TACSTD2$","|^OR\\d+","|^TAS2R\\d+")
  
  deg_tbl <- deg %>%                             # topTable の出力を想定
    rownames_to_column("EntrezID") %>%             
    mutate(Symbol = mapIds(org.Hs.eg.db,
                           keys      = EntrezID,
                           column    = "SYMBOL",
                           keytype   = "ENTREZID",
                           multiVals = "first")) %>% 
    filter(!grepl(skip_pattern, Symbol) & !is.na(Symbol)) %>% 
    dplyr::select(Symbol, logFC, adj.P.Val)         # ← ここを明示
  
  ## 5. STRING ID マッピング
  map_tbl <- str_db$map(deg_tbl, "Symbol", removeUnmappedRows = TRUE)
  cat(sprintf("🧬  STRING mapping: %d / %d genes hit\n",
              nrow(map_tbl), nrow(deg_tbl)))
  
  ## 6. ネットワーク構築 + 属性付与
  g_net <- str_db$get_subnetwork(map_tbl$STRING_id)
  
  V(g_net)$Symbol <- map_tbl$Symbol[match(V(g_net)$name, map_tbl$STRING_id)]
  V(g_net)$logFC  <- deg_tbl$logFC [match(V(g_net)$Symbol, deg_tbl$Symbol)]
  
  deg_vec <- igraph::degree(g_net)  # ← 衝突を回避
  k      <- min(30, igraph::vcount(g_net))
  hub_id <- names(sort(deg_vec, decreasing = TRUE))[seq_len(k)]
  V(g_net)$is_hub <- V(g_net)$name %in% hub_id
  
  ## 念のため：combined_score の保険（環境差対策）
  if (!"combined_score" %in% igraph::edge_attr_names(g_net) &&
      "weight" %in% igraph::edge_attr_names(g_net)) {
    E(g_net)$combined_score <- E(g_net)$weight
  }
}

## 7. 可視化
pal_fun <- colorRampPalette(c("navy","white","firebrick"))
ggraph(g_net, layout = "fr") +
  geom_edge_link(aes(width = combined_score / 1000),
                 colour = "grey70", alpha = 0.3) +
  scale_edge_width(range = c(0.2, 1.5),          # 幅調整
                   name  = "score") +  # ← ここで好きな文字列
  geom_node_point(aes(fill = logFC, size = is_hub),
                  shape = 21, colour = "black") +
  geom_node_text(aes(label = Symbol), size = 2,
                 repel = TRUE, check_overlap = TRUE) +
  scale_fill_gradientn(colours = pal_fun(100), name = "logFC") +
  scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2), guide = "none") +
  theme_void(base_size = 12) +
  labs(title = sprintf("STRING PPI  |  |logFC| ≥ 0.30  |  score ≥ 0.5  |  %d genes",
                       vcount(g_net)))



# connectionのある遺伝子のみラベル表示
{
  library(igraph)
  library(ggraph)
  library(ggrepel)
  
  ## 1) 連結成分のサイズ
  comp      <- igraph::components(g_net)
  comp_size <- comp$csize[comp$membership]
  
  ## 2) ラベル表示フラグ
  igraph::V(g_net)$show_lab <- comp_size >= 3
  igraph::V(g_net)$label    <- ifelse(igraph::V(g_net)$show_lab,
                                      igraph::V(g_net)$Symbol, "")
  
  ## 参考: ハブ判定は igraph::degree を明示
  hub_thr <- 3
  igraph::V(g_net)$is_hub <- igraph::degree(g_net) >= hub_thr
  
  ## 念のため: combined_score が無い環境への保険
  if (!"combined_score" %in% igraph::edge_attr_names(g_net)) {
    if ("weight" %in% igraph::edge_attr_names(g_net)) {
      igraph::E(g_net)$combined_score <- igraph::E(g_net)$weight
    } else {
      igraph::E(g_net)$combined_score <- 500
    }
  }
  
  pal_fun <- colorRampPalette(c("navy","white","firebrick"))
  
  ggraph(g_net, layout = "fr") +
    geom_edge_link(aes(width = combined_score/1000),
                   colour = "grey70", alpha = .3) +
    scale_edge_width(range = c(0.2, 1.5), name = "score") +
    geom_node_point(aes(fill = logFC, size = is_hub),
                    shape = 21, colour = "black") +
    scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2), guide = "none") +
    geom_node_text(
      data = function(d) d[d$show_lab, ],
      aes(label = label),
      repel = TRUE,
      max.overlaps = Inf,
      box.padding  = .4,
      point.padding = .3,
      segment.size = .2,
      segment.alpha = .5,
      size = 2.5
    ) +
    scale_fill_gradientn(colours = pal_fun(100), name = "logFC") +
    theme_void(base_size = 12) +
    labs(title = sprintf(
      "STRING PPI  |  |logFC| ≥ 0.30  |  score ≥ 0.5  |  %d genes",
      igraph::vcount(g_net)))
}




# connectionの無い遺伝子を全て除去
{
  library(igraph)
  
  ## 0) 次数ベクトル（明示的に igraph::degree）
  deg_vec  <- igraph::degree(g_net, mode = "all")
  iso_vtx  <- igraph::V(g_net)[deg_vec == 0]
  g_net_sub <- igraph::delete_vertices(g_net, iso_vtx)
  
  deg_vec_sub <- igraph::degree(g_net_sub, mode = "all")
  hub_thr <- 2
  is_hub  <- deg_vec_sub >= hub_thr
  
  igraph::V(g_net_sub)$is_hub <- is_hub
  igraph::V(g_net_sub)$label  <- ifelse(is_hub, igraph::V(g_net_sub)$Symbol, "")
  
  library(ggraph)
  pal_fun <- colorRampPalette(c("navy", "white", "firebrick"))
  
  ## 念のため: combined_score の保険
  if (!"combined_score" %in% igraph::edge_attr_names(g_net_sub)) {
    if ("weight" %in% igraph::edge_attr_names(g_net_sub)) {
      igraph::E(g_net_sub)$combined_score <- igraph::E(g_net_sub)$weight
    } else {
      igraph::E(g_net_sub)$combined_score <- 500
    }
  }
  
  ggraph(g_net_sub, layout = "fr") +
    geom_edge_link(aes(width = combined_score / 1000),
                   colour = "grey70", alpha = .3) +
    scale_edge_width(range = c(0.2, 1.5), name = "score") +
    geom_node_point(aes(fill = logFC, size = is_hub),
                    shape = 21, colour = "black") +
    scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2), guide = "none") +
    geom_node_text(
      data = function(d) d[d$is_hub, ],
      aes(label = label),
      repel = TRUE, size = 3.5,
      box.padding = .35, point.padding = .3,
      segment.size = .2, segment.alpha = .5
    ) +
    scale_fill_gradientn(colours = pal_fun(100), name = "logFC") +
    theme_void(base_size = 12) +
    labs(title = sprintf(
      "STRING PPI  |  |logFC| ≥ 0.30  |  score ≥ 0.5  |  connected genes: %d",
      igraph::vcount(g_net_sub)))
  
  
}



# ノード数
vcount(g_net)   # → 231

# エッジ数
ecount(g_net)   # → 328

# 1. DEG テーブルに存在するか
deg_tbl |> dplyr::filter(Symbol == "CIITA")

# 2. STRING へ変換されたか
map_tbl |> dplyr::filter(Symbol == "ACADS")

# ハブ相関チェック
corr <- cor(mat_hub_z)
pheatmap(corr, clustering_distance_rows="euclidean",
         main="Sample–sample correlation")


# CGB3, SERPINA3: 遺伝子名重複チェック
library(STRINGdb)
db <- STRINGdb$new(version="11.5", species=9606, score_threshold=500)
ids <- db$map(data.frame(gene="CGB3"), "gene", removeUnmappedRows = TRUE)
ids$STRING_id    # ← 複数 ID が返ることを確認
# [1] "9606.ENSP00000470813" "9606.ENSP00000301408" "9606.ENSP00000349954"

ids_2 <- db$map(data.frame(gene="SERPINA3"), "gene", removeUnmappedRows = TRUE)
ids_2$STRING_id    # ← 複数 ID が返ることを確認
# [1] "9606.ENSP00000376793" "9606.ENSP00000450540"
  

# C-2.  Hub gene heatmap ---------------------------------------------------------
{
  hub_ids    <- names(sort(igraph::degree(g_net), decreasing = TRUE))[1:102]
  hub_syms   <- V(g_net)$Symbol[match(hub_ids, V(g_net)$name)] |> na.omit() |> unique()
  hub_syms   <- intersect(hub_syms, rownames(expr_sym))           # expr に存在するものだけ
  
  mat_hub_z  <- t(scale(t(expr_sym[hub_syms, ])))
  mat_hub_z  <- mat_hub_z[rowSums(is.na(mat_hub_z)) == 0, ]
  
  ## 1. アノテーションデータフレームを拡張
  # 既存 ann_col は 1 列 (Group) でした
  ann_col <- data.frame(
    Group  = Group,                     # 既存
    Cohort = factor(Cohort,             # ★追加
                    levels = c("C1","C2","C3","C4"))
  )
  rownames(ann_col) <- colnames(mat_hub_z)
  
  
  ## 2. 色指定も 2 列分用意
  ann_cols <- list(
    Group  = c(Early = "#2C7BB6", Late = "#D7191C"),
    Cohort = c(C1 = "black",C2 = "grey", C3 = "#4DAF4A", C4 = "purple")   # ★C3 を緑で強調
  )
  
  
  ## 3. pheatmap 呼び出し
  pheatmap(
    mat_hub_z,
    scale              = "column",             # ★列方向で Z スケール
    color              = colorRampPalette(c("navy","white","firebrick3"))(50),
    annotation_col     = ann_col,
    annotation_colors  = ann_cols,
    cluster_rows       = TRUE,
    cluster_cols       = TRUE,
    show_rownames      = TRUE,
    fontsize_row       = 4,
    border_color    = TRUE,
    main = sprintf("STRING hub genes | Top %d genes by degree × %d samples",
                   nrow(mat_hub_z), ncol(mat_hub_z))
  )
}

## ④ 重複チェック
stopifnot(!any(duplicated(rownames(mat_hub_z))))
stopifnot(!any(duplicated(colnames(mat_hub_z))))

# C-3.  Hub-gene PCA -----------------------------------------------------------------
{
  pc <- prcomp(t(mat_hub_z))
  
  ## 1. aesthetic ベクトル
  col_vec <- c(Early = "#2C7BB6",  Late = "#D7191C")[Group]      # 青／赤
  pch_vec <- c(C1 = 17, C2 = 7, C3 = 5, C4 = 3)[Cohort]                # ■ ▲ □
  
  ## 2. プロット
  plot(pc$x[, 1:2],
       col = col_vec,
       pch = pch_vec,
       xlab = sprintf("PC1 (%.1f%%)", 100 * pc$sdev[1]^2 / sum(pc$sdev^2)),
       ylab = sprintf("PC2 (%.1f%%)", 100 * pc$sdev[2]^2 / sum(pc$sdev^2)))
  
  ## 3. 凡例
  legend("topleft",  legend = levels(Group),
         col = c("#2C7BB6", "#D7191C"), pch = 17, title = "Group", bty = "n")
  legend("topright", legend = levels(Cohort),
         col = "black",          pch = c(17, 7, 3, 5),
         title = "Cohort", bty = "n")
}

# C-4.  Hub-gene RLE -----------------------------------------------------------------
## ─ 1.  ライブラリ
if (!requireNamespace("EDASeq", quietly = TRUE))
  BiocManager::install("EDASeq", ask = FALSE, update = FALSE)
library(EDASeq)

## ─ 2.  RLE 値の計算
{
  rle_mat <- log(expr_sym)                           # ① すでに log2 ならこの行は不要
  rle_mat <- rle_mat - rowMedians(rle_mat, na.rm=TRUE)  # ② 各行の中央値を引く
  
  ## ─ 3.  可視化 (boxplot
  boxplot(as.data.frame(rle_mat),
          las = 2,                                 # 縦軸ラベル横向き
          outline = FALSE,
          col = c(C1="black","grey40","#4DAF4A","purple")[Cohort],  # Cohort 着色
          ylab = "RLE (log2, gene-median centred)",
          main = "Relative Log Expression (RLE)")
  
  legend("topright", legend = levels(Cohort),
         fill = c("black","grey40","#4DAF4A","purple"), bty="n")
}

##### D: --------------------------------------------------------ORA Analysis----------------------------------------------------------------
# D-1. KEGG / Reactome / GO‑BP-------------------------------------------------------
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(stringr)
  # enrichplot は clusterProfiler の依存として自動ロード済み
})


## 0. 前提チェック
{
  stopifnot(exists("combat"), exists("deg"))
  
  ## 1. 背景遺伝子 & 解析対象遺伝子
  bg_ids  <- unique(sub("_.*", "", rownames(combat)))    # 22k 強
  deg_ids <- rownames(deg)                             # 162 Entrez
  
  ## 2. KEGG ORA
  kegg_res <- enrichKEGG(gene         = deg_ids,
                         universe     = bg_ids,
                         organism     = "hsa",
                         keyType      = "ncbi-geneid",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.10)
  
  library(dplyr)
  
  ## ── 病名を含む KEGG 経路を取り除く ──
  kegg_clean <- kegg_res@result %>% 
    filter(!grepl("^hsa05|disease", ID))      # “05xxx” は疾病カテゴリー
  
  ## enrichplot 系で描画する場合は、結果を再ラップ
  kegg_clean <- new("enrichResult", kegg_res, result = kegg_clean)
  
  ## 3. Reactome ORA
  react_res <- enrichPathway(gene         = deg_ids,
                             universe     = bg_ids,
                             organism     = "human",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.10)
  
  ## 4. GO‑BP ORA（冗長削減）
  go_res <- enrichGO(gene          = deg_ids,
                     universe      = bg_ids,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.10) |>
    simplify(cutoff = 0.5, by = "p.adjust", select_fun = min)
  
  
  ## 5. 可視化用ユーティリティ
  pub_dotplot <- function(eres, n = 15, fdr_cut = 0.05,
                          title = "", wrap_len = 35) {
    
    df <- as.data.frame(eres@result)
    df <- df[df$p.adjust <= fdr_cut, ]
    if (!nrow(df)) {
      message("＜", title, "＞ 有意経路なし (FDR <", fdr_cut, ")")
      return(invisible(NULL))
    }
    df <- df[order(df$p.adjust), ][1:min(n, nrow(df)), ]
    df$GeneRatio <- sapply(strsplit(df$GeneRatio, "/"),
                           \(x) as.numeric(x[1])/as.numeric(x[2]))
    df$Description <- factor(str_wrap(df$Description, wrap_len),
                             levels = rev(str_wrap(df$Description, wrap_len)))
    
    ggplot(df, aes(x = GeneRatio, y = Description)) +
      geom_point(aes(size = Count, colour = p.adjust)) +
      scale_colour_gradient(low = "#b2182b", high = "#2166ac", name = "FDR") +
      scale_size(range = c(3, 8)) +
      labs(title = title, x = "Gene ratio", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = .5),
            axis.text.y = element_text(size = 9))
  }
  
  pub_barplot <- function(eres, n = 15, fdr_cut = 0.05,
                          title = "", wrap_len = 40) {
    
    df <- as.data.frame(eres@result)
    df <- df[df$p.adjust <= fdr_cut, ]
    if (!nrow(df)) {
      message("＜", title, "＞ 有意経路なし (FDR <", fdr_cut, ")")
      return(invisible(NULL))
    }
    df <- df[order(df$p.adjust), ][1:min(n, nrow(df)), ]
    df$GeneRatio <- sapply(strsplit(df$GeneRatio, "/"),
                           \(x) as.numeric(x[1])/as.numeric(x[2]))
    df$Description <- factor(str_wrap(df$Description, wrap_len),
                             levels = rev(str_wrap(df$Description, wrap_len)))
    
    ggplot(df, aes(x = GeneRatio, y = Description)) +
      geom_col(aes(fill = p.adjust), width = .7) +
      scale_fill_gradient(low = "#b2182b", high = "#2166ac", name = "FDR") +
      labs(title = title, x = "Gene ratio", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = .5),
            axis.text.y = element_text(size = 9),
            panel.grid.major.y = element_blank())
  }
  
}


# D-2. 図の出力例---------------------------------
pub_dotplot(react_res, n = 20, fdr_cut = 0.05,
            title = "Reactome enrichment")
pub_barplot(react_res, n = 20, fdr_cut = 0.10,
            title = "Reactome enrichment")
# 7yr,8yr では有意経路無し


pub_dotplot(go_res,    n = 20, fdr_cut = 0.05,
            title = "GO-BP enrichment")
pub_barplot(go_res,    n = 20, fdr_cut = 0.05,
            title = "GO-BP enrichment")


pub_dotplot(kegg_clean,  n = 20, fdr_cut = 0.05,
            title = "KEGG enrichment")
pub_barplot(kegg_clean,  n = 20, fdr_cut = 0.20,
            title = "KEGG enrichment")
# 8yr では有意経路無し


##### E: ----------------------------------------------------------------GSEA--------------------------------------------------------------------
# E-1. GSEA (Reactome) Analysis---------------------------------------------------------
suppressPackageStartupMessages({
  library(limma)
  library(clusterProfiler)
  library(BiocParallel)
  library(ggplot2)
  library(ggridges)
})

# 0. 解析パラメータ
{
  LFC_CUTOFF     <- 0.30     # treat() の閾値
  P_CUTOFF       <- 0.05     # gseaPathway の FDR 閾値
  SHOW_TOP       <- 20       # ridgeplot に表示する経路数
  X_BREAKS       <- seq(-2, 2, 1)
  
  # 1. モデルフィット & ランクベクトル
  designGLM <- model.matrix(~ 0 + Group)
  fit.lfc   <- v |>
    lmFit(designGLM) |>
    contrasts.fit(cbind(Late_vs_Early = c(-1, 1))) |>
    eBayes(trend = TRUE) |>
    treat(lfc = LFC_CUTOFF)
  
  # EntrezID を取り出して geneList を作成
  fc_vec <- setNames(fit.lfc$coefficients,
                     sub("_.*", "", rownames(fit.lfc$coefficients)))
  geneList <- sort(fc_vec, decreasing = TRUE)
  
  # 2. BiocParallel を安全に設定
  register(SerialParam())      # ← 安定運用。並列化するなら SnowParam() へ変更
  # 例: register(SnowParam(workers = 4, type = "SOCK", progressbar = TRUE))
  
  #3. GSEA 実行
  set.seed(1234)
  gsea_react <- gsePathway(
    geneList     = geneList,
    organism     = "human",
    pvalueCutoff = P_CUTOFF,   # 有意経路が少なければ 1 にして上位表示
    maxGSSize    = 3000,
    minGSSize    = 10,
    eps          = 0           # fgseaMultilevel を使用（高速・高精度）
  )
}


#4. 可視化関数
{
  RIDGE_LABEL_SZ <- 10        # y 軸ラベル文字サイズ
  ridgeplot2 <- function(gsea_obj,
                         showCategory      = 20,
                         pathway_text_size = 10,
                         x_breaks          = seq(-2, 2, 1))
  {
    enrichplot::ridgeplot(gsea_obj,
                          showCategory = showCategory,
                          fill         = "p.adjust") +
      scale_fill_gradient(low  = "#b2182b",
                          high = "#2166ac",
                          name = "FDR") +
      scale_x_continuous(name   = "Running enrichment score",
                         breaks = x_breaks,
                         labels = sprintf("%.0f", x_breaks)) +
      theme(axis.text.y  = element_text(size = pathway_text_size),
            axis.title.y = element_blank(),
            legend.position = "right",
            plot.title = element_text(hjust = 0.5))
  }
  
  # 5. Ridgeplot 描画
  {p_ridge <- ridgeplot2(
    gsea_obj          = gsea_react,
    showCategory      = SHOW_TOP,
    pathway_text_size = RIDGE_LABEL_SZ,   # ← y 軸ラベル（経路名）のサイズ
    x_breaks          = X_BREAKS
  ) +
      ggtitle(sprintf("Reactome enrichment (|log2FC| > %.2f)", LFC_CUTOFF)) +
      theme(
        plot.title  = element_text(hjust = 0.5, size = 14,  margin = margin(b = 6)),  # タイトル
        axis.text.x = element_text(size = 10),   # x 軸目盛
        axis.title.x= element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title= element_text(size = 10)
      )
    print(p_ridge)
  }
}


# E-2. ヒットした全経路をコンソール表示する-----------------------------------
{
  ## ❶ まず FDR<0.05 でフィルタリング
  tbl_short <- as.data.frame(gsea_react) |>
    dplyr::filter(p.adjust < 0.05) |>
    dplyr::select(Description, NES, p.adjust) |>
    dplyr::arrange(p.adjust)
  
  ## ❷ 数値列だけを有効桁 3 桁に丸めて表示
  tbl_print <- tbl_short |>
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),
        ~ signif(.x, 3)          # ← 数値列にだけ適用
      )
    )
  sig3 <- function(x, sig = 3){
    vapply(x, function(v){
      if (v == 0) return("0")
      ## 桁数 = sig - 1 - floor(log10(v))
      dec <- max(sig - 1 - floor(log10(abs(v))), 0)
      formatC(v, digits = dec, format = "f", drop0trailing = FALSE)
    }, FUN.VALUE = "")
  }
  tbl_print <- tbl_short %>%
    dplyr::mutate(p.adjust = sig3(p.adjust))
  print(tbl_print, row.names = FALSE)
  
}



# supplementary exploratory analysis with past research data--------------------

## 1) 線維面積率/筋線維面積率 vs MFD の非線形モデル
cor(MFD, MFA, method = "spearman")
qplot(MFD, MFA) +
  geom_smooth(span = 1.5, color = "red", se = F)

cor(MFD, CFA, method = "spearman")
qplot(MFD, CFA)+
  geom_smooth(span = 1, color = "red", se = F)

## 2) 部分スピアマン（年齢補正）
library(ppcor)
pcor.test(MFD, MFA, ABx, method = "spearman")
pcor.test(MFD, CFA, ABx, method = "spearman")

## 3) piecewise 回帰例
library(segmented)
lin_mod <- lm(MFA ~ MFD)
seg_mod <- segmented(lin_mod, seg.Z = ~ MFD, psi = list(MFD = 400))
summary(seg_mod)      # しきい値 (break‑point) と両側傾きを取得






# 1. DEG テーブルに存在するか (再掲) ---------------------------------------------------
deg_tbl |> dplyr::filter(Symbol == "CIITA")

# 2. STRING へ変換されたか
map_tbl |> dplyr::filter(Symbol == "Myogenin")





# 3. 免疫チェックポイント：興味遺伝子の検証-------------------------------------

# ① Entrez → Symbol の対応表を作る
id2symbol <- mapIds(org.Hs.eg.db,
                    keys       = rownames(fit$coefficients),
                    column     = "SYMBOL",
                    keytype    = "ENTREZID",
                    multiVals  = "first")

# ② limma の統計テーブルを取得して Symbol を付ける
tt <- topTable(fit,
               coef   = "Late_vs_Early",
               number = Inf,
               sort.by= "none") %>%
  rownames_to_column("ENTREZID") %>%
  mutate(Symbol = id2symbol[ENTREZID])     # 行ごとにシンボル付与


## 1. gsea_react@result から該当行を拾う
row_ifn <- gsea_react@result %>%
  filter(Description == "Assembly of collagen fibrils and other multimeric structures")

## 2. core_enrichment 列（EntrezID が "/" 区切り）を分割
lead_entrez <- str_split(row_ifn$core_enrichment, pattern = "/", simplify = TRUE)
lead_entrez <- lead_entrez[lead_entrez != ""]           # 空文字除去

## 3. Entrez → Symbol 変換
lead_symbol <- mapIds(org.Hs.eg.db,
                      keys      = lead_entrez,
                      column    = "SYMBOL",
                      keytype   = "ENTREZID",
                      multiVals = "first")

## 4. limma の統計テーブルから logFC と p 値を付ける
tt <- topTable(fit,
               coef = "Late_vs_Early",
               number = Inf, sort.by = "none") %>%
  rownames_to_column("ENTREZID") %>%
  mutate(Symbol = mapIds(org.Hs.eg.db,
                         keys = ENTREZID,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first"))

lead_tbl <- tt %>%
  filter(ENTREZID %in% lead_entrez) %>%
  select(Symbol, logFC, P.Value, adj.P.Val) %>%
  arrange(P.Value)

lead_tbl %>%
  tibble::as_tibble() %>%
  print(n = 100)
   # 先頭 20 行を確認
# A tibble: 12 × 4
# Symbol   logFC   adj.P.Val
# <chr>    <dbl>     <dbl>
# 1  HLA-DQA1 2.36   1.73e-10
# 2  HLA-DQB1 1.20   5.02e- 9
# 3  TRIM26   0.370  6.61e- 7
# 4  HLA-DRB1 1.15   2.02e- 5
# 5  HLA-DRB5 1.46   4.66e- 5
# 6  HLA-DRA  0.597  6.84e- 5
# 7  HLA-DPA1 0.745  7.26e- 5
# 8  HLA-DRB3 1.77   9.54e- 5
# 9  GBP3     0.492  5.55e- 4
# 10 HLA-E    0.338  5.97e- 4
# 11 MID1     0.434  1.54e- 3
# 12 TRIM34   0.369  4.73e- 3

ifn_idx <- which(gsea_react@result$Description == "Interferon gamma signaling")
gene_set <- gsea_react@result$core_enrichment[ifn_idx] %>%
  str_split("/", simplify = TRUE) %>% c()

camera_res <- camera(combat,
                     index  = list(ifn = gene_set),
                     design = model.matrix(~ Group)); camera_res
# NGenes Direction              PValue
# ifn     12        Up 0.00000000001333338


pd1_idx <- which(gsea_react@result$Description == "Co-inhibition by PD-1")
lead_pd1 <- str_split(gsea_react@result$core_enrichment[pd1_idx], "/",
                      simplify = TRUE)[1, ]
lead_pd1 <- lead_pd1[lead_pd1 != ""]

cat("\n### Co‑inhibition by PD‑1 — leading‑edge\n")
print(tt %>%
        filter(ENTREZID %in% lead_pd1) %>%
        select(Symbol, logFC, P.Value, adj.P.Val) %>%
        arrange(P.Value),
      row.names = FALSE)

# Symbol     logFC          adj.P.Val
# HLA-DQA1 2.3609192 0.0000000001729528
# HLA-DQB1 1.1974251  0.0000000050249192
# HLA-DRB1 1.1491410  0.0000201539947374
# HLA-DRB5 1.4614403  0.0000466135719200
# HLA-DRA  0.5967247  0.0000684066463499
# HLA-DPA1 0.7450210  0.0000726298343722
# HLA-DRB3 1.7741938  0.0000953766937888

camera_pd1 <- camera(combat,
                     index  = list(PD1 = lead_pd1),
                     design = model.matrix(~ Group))
camera_pd1$FDR <- p.adjust(camera_pd1$PValue, "BH")

cat("\n### Co‑inhibition by PD‑1 — camera\n")
print(camera_pd1)
# NGenes Direction            PValue
# PD1      7        Up 0.000000006274989

# PPARgammaについては
# Symbol    logFC     adj.P.Val
# 1 SCD    3.12          0.0000359
# 2 ADIPOQ 2.68          0.000531 
# 3 PLIN1  2.73          0.00434  
# 4 DGAT2  1.92          0.00459  
# 5 THRSP  1.03          0.0222   
# 6 CEBPA  0.824         0.0248   
# 7 FABP4  0.876         0.0357
#-------------------------------
# 8 PPARG  0.540         0.192    
# 9 H2BC1  0.647         0.557    
# 10 GPAM   0.633        0.759    
# 11 AGPAT2 0.491        0.778    
# 12 ACSS3  0.372         1        
# 13 H3C11  0.375         1        
# 14 H3-3B  0.334         1        
# 15 H4C8   0.330         1      


# NGenes Direction              PValue
# PPARG     15        Up 0.00000000002426011

## ------ 1. PD‑1 経路 ID を安全に取得
pd1_row <- gsea_react@result %>%
  filter(grepl("PD-1", Description, fixed = TRUE))

if (nrow(pd1_row) == 0) {
  stop("Reactome GSEA 結果に 'PD‑1' を含む経路が見つかりません。")
}

pd1_id <- pd1_row$ID[1]                       # 例: "R-HSA-389948"
pd1_entrez <- gsea_react@geneSets[[pd1_id]]

if (length(pd1_entrez) == 0) {
  stop("PD‑1 経路の遺伝子セットが空です。")
}

## ------ 2. limma 統計テーブルと突合
# 'tt' は既に存在 (ENTREZID 列つき). 無ければ再作成:
#   tt <- topTable(... ) |> rownames_to_column("ENTREZID") |> mutate(Symbol = id2sym[ENTREZID])

pd1_tbl <- tibble(ENTREZID = pd1_entrez) %>%
  left_join(tt, by = "ENTREZID") %>%       # 統計値付加
  mutate(Symbol = coalesce(Symbol,
                           mapIds(org.Hs.eg.db,
                                  keys      = ENTREZID,
                                  keytype   = "ENTREZID",
                                  column    = "SYMBOL",
                                  multiVals = "first"))) %>%
  select(Symbol, logFC, P.Value, adj.P.Val)

## ------ 3. 結果表示
print(pd1_tbl %>% arrange(P.Value, Symbol), n = Inf, row.names = FALSE)
# Symbol      logFC    P.Value     adj.P.Val
#   1 HLA-DQA1  2.36     6.06 e-13  0.0000000109
# 2 HLA-DQB1  1.20     2.35 e- 9  0.00000544  
# 3 HLA-DRB5  1.46     3.19 e- 5  0.00493     
# 4 HLA-DRB1  1.15     3.70 e- 5  0.00552     
# 5 HLA-DRB3  1.77     3.82 e- 5  0.00565     
# 6 HLA-DPA1  0.745    1.12 e- 3  0.0696      
# 7 HLA-DRA   0.597    4.57 e- 3  0.186       
# 8 CD3G      0.154    8.32 e- 1  1           
# 9 PDCD1LG2  0.210    8.52 e- 1  1           
# 10 CD3E      0.212    8.74 e- 1  1           
# 11 HLA-DRB4  0.0961   9.05 e- 1  1           
# 12 HLA-DPB1  0.0791   9.45 e- 1  1           
# 13 PTPN11   -0.136    9.54 e- 1  1           
# 14 HLA-DQA2  0.0158   9.83 e- 1  1           
# 15 CD247    -0.134    9.87 e- 1  1           
# 16 PDCD1    -0.0565   9.93 e- 1  1           
# 17 CD4       0.139    9.96 e- 1  1           
# 18 PTPN6     0.0789   9.99 e- 1  1           
# 19 CD274    -0.0281   9.99 e- 1  1           
# 20 CD3D      0.0187   1.000e+ 0  1           
# 21 LCK       0.0219   1.000e+ 0  1           
# 22 CSK       0.00280  1.000e+ 0  1


# バージョン情報----------------------------------------------------------------
##  ログ取得スクリプト  —— 主要ライブラリのバージョンを一覧表示

{
  pkgs <- c(
    # 解析コア
    "clusterProfiler", "ReactomePA", "org.Hs.eg.db", "reactome.db", "fgsea","STRINGdb",
    # 可視化・補助
    "enrichplot", "ggplot2", "ggridges", "stringr", "dplyr",
    # データ処理・依存
    "Biobase", "AnnotationDbi", "GOSemSim", "qsmooth", "impute"
  )
  
  ver <- sapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE))
      as.character(packageVersion(p)) else NA
  })
  
  version_df <- data.frame(Package = pkgs, Version = ver, row.names = NULL)
  
  print(version_df, right = FALSE, row.names = FALSE)
}

# 自宅 (2025.07.26)
# Package         Version
# clusterProfiler 4.17.0 
# ReactomePA      1.53.0 
# org.Hs.eg.db    3.21.0 
# reactome.db     1.92.0 
# fgsea           1.35.6 
# STRINGdb        2.21.0
# enrichplot      1.29.2 
# ggplot2         3.5.2  
# ggridges        0.5.6  
# stringr         1.5.1  
# dplyr           1.1.4  
# Biobase         2.69.0 
# AnnotationDbi   1.71.0 
# GOSemSim        2.35.0
# qsmooth         1.25.0
# impute          1.83.0 

sessionInfo()
# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26100)

# Matrix products: default
# LAPACK version 3.12.1

# locale:
#   [1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8   
# [3] LC_MONETARY=Japanese_Japan.utf8 LC_NUMERIC=C                   
# [5] LC_TIME=Japanese_Japan.utf8    

# time zone: Etc/GMT-9
# tzcode source: internal

# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] ggridges_0.5.6         enrichplot_1.29.2      XVector_0.49.0         ggrepel_0.9.6         
# [5] matrixStats_1.5.0      qsmooth_1.25.0         preprocessCore_1.71.2  hgu133b.db_3.13.0     
# [9] hgu133a.db_3.13.0      hgu133plus2.db_3.13.0  GEOquery_2.77.1        RColorBrewer_1.1-3    
# [13] digest_0.6.37          STRINGdb_2.21.0        ReactomePA_1.53.0      clusterProfiler_4.17.0
# [17] org.Hs.eg.db_3.21.0    AnnotationDbi_1.71.0   IRanges_2.43.0         S4Vectors_0.47.0      
# [21] Biobase_2.69.0         BiocGenerics_0.55.0    generics_0.1.4         impute_1.83.0         
# [25] sva_3.57.0             BiocParallel_1.43.4    genefilter_1.91.0      mgcv_1.9-3            
# [29] nlme_3.1-168           limma_3.65.1           ggraph_2.2.1           igraph_2.1.4          
# [33] pheatmap_1.0.13        lubridate_1.9.4        forcats_1.0.0          stringr_1.5.1         
# [37] dplyr_1.1.4            purrr_1.1.0            readr_2.1.5            tidyr_1.3.1           
# [41] tibble_3.3.0           ggplot2_3.5.2          tidyverse_2.0.0       

# loaded via a namespace (and not attached):
#   [1] fs_1.6.6                    bitops_1.0-9                httr_1.4.7                 
# [4] tools_4.5.1                 backports_1.5.0             R6_2.6.1                   
# [7] lazyeval_0.2.2              hgu95e.db_3.13.0            withr_3.0.2                
# [10] graphite_1.55.0             prettyunits_1.2.0           gridExtra_2.3              
# [13] cli_3.6.5                   labeling_0.4.3              Rsamtools_2.25.1           
# [16] yulab.utils_0.2.0           gson_0.1.0                  foreign_0.8-90             
# [19] DOSE_4.3.0                  R.utils_2.13.0              hgu95av2.db_3.13.0         
# [22] hgu95b.db_3.13.0            rentrez_1.2.4               plotrix_3.8-4              
# [25] rstudioapi_0.17.1           RSQLite_2.4.2               BiocIO_1.19.0              
# [28] gridGraphics_0.5-1          hwriter_1.3.2.1             gtools_3.9.5               
# [31] GO.db_3.21.0                interp_1.1-6                Matrix_1.7-3               
# [34] abind_1.4-8                 R.methodsS3_1.8.2           lifecycle_1.0.4            
# [37] yaml_2.3.10                 edgeR_4.7.3                 SummarizedExperiment_1.39.1
# [40] hgu95c.db_3.13.0            gplots_3.2.0                qvalue_2.41.0              
# [43] SparseArray_1.9.1           BiocFileCache_2.99.5        grid_4.5.1                 
# [46] blob_1.2.4                  pwalign_1.5.0               crayon_1.5.3               
# [49] ggtangle_0.0.7              lattice_0.22-7              cowplot_1.2.0              
# [52] GenomicFeatures_1.61.5      annotate_1.87.0             KEGGREST_1.49.1            
# [55] EDASeq_2.43.0               pillar_1.11.0               knitr_1.50                 
# [58] fgsea_1.35.6                GenomicRanges_1.61.1        rjson_0.2.23               
# [61] codetools_0.2-20            fastmatch_1.1-6             glue_1.8.0                 
# [64] ShortRead_1.67.0            ggfun_0.2.0                 data.table_1.17.8          
# [67] vctrs_0.6.5                 png_0.1-8                   treeio_1.33.0              
# [70] gtable_0.3.6                gsubfn_0.7                  cachem_1.1.0               
# [73] aroma.light_3.39.0          xfun_0.52                   S4Arrays_1.9.1             
# [76] tidygraph_1.3.1             Seqinfo_0.99.1              survival_3.8-3             
# [79] statmod_1.5.0               ggtree_3.17.1               bit64_4.6.0-1              
# [82] progress_1.2.3              filelock_1.0.3              GenomeInfoDb_1.45.8        
# [85] KernSmooth_2.23-26          rpart_4.1.24                colorspace_2.1-1           
# [88] DBI_1.2.3                   Hmisc_5.2-3                 nnet_7.3-20                
# [91] tidyselect_1.2.1            bit_4.6.0                   compiler_4.5.1             
# [94] curl_6.4.0                  chron_2.3-62                httr2_1.2.1                
# [97] graph_1.87.0                htmlTable_2.4.3             xml2_1.3.8                 
# [100] DelayedArray_0.35.2         rtracklayer_1.69.1          checkmate_2.3.2            
# [103] scales_1.4.0                caTools_1.18.3              rappdirs_0.3.3             
# [106] rmarkdown_2.29              jpeg_0.1-11                 htmltools_0.5.8.1          
# [109] pkgconfig_2.0.3             base64enc_0.1-3             MatrixGenerics_1.21.0      
# [112] dbplyr_2.5.0                fastmap_1.2.0               rlang_1.1.6                
# [115] htmlwidgets_1.6.4           UCSC.utils_1.5.0            farver_2.1.2               
# [118] jsonlite_2.0.0              GOSemSim_2.35.0             R.oo_1.27.1                
# [121] RCurl_1.98-1.17             magrittr_2.0.3              Formula_1.2-5              
# [124] ggplotify_0.1.2             patchwork_1.3.1             Rcpp_1.1.0                 
# [127] ape_5.8-1                   viridis_0.6.5               proto_1.0.0                
# [130] sqldf_0.4-11                stringi_1.8.7               MASS_7.3-65                
# [133] plyr_1.8.9                  parallel_4.5.1              deldir_2.0-4               
# [136] Biostrings_2.77.2           graphlayouts_1.2.2          splines_4.5.1              
# [139] hash_2.2.6.3                hms_1.1.3                   locfit_1.5-9.12            
# [142] reshape2_1.4.4              biomaRt_2.65.0              XML_3.99-0.18              
# [145] evaluate_1.0.4              latticeExtra_0.6-30         BiocManager_1.30.26        
# [148] tzdb_0.5.0                  tweenr_2.0.3                polyclip_1.10-7            
# [151] conflicted_1.2.0            ggforce_0.5.0               xtable_1.8-4               
# [154] restfulr_0.0.16             hgu95d.db_3.13.0            reactome.db_1.92.0         
# [157] tidytree_0.4.6              viridisLite_0.4.2           snow_0.4-4                 
# [160] aplot_0.2.8                 GenomicAlignments_1.45.1    memoise_2.0.1              
# [163] cluster_2.1.8.1             timechange_0.3.0 




## 保存：現在のワークスペースまるごと-------------------------------------------
getwd()         # 例: "C:/DMD_project"
# 保存先フォルダを指定（例：プロジェクト直下の /cache）
#dir.create("cache", showWarnings = FALSE)
#save.image(file = "cache/workspace_2025.0807.RData")


## 復元
load("cache/cache/workspace_2025.0807.RData")
# 読み込むと、その時点のオブジェクトが Environment に再現されます

deg_tbl |> dplyr::filter(Symbol == "NRF2")

# 2. STRING へ変換されたか
map_tbl |> dplyr::filter(Symbol == "Myogenin")