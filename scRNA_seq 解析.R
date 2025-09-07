## ============================================================
## 0) 環境セットアップ（R 4.5.1 で確認）
## ============================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("glmGamPoi","edgeR","scDblFinder"), ask = FALSE, update = FALSE)

pkgs <- c("Seurat","harmony","Matrix","lme4","remotes")
install.packages(setdiff(pkgs, rownames(installed.packages())), Ncpus = 2)

# Rcpp 系は先に
if (!"Rcpp" %in% rownames(installed.packages())) install.packages("Rcpp", type = "binary")
if (!"RcppArmadillo" %in% rownames(installed.packages())) install.packages("RcppArmadillo", type = "binary")

update.packages(ask = FALSE, checkBuilt = TRUE)

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(edgeR); library(lme4); library(Matrix)
  library(SingleCellExperiment); library(scDblFinder)
})

set.seed(123)

## ============================================================
## 1) ロギングとオプション
## ============================================================
data_root <- "C:/DMD_project/exprdata"
dir.create(data_root, showWarnings = FALSE, recursive = TRUE)

logfile <- file.path(data_root, "nichenet_run_mhc_high.log")
sink(logfile, split = TRUE)
.on.exit_sinks <- function(){ try(sink(NULL), silent = TRUE) }
on.exit(.on.exit_sinks(), add = TRUE)

opt_old_warn  <- getOption("warn")
opt_old_error <- getOption("error")
options(
  warn = 1,
  error = function(){
    cat("\n!! 未捕捉エラー発生。ログを確認してください。\n")
    traceback(2)
    .on.exit_sinks()
  }
)
message("[START] ", Sys.time())
options(bitmapType = "cairo")   # Windows の描画デバイス対策

## ============================================================
## 2) 10X 互換 mtx 読み込み（基本形は変更しない）
## ============================================================
read_one <- function(prefix, id){
  mtx   <- file.path(data_root, paste0(prefix, "_matrix.mtx.gz"))
  feat  <- file.path(data_root, paste0(prefix, "_features.tsv.gz"))
  bc    <- file.path(data_root, paste0(prefix, "_barcodes.tsv.gz"))
  m     <- ReadMtx(mtx, bc, feat)
  obj   <- CreateSeuratObject(m, project = id)
  obj$orig.ident <- id
  obj
}

objs <- list(
  wt1   = read_one("GSM6596509_wtNSG01",   "wt1"),
  wt2   = read_one("GSM6596510_wtNSG02",   "wt2"),
  mdx1  = read_one("GSM6596511_mdxNSG01",  "mdx1"),
  mdx2  = read_one("GSM6596512_mdxNSG02",  "mdx2"),
  d2_1  = read_one("GSM6596513_mdxD2NSG01","d2_1"),
  d2_2  = read_one("GSM6596514_mdxD2NSG02","d2_2")
)

mhc_core <- c("H2-Aa","H2-Ab1","H2-Eb1","Ciita")

## ============================================================
## 3) QC（ミト％・UMI・検出遺伝子数）
## ============================================================
objs <- lapply(objs, function(o){
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^mt-")
  subset(
    o,
    subset = nFeature_RNA > 200 &
      nFeature_RNA < 8000 &
      nCount_RNA   < 20000 &
      percent.mt   < 10
  )
})

## ============================================================
## 3.5) ダブルレット除去（scDblFinder; 各サンプルごと）
## ============================================================
## --- ダブルレット除去の前後で細胞数を計測し、除外数を算出 ---
# 除去前の細胞数（QC後）
pre_n  <- sapply(objs, function(o) length(Cells(o)))   # 除去“前”の細胞数を保存

# ダブルレット除去（既存コード）
objs <- lapply(objs, run_scDblFinder)

post_n <- sapply(objs, function(o) length(Cells(o)))   # 除去“後”の細胞数
rm_summ <- data.frame(
  sample        = names(pre_n),
  n_before      = as.integer(pre_n),
  n_after       = as.integer(post_n),
  n_removed     = as.integer(pre_n - post_n),
  frac_removed  = round(100 * (pre_n - post_n) / pmax(pre_n, 1), 2)  # %
)
print(rm_summ)
#      sample n_before n_after n_removed frac_removed
# wt1     wt1     3200    3067       133         4.16
# wt2     wt2     1807    1730        77         4.26
# mdx1   mdx1     5229    4869       360         6.88
# mdx2   mdx2     1746    1664        82         4.70
# d2_1   d2_1     7326    6735       591         8.07
# d2_2   d2_2     4959    4613       346         6.98

## ============================================================
## 4) SCTransform（各サンプル）と軽量統合（rpca）
## ============================================================
objs <- lapply(
  objs,
  SCTransform,
  vars.to.regress = "percent.mt",
  method = if ("glmGamPoi" %in% rownames(installed.packages())) "glmGamPoi" else "poisson",
  verbose = TRUE
)

var.features <- SelectIntegrationFeatures(objs, nfeatures = 5000)

# 追加遺伝子は「全サンプル共通」に限定（←ここが重要）
common_all   <- Reduce(intersect, lapply(objs, function(o) rownames(o)))
features_use <- unique(c(var.features, intersect(mhc_core, common_all)))

# 念のための整合性チェック
stopifnot(all(sapply(objs, function(o) all(features_use %in% rownames(o)))))

objs <- PrepSCTIntegration(objs, anchor.features = features_use)

# rpca 前処理
objs <- lapply(objs, RunPCA, features = features_use, npcs = 30, verbose = FALSE)

anchors_sc <- FindIntegrationAnchors(
  object.list          = objs,
  normalization.method = "SCT",
  anchor.features      = features_use,
  reduction            = "rpca",
  dims                 = 1:30
)

muscle_sc <- IntegrateData(
  anchorset            = anchors_sc,
  normalization.method = "SCT",
  dims                 = 1:30,
  features.to.integrate= features_use
)

rm(objs, anchors_sc); gc()
#              used    (Mb)  gc trigger     (Mb)    max used     (Mb)
# Ncells   26699446  1426.0    43092531   2301.4    43092531   2301.4
# Vcells 8817841601 67274.8 17284033868 131866.8 14133710207 107831.7
## ============================================================
## 5) 次元削減・クラスタリング
## ============================================================
DefaultAssay(muscle_sc) <- "integrated"
muscle_sc <- RunPCA(muscle_sc, npcs = 50, verbose = TRUE)
muscle_sc <- RunUMAP(muscle_sc, dims = 1:30)
muscle_sc <- FindNeighbors(muscle_sc, dims = 1:30)
muscle_sc <- FindClusters(muscle_sc, resolution = 0.4)
# Maximum modularity in 10 random starts: 0.9466
# Number of communities: 23

## 表示用メタ情報
muscle_sc$geno <- factor(sub("^(wt|mdx|d2).*","\\1", muscle_sc$orig.ident),
                         levels = c("wt","mdx","d2"))

## まだアノテーションがない場合はクラスタ番号を仮の celltype に
if (!"celltype" %in% colnames(muscle_sc@meta.data)) {
  muscle_sc$celltype <- Idents(muscle_sc)
}

## 保存
#rds_dir <- "C:/DMD_project/rds"
#dir.create(rds_dir, showWarnings = FALSE)
#saveRDS(muscle_sc, file.path(rds_dir, "muscle_sc_proto.rds"), compress = "gzip")

muscle_proto <- readRDS("C:/DMD_project/rds/muscle_sc_proto.rds")


## 最低限の描画（デバイス不具合回避のため print のみ）
print(DimPlot(muscle_proto, group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(DimPlot(muscle_proto, group.by = "orig.ident"))
print(DimPlot(muscle_proto, group.by = "celltype", split.by = "geno"))

## --- WT 中のクラスタ22の個数を数える（最短） ---
meta <- muscle_proto@meta.data
geno_col <- if ("geno" %in% colnames(meta)) "geno" else "genotype"
wt_label <- if ("wt" %in% meta[[geno_col]]) "wt" else "WT"

n_22_WT <- sum(as.character(meta$seurat_clusters) == "22" & meta[[geno_col]] == wt_label)
n_22_WT
# [1] 0

table( meta[[geno_col]][ as.character(meta$seurat_clusters) == "22" ] )
#  wt mdx  d2 
#   0   7  86


muscle_origin <- muscle_proto

# 読み込み（割り当て検索用）
#muscle_sc <- readRDS("C:/DMD_project/rds/muscle_sc_proto.rds")
#muscle_origin <- readRDS("C:/DMD_project/rds/muscle_sc_proto.rds")
muscle_origin$cluster <- Idents(muscle_origin) 
DimPlot(muscle_origin, reduction = "umap", label = TRUE)

muscle_sc$cluster <- Idents(muscle_sc)        # 後で参照する用
DimPlot(muscle_sc, reduction = "umap", label = TRUE)

## --- バッチ混合の指標として iLISI を確認（目安 2 以上)
# PCA で iLISI（30次元を使用）
X <- Embeddings(muscle_sc, "pca")[, 1:30, drop = FALSE]
ilisi <- compute_lisi(X, muscle_sc@meta.data, label_colnames = "orig.ident")
summary(ilisi$orig.ident)  # 2.447


## --- 前提: Seurat オブジェクト muscle_sc がメモリにある ---
stopifnot("Seurat" %in% class(muscle_origin))
stopifnot("Seurat" %in% class(muscle_sc))

## 1) アッセイとUMAPの存在を整える ---------------------------------
assay_names <- names(Assays(muscle_origin))
assay_names <- names(Assays(muscle_origin))

DefaultAssay(muscle_origin) <- if ("SCT" %in% assay_names) "SCT" else if ("RNA" %in% assay_names) "RNA" else assay_names[1]
DefaultAssay(muscle_origin) <- if ("SCT" %in% assay_names) "SCT" else if ("RNA" %in% assay_names) "RNA" else assay_names[1]

# まずはSCTを既定に
DefaultAssay(muscle_origin) <- "SCT"
DefaultAssay(muscle_origin) <- "SCT"

# SCT内で使えるレイヤーを自動選択（data が無ければ counts）
layer_use <- intersect(c("data","counts","scale.data"), Layers(muscle_origin[["SCT"]]))[1]
layer_use <- intersect(c("data","counts","scale.data"), Layers(muscle_origin[["SCT"]]))[1]

if (!"umap" %in% names(muscle_origin@reductions)) {
  suppressMessages({
    muscle_sc <- RunPCA(muscle_origin, npcs = 50, verbose = FALSE)
    muscle_sc <- RunUMAP(muscle_origin, dims = 1:30, verbose = FALSE)
  })
}

if (!"umap" %in% names(muscle_origin@reductions)) {
  suppressMessages({
    muscle_origin <- RunPCA(muscle_origin, npcs = 50, verbose = FALSE)
    muscle_origin <- RunUMAP(muscle_origin, dims = 1:30, verbose = FALSE)
  })
}

## 2) ユーティリティ: 遺伝子の存在確認とブレンド描画 ---------------
check_genes <- function(obj, genes) {
  miss <- setdiff(genes, rownames(obj))
  if (length(miss)) {
    warning(sprintf("オブジェクトに存在しない遺伝子: %s",
                    paste(miss, collapse = ", ")))
  }
  invisible(setdiff(genes, miss))
}

blend2 <- function(obj, g1, g2, title = NULL) {
  gg <- FeaturePlot(
    obj, features = c(g1, g2),
    blend = TRUE, order = TRUE, max.cutoff = "q95"
  )
  if (!is.null(title)) gg <- gg + ggtitle(title)
  gg
}

## 3) クラスタ図（基準用） ------------------------------------------
p_clusters <- DimPlot(muscle_origin, reduction = "umap", label = TRUE) +
  ggtitle("Clusters")

## 4) 目的集団の2遺伝子ブレンド可視化 -------------------------------
# Schwann
check_genes(muscle_origin, c("Egr2","Ngfr"))
p_schwann <- blend2(muscle_origin, "Egr2","Ngfr",  "Schwann")


# Neutrophil
check_genes(muscle_origin, c("S100a9","Cxcr2"))
p_neutro  <- blend2(muscle_origin, "S100a9", "Cxcr2",   "Neutrophil")

# Mast
check_genes(muscle_origin, c("Cma1","Hpgds"))
p_mast  <- blend2(muscle_origin, "Cma1",   "Hpgds",    "Mast")
p_mast

# Dendritic
check_genes(muscle_origin, c("Itgax","Flt3",))
p_dc      <- blend2(muscle_origin, "Itgax",   "Flt3",    "DC")
p_dc

p_dc2      <- blend2(muscle_origin, "Cd74",   "Itgax",    "DC")
p_dc2

# Immuno-stromal
ImStromal1 <-  blend2(muscle_origin, "Ly6a","Pdgfra","IS1")
ImStromal <-  blend2(muscle_origin, "Lum","Cd53","Immune-Stromal")
ImStromal1/ImStromal

# Platelet
check_genes(muscle_origin, c("Pf4","Plek"))
p_Plt      <- blend2(muscle_origin,
                     "Pf4",
                     "Plek",
                     "Platelet")
p_Plt2      <- blend2(muscle_origin,
                     "Itga2b",
                     "Nbeal2",
                     "Platelet")
p_Plt_3 <- blend2(muscle_origin,
                      "Cd74",
                      "Itgax",
                      "Platelet_negative")
p_Plt_4 <- blend2(muscle_origin,
                  "Csf1r",
                  "Ly6c2",
                  "Platelet_negative")
p_Plt/p_Plt2/p_Plt_3/p_Plt_4

# RBC
check_genes(muscle_origin, c("Hbb-bt","Hba-a1"))
p_RBC      <- blend2(muscle_origin, "Hbb-bt",   "Hba-a1",    "Erythrocyte")
p_RBC

# Myocyte
check_genes(muscle_origin, c("Myl1","Tnni2"))
p_myocyte      <- blend2(muscle_origin, "Myl1",   "Tnni2",    "Myocyte")

# Satellite
check_genes(muscle_origin, c("Pax7","Myf5"))
p_St      <- blend2(muscle_origin, "Pax7",   "Myf5",    "Satellite")

# Pericyte
check_genes(muscle_origin, c("Rgs5","Abcc9"))
p_Pericy      <- blend2(muscle_origin, "Rgs5",   "Abcc9",    "Pericyte")

# Tendon
check_genes(muscle_origin, c("Scx","Tnmd"))
p_Tendo      <- blend2(muscle_origin, "Scx",   "Tnmd",    "Tendon")

# EC
check_genes(muscle_origin, c("Pecam1","Cdh5"))
p_EC      <- blend2(muscle_origin, "Pecam1",   "Cdh5",    "EC")

# Macro
check_genes(muscle_origin, c("Ptprc","Cd68"))
p_Macro      <- blend2(muscle_origin, "Ptprc",   "Cd68",    "Macro")
p_Macro

# MHC-II
check_genes(muscle_origin, c("H2-Aa","H2-Ab1"))
p_MHCII      <- blend2(muscle_origin, "H2-Aa",   "H2-Ab1",    "MHC-II")

# Stromal_1
check_genes(muscle_origin, c("Pdgfra","Meox1"))
p_Stromal_1      <- blend2(muscle_origin, "Pdgfra",   "Meox1",    "Stromal_1")

# Stromal_2
check_genes(muscle_origin, c("Pi16","Postn"))
p_Stromal_2      <- blend2(muscle_origin, "Pi16",   "Postn",    "Stromal_2")

# Stromal_3
check_genes(muscle_origin, c("Vcan","Ly6a"))
p_Stromal_3      <- blend2(muscle_origin, "Vcan","Ly6a", "Stromal_3")


## 5) まとめて表示（patchwork） ---------------------------------------
library(patchwork) 

(p_Macro | p_MHCII)  / (p_dc | ImStromal) / (p_schwann | p_neutro) / (p_Plt | p_RBC)
(p_St | p_myocyte) / (p_Pericy | p_Tendo) / (p_EC | p_Stromal_1) / (p_Stromal_2 | p_Stromal_3)

p_schwann # 21
p_mast # 9の一部
p_dc # 17,22
p_Macro
p_MHCII
p_neutro # 16
p_Plt
p_RBC # 12
p_St # 19
p_Pericy # 20
p_Tendo # 13
p_EC # 3
p_FAP
p_Stromal_1
p_Stromal_2
p_Stromal_3
p_myocyte # 11
p_Macro / p_MHCII
# MHC-high macrophage: 2,17
# Immune-stromal(CD53+): 6
# MHC-low macrophage: 0,4,8
# Schwann cell: 21
# Satellite cell: 19
# Myocyte/Myonucleus: 11
# Erythrocyte: 12
# Neutrophil: 16
# Endothelial cell(EC): 3
# Pericyte/Vascular SMC: 20
# Tenocyte/Tendon fibroblast: 13
# Dendritic 22
# Stromal: 1,5,7,10,14,15
# Mast cell: 9
# Platelet: 18


muscle_dotplot <- muscle_origin

Assays(muscle_sc)            # 例: "RNA" "SCT" "integrated"
Assays(muscle_dotplot)            # 例: "RNA" "SCT" "integrated"

DefaultAssay(muscle_sc)      # 何になっているか
DefaultAssay(muscle_dotplot)      # 何になっているか

Layers(muscle_sc[["SCT"]])   # SCT assay 内の layer 一覧
Layers(muscle_dotplot[["SCT"]])   # SCT assay 内の layer 一覧

Layers(muscle_sc[["RNA"]])   # RNA assay (あれば)
Layers(muscle_dotplot[["RNA"]])   # RNA assay (あれば)

# 適切な assay に合わせて
DefaultAssay(muscle_sc) <- "SCT"
DefaultAssay(muscle_dotplot) <- "SCT"

## 2-2 ラベリング  細胞型ラベル付与--------------------------------------------------
cluster_cell <- setNames(
  c("MHC-low macrophage", # 0 low
    "Stromal", # 1
    "MHC-high macrophage", # 2 high
    "Endothelial cell(EC)", # 3
    "MHC-low macrophage", # 4 low
    "Stromal", # 5
    "Immune-Stromal(Cd53+)", # 6 low
    "Stromal", # 7
    "MHC-low macrophage", # 8 low
    "Immune-Stromal(Cd53+)", # 9
    "Stromal", # 10
    "Myocyte/Myonucleus", # 11
    "Erythrocyte", # 12
    "Tenocyte/Tendon fibroblast", # 13
    "Stromal", # 14
    "Stromal",  # 15
    "Neutrophil",  # 16
    "Dendritic cell(DC)", # 17 high
    "MHC-low macrophage", # 18
    "Satellite cell", # 19
    "Pericyte/Vascular SMC", # 20
    "Schwann cell", # 21
    "Dendritic cell(DC)"),# 22
  0:22
)


Idents(muscle_sc)   <- "orig_cluster"
Idents(muscle_dotplot)   <- "orig_cluster"

Idents(muscle_sc) <- "seurat_clusters"   # 念押し
Idents(muscle_dotplot) <- "seurat_clusters"   # 念押し

muscle_sc <- RenameIdents(muscle_sc, cluster_cell)
muscle_dotplot <- RenameIdents(muscle_dotplot, cluster_cell)

muscle_sc$celltype <- Idents(muscle_sc)  # 固定
muscle_dotplot$celltype <- Idents(muscle_dotplot)  # 固定

DefaultAssay(muscle_sc) <- "RNA"   # 表示は生カウント基準
DefaultAssay(muscle_dotplot) <- "RNA"   # 表示は生カウント基準

DimPlot(muscle_dotplot, reduction = "umap", label = TRUE) +
  ggtitle("UMAP with assigned cell types")


# 10X ライブラリ名 → genotype へ写像
geno_map <- c(
  wt1  = "WT",   wt2  = "WT",
  mdx1 = "mdx",  mdx2 = "mdx",
  d2_1 = "mdxD2", d2_2 = "mdxD2"
)

# orig.ident が入っていることを確認
table(muscle_sc$orig.ident)
table(muscle_dotplot$orig.ident)

# メタデータに追加（上書き可）
muscle_sc$genotype <- base::unname(geno_map[muscle_sc$orig.ident])
muscle_dotplot$genotype <- base::unname(geno_map[muscle_dotplot$orig.ident])

# 因子にして順序を固定（任意）
muscle_sc$genotype <- factor(muscle_sc$genotype,
                             levels = c("WT","mdx","mdxD2"))
muscle_dotplot$genotype <- factor(muscle_dotplot$genotype,
                                   levels = c("WT","mdx","mdxD2"))

Idents(muscle_sc) <- "celltype"   # ← 念のため
Idents(muscle_dotplot) <- "celltype"   # ← 念のため

library(ggplot2)
dot_obj <- DotPlot(
  muscle_dotplot,
  features = list(
    "MHC-II core" = c("H2-Aa","H2-Ab1","H2-Eb1", "Ciita"),
    "Macrophage"  = c("Ptprc","Cd53","Cd68","Ctss"),
    "DC" = c("Itgax","Flt3"),
    "Neutro" = c("Lcn2","Cxcr2"),
    "Stromal" = c("Pdgfra","Ly6a"),
    "FAP" = c("Col15a1","Dpp4","Pi16","Postn"),
    "Fat/Adipogenesis" = c("Ptgs2","Cebpb","Scd1","Fabp3","Fgfr4"),
    "ECM" = c("Fn1","Lum","Thbs2"),
    "Fibroblast" = c("Fbln1","Fbln2","Vcan"),
    #"IFNg" = c("Irf1","Stat1"),
    #"Co-stim" = c("Cd274","Cd86","Cd38","Tnfsf9","Cd83"),
    #"FIM" = c("Lpl","C1qa","C1qb","C1qc"),
    "Myocyte"  = c("Tnni2","Des"),
    "Satellite"  = c("Pax7","Myf5"),
    "EC" = c("Pecam1","Cdh5"),
    "Pericyte" = c("Rgs5","Abcc9"),
    "Tendon" = c("Scx","Tnmd")
  ),
  group.by  = "celltype",
  split.by  = "genotype",
  assay     = "RNA",
  cols      = c("lightgrey","mediumpurple","darkorchid4"),
  dot.scale = 5,            # ドット径の上限
  dot.min   = 0.005          # 1% 未満は非表示
) +
  RotatedAxis() +
  theme(                                   # ← ここで文字サイズなどを調整
    axis.text.x  = element_text(size = 8), # 遺伝子名
    axis.text.y  = element_text(size = 8), # 行ラベル
    strip.text.x = element_text(size = 10)  # “MHC‑II core” などの列見出し
  ); print(dot_obj)


## 保存
#rds_dir <- "C:/DMD_project/rds"
#dir.create(rds_dir, showWarnings = FALSE)
#saveRDS(muscle_sc, file.path(rds_dir, "muscle_sc_after_dotplot.rds"), compress = "gzip")

#muscle_sc <- readRDS("C:/DMD_project/rds/muscle_sc_after_dotplot.rds")





## ==== DMD scRNA-seq: Doublet（layers単位）→ RBC+doublet 除外 → 擬似バルク → GLM / 共変量調整付きMHC-II検定 ====
suppressPackageStartupMessages({
  library(Seurat)            # v5 前提
  library(SeuratObject)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(edgeR)
  library(limma)
})

set.seed(123)
options(stringsAsFactors = FALSE)

## 0) 入力 ---------------------------------------------------------------
##  事前に Seurat v5 オブジェクト（RNA assay / counts.* layersあり）を用意
muscle_sc <- readRDS("C:/DMD_project/rds/muscle_sc_after_dotplot.rds")
stopifnot(inherits(muscle_sc, "Seurat"))
stopifnot("RNA" %in% names(muscle_sc@assays))
DefaultAssay(muscle_sc) <- "RNA"

## 1) counts.* レイヤー列挙（JoinLayersしない）--------------------------
layers_all <- Layers(muscle_sc[["RNA"]])
layers <- grep("^counts(\\.|$)", layers_all, value = TRUE)
layers <- setdiff(layers, "counts")   # ベース"counts"は除く（生サンプルごとの層に限定）
stopifnot(length(layers) > 0)

## 2) 出力ベクトル（グローバル名で初期化）
all_cells <- colnames(muscle_sc)
dbl_class <- setNames(rep(NA_character_, length(all_cells)), all_cells)
dbl_score <- setNames(rep(NA_real_,      length(all_cells)), all_cells)

## 3) 各レイヤーで scDblFinder → Cells(..., layer=) で厳密対応 ----------
summary_C <- data.frame(layer=character(), cells=integer(), dbl=integer(), ratio=numeric())
for (lay in layers) {
  mat <- muscle_sc[["RNA"]]@layers[[lay]]
  if (is.null(mat) || ncol(mat) == 0) {
    message(sprintf("[skip] %s: empty", lay)); next
  }
  cells_lay <- Cells(muscle_sc, layer = lay)  # グローバル列名（matの列順と一致）
  stopifnot(length(cells_lay) == ncol(mat))
  
  sce <- SingleCellExperiment(list(counts = mat))
  set.seed(123)
  sce <- scDblFinder(sce)  # 既定でサンプル内の局所密度・ライブラリサイズなどを用いて推定
  
  calls  <- as.character(colData(sce)$scDblFinder.class); names(calls)  <- cells_lay
  scores <- as.numeric(  colData(sce)$scDblFinder.score ); names(scores) <- cells_lay
  
  dbl_class[cells_lay] <- calls[cells_lay]
  dbl_score[cells_lay] <- scores[cells_lay]
  
  summary_C <- rbind(summary_C, data.frame(
    layer = lay, cells = length(cells_lay),
    dbl   = sum(calls == "doublet"),
    ratio = mean(calls == "doublet")
  ))
}

## 4) メタデータへ反映（未割当は singlet として保守的に補完）-----------
dbl_class[is.na(dbl_class)] <- "singlet"
muscle_sc$scDblFinder.class <- dbl_class
muscle_sc$scDblFinder.score <- dbl_score

cat(sprintf("\n[scDblFinder] Overall doublet = %.2f%%\n",
            100 * mean(muscle_sc$scDblFinder.class == "doublet", na.rm = TRUE)))
# [scDblFinder] Overall doublet = 3.21%

print(summary_C)
#         layer cells dbl      ratio
# 1  counts.wt1  3067  53 0.01728073
# 2  counts.wt2  1730  49 0.02832370
# 3 counts.mdx1  4869 164 0.03368248
# 4 counts.mdx2  1664  70 0.04206731
# 5 counts.d2_1  6735 176 0.02613215
# 6 counts.d2_2  4613 216 0.04682419

print(addmargins(table(muscle_sc$orig.ident, muscle_sc$scDblFinder.class, useNA = "ifany")))

## 5) RBC と doublet を同時除外 -----------------------------------------
stopifnot("celltype" %in% colnames(muscle_sc@meta.data))
n0   <- ncol(muscle_sc)
nRBC <- sum(muscle_sc$celltype == "Erythrocyte", na.rm = TRUE)
nDBL <- sum(muscle_sc$scDblFinder.class == "doublet", na.rm = TRUE)
message(sprintf("Before filter: %d cells (Erythrocyte=%d, doublet=%d)", n0, nRBC, nDBL))
# Before filter: 22678 cells (Erythrocyte=672, doublet=728)

muscle_sc <- subset(
  muscle_sc,
  subset = (celltype != "Erythrocyte") & (scDblFinder.class != "doublet")
)
muscle_sc$celltype <- droplevels(muscle_sc$celltype)
Idents(muscle_sc)  <- "celltype"

n1 <- ncol(muscle_sc)
message(sprintf("After Erythrocyte+doublet removal: %d cells kept (removed %d)", n1, n0 - n1))
# After Erythrocyte+doublet removal: 21294 cells kept (removed 1384)

print(addmargins(table(muscle_sc$celltype, muscle_sc$orig.ident)))
#                               d2_1  d2_2  mdx1  mdx2   wt1   wt2   Sum
# MHC-low macrophage          2490  1792  1327    83   290    98  6080
# Stromal                     1048  1196  1632   588  1409   922  6795
# MHC-high macrophage          555   390   607   122   133    50  1857
# Endothelial cell(EC)         245   129   304   190   574   244  1686
# Immune-Stromal(Cd53+)       1376   183   120   268    26     4  1977
# Myocyte/Myonucleus            37    28   161    89   168   225   708
# Tenocyte/Tendon fibroblast    71   208   165    38    84    36   602
# Neutrophil                   144   139   137     8    50     4   482
# Dendritic cell(DC)           131   210   106    47    35    12   541
# Satellite cell                23     9    86    57    76    21   272
# Pericyte/Vascular SMC         14    34    34    24    48    29   183
# Schwann cell                  16     8    10    16    41    20   111
# Sum                         6150  4326  4689  1530  2934  1665 21294

## 6) 擬似バルク（celltype × orig.ident）---------------------------------
##   counts合算（RNA/slot=counts）を取得。結果は行=遺伝子、列=celltype_sample
agg <- AggregateExpression(
  muscle_sc,
  group.by = c("celltype","orig.ident"),
  assays   = "RNA",
  slot     = "counts"
)
pseudo_sc <- agg$RNA
cat("pseudo_sc dim = ", paste(dim(pseudo_sc), collapse=" x "), "\n")
# pseudo_sc dim =  31053 x 72

stopifnot(is.matrix(pseudo_sc) || inherits(pseudo_sc, "dgCMatrix"))

## 7) 関心遺伝子パネル ---------------------------------------------------
gene_target_sc <- unique(intersect(
  c("H2-Aa","H2-Ab1","H2-Eb1","Ciita",      # MHC-II core
    "Ptprc","Cd68","Mrc1","Ctss","Ccl2",    # Mac
    "Col15a1","Dpp4","Pi16","Sfrp2","Postn",# FAP/ECM
    "Ptgs2","Cebpb","Scd1","Fabp3","Fgfr4",
    "Fn1","Lum","Thbs2","Fbln1","Fbln2","Vcan",
    "Irf1","Stat1","Cd274","Cd86","Cd38","Tnfsf9","Cd83",
    "Lpl","C1qa","C1qb","C1qc",
    "Myog","Des","Pax7","Myf5",
    "Pecam1","Cdh5","Pdgfra","Meox1","Ly6a",
    "Rgs5","Abcc9","Scx","Tnmd"),
  rownames(pseudo_sc)
))
stopifnot(length(gene_target_sc) > 0)

## 8) 細胞型ごとの edgeR GLM ---------------------------------------------
esc <- function(s) gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", s))  # 正規表現エスケープ

edgeR_by_celltype <- function(ct) {
  cols <- colnames(pseudo_sc)[startsWith(colnames(pseudo_sc), paste0(ct, "_"))]
  if (!length(cols)) return(NULL)
  
  labs <- sub(".*_", "", cols)
  grp  <- factor(sub("^(wt|mdx|d2).*","\\1", labs), levels=c("wt","mdx","d2"))
  if (length(unique(grp)) < 2) return(NULL)  # 群が1つのみは不可
  
  y <- DGEList(counts = as.matrix(pseudo_sc[, cols, drop=FALSE]))
  y <- calcNormFactors(y, method="TMMwsp")
  
  X <- model.matrix(~0 + grp)   # wt, mdx, d2 のダミー
  y <- estimateDisp(y, X, robust=TRUE)
  fit <- glmQLFit(y, X, robust=TRUE)
  
  contr <- makeContrasts(
    mdx_vs_wt = grpmdx - grpwt,
    d2_vs_mdx = grpd2  - grpmdx,
    d2_vs_wt  = grpd2  - grpwt, levels = X
  )
  
  list(
    mdx_vs_wt = topTags(glmQLFTest(fit, contrast = contr[,"mdx_vs_wt"]), n=Inf),
    d2_vs_mdx = topTags(glmQLFTest(fit, contrast = contr[,"d2_vs_mdx"]), n=Inf),
    d2_vs_wt  = topTags(glmQLFTest(fit, contrast = contr[,"d2_vs_wt"]),  n=Inf)
  )
}

ct_levels <- levels(muscle_sc$celltype)
edgeR_out <- setNames(lapply(ct_levels, edgeR_by_celltype), ct_levels)
edgeR_out <- Filter(Negate(is.null), edgeR_out)

## 9) ヒット数（全遺伝子 / 関心遺伝子）-----------------------------------
.safe_nrow <- function(x) if (is.null(x) || !is.data.frame(x)) 0L else nrow(x)
collect_hits <- function(lst, only_targets = FALSE) {
  lapply(lst, function(tts){
    lapply(tts, function(tt){
      if (is.null(tt)) return(NULL)
      tab <- tt$table
      if (only_targets) tab <- tab[rownames(tab) %in% gene_target_sc, , drop=FALSE]
      tab[tab$FDR < 0.05, , drop=FALSE]
    })
  })
}
hit_count <- function(sig){
  do.call(rbind, lapply(names(sig), function(ct){
    x <- sig[[ct]]
    data.frame(
      celltype    = ct,
      mdx_vs_wt_n = .safe_nrow(x$mdx_vs_wt),
      d2_vs_mdx_n = .safe_nrow(x$d2_vs_mdx),
      d2_vs_wt_n  = if ("d2_vs_wt" %in% names(x)) .safe_nrow(x$d2_vs_wt) else NA_integer_
    )}))
}

sig_full   <- collect_hits(edgeR_out, only_targets = FALSE)
sig_target <- collect_hits(edgeR_out, only_targets = TRUE)

cat("\n[Hit count: 全遺伝子]\n"); print(hit_count(sig_full),   row.names=FALSE)
cat("\n[Hit count: 関心遺伝子のみ]\n"); print(hit_count(sig_target), row.names=FALSE)


## ---- 3) 実行例 ------------------------------------------------------------
sig_target[["Endothelial cell(EC)"]][["d2_vs_mdx"]]
sig_target[["Endothelial cell(EC)"]][["mdx_vs_wt"]]
sig_target[["Endothelial cell(EC)"]][["d2_vs_wt"]]



## ---- 0) 前提 -------------------------------------------------------------
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

stopifnot(exists("pseudo_sc"))  # AggregateExpression の RNA-counts 合成結果（genes x samples）
esc <- function(s) gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", s))

## ---- 1) 代表遺伝子セット（必要に応じて差し替え） --------------------------
gene_sets <- list(
  ECM_org      = c("Fn1","Fbn2","Fbln1","Col1a1","Col3a1","Col4a1","Col6a1","Col8a1",
                   "Col14a1","Col16a1","Plod2","Postn","Vcan","Adam12"),
  Platelet_act = c("Cd63","Plek","Pf4","Selp","Fcer1g","Itgb1","Vwf","Col1a1","Fn1"),
  MHCII        = c("H2-Aa","H2-Ab1","H2-Eb1","Ciita"),
  
  # --- 脂質関連を分解 ---
  Lipo_Synthesis      = c("Acaca","Fasn","Scd1","Elovl6","Me1","Acly","Acss2",
                          "Gpam","Agpat2","Lpin1","Dgat1","Dgat2","Srebf1","Pparg","Cebpa","Cebpb"),
  FA_Uptake_Transport = c("Cd36","Lpl","Gpihbp1","Fabp3","Fabp4","Fabp5","Ldlr","Scarb1","Mfsd2a"),
  FAO_mito            = c("Cpt1a","Cpt1b","Acadvl","Acadm","Hadha","Hadhb","Echs1","Etfdh","Etfb","Acaa2"),
  FAO_perox           = c("Acox1","Ehhadh","Hsd17b4","Abcd3","Acot1","Pecr","Pex11a"),
  
  # 既存の代謝系（先生のまま）
  OxPhos          = c("Ndufa1","Ndufa2","Ndufb5","Cox4i1","Atp5f1a","Atp5f1b","Atp5mc1","Uqcrc1","Uqcrq",
                      "Ndufs2","Ndufb8","Cox5a","Cox6a1","Atp5me","Atp5md","Uqcrc2","Uqcrfs1"),
  TCA_cycle       = c("Cs","Aco2","Idh3a","Idh3b","Ogdh","Dlst","Sucb1","Sdha","Sdhb","Fh1","Mdh2","Pck2","Mpc1","Mpc2","Pdhb","Dlat"),
  Mito_Biogenesis = c("Ppargc1a","Ppargc1b","Nrf1","Nrf2","Tfam","Tfb1m","Tfb2m","Tomm20","Timm23","Polg")
)

## ---- 2) 共通ユーティリティ -------------------------------------------------
.get_cols_by_ct <- function(ct){
  pfx  <- paste0(ct, "_")
  cols <- colnames(pseudo_sc)[startsWith(colnames(pseudo_sc), pfx)]
  if (!length(cols)) stop(sprintf("列が見つかりません：celltype='%s'", ct))
  cols
}


.get_grp_from_labels <- function(labs){
  # labs 例: "wt1", "mdx_2", "d2-1" などを許容
  grp <- sub("^(wt|mdx|d2).*", "\\1", labs)
  factor(grp, levels = c("wt","mdx","d2"))
}

.build_design_nocov <- function(y, labs){
  grp <- .get_grp_from_labels(labs)
  if (all(table(grp) == 0)) stop("群が認識できません（列名の末尾に wt/mdx/d2 を含めてください）")
  model.matrix(~0 + grp)
}

.make_contrasts3 <- function(X){
  makeContrasts(
    mdx_vs_wt = grpmdx - grpwt,
    d2_vs_mdx = grpd2  - grpmdx,
    d2_vs_wt  = grpd2  - grpwt, levels = X
  )
}

.index_sets <- function(gene_sets, genes_in_matrix){
  idx_list <- lapply(gene_sets, function(gs){
    ix <- match(intersect(gs, genes_in_matrix), genes_in_matrix)
    ix <- ix[!is.na(ix)]
    if (length(ix) >= 3) ix else integer(0)
  })
  idx_list[sapply(idx_list, length) >= 3]
}

## ---- 3) camera/roast を一括実行 ----------------------------------
camera_all_sets <- function(ct,
                            contrast = c("d2_vs_mdx","mdx_vs_wt","d2_vs_wt"),
                            gene_sets){
  contrast <- match.arg(contrast)
  
  ## -- 細胞型の擬似バルクを抽出
  cols <- .get_cols_by_ct(ct)
  pb   <- pseudo_sc[, cols, drop = FALSE]
  
  ## -- ラベルと DGE オブジェクト
  labs <- sub(".*_", "", cols)   # 末尾トークン（例: wt1 / mdx2 / d2_1）
  y    <- DGEList(counts = as.matrix(pb))
  y    <- calcNormFactors(y, method = "TMMwsp")
  
  ## -- デザイン
  X <- .build_design_nocov(y, labs)
  # 注意：voom は design を必要とします。camera は voom オブジェクトを受け取ります。
  v <- voom(y, X, plot = FALSE)
  
  ## -- コントラスト
  contr <- .make_contrasts3(X)
  K <- switch(contrast,
              mdx_vs_wt = contr[,"mdx_vs_wt"],
              d2_vs_mdx = contr[,"d2_vs_mdx"],
              d2_vs_wt  = contr[,"d2_vs_wt"])
  
  ## -- 遺伝子セット index（>=3 遺伝子のみ）
  idx_list <- .index_sets(gene_sets, rownames(v$E))
  if (!length(idx_list)) stop("有効な遺伝子セットがありません（>=3 遺伝子が必要）。")
  
  ## -- camera / roast（回転乱数の再現性確保）
  set.seed(1)
  cam <- camera(v, index = idx_list, design = X, contrast = K)
  cam$FDR <- p.adjust(cam$PValue, "BH")
  
  rs <- lapply(idx_list, function(ix)
    roast(v, index = ix, design = X, contrast = K, nrot = 9999))
  
  list(celltype  = ct,
       contrast  = contrast,
       camera    = cam[order(cam$PValue), , drop = FALSE],
       roast     = rs,
       n_samples = table(.get_grp_from_labels(labs)),
       design_cols = colnames(X))
}

## ---- 4) 遺伝子セットの要約出力ヘルパ --------------------------------------
## 有意（FDR < alpha）のセットのみを表示する版に差し替え
print_top_sets <- function(res, n = NULL, alpha = 0.05, sort_by = c("FDR", "PValue")) {
  sort_by <- match.arg(sort_by)
  tab <- res$camera
  if (!"FDR" %in% names(tab)) {
    tab$FDR <- p.adjust(tab$PValue, method = "BH")
  }
  
  sig <- tab[tab$FDR < alpha, , drop = FALSE]
  if (nrow(sig) > 0L) {
    sig <- sig[order(sig[[sort_by]], sig$PValue), , drop = FALSE]
  }
  
  cat(sprintf("\n[%s | %s] camera 有意セット（FDR < %.3f）: %d件\n",
              res$celltype, res$contrast, alpha, nrow(sig)))
  
  if (nrow(sig) == 0L) {
    cat("  有意な遺伝子セットはありません。\n")
  } else {
    print(sig)
  }
  
  cat("\nサンプル数（wt/mdx/d2）:\n")
  print(res$n_samples)
  invisible(res)
}


## ---- 5) 使い方例（必要に応じて実行） ---------------------------------------
## d2_vs_mdx を検定
res_ec <- camera_all_sets("Endothelial cell(EC)", "d2_vs_mdx", gene_sets)
print_top_sets(res_ec, alpha = 0.05)  # 既定は 0.05
#                     NGenes Direction                         FDR
# Platelet_act             9        Up  0.000000000000000000001172269
# ECM_org                 14        Up  0.000000000000000618550412368
# FA_Uptake_Transport      9      Down  0.000299071754007749395290472
# MHCII                    4        Up  0.000763714781201935903436062
# FAO_mito                10        Up  0.041719318893955359417446971


res_st <- camera_all_sets("Stromal", "d2_vs_mdx", gene_sets)
print_top_sets(res_st, alpha = 0.05)
#                     NGenes Direction                         FDR
# ECM_org                 14        Up  0.000000000000000000002813363
# Platelet_act             9        Up  0.000000000000015125974029502
# MHCII                    4        Up  0.000026313914253506312773256
# FAO_mito                10        Up  0.003108951670574953626835502
# FA_Uptake_Transport      9        Up  0.020862924890370077718459996


res_satellite <- camera_all_sets("Satellite cell", "d2_vs_mdx", gene_sets)
print_top_sets(res_satellite, alpha = 0.05)
#                     NGenes Direction       FDR
# FA_Uptake_Transport      9      Down  0.0008503544
# Platelet_act             9        Up  0.0008503544
# Lipo_Synthesis          16      Down  0.0013746888
# Mito_Biogenesis          9      Down  0.0019878389
# TCA_cycle               15      Down  0.0042266178
# FAO_perox                7      Down  0.0284728427
# FAO_mito                10      Down  0.0441618366


res_Myo <- camera_all_sets("Myocyte/Myonucleus", "d2_vs_mdx", gene_sets)
print_top_sets(res_Myo, alpha = 0.05)
#                     NGenes Direction                             FDR
# TCA_cycle               15      Down  0.0000000000000000000000003262957
# OxPhos                  12      Down  0.0000000000000000000000007060666
# FAO_mito                10      Down  0.0000001228418553653926617152264
# ECM_org                 14        Up  0.0010810582945751461524108716361
# Mito_Biogenesis          9      Down  0.0014749817065457249132265360458
# Platelet_act             9        Up  0.0014749817065457249132265360458
# FAO_perox                7      Down  0.0020068826734276363349396277158
# Lipo_Synthesis          16      Down  0.0089570503044997962921458523056
# FA_Uptake_Transport      9      Down  0.0379222898914119699309566158263


res_CD53 <- camera_all_sets("Immune-Stromal(Cd53+)", "d2_vs_mdx", gene_sets)
print_top_sets(res_CD53, alpha = 0.05)
#                 NGenes Direction       FDR
# OxPhos              12        Up  1.160267e-08
# TCA_cycle           15        Up  2.489921e-06
# Platelet_act         9        Up  1.596754e-04
# FAO_mito            10        Up  1.562747e-03
# MHCII                4        Up  7.714488e-03
# ECM_org             14        Up  1.587971e-02
# FAO_perox            7        Up  1.763567e-02
# Mito_Biogenesis      9        Up  1.842259e-02


res_DC <- camera_all_sets("Dendritic cell(DC)", "d2_vs_mdx", gene_sets)
print_top_sets(res_DC, alpha = 0.05)
#                     NGenes Direction               FDR
# ECM_org                 14        Up  0.00000000000100843
# MHCII                    4      Down  0.00012466696860586
# TCA_cycle               15        Up  0.00038725561827587
# OxPhos                  12        Up  0.00645226525953102
# FA_Uptake_Transport      9      Down  0.01437414027789251
# Mito_Biogenesis          9        Up  0.01437414027789251
# FAO_mito                10        Up  0.01437414027789251
# Platelet_act             9        Up  0.01437697544269053


res_High <- camera_all_sets("MHC-high macrophage", "d2_vs_mdx", gene_sets)
print_top_sets(res_High, alpha = 0.05)
#                     NGenes Direction                          FDR
# ECM_org                 14        Up  0.0000000000000000000005564824
# FA_Uptake_Transport      9      Down  0.0000000019965528887113431281
# Platelet_act             9        Up  0.0004036007544234108393096938
# FAO_mito                10        Up  0.0004036007544234108393096938
# Lipo_Synthesis          16      Down  0.0194314816257071670824174703
# MHCII                    4      Down  0.0366327865794700349710844023


res_Low <- camera_all_sets("MHC-low macrophage", "d2_vs_mdx", gene_sets)
print_top_sets(res_Low, alpha = 0.05)
#              NGenes Direction                        FDR
# MHCII             4      Down  0.00000000000000000001820953
# ECM_org          14        Up  0.00000000000001320152470235
# OxPhos           12        Up  0.00000582855050288573875049
# TCA_cycle        15        Up  0.00008249127537354697012675
# FAO_mito         10        Up  0.00008249127537354697012675
# Platelet_act      9        Up  0.00122406352979519106836948


res_Neu <- camera_all_sets("Neutrophil", "d2_vs_mdx", gene_sets)
print_top_sets(res_Neu, alpha = 0.05)
#                     NGenes Direction                          FDR
# Platelet_act             9        Up  0.0000000000000000000003620995
# ECM_org                 14        Up  0.0000000000000000000003620995
# MHCII                    4        Up  0.0000000000000000101501969115
# FA_Uptake_Transport      9        Up  0.0022632851880849016827645936


res_Peri <- camera_all_sets("Pericyte/Vascular SMC", "d2_vs_mdx", gene_sets)
print_top_sets(res_Peri, alpha = 0.05)
#                     NGenes Direction                     FDR
# Platelet_act             9        Up  0.00000000000000001624227
# ECM_org                 14        Up  0.00000000000000023127067
# FA_Uptake_Transport      9        Up  0.00005855833823792290364
# FAO_mito                10        Up  0.00013794826963615022786
# OxPhos                  12        Up  0.00742452651312596172833
# Mito_Biogenesis          9      Down  0.02877214841419776722442


res_Sch <- camera_all_sets("Schwann cell", "d2_vs_mdx", gene_sets)
print_top_sets(res_Sch, alpha = 0.05)
#                NGenes Direction         FDR
# Platelet_act        9        Up  0.00001171314
# ECM_org            14      Down  0.00002527173
# Lipo_Synthesis     16      Down  0.00022611242


res_Teno <- camera_all_sets("Tenocyte/Tendon fibroblast", "d2_vs_mdx", gene_sets)
print_top_sets(res_Teno, alpha = 0.05)
#                     NGenes Direction                            FDR
# ECM_org                 14        Up  0.000000000000000000000364033
# Platelet_act             9        Up  0.000000000000005396598576695
# OxPhos                  12        Up  0.000854426868896658326529847
# FAO_mito                10        Up  0.000854426868896658326529847
# MHCII                    4        Up  0.001508411467867410262269634
# TCA_cycle               15        Up  0.001508411467867410262269634
# FA_Uptake_Transport      9        Up  0.004698316744163056821181890


## ---------- ①-0 準備 ----------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr)
})
dir.create("C:/DMD_project/figs", showWarnings = FALSE, recursive = TRUE)

stopifnot(exists("edgeR_out"), exists("pseudo_sc"))
mhc_core <- c("H2-Aa","H2-Ab1","H2-Eb1","Ciita")

## ---------- ①-A 代表遺伝子パネル：ドット・ヒートマップ ----------
# edgeR_out -> tidy
ed_list <- lapply(names(edgeR_out), function(ct){
  lapply(names(edgeR_out[[ct]]), function(contr){
    tt <- edgeR_out[[ct]][[contr]]$table
    if (is.null(tt) || nrow(tt)==0) return(NULL)
    transform(as.data.frame(tt), gene = rownames(tt),
              celltype = ct, contrast = contr, row.names = NULL)
  }) |> dplyr::bind_rows()
}) |> dplyr::bind_rows()

panel_genes <- unique(c(
  mhc_core,                       # MHC-II core
  c("Fn1","Lum","Thbs2","Col1a1","Col15a1","Postn","Vcan"),  # ECM/線維化
  c("Lpl","Spp1","C1qa","C1qb","C1qc","Cd53","Ctss"),        # FIM/免疫
  c("Tnni2","Des","Pax7","Myf5"),                            # 筋/サテライト
  c("Pecam1","Cdh5","Pdgfra","Pi16","Rgs5","Abcc9","Scx","Tnmd") # EC/間質/周皮/Tendon
))

df_panel <- ed_list |>
  dplyr::filter(contrast == "d2_vs_mdx", gene %in% panel_genes)

p_panel <- ggplot(df_panel, aes(x = gene, y = celltype)) +
  geom_point(aes(size = pmin(-log10(FDR), 10), color = logFC)) +
  scale_color_gradient2(low = "steelblue", mid = "grey95", high = "firebrick", midpoint = 0) +
  scale_size(range = c(0.5, 5)) +
  labs(title = "Pseudo Bulk（d2_vs_mdx）：The effect size and statistical significance of representative genes",
       x = NULL, y = NULL, color = "logFC", size = "-log10(FDR)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
p_panel

## ---------- ①-B camera 結果：遺伝子セット×細胞型タイル ----------
stopifnot(exists("camera_all_sets"), exists("gene_sets"))

cts_use <- intersect(names(edgeR_out), c(
  "MHC-low macrophage","MHC-high macrophage","Dendritic cell(DC)","Neutrophil",
  "Stromal","Immune-Stromal(Cd53+)","Endothelial cell(EC)","Pericyte/Vascular SMC","Tenocyte/Tendon fibroblast",
  "Schwann cell","Satellite cell","Myocyte/Myonucleus"
))

cam_tab <- lapply(cts_use, function(ct){
  res <- camera_all_sets(ct, "d2_vs_mdx", gene_sets)
  if (is.null(res) || is.null(res$camera)) return(NULL)
  x <- res$camera
  x$celltype <- ct
  x$geneset  <- rownames(x)
  x
}) |> dplyr::bind_rows()

df_cam <- cam_tab |>
  mutate(dir_num = ifelse(Direction == "Up", 1, -1),
         signed  = dir_num * pmax(-log10(FDR), 0),
         geneset = factor(geneset, levels = unique(geneset))) |>
  filter(FDR < 0.05)  # 有意のみを色付け

p_cam <- ggplot(df_cam, aes(x = geneset, y = celltype, fill = signed)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0,
                       name = "Direction ×\n -log10(FDR)") +
  labs(title = "Pseudo Bulk（d2_vs_mdx）：Statistical significance of gene sets（camera）", x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p_cam

## ---------- ①-C MHC-II コアの火山図（細胞型ファセット） ----------
lab_genes <- mhc_core
df_mhc <- ed_list |>
  dplyr::filter(contrast == "d2_vs_mdx")

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
base_vol <- ggplot(df_mhc, aes(x = logFC, y = -log10(FDR))) +
  geom_point(alpha = 0.30, size = 0.8) +
  geom_point(data = subset(df_mhc, gene %in% lab_genes),
             color = "firebrick", size = 1.5) +
  labs(title = "Pseudo Bulk（d2_vs_mdx）：Volcano Plot: MHC-II Core",
       x = "log2 fold-change", y = "-log10(FDR)") +
  facet_wrap(~ celltype, scales = "free_y") +
  theme_bw(base_size = 9)

p_vol <- if (has_ggrepel) {
  base_vol + ggrepel::geom_text_repel(
    data = subset(df_mhc, gene %in% lab_genes),
    aes(label = gene), size = 3, max.overlaps = Inf)
} else base_vol; p_vol


## 4-3 一通りの解析が完了したら-------------------------------------------------
#saveRDS(object   = muscle_sc,file = "C:/DMD_project/rds/muscle_sc_final.rds",compress = "gzip")











































## ======================================================================
## NicheNet × Seurat v5  再現可能フルスクリプト（最新安定 / 2025-08-14）
## 目的：
##  - 環境依存を避けて “必ず走る”
##  - 途中で止まっても “なぜ止まったか” がログに残る
##  - 最新 CRAN/Bioc + nichenetr v2.2.0 を子プロセスで導入
##  - Seurat v5 Assay5/layer 差分に対応（冪等）
##  - aggregate ラッパーで失敗時は公開APIへ自動フォールバック
## ======================================================================
## ---- 0. 設定：パス（必要に応じて調整） --------------------------------------
logfile <- "C:/DMD_project/exprdata/nichenet_run_mhc_high.log"
rds_dir <- "C:/DMD_project/exprdata"
paths <- list(
  ligand_target_matrix = file.path(rds_dir, "ligand_target_matrix_nsga2r_final_mouse.rds"),
  lr_network           = file.path(rds_dir, "lr_network_mouse_21122021.rds"),
  weighted_networks    = file.path(rds_dir, "weighted_networks_nsga2r_final_mouse.rds"),
  seurat_obj           = "C:/DMD_project/rds/muscle_sc_final.rds"
)

## ---- 1. ログ：安全な開始／終了ラッパ ----------------------------------------
start_log <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt", encoding = "UTF-8")
  sink(con, split = TRUE)              # 標準出力
  sink(con, type = "message")          # メッセージ
  assign(".__log_con", con, envir = .GlobalEnv)
  invisible(TRUE)
}
stop_log <- function() {
  ## message用 sink を1段だけ閉じる（開いていれば）
  try(sink(type = "message"), silent = TRUE)
  ## 出力用 sink を全部閉じる
  while (sink.number() > 0) { try(sink(), silent = TRUE) }
  ## コネクションを閉じる
  if (exists(".__log_con", envir = .GlobalEnv)) {
    con <- get(".__log_con", envir = .GlobalEnv)
    try(close(con), silent = TRUE)
    rm(".__log_con", envir = .GlobalEnv)
  }
  invisible(TRUE)
}

## 0) まっさらなセッション推奨。開きっぱなしの sink があれば閉じる
while (sink.number() > 0) sink()
## ---- 2. セクション実行ユーティリティ（ログ一体） ------------------------------
say <- function(...) { cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n") }
section <- function(title) { cat("\n", paste0("---- ", title, " ----"), "\n", sep="") }
run_step <- function(title, expr) {
  section(title)
  tryCatch(force(expr),
           error = function(e) {
             say("ERROR in ", title, ":", conditionMessage(e))
             cat("Calls:\n"); print(sys.calls())
             stop("Step failed: ", title, call. = FALSE)
           },
           warning = function(w) {
             say("WARN in ", title, ":", conditionMessage(w))
             invokeRestart("muffleWarning")
           }
  )
}

## ---- 3. ログ開始 --------------------------------------------------------------
start_log(logfile)

## —— エラーを必ず表示（終了処理なし・最短）——
options(show.error.messages = TRUE, warn = 1)

options(error = function(){
  ## ログの sink に奪われていても必ずコンソールに出すため、まず“消音”を解除
  try({ while (sink.number(type = "message") > 0) sink(type = "message") }, silent = TRUE)
  try({ while (sink.number() > 0) sink() }, silent = TRUE)
  
  cat("\n[ERROR] ", geterrmessage(), "\n", sep = "")
  utils::traceback()   # スタックを表示
  flush.console()
})

say("[START]")



## 3) ライブラリ読み込み（遅延ロード対策で suppress）
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(nichenetr)
  library(cli)
})



cli::cli_h1("Environment check")
cli::cli_inform(c(
  "i" = paste("Seurat", as.character(packageVersion("Seurat"))),
  "i" = paste("SeuratObject", as.character(packageVersion("SeuratObject"))),
  "i" = paste("nichenetr", as.character(packageVersion("nichenetr")))
))

## 4) 入力ファイル（先生のパスに合わせています）
rds_dir <- "C:/DMD_project/exprdata"
paths <- list(
  ligand_target_matrix = file.path(rds_dir, "ligand_target_matrix_nsga2r_final_mouse.rds"),
  lr_network           = file.path(rds_dir, "lr_network_mouse_21122021.rds"),
  weighted_networks    = file.path(rds_dir, "weighted_networks_nsga2r_final_mouse.rds"),
  seurat_obj           = "C:/DMD_project/rds/muscle_sc_final.rds"
)
read_if <- function(p) {
  if (!file.exists(p)) stop("ファイルが見つかりません: ", p, call. = FALSE)
  readRDS(p)
}

cli::cli_h1("Load data")
muscle_sc            <- read_if(paths$seurat_obj)
ligand_target_matrix <- read_if(paths$ligand_target_matrix)
lr_network           <- read_if(paths$lr_network)
weighted_networks    <- read_if(paths$weighted_networks)

cli::cli_inform(c(
  "v" = paste("Assay RNA layers:", paste(Layers(muscle_sc[["RNA"]]), collapse = ", ")),
  "v" = paste("Cells:", ncol(muscle_sc), "| Genes:", nrow(muscle_sc))
))

## 5) 種（human/mouse）の整合性チェックと必要な場合の置換
is_mouse_object <- "Ptprc" %in% rownames(muscle_sc)
has_hla_in_ltm  <- any(grepl("^HLA-", rownames(ligand_target_matrix)))
if (is_mouse_object && has_hla_in_ltm) {
  cli::cli_inform(c(">" = "Human→Mouse 記号変換を実施（HLA-* → H2-* など）"))
  human2mouse <- nichenetr::convert_human_to_mouse_symbols
  rownames(ligand_target_matrix) <- human2mouse(rownames(ligand_target_matrix))
  colnames(ligand_target_matrix) <- human2mouse(colnames(ligand_target_matrix))
  lr_network$from <- human2mouse(lr_network$from)
  lr_network$to   <- human2mouse(lr_network$to)
  weighted_networks$from <- human2mouse(weighted_networks$from)
  weighted_networks$to   <- human2mouse(weighted_networks$to)
}

## 6) Seurat v5 の layer を安全に整える（Assay5 で "data" が無い時に作成）
DefaultAssay(muscle_sc) <- "RNA"
## 複数の counts.* があれば単一 "counts" に統合（JoinLayers）
if (sum(grepl("^counts", Layers(muscle_sc[["RNA"]]))) > 1) {
  muscle_sc <- JoinLayers(muscle_sc)  # Seurat v5 推奨の layer 統合
}
## "data" layer が無ければ NormalizeData で生成
if (!"data" %in% Layers(muscle_sc[["RNA"]])) {
  muscle_sc <- NormalizeData(muscle_sc, assay = "RNA",
                             normalization.method = "LogNormalize", verbose = FALSE)
}
## デフォルト表示/取得に使う layer を "data" に
DefaultLayer(muscle_sc[["RNA"]]) <- "data"

## 7) 解析パラメータ
Idents(muscle_sc) <- "celltype"
receiver_oi <- "Myocyte/Myonucleus"
sender_oi   <- c("MHC-high macrophage","MHC-low macrophage","Stromal")
# candidate
# "MHC-low macrophage", # 0 low
# "MHC-high macrophage", # 2 high
# "Endothelial cell(EC)", # 3
# "Stromal", # 5
# "Immune-Stromal(Cd53+)", # 6 low
# "Myocyte/Myonucleus", # 11
# "Tenocyte/Tendon fibroblast", # 13
# "Neutrophil",  # 16
# "Dendritic cell(DC)", # 17 high
# "Satellite cell", # 19
# "Pericyte/Vascular SMC", # 20
# "Schwann cell", # 21
condition_col       <- "genotype"
condition_reference <- "mdx"
condition_oi        <- "mdxD2"

## 送受信セル数の確認（0 なら即停止）
recv_cells <- WhichCells(muscle_sc, idents = receiver_oi)
send_cells <- WhichCells(muscle_sc, idents = sender_oi)
stopifnot(length(recv_cells) > 0, length(send_cells) > 0)
cli::cli_inform(c(
  "✔" = paste("Receiver:", length(recv_cells), "cells"),
  "✔" = paste("Sender  :", length(send_cells), "cells")
))

# Receiver: 1805 cells
# Sender : 9151 cells


## 8) NicheNet 本体：aggregate ラッパーを素直に呼ぶ
##    （assay_oi = 'RNA' を明示。expression_pct は先生の既定 0.1）
ligand_target_matrix <- as.matrix(ligand_target_matrix)

cli::cli_h1("Run NicheNet (aggregate wrapper)")
nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj           = muscle_sc,
  receiver             = receiver_oi,
  condition_colname    = condition_col,
  condition_oi         = condition_oi,
  condition_reference  = condition_reference,
  sender               = sender_oi,
  ligand_target_matrix = ligand_target_matrix,
  lr_network           = lr_network,
  weighted_networks    = weighted_networks,
  assay_oi             = "RNA",
  expression_pct       = 0.10,
  verbose              = TRUE
)

## 9) 最低限の検証と可視化（ゼロ列リガンド等は除外して描画）
stopifnot(!is.null(nichenet_output$ligand_activities),
          nrow(nichenet_output$ligand_activities) > 0)

## 全ゼロのリガンド列があれば除外（まれにリソースと記号が噛み合わない時の保険）
zero_ligands <- colnames(ligand_target_matrix)[colSums(ligand_target_matrix) == 0]
if (length(zero_ligands)) {
  nichenet_output$ligand_activities <- subset(
    nichenet_output$ligand_activities,
    !(test_ligand %in% zero_ligands)
  )
  if (!is.null(nichenet_output$ligand_target_links)) {
    nichenet_output$ligand_target_links <- subset(
      nichenet_output$ligand_target_links,
      !(ligand %in% zero_ligands)
    )
  }
}

## 推奨の複合ヒートマップ（v2 系の標準プロット）
if (!is.null(nichenet_output$ligand_activity_target_heatmap)) {
  print(nichenet_output$ligand_activity_target_heatmap)
}



suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# 横並び 2 パネル（Expression・LFC）
nnk_expr_lfc_two_panel <- function(
    seurat_obj,
    senders,                       # 例: c("Pericyte/Vascular SMC","Satellite cell")
    cond_col, cond_oi, cond_ref,   # 条件列と比較（oi - ref）
    ligands_to_plot = NULL,        # 例: c("Col18a1","Cd47","Calm1",...)
    nichenet_output = NULL,        # あれば top ligands を自動抽出に使用
    assay = "RNA",
    expr_slot = c("data","counts","scale.data")[1],
    pct_threshold = 0,             # 0 だと全ての点を表示（非発現は size=0）
    lfc_eps = 1e-9,                # LFC のゼロ割回避
    dot_max_size = 6,
    lfc_limits = c(-3, 3),
    title_left  = "Expression in Sender",
    title_right = "LFC in Sender"
){
  stopifnot(inherits(seurat_obj, "Seurat"))
  DefaultAssay(seurat_obj) <- assay
  Idents(seurat_obj) <- "celltype"
  
  # 1) リガンド集合（引数 > nichenet_output > 例示トップ10 の順で決定）
  if (is.null(ligands_to_plot)) {
    if (!is.null(nichenet_output) && !is.null(nichenet_output$ligand_activities)) {
      nnk_la <- nichenet_output$ligand_activities
      nnk_score <- if ("aupr" %in% names(nnk_la)) "aupr" else if ("pearson" %in% names(nnk_la)) "pearson" else names(nnk_la)[2]
      nnk_lig_col <- intersect(c("test_ligand","ligand"), names(nnk_la))[1]
      ligands_to_plot <- head(nnk_la[order(-nnk_la[[nnk_score]]), nnk_lig_col, drop=TRUE], 10)
    } else {
      stop("ligands_to_plot を指定するか、nichenet_output に ligand_activities を含めてください。")
    }
  }
  ligands_to_plot <- unique(ligands_to_plot)
  
  # 2) 行列取得（data スロットが無ければ counts を使用）
  nnk_mat <- tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = expr_slot),
    error = function(e) GetAssayData(seurat_obj, assay = assay, slot = expr_slot)
  )
  nnk_mat <- nnk_mat[intersect(rownames(nnk_mat), ligands_to_plot), , drop = FALSE]
  
  # メタデータ
  nnk_meta <- seurat_obj[[]]
  stopifnot(all(c(cond_col) %in% colnames(nnk_meta)))
  
  # 3) Expression in Sender（各 sender × ligand：発現率％・平均発現）
  nnk_expr_list <- lapply(senders, function(sdr){
    nnk_cells <- WhichCells(seurat_obj, idents = sdr)
    if (!length(nnk_cells)) return(NULL)
    sub <- nnk_mat[, nnk_cells, drop = FALSE]
    data.frame(
      sender  = sdr,
      ligand  = rownames(sub),
      pct     = Matrix::rowSums(sub > 0) / ncol(sub) * 100,
      avg_exp = as.numeric(Matrix::rowMeans(sub)),
      stringsAsFactors = FALSE
    )
  })
  nnk_expr_df <- bind_rows(nnk_expr_list)
  if (!nrow(nnk_expr_df)) stop("指定 sender に該当する細胞が見つかりません。")
  
  # 4) LFC in Sender（各 sender × ligand：log2(oi/ref)）
  nnk_lfc_list <- lapply(senders, function(sdr){
    all_cells <- WhichCells(seurat_obj, idents = sdr)
    if (!length(all_cells)) return(NULL)
    oi_cells  <- intersect(all_cells, rownames(nnk_meta)[nnk_meta[[cond_col]] == cond_oi])
    ref_cells <- intersect(all_cells, rownames(nnk_meta)[nnk_meta[[cond_col]] == cond_ref])
    if (length(oi_cells) == 0 || length(ref_cells) == 0) return(NULL)
    
    sub_oi  <- nnk_mat[, oi_cells,  drop = FALSE]
    sub_ref <- nnk_mat[, ref_cells, drop = FALSE]
    
    m_oi  <- Matrix::rowMeans(sub_oi)
    m_ref <- Matrix::rowMeans(sub_ref)
    lfc   <- log2((m_oi + lfc_eps) / (m_ref + lfc_eps))
    
    data.frame(
      sender = sdr,
      ligand = names(lfc),
      lfc    = as.numeric(lfc),
      stringsAsFactors = FALSE
    )
  })
  nnk_lfc_df <- bind_rows(nnk_lfc_list)
  
  # 5) 因子順序（横＝リガンド、縦＝sender）
  nnk_expr_df$ligand <- factor(nnk_expr_df$ligand, levels = rev(ligands_to_plot))
  nnk_lfc_df$ligand  <- factor(nnk_lfc_df$ligand,  levels = rev(ligands_to_plot))
  nnk_expr_df$sender <- factor(nnk_expr_df$sender, levels = senders)
  nnk_lfc_df$sender  <- factor(nnk_lfc_df$sender,  levels = senders)
  
  # 6) プロット作成
  nnk_expr_plot <- ggplot(nnk_expr_df, aes(x = ligand, y = sender)) +
    geom_point(aes(size = pmax(pct, pct_threshold), fill = avg_exp),
               shape = 21, colour = "grey30", stroke = 0.25) +
    scale_size_area(max_size = dot_max_size, breaks = c(0,25,50,75,100),
                    name = "Percent expressed (%)") +
    scale_fill_distiller(palette = "RdYlBu", direction = -1,
                         name = "Average expression") +
    labs(title = title_left, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))+
    coord_flip()
  
  nnk_lfc_plot <- ggplot(nnk_lfc_df, aes(x = ligand, y = sender, fill = lfc)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#2166AC", high = "#B2182B", mid = "white",
                         limits = lfc_limits, oob = scales::squish,
                         name = "LFC (oi - ref)") +
    labs(title = title_right, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))+
    coord_flip()
  
  # 7) 横並びで返す
  nnk_combined <- nnk_expr_plot + nnk_lfc_plot +
    plot_layout(ncol = 2, widths = c(1, 1))
  
  return(list(
    expr_df   = nnk_expr_df,
    lfc_df    = nnk_lfc_df,
    expr_plot = nnk_expr_plot,
    lfc_plot  = nnk_lfc_plot,
    combined  = nnk_combined
  ))
}



# 例: Pericyte/Vascular SMC と Satellite を比較し、v2 の標準パネルと同じ上位 9 リガンドを表示
## ===== 全ヒット（top9廃止）版：呼び出しブロックをこれに置換 =====

## NicheNetのligand_activitiesから“全ヒット”を作る
la <- nichenet_output$ligand_activities
lig_col <- intersect(c("test_ligand","ligand","ligand_gene","ligand_symbol"), colnames(la))[1]
if (is.na(lig_col) || is.null(lig_col)) stop("nichenet_output$ligand_activities にリガンド列が見つかりません。")

if ("hit" %in% colnames(la)) {
  ligands_all <- unique(na.omit(as.character(la[[lig_col]][la$hit %in% c(TRUE, 1)])))
} else if ("pearson" %in% colnames(la)) {
  ligands_all <- unique(na.omit(as.character(la[[lig_col]][la$pearson > 0])))
} else {
  ligands_all <- unique(na.omit(as.character(la[[lig_col]])))
}

## データ内に実在する遺伝子だけに限定（Seurat v4/v5両対応）
nnk_mat_tmp <- tryCatch(
  GetAssayData(muscle_sc, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(muscle_sc, assay = "RNA", slot = "data")
)
ligands_all <- intersect(ligands_all, rownames(nnk_mat_tmp))
if (length(ligands_all) == 0) stop("全ヒット候補に、データ内のリガンドがありません。")

## 上位9のベクトル指定をやめ、全ヒットを渡す
nnk_panels <- nnk_expr_lfc_two_panel(
  seurat_obj      = muscle_sc,
  senders         = c("MHC-high macrophage","MHC-low macrophage","Stromal"),
  cond_col        = "genotype",
  cond_oi         = "mdxD2",
  cond_ref        = "mdx",
  ligands_to_plot = ligands_all,   # ← ここが“全ヒット”
  nichenet_output = NULL           # ← すでに抽出済みなので不要
)
# candidate
# "MHC-low macrophage", # 0 low
# "MHC-high macrophage", # 2 high
# "Endothelial cell(EC)", # 3
# "Stromal", # 5
# "Immune-Stromal(Cd53+)", # 6 low
# "Myocyte/Myonucleus", # 11
# "Tenocyte/Tendon fibroblast", # 13
# "Neutrophil",  # 16
# "Dendritic cell(DC)", # 17 high
# "Satellite cell", # 19
# "Pericyte/Vascular SMC", # 20
# "Schwann cell", # 21
## 表示
print(nnk_panels$expr_plot)
print(nnk_panels$lfc_plot)
print(nnk_panels$combined)


# それぞれを表示
print(nnk_panels$expr_plot)
print(nnk_panels$lfc_plot)

# 横並び完成図
print(nnk_panels$combined)

# 保存する場合
# ggsave("expr_lfc_two_panel.pdf", nnk_panels$combined, width = 9, height = 4.5, useDingbats = FALSE)





## ---- 15. 簡易サマリ出力 -------------------------------------------------------
say("diag class(nichenet_output)=", paste(class(nichenet_output), collapse=", "),
    " | typeof=", typeof(nichenet_output),
    " | names(head)=",
    paste(utils::head(names(nichenet_output), 6), collapse=", "))

section <- function(title) cat("\n", paste0("---- ", title, " ----"), "\n", sep = "")

run_step <- function(title, expr) {
  section(title)
  exp <- substitute(expr)                     # 式を捕まえる
  tryCatch(
    withCallingHandlers(
      eval(exp, envir = parent.frame()),      # ★ 親フレームで評価
      warning = function(w) {                 # 警告はログに出すが止めない
        say("WARN in ", title, ": ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {                     # エラーは詳細を出して止める
      say("ERROR in ", title, ": ", conditionMessage(e))
      cat("Calls:\n"); print(sys.calls())
      stop("Step failed: ", title, call. = FALSE)
    }
  )
}

## ======================================================================
## ---- 16. Sanity report (numbers you can trust) -------------------------------
expr_pct_safe <- get0(".expr_pct", ifnotfound = 0.10)

run_step("Sanity report", {
  la <- nichenet_output$ligand_activities
  if (is.null(la) || nrow(la) == 0L) stop("ligand_activities が空です。")
  
  total_ligands   <- nrow(la)
  pos_pearson     <- sum(la$pearson > 0, na.rm = TRUE)
  top10           <- head(la[order(-la$pearson), c("test_ligand","pearson")], 10)
  say("Total candidate ligands: ", total_ligands, " | pearson>0: ", pos_pearson)
  say("Top10 ligands: ", paste(top10$test_ligand, collapse = ", "))
  
  ## LR ネットワークにおける裏付け（受け手/送信側の発現と合致するエッジ数）
  ## 受け手・送信側の“発現遺伝子”は Preflight でのロジックを再利用
  Seurat::Idents(muscle_sc) <- "celltype"
  recv_cells <- Seurat::WhichCells(muscle_sc, idents = receiver_oi)
  
  get_expressed_genes_v5 <- function(object, cells, assay = "RNA", pct = .expr_pct) {
    m <- SeuratObject::GetAssayData(object = object, assay = assay)
    sub <- m[, cells, drop = FALSE]
    frac <- Matrix::rowSums(sub > 0) / ncol(sub)
    rownames(m)[which(frac >= pct)]
  }
  expr_recv   <- get_expressed_genes_v5(muscle_sc, recv_cells, assay = "RNA", pct = .expr_pct)
  expr_sender <- unique(unlist(lapply(sender_oi, function(s)
    get_expressed_genes_v5(muscle_sc, Seurat::WhichCells(muscle_sc, idents = s),
                           assay = "RNA", pct = .expr_pct)
  )))
  
  lr_ok <- subset(lr_network, from %in% la$test_ligand & to %in% expr_recv)
  lr_ok <- subset(lr_ok, from %in% expr_sender)
  say("LR edges supported by expression: ", nrow(lr_ok))
  
  ## LTM の“空列”に該当するリガンドの混入チェック（ゼロの再検証）
  zero_ligands <- colnames(ligand_target_matrix)[colSums(ligand_target_matrix) == 0]
  n_zero <- sum(la$test_ligand %in% zero_ligands)
  if (n_zero > 0L) say("Ligands with zero LTM column detected: ", n_zero, " (they will be ignored downstream)")
})


#---- Sanity report ----
# 2025-08-14 10:15:40 | Total candidate ligands:  215  | pearson>0:  190 
# 2025-08-14 10:15:40 | Top10 ligands:  Prg4, Rarres2, Clec11a, Saraf, Tgfb2, Timp1, Hspg2, Jam2, Efnb1, Icam1 
# 2025-08-14 10:15:41 | LR edges supported by expression:  13

# ---- Sanity report ---- EC and Peri -> Myocyte
# 2025-08-15 11:03:59 | Total candidate ligands:  14  | pearson>0:  14 
# 2025-08-15 11:03:59 | Top10 ligands:  Calm1, Col18a1, Timp1, Calm2, Adm, Lamb2, Calm3, Ptgfrn, Angptl4, Cd47 
# 2025-08-15 11:04:00 | LR edges supported by expression:  14

# ---- Sanity report ---- Mac and Stromal -> Myocyte
# 2025-08-15 11:37:40 | Total candidate ligands:  15  | pearson>0:  15 
# 2025-08-15 11:37:40 | Top10 ligands:  Calm1, Col18a1, Timp1, Calm2, Adm, Lamb2, Calm3, Ptgfrn, Angptl4, Cd47 
# 2025-08-15 11:37:41 | LR edges supported by expression:  15

# ---- Sanity report ----Immune-stromal and EC -> myocyte
# 2025-08-15 11:59:17 | Total candidate ligands:  13  | pearson>0:  13 
# 2025-08-15 11:59:17 | Top10 ligands:  Calm1, Col18a1, Timp1, Calm2, Adm, Lamb2, Calm3, Ptgfrn, Cd47, App 
# 2025-08-15 11:59:18 | LR edges supported by expression:  13

# ---- Sanity report ----Peri and Satellite
# 2025-08-20 13:40:59 | Total candidate ligands:  15  | pearson>0:  15 
# 2025-08-20 13:40:59 | Top10 ligands:  Calm1, Col18a1, Timp1, Calm2, Adm, Lamb2, Calm3, Ptgfrn, Angptl4, Cd47 
# 2025-08-20 13:41:00 | LR edges supported by expression:  15

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(Matrix); library(circlize)
})

# 送受信クラスタ名（任意に変更可）
receiver_oi <- "Myocyte/Myonucleus"
sender_sat  <- c("Pericyte/Vascular SMC")
sender_str  <- "Satellite cell"
# candidate
# "MHC-low macrophage", # 0 low
# "MHC-high macrophage", # 2 high
# "Endothelial cell(EC)", # 3
# "Stromal", # 5
# "Immune-Stromal(Cd53+)", # 6 low
# "Myocyte/Myonucleus", # 11
# "Tenocyte/Tendon fibroblast", # 13
# "Neutrophil",  # 16
# "Dendritic cell(DC)", # 17 high
# "Satellite cell", # 19
# "Pericyte/Vascular SMC", # 20
# "Schwann cell", # 21

# 1) クラスタ名と細胞数を確認（0 ならここで止める）
lv <- levels(Idents(muscle_sc))
cat("Idents levels:\n"); print(lv)
nc <- table(Idents(muscle_sc))
cat("Cell counts:\n"); print(nc[c(sender_sat, sender_str, receiver_oi)])
stopifnot(length(WhichCells(muscle_sc, idents = receiver_oi)) > 0)
stopifnot(length(WhichCells(muscle_sc, idents = sender_sat )) > 0)
stopifnot(length(WhichCells(muscle_sc, idents = sender_str )) > 0)

# 2) 描画デバイスを必ず開く（Windows例）
if (identical(names(dev.cur()), "null device")) windows(width = 13, height = 8)


# 0) 発現遺伝子の抽出（Seurat v5 を想定：layer="data" を優先）
pct_expr <- 0.10
get_genes <- function(idents, pct = 0.10){
  cells <- WhichCells(muscle_sc, idents = idents)
  m <- tryCatch(
    SeuratObject::GetAssayData(muscle_sc, assay = "RNA", layer = "data"),
    error = function(e) SeuratObject::GetAssayData(muscle_sc, assay = "RNA")
  )
  sub <- m[, cells, drop = FALSE]
  rownames(m)[ Matrix::rowSums(sub > 0) / ncol(sub) >= pct ]
}
genes_sat  <- get_genes(sender_sat,  pct_expr)
genes_str  <- get_genes(sender_str,  pct_expr)
genes_recv <- get_genes(receiver_oi, pct_expr)

# 1) LTM の向きを統一（行=ligand / 列=target）
mat <- as.matrix(ligand_target_matrix)
lig_act <- nichenet_output$ligand_activities
stopifnot(is.data.frame(lig_act), nrow(lig_act) > 0)
lig_col   <- intersect(c("test_ligand","ligand"), names(lig_act))[1]
score_col <- if ("aupr" %in% names(lig_act)) "aupr" else "pearson"
lig_ranked <- lig_act[order(-lig_act[[score_col]]), , drop = FALSE]
ligs  <- lig_ranked[[lig_col]]
if (sum(ligs %in% rownames(mat)) < sum(ligs %in% colnames(mat))) mat <- t(mat)

# 2) 候補集合
lig_all <- rownames(mat); tgt_all <- colnames(mat)
tgt_key <- intersect(c("top_targets","prioritized_targets"), names(nichenet_output))[1]
targets_all <- if (length(tgt_key)) unique(unlist(nichenet_output[[tgt_key]])) else tgt_all
lig_ok <- intersect(ligs, lig_all)
tgt_ok <- intersect(targets_all, tgt_all)
stopifnot(length(lig_ok) > 0, length(tgt_ok) > 0)

# 3) 各リガンド上位ターゲット（N_per_lig まで）
N_per_lig <- 200
get_top_links <- function(mat, lig_set, tgt_set, n = 200){
  do.call(rbind, lapply(lig_set, function(lg){
    sc <- mat[lg, tgt_set]; ord <- order(sc, decreasing = TRUE)[seq_len(min(n, length(sc)))]
    data.frame(ligand = lg, target = tgt_set[ord],
               regulatory_potential = as.numeric(sc[ord]), stringsAsFactors = FALSE)
  }))
}
links_df <- get_top_links(mat, lig_ok, tgt_ok, N_per_lig)
stopifnot(nrow(links_df) > 0)

# 4) 起源グループ付与（Macro / Stromal / Both / Others）
lig_sat <- intersect(lig_all, genes_sat)
lig_str <- intersect(lig_all, genes_str)
links_df$group <- dplyr::case_when(
  links_df$ligand %in% lig_sat & links_df$ligand %in% lig_str ~ "Both",
  links_df$ligand %in% lig_sat                                 ~ "Pericyte",
  links_df$ligand %in% lig_str                                 ~ "Satellite",
  TRUE                                                         ~ "Others"
)
cat("Links by origin:\n"); print(table(links_df$group))



# 5) 受容体裏付け（receiver で発現する受容体の割合）
rec_frac_by_lig <- links_df %>% distinct(ligand) %>% rowwise() %>% mutate(frac = {
  recs <- unique(lr_network$to[lr_network$from == ligand])
  recs <- intersect(recs, rownames(muscle_sc[["RNA"]]))
  if (!length(recs)) 0 else mean(recs %in% genes_recv)
}) %>% ungroup() %>% { setNames(.$frac, .$ligand) }

# 6) 合成スコア → 上位 N_global を採用
rescale01 <- function(x){ x <- as.numeric(x); rng <- range(x, finite = TRUE); 
if (!is.finite(rng[2]-rng[1]) || (rng[2]-rng[1]==0)) return(rep(0,length(x)))
(x - rng[1])/(rng[2]-rng[1])
}
N_global <- 300
lig_act_norm <- setNames(rescale01(lig_ranked[[score_col]]), lig_ranked[[lig_col]])
plot_df <- links_df %>%
  mutate(
    rp_n   = rescale01(regulatory_potential),
    act_n  = ifelse(is.na(lig_act_norm[ligand]), 0, lig_act_norm[ligand]),
    rec_n  = ifelse(is.na(rec_frac_by_lig[ligand]), 0, rec_frac_by_lig[ligand]),
    score  = 0.55*rp_n + 0.30*act_n + 0.15*rec_n,
    ligandL = paste0("L:", ligand),
    targetT = paste0("T:", target)
  ) %>%
  arrange(desc(score)) %>%
  slice_head(n = N_global)
stopifnot(all(c("ligandL","targetT","group","score") %in% names(plot_df)))

# 7) セクタと色の準備（ここが今回の“透明化”の山場）
lig_sectors <- unique(plot_df$ligandL)
tgt_sectors <- unique(plot_df$targetT)
sectors <- c(lig_sectors, tgt_sectors)
grid_col <- c(
  setNames(rep("#2874c5", length(lig_sectors)), lig_sectors),  # Sender=blue
  setNames(rep("#d73027", length(tgt_sectors)), tgt_sectors)   # Receiver=red
)

# リンクの色：キーは必ず plot_df$group と一致させる
col_map <- c(Pericyte="#2874c5", Satellite="#45a045", Both="#7fff00", Others="#999999")
bad <- setdiff(unique(as.character(plot_df$group)), names(col_map))
stopifnot(length(bad) == 0)  # 不一致があればここで止める

cols_base <- col_map[match(as.character(plot_df$group), names(col_map))]
cols_base[is.na(cols_base)] <- "#999999"
edge_col  <- grDevices::adjustcolor(cols_base, alpha.f = 0.45)
stopifnot(length(edge_col) == nrow(plot_df))

# 線の太さ
link_lwd <- scales::rescale(plot_df$score, to = c(1.2, 6.0))


# レイアウト
op <- par(no.readonly = TRUE); on.exit({ par(op); circos.clear() }, add = TRUE)
layout(matrix(c(1,2), 1), widths = c(0.76, 0.24)); par(oma = c(0,0,3.4,0))

# 左：サークル本体
par(xpd = NA, mar = c(2,2,1,0.2))
circlize::circos.clear()
nL <- length(lig_sectors); nT <- length(tgt_sectors)
gap_after <- c(rep(1, nL-1), 6, rep(1, nT-1), 6)

# データと order の不一致を吸収
cats_x <- unique(c(plot_df$ligandL, plot_df$targetT))
if (!setequal(sectors, cats_x)) {
  sectors  <- intersect(sectors, cats_x)
  grid_col <- grid_col[names(grid_col) %in% sectors]
  nL <- sum(startsWith(sectors, "L:")); nT <- sum(startsWith(sectors, "T:"))
  gap_after <- c(rep(1, nL-1), 6, rep(1, nT-1), 6)
}

circlize::circos.par(
  start.degree = 90, gap.after = gap_after, track.margin = c(0.003, 0.003),
  cell.padding = c(0,0,0,0), canvas.xlim = c(-0.8,0.8), canvas.ylim = c(-0.8,0.8),
  points.overflow.warning = FALSE
)

df_links <- data.frame(from = plot_df$ligandL, to = plot_df$targetT, value = plot_df$score)
circlize::chordDiagram(
  x = df_links, order = sectors, grid.col = grid_col,
  col = edge_col, link.lwd = link_lwd, directional = 1,
  link.arr.type = "big.arrow", link.sort = TRUE, link.decreasing = FALSE,
  link.border = NA, annotationTrack = "grid", preAllocateTracks = 1
)

# 外周リング
RING_ALPHA <- 0.55
circlize::circos.trackPlotRegion(
  ylim = c(0,1), track.height = 0.022, bg.border = NA,
  panel.fun = function(x,y){
    nm <- circlize::get.cell.meta.data("sector.index")
    fill <- if (startsWith(nm,"L:")) "#2874c5" else "#d73027"
    circlize::circos.rect(circlize::CELL_META$xlim[1], 0, circlize::CELL_META$xlim[2], 1,
                          col = grDevices::adjustcolor(fill, alpha.f = RING_ALPHA), border = NA)
  }
)

# ラベル（送信/受信で色を変える）
keep_lig <- head(lig_ranked[[lig_col]], 200)
keep_tgt <- plot_df %>% group_by(target) %>% summarise(s = sum(score), .groups="drop") %>%
  arrange(desc(s)) %>% slice_head(n = 80) %>% pull(target)
keep_ligL <- paste0("L:", keep_lig); keep_tgtT <- paste0("T:", keep_tgt)

LAB_CEX <- 0.60; LABEL_OFFSET <- 0.10
circlize::circos.trackPlotRegion(
  track.index = 1, bg.border = NA,
  panel.fun = function(x,y){
    nm  <- circlize::get.cell.meta.data("sector.index")
    lab <- sub("^[LT]:","", nm)
    x0  <- circlize::CELL_META$xcenter
    y0  <- circlize::CELL_META$ylim[1] - circlize::mm_y(LABEL_OFFSET)
    if (startsWith(nm,"L:") && nm %in% keep_ligL)
      circlize::circos.text(x0,y0,lab,facing="clockwise",niceFacing=TRUE,adj=c(0,0.5),cex=LAB_CEX,col="#2874c5")
    if (startsWith(nm,"T:") && nm %in% keep_tgtT)
      circlize::circos.text(x0,y0,lab,facing="clockwise",niceFacing=TRUE,adj=c(0,0.5),cex=LAB_CEX,col="#d73027")
  }
)

# 右：凡例とタイトル
par(mar = c(2,1.5,1,2)); plot.new(); par(usr = c(0,1,0,1))
legend(0.05,0.95, legend = c("Pericyte_Vascular","Satellite","Both","Others"),
       fill = c("#2874c5","#45a045","#7fff00","#999999"),
       border = NA, bty = "n", cex = 1.05, title = "Link color (origin)")
legend(0.05,0.72, legend = c("Sender sectors","Receiver sectors"),
       fill = c("#2874c5","#d73027"), border = NA, bty = "n", cex = 1.05, title = "Sector color")
text(0.05, 0.56, "Composite score", adj = 0, font = 2, cex = 1.10)
yy <- seq(0.50, 0.35, length.out = 8); ww <- seq(1.2, 6.0, length.out = 8)
for (i in seq_along(yy)) segments(0.05, yy[i], 0.55, yy[i], lwd = ww[i], col = "#000000", lend = 1)
mtext("mdxD2 vs mdx: Pericyte_Vascular/Satellite \u2192 Myocyte/Myonucleus  |  Composite‑weighted chord diagram",
      outer = TRUE, side = 3, line = 1.4, font = 2, cex = 1.35)

# 1) 色分けグループの本当の内訳（描画対象のみ）
table(plot_df$group)

# 2) 合成スコアの群別寄与（どの群が太いリンクを占めるか）
tapply(plot_df$score, plot_df$group, sum)

# 3) 色マップのキー不一致ゼロ確認（透明化の典型原因）
setdiff(unique(as.character(plot_df$group)), names(col_map))

c("Il6","Tgfb1","Timp1","Hbegf","Vegfa","Icam1") -> check_ligs
subset(plot_df, ligand %in% check_ligs)[, c("ligand","group")] |> unique()



## 推奨の複合ヒートマップ（v2 系の標準プロット）
if (!is.null(nichenet_output$ligand_activity_target_heatmap)) {
  print(nichenet_output$ligand_activity_target_heatmap)
}




# ---- sender 組合せを一括で走らせて、各ケースの上位リガンドと妥当性を比較 ----
screen_senders <- function(receiver, sender_sets, cond_col, cond_oi, cond_ref,
                           assay="RNA", expr_pct=.expr_pct %||% 0.10, top_n=10) {
  stopifnot(exists("muscle_sc"), exists("ligand_target_matrix"),
            exists("lr_network"), exists("weighted_networks"))
  DefaultAssay(muscle_sc) <- assay
  out <- lapply(sender_sets, function(senders) {
    res <- tryCatch(
      nichenetr::nichenet_seuratobj_aggregate(
        seurat_obj = muscle_sc, receiver = receiver,
        condition_colname = cond_col, condition_oi = cond_oi, condition_reference = cond_ref,
        sender = senders,
        ligand_target_matrix = as.matrix(ligand_target_matrix),
        lr_network = lr_network, weighted_networks = weighted_networks,
        assay_oi = assay, expression_pct = expr_pct, verbose = FALSE
      ), error = function(e) NULL
    )
    if (is.null(res) || is.null(res$ligand_activities)) return(NULL)
    la <- res$ligand_activities
    score <- if ("aupr" %in% names(la)) "aupr" else "pearson"
    top <- head(la[order(-la[[score]]), c("test_ligand", score)], top_n)
    data.frame(sender_combo = paste(senders, collapse="+"),
               top_ligand = top$test_ligand, score = top[[score]],
               stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(out)
}

# 例：Myocyte をレシーバーに sender 候補を比較
# "MHC-low macrophage", # 0 low
# "MHC-high macrophage", # 2 high
# "Endothelial cell(EC)", # 3
# "Stromal", # 5
# "Immune-Stromal(Cd53+)", # 6 low
# "Myocyte/Myonucleus", # 11
# "Tenocyte/Tendon fibroblast", # 13
# "Neutrophil",  # 16
# "Dendritic cell(DC)", # 17 high
# "Satellite cell", # 19
# "Pericyte/Vascular SMC", # 20
# "Schwann cell", # 21
sender_grid <- list(
  c("Stromal"), c("MHC-high macrophage"),
  c("Endothelial cell(EC)"), c("Pericyte/Vascular SMC"),
  c("Stromal","MHC-high macrophage"),
  c("Stromal","Pericyte/Vascular SMC"),
  c("Endothelial cell(EC)","Tenocyte/Tendon fibroblast"),
  c("Endothelial cell(EC)","Immune-Stromal(Cd53+)"),
  c("Stromal","Tenocyte/Tendon fibroblast"),
  c("Stromal","Immune-Stromal(Cd53+)"),
  c("Stromal","MHC-high macrophage","Endothelial cell(EC)")
)
res_grid_myo <- screen_senders(
  receiver = "Myocyte/Myonucleus",
  sender_sets = sender_grid,
  cond_col = "genotype", cond_oi = "mdxD2", cond_ref = "mdx"
)
print(res_grid_myo)

res_grid_satellite <- screen_senders(
  receiver = "Satellite cell",
  sender_sets = sender_grid,
  cond_col = "genotype", cond_oi = "mdxD2", cond_ref = "mdx"
)
print(res_grid_satellite)


## === 追加ユーティリティ：サークル図の色（送信グループ）出現サマリ ===
## 目的：
##  - chord を描く前に、各色グループ（例：Macrophage / Stromal / Both / Shared）ごとの
##    本数・スコア比率を数値で確認する
##  - 「青リンクが 0」「緑が支配的」などを即時に判定

suppressPackageStartupMessages({library(dplyr); library(Matrix)})

chord_links_and_stats <- function(
    nichenet_output, muscle_sc, receiver_oi, senders,
    ligand_target_matrix, lr_network,
    expr_pct = 0.10, N_per_lig = 100, N_global = 300,
    group_labels = NULL   # 例: c("Macrophage","Stromal") / c("Endothelial","Peri_vascular")
){
  stopifnot(is.list(nichenet_output), inherits(muscle_sc, "Seurat"))
  Idents(muscle_sc) <- "celltype"
  
  # 送信クラスターごとの「発現遺伝子」
  get_genes <- function(idents, pct = 0.10){
    cells <- Seurat::WhichCells(muscle_sc, idents = idents)
    m <- tryCatch(
      SeuratObject::GetAssayData(muscle_sc, assay = "RNA", layer = "data"),
      error = function(e) SeuratObject::GetAssayData(muscle_sc, assay = "RNA")
    )
    sub <- m[, cells, drop = FALSE]
    rownames(m)[ Matrix::rowSums(sub > 0) / ncol(sub) >= pct ]
  }
  genes_list <- lapply(senders, get_genes, pct = expr_pct)
  names(genes_list) <- if (is.null(group_labels)) senders else group_labels
  
  # 受け手の発現遺伝子（受容体裏付け用）
  recv_cells <- Seurat::WhichCells(muscle_sc, idents = receiver_oi)
  m_all <- tryCatch(
    SeuratObject::GetAssayData(muscle_sc, assay = "RNA", layer = "data"),
    error = function(e) SeuratObject::GetAssayData(muscle_sc, assay = "RNA")
  )
  sub <- m_all[, recv_cells, drop = FALSE]
  genes_recv <- rownames(m_all)[ Matrix::rowSums(sub > 0) / ncol(sub) >= expr_pct ]
  
  # LTM の向き調整（行=ligand / 列=target になるように）
  mat <- as.matrix(ligand_target_matrix)
  lig_act <- nichenet_output$ligand_activities
  lig_col   <- intersect(c("test_ligand","ligand"), names(lig_act))[1]
  score_col <- if ("aupr" %in% names(lig_act)) "aupr" else "pearson"
  lig_ranked <- lig_act[order(-lig_act[[score_col]]), , drop = FALSE]
  ligs  <- lig_ranked[[lig_col]]
  if (sum(ligs %in% rownames(mat)) < sum(ligs %in% colnames(mat))) mat <- t(mat)
  
  lig_all <- rownames(mat); tgt_all <- colnames(mat)
  tgt_key <- intersect(c("top_targets","prioritized_targets"), names(nichenet_output))[1]
  targets_all <- if (length(tgt_key)) unique(unlist(nichenet_output[[tgt_key]])) else tgt_all
  lig_ok <- intersect(ligs, lig_all)
  tgt_ok <- intersect(targets_all, tgt_all)
  
  # 各リガンドの上位ターゲット（規定 N_per_lig 件）
  get_top_links <- function(mat, lig_set, tgt_set, n = 200){
    do.call(rbind, lapply(lig_set, function(lg){
      sc <- mat[lg, tgt_set]; ord <- order(sc, decreasing = TRUE)[seq_len(min(n, length(sc)))]
      data.frame(ligand = lg, target = tgt_set[ord],
                 regulatory_potential = as.numeric(sc[ord]), stringsAsFactors = FALSE)
    }))
  }
  links_df <- get_top_links(mat, lig_ok, tgt_ok, N_per_lig)
  
  # グループ付与：各リガンドがどの sender 群で「発現」しているか
  #  - 1群だけ → その群名
  #  - 2群以上 → "Shared"（= Both の一般化）
  #  - 0群     → "Not_in_senders"（凡例には出さないが統計には残す）
  membership <- sapply(links_df$ligand, function(lg){
    hit <- names(Filter(function(gs) lg %in% gs, genes_list))
    if (length(hit) == 0L) "Not_in_senders" else if (length(hit) == 1L) hit else "Shared"
  })
  links_df$group <- base::unname(membership)
  
  # 合成スコア（chord と同じ定義：RP 55%, LigAct 30%, Receptor裏付け 15%）
  rescale01 <- function(x){ x <- as.numeric(x); rng <- range(x, finite = TRUE);
  if (!is.finite(rng[2]-rng[1]) || (rng[2]-rng[1]==0)) return(rep(0,length(x)))
  (x-rng[1])/(rng[2]-rng[1]) }
  # 受容体裏付け（受け手で発現する受容体の割合）
  rec_frac_by_lig <- links_df %>% dplyr::distinct(ligand) %>% dplyr::rowwise() %>%
    dplyr::mutate(frac = {
      recs <- unique(lr_network$to[lr_network$from == ligand])
      recs <- intersect(recs, rownames(muscle_sc[["RNA"]]))
      if (!length(recs)) 0 else mean(recs %in% genes_recv)
    }) %>% dplyr::ungroup() %>% { stats::setNames(.$frac, .$ligand) }
  
  lig_norm <- stats::setNames(rescale01(lig_ranked[[score_col]]), lig_ranked[[lig_col]])
  plot_df <- links_df %>%
    mutate(rp_n   = rescale01(regulatory_potential),
           act_n  = ifelse(is.na(lig_norm[ligand]), 0, lig_norm[ligand]),
           rec_n  = ifelse(is.na(rec_frac_by_lig[ligand]), 0, rec_frac_by_lig[ligand]),
           score  = 0.55*rp_n + 0.30*act_n + 0.15*rec_n,
           ligandL = paste0("L:", ligand),
           targetT = paste0("T:", target)) %>%
    arrange(desc(score)) %>%
    slice_head(n = N_global)
  
  # グループ別サマリ（本数・合成スコアの総和と比率）
  sum_df <- plot_df %>%
    count(group, name = "n_edges") %>%
    mutate(score_sum  = tapply(plot_df$score, plot_df$group, sum)[as.character(group)],
           score_frac = score_sum / sum(score_sum)) %>%
    arrange(desc(score_sum))
  sum_df[is.na(sum_df)] <- 0
  
  present <- sum_df$group[sum_df$n_edges > 0]
  wanted  <- c(names(genes_list), "Shared")   # 期待する凡例キー
  miss    <- setdiff(wanted, present)
  
  list(plot_df = plot_df,
       summary = sum_df,
       present_groups = present,
       missing_groups = miss)
}



# "MHC-low macrophage", # 0 low
# "MHC-high macrophage", # 2 high
# "Endothelial cell(EC)", # 3
# "Stromal", # 5
# "Immune-Stromal(Cd53+)", # 6 low
# "Myocyte/Myonucleus", # 11
# "Tenocyte/Tendon fibroblast", # 13
# "Neutrophil",  # 16
# "Dendritic cell(DC)", # 17 high
# "Satellite cell", # 19
# "Pericyte/Vascular SMC", # 20
# "Schwann cell", # 21
{
  stats_MS <- chord_links_and_stats(
    nichenet_output, muscle_sc,
    receiver_oi = "Myocyte/Myonucleus",
    senders = c("Pericyte/Vascular SMC","Satellite cell"),
    group_labels = c("Pericyte/Vascular SMC","Satellite cell"),
    ligand_target_matrix = ligand_target_matrix,
    lr_network = lr_network,
    expr_pct = 0.10, N_per_lig = 100, N_global = 300
  )
  
  print(stats_MS$summary)
  if (length(stats_MS$missing_groups)) {
    message("欠落している色グループ: ", paste(stats_MS$missing_groups, collapse = ", "))
  } else {
    message("すべての色グループが出現しています。")
  }
}

### 1: Myocyte
# Peri/EC
# 1                Shared      70 14.597953  0.7850038
# 2 Pericyte/Vascular SMC      14  2.657834  0.1429248
# 3  Endothelial cell(EC)       7  1.340242  0.0720714

# Stromal/Macro
# 1  Shared      56  11.74200  0.6314251
# 2 Stromal      35   6.85403  0.3685749

# EC/Immune-Stromal
# 1  Endothelial cell(EC)      42  8.370836  0.4501411
# 2                Shared      35  7.567358  0.4069341
# 3 Immune-Stromal(Cd53+)      14  2.657834  0.1429248




### 2: Satellite cel
# 1                Shared      70 15.8002258 0.79805967
# 2 Pericyte/Vascular SMC      14  3.1828336 0.16076297
# 3  Endothelial cell(EC)       7  0.8152418 0.04117736

# Stromal/Macro
# 1  Shared      56  12.57927  0.6353713
# 2 Stromal      35   7.21903  0.3646287

# EC/Immune-Stromal
# 1  Endothelial cell(EC)      42  8.570836  0.4329077
# 2                Shared      35  8.044631  0.4063294
# 3 Immune-Stromal(Cd53+)      14  3.182834  0.1607630




## 構成比の群間比較（propeller）
# BiocManager が無ければ導入
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Bioconductor 版 speckle をインストール
BiocManager::install("speckle")

library(speckle)
packageVersion("speckle")

res_macro1 <- propeller(clusters = muscle_sc$celltype,
                 sample   = muscle_sc$orig.ident,
                 group    = muscle_sc$genotype,   # wt / mdx / d2
                 transform = "logit")

# まず列名を確認
names(res_macro1)

# MHC-high macrophage の行だけ抜き出す（行名にマッチさせる）
res_macro1[grep("MHC-high macrophage", rownames(res_macro1)),
           c("PropMean.WT","PropMean.mdx","PropMean.mdxD2",
             "Fstatistic","P.Value","FDR")]
#                     PropMean.WT PropMean.mdx PropMean.mdxD2 Fstatistic   P.Value       FDR
# MHC-high macrophage  0.03768032    0.1045952     0.09019823   1.553551 0.2114956 0.2307224


res_macro1[grep("MHC-low macrophage", rownames(res_macro1)),
           c("PropMean.WT","PropMean.mdx","PropMean.mdxD2",
             "Fstatistic","P.Value","FDR")]
#                    PropMean.WT PropMean.mdx PropMean.mdxD2 Fstatistic     P.Value        FDR
# MHC-low macrophage  0.07885002    0.1686256      0.4095588   5.368641 0.004660458 0.01398138



































## ==== Pseudobulk camera + fgsea pipeline (Seurat v5) =======================
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(SingleCellExperiment); library(scDblFinder)
  library(edgeR); library(limma)
  library(fgsea)         # GSEA
  suppressWarnings(library(msigdbr))  # オフライン時は不要
  library(data.table)
})

set.seed(123); options(stringsAsFactors = FALSE)

## 0) 入力（Seurat v5 オブジェクト：RNA assay / counts.* layers 必須）
obj_path <- "C:/DMD_project/rds/muscle_sc_after_dotplot.rds"
muscle_sc <- readRDS(obj_path)
stopifnot(inherits(muscle_sc, "Seurat"))
DefaultAssay(muscle_sc) <- "RNA"

## 1) scDblFinder を layer 単位で実行（先生の方針を踏襲）
layers_all <- Layers(muscle_sc[["RNA"]])
layers <- setdiff(grep("^counts(\\.|$)", layers_all, value = TRUE), "counts")
dbl_class <- setNames(rep(NA_character_, ncol(muscle_sc)), colnames(muscle_sc))
dbl_score <- setNames(rep(NA_real_,      ncol(muscle_sc)), colnames(muscle_sc))

for (lay in layers) {
  mat <- muscle_sc[["RNA"]]@layers[[lay]]
  if (is.null(mat) || ncol(mat) == 0) next
  cells_lay <- Cells(muscle_sc, layer = lay)
  sce <- SingleCellExperiment(list(counts = mat))
  sce <- scDblFinder(sce)
  dbl_class[cells_lay] <- as.character(colData(sce)$scDblFinder.class)
  dbl_score[cells_lay] <- as.numeric(  colData(sce)$scDblFinder.score)
}
dbl_class[is.na(dbl_class)] <- "singlet"
muscle_sc$scDblFinder.class <- dbl_class
muscle_sc$scDblFinder.score <- dbl_score

## 2) RBC + doublet 除外
stopifnot("celltype" %in% colnames(muscle_sc@meta.data))
muscle_sc <- subset(muscle_sc,
                    subset = (celltype != "Erythrocyte") &
                      (scDblFinder.class != "doublet"))
muscle_sc$celltype <- droplevels(muscle_sc$celltype)
Idents(muscle_sc) <- "celltype"

## 3) 擬似バルク（celltype × orig.ident）
agg <- AggregateExpression(
  muscle_sc, group.by = c("celltype","orig.ident"),
  assays = "RNA", slot = "counts"
)
pb_counts <- agg$RNA       # genes x (celltype_sample) の行列
stopifnot(is.matrix(pb_counts) || inherits(pb_counts, "dgCMatrix"))

## 4) ユーティリティ
.parse_group <- function(labels) {
  # ラベル末尾の wt / mdx / d2 を群として抽出（例: "Myocyte_wt1" → "wt"）
  grp <- sub("^(wt|mdx|d2).*", "\\1", labels)
  factor(grp, levels = c("wt","mdx","d2"))
}
.get_cols_by_ct <- function(ct) {
  cols <- colnames(pb_counts)[startsWith(colnames(pb_counts), paste0(ct, "_"))]
  if (!length(cols)) stop(sprintf("列がありません：%s", ct))
  cols
}

## ==== Pathway utility: msigdbr / GMT の揺らぎを吸収して名前付きリストで返す ====
suppressPackageStartupMessages({
  library(msigdbr)    # オンライン時に使用（なければフォールバック）
  library(fgsea)      # gmtPathways を使う
  library(edgeR); library(limma); library(data.table)
})

get_pathways_safe <- function(collections = c("REACTOME","GO:BP","CP:KEGG"),
                              species = "Mus musculus",
                              min_size = 10, max_size = 500,
                              gmt_files = list()) {
  ## collections は "REACTOME", "GO:BP", "HALLMARK", "CP:KEGG" などの文字列
  collections <- toupper(collections)
  
  build_list_from_tibble <- function(tb) {
    stopifnot(all(c("gs_name","gene_symbol") %in% names(tb)))
    split(tb$gene_symbol, tb$gs_name)
  }
  
  out_list <- list()
  
  for (coll in collections) {
    ## 1) まず msigdbr を試す
    msig_ok <- FALSE; lst <- NULL
    cat(sprintf("[get_pathways] try msigdbr for %s ...\n", coll))
    tb <- try({
      ## coll を C2/C5 へ写像
      if (coll == "REACTOME") {
        # db_species は新しめのmsigdbrで有効。古い環境では無視される→tryCatchで二段構え
        tryCatch(
          msigdbr::msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME",
                           db_species = "MM"),
          error = function(e) msigdbr::msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME")
        )
      } else if (coll %in% c("GO:BP","GOBP","GO_BP")) {
        tryCatch(
          msigdbr::msigdbr(species = species, category = "C5", subcategory = "GO:BP",
                           db_species = "MM"),
          error = function(e) msigdbr::msigdbr(species = species, category = "C5", subcategory = "GO:BP")
        )
      } else if (coll %in% c("CP:KEGG","KEGG")) {
        # ライセンス都合でMSigDB側に無いことがある→失敗したらGMTへ
        tryCatch(
          msigdbr::msigdbr(species = species, category = "C2", subcategory = "CP:KEGG",
                           db_species = "MM"),
          error = function(e) NULL
        )
      } else if (coll == "HALLMARK") {
        tryCatch(
          msigdbr::msigdbr(species = species, category = "H",
                           db_species = "MM"),
          error = function(e) msigdbr::msigdbr(species = species, category = "H")
        )
      } else {
        stop(sprintf("Unsupported collection: %s", coll))
      }
    }, silent = TRUE)
    
    if (inherits(tb, "data.frame") && nrow(tb) > 0L) {
      lst <- build_list_from_tibble(tb)
      msig_ok <- TRUE
    }
    
    ## 2) 失敗 or 空なら GMT にフォールバック
    if (!msig_ok) {
      gmt_path <- gmt_files[[coll]]
      if (is.null(gmt_path)) {
        cat(sprintf("[get_pathways] msigdbr failed and no GMT for %s. Skip.\n", coll))
        next
      }
      cat(sprintf("[get_pathways] use GMT for %s : %s\n", coll, gmt_path))
      lst <- fgsea::gmtPathways(gmt_path)  # これは最初から「名前付きリスト」を返す
    }
    
    ## 3) サイズでフィルタし、重複を除去
    lst <- lst[!duplicated(names(lst))]
    keep <- vapply(lst, length, integer(1))
    lst <- lst[keep >= min_size & keep <= max_size]
    cat(sprintf("[get_pathways] %s -> %d gene sets kept (size %d–%d)\n",
                coll, length(lst), min_size, max_size))
    out_list <- c(out_list, lst)
  }
  
  out_list
}

## ==== 擬似バルク: celltypeごとの limma-voom + camera & fgsea（先生の設計踏襲） ====
## 前提: pb_counts は genes x samples のカウント行列、
##       列名は "Celltype_sampleLabel" 形式（例: "Myocyte/Myonucleus_mdx1"）
## ==== 前提：pb_counts は genes x samples の擬似バルク行列 ====
## 必要パッケージ
suppressPackageStartupMessages({
  library(edgeR); library(limma)
  library(fgsea); library(data.table)
  library(BiocParallel)
})

## 並列は止めておく（環境依存の警告・エラーを封じる）
BiocParallel::register(BiocParallel::SerialParam())
data.table::setDTthreads(1)

## 先生の前提を再掲
stopifnot(exists("pb_counts"))  # genes x samples
.get_cols_by_ct <- function(ct) {
  pfx <- paste0(ct, "_")
  cols <- colnames(pb_counts)[startsWith(colnames(pb_counts), pfx)]
  if (!length(cols)) stop(sprintf("No columns for celltype='%s'", ct))
  cols
}
.parse_group <- function(labels) factor(sub("^(wt|mdx|d2).*","\\1", labels), levels = c("wt","mdx","d2"))

run_ct_tests_safe <- function(ct,
                              contrast = c("mdx_vs_wt","d2_vs_mdx","d2_vs_wt"),
                              pathways,                # named list: pathway -> gene symbols
                              min_gs = 10, max_gs = 500,
                              prefer_multilevel = TRUE) {
  
  contrast <- match.arg(contrast)
  cols <- .get_cols_by_ct(ct)
  mat  <- as.matrix(pb_counts[, cols, drop = FALSE])
  
  labels <- sub(".*_", "", cols)
  grp    <- .parse_group(labels)
  
  y <- DGEList(counts = mat)
  keep <- filterByExpr(y, group = grp)
  y <- y[keep,, keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMMwsp")
  
  X <- model.matrix(~0 + grp)
  v <- voomWithQualityWeights(y, X, plot = FALSE)
  
  vfit <- lmFit(v, X)
  Kall <- makeContrasts(
    mdx_vs_wt = grpmdx - grpwt,
    d2_vs_mdx = grpd2  - grpmdx,
    d2_vs_wt  = grpd2  - grpwt, levels = X
  )
  K <- switch(contrast,
              mdx_vs_wt = Kall[, "mdx_vs_wt"],
              d2_vs_mdx = Kall[, "d2_vs_mdx"],
              d2_vs_wt  = Kall[, "d2_vs_wt"])
  vfit <- contrasts.fit(vfit, K)
  vfit <- eBayes(vfit)
  
  genes_in <- rownames(vfit)
  
  ## ---- ここが肝: camera と fgsea で渡す形式を分ける ----
  ## fgsea 用：遺伝子名リスト（サイズ制約を同時にかける）
  pw_list <- lapply(pathways, function(gs) intersect(gs, genes_in))
  sz <- vapply(pw_list, length, 1L)
  pw_list <- pw_list[sz >= min_gs & sz <= max_gs]
  
  ## camera 用：上と同じ集合を「行インデックス」に変換
  idx <- lapply(pw_list, function(gs) match(gs, genes_in))
  
  ## ---- camera（相関補正あり） ----
  cam <- camera(v, index = idx, design = X, contrast = K)
  cam$FDR <- p.adjust(cam$PValue, "BH")
  cam <- cam[order(cam$FDR, cam$PValue), , drop = FALSE]
  
  ## ---- fgsea（Multilevel → Simple フォールバック） ----
  ranks <- setNames(vfit$t[,1], genes_in)
  fg <- tryCatch({
    if (prefer_multilevel) {
      suppressWarnings(fgseaMultilevel(pathways = pw_list, stats = ranks,
                                       minSize = min_gs, maxSize = max_gs))
    } else {
      stop("skip")
    }
  }, error = function(e) {
    fgseaSimple(pathways = pw_list, stats = ranks,
                nperm = 20000, minSize = min_gs, maxSize = max_gs)
  })
  fg <- as.data.table(fg)[order(padj, pval)]
  
  list(
    celltype  = ct, contrast = contrast,
    n_samples = table(grp),
    topTable  = topTable(vfit, number = Inf, sort.by = "P"),
    camera    = cam,
    fgsea     = fg,
    sets_kept = names(pw_list)
  )
}



## ==== 実行例（Myocyte/Satellite, d2_vs_mdx） ====
## 1) 経路集合の用意：オンラインなら msigdbr、オフラインなら gmt を指定
pathways_all <- get_pathways_safe(
  collections = c("REACTOME","GO:BP"),           # 必要なら "CP:KEGG","HALLMARK" を追加
  species     = "Mus musculus",
  min_size    = 10, max_size = 500,
  gmt_files   = list(
    # 先生のローカルに置いたGMTがあればパスを指定（無ければNULLのままでOK）
    # "REACTOME" = "data/c2.cp.reactome.mouse.symbols.gmt",
    # "GO:BP"    = "data/c5.go.bp.mouse.symbols.gmt",
    # "CP:KEGG"  = "data/c2.cp.kegg.symbols.gmt"
    # 例示なのでコメントアウト
  )
)

## パスウェイ集合（先生の get_pathways_safe の出力をそのまま渡せます）
## pathways_all は “名前付きリスト: パスウェイ名 -> 遺伝子名ベクター”
stopifnot(exists("pathways_all"))

## まずは Myocyte / d2_vs_mdx
res_myocyte <- run_ct_tests_safe("Myocyte/Myonucleus", "d2_vs_mdx", pathways_all)

## フォーカス名は、実際に残っている集合に合わせてからフィルタ
focus_names <- grep(paste(c(
  "OXIDATIVE_PHOSPHORYLATION", "RESPIRATORY_ELECTRON_TRANSPORT",
  "TCA", "FATTY_ACID.*OXIDATION", "MITOCHONDRIAL",
  "EXTRACELLULAR_MATRIX_ORGANIZATION", "COLLAGEN", "MATRIX_METALLOPROTEINASE",
  "PLATELET.*DEGRANULATION", "TGF.*BETA", "ADIPOGENESIS", "PPAR"
), collapse="|"), names(pathways_all), value = TRUE)

focus_kept <- intersect(focus_names, res_myocyte$sets_kept)

## camera と fgsea の整合チェック
res_myocyte$camera[rownames(res_myocyte$camera) %in% focus_kept, ][1:21, ]
#                                                                       NGenes Direction         FDR
# GOBP_OXIDATIVE_PHOSPHORYLATION                                           115      Down  0.00000000000000424918
# REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT          183      Down  0.00000000000038058772
# REACTOME_RESPIRATORY_ELECTRON_TRANSPORT                                  120      Down  0.00000000000103919154
# GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE                  47      Down  0.00000000000678359936
# REACTOME_COLLAGEN_FORMATION                                               15        Up  0.00000005638168050723
# REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES     12        Up  0.00000016457457693875
# REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES                      11        Up  0.00000158611099339044
# REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION                                47        Up  0.00000191172142942930
# GOBP_COLLAGEN_FIBRIL_ORGANIZATION                                         10        Up  0.00000416549211786893
# GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_UBIQUINOL_TO_CYTOCHROME_C           12      Down  0.00001296767256598388
# REACTOME_COLLAGEN_DEGRADATION                                             11        Up  0.00001507502563847225
# REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION                                65      Down  0.00013973296894387145
# GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN              18      Down  0.00015084533399006370
# REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE                                      28      Down  0.00024982087320073171
# GOBP_COLLAGEN_METABOLIC_PROCESS                                           22        Up  0.00031306588111104335
# REACTOME_MITOCHONDRIAL_BIOGENESIS                                         46      Down  0.00054825620519424788
# GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY                     67      Down  0.00237718667779414361
# GOBP_INNER_MITOCHONDRIAL_MEMBRANE_ORGANIZATION                            24      Down  0.01432916502223899671
# REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT                                     50      Down  0.01626450318504392770
# REACTOME_MATURATION_OF_TCA_ENZYMES_AND_REGULATION_OF_TCA_CYCLE            17      Down  0.03164560367175046213
# GOBP_MITOCHONDRIAL_GENE_EXPRESSION                                        86      Down  0.04794620212505335682


res_myocyte$fgsea [pathway %in% focus_kept][order(padj)][1:22]
#                                                                   pathway                                       
# <char>                                      <num>
# 1:       REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT    -3.284821
# 2:                                        GOBP_OXIDATIVE_PHOSPHORYLATION    -3.505842
# 3:                               REACTOME_RESPIRATORY_ELECTRON_TRANSPORT    -3.406940
# 4:              GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE    -3.270048
# 5:                            REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION     2.397778 Cd47, Col1a1, Col1a2, Fn1, Spp1, Timp1
# 6:                            REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION    -2.402719
# 7:          GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN    -2.492191
# 8: REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES     2.203094
# 9:                                       GOBP_COLLAGEN_METABOLIC_PROCESS     2.241200    
# 10:                 GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY   -2.219578
# 11:       GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_UBIQUINOL_TO_CYTOCHROME_C   -2.348689
# 12:                                     REACTOME_MITOCHONDRIAL_BIOGENESIS   -2.307869
# 13:                                           REACTOME_COLLAGEN_FORMATION    2.196253    
# 14:                                         REACTOME_COLLAGEN_DEGRADATION    2.102196    
# 15:                                  REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE   -2.267337
# 16:                  REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES    1.994242    
# 17:                                     GOBP_COLLAGEN_FIBRIL_ORGANIZATION    1.947975    
# 18:                        GOBP_INNER_MITOCHONDRIAL_MEMBRANE_ORGANIZATION   -2.100871
# 19:                                 REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT   -1.906761
# 20:                                    GOBP_MITOCHONDRIAL_GENE_EXPRESSION   -1.657863
# 21:                                    REACTOME_MITOCHONDRIAL_TRANSLATION   -1.667341
# 22:                                        GOBP_MITOCHONDRIAL_TRANSLATION   -1.608181

library(data.table); library(dplyr)

ligands_interest <- c("Spp1","Fn1","Lamb2","Col1a1","Col1a2","Col18a1","Timp1",
                      "Cd47")
lig_low <- tolower(ligands_interest)

fg <- as.data.table(res_myocyte$fgsea)

## ① NicheNet で検証したい経路だけに限定（どちらか一方を使う）
target_ids <- unique(hit_table$fg_best)            # ← hit_table の “fg_best” に絞る
# または:
# target_ids <- intersect(focus_kept, fg$pathway)  # ← focus_kept に絞る

## ② その中で「有意（padj≤0.05）かつ NES>0」だけを確認
fg_pos <- fg[padj <= 0.05 & NES > 0 & pathway %in% target_ids]

## ③ leading‑edge にリガンドが含まれているかを機械的に確認
le_long <- fg_pos[, .(gene = as.character(unlist(leadingEdge))),
                  by = .(pathway, NES, padj)]
le_long[, gene_low := tolower(gene)]

hits_long <- le_long[gene_low %chin% lig_low]
hits_long[, ligand := ligands_interest[match(gene_low, lig_low)]]

## ④ 経路ごとのリガンドヒット要約（padj 昇順）
hits_by_pathway_focus <- hits_long[, .(
  n_found = uniqueN(ligand),
  found_ligands = paste(sort(unique(ligand)), collapse = ", ")
), by = .(pathway, NES, padj)][order(padj)]

print(hits_by_pathway_focus)

## ⑤ “余計な経路” が混ざっていないかの検算（空なら OK）
setdiff(unique(hits_by_pathway_focus$pathway), target_ids)



## ============== 文字整形 & 可視化：衝突回避の安全版 ==========================
suppressPackageStartupMessages({
  library(ggplot2); library(ggridges); library(dplyr); library(stringr); library(data.table)
})

## 1) fgsea 上位20（先生が貼ってくださった表と同じ集合・順序）をそのまま使う
ids_to_show <- as.character(
  res_myocyte$fgsea[pathway %in% focus_kept][order(padj)][1:20, ]$pathway
)

## 2) ランク（limma の t 値）— camera/fgsea と同じもの
ranks_named <- setNames(res_myocyte$topTable$t, rownames(res_myocyte$topTable))
ranks_named <- sort(ranks_named[!is.na(ranks_named)], decreasing = TRUE)

## 3) running ES を“再構成”する関数（p=1；fgsea と同じ重み）
calc_running_es <- function(stats, genes, gsea_p = 1) {
  stats <- stats[!is.na(stats)]
  stats <- sort(stats, decreasing = TRUE)
  in_set <- names(stats) %in% genes
  Nh <- sum(in_set); N <- length(stats); Nm <- N - Nh
  if (Nh == 0L || Nh == N) return(rep(0, N))
  w <- abs(stats)^gsea_p
  Phit  <- cumsum(ifelse(in_set, w, 0)) / sum(w[in_set])
  Pmiss <- cumsum(ifelse(in_set, 0, 1)) / Nm
  Phit - Pmiss
}

## 4) ridge 用データ（各経路の running ES を抽出；描画軽量化のため間引き）
build_ridge_df <- function(stats, gene_sets, ids, sample_each = 1500) {
  ids <- intersect(ids, names(gene_sets))
  out <- vector("list", length(ids)); k <- 0L
  for (id in ids) {
    es <- calc_running_es(stats, intersect(gene_sets[[id]], names(stats)), gsea_p = 1)
    if (length(es) > sample_each) {
      take <- unique(round(seq(1, length(es), length.out = sample_each)))
      es <- es[take]
    }
    k <- k + 1L
    out[[k]] <- data.frame(pathway = id, es = es, stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(out)
}

ridge_df <- build_ridge_df(ranks_named, pathways_all, ids_to_show, sample_each = 1500)

## 5) 見やすいラベルと、色に使う“先生の数値そのもの”
format_pathway_labels_clean <- function(pathways, wrap_width = 36){
  lab <- gsub("^REACTOME_|^GOBP_", "", pathways)
  lab <- gsub("_", " ", lab)
  lab <- stringr::str_to_sentence(lab)
  lab <- gsub("Tca","TCA",lab,fixed=TRUE)
  lab <- gsub("Ppar","PPAR",lab,fixed=TRUE)
  lab <- gsub("Ecm","ECM",lab,fixed=TRUE)
  lab <- gsub("Mmp","MMP",lab,fixed=TRUE)
  stringr::str_wrap(lab, width = wrap_width)
}
pretty_map <- setNames(format_pathway_labels_clean(unique(ridge_df$pathway), 36),
                       unique(ridge_df$pathway))
ridge_df$label <- factor(pretty_map[ridge_df$pathway], levels = pretty_map[ids_to_show])

camera_fdr_map <- setNames(res_myocyte$camera$FDR, rownames(res_myocyte$camera))
fgsea_padj_map <- setNames(res_myocyte$fgsea$padj,  res_myocyte$fgsea$pathway)

ridge_df$camera_fdr <- camera_fdr_map[ridge_df$pathway]
ridge_df$fgsea_padj <- fgsea_padj_map[ridge_df$pathway]

## 6) 描画関数（“色=先生の FDR/padj そのもの”、曲線=既存ランクからの再構成）
plot_ridge_constantfill <- function(df, fill_col, title_text, legend_title) {
  ggplot(df, aes(x = es, y = label, group = label, fill = !!as.name(fill_col))) +
    ggridges::geom_density_ridges(scale = 1.2, rel_min_height = 0.01,
                                  size = 0.2, color = "white") +
    scale_x_continuous(name = "Running enrichment score", breaks = seq(-2, 2, 1)) +
    scale_fill_gradient(low = "#b2182b", high = "#2166ac", name = legend_title) +
    labs(y = NULL, title = title_text) +
    theme_ridges(font_size = 11) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
}

plot_ridge_camera <- plot_ridge_constantfill(
  ridge_df, "camera_fdr",
  "Ridge plot（camera／color = camera FDR；curve = rank）", "camera FDR"
)
plot_ridge_fgsea <- plot_ridge_constantfill(
  ridge_df, "fgsea_padj",
  "Ridge plot（fgsea／color = fgsea padj；curve = rank）", "fgsea padj"
)

## 7) —— 画面に順に描画 ——
print(plot_ridge_camera)
print(plot_ridge_fgsea)

## 8) 〈任意〉整合の数式的チェック：fgsea の ES と“再構成 ES”が一致することの確認
check_es <- sapply(ids_to_show, function(id){
  es <- calc_running_es(ranks_named, intersect(pathways_all[[id]], names(ranks_named)))
  es_star <- if (abs(max(es)) >= abs(min(es))) max(es) else min(es)  # fgsea の ES と同じ定義
  c(reconstructed_ES = es_star,
    fgsea_ES = res_myocyte$fgsea[pathway == id, ES][1])
})
print(t(check_es))  # 2列（reconstructed_ES と fgsea_ES）がほぼ一致するはずです


## ========== FIX: NicheNet 整合・経路ヒットの要約（再計算なし） ==================
check_ligand_pathway_hits <- function(res_list,
                                      ligand_groups = list(
                                        ECM      = c("Spp1","Fn1","Lamb2","Col1a1","Col1a2","Col18a1","Timp1","Cd47"),
                                        Platelet = c("Spp1","Fn1","Timp1")
                                      ),
                                      pathway_patterns = list(
                                        ECM_org   = "EXTRACELLULAR_MATRIX_ORGANIZATION",
                                        Col_form  = "COLLAGEN_FORMATION",
                                        Col_asm   = "ASSEMBLY_OF_COLLAGEN_FIBRILS",
                                        Col_meta  = "COLLAGEN_METABOLIC_PROCESS",
                                        Col_bios  = "COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES",
                                        Col_deg   = "COLLAGEN_DEGRADATION",
                                        Platelet  = "PLATELET_(DE|AC)GRANULATION|PLATELET_ACTIVATION"
                                      ),
                                      fdr_cut  = 0.05,
                                      padj_cut = 0.05) {
  ## camera / fgsea を data.frame として扱い、pathway 列を明示的に用意
  camera_df <- as.data.frame(res_list$camera, stringsAsFactors = FALSE)
  camera_df$pathway <- rownames(res_list$camera)
  rownames(camera_df) <- NULL
  
  fgsea_df <- as.data.frame(res_list$fgsea, stringsAsFactors = FALSE)
  fgsea_df$pathway <- as.character(fgsea_df$pathway)
  
  pull_one <- function(pat) {
    cam_rows <- camera_df[ base::grepl(pat, camera_df$pathway, ignore.case = TRUE), , drop = FALSE]
    fg_rows  <-  fgsea_df[ base::grepl(pat,  fgsea_df$pathway,  ignore.case = TRUE), , drop = FALSE]
    
    if (!nrow(cam_rows) && !nrow(fg_rows)) return(NULL)
    
    cam_FDR_min   <- if (nrow(cam_rows)) min(cam_rows$FDR, na.rm = TRUE) else NA_real_
    cam_dir_atmin <- if (nrow(cam_rows)) cam_rows$Direction[ which.min(cam_rows$FDR) ] else NA_character_
    cam_best_name <- if (nrow(cam_rows)) cam_rows$pathway[   which.min(cam_rows$FDR) ] else NA_character_
    cam_hit_up    <- nrow(cam_rows) > 0 && any(cam_rows$FDR <= fdr_cut & cam_rows$Direction == "Up", na.rm = TRUE)
    
    fg_padj_min   <- if (nrow(fg_rows)) min(fg_rows$padj, na.rm = TRUE) else NA_real_
    fg_idx_min    <- if (nrow(fg_rows)) which.min(fg_rows$padj) else integer(0)
    fg_NES_atmin  <- if (length(fg_idx_min)) fg_rows$NES[fg_idx_min] else NA_real_
    fg_best_name  <- if (length(fg_idx_min)) fg_rows$pathway[fg_idx_min] else NA_character_
    fg_hit_up     <- nrow(fg_rows) > 0 && any(fg_rows$padj <= padj_cut & fg_rows$NES > 0, na.rm = TRUE)
    
    data.frame(
      pattern      = pat,
      cam_hit_up   = cam_hit_up,
      cam_FDR_min  = cam_FDR_min,
      cam_best     = cam_best_name,
      cam_dir_best = cam_dir_atmin,
      fg_hit_up    = fg_hit_up,
      fg_padj_min  = fg_padj_min,
      fg_NES_best  = fg_NES_atmin,
      fg_best      = fg_best_name,
      ligands_note = paste(unique(unlist(ligand_groups)), collapse = ", "),
      stringsAsFactors = FALSE
    )
  }
  
  rows <- lapply(pathway_patterns, pull_one)
  hits <- do.call(rbind, rows)
  hits[order(hits$cam_FDR_min, hits$fg_padj_min), ]
}
## 実行例：
hit_table <- check_ligand_pathway_hits(res_myocyte)
print(hit_table)






suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(stringr)})

## --- pretty label (英語・略語整形) -------------------------------------------
format_pathway_id <- function(ids){
  labs <- gsub("^REACTOME_|^GOBP_", "", ids)
  labs <- gsub("_", " ", labs, fixed = TRUE)
  labs <- str_to_sentence(labs)
  labs <- gsub("Ecm", "ECM", labs, fixed = TRUE)
  labs <- gsub("Tca", "TCA", labs, fixed = TRUE)
  labs <- gsub("Ppar","PPAR", labs, fixed = TRUE)
  labs <- gsub("Mmp","MMP", labs, fixed = TRUE)
  labs
}

## --- fgsea の最良ヒットだけを利用（再計算なし） -------------------------------
fg_bar_df <- hit_table %>%
  transmute(pathway_id = fg_best,
            FDR       = fg_padj_min,
            NES       = fg_NES_best) %>%
  filter(!is.na(pathway_id), is.finite(FDR), is.finite(NES)) %>%
  distinct(pathway_id, .keep_all = TRUE) %>%                   # 念のため重複除去
  arrange(FDR)

## ラベルと並び（FDR 降順 = -log10(FDR) が大きい順）
fg_bar_df$pathway_label <- factor(
  format_pathway_id(fg_bar_df$pathway_id),
  levels = format_pathway_id(fg_bar_df$pathway_id[order(-fg_bar_df$FDR)])
)


## カラースケールは NES の正負を見せる二色グラデーション
nes_lim <- max(abs(fg_bar_df$NES))
p_fgsea_bar <- ggplot(fg_bar_df, aes(x = -log10(FDR), y = pathway_label, fill = NES)) +
  geom_col(width = 0.7) +
  scale_fill_gradient2(limits = c(1.8, nes_lim),
                       low = "white", mid = "#2166AC", high = "#B2182B",
                       midpoint = 2.1, name = "NES") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  labs(title = "fgsea: pathway hits (Myocyte mdxd2 vs mdx)",
       x = "-log10(FDR) (adjusted p)", y = "Pathway") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 9),
        legend.position = "right")

print(p_fgsea_bar)
