## Script 2: scRNA-seq analysis

## 0) Environment setup (validated on R 4.5.1)----------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("glmGamPoi","edgeR","scDblFinder"), ask = FALSE, update = FALSE)

pkgs <- c("Seurat","harmony","Matrix","lme4","remotes")
install.packages(setdiff(pkgs, rownames(installed.packages())), Ncpus = 2)

# Install Rcpp-related packages first
if (!"Rcpp" %in% rownames(installed.packages())) install.packages("Rcpp", type = "binary")
if (!"RcppArmadillo" %in% rownames(installed.packages())) install.packages("RcppArmadillo", type = "binary")

update.packages(ask = FALSE, checkBuilt = TRUE)

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(edgeR); library(lme4); library(Matrix)
  library(SingleCellExperiment); library(scDblFinder)
})

set.seed(123)


## 1) Logging and options-------------------------------------------------------

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
    cat("\n!! Unhandled error occurred. Please check the log.\n")
    traceback(2)
    .on.exit_sinks()
  }
)
message("[START] ", Sys.time())
options(bitmapType = "cairo")   # Workaround for the Windows graphics device


## 2) Read 10X-compatible mtx files (core structure preserved)------------------

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


## 3) QC (mitochondrial percentage, UMI count, and detected gene count)---------

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


## 4) Doublet removal (scDblFinder; per sample)-------------------------------

## --- Count cells before and after doublet removal and calculate the number removed ---
# Cell counts before removal (after QC)
pre_n  <- sapply(objs, function(o) length(Cells(o)))   # Store the cell counts before removal

# Doublet removal (existing code)
objs <- lapply(objs, run_scDblFinder)

post_n <- sapply(objs, function(o) length(Cells(o)))   # Cell counts after removal
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


## 5) SCTransform (per sample) and lightweight integration (rpca)---------------

objs <- lapply(
  objs,
  SCTransform,
  vars.to.regress = "percent.mt",
  method = if ("glmGamPoi" %in% rownames(installed.packages())) "glmGamPoi" else "poisson",
  verbose = TRUE
)

var.features <- SelectIntegrationFeatures(objs, nfeatures = 5000)

# Restrict additional genes to those shared across all samples
common_all   <- Reduce(intersect, lapply(objs, function(o) rownames(o)))
features_use <- unique(c(var.features, intersect(mhc_core, common_all)))

# Consistency check
stopifnot(all(sapply(objs, function(o) all(features_use %in% rownames(o)))))

objs <- PrepSCTIntegration(objs, anchor.features = features_use)

# rpca Preprocessing
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


## 6) Dimensional reduction and clustering--------------------------------------

DefaultAssay(muscle_sc) <- "integrated"
muscle_sc <- RunPCA(muscle_sc, npcs = 50, verbose = TRUE)
muscle_sc <- RunUMAP(muscle_sc, dims = 1:30)
muscle_sc <- FindNeighbors(muscle_sc, dims = 1:30)
muscle_sc <- FindClusters(muscle_sc, resolution = 0.4)
# Maximum modularity in 10 random starts: 0.9466
# Number of communities: 23

## Metadata used for visualization
muscle_sc$geno <- factor(sub("^(wt|mdx|d2).*","\\1", muscle_sc$orig.ident),
                         levels = c("wt","mdx","d2"))

## If annotation is not yet available, use cluster IDs as provisional celltype labels
if (!"celltype" %in% colnames(muscle_sc@meta.data)) {
  muscle_sc$celltype <- Idents(muscle_sc)
}

## Save
#rds_dir <- "C:/DMD_project/rds"
#dir.create(rds_dir, showWarnings = FALSE)
#saveRDS(muscle_sc, file.path(rds_dir, "muscle_sc_proto.rds"), compress = "gzip")

muscle_proto <- readRDS("C:/DMD_project/rds/muscle_sc_proto.rds")


## Minimal plotting (print only, to avoid device-related issues)
print(DimPlot(muscle_proto, group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(DimPlot(muscle_proto, group.by = "orig.ident"))
print(DimPlot(muscle_proto, group.by = "celltype", split.by = "geno"))

## --- Count the number of cluster 22 cells in WT (minimal code) ---
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

# Reload for assignment lookup
#muscle_sc <- readRDS("C:/DMD_project/rds/muscle_sc_proto.rds")
#muscle_origin <- readRDS("C:/DMD_project/rds/muscle_sc_proto.rds")
muscle_origin$cluster <- Idents(muscle_origin) 
DimPlot(muscle_origin, reduction = "umap", label = TRUE)

muscle_sc$cluster <- Idents(muscle_sc)        # Retained for downstream reference
DimPlot(muscle_sc, reduction = "umap", label = TRUE)

## --- Check iLISI as a metric of batch mixing
# Compute iLISI in PCA space (30 dimensions)
X <- Embeddings(muscle_sc, "pca")[, 1:30, drop = FALSE]
ilisi <- compute_lisi(X, muscle_sc@meta.data, label_colnames = "orig.ident")
summary(ilisi$orig.ident)  # 2.447


## --- Prerequisite: Seurat object muscle_sc is present in memory ---
stopifnot("Seurat" %in% class(muscle_origin))
stopifnot("Seurat" %in% class(muscle_sc))

## 7) Ensure that the assay and UMAP are available ---------------------------------
assay_names <- names(Assays(muscle_origin))

DefaultAssay(muscle_origin) <- if ("SCT" %in% assay_names) "SCT" else if ("RNA" %in% assay_names) "RNA" else assay_names[1]

# Set SCT as the default assay first
DefaultAssay(muscle_origin) <- "SCT"

# Automatically select an available layer in SCT (use counts if data is absent)
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

## 8) Utilities: check gene availability and generate blended plots ---------------
check_genes <- function(obj, genes) {
  miss <- setdiff(genes, rownames(obj))
  if (length(miss)) {
    warning(sprintf("Genes not present in the object: %s",
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

## 9) Cluster plot (reference view) ------------------------------------------
p_clusters <- DimPlot(muscle_origin, reduction = "umap", label = TRUE) +
  ggtitle("Clusters")

## 10) Two-gene blended visualization for target populations -------------------------------
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


## 11) Display summary panels with patchwork ---------------------------------------
library(patchwork) 

(p_Macro | p_MHCII)  / (p_dc | ImStromal) / (p_schwann | p_neutro) / (p_Plt | p_RBC)
(p_St | p_myocyte) / (p_Pericy | p_Tendo) / (p_EC | p_Stromal_1) / (p_Stromal_2 | p_Stromal_3)

p_schwann # 21
p_mast # 9partial signal
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

Assays(muscle_sc)            # example: "RNA" "SCT" "integrated"
Assays(muscle_dotplot)            # example: "RNA" "SCT" "integrated"

DefaultAssay(muscle_sc)      # current value
DefaultAssay(muscle_dotplot)      # current value

Layers(muscle_sc[["SCT"]])   # layers present in the SCT assay
Layers(muscle_dotplot[["SCT"]])   # layers present in the SCT assay

Layers(muscle_sc[["RNA"]])   # RNA assay (if present)
Layers(muscle_dotplot[["RNA"]])   # RNA assay (if present)

# Match to the appropriate assay
DefaultAssay(muscle_sc) <- "SCT"
DefaultAssay(muscle_dotplot) <- "SCT"

## 12) Labeling: assign cell-type labels--------------------------------------------------
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

Idents(muscle_sc) <- "seurat_clusters"   # for safety
Idents(muscle_dotplot) <- "seurat_clusters"   # for safety

muscle_sc <- RenameIdents(muscle_sc, cluster_cell)
muscle_dotplot <- RenameIdents(muscle_dotplot, cluster_cell)

muscle_sc$celltype <- Idents(muscle_sc)  # store fixed celltype labels in metadata
muscle_dotplot$celltype <- Idents(muscle_dotplot)  # store fixed celltype labels in metadata

DefaultAssay(muscle_sc) <- "RNA"   # Use raw-count-based display
DefaultAssay(muscle_dotplot) <- "RNA"   # Use raw-count-based display

DimPlot(muscle_dotplot, reduction = "umap", label = TRUE) +
  ggtitle("UMAP with assigned cell types")


# Map 10X library names to genotype
geno_map <- c(
  wt1  = "WT",   wt2  = "WT",
  mdx1 = "mdx",  mdx2 = "mdx",
  d2_1 = "mdxD2", d2_2 = "mdxD2"
)

# Confirm that orig.ident is present
table(muscle_sc$orig.ident)
table(muscle_dotplot$orig.ident)

# Add to metadata (overwrite allowed)
muscle_sc$genotype <- base::unname(geno_map[muscle_sc$orig.ident])
muscle_dotplot$genotype <- base::unname(geno_map[muscle_dotplot$orig.ident])

# Convert to a factor and fix the ordering (optional)
muscle_sc$genotype <- factor(muscle_sc$genotype,
                             levels = c("WT","mdx","mdxD2"))
muscle_dotplot$genotype <- factor(muscle_dotplot$genotype,
                                  levels = c("WT","mdx","mdxD2"))

Idents(muscle_sc) <- "celltype"   # for safety
Idents(muscle_dotplot) <- "celltype"   # for safety

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
  dot.scale = 5,            # maximum dot size
  dot.min   = 0.005          # hide values below 1%
) +
  RotatedAxis() +
  theme(                                   # ← adjust font sizes here
    axis.text.x  = element_text(size = 8), # gene labels
    axis.text.y  = element_text(size = 8), # row labels
    strip.text.x = element_text(size = 10)  # "MHC‑II core" column headers such as "MHC-II core"
  ); print(dot_obj)


## Save
#rds_dir <- "C:/DMD_project/rds"
#dir.create(rds_dir, showWarnings = FALSE)
#saveRDS(muscle_sc, file.path(rds_dir, "muscle_sc_after_dotplot.rds"), compress = "gzip")

## Road
#muscle_sc <- readRDS("C:/DMD_project/rds/muscle_sc_after_dotplot.rds")





## ==== DMD scRNA-seq: layer-wise doublet detection -> RBC + doublet exclusion -> pseudobulk -> GLM / covariate-adjusted MHC-II testing ====
suppressPackageStartupMessages({
  library(Seurat)            # Seurat v5 assumed
  library(SeuratObject)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(edgeR)
  library(limma)
})

set.seed(123)
options(stringsAsFactors = FALSE)

## 0) Input ---------------------------------------------------------------
##  Prepare a Seurat v5 object in advance (RNA assay with counts.* layers)
muscle_sc <- readRDS("C:/DMD_project/rds/muscle_sc_after_dotplot.rds")
stopifnot(inherits(muscle_sc, "Seurat"))
stopifnot("RNA" %in% names(muscle_sc@assays))
DefaultAssay(muscle_sc) <- "RNA"

## 1) Enumerate counts.* layers (do not use JoinLayers)--------------------------
layers_all <- Layers(muscle_sc[["RNA"]])
layers <- grep("^counts(\\.|$)", layers_all, value = TRUE)
layers <- setdiff(layers, "counts")   # Exclude the base "counts" layer (restrict to sample-specific layers)
stopifnot(length(layers) > 0)

## 2) Output vectors (initialize with global cell names)
all_cells <- colnames(muscle_sc)
dbl_class <- setNames(rep(NA_character_, length(all_cells)), all_cells)
dbl_score <- setNames(rep(NA_real_,      length(all_cells)), all_cells)

## 2) Run scDblFinder for each layer and align strictly using Cells(..., layer=) ----------
summary_C <- data.frame(layer=character(), cells=integer(), dbl=integer(), ratio=numeric())
for (lay in layers) {
  mat <- muscle_sc[["RNA"]]@layers[[lay]]
  if (is.null(mat) || ncol(mat) == 0) {
    message(sprintf("[skip] %s: empty", lay)); next
  }
  cells_lay <- Cells(muscle_sc, layer = lay)  # Global cell names (matching the matrix column order)
  stopifnot(length(cells_lay) == ncol(mat))
  
  sce <- SingleCellExperiment(list(counts = mat))
  set.seed(123)
  sce <- scDblFinder(sce)  # By default, estimation uses within-sample local density, library size, and related features
  
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

## 3) Write back to metadata (conservatively fill unassigned cells as singlet)-----------
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
#      doublet singlet   Sum
# d2_1     176    6559  6735
# d2_2     216    4397  4613
# mdx1     164    4705  4869
# mdx2      70    1594  1664
# wt1       53    3014  3067
# wt2       49    1681  1730
# Sum      728   21950 22678

## 4) Exclude RBCs and doublets simultaneously -----------------------------------------
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

## 5) Pseudobulk (celltype × orig.ident)---------------------------------
##   Aggregate counts (RNA/slot=counts). Output is genes in rows and celltype_sample in columns
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

## 6) Target gene panel ---------------------------------------------------
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

## 7) edgeR GLM for each cell type ---------------------------------------------
esc <- function(s) gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", s))  # Regular-expression escape

edgeR_by_celltype <- function(ct) {
  cols <- colnames(pseudo_sc)[startsWith(colnames(pseudo_sc), paste0(ct, "_"))]
  if (!length(cols)) return(NULL)
  
  labs <- sub(".*_", "", cols)
  grp  <- factor(sub("^(wt|mdx|d2).*","\\1", labs), levels=c("wt","mdx","d2"))
  if (length(unique(grp)) < 2) return(NULL)  # A single-group design is not valid
  
  y <- DGEList(counts = as.matrix(pseudo_sc[, cols, drop=FALSE]))
  y <- calcNormFactors(y, method="TMMwsp")
  
  X <- model.matrix(~0 + grp)   # dummy variables for wt, mdx, and d2
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

## 8) Hit counts (all genes / target genes)-----------------------------------
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

cat("\n[Hit count: All genes]\n"); print(hit_count(sig_full),   row.names=FALSE)
# [Hit count: All genes]
#           celltype     mdx_vs_wt_n  d2_vs_mdx_n   d2_vs_wt_n
# MHC-low macrophage         1040         588       2272
# Stromal                     664        1390       3068
# MHC-high macrophage         417         314        763
# Endothelial cell(EC)        184         406       1077
# Immune-Stromal(Cd53+)         2          15         93
# Myocyte/Myonucleus           93          75        239
# Tenocyte/Tendon fibroblast   85         215        483
# Neutrophil                   13           4         14
# Dendritic cell(DC)           70         145        300
# Satellite cell               87          76        325
# Pericyte/Vascular SMC        20          43        113
# Schwann cell                 43          50        138

cat("\n[Hit count: Target genes only]\n"); print(hit_count(sig_target), row.names=FALSE)
# [Hit count: Target genes only]
#       celltype        mdx_vs_wt_n  d2_vs_mdx_n  d2_vs_wt_n
# MHC-low macrophage          15          17         22
# Stromal                     14          25         35
# MHC-high macrophage         10           9         12
# Endothelial cell(EC)         7          12         21
# Immune-Stromal(Cd53+)        0           0          4
# Myocyte/Myonucleus           8           3         12
# Tenocyte/Tendon fibroblast   2          11         14
# Neutrophil                   0           0          0
# Dendritic cell(DC)           3           2          8
# Satellite cell               6           2          9
# Pericyte/Vascular SMC        0           1          3
# Schwann cell                 2           1          3

# Execution examples
sig_target[["Endothelial cell(EC)"]][["d2_vs_mdx"]]
sig_target[["Endothelial cell(EC)"]][["mdx_vs_wt"]]
sig_target[["Endothelial cell(EC)"]][["d2_vs_wt"]]



## ---- 0) Prerequisites -------------------------------------------------------------
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

stopifnot(exists("pseudo_sc"))  # Aggregated RNA-count matrix from AggregateExpression (genes x samples)
esc <- function(s) gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", s))

## ---- 1) Representative gene sets (replace as needed) --------------------------
gene_sets <- list(
  ECM_org      = c("Fn1","Fbn2","Fbln1","Col1a1","Col3a1","Col4a1","Col6a1","Col8a1",
                   "Col14a1","Col16a1","Plod2","Postn","Vcan","Adam12"),
  Platelet_act = c("Cd63","Plek","Pf4","Selp","Fcer1g","Itgb1","Vwf","Col1a1","Fn1"),
  MHCII        = c("H2-Aa","H2-Ab1","H2-Eb1","Ciita"),
  
  # --- Split lipid-related programs into subcategories ---
  Lipo_Synthesis      = c("Acaca","Fasn","Scd1","Elovl6","Me1","Acly","Acss2",
                          "Gpam","Agpat2","Lpin1","Dgat1","Dgat2","Srebf1","Pparg","Cebpa","Cebpb"),
  FA_Uptake_Transport = c("Cd36","Lpl","Gpihbp1","Fabp3","Fabp4","Fabp5","Ldlr","Scarb1","Mfsd2a"),
  FAO_mito            = c("Cpt1a","Cpt1b","Acadvl","Acadm","Hadha","Hadhb","Echs1","Etfdh","Etfb","Acaa2"),
  FAO_perox           = c("Acox1","Ehhadh","Hsd17b4","Abcd3","Acot1","Pecr","Pex11a"),
  
  # Existing metabolic programs
  OxPhos          = c("Ndufa1","Ndufa2","Ndufb5","Cox4i1","Atp5f1a","Atp5f1b","Atp5mc1","Uqcrc1","Uqcrq",
                      "Ndufs2","Ndufb8","Cox5a","Cox6a1","Atp5me","Atp5md","Uqcrc2","Uqcrfs1"),
  TCA_cycle       = c("Cs","Aco2","Idh3a","Idh3b","Ogdh","Dlst","Sucb1","Sdha","Sdhb","Fh1","Mdh2","Pck2","Mpc1","Mpc2","Pdhb","Dlat"),
  Mito_Biogenesis = c("Ppargc1a","Ppargc1b","Nrf1","Nrf2","Tfam","Tfb1m","Tfb2m","Tomm20","Timm23","Polg")
)

## ---- 2) Shared utilities -------------------------------------------------
.get_cols_by_ct <- function(ct){
  pfx  <- paste0(ct, "_")
  cols <- colnames(pseudo_sc)[startsWith(colnames(pseudo_sc), pfx)]
  if (!length(cols)) stop(sprintf("No columns found: celltype='%s'", ct))
  cols
}


.get_grp_from_labels <- function(labs){
  # Example labels: allow formats such as "wt1", "mdx_2", and "d2-1"
  grp <- sub("^(wt|mdx|d2).*", "\\1", labs)
  factor(grp, levels = c("wt","mdx","d2"))
}

.build_design_nocov <- function(y, labs){
  grp <- .get_grp_from_labels(labs)
  if (all(table(grp) == 0)) stop("Group labels could not be recognized (column names must end with wt/mdx/d2)")
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

## ---- 3) Run camera/roast in batch ----------------------------------
camera_all_sets <- function(ct,
                            contrast = c("d2_vs_mdx","mdx_vs_wt","d2_vs_wt"),
                            gene_sets){
  contrast <- match.arg(contrast)
  
  ## -- Extract pseudobulk counts for the cell type
  cols <- .get_cols_by_ct(ct)
  pb   <- pseudo_sc[, cols, drop = FALSE]
  
  ## -- Labels and DGE object
  labs <- sub(".*_", "", cols)   # terminal tokens (e.g., wt1 / mdx2 / d2_1)
  y    <- DGEList(counts = as.matrix(pb))
  y    <- calcNormFactors(y, method = "TMMwsp")
  
  ## -- Design
  X <- .build_design_nocov(y, labs)
  # Note: voom requires a design matrix, and camera takes the voom object as input.
  v <- voom(y, X, plot = FALSE)
  
  ## -- Contrasts
  contr <- .make_contrasts3(X)
  K <- switch(contrast,
              mdx_vs_wt = contr[,"mdx_vs_wt"],
              d2_vs_mdx = contr[,"d2_vs_mdx"],
              d2_vs_wt  = contr[,"d2_vs_wt"])
  
  ## -- Gene-set index (only sets with >=3 genes)
  idx_list <- .index_sets(gene_sets, rownames(v$E))
  if (!length(idx_list)) stop("No valid gene sets are available (>=3 genes required).")
  
  ## -- camera / roast (with reproducible random rotations)
  set.seed(123)
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

## ---- 4) Helper for summarized gene-set output --------------------------------------
## Display only significant sets (FDR < alpha)
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
  
  cat(sprintf("\n[%s | %s] camera-significant sets(FDR < %.3f): %dhits\n",
              res$celltype, res$contrast, alpha, nrow(sig)))
  
  if (nrow(sig) == 0L) {
    cat("  No significant gene sets were detected.\n")
  } else {
    print(sig)
  }
  
  cat("\nSample counts (wt/mdx/d2):\n")
  print(res$n_samples)
  invisible(res)
}


## ---- 5) Usage examples (run as needed) ---------------------------------------
## d2_vs_mdx test
res_ec <- camera_all_sets("Endothelial cell(EC)", "d2_vs_mdx", gene_sets)
print_top_sets(res_ec, alpha = 0.05)  # default = 0.05
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


## ---------- ①-0 Preparation ----------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr)
})
dir.create("C:/DMD_project/figs", showWarnings = FALSE, recursive = TRUE)

stopifnot(exists("edgeR_out"), exists("pseudo_sc"))
mhc_core <- c("H2-Aa","H2-Ab1","H2-Eb1","Ciita")

## ---------- ①-A Representative gene panel: dot plot and heatmap ----------
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
  c("Fn1","Lum","Thbs2","Col1a1","Col15a1","Postn","Vcan"),  # ECM/fibrosis
  c("Lpl","Spp1","C1qa","C1qb","C1qc","Cd53","Ctss"),        # FIM/immune
  c("Tnni2","Des","Pax7","Myf5"),                            # myogenic/satellite
  c("Pecam1","Cdh5","Pdgfra","Pi16","Rgs5","Abcc9","Scx","Tnmd") # EC/stromal/pericyte/tendon
))

df_panel <- ed_list |>
  dplyr::filter(contrast == "d2_vs_mdx", gene %in% panel_genes)

p_panel <- ggplot(df_panel, aes(x = gene, y = celltype)) +
  geom_point(aes(size = pmin(-log10(FDR), 10), color = logFC)) +
  scale_color_gradient2(low = "steelblue", mid = "grey95", high = "firebrick", midpoint = 0) +
  scale_size(range = c(0.5, 5)) +
  labs(title = "Pseudo Bulk(d2_vs_mdx):The effect size and statistical significance of representative genes",
       x = NULL, y = NULL, color = "logFC", size = "-log10(FDR)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))


## ---------- ①-B camera results: gene-set by cell-type tiles ----------
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
  filter(FDR < 0.05)  # Color only significant results

p_cam <- ggplot(df_cam, aes(x = geneset, y = celltype, fill = signed)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0,
                       name = "Direction ×\n -log10(FDR)") +
  labs(title = "Pseudo Bulk(d2_vs_mdx):Statistical significance of gene sets(camera)", x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


## ---------- ①-C Volcano plot of the MHC-II core (faceted by cell type) ----------
lab_genes <- mhc_core
df_mhc <- ed_list |>
  dplyr::filter(contrast == "d2_vs_mdx")

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
base_vol <- ggplot(df_mhc, aes(x = logFC, y = -log10(FDR))) +
  geom_point(alpha = 0.30, size = 0.8) +
  geom_point(data = subset(df_mhc, gene %in% lab_genes),
             color = "firebrick", size = 1.5) +
  labs(title = "Pseudo Bulk(d2_vs_mdx):Volcano Plot: MHC-II Core",
       x = "log2 fold-change", y = "-log10(FDR)") +
  facet_wrap(~ celltype, scales = "free_y") +
  theme_bw(base_size = 9)

p_vol <- if (has_ggrepel) {
  base_vol + ggrepel::geom_text_repel(
    data = subset(df_mhc, gene %in% lab_genes),
    aes(label = gene), size = 3, max.overlaps = Inf)
} else base_vol; p_vol



## ==== Pseudobulk camera + fgsea pipeline (Seurat v5) =======================
install.packages(c("fgsea","msigdbr","data.table"), dependencies = TRUE)
# Use Bioconductor only if fgsea cannot be installed from CRAN
if (!requireNamespace("fgsea", quietly = TRUE)) BiocManager::install("fgsea")

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(SingleCellExperiment); library(scDblFinder)
  library(edgeR); library(limma)
  library(fgsea)         # GSEA
  suppressWarnings(library(msigdbr))  # Not required in offline settings
  library(data.table)
})

set.seed(123); options(stringsAsFactors = FALSE)

## 0) Input (Seurat v5 object; RNA assay with counts.* layers required)
obj_path <- "C:/DMD_project/rds/muscle_sc_after_dotplot.rds"
muscle_sc <- readRDS(obj_path)
stopifnot(inherits(muscle_sc, "Seurat"))
DefaultAssay(muscle_sc) <- "RNA"

## 1) Run scDblFinder for each layer
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

## 2) RBC + doublet exclusion
stopifnot("celltype" %in% colnames(muscle_sc@meta.data))
muscle_sc <- subset(muscle_sc,
                    subset = (celltype != "Erythrocyte") &
                      (scDblFinder.class != "doublet"))
muscle_sc$celltype <- droplevels(muscle_sc$celltype)
Idents(muscle_sc) <- "celltype"

## 3) Pseudobulk (celltype × orig.ident)
agg <- AggregateExpression(
  muscle_sc, group.by = c("celltype","orig.ident"),
  assays = "RNA", slot = "counts"
)
pb_counts <- agg$RNA       # genes x (celltype_sample)  matrix
stopifnot(is.matrix(pb_counts) || inherits(pb_counts, "dgCMatrix"))

## 4) Utilities
.parse_group <- function(labels) {
  # Extract group labels from the suffix wt / mdx / d2 (e.g., "Myocyte_wt1" -> "wt")
  grp <- sub("^(wt|mdx|d2).*", "\\1", labels)
  factor(grp, levels = c("wt","mdx","d2"))
}
.get_cols_by_ct <- function(ct) {
  cols <- colnames(pb_counts)[startsWith(colnames(pb_counts), paste0(ct, "_"))]
  if (!length(cols)) stop(sprintf("No columns found: %s", ct))
  cols
}

## ==== Pathway utility: absorb msigdbr / GMT variability and return a named list ====
suppressPackageStartupMessages({
  library(msigdbr)    # Use in online settings (with fallback if unavailable)
  library(fgsea)      # Use gmtPathways
  library(edgeR); library(limma); library(data.table)
})

get_pathways_safe <- function(collections = c("REACTOME","GO:BP","CP:KEGG"),
                              species = "Mus musculus",
                              min_size = 10, max_size = 500,
                              gmt_files = list()) {
  ## collections include "REACTOME", "GO:BP", "HALLMARK", "CP:KEGG" and related collection names
  collections <- toupper(collections)
  
  build_list_from_tibble <- function(tb) {
    stopifnot(all(c("gs_name","gene_symbol") %in% names(tb)))
    split(tb$gene_symbol, tb$gs_name)
  }
  
  out_list <- list()
  
  for (coll in collections) {
    ## 1) First try msigdbr
    msig_ok <- FALSE; lst <- NULL
    cat(sprintf("[get_pathways] try msigdbr for %s ...\n", coll))
    tb <- try({
      ## Map coll to C2/C5
      if (coll == "REACTOME") {
        # db_species is used in newer msigdbr versions; older versions ignore it, so a two-step tryCatch is used
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
        # Some sets may be absent from MSigDB for licensing reasons; fall back to GMT if needed
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
    
    ## 2) Fall back to GMT if retrieval fails or returns empty
    if (!msig_ok) {
      gmt_path <- gmt_files[[coll]]
      if (is.null(gmt_path)) {
        cat(sprintf("[get_pathways] msigdbr failed and no GMT for %s. Skip.\n", coll))
        next
      }
      cat(sprintf("[get_pathways] use GMT for %s : %s\n", coll, gmt_path))
      lst <- fgsea::gmtPathways(gmt_path)  # This already returns a named list
    }
    
    ## 3) Filter by size and remove duplicates
    lst <- lst[!duplicated(names(lst))]
    keep <- vapply(lst, length, integer(1))
    lst <- lst[keep >= min_size & keep <= max_size]
    cat(sprintf("[get_pathways] %s -> %d gene sets kept (size %d–%d)\n",
                coll, length(lst), min_size, max_size))
    out_list <- c(out_list, lst)
  }
  
  out_list
}

## ==== Pseudobulk: celltype-wise limma-voom + camera & fgsea ====
## Prerequisite: pb_counts is a genes x samples count matrix,
##       column names must follow the "Celltype_sampleLabel" format (e.g., "Myocyte/Myonucleus_mdx1")
## ==== Prerequisite: pb_counts is a genes x samples pseudobulk matrix ====
## Required packages
suppressPackageStartupMessages({
  library(edgeR); library(limma)
  library(fgsea); library(data.table)
  library(BiocParallel)
})

## Disable parallelization to avoid environment-dependent warnings and errors
BiocParallel::register(BiocParallel::SerialParam())
data.table::setDTthreads(1)

## Restate prerequisites
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
  
  ## ---- Use separate formats for camera and fgsea ----
  ## For fgsea: gene-symbol list with size constraints applied
  pw_list <- lapply(pathways, function(gs) intersect(gs, genes_in))
  sz <- vapply(pw_list, length, 1L)
  pw_list <- pw_list[sz >= min_gs & sz <= max_gs]
  
  ## For camera: convert the same sets to row indices
  idx <- lapply(pw_list, function(gs) match(gs, genes_in))
  
  ## ---- camera (with inter-gene correlation adjustment) ----
  cam <- camera(v, index = idx, design = X, contrast = K)
  cam$FDR <- p.adjust(cam$PValue, "BH")
  cam <- cam[order(cam$FDR, cam$PValue), , drop = FALSE]
  
  ## ---- fgsea (Multilevel with fallback to Simple) ----
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


## ==== Execution examples(Myocyte/Satellite, d2_vs_mdx) ====
## 1) Prepare pathway collections: use msigdbr online or specify GMT files offline
pathways_all <- get_pathways_safe(
  collections = c("REACTOME","GO:BP"),           # Add "CP:KEGG" or "HALLMARK" if needed
  species     = "Mus musculus",
  min_size    = 10, max_size = 500,
  gmt_files   = list(
    # If local GMT files are available, specify the paths here (otherwise leave as NULL)
    # "REACTOME" = "data/c2.cp.reactome.mouse.symbols.gmt",
    # "GO:BP"    = "data/c5.go.bp.mouse.symbols.gmt",
    # "CP:KEGG"  = "data/c2.cp.kegg.symbols.gmt"
    # Left commented out as examples
  )
)

## Pathway collection (pass the output of get_pathways_safe directly)
## pathways_all is a named list: pathway name -> gene-symbol vector
stopifnot(exists("pathways_all"))

## Start with Myocyte / d2_vs_mdx
res_myocyte <- run_ct_tests_safe("Myocyte/Myonucleus", "d2_vs_mdx", pathways_all)

## Filter focus names after matching them to the sets that were actually retained
focus_names <- grep(paste(c(
  "OXIDATIVE_PHOSPHORYLATION", "RESPIRATORY_ELECTRON_TRANSPORT",
  "TCA", "FATTY_ACID.*OXIDATION", "MITOCHONDRIAL",
  "EXTRACELLULAR_MATRIX_ORGANIZATION", "COLLAGEN", "MATRIX_METALLOPROTEINASE",
  "PLATELET.*DEGRANULATION", "TGF.*BETA", "ADIPOGENESIS", "PPAR"
), collapse="|"), names(pathways_all), value = TRUE)

focus_kept <- intersect(focus_names, res_myocyte$sets_kept)

## Consistency check between camera and fgsea
res_myocyte$camera[rownames(res_myocyte$camera) %in% focus_kept, ][1:22, ]
#                                                                        NGenes  Direction         FDR
#1 GOBP_OXIDATIVE_PHOSPHORYLATION                                           114     Down  0.00000000000003156473
#2 REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT          182     Down  0.00000000000108061961
#3 REACTOME_RESPIRATORY_ELECTRON_TRANSPORT                                  120     Down  0.00000000000605963315
#4 GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE                  47     Down  0.00000000006654059154
#5 REACTOME_COLLAGEN_FORMATION                                               15       Up  0.00000004591237902694
#6 REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES     12       Up  0.00000014562924397284
#7 REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION                                46       Up  0.00000056750552993768
#8 REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES                      11       Up  0.00000378822332895240
#9 GOBP_COLLAGEN_FIBRIL_ORGANIZATION                                         10       Up  0.00000578729254059024
#10 REACTOME_COLLAGEN_DEGRADATION                                            11       Up  0.00001542674307388481
#11 GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_UBIQUINOL_TO_CYTOCHROME_C          12     Down  0.00002002958620514243
#12 GOBP_COLLAGEN_METABOLIC_PROCESS                                          22       Up  0.00007463809463265978
#13 GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN             18     Down  0.00021503451153735673
#14 REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION                               65     Down  0.00021949123698015828
#15 REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE                                     28     Down  0.00026208064434801368
#16 REACTOME_MITOCHONDRIAL_BIOGENESIS                                        45     Down  0.00026629149942877599
#17 GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY                    67     Down  0.00210054022750523527
#18 GOBP_INNER_MITOCHONDRIAL_MEMBRANE_ORGANIZATION                           24     Down  0.01338857659084882572
#19 REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT                                    50     Down  0.01457478183833985086
#20 REACTOME_MATURATION_OF_TCA_ENZYMES_AND_REGULATION_OF_TCA_CYCLE           17     Down  0.02600637811184790260
#21 GOBP_MITOCHONDRIAL_GENE_EXPRESSION                                       86     Down  0.04252474814994626484
#22 REACTOME_MITOCHONDRIAL_TRANSLATION                                       71     Down  0.04908644421550474590


res_myocyte$fgsea [pathway %in% focus_kept][order(padj)][1:23]
#                                                                   pathway            FDR                           NES                           
# 1:       REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT   0.00000000000000000000000000000178101 -3.280589   
# 2:                                        GOBP_OXIDATIVE_PHOSPHORYLATION   0.00000000000000000000000000015473566 -3.445820   
# 3:                               REACTOME_RESPIRATORY_ELECTRON_TRANSPORT   0.00000000000000000000000039617273774 -3.325534   
# 4:              GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE   0.00000000000000125674580941039857999 -3.163266   
# 5:                            REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION   0.00000148437541074624543063417703159  2.293537 
# 6:                            REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION   0.00000492683984234174676512807433681 -2.303110  
# 7:                 GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY   0.00002839651545089853749460895304679 -2.205457 
# 8:                                     REACTOME_MITOCHONDRIAL_BIOGENESIS   0.00004695740292645269865941248799146 -2.334163 
# 9:       GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_UBIQUINOL_TO_CYTOCHROME_C   0.00006410354196609753540687243189566 -2.326744 
# 10:                                       GOBP_COLLAGEN_METABOLIC_PROCESS  0.00008455811161698221687342547081556  2.172928  
# 11:          GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN  0.00009648295402232930485439399426184 -2.312916  
# 12: REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES  0.00028301020418923484526368961056164  2.052750  
# 13:                                  REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE  0.00031412920434968681724502426177992 -2.320893  
# 14:                                           REACTOME_COLLAGEN_FORMATION  0.00032219237751760540066497950917324  2.051906  
# 15:                                         REACTOME_COLLAGEN_DEGRADATION  0.00055955631797652987963892501355190  2.042714  
# 16:                        GOBP_INNER_MITOCHONDRIAL_MEMBRANE_ORGANIZATION  0.00347929092930622475485225031377468 -2.120391  
# 17:                  REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES  0.00440822069480227562987550982143148  1.909381  
# 18:                                     GOBP_COLLAGEN_FIBRIL_ORGANIZATION  0.00737205809345011561090066365409257  1.835492  
# 19:                                 REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT  0.00818580698123624710127455017527609 -1.820759  
# 20:                                    REACTOME_MITOCHONDRIAL_TRANSLATION  0.02054847378703995863791220699567930 -1.602833  
# 21:                                    GOBP_MITOCHONDRIAL_GENE_EXPRESSION  0.02071480022600939024934874055361433 -1.569344  
# 22:                                        GOBP_MITOCHONDRIAL_TRANSLATION  0.03913632692061281420636831285264634 -1.517916  
# 23:        REACTOME_MATURATION_OF_TCA_ENZYMES_AND_REGULATION_OF_TCA_CYCLE  0.04967459051662621860590007827340742 -1.677529  



## ==== String formatting and visualization: conflict-safe version ==========================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(stringr); library(data.table)
})

format_pathway_labels_unique <- function(pathways, wrap_width = 36){
  raw <- pathways
  src <- ifelse(startsWith(raw, "REACTOME_"), " [Reactome]",
                ifelse(startsWith(raw, "GOBP_"),     " [GO]", ""))
  
  lab <- gsub("^REACTOME_|^GOBP_", "", raw)
  lab <- gsub("_", " ", lab)
  lab <- stringr::str_to_sentence(lab)
  lab <- gsub("Tca","TCA", lab, fixed=TRUE)
  lab <- gsub("Ppar","PPAR",lab, fixed=TRUE)
  lab <- gsub("Ecm","ECM",  lab, fixed=TRUE)
  lab <- gsub("Mmp","MMP",  lab, fixed=TRUE)
  
  lab <- paste0(lab, src)
  lab <- make.unique(lab, sep = " ")        # ★ Resolve duplicate labels completely here
  lab <- stringr::str_wrap(lab, width = wrap_width)
  stats::setNames(lab, raw)
}

## 2) Reconstruct the running ES (using the same p=1 setting as fgsea)
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

## 3) Construct ridge-plot data (with downsampling)
build_ridge_df <- function(stats, gene_sets, ids, sample_each = 1500) {
  ids <- intersect(ids, names(gene_sets))
  out <- vector("list", length(ids)); k <- 0L
  for (id in ids) {
    genes <- intersect(gene_sets[[id]], names(stats))
    es <- calc_running_es(stats, genes, gsea_p = 1)
    if (length(es) > sample_each) {
      take <- unique(round(seq(1, length(es), length.out = sample_each)))
      es <- es[take]
    }
    k <- k + 1L
    out[[k]] <- data.frame(pathway = id, es = es, stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(out)
}

## 4) Visualization utility (assuming ggplot2 >= 3.4: use linewidth for line width)
plot_ridge_constantfill <- function(df, pretty_map, ids_to_show,
                                    fill_col, title_text, legend_title) {
  df$label <- factor(pretty_map[df$pathway],
                     levels = unique(pretty_map[ids_to_show]))  # ★ avoid duplicated factor levels
  ggplot2::ggplot(df, ggplot2::aes(x = es, y = label, group = label,
                                   fill = !!as.name(fill_col))) +
    ggridges::geom_density_ridges(scale = 1.2, rel_min_height = 0.01,
                                  linewidth = 0.2, color = "white") +  # ★ linewidth
    ggplot2::scale_x_continuous(name = "Running enrichment score",
                                limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    ggplot2::scale_fill_gradient(low = "#b2182b", high = "#2166ac",
                                 name = legend_title) +
    ggplot2::labs(y = NULL, title = title_text) +
    ggridges::theme_ridges(font_size = 11) +                           # ★ use explicit namespace qualification
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "right")
}


## 5) Example plot using the top 23 pathways (21 for camera)
ids_to_show <- as.character(
  res_myocyte$fgsea[pathway %in% focus_kept][order(padj)][1:23, ]$pathway
)

ranks_named <- setNames(res_myocyte$topTable$t, rownames(res_myocyte$topTable))
ranks_named <- sort(ranks_named[!is.na(ranks_named)], decreasing = TRUE)

ridge_df <- build_ridge_df(ranks_named, pathways_all, ids_to_show, sample_each = 1500)

pretty_map <- format_pathway_labels_unique(unique(ridge_df$pathway), wrap_width = 36)

camera_fdr_map <- setNames(res_myocyte$camera$FDR, rownames(res_myocyte$camera))
fgsea_padj_map <- setNames(res_myocyte$fgsea$padj,  res_myocyte$fgsea$pathway)
ridge_df$camera_fdr <- camera_fdr_map[ridge_df$pathway]
ridge_df$fgsea_padj <- fgsea_padj_map[ridge_df$pathway]

plot_ridge_camera <- plot_ridge_constantfill(
  ridge_df, pretty_map, ids_to_show,
  fill_col = "camera_fdr",
  title_text = "Ridge plot(camera/color = camera FDR；curve = rank)",
  legend_title = "camera FDR"
)

plot_ridge_fgsea <- plot_ridge_constantfill(
  ridge_df, pretty_map, ids_to_show,
  fill_col = "fgsea_padj",
  title_text = "Ridge plot(fgsea/color = fgsea padj；curve = rank)",
  legend_title = "fgsea padj"
)

# print(plot_ridge_camera)
print(plot_ridge_fgsea)

## 8) [Optional] Mathematical consistency check: verify that the fgsea ES matches the reconstructed ES
check_es <- sapply(ids_to_show, function(id){
  es <- calc_running_es(ranks_named, intersect(pathways_all[[id]], names(ranks_named)))
  es_star <- if (abs(max(es)) >= abs(min(es))) max(es) else min(es)  # same definition as the fgsea ES
  c(reconstructed_ES = es_star,
    fgsea_ES = res_myocyte$fgsea[pathway == id, ES][1])
})
print(t(check_es))  # The two columns (reconstructed_ES and fgsea_ES) should be nearly identical
# reconstructed_ES   fgsea_ES
# REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT             -0.6506216 -0.6506216
# GOBP_OXIDATIVE_PHOSPHORYLATION                                              -0.7364245 -0.7364245
# REACTOME_RESPIRATORY_ELECTRON_TRANSPORT                                     -0.7000533 -0.7000533
# GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE                    -0.7968745 -0.7968745
# REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION                                   0.6971782  0.6971782
# REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION                                  -0.5481986 -0.5481986
# GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY                       -0.5209188 -0.5209188
# GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_UBIQUINOL_TO_CYTOCHROME_C             -0.8532404 -0.8532404
# GOBP_COLLAGEN_METABOLIC_PROCESS                                              0.7692346  0.7692346
# REACTOME_MITOCHONDRIAL_BIOGENESIS                                           -0.5874281 -0.5874281
# GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN                -0.7628066 -0.7628066
# REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES        0.8480716  0.8480716
# REACTOME_COLLAGEN_FORMATION                                                  0.8009356  0.8009356
# REACTOME_COLLAGEN_DEGRADATION                                                0.8588034  0.8588034
# REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE                                        -0.6483891 -0.6483891
# GOBP_INNER_MITOCHONDRIAL_MEMBRANE_ORGANIZATION                              -0.6223046 -0.6223046
# REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES                         0.8027470  0.8027470
# REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT                                       -0.4563661 -0.4563661
# GOBP_COLLAGEN_FIBRIL_ORGANIZATION                                            0.7924162  0.7924162
# GOBP_MITOCHONDRIAL_GENE_EXPRESSION                                          -0.3542381 -0.3542381
# REACTOME_MITOCHONDRIAL_TRANSLATION                                          -0.3735884 -0.3735884
# GOBP_MITOCHONDRIAL_TRANSLATION                                              -0.3435387 -0.3435387
# REACTOME_MATURATION_OF_TCA_ENZYMES_AND_REGULATION_OF_TCA_CYCLE              -0.5600417 -0.5600417




## Myofiber-type (slow/fast) analysis: covariate adjustment plus within-type stratification (revised version)-------------
## Prerequisite: muscle_sc (after RBC/doublet removal and annotation) and pb_counts are available

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(edgeR);  library(limma); library(dplyr); library(stringr)
})

set.seed(123); options(stringsAsFactors = FALSE)

## 0) Prerequisite checks and utilities --------
stopifnot(exists("muscle_sc"), inherits(muscle_sc, "Seurat"))
stopifnot(exists("pb_counts"))

if (!"celltype" %in% colnames(muscle_sc@meta.data))
  stop("muscle_sc does not contain 'celltype' (check the annotation)")

if (!exists(".parse_group")) {
  .parse_group <- function(labels)
    factor(sub("^(wt|mdx|d2).*","\\1", labels), levels = c("wt","mdx","d2"))
}

.get_cols_by_prefix <- function(prefix, mat) {
  ## Example: for "Myocyte/Myonucleus", use ^Myocyte/Myonucleus[_:.]
  pat <- paste0("^", gsub("([][\\^$.|?*+(){}\\\\])","\\\\\\1", prefix), "([_.:])")
  grep(pat, colnames(mat), value = TRUE, perl = TRUE)
}
.keep_present <- function(genes, obj) {
  g <- intersect(genes, rownames(obj))
  if (!length(g)) warning("Specified genes were not found: ", paste(head(genes, 5), "..."))
  g
}
v_add <- function(x, nm){
  v <- setNames(rep(NA_real_, ncol(muscle_sc)), colnames(muscle_sc))
  v[colnames(myo_obj)] <- x
  assign("muscle_sc", AddMetaData(muscle_sc, v, col.name = nm), inherits = TRUE)
}

## 1) Myocyte slow/fast module scores -> fiber_index --------
DefaultAssay(muscle_sc) <- "RNA"
if (!"data" %in% Layers(muscle_sc[["RNA"]])) {
  muscle_sc <- NormalizeData(muscle_sc, assay = "RNA",
                             normalization.method = "LogNormalize", verbose = FALSE)
}
DefaultLayer(muscle_sc[["RNA"]]) <- "data"

cells_myo <- WhichCells(muscle_sc, idents = "Myocyte/Myonucleus")
stopifnot(length(cells_myo) > 0)

slow_genes_raw <- c("Myh7","Myl2","Tnni1","Tnnt1","Mb")
fast_genes_raw <- c("Myh1","Myh2","Myh4","Myh13","Tnni2","Tnnt3","Actn3")
slow_genes <- .keep_present(slow_genes_raw, muscle_sc)
fast_genes <- .keep_present(fast_genes_raw, muscle_sc)

myo_obj <- subset(muscle_sc, cells = cells_myo)
myo_obj <- AddModuleScore(myo_obj, features = list(slow_genes), name = "SLOW")
myo_obj <- AddModuleScore(myo_obj, features = list(fast_genes), name = "FAST")
myo_obj$slow_score  <- myo_obj$SLOW1
myo_obj$fast_score  <- myo_obj$FAST1
myo_obj$fiber_index <- myo_obj$slow_score - myo_obj$fast_score  # slow-dominant (+) / fast-dominant (-)

## Write back to the parent object
v_add(myo_obj$slow_score,  "slow_score")
v_add(myo_obj$fast_score,  "fast_score")
v_add(myo_obj$fiber_index, "fiber_index")

## 2) Construct sample-level covariates and a group-orthogonalized index --------
cov_df <- myo_obj@meta.data %>%
  mutate(sample = gsub("-", "_", orig.ident)) %>%
  group_by(sample) %>%
  summarise(
    slow_mean  = mean(slow_score,  na.rm = TRUE),
    fast_mean  = mean(fast_score,  na.rm = TRUE),
    index_mean = mean(fiber_index, na.rm = TRUE),
    n_nuclei   = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    slow_c  = as.numeric(scale(slow_mean,  scale = FALSE)),
    fast_c  = as.numeric(scale(fast_mean,  scale = FALSE)),
    index_c = as.numeric(scale(index_mean, scale = FALSE))
  )
print(cov_df)
#  sample  slow_mean fast_mean index_mean n_nuclei    slow_c  fast_c index_c
# 1 d2_1      0.0959      1.40     -1.30        37 -0.0545   -0.111   0.0561
# 2 d2_2      0.114       1.02     -0.904       28 -0.0365   -0.489   0.452 
# 3 mdx1      0.168       1.59     -1.42       160  0.0180    0.0791 -0.0611
# 4 mdx2      0.188       1.56     -1.37        89  0.0380    0.0544 -0.0164
# 5 wt1       0.185       1.91     -1.73       168  0.0348    0.406  -0.372 
# 6 wt2       0.151       1.57     -1.42       225  0.000129  0.0597 -0.0596

## 3) Myocyte pseudobulk -> voom/limma design (with group-orthogonalized index) --------
pb_cols_myo <- .get_cols_by_prefix("Myocyte/Myonucleus", pb_counts)
pb_myo <- as.matrix(pb_counts[, pb_cols_myo, drop = FALSE])

## Extract sample labels from column names and standardize the format
labs <- sub("^[^_:]+[_:]", "", colnames(pb_myo))    # "Myocyte/Myonucleus_d2-1" → "d2-1"
labs <- gsub("-", "_", labs)
grp  <- .parse_group(labs)

## Orthogonalize the index to group (absorb group-mean differences)
idx_raw   <- cov_df$index_c[ match(labs, cov_df$sample) ]
idx_raw[is.na(idx_raw)] <- 0
idx_resid <- resid(lm(idx_raw ~ 0 + grp))
idx_resid <- as.numeric(scale(idx_resid, scale = FALSE))

## voom/limma
y <- DGEList(counts = pb_myo)
keep <- filterByExpr(y, group = grp)
y <- y[keep,, keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMMwsp")

design_ortho <- cbind(model.matrix(~0 + grp), index_resid = idx_resid)
stopifnot(nrow(design_ortho) == ncol(y))

v    <- voomWithQualityWeights(y, design_ortho, plot = FALSE)
fit0 <- lmFit(v, design_ortho)
K    <- makeContrasts(
  mdx_vs_wt = grpmdx - grpwt,
  d2_vs_mdx = grpd2  - grpmdx,
  d2_vs_wt  = grpd2  - grpwt, levels = design_ortho
)
fit2 <- eBayes(contrasts.fit(fit0, K))
genes_in <- rownames(fit2)

## 4) Create consensus MSigDB sets (1:1 correspondence with the ridge plot) --------
msig_ids <- list(
  OxPhos_Ridge = c(
    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
    "GO_OXIDATIVE_PHOSPHORYLATION",
    "GO_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE",
    "GO_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN"
  ),
  TCA_Ridge = c(
    "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
    "GO_TRICARBOXYLIC_ACID_CYCLE",
    "REACTOME_MATURATION_OF_TCA_ENZYMES_AND_REGULATION_OF_TCA_CYCLE"
  ),
  Mito_Transl_Import_Ridge = c(
    "REACTOME_MITOCHONDRIAL_TRANSLATION",
    "REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT",
    "GO_MITOCHONDRIAL_GENE_EXPRESSION",
    "GO_INNER_MITOCHONDRIAL_MEMBRANE_ORGANIZATION"
  ),
  ECM_Collagen_Ridge = c(
    "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
    "REACTOME_COLLAGEN_FORMATION",
    "GO_COLLAGEN_FIBRIL_ORGANIZATION",
    "REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES",
    "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES",
    "REACTOME_COLLAGEN_DEGRADATION",
    "GO_COLLAGEN_METABOLIC_PROCESS"
  )
)

## Prefer msigdbr; otherwise use gmtPathways (fgsea required)
gmt_file <- NULL  # example: "resources/msigdb.v2023.2.Mm.symbols.gmt"
.pull_msig <- function(id_vec, species = "Mus musculus", gmt_file = NULL) {
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    tb <- msigdbr::msigdbr(species = species)
    tb <- tb[tb$gs_name %in% id_vec, c("gs_name", "gene_symbol")]
    split(tb$gene_symbol, tb$gs_name)
  } else {
    if (is.null(gmt_file))
      stop("msigdbr is not available; please provide a GMT path in gmt_file.")
    if (!requireNamespace("fgsea", quietly = TRUE))
      stop("GMT import requires fgsea (install fgsea or use msigdbr).")
    lst <- fgsea::gmtPathways(gmt_file)
    keep <- intersect(names(lst), id_vec)
    setNames(lst[keep], keep)
  }
}

Gene <- function(...) unique(unlist(list(...), use.names = FALSE))
gene_sets_core <- lapply(msig_ids, function(ids) {
  raw <- .pull_msig(ids, species = "Mus musculus", gmt_file = gmt_file)
  Gene(unlist(unname(raw), use.names = FALSE))
})

## Control set size (exclude extremely large or small sets)
min_size <- 10; max_size <- 2000
gene_sets_core <- Filter(function(gs) length(gs) >= min_size && length(gs) <= max_size, gene_sets_core)
cat("[gene_sets_core sizes]\n"); print(vapply(gene_sets_core, length, 1L))

## Map to camera indices (row indices)
idx_map <- lapply(gene_sets_core, function(gs){
  ix <- match(intersect(gs, genes_in), genes_in)
  ix[!is.na(ix)]
})
sizes <- vapply(idx_map, length, 1L)
if (any(sizes < 3L)) warning("Some sets have size < 3: ", paste(names(sizes)[sizes<3L], collapse = ", "))

## 5) Run camera / fgsea (using d2_vs_mdx as an example) --------
camera_with_new_sets <- function(contrast = c("d2_vs_mdx","mdx_vs_wt","d2_vs_wt")) {
  contrast <- match.arg(contrast)
  cam <- camera(v, index = idx_map, design = design_ortho, contrast = K[, contrast])
  cam$FDR <- p.adjust(cam$PValue, "BH")
  cam[order(cam$FDR, cam$PValue), , drop = FALSE]
}

cat("\n[Myocyte | orthogonalized-index adjustment | camera] d2_vs_mdx (reported with the updated consensus labels)\n")
print(camera_with_new_sets("d2_vs_mdx"))
#                          NGenes Direction           FDR
# OxPhos_Ridge                120      Down  8.833408e-14
# ECM_Collagen_Ridge           46        Up  1.201938e-09
# TCA_Ridge                    28      Down  2.086541e-05
# Mito_Transl_Import_Ridge    121      Down  6.051588e-04

# 4) fgsea(only when fgsea is available; pathway names are shown with the new labels)
if (requireNamespace("fgsea", quietly = TRUE)) {
  suppressPackageStartupMessages(library(data.table))
  ranks <- setNames(fit2$t[, "d2_vs_mdx"], rownames(fit2$t))
  pw_for_fgsea <- lapply(gene_sets_core, function(gs) intersect(gs, names(ranks)))
  fg <- tryCatch(
    fgsea::fgseaMultilevel(pathways = pw_for_fgsea, stats = ranks, minSize = 5, maxSize = 500),
    error = function(e)
      fgsea::fgseaSimple(pathways = pw_for_fgsea, stats = ranks, nperm = 20000, minSize = 5, maxSize = 500)
  )
  fg <- data.table::as.data.table(fg)[order(padj, pval)]
  cat("\n[Myocyte | orthogonalized-index adjustment | fgsea] d2_vs_mdx (reported with the updated consensus labels)\n")
  print(fg)
} else {
  message("[info] fgsea is not installed; running camera only")
}
# [Myocyte | orthogonalized-index adjustment | fgsea] d2_vs_mdx (reported with the updated consensus labels)
#                     pathway          padj   log2err         ES       NES  size  leadingEdge
# 1:             OxPhos_Ridge  2.026564e-25 1.3267161 -0.6902828 -3.253270   120 Uqcr11, ....
# 2:       ECM_Collagen_Ridge  3.184848e-09 0.7881868  0.7194894  2.301504    46 Fn1, Act....
# 3: Mito_Transl_Import_Ridge  1.842620e-05 0.5756103 -0.3885247 -1.840631   121 Atp5g1, ....
# 4:                TCA_Ridge  1.842620e-05 0.5756103 -0.6539155 -2.311371    28 Idh2, Md....


## Define strata (30th/70th percentiles of fiber_index within Myocyte)
q_low  <- quantile(myo_obj$fiber_index, 0.30, na.rm = TRUE)
q_high <- quantile(myo_obj$fiber_index, 0.70, na.rm = TRUE)
myo_obj$fiber_stratum <- ifelse(myo_obj$fiber_index <= q_low, "fast_like",
                                ifelse(myo_obj$fiber_index >= q_high, "slow_like", "middle"))
myo_obj$fiber_stratum <- NA_character_
for (s in unique(myo_obj$orig.ident)) {
  ix <- which(myo_obj$orig.ident == s)
  ql <- quantile(myo_obj$fiber_index[ix], 0.30, na.rm=TRUE)
  qh <- quantile(myo_obj$fiber_index[ix], 0.70, na.rm=TRUE)
  myo_obj$fiber_stratum[ix] <- ifelse(myo_obj$fiber_index[ix] <= ql, "fast_like",
                                      ifelse(myo_obj$fiber_index[ix] >= qh, "slow_like", "middle"))
}

## Exclude middle stratum and normalize sample IDs (hyphen -> underscore)
myo_sub <- subset(myo_obj, subset = fiber_stratum %in% c("slow_like","fast_like"))
myo_sub$sample_id <- gsub("-", "_", myo_sub$orig.ident)

agg_strata <- AggregateExpression(
  myo_sub, group.by = c("fiber_stratum","sample_id"),
  assays = "RNA", slot = "counts"
)
pb_strata <- agg_strata$RNA

## Fully standardize column names: "MyofiberSlow_<sample>" and "MyofiberFast_<sample>"
nm <- colnames(pb_strata)
nm <- gsub("-", "_", nm, fixed = TRUE)
nm <- sub("^slow[_:.]?like[_:.]", "MyofiberSlow_", nm, perl = TRUE)
nm <- sub("^fast[_:.]?like[_:.]", "MyofiberFast_", nm, perl = TRUE)
colnames(pb_strata) <- nm

.get_cols_by_prefix2 <- function(prefix, mat){
  pat <- paste0("^", gsub("([][\\^$.|?*+(){}\\\\])","\\\\\\1", prefix), "([_.:])")
  grep(pat, colnames(mat), value = TRUE, perl = TRUE)
}
stopifnot(length(.get_cols_by_prefix2("MyofiberSlow", pb_strata)) > 0)
stopifnot(length(.get_cols_by_prefix2("MyofiberFast", pb_strata)) > 0)

.run_pb_tests_matrix <- function(pb_mat, prefix,
                                 contrast = c("d2_vs_mdx","mdx_vs_wt","d2_vs_wt"),
                                 gene_sets = gene_sets_core) {
  contrast <- match.arg(contrast)
  cols <- .get_cols_by_prefix2(prefix, pb_mat); stopifnot(length(cols)>0)
  mat  <- as.matrix(pb_mat[, cols, drop = FALSE])
  labs <- sub("^[^_:]+[_:]", "", cols)  # "MyofiberSlow_d2_1" → "d2_1"
  grp  <- .parse_group(labs)
  
  y <- DGEList(counts = mat)
  keep <- filterByExpr(y, group = grp)
  y <- y[keep,, keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMMwsp")
  
  X  <- model.matrix(~0 + grp)
  v  <- voomWithQualityWeights(y, X, plot = FALSE)
  fit<- lmFit(v, X)
  Kk <- makeContrasts(
    mdx_vs_wt = grpmdx - grpwt,
    d2_vs_mdx = grpd2  - grpmdx,
    d2_vs_wt  = grpd2  - grpwt, levels = X
  )
  fitc <- eBayes(contrasts.fit(fit, Kk))
  
  genes_in2 <- rownames(fitc)
  idx2 <- lapply(gene_sets, function(gs){ ix <- match(intersect(gs, genes_in2), genes_in2); ix[!is.na(ix)] })
  idx2 <- idx2[sapply(idx2, length) >= 3]
  
  cam <- camera(v, index = idx2, design = X, contrast = Kk[, contrast])
  cam$FDR <- p.adjust(cam$PValue, "BH")
  cam[order(cam$FDR, cam$PValue), , drop = FALSE]
}

cat("\n[MyofiberSlow | camera] d2_vs_mdx\n")
print(.run_pb_tests_matrix(pb_strata, "MyofiberSlow", "d2_vs_mdx")[1:4, ])
#                          NGenes Direction           FDR
# OxPhos_Ridge                100      Down  3.398919e-17
# ECM_Collagen_Ridge           22        Up  6.883917e-05
# TCA_Ridge                    22      Down  7.978319e-04
# Mito_Transl_Import_Ridge     83      Down  1.091579e-02

cat("\n[MyofiberFast | camera] d2_vs_mdx\n")
print(.run_pb_tests_matrix(pb_strata, "MyofiberFast",  "d2_vs_mdx")[1:4, ])
#                          NGenes Direction           FDR
# ECM_Collagen_Ridge           18        Up 2.149899e-07
# OxPhos_Ridge                 88      Down 2.248188e-04
# Mito_Transl_Import_Ridge     58      Down 4.781476e-02
# TCA_Ridge                    15      Down 6.287234e-02



## Prerequisite: myo_obj contains slow_score / fast_score / fiber_index / orig.ident

## 1) Hard call (classify slow_score >= fast_score as Type-I-like)
myo_obj$myofiber_type_hard <- ifelse(
  myo_obj$slow_score >= myo_obj$fast_score, "TypeI_like", "TypeII_like"
)

# Total counts
cat("\n[Hard call] Total counts\n")
print(table(myo_obj$myofiber_type_hard))
# TypeI_like TypeII_like 
# 59         648

# By sample
cat("\n[Hard call] Counts by sample (orig.ident)\n")
print(as.data.frame.matrix(table(myo_obj$orig.ident, myo_obj$myofiber_type_hard)))

#     TypeI_like TypeII_like
#d2_1          3          34
#d2_2          5          23
#mdx1         18         142
#mdx2         11          78
#wt1           2         166
#wt2          20         205

# By genotype (if available)
if ("genotype" %in% colnames(myo_obj@meta.data)) {
  cat("\n[Hard call] Counts by genotype\n")
  print(as.data.frame.matrix(table(myo_obj$genotype, myo_obj$myofiber_type_hard)))
}
# [Hard call] Counts by genotype
#       TypeI_like TypeII_like
#WT            22         371
#mdx           29         220
#mdxD2          8          57

## 2) Quantile call (split within each sample at the median: TypeI_like / TypeII_like)
myo_obj$myofiber_type_q50 <- NA_character_
for (s in unique(myo_obj$orig.ident)) {
  ix <- which(myo_obj$orig.ident == s)
  if (length(ix) == 0) next
  cut <- stats::quantile(myo_obj$fiber_index[ix], 0.50, na.rm = TRUE)  # median
  myo_obj$myofiber_type_q50[ix] <- ifelse(myo_obj$fiber_index[ix] >= cut, "TypeI_like", "TypeII_like")
}

cat("\n[Quantile (within-sample median) call] Total counts\n")
print(table(myo_obj$myofiber_type_q50))
# TypeI_like TypeII_like 
#   355         352

cat("\n[Quantile call] Counts by sample (orig.ident)\n")
print(as.data.frame.matrix(table(myo_obj$orig.ident, myo_obj$myofiber_type_q50)))
#      TypeI_like TypeII_like
# d2_1         19          18
# d2_2         14          14
# mdx1         80          80
# mdx2         45          44
# wt1          84          84
# wt2         113         112


if ("genotype" %in% colnames(myo_obj@meta.data)) {
  cat("\n[Quantile call] Counts by genotype\n")
  print(as.data.frame.matrix(table(myo_obj$genotype, myo_obj$myofiber_type_q50)))
}
# [Quantile call] Counts by genotype
#       TypeI_like TypeII_like
#WT           197         196
#mdx          125         124
#mdxD2         33          32

## (Optional) More stringent classification using only the 30/70% tails as confident calls (excluding the middle stratum)
myo_obj$myofiber_type_q30_70 <- NA_character_
for (s in unique(myo_obj$orig.ident)) {
  ix <- which(myo_obj$orig.ident == s)
  if (length(ix) == 0) next
  ql <- stats::quantile(myo_obj$fiber_index[ix], 0.30, na.rm = TRUE)
  qh <- stats::quantile(myo_obj$fiber_index[ix], 0.70, na.rm = TRUE)
  myo_obj$myofiber_type_q30_70[ix] <- ifelse(
    myo_obj$fiber_index[ix] >= qh, "TypeI_like",
    ifelse(myo_obj$fiber_index[ix] <= ql, "TypeII_like", "middle")
  )
}

cat("\n[Quantile 30/70] Confident TypeI/II only (middle excluded) — Total counts\n")
print(table(myo_obj$myofiber_type_q30_70))
# middle  TypeI_like TypeII_like 
#  279         214         214


cat("\n[Quantile 30/70] by sample (orig.ident)\n")
print(as.data.frame.matrix(table(myo_obj$orig.ident, myo_obj$myofiber_type_q30_70)))
# [Quantile 30/70] by sample (orig.ident)
#      middle TypeI_like TypeII_like
#d2_1     15         11          11
#d2_2     10          9           9
#mdx1     64         48          48
#mdx2     35         27          27
#wt1      66         51          51
#wt2      89         68          68



suppressPackageStartupMessages(library(mclust)); library(dplyr)
stopifnot(exists("myo_obj"), "fiber_index" %in% colnames(myo_obj@meta.data))
meta <- myo_obj@meta.data
meta$group <- if("genotype"%in%names(myo_obj@meta.data)) sub("^mdx.?d2$","d2", sub("^mdx$","mdx", sub("^wt$","wt", tolower(meta$genotype)))) else sub("^(wt|mdx|d2).*","\\1", tolower(meta$orig.ident))
x_wt <- na.omit(meta$fiber_index[meta$group=="wt"]); fit <- Mclust(x_wt, G=2, verbose=FALSE)

pi <- fit$parameters$pro; mu <- fit$parameters$mean
s2 <- fit$parameters$variance$sigmasq; if(length(s2)==1) s2 <- rep(s2,2)
f  <- function(z) pi[1]*dnorm(z,mu[1],sqrt(s2[1])) - pi[2]*dnorm(z,mu[2],sqrt(s2[2]))
rng <- range(mu); t  <- if (sign(f(rng[1]))==sign(f(rng[2]))) mean(mu) else uniroot(f, rng)$root
hi <- which.max(mu); pred <- predict(fit, newdata = meta$fiber_index)

meta$TypeI_hard <- !is.na(meta$fiber_index) & (meta$fiber_index >= t)
meta$TypeI_soft <- pred$z[, hi]

summ <- meta %>%
  dplyr::filter(group %in% c("wt","mdx","d2")) %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n=dplyr::n(),
                   hard=sum(TypeI_hard,na.rm=TRUE),
                   frac_hard=hard/n,
                   frac_soft=mean(TypeI_soft,na.rm=TRUE),
                   .groups="drop") %>%
  dplyr::mutate(se_hard=sqrt(frac_hard*(1-frac_hard)/n),
                ci_lo=pmax(0, frac_hard-1.96*se_hard),
                ci_hi=pmin(1, frac_hard+1.96*se_hard))
print(list(threshold_t=t, summary=summ))
# A tibble: 3 × 7
#  group     n  frac_hard frac_soft     se ci_lo ci_hi
# <chr> <int>     <dbl>     <dbl>  <dbl> <dbl> <dbl>
# 1 d2       65     0.308    0.0467 0.0572 0.195 0.420
# 2 mdx     249     0.197    0.0356 0.0252 0.147 0.246
# 3 wt      393     0.201    0.0205 0.0202 0.161 0.241

ref <- summ[summ$group=="wt", c("frac_hard","n")]; vref <- ref$frac_hard*(1-ref$frac_hard)/ref$n
cmp <- summ %>% dplyr::filter(group!="wt") %>%
  dplyr::mutate(diff_hard=frac_hard-ref$frac_hard,
                se_diff=sqrt(frac_hard*(1-frac_hard)/n + vref),
                diff_lo=diff_hard-1.96*se_diff, diff_hi=diff_hard+1.96*se_diff) %>%
  dplyr::select(group, n, diff_hard, diff_lo, diff_hi)
print(cmp)  # Both mdx and d2 are reported
# A tibble: 2 × 5
# group        n     diff diff_lo diff_hi
# <chr>    <int>    <dbl>   <dbl>   <dbl>
#   1 d2       65   0.107   -0.0123  0.226 
#   2 mdx     249  -0.00423 -0.0675  0.0591

suppressPackageStartupMessages({library(mclust);library(dplyr)})
stopifnot(exists("myo_obj"), "fiber_index"%in%colnames(myo_obj@meta.data))
md <- myo_obj@meta.data
grp <- if ("genotype" %in% names(md)) {
  g <- tolower(as.character(md$genotype))
  sub("^mdx.?d2$","d2", sub("^mdx$","mdx", sub("^wt$","wt", g)))
} else {
  sub("^(wt|mdx|d2).*", "\\1", tolower(as.character(md$orig.ident)))
}
x  <- md$fiber_index; x_wt <- na.omit(x[grp=="wt"])

## A) Mixture-model crossing point (posterior = 0.5)
fit <- Mclust(x_wt, G=2, verbose=FALSE); pi <- fit$parameters$pro; mu <- fit$parameters$mean
s2 <- fit$parameters$variance$sigmasq; if(length(s2)==1) s2 <- rep(s2,2)
f  <- function(z) pi[1]*dnorm(z,mu[1],sqrt(s2[1])) - pi[2]*dnorm(z,mu[2],sqrt(s2[2]))
tA <- {rng<-range(mu); if(sign(f(rng[1]))==sign(f(rng[2]))) mean(mu) else uniroot(f, rng)$root}

## B) Otsu threshold (maximize separation in the histogram)
otsu <- function(v, nb=256){v<-v[is.finite(v)]; h<-hist(v, nb, plot=FALSE); p<-h$counts/sum(h$counts); w<-cumsum(p); m<-cumsum(p*h$mids)
mu<-m[length(m)]; s2<- (mu*w - m)^2 /(w*(1-w)+1e-12); thr <- h$mids[which.max(replace(s2, !is.finite(s2), -Inf))]; thr}
tB <- otsu(x_wt)

## C) Archetype ROC (use the lower and upper 20% of WT as anchors)
q <- quantile(x_wt, c(.2,.8), na.rm=TRUE); lab <- x_wt; y <- ifelse(lab>=q[2],1, ifelse(lab<=q[1],0, NA)); z <- lab[!is.na(y)]; y <- y[!is.na(y)]
cand <- sort(unique(z)); J <- sapply(cand, function(c) {tp<-mean(y[z>=c]==1); tn<-mean(y[z<c]==0); tp+tn-1})
tC <- cand[which.max(J)]

## Final threshold t* (median) and posterior-based soft assignment
t_star <- median(c(tA,tB,tC), na.rm=TRUE)
post <- tryCatch({pz<-predict(fit, newdata = x)$z; pz[,which.max(mu)]}, error=function(e) rep(NA_real_, length(x)))

## Summary: Type-I-like nuclear fractions (hard/soft) with 95% Wald confidence intervals
summ <- data.frame(group=grp, x=x, pI=post) %>%
  filter(group%in%c("wt","mdx","d2") & is.finite(x)) %>%
  group_by(group) %>%
  summarise(n=n(),
            frac_hard=mean(x>=t_star),
            frac_soft=mean(pI, na.rm=TRUE),
            se=sqrt(frac_hard*(1-frac_hard)/n), .groups="drop") %>%
  mutate(ci_lo=pmax(0, frac_hard-1.96*se), ci_hi=pmin(1, frac_hard+1.96*se))

## Difference from WT (hard call) with 95% confidence intervals
ref <- summ[summ$group=="wt", c("frac_hard","n")]; vref <- with(ref, frac_hard*(1-frac_hard)/n)
cmp <- summ %>% filter(group!="wt") %>%
  mutate(diff=frac_hard-ref$frac_hard, se_diff=sqrt(frac_hard*(1-frac_hard)/n + vref),
         diff_lo=diff-1.96*se_diff, diff_hi=diff+1.96*se_diff) %>%
  select(group, n, diff, diff_lo, diff_hi)

print(list(thresholds=list(tA=tA,tB=tB,tC=tC,t_star=t_star),
           summary=summ %>% arrange(match(group,c("wt","mdx","d2"))),
           diff_vs_WT=cmp))

# group       n frac_hard frac_soft     se ci_lo ci_hi
# 1 wt      393     0.201    0.0205 0.0202 0.161 0.241
# 2 mdx     249     0.197    0.0356 0.0252 0.147 0.246
# 3 d2       65     0.308    0.0467 0.0572 0.195 0.420

# group      n     diff diff_lo diff_hi
# 1 d2       65  0.107   -0.0123  0.226 
# 2 mdx     249 -0.00423 -0.0675  0.0591

## After the main analysis block is completed
#saveRDS(object   = muscle_sc,file = "C:/DMD_project/rds/muscle_sc_final.rds",compress = "gzip")



                 
## Ortholog projection analysis------------------------------------------------------------

## 0) User settings -------------------------------------------------------------

proj_root_xsp <- "C:/DMD_project"

cache_dir_bulk_xsp <- file.path(proj_root_xsp, "cache_bulk_revision")
rds_dir_mouse_xsp  <- file.path(proj_root_xsp, "rds")
res_dir_xsp        <- file.path(proj_root_xsp, "results_cross_species_concordance")
fig_dir_xsp        <- file.path(res_dir_xsp, "figures")

dir.create(res_dir_xsp, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir_xsp, showWarnings = FALSE, recursive = TRUE)

# Human inputs generated by Script 1.
human_stage6_path_xsp <- file.path(cache_dir_bulk_xsp, "stage6_combat_only_bulk.rds")

# Optional fallback if the stage6 cache is unavailable.
# The stage6 cache is preferred because it preserves the exact primary DEG object.
additional_file2_candidates_xsp <- c(
  file.path(proj_root_xsp, "Additional_file_2_Primary_DEG_and_sensitivity_analysis_table.xlsx"),
  file.path(proj_root_xsp, "Additional_file_2_Primary_DEG_and_sensitivity_analysis_table(4).xlsx"),
  file.path(proj_root_xsp, "results_bulk_revision", "Additional_file_2_Primary_DEG_and_sensitivity_analysis_table.xlsx")
)

# TSV fallback generated by Script 1 when write_deg_table_bulk <- TRUE.
human_deg_tsv_candidates_xsp <- c(
  file.path(proj_root_xsp, "results_bulk_revision", "deg_bulk.tsv"),
  file.path(proj_root_xsp, "deg_bulk.tsv")
)

use_additional_file2_if_stage6_missing_xsp <- TRUE
use_deg_tsv_if_stage6_missing_xsp <- TRUE

# The script always writes TSV and RDS outputs. If openxlsx is already installed,
# it also writes an XLSX workbook. openxlsx is intentionally optional so that
# package-download problems do not block the concordance analysis.
write_xlsx_if_openxlsx_available_xsp <- TRUE

# Mouse inputs generated by Script 2.
mouse_rds_candidates_xsp <- c(
  file.path(rds_dir_mouse_xsp, "muscle_sc_final.rds"),
  file.path(rds_dir_mouse_xsp, "muscle_sc_after_dotplot.rds"),
  file.path(rds_dir_mouse_xsp, "muscle_sc_proto.rds")
)

# Whether to run layer-wise scDblFinder when scDblFinder.class is absent.
# If the loaded object is muscle_sc_final.rds after doublet/RBC exclusion, this is usually unnecessary.
run_scdblfinder_if_missing_xsp <- TRUE
rerun_scdblfinder_even_if_present_xsp <- FALSE

# Analysis parameters.
seed_xsp <- 123
human_deg_fdr_cut_xsp <- 0.05
human_abs_logfc_min_xsp <- 0
ortholog_min_support_xsp <- 2
fgsea_min_size_xsp <- 10
fgsea_max_size_xsp <- 1000
camera_min_size_xsp <- 5
contrast_name_xsp <- "mdxD2_vs_mdx"

# Require at least this many pseudobulk samples per contrasted mouse genotype.
min_samples_per_contrast_group_xsp <- 2

# Package installation switch. Keep FALSE when running in a curated local R environment.
auto_install_packages_xsp <- FALSE



## 1) Package management and general utilities ----------------------------------

cran_pkgs_xsp <- c(
  "dplyr", "tibble", "tidyr", "stringr", "ggplot2", "data.table",
  "Matrix", "babelgene", "scales"
)

bioc_pkgs_xsp <- c(
  "Seurat", "SeuratObject", "SingleCellExperiment", "scDblFinder",
  "edgeR", "limma", "fgsea"
)

optional_cran_pkgs_xsp <- c("openxlsx")

if (isTRUE(auto_install_packages_xsp)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  missing_cran_xsp <- cran_pkgs_xsp[!vapply(cran_pkgs_xsp, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_cran_xsp) > 0) {
    install.packages(missing_cran_xsp, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
  missing_bioc_xsp <- bioc_pkgs_xsp[!vapply(bioc_pkgs_xsp, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc_xsp) > 0) {
    BiocManager::install(missing_bioc_xsp, ask = FALSE, update = FALSE)
  }
}

required_pkgs_xsp <- unique(c(cran_pkgs_xsp, bioc_pkgs_xsp))
missing_pkgs_xsp <- required_pkgs_xsp[!vapply(required_pkgs_xsp, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs_xsp) > 0) {
  stop(
    "Missing packages for cross-species concordance analysis:\n  ",
    paste(missing_pkgs_xsp, collapse = "\n  "),
    "\nInstall them first or set auto_install_packages_xsp <- TRUE."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(data.table)
  library(Matrix)
  library(babelgene)
  library(scales)
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(edgeR)
  library(limma)
  library(fgsea)
})

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(seed_xsp)

ts_xsp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
say_xsp <- function(...) message(ts_xsp(), " | ", ...)

first_existing_xsp <- function(paths_xsp) {
  hit_xsp <- paths_xsp[file.exists(paths_xsp)]
  if (length(hit_xsp) == 0) return(NULL)
  hit_xsp[1]
}

first_col_xsp <- function(df_xsp, candidates_xsp) {
  hit_xsp <- candidates_xsp[candidates_xsp %in% colnames(df_xsp)]
  if (length(hit_xsp) == 0) return(NULL)
  hit_xsp[1]
}

write_tsv_xsp <- function(df_xsp, file_xsp) {
  utils::write.table(
    df_xsp,
    file = file_xsp,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

safe_neglog10_xsp <- function(p_xsp) {
  -log10(pmax(as.numeric(p_xsp), .Machine$double.xmin))
}

escape_regex_xsp <- function(x_xsp) {
  gsub("([][\\^$.|?*+(){}\\\\])", "\\\\\\1", x_xsp)
}

collapse_rank_by_gene_xsp <- function(stats_xsp) {
  stats_xsp <- stats_xsp[is.finite(stats_xsp) & !is.na(names(stats_xsp)) & names(stats_xsp) != ""]
  if (length(stats_xsp) == 0) return(stats_xsp)
  tmp_xsp <- data.frame(
    gene = names(stats_xsp),
    stat = as.numeric(stats_xsp),
    stringsAsFactors = FALSE
  )
  tmp_xsp <- tmp_xsp |>
    dplyr::group_by(gene) |>
    dplyr::arrange(dplyr::desc(abs(stat)), .by_group = TRUE) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  out_xsp <- tmp_xsp$stat
  names(out_xsp) <- tmp_xsp$gene
  sort(out_xsp, decreasing = TRUE)
}

standardize_mouse_sample_xsp <- function(x_xsp) {
  x_xsp <- as.character(x_xsp)
  x_xsp <- gsub("-", "_", x_xsp)
  x_xsp <- gsub("^mdxD2", "d2", x_xsp, ignore.case = TRUE)
  x_xsp <- gsub("^mdxd2", "d2", x_xsp, ignore.case = TRUE)
  x_xsp
}

infer_mouse_genotype_short_xsp <- function(sample_or_genotype_xsp) {
  s_xsp <- tolower(as.character(sample_or_genotype_xsp))
  s_xsp <- gsub("-", "_", s_xsp)
  out_xsp <- rep(NA_character_, length(s_xsp))
  out_xsp[grepl("^(wt|wild)", s_xsp)] <- "wt"
  out_xsp[grepl("^(d2|mdxd2|mdx_d2|mdx\\.d2|mdx-d2)", s_xsp)] <- "d2"
  out_xsp[is.na(out_xsp) & grepl("^mdx", s_xsp)] <- "mdx"
  factor(out_xsp, levels = c("wt", "mdx", "d2"))
}

pretty_signature_xsp <- function(x_xsp) {
  dplyr::case_when(
    x_xsp == "Human_Late_up_DEG_orthologs"  ~ "Human Late-up DEG orthologs",
    x_xsp == "Human_Early_up_DEG_orthologs" ~ "Human Early-up DEG orthologs",
    TRUE ~ x_xsp
  )
}



## 2) Load human DEG signature and convert to mouse orthologs -------------------

load_human_deg_from_stage6_xsp <- function(stage6_path_xsp) {
  if (!file.exists(stage6_path_xsp)) return(NULL)
  say_xsp("Loading human bulk cache: ", stage6_path_xsp)
  st_xsp <- readRDS(stage6_path_xsp)
  if (!"deg_bulk" %in% names(st_xsp)) {
    stop("stage6_combat_only_bulk.rds was found, but it does not contain deg_bulk.")
  }
  deg_xsp <- as.data.frame(st_xsp$deg_bulk)
  deg_xsp$EntrezID <- sub("_.*", "", rownames(deg_xsp))
  deg_xsp$Source <- "stage6_combat_only_bulk.rds"
  deg_xsp
}

load_human_deg_from_tsv_xsp <- function(candidates_xsp) {
  path_xsp <- first_existing_xsp(candidates_xsp)
  if (is.null(path_xsp)) return(NULL)
  say_xsp("Using DEG TSV fallback for human DEG table: ", path_xsp)
  deg_xsp <- utils::read.delim(path_xsp, stringsAsFactors = FALSE, check.names = FALSE)
  
  entrez_col_xsp <- first_col_xsp(deg_xsp, c("EntrezID", "EntrezID_bulk", "ENTREZID", "GeneID"))
  if (is.null(entrez_col_xsp)) {
    stop("The DEG TSV fallback does not contain an EntrezID-like column.")
  }
  deg_xsp$EntrezID <- as.character(deg_xsp[[entrez_col_xsp]])
  deg_xsp$Source <- basename(path_xsp)
  rownames(deg_xsp) <- make.unique(deg_xsp$EntrezID)
  deg_xsp
}

load_human_deg_from_additional_file2_xsp <- function(candidates_xsp) {
  path_xsp <- first_existing_xsp(candidates_xsp)
  if (is.null(path_xsp)) return(NULL)
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    say_xsp("Additional file 2 was found, but openxlsx is not installed; skipping XLSX fallback.")
    return(NULL)
  }
  say_xsp("Using Additional file 2 fallback for human DEG table: ", path_xsp)
  deg_xsp <- openxlsx::read.xlsx(path_xsp, sheet = "All_DEGs")
  if (!"EntrezID" %in% colnames(deg_xsp)) {
    stop("The All_DEGs sheet does not contain an EntrezID column.")
  }
  if ("FDR" %in% colnames(deg_xsp) && !"adj.P.Val" %in% colnames(deg_xsp)) {
    deg_xsp$adj.P.Val <- deg_xsp$FDR
  }
  if ("Pvalue" %in% colnames(deg_xsp) && !"P.Value" %in% colnames(deg_xsp)) {
    deg_xsp$P.Value <- deg_xsp$Pvalue
  }
  rownames(deg_xsp) <- as.character(deg_xsp$EntrezID)
  deg_xsp$Source <- "Additional_file_2_All_DEGs"
  deg_xsp
}

standardize_human_deg_table_xsp <- function(deg_xsp) {
  logfc_col_xsp <- first_col_xsp(deg_xsp, c("logFC", "LFC", "log2FC"))
  fdr_col_xsp   <- first_col_xsp(deg_xsp, c("adj.P.Val", "FDR", "padj", "p.adjust", "qvalue"))
  p_col_xsp     <- first_col_xsp(deg_xsp, c("P.Value", "Pvalue", "pvalue", "P"))
  entrez_col_xsp <- first_col_xsp(deg_xsp, c("EntrezID", "EntrezID_bulk", "ENTREZID", "GeneID"))
  symbol_col_xsp <- first_col_xsp(deg_xsp, c("SYMBOL", "Symbol", "GeneSymbol", "Gene_SYMBOL", "SYMBOL_bulk"))
  
  if (is.null(entrez_col_xsp)) stop("No EntrezID-like column was found in the human DEG table.")
  if (is.null(logfc_col_xsp)) stop("No logFC-like column was found in the human DEG table.")
  if (is.null(fdr_col_xsp))   stop("No FDR-like column was found in the human DEG table.")
  
  out_xsp <- deg_xsp |>
    dplyr::mutate(
      EntrezID = as.character(.data[[entrez_col_xsp]]),
      logFC = as.numeric(.data[[logfc_col_xsp]]),
      FDR = as.numeric(.data[[fdr_col_xsp]]),
      Pvalue = if (!is.null(p_col_xsp)) as.numeric(.data[[p_col_xsp]]) else NA_real_,
      Human_SYMBOL_input = if (!is.null(symbol_col_xsp)) as.character(.data[[symbol_col_xsp]]) else NA_character_
    ) |>
    dplyr::filter(!is.na(EntrezID), EntrezID != "", !is.na(logFC), !is.na(FDR)) |>
    dplyr::filter(FDR < human_deg_fdr_cut_xsp, abs(logFC) >= human_abs_logfc_min_xsp) |>
    dplyr::mutate(
      Human_direction = dplyr::if_else(logFC > 0, "Late_up", "Early_up"),
      Signature = dplyr::if_else(
        logFC > 0,
        "Human_Late_up_DEG_orthologs",
        "Human_Early_up_DEG_orthologs"
      )
    )
  
  out_xsp
}

map_human_deg_entrez_to_mouse_orthologs_xsp <- function(human_deg_xsp, min_support_xsp = 2) {
  human_entrez_xsp <- unique(as.character(human_deg_xsp$EntrezID))
  human_entrez_xsp <- human_entrez_xsp[!is.na(human_entrez_xsp) & human_entrez_xsp != ""]
  if (length(human_entrez_xsp) == 0) stop("No human Entrez IDs were available for ortholog mapping.")
  
  say_xsp("Mapping human Entrez IDs to mouse orthologs with babelgene: ", length(human_entrez_xsp), " Entrez IDs")
  
  ort_xsp <- tryCatch(
    babelgene::orthologs(
      genes = human_entrez_xsp,
      species = "mouse",
      human = TRUE,
      min_support = min_support_xsp,
      top = TRUE
    ),
    error = function(e1) {
      say_xsp("babelgene call with min_support/top failed; retrying with the minimal argument set.")
      tryCatch(
        babelgene::orthologs(genes = human_entrez_xsp, species = "mouse", human = TRUE),
        error = function(e2) stop("babelgene ortholog mapping failed: ", conditionMessage(e2))
      )
    }
  )
  
  if (!is.data.frame(ort_xsp) || nrow(ort_xsp) == 0) {
    stop("babelgene returned no ortholog mappings.")
  }
  
  human_entrez_col_xsp <- first_col_xsp(ort_xsp, c("human_entrez", "human_entrez_id", "human_ncbi_gene_id"))
  human_symbol_col_xsp <- first_col_xsp(ort_xsp, c("human_symbol", "human_gene", "human_gene_symbol", "input", "query"))
  mouse_col_xsp <- first_col_xsp(ort_xsp, c("symbol", "mouse_symbol", "gene_symbol", "ortholog_symbol"))
  mouse_entrez_col_xsp <- first_col_xsp(ort_xsp, c("entrez", "mouse_entrez", "mouse_entrez_id"))
  support_col_xsp <- first_col_xsp(ort_xsp, c("support_n", "support", "n_support", "num_ortholog_sources"))
  
  if (is.null(human_entrez_col_xsp)) {
    stop("Could not identify the human Entrez column in the babelgene output. Columns: ", paste(colnames(ort_xsp), collapse = ", "))
  }
  if (is.null(mouse_col_xsp)) {
    stop("Could not identify the mouse-symbol column in the babelgene output. Columns: ", paste(colnames(ort_xsp), collapse = ", "))
  }
  
  out_xsp <- data.frame(
    EntrezID = as.character(ort_xsp[[human_entrez_col_xsp]]),
    Human_SYMBOL_ortholog = if (!is.null(human_symbol_col_xsp)) as.character(ort_xsp[[human_symbol_col_xsp]]) else NA_character_,
    Mouse_SYMBOL = as.character(ort_xsp[[mouse_col_xsp]]),
    Mouse_EntrezID = if (!is.null(mouse_entrez_col_xsp)) as.character(ort_xsp[[mouse_entrez_col_xsp]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(support_col_xsp)) {
    support_numeric_xsp <- suppressWarnings(as.numeric(ort_xsp[[support_col_xsp]]))
    out_xsp$Ortholog_support <- support_numeric_xsp
    if (any(is.finite(support_numeric_xsp))) {
      out_xsp <- out_xsp[is.na(out_xsp$Ortholog_support) | out_xsp$Ortholog_support >= min_support_xsp, , drop = FALSE]
    }
  } else {
    out_xsp$Ortholog_support <- NA_real_
  }
  
  out_xsp |>
    dplyr::filter(!is.na(EntrezID), EntrezID != "", !is.na(Mouse_SYMBOL), Mouse_SYMBOL != "") |>
    dplyr::distinct(EntrezID, Mouse_SYMBOL, .keep_all = TRUE)
}

human_deg_raw_xsp <- load_human_deg_from_stage6_xsp(human_stage6_path_xsp)
if (is.null(human_deg_raw_xsp) && isTRUE(use_deg_tsv_if_stage6_missing_xsp)) {
  human_deg_raw_xsp <- load_human_deg_from_tsv_xsp(human_deg_tsv_candidates_xsp)
}
if (is.null(human_deg_raw_xsp)) {
  if (!isTRUE(use_additional_file2_if_stage6_missing_xsp)) {
    stop("Human stage6 cache was not found: ", human_stage6_path_xsp)
  }
  human_deg_raw_xsp <- load_human_deg_from_additional_file2_xsp(additional_file2_candidates_xsp)
}
if (is.null(human_deg_raw_xsp)) {
  stop("No human DEG input was found. Expected stage6 cache, DEG TSV fallback, or Additional file 2 fallback.")
}

human_deg_xsp <- standardize_human_deg_table_xsp(human_deg_raw_xsp)
say_xsp("Human DEGs retained for signature construction: ", nrow(human_deg_xsp))
say_xsp("  Late-up:  ", sum(human_deg_xsp$Human_direction == "Late_up"))
say_xsp("  Early-up: ", sum(human_deg_xsp$Human_direction == "Early_up"))

ortholog_tbl_xsp <- map_human_deg_entrez_to_mouse_orthologs_xsp(
  human_deg_xsp = human_deg_xsp,
  min_support_xsp = ortholog_min_support_xsp
)

human_to_mouse_sig_xsp <- human_deg_xsp |>
  dplyr::left_join(ortholog_tbl_xsp, by = "EntrezID") |>
  dplyr::mutate(
    Human_SYMBOL = dplyr::coalesce(Human_SYMBOL_input, Human_SYMBOL_ortholog, EntrezID)
  ) |>
  dplyr::filter(!is.na(Mouse_SYMBOL), Mouse_SYMBOL != "") |>
  dplyr::distinct(Signature, Human_direction, EntrezID, Human_SYMBOL, Mouse_SYMBOL, .keep_all = TRUE)

gene_sets_mouse_xsp <- split(human_to_mouse_sig_xsp$Mouse_SYMBOL, human_to_mouse_sig_xsp$Signature)
gene_sets_mouse_xsp <- lapply(gene_sets_mouse_xsp, function(x_xsp) sort(unique(x_xsp[!is.na(x_xsp) & x_xsp != ""])))

required_signatures_xsp <- c("Human_Late_up_DEG_orthologs", "Human_Early_up_DEG_orthologs")
if (!all(required_signatures_xsp %in% names(gene_sets_mouse_xsp))) {
  stop("Both Late-up and Early-up human ortholog signatures were not constructed.")
}
for (nm_xsp in names(gene_sets_mouse_xsp)) {
  say_xsp(nm_xsp, ": ", length(gene_sets_mouse_xsp[[nm_xsp]]), " mouse ortholog genes")
  if (length(gene_sets_mouse_xsp[[nm_xsp]]) < fgsea_min_size_xsp) {
    stop("Signature too small after ortholog mapping: ", nm_xsp)
  }
}



## 3) Load mouse Seurat object and generate cell-class pseudobulk counts --------

run_layerwise_scdblfinder_xsp <- function(obj_xsp) {
  if (!inherits(obj_xsp, "Seurat")) stop("Input object is not a Seurat object.")
  if (!"RNA" %in% names(obj_xsp@assays)) stop("The Seurat object does not contain an RNA assay.")
  Seurat::DefaultAssay(obj_xsp) <- "RNA"
  
  layers_all_xsp <- tryCatch(SeuratObject::Layers(obj_xsp[["RNA"]]), error = function(e) character(0))
  count_layers_xsp <- grep("^counts(\\.|$)", layers_all_xsp, value = TRUE)
  count_layers_xsp <- setdiff(count_layers_xsp, "counts")
  
  if (length(count_layers_xsp) == 0) {
    say_xsp("No sample-specific counts.* layers were found. scDblFinder layer-wise rerun is skipped.")
    return(obj_xsp)
  }
  
  all_cells_xsp <- colnames(obj_xsp)
  dbl_class_xsp <- setNames(rep(NA_character_, length(all_cells_xsp)), all_cells_xsp)
  dbl_score_xsp <- setNames(rep(NA_real_, length(all_cells_xsp)), all_cells_xsp)
  summary_list_xsp <- list()
  
  say_xsp("Running layer-wise scDblFinder on ", length(count_layers_xsp), " count layers")
  for (lay_xsp in count_layers_xsp) {
    mat_xsp <- obj_xsp[["RNA"]]@layers[[lay_xsp]]
    if (is.null(mat_xsp) || ncol(mat_xsp) == 0) next
    cells_lay_xsp <- tryCatch(SeuratObject::Cells(obj_xsp, layer = lay_xsp), error = function(e) colnames(mat_xsp))
    if (length(cells_lay_xsp) != ncol(mat_xsp)) {
      warning("Layer cell names do not match matrix columns for layer: ", lay_xsp)
      next
    }
    set.seed(seed_xsp)
    sce_xsp <- SingleCellExperiment::SingleCellExperiment(list(counts = mat_xsp))
    sce_xsp <- scDblFinder::scDblFinder(sce_xsp)
    calls_xsp <- as.character(SummarizedExperiment::colData(sce_xsp)$scDblFinder.class)
    scores_xsp <- as.numeric(SummarizedExperiment::colData(sce_xsp)$scDblFinder.score)
    names(calls_xsp) <- cells_lay_xsp
    names(scores_xsp) <- cells_lay_xsp
    dbl_class_xsp[cells_lay_xsp] <- calls_xsp[cells_lay_xsp]
    dbl_score_xsp[cells_lay_xsp] <- scores_xsp[cells_lay_xsp]
    summary_list_xsp[[lay_xsp]] <- data.frame(
      layer = lay_xsp,
      cells = length(cells_lay_xsp),
      doublets = sum(calls_xsp == "doublet", na.rm = TRUE),
      doublet_ratio = mean(calls_xsp == "doublet", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  dbl_class_xsp[is.na(dbl_class_xsp)] <- "singlet"
  obj_xsp$scDblFinder.class <- dbl_class_xsp[colnames(obj_xsp)]
  obj_xsp$scDblFinder.score <- dbl_score_xsp[colnames(obj_xsp)]
  
  summary_df_xsp <- dplyr::bind_rows(summary_list_xsp)
  if (nrow(summary_df_xsp) > 0) {
    write_tsv_xsp(summary_df_xsp, file.path(res_dir_xsp, "mouse_layerwise_scdblfinder_summary.tsv"))
  }
  obj_xsp
}

join_rna_layers_if_needed_xsp <- function(obj_xsp) {
  Seurat::DefaultAssay(obj_xsp) <- "RNA"
  layers_all_xsp <- tryCatch(SeuratObject::Layers(obj_xsp[["RNA"]]), error = function(e) character(0))
  count_layers_xsp <- grep("^counts", layers_all_xsp, value = TRUE)
  if (length(count_layers_xsp) > 1) {
    say_xsp("Joining RNA layers before pseudobulk count extraction.")
    if (exists("JoinLayers", where = asNamespace("Seurat"), inherits = FALSE)) {
      obj_xsp <- Seurat::JoinLayers(obj_xsp, assay = "RNA")
    } else if (exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
      obj_xsp <- SeuratObject::JoinLayers(obj_xsp, assay = "RNA")
    } else {
      warning("JoinLayers was not found in Seurat or SeuratObject; count extraction may fail if multiple layers remain.")
    }
  }
  obj_xsp
}

get_rna_counts_xsp <- function(obj_xsp) {
  Seurat::DefaultAssay(obj_xsp) <- "RNA"
  counts_xsp <- tryCatch(
    Seurat::GetAssayData(obj_xsp, assay = "RNA", layer = "counts"),
    error = function(e1) {
      tryCatch(
        Seurat::GetAssayData(obj_xsp, assay = "RNA", slot = "counts"),
        error = function(e2) stop("Could not extract RNA counts from the Seurat object: ", conditionMessage(e2))
      )
    }
  )
  counts_xsp
}

make_pseudobulk_counts_xsp <- function(obj_xsp) {
  if (!"celltype" %in% colnames(obj_xsp@meta.data)) {
    stop("The mouse Seurat object must contain a 'celltype' metadata column.")
  }
  if (!"orig.ident" %in% colnames(obj_xsp@meta.data)) {
    stop("The mouse Seurat object must contain an 'orig.ident' metadata column.")
  }
  
  if (!"genotype" %in% colnames(obj_xsp@meta.data)) {
    say_xsp("No genotype metadata column was found; inferring genotype from orig.ident.")
    obj_xsp$genotype <- dplyr::case_when(
      as.character(infer_mouse_genotype_short_xsp(obj_xsp$orig.ident)) == "wt"  ~ "WT",
      as.character(infer_mouse_genotype_short_xsp(obj_xsp$orig.ident)) == "mdx" ~ "mdx",
      as.character(infer_mouse_genotype_short_xsp(obj_xsp$orig.ident)) == "d2"  ~ "mdxD2",
      TRUE ~ NA_character_
    )
  }
  
  if ((rerun_scdblfinder_even_if_present_xsp || !"scDblFinder.class" %in% colnames(obj_xsp@meta.data)) &&
      isTRUE(run_scdblfinder_if_missing_xsp)) {
    obj_xsp <- run_layerwise_scdblfinder_xsp(obj_xsp)
  }
  
  n_before_xsp <- ncol(obj_xsp)
  keep_xsp <- rep(TRUE, n_before_xsp)
  names(keep_xsp) <- colnames(obj_xsp)
  
  if ("celltype" %in% colnames(obj_xsp@meta.data)) {
    keep_xsp <- keep_xsp & !(as.character(obj_xsp$celltype) %in% c("Erythrocyte", "erythrocyte", "RBC", "Red blood cell"))
  }
  if ("scDblFinder.class" %in% colnames(obj_xsp@meta.data)) {
    keep_xsp <- keep_xsp & !(as.character(obj_xsp$scDblFinder.class) == "doublet")
  }
  
  obj_xsp <- subset(obj_xsp, cells = names(keep_xsp)[keep_xsp])
  obj_xsp$celltype <- droplevels(factor(obj_xsp$celltype))
  Seurat::Idents(obj_xsp) <- "celltype"
  n_after_xsp <- ncol(obj_xsp)
  say_xsp("Mouse cells retained for pseudobulk: ", n_after_xsp, " / ", n_before_xsp)
  
  obj_xsp <- join_rna_layers_if_needed_xsp(obj_xsp)
  counts_xsp <- get_rna_counts_xsp(obj_xsp)
  meta_xsp <- obj_xsp@meta.data[colnames(counts_xsp), , drop = FALSE]
  
  group_key_xsp <- paste(as.character(meta_xsp$celltype), as.character(meta_xsp$orig.ident), sep = "___")
  group_fac_xsp <- factor(group_key_xsp, levels = sort(unique(group_key_xsp)))
  design_agg_xsp <- Matrix::sparseMatrix(
    i = seq_along(group_fac_xsp),
    j = as.integer(group_fac_xsp),
    x = 1,
    dims = c(length(group_fac_xsp), nlevels(group_fac_xsp)),
    dimnames = list(colnames(counts_xsp), levels(group_fac_xsp))
  )
  
  pb_counts_xsp <- counts_xsp %*% design_agg_xsp
  colnames(pb_counts_xsp) <- levels(group_fac_xsp)
  
  pb_meta_xsp <- data.frame(
    pb_col = colnames(pb_counts_xsp),
    stringsAsFactors = FALSE
  )
  split_xsp <- strsplit(pb_meta_xsp$pb_col, "___", fixed = TRUE)
  pb_meta_xsp$celltype <- vapply(split_xsp, function(x) x[1], character(1))
  pb_meta_xsp$sample <- vapply(split_xsp, function(x) x[2], character(1))
  pb_meta_xsp$sample_std <- standardize_mouse_sample_xsp(pb_meta_xsp$sample)
  
  sample_genotype_map_xsp <- meta_xsp |>
    dplyr::mutate(sample = as.character(orig.ident), genotype = as.character(genotype)) |>
    dplyr::distinct(sample, genotype)
  
  pb_meta_xsp <- pb_meta_xsp |>
    dplyr::left_join(sample_genotype_map_xsp, by = "sample") |>
    dplyr::mutate(
      genotype_short = infer_mouse_genotype_short_xsp(dplyr::if_else(!is.na(genotype), genotype, sample_std))
    )
  
  cell_count_df_xsp <- as.data.frame.matrix(table(meta_xsp$celltype, meta_xsp$orig.ident))
  cell_count_df_xsp <- tibble::rownames_to_column(cell_count_df_xsp, "celltype")
  cell_count_df_xsp$Total <- rowSums(cell_count_df_xsp[, setdiff(colnames(cell_count_df_xsp), "celltype"), drop = FALSE])
  write_tsv_xsp(cell_count_df_xsp, file.path(res_dir_xsp, "mouse_cell_counts_used_for_cross_species.tsv"))
  
  list(
    seurat = obj_xsp,
    pb_counts = pb_counts_xsp,
    pb_meta = pb_meta_xsp,
    cell_counts = cell_count_df_xsp
  )
}

mouse_rds_path_xsp <- first_existing_xsp(mouse_rds_candidates_xsp)
if (is.null(mouse_rds_path_xsp)) {
  stop("No mouse Seurat RDS was found. Expected one of:\n  ", paste(mouse_rds_candidates_xsp, collapse = "\n  "))
}

say_xsp("Loading mouse Seurat object: ", mouse_rds_path_xsp)
mouse_obj_xsp <- readRDS(mouse_rds_path_xsp)
if (!inherits(mouse_obj_xsp, "Seurat")) stop("The loaded mouse RDS is not a Seurat object.")

mouse_pb_xsp <- make_pseudobulk_counts_xsp(mouse_obj_xsp)
pb_counts_xsp <- mouse_pb_xsp$pb_counts
pb_meta_xsp <- mouse_pb_xsp$pb_meta

say_xsp("Mouse pseudobulk matrix: ", nrow(pb_counts_xsp), " genes x ", ncol(pb_counts_xsp), " celltype-sample columns")



## 4) Mouse cell-class pseudobulk differential expression and signature projection ------------

celltype_order_preferred_xsp <- c(
  "Myocyte/Myonucleus",
  "MHC-low macrophage",
  "MHC-high macrophage",
  "Dendritic cell(DC)",
  "Neutrophil",
  "Stromal",
  "Immune-Stromal(Cd53+)",
  "Endothelial cell(EC)",
  "Pericyte/Vascular SMC",
  "Tenocyte/Tendon fibroblast",
  "Satellite cell",
  "Schwann cell"
)

available_celltypes_xsp <- sort(unique(pb_meta_xsp$celltype))
celltypes_to_run_xsp <- c(
  intersect(celltype_order_preferred_xsp, available_celltypes_xsp),
  setdiff(available_celltypes_xsp, celltype_order_preferred_xsp)
)

run_mouse_celltype_de_xsp <- function(ct_xsp, pb_counts_xsp, pb_meta_xsp) {
  cols_xsp <- pb_meta_xsp |>
    dplyr::filter(celltype == ct_xsp, genotype_short %in% c("wt", "mdx", "d2")) |>
    dplyr::pull(pb_col)
  
  if (length(cols_xsp) == 0) return(NULL)
  meta_ct_xsp <- pb_meta_xsp[match(cols_xsp, pb_meta_xsp$pb_col), , drop = FALSE]
  grp_xsp <- factor(as.character(meta_ct_xsp$genotype_short), levels = c("wt", "mdx", "d2"))
  
  n_by_grp_xsp <- table(grp_xsp)
  if (any(n_by_grp_xsp[c("mdx", "d2")] < min_samples_per_contrast_group_xsp)) {
    say_xsp("Skipping ", ct_xsp, ": insufficient mdx or mdxD2 pseudobulk samples.")
    return(NULL)
  }
  
  mat_xsp <- as.matrix(pb_counts_xsp[, cols_xsp, drop = FALSE])
  y_xsp <- edgeR::DGEList(counts = mat_xsp)
  keep_gene_xsp <- edgeR::filterByExpr(y_xsp, group = grp_xsp)
  if (sum(keep_gene_xsp) < 50) {
    say_xsp("Skipping ", ct_xsp, ": fewer than 50 genes after filterByExpr.")
    return(NULL)
  }
  y_xsp <- y_xsp[keep_gene_xsp, , keep.lib.sizes = FALSE]
  y_xsp <- edgeR::calcNormFactors(y_xsp, method = "TMMwsp")
  
  design_xsp <- stats::model.matrix(~ 0 + grp_xsp)
  colnames(design_xsp) <- levels(grp_xsp)
  
  v_xsp <- limma::voomWithQualityWeights(y_xsp, design_xsp, plot = FALSE)
  fit0_xsp <- limma::lmFit(v_xsp, design_xsp)
  contr_xsp <- limma::makeContrasts(
    mdx_vs_wt = mdx - wt,
    d2_vs_mdx = d2 - mdx,
    d2_vs_wt = d2 - wt,
    levels = design_xsp
  )
  fit_xsp <- limma::contrasts.fit(fit0_xsp, contr_xsp)
  fit_xsp <- limma::eBayes(fit_xsp)
  
  top_xsp <- limma::topTable(fit_xsp, coef = "d2_vs_mdx", number = Inf, sort.by = "none")
  top_xsp$gene <- rownames(top_xsp)
  
  stats_xsp <- fit_xsp$t[, "d2_vs_mdx"]
  names(stats_xsp) <- rownames(fit_xsp$t)
  stats_xsp <- collapse_rank_by_gene_xsp(stats_xsp)
  
  list(
    celltype = ct_xsp,
    meta = meta_ct_xsp,
    group_counts = as.data.frame(n_by_grp_xsp),
    y = y_xsp,
    design = design_xsp,
    voom = v_xsp,
    contrast_vector = contr_xsp[, "d2_vs_mdx"],
    fit = fit_xsp,
    topTable = top_xsp,
    rank_stats = stats_xsp
  )
}

run_fgsea_projection_xsp <- function(de_obj_xsp, gene_sets_mouse_xsp) {
  stats_xsp <- de_obj_xsp$rank_stats
  pathways_xsp <- lapply(gene_sets_mouse_xsp, function(gs_xsp) intersect(gs_xsp, names(stats_xsp)))
  sz_xsp <- vapply(pathways_xsp, length, integer(1))
  pathways_xsp <- pathways_xsp[sz_xsp >= fgsea_min_size_xsp & sz_xsp <= fgsea_max_size_xsp]
  if (length(pathways_xsp) == 0) return(data.frame())
  
  fg_xsp <- tryCatch(
    fgsea::fgseaMultilevel(
      pathways = pathways_xsp,
      stats = stats_xsp,
      minSize = fgsea_min_size_xsp,
      maxSize = fgsea_max_size_xsp,
      nproc = 1
    ),
    error = function(e) {
      say_xsp("fgseaMultilevel failed for ", de_obj_xsp$celltype, "; falling back to fgseaSimple.")
      fgsea::fgseaSimple(
        pathways = pathways_xsp,
        stats = stats_xsp,
        nperm = 20000,
        minSize = fgsea_min_size_xsp,
        maxSize = fgsea_max_size_xsp,
        nproc = 1
      )
    }
  )
  
  fg_xsp <- data.table::as.data.table(fg_xsp)
  if (nrow(fg_xsp) == 0) return(data.frame())
  fg_xsp <- fg_xsp[order(padj, pval)]
  fg_xsp$celltype <- de_obj_xsp$celltype
  fg_xsp$leadingEdge <- vapply(fg_xsp$leadingEdge, paste, character(1), collapse = ";")
  as.data.frame(fg_xsp)
}

run_camera_projection_xsp <- function(de_obj_xsp, gene_sets_mouse_xsp) {
  genes_in_xsp <- rownames(de_obj_xsp$voom$E)
  idx_xsp <- lapply(gene_sets_mouse_xsp, function(gs_xsp) {
    ix_xsp <- match(intersect(gs_xsp, genes_in_xsp), genes_in_xsp)
    ix_xsp[!is.na(ix_xsp)]
  })
  idx_xsp <- idx_xsp[vapply(idx_xsp, length, integer(1)) >= camera_min_size_xsp]
  if (length(idx_xsp) == 0) return(data.frame())
  
  cam_xsp <- limma::camera(
    de_obj_xsp$voom,
    index = idx_xsp,
    design = de_obj_xsp$design,
    contrast = de_obj_xsp$contrast_vector
  )
  cam_xsp$FDR_local <- p.adjust(cam_xsp$PValue, method = "BH")
  cam_xsp$signature <- rownames(cam_xsp)
  cam_xsp$celltype <- de_obj_xsp$celltype
  cam_xsp <- cam_xsp[order(cam_xsp$FDR_local, cam_xsp$PValue), , drop = FALSE]
  rownames(cam_xsp) <- NULL
  as.data.frame(cam_xsp)
}

say_xsp("Running cell-class pseudobulk DE and cross-species signature projection.")

de_objects_xsp <- list()
fgsea_all_xsp <- list()
camera_all_xsp <- list()
top_tables_xsp <- list()

for (ct_xsp in celltypes_to_run_xsp) {
  say_xsp("Cell class: ", ct_xsp)
  de_ct_xsp <- run_mouse_celltype_de_xsp(ct_xsp, pb_counts_xsp, pb_meta_xsp)
  if (is.null(de_ct_xsp)) next
  de_objects_xsp[[ct_xsp]] <- de_ct_xsp
  
  fgsea_ct_xsp <- run_fgsea_projection_xsp(de_ct_xsp, gene_sets_mouse_xsp)
  if (nrow(fgsea_ct_xsp) > 0) fgsea_all_xsp[[ct_xsp]] <- fgsea_ct_xsp
  
  camera_ct_xsp <- run_camera_projection_xsp(de_ct_xsp, gene_sets_mouse_xsp)
  if (nrow(camera_ct_xsp) > 0) camera_all_xsp[[ct_xsp]] <- camera_ct_xsp
  
  tt_xsp <- de_ct_xsp$topTable |>
    dplyr::mutate(celltype = ct_xsp) |>
    dplyr::select(celltype, gene, dplyr::everything())
  top_tables_xsp[[ct_xsp]] <- tt_xsp
}

if (length(de_objects_xsp) == 0) {
  stop("No cell class could be analysed. Check pseudobulk sample labels and genotype inference.")
}

fgsea_df_xsp <- dplyr::bind_rows(fgsea_all_xsp)
camera_df_xsp <- dplyr::bind_rows(camera_all_xsp)
mouse_de_gene_tables_xsp <- dplyr::bind_rows(top_tables_xsp)

if (nrow(fgsea_df_xsp) == 0) stop("fgsea returned no results for the human signatures.")
if (nrow(camera_df_xsp) == 0) say_xsp("camera returned no results; fgsea results will still be exported.")

fgsea_df_xsp <- fgsea_df_xsp |>
  dplyr::rename(signature = pathway) |>
  dplyr::mutate(
    padj_global = p.adjust(pval, method = "BH"),
    signature_label = pretty_signature_xsp(signature),
    expected_sign = dplyr::case_when(
      signature == "Human_Late_up_DEG_orthologs" ~ 1,
      signature == "Human_Early_up_DEG_orthologs" ~ -1,
      TRUE ~ NA_real_
    ),
    concordant_direction = sign(NES) == expected_sign,
    signed_neglog10_padj_global = sign(NES) * safe_neglog10_xsp(padj_global),
    signed_neglog10_padj_local = sign(NES) * safe_neglog10_xsp(padj),
    contrast = contrast_name_xsp
  ) |>
  dplyr::select(
    celltype, contrast, signature, signature_label,
    size, ES, NES, pval, padj, padj_global,
    expected_sign, concordant_direction,
    signed_neglog10_padj_local, signed_neglog10_padj_global,
    leadingEdge, dplyr::everything()
  )

if (nrow(camera_df_xsp) > 0) {
  camera_df_xsp <- camera_df_xsp |>
    dplyr::mutate(
      FDR_global = p.adjust(PValue, method = "BH"),
      signature_label = pretty_signature_xsp(signature),
      expected_direction = dplyr::case_when(
        signature == "Human_Late_up_DEG_orthologs" ~ "Up",
        signature == "Human_Early_up_DEG_orthologs" ~ "Down",
        TRUE ~ NA_character_
      ),
      concordant_direction = Direction == expected_direction,
      signed_neglog10_FDR_global = dplyr::case_when(
        Direction == "Up" ~ safe_neglog10_xsp(FDR_global),
        Direction == "Down" ~ -safe_neglog10_xsp(FDR_global),
        TRUE ~ NA_real_
      ),
      contrast = contrast_name_xsp
    ) |>
    dplyr::select(
      celltype, contrast, signature, signature_label,
      NGenes, Direction, PValue, FDR_local, FDR_global,
      expected_direction, concordant_direction,
      signed_neglog10_FDR_global, dplyr::everything()
    )
}

# Cell-class summary score:
#   Positive score means Late-up human orthologs enrich toward mdxD2-up genes and
#   Early-up human orthologs enrich toward mdxD2-down genes.
score_wide_xsp <- fgsea_df_xsp |>
  dplyr::select(celltype, signature, signed_neglog10_padj_global, NES, padj_global, concordant_direction) |>
  tidyr::pivot_wider(
    names_from = signature,
    values_from = c(signed_neglog10_padj_global, NES, padj_global, concordant_direction),
    names_sep = "__"
  )

late_signed_col_xsp <- "signed_neglog10_padj_global__Human_Late_up_DEG_orthologs"
early_signed_col_xsp <- "signed_neglog10_padj_global__Human_Early_up_DEG_orthologs"

summary_score_xsp <- score_wide_xsp |>
  dplyr::mutate(
    phase_concordance_score = .data[[late_signed_col_xsp]] - .data[[early_signed_col_xsp]],
    late_component = .data[[late_signed_col_xsp]],
    early_component = -.data[[early_signed_col_xsp]]
  ) |>
  dplyr::arrange(dplyr::desc(phase_concordance_score))

celltype_plot_order_xsp <- summary_score_xsp$celltype
fgsea_df_xsp$celltype_ordered <- factor(fgsea_df_xsp$celltype, levels = rev(celltype_plot_order_xsp))
if (nrow(camera_df_xsp) > 0) {
  camera_df_xsp$celltype_ordered <- factor(camera_df_xsp$celltype, levels = rev(celltype_plot_order_xsp))
}



## 5) Export tables -------------------------------------------------------------

say_xsp("Writing cross-species concordance tables.")

human_signature_summary_xsp <- human_to_mouse_sig_xsp |>
  dplyr::select(Signature, Human_direction, EntrezID, Human_SYMBOL, Mouse_SYMBOL, logFC, FDR, Pvalue, Ortholog_support, Source) |>
  dplyr::arrange(Signature, Human_SYMBOL, Mouse_SYMBOL)

signature_size_summary_xsp <- human_signature_summary_xsp |>
  dplyr::group_by(Signature, Human_direction) |>
  dplyr::summarise(
    human_deg_symbols_n = dplyr::n_distinct(Human_SYMBOL),
    mouse_ortholog_symbols_n = dplyr::n_distinct(Mouse_SYMBOL),
    .groups = "drop"
  )

write_tsv_xsp(signature_size_summary_xsp, file.path(res_dir_xsp, "human_signature_size_summary.tsv"))
write_tsv_xsp(human_signature_summary_xsp, file.path(res_dir_xsp, "human_DEG_to_mouse_ortholog_signature_table.tsv"))
write_tsv_xsp(fgsea_df_xsp, file.path(res_dir_xsp, "cross_species_fgsea_by_cellclass.tsv"))
if (nrow(camera_df_xsp) > 0) write_tsv_xsp(camera_df_xsp, file.path(res_dir_xsp, "cross_species_camera_by_cellclass.tsv"))
write_tsv_xsp(summary_score_xsp, file.path(res_dir_xsp, "cross_species_cellclass_concordance_summary.tsv"))
write_tsv_xsp(mouse_de_gene_tables_xsp, file.path(res_dir_xsp, "mouse_cellclass_limma_voom_d2_vs_mdx_gene_tables.tsv"))
write_tsv_xsp(pb_meta_xsp, file.path(res_dir_xsp, "mouse_pseudobulk_column_metadata.tsv"))

# Compact leading-edge table.
leading_edge_xsp <- fgsea_df_xsp |>
  dplyr::select(celltype, signature, signature_label, NES, padj_global, leadingEdge) |>
  tidyr::separate_rows(leadingEdge, sep = ";") |>
  dplyr::rename(leading_edge_gene = leadingEdge) |>
  dplyr::filter(!is.na(leading_edge_gene), leading_edge_gene != "")
write_tsv_xsp(leading_edge_xsp, file.path(res_dir_xsp, "cross_species_fgsea_leading_edge_genes.tsv"))

# XLSX workbook for manuscript supplement. This is optional; TSV files above are
# the primary machine-readable outputs and are always written.
xlsx_output_path_xsp <- file.path(res_dir_xsp, "Additional_file_5_cross_species_concordance_outputs.xlsx")
if (isTRUE(write_xlsx_if_openxlsx_available_xsp) && requireNamespace("openxlsx", quietly = TRUE)) {
  wb_xsp <- openxlsx::createWorkbook()
  add_sheet_xsp <- function(wb_xsp, sheet_name_xsp, df_xsp) {
    sheet_name_xsp <- substr(sheet_name_xsp, 1, 31)
    openxlsx::addWorksheet(wb_xsp, sheet_name_xsp)
    openxlsx::writeData(wb_xsp, sheet_name_xsp, df_xsp)
    openxlsx::freezePane(wb_xsp, sheet_name_xsp, firstActiveRow = 2)
    openxlsx::setColWidths(wb_xsp, sheet_name_xsp, cols = 1:ncol(df_xsp), widths = "auto")
  }
  
  add_sheet_xsp(wb_xsp, "Signature_sizes", signature_size_summary_xsp)
  add_sheet_xsp(wb_xsp, "Human_to_mouse_signature", human_signature_summary_xsp)
  add_sheet_xsp(wb_xsp, "fgsea_by_cellclass", fgsea_df_xsp)
  if (nrow(camera_df_xsp) > 0) add_sheet_xsp(wb_xsp, "camera_by_cellclass", camera_df_xsp)
  add_sheet_xsp(wb_xsp, "Concordance_summary", summary_score_xsp)
  add_sheet_xsp(wb_xsp, "Leading_edges", leading_edge_xsp)
  add_sheet_xsp(wb_xsp, "Mouse_pb_metadata", pb_meta_xsp)
  add_sheet_xsp(wb_xsp, "Mouse_cell_counts", mouse_pb_xsp$cell_counts)
  openxlsx::saveWorkbook(wb_xsp, xlsx_output_path_xsp, overwrite = TRUE)
  say_xsp("Optional XLSX workbook written: ", xlsx_output_path_xsp)
} else {
  say_xsp("openxlsx is not installed or XLSX output is disabled; skipping optional XLSX workbook.")
}

# Save reusable R object.
saveRDS(
  list(
    settings = list(
      proj_root_xsp = proj_root_xsp,
      human_stage6_path_xsp = human_stage6_path_xsp,
      mouse_rds_path_xsp = mouse_rds_path_xsp,
      human_deg_fdr_cut_xsp = human_deg_fdr_cut_xsp,
      human_abs_logfc_min_xsp = human_abs_logfc_min_xsp,
      ortholog_min_support_xsp = ortholog_min_support_xsp,
      fgsea_min_size_xsp = fgsea_min_size_xsp,
      fgsea_max_size_xsp = fgsea_max_size_xsp,
      camera_min_size_xsp = camera_min_size_xsp
    ),
    human_deg = human_deg_xsp,
    ortholog_table = ortholog_tbl_xsp,
    human_signature_summary = human_signature_summary_xsp,
    gene_sets_mouse = gene_sets_mouse_xsp,
    pb_meta = pb_meta_xsp,
    fgsea = fgsea_df_xsp,
    camera = camera_df_xsp,
    summary = summary_score_xsp,
    sessionInfo = utils::sessionInfo()
  ),
  file = file.path(res_dir_xsp, "cross_species_concordance_analysis_objects.rds"),
  compress = "xz",
  version = 3
)



## 6) Figures -------------------------------------------------------------------

say_xsp("Generating figures.")

fgsea_plot_df_xsp <- fgsea_df_xsp |>
  dplyr::mutate(
    neglog10_FDR = pmin(safe_neglog10_xsp(padj_global), 25),
    significance_label = dplyr::if_else(padj_global < 0.05, "FDR < 0.05", "NS"),
    signature_label = factor(
      signature_label,
      levels = c("Human Early-up DEG orthologs", "Human Late-up DEG orthologs")
    )
  )

p_dot_xsp <- ggplot2::ggplot(
  fgsea_plot_df_xsp,
  ggplot2::aes(x = signature_label, y = celltype_ordered)
) +
  ggplot2::geom_point(
    ggplot2::aes(size = neglog10_FDR, colour = NES, shape = significance_label),
    alpha = 0.95
  ) +
  ggplot2::scale_colour_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "fgsea NES"
  ) +
  ggplot2::scale_size_continuous(
    range = c(2, 8),
    name = expression(-log[10]~global~FDR)
  ) +
  ggplot2::scale_shape_manual(values = c("FDR < 0.05" = 16, "NS" = 1), name = NULL) +
  ggplot2::labs(
    title = "Cross-species concordance: human age-phase DEG signatures in mouse cell-class pseudobulk ranks",
    subtitle = "Mouse ranking metric: limma-voom t-statistic for mdxD2 versus mdx within each cell class",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )

summary_plot_df_xsp <- summary_score_xsp |>
  dplyr::mutate(celltype_ordered = factor(celltype, levels = rev(celltype_plot_order_xsp)))

p_score_xsp <- ggplot2::ggplot(
  summary_plot_df_xsp,
  ggplot2::aes(x = phase_concordance_score, y = celltype_ordered)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::labs(
    title = "Net phase-concordance score by mouse cell class",
    subtitle = "Score = signed Late-up enrichment minus signed Early-up enrichment; higher values indicate stronger human Late/Early concordance",
    x = "Net phase-concordance score",
    y = NULL
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

fgsea_wide_for_scatter_xsp <- fgsea_df_xsp |>
  dplyr::select(celltype, signature, NES, padj_global) |>
  tidyr::pivot_wider(names_from = signature, values_from = c(NES, padj_global), names_sep = "__")

p_scatter_xsp <- ggplot2::ggplot(
  fgsea_wide_for_scatter_xsp,
  ggplot2::aes(
    x = NES__Human_Early_up_DEG_orthologs,
    y = NES__Human_Late_up_DEG_orthologs,
    label = celltype
  )
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(check_overlap = TRUE, nudge_y = 0.05, size = 3) +
  ggplot2::labs(
    title = "Human Early-up versus Late-up signature NES by mouse cell class",
    x = "NES for human Early-up DEG orthologs",
    y = "NES for human Late-up DEG orthologs"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

fig_pdf_xsp <- file.path(fig_dir_xsp, "cross_species_concordance_fgsea_dotplot.pdf")
fig_png_xsp <- file.path(fig_dir_xsp, "cross_species_concordance_fgsea_dotplot.png")
score_pdf_xsp <- file.path(fig_dir_xsp, "cross_species_concordance_net_score.pdf")
score_png_xsp <- file.path(fig_dir_xsp, "cross_species_concordance_net_score.png")
scatter_pdf_xsp <- file.path(fig_dir_xsp, "cross_species_concordance_NES_scatter.pdf")
scatter_png_xsp <- file.path(fig_dir_xsp, "cross_species_concordance_NES_scatter.png")

# Sizes chosen for supplementary figure readability.
ggplot2::ggsave(fig_pdf_xsp, p_dot_xsp, width = 8.5, height = 6.8, units = "in")
ggplot2::ggsave(fig_png_xsp, p_dot_xsp, width = 8.5, height = 6.8, units = "in", dpi = 600)
ggplot2::ggsave(score_pdf_xsp, p_score_xsp, width = 8.0, height = 6.0, units = "in")
ggplot2::ggsave(score_png_xsp, p_score_xsp, width = 8.0, height = 6.0, units = "in", dpi = 600)
ggplot2::ggsave(scatter_pdf_xsp, p_scatter_xsp, width = 7.0, height = 5.5, units = "in")
ggplot2::ggsave(scatter_png_xsp, p_scatter_xsp, width = 7.0, height = 5.5, units = "in", dpi = 600)

# Optional camera dot plot.
if (nrow(camera_df_xsp) > 0) {
  camera_plot_df_xsp <- camera_df_xsp |>
    dplyr::mutate(
      neglog10_FDR = pmin(safe_neglog10_xsp(FDR_global), 25),
      Direction_numeric = dplyr::case_when(Direction == "Up" ~ 1, Direction == "Down" ~ -1, TRUE ~ 0),
      signed_FDR_metric = Direction_numeric * safe_neglog10_xsp(FDR_global),
      significance_label = dplyr::if_else(FDR_global < 0.05, "FDR < 0.05", "NS"),
      signature_label = factor(
        signature_label,
        levels = c("Human Early-up DEG orthologs", "Human Late-up DEG orthologs")
      )
    )
  
  p_camera_xsp <- ggplot2::ggplot(
    camera_plot_df_xsp,
    ggplot2::aes(x = signature_label, y = celltype_ordered)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = neglog10_FDR, colour = signed_FDR_metric, shape = significance_label),
      alpha = 0.95
    ) +
    ggplot2::scale_colour_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = "Direction × -log10(FDR)"
    ) +
    ggplot2::scale_size_continuous(range = c(2, 8), name = expression(-log[10]~global~FDR)) +
    ggplot2::scale_shape_manual(values = c("FDR < 0.05" = 16, "NS" = 1), name = NULL) +
    ggplot2::labs(
      title = "Cross-species concordance: camera competitive gene-set test",
      subtitle = "Mouse contrast: mdxD2 versus mdx within each cell class",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  ggplot2::ggsave(file.path(fig_dir_xsp, "cross_species_concordance_camera_dotplot.pdf"), p_camera_xsp, width = 8.5, height = 6.8, units = "in")
  ggplot2::ggsave(file.path(fig_dir_xsp, "cross_species_concordance_camera_dotplot.png"), p_camera_xsp, width = 8.5, height = 6.8, units = "in", dpi = 600)
}



## 7) Session information and console summary -----------------------------------

session_txt_xsp <- capture.output(utils::sessionInfo())
writeLines(session_txt_xsp, con = file.path(res_dir_xsp, "sessionInfo_cross_species_concordance.txt"))

say_xsp("Done. Main outputs:")
if (file.exists(xlsx_output_path_xsp)) say_xsp("  ", xlsx_output_path_xsp)
say_xsp("  ", file.path(res_dir_xsp, "cross_species_fgsea_by_cellclass.tsv"))
say_xsp("  ", file.path(res_dir_xsp, "cross_species_cellclass_concordance_summary.tsv"))
say_xsp("  ", file.path(fig_dir_xsp, "cross_species_concordance_fgsea_dotplot.pdf"))
say_xsp("  ", file.path(fig_dir_xsp, "cross_species_concordance_net_score.pdf"))

print(signature_size_summary_xsp)
print(summary_score_xsp)


## 8) Option -----------------------------------------------------------------

##──────────────────────────────────────────────────────────────────────────────
## Purpose:
##   Display only the Myocyte/Myonucleus result from the exploratory
##   human-to-mouse ortholog projection.
##
## Required input files:
##   C:/DMD_project/results_cross_species_concordance/
##     cross_species_fgsea_by_cellclass.tsv
##     human_signature_size_summary.tsv
##
## This script does NOT rerun the analysis.
## It only reads the exported TSV files and prints a ggplot object.
##──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})



## ---- 0) Paths ---------------------------------------------------------------------

res_dir_fig5a <- "C:/DMD_project/results_cross_species_concordance"
fig_dir_fig5a <- file.path(res_dir_fig5a, "figures")
dir.create(fig_dir_fig5a, showWarnings = FALSE, recursive = TRUE)

fgsea_file_fig5a <- file.path(res_dir_fig5a, "cross_species_fgsea_by_cellclass.tsv")
size_file_fig5a  <- file.path(res_dir_fig5a, "human_signature_size_summary.tsv")

if (!file.exists(fgsea_file_fig5a)) {
  stop("Missing file: ", fgsea_file_fig5a)
}
if (!file.exists(size_file_fig5a)) {
  stop("Missing file: ", size_file_fig5a)
}



## ---- 1) Read exported results ----------------------------------------------------

fgsea_fig5a <- read.delim(
  fgsea_file_fig5a,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sigsize_fig5a <- read.delim(
  size_file_fig5a,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)



## ---- 2) Print significant hits from the full cell-class screen --------------------
## This is useful for checking that the main panel is not cherry-picked.

sig_hits_fig5a <- fgsea_fig5a %>%
  filter(padj_global < 0.05) %>%
  arrange(padj_global) %>%
  select(
    celltype,
    signature_label,
    size,
    NES,
    pval,
    padj,
    padj_global,
    concordant_direction
  )

cat("\n[Global FDR-significant hits from the full cell-class screen]\n")
if (nrow(sig_hits_fig5a) == 0) {
  cat("No globally significant hit was detected.\n")
} else {
  print(sig_hits_fig5a, row.names = FALSE, digits = 3)
}



## ---- 3) Prepare Myocyte/Myonucleus-only plotting table ----------------------------

myo_fig5a <- fgsea_fig5a %>%
  filter(celltype == "Myocyte/Myonucleus") %>%
  mutate(
    signature_short = case_when(
      signature == "Human_Early_up_DEG_orthologs" ~ "Human Early-up DEG orthologs",
      signature == "Human_Late_up_DEG_orthologs"  ~ "Human Late-up DEG orthologs",
      TRUE ~ signature_label
    ),
    signature_short = factor(
      signature_short,
      levels = c("Human Early-up DEG orthologs",
                 "Human Late-up DEG orthologs")
    ),
    biological_readout = case_when(
      signature == "Human_Early_up_DEG_orthologs" ~ "Early programme lost in mdxD2",
      signature == "Human_Late_up_DEG_orthologs" ~ "Late programme directional gain",
      TRUE ~ ""
    ),
    significance_label = ifelse(
      padj_global < 0.05,
      "Global FDR < 0.05",
      "Not significant"
    ),
    fdr_text = ifelse(
      padj_global < 0.001,
      "global FDR < 0.001",
      paste0("global FDR = ", formatC(padj_global, format = "fg", digits = 2))
    ),
    nes_text = paste0("NES = ", sprintf("%.2f", NES)),
    n_text = paste0("genes in rank n = ", size),
    label_text = paste0(
      biological_readout, "\n",
      nes_text, "\n",
      fdr_text, "\n",
      n_text
    )
  )

if (nrow(myo_fig5a) != 2) {
  stop("Expected two Myocyte/Myonucleus rows: Early-up and Late-up signatures.")
}

## label position for horizontal bars
x_pad_fig5a <- 0.22
myo_fig5a <- myo_fig5a %>%
  mutate(
    label_x = ifelse(NES >= 0, NES + x_pad_fig5a, NES - x_pad_fig5a),
    label_hjust = ifelse(NES >= 0, 0, 1)
  )

## Add total ortholog-set sizes for the caption-like note.
size_note_fig5a <- sigsize_fig5a %>%
  mutate(
    label = case_when(
      Signature == "Human_Early_up_DEG_orthologs" ~
        paste0("Early-up set: ", mouse_ortholog_symbols_n, " mouse orthologs"),
      Signature == "Human_Late_up_DEG_orthologs" ~
        paste0("Late-up set: ", mouse_ortholog_symbols_n, " mouse orthologs"),
      TRUE ~ paste0(Signature, ": ", mouse_ortholog_symbols_n)
    )
  ) %>%
  pull(label) %>%
  paste(collapse = "; ")

x_lim_fig5a <- max(abs(myo_fig5a$NES), na.rm = TRUE) * 1.65

## ---- 4) Plot ----------------------------------------------------------------------
## Horizontal version for a wide Figure 5A panel

p_fig5a <- ggplot(
  myo_fig5a,
  aes(y = signature_short, x = NES)
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    linetype = "dashed",
    colour = "grey45"
  ) +
  geom_col(
    aes(fill = signature),
    width = 0.55,
    colour = "grey25",
    linewidth = 0.35
  ) +
  geom_point(
    aes(shape = significance_label),
    size = 3.0,
    colour = "black"
  ) +
  geom_text(
    aes(x = label_x, label = label_text, hjust = label_hjust),
    size = 3.3,
    lineheight = 1.02
  ) +
  scale_fill_manual(
    values = c(
      "Human_Early_up_DEG_orthologs" = "#2C7BB6",
      "Human_Late_up_DEG_orthologs"  = "#D7191C"
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("Global FDR < 0.05" = 16, "Not significant" = 1),
    name = NULL
  ) +
  coord_cartesian(
    xlim = c(-x_lim_fig5a, x_lim_fig5a),
    clip = "off"
  ) +
  labs(
    title = "Myocyte/Myonucleus-resolved ortholog projection",
    subtitle = "Human bulk age-phase DEG signatures projected onto the mouse mdxD2-versus-mdx Myocyte/Myonucleus ranking",
    x = "fgsea normalized enrichment score",
    y = NULL,
    caption = paste0(
      "Negative NES indicates enrichment among mdxD2-down genes; positive NES indicates enrichment among mdxD2-up genes.\n",
      "Global FDR was adjusted across all tested cell class-signature combinations. ",
      size_note_fig5a
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    plot.caption = element_text(size = 8, hjust = 0),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 28, b = 10, l = 10)
  )

## ---- 5) Display in the R plotting window ------------------------------------------

print(p_fig5a)
                 



### Session information---------------------------------------------------------

sessionInfo()
# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)

# Matrix products: default
# LAPACK version 3.12.1

# locale:
#   [1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8    LC_MONETARY=Japanese_Japan.utf8 LC_NUMERIC=C                   
# [5] LC_TIME=Japanese_Japan.utf8    

# time zone: Etc/GMT-9
# tzcode source: internal

# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] cli_3.6.5                   nichenetr_2.2.0             edgeR_4.6.3                 scDblFinder_1.22.0          SingleCellExperiment_1.30.1
# [6] SummarizedExperiment_1.38.1 GenomicRanges_1.60.0        GenomeInfoDb_1.44.1         MatrixGenerics_1.20.0       Matrix_1.7-3               
# [11] Seurat_5.3.0                SeuratObject_5.1.0          sp_2.2-0                    matrixStats_1.5.0           hgu133b.db_3.13.0          
# [16] hgu133a.db_3.13.0           hgu133plus2.db_3.13.0       ggraph_2.2.1                igraph_2.1.4                pheatmap_1.0.13            
# [21] STRINGdb_2.20.0             ReactomePA_1.52.0           clusterProfiler_4.16.0      org.Hs.eg.db_3.21.0         AnnotationDbi_1.70.0       
# [26] IRanges_2.42.0              S4Vectors_0.46.0            Biobase_2.68.0              BiocGenerics_0.54.0         generics_0.1.4             
# [31] impute_1.82.0               sva_3.56.0                  BiocParallel_1.42.1         genefilter_1.90.0           mgcv_1.9-3                 
# [36] nlme_3.1-168                limma_3.64.3                lubridate_1.9.4             forcats_1.0.0               stringr_1.5.2              
# [41] dplyr_1.1.4                 purrr_1.1.0                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.3.0               
# [46] ggplot2_4.0.0               tidyverse_2.0.0            

# loaded via a namespace (and not attached):
#   [1] hash_2.2.6.3             graph_1.86.0             ica_1.0-3                plotly_4.12.0            Formula_1.2-5           
# [6] scater_1.36.0            MBESS_4.9.41             tidyselect_1.2.1         bit_4.6.0                doParallel_1.0.17       
# [11] BWStest_0.2.3            clue_0.3-66              lattice_0.22-7           rjson_0.2.23             bridgesampling_1.2-1    
# [16] blob_1.3.0               S4Arrays_1.8.1           parallel_4.5.1           caret_7.0-1              png_0.1-8               
# [21] plotrix_3.8-4            ggplotify_0.1.2          goftest_1.2-3            BiocIO_1.18.0            bluster_1.18.0          
# [26] BiocNeighbors_2.2.0      uwot_0.2.3               shadowtext_0.1.6         curl_7.0.0               mime_0.13               
# [31] evaluate_1.0.5           tidytree_0.4.6           gsubfn_0.7               ComplexHeatmap_2.24.1    stringi_1.8.7           
# [36] pROC_1.19.0.1            backports_1.5.0          PMCMRplus_1.9.12         XML_3.99-0.18            httpuv_1.6.16           
# [41] magrittr_2.0.3           rappdirs_0.3.3           splines_4.5.1            prodlim_2025.04.28       ggbeeswarm_0.7.2        
# [46] sctransform_0.4.2        effsize_0.8.1            DBI_1.3.0                reactome.db_1.92.0       withr_3.0.2             
# [51] reformulas_0.4.1         class_7.3-23             systemfonts_1.3.1        xgboost_1.7.11.1         enrichplot_1.28.4       
# [56] lmtest_0.9-40            ggnewscale_0.5.2         tidygraph_1.3.1          rtracklayer_1.68.0       BiocManager_1.30.27     
# [61] htmlwidgets_1.6.4        fs_1.6.6                 SuppDists_1.1-9.9        ggrepel_0.9.6            labeling_0.4.3          
# [66] SparseArray_1.8.1        annotate_1.86.1          reticulate_1.43.0        zoo_1.8-14               XVector_0.48.0          
# [71] knitr_1.51               UCSC.utils_1.4.0         timechange_0.3.0         foreach_1.5.2            patchwork_1.3.2         
# [76] caTools_1.18.3           visNetwork_2.1.4         grid_4.5.1               data.table_1.17.8        timeDate_4041.110       
# [81] ggtree_3.16.3            R.oo_1.27.1              ggiraph_0.9.1            RSpectra_0.16-2          irlba_2.3.5.1           
# [86] DiagrammeR_1.0.11        fastDummies_1.7.5        gridGraphics_0.5-1       yaml_2.3.10              lazyeval_0.2.2          
# [91] conflicted_1.2.0         survival_3.8-3           scattermore_1.2          crayon_1.5.3             tensorA_0.36.2.1        
# [96] RcppAnnoy_0.0.22         RColorBrewer_1.1-3       progressr_0.18.0         tweenr_2.0.3             later_1.4.2             
# [101] ggridges_0.5.7           codetools_0.2-20         base64enc_0.1-3          GlobalOptions_0.1.3      KEGGREST_1.48.1         
# [106] Rtsne_0.17               shape_1.4.6.1            estimability_1.5.1       Rsamtools_2.24.0         sqldf_0.4-11            
# [111] foreign_0.8-90           pkgconfig_2.0.3          spatstat.univar_3.1-4    ggpubr_0.6.3             GenomicAlignments_1.44.0
# [116] aplot_0.2.8              spatstat.sparse_3.1-0    ape_5.8-1                viridisLite_0.4.2        xtable_1.8-4            
# [121] car_3.1-3                plyr_1.8.9               httr_1.4.8               rbibutils_2.3            tools_4.5.1             
# [126] globals_0.19.1           brms_2.23.0              hardhat_1.4.2            bayesplot_1.15.0         beeswarm_0.4.0          
# [131] htmlTable_2.4.3          broom_1.0.12             checkmate_2.3.3          loo_2.9.0                lme4_1.1-37             
# [136] assertthat_0.2.1         digest_0.6.37            farver_2.1.2             tzdb_0.5.0               reshape2_1.4.4          
# [141] ModelMetrics_1.2.2.2     yulab.utils_0.2.4        viridis_0.6.5            rpart_4.1.24             glue_1.8.0              
# [146] cachem_1.1.0             polyclip_1.10-7          Hmisc_5.2-3              Biostrings_2.76.0        mvtnorm_1.3-3           
# [151] presto_1.0.0             proto_1.0.0              parallelly_1.45.1        statmod_1.5.0            RcppHNSW_0.6.0          
# [156] ScaledMatrix_1.16.0      carData_3.0-5            minqa_1.2.8              pbapply_1.7-4            pwr_1.3-0               
# [161] spam_2.11-1              gson_0.1.0               dqrng_0.4.1              utf8_1.2.6               gower_1.0.2             
# [166] graphlayouts_1.2.2       gtools_3.9.5             ggsignif_0.6.4           gridExtra_2.3            shiny_1.11.1            
# [171] lava_1.8.2               GenomeInfoDbData_1.2.14  R.utils_2.13.0           RCurl_1.98-1.17          memoise_2.0.1           
# [176] rmarkdown_2.30           scales_1.4.0             R.methodsS3_1.8.2        future_1.70.0            RANN_2.6.2              
# [181] spatstat.data_3.1-8      rstudioapi_0.18.0        cluster_2.1.8.1          msigdbr_26.1.0           rstantools_2.6.0        
# [186] spatstat.utils_3.1-5     hms_1.1.4                fitdistrplus_1.2-6       fdrtool_1.2.18           cowplot_1.2.0           
# [191] colorspace_2.1-2         rlang_1.1.6              ipred_0.9-15             dotCall64_1.2            scuttle_1.18.0          
# [196] ggforce_0.5.0            circlize_0.4.17          ggtangle_0.0.7           xfun_0.53                multcompView_0.1-11     
# [201] coda_0.19-4.1            e1071_1.7-16             TH.data_1.1-4            posterior_1.6.1          recipes_1.3.1           
# [206] iterators_1.0.14         emmeans_2.0.2            abind_1.4-8              randomForest_4.7-1.2     GOSemSim_2.34.0         
# [211] treeio_1.32.0            gmp_0.7-5                Rdpack_2.6.4             bitops_1.0-9             promises_1.3.3          
# [216] RSQLite_2.4.3            qvalue_2.40.0            sandwich_3.1-1           fgsea_1.34.2             DelayedArray_0.34.1     
# [221] proxy_0.4-27             Rmpfr_1.1-1              GO.db_3.21.0             compiler_4.5.1           beachmat_2.24.0         
# [226] boot_1.3-32              distributional_0.7.0     graphite_1.54.0          listenv_0.10.1           Rcpp_1.1.0              
# [231] BiocSingular_1.24.0      tensor_1.5.1             MASS_7.3-65              kSamples_1.2-12          uuid_1.2-1              
# [236] babelgene_22.9           spatstat.random_3.4-1    R6_2.6.1                 fastmap_1.2.0            multcomp_1.4-28         
# [241] fastmatch_1.1-6          rstatix_0.7.3            vipor_0.4.7              ROCR_1.0-12              rsvd_1.0.5              
# [246] nnet_7.3-20              gtable_0.3.6             KernSmooth_2.23-26       miniUI_0.1.2             deldir_2.0-4            
# [251] htmltools_0.5.8.1        RcppParallel_5.1.11-1    bit64_4.6.0-1            spatstat.explore_3.5-2   lifecycle_1.0.5         
# [256] S7_0.2.0                 Brobdingnag_1.2-9        restfulr_0.0.16          nloptr_2.2.1             vctrs_0.6.5             
# [261] spatstat.geom_3.5-0      DOSE_4.2.0               scran_1.36.0             ggfun_0.2.0              future.apply_1.20.2     
# [266] pillar_1.11.1            gplots_3.2.0             metapod_1.16.0           locfit_1.5-9.12          otel_0.2.0              
# [271] jsonlite_2.0.0           chron_2.3-62             GetoptLong_1.1.0 



## Save the entire current workspace-------------------------------------------
getwd()         # example: "C:/DMD_project"
# Specify the save directory (e.g., /cache under the project root)
#save_dir <- "C:/DMD_project/cache"
#dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
#save.image(file = file.path(save_dir, "workspace_2026-04-24.RData"))
#savehistory(file = file.path(save_dir, "workspace_2026-04-24.Rhistory"))

#load("cache/workspace_2026-04-24.RData")
