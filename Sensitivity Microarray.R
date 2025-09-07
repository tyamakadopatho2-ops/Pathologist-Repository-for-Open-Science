##### E: ------------------------------------------------Sensitivity Analysis-------------------------------------------------------
# 年齢カットオフを 4/6/7歳で置換し、DEG 数推移と GSEA NES の相関を表で提示
###############################################################################
## 0.  年齢 6 歳カットオフのグループ因子
###############################################################################
Group_sense <- factor(ifelse(AgeL[PID] < 6, "Early", "Late"),
                      levels = c("Early","Late"))

###############################################################################
## 1.  ComBat を “6歳版” で再実行           （expr_qs は再利用）
###############################################################################
library(sva)

combat_sense <- ComBat(expr_qs,
                       batch     = Cohort,
                       mod       = model.matrix(~ Group_sense),
                       ref.batch = "C1",
                       par.prior = TRUE,
                       mean.only = FALSE)

## さらに limma::removeBatchEffect で微調整
library(limma)
combat_sense <- removeBatchEffect(combat_sense,
                                  batch  = Cohort,
                                  design = model.matrix(~ Group_sense))

###############################################################################
## 2.  limma vooma‑treat パイプラインを “_sense” 接尾辞で実行
###############################################################################
designGLM_sense <- model.matrix(~ 0 + Group_sense)
v_sense   <- vooma(combat_sense, designGLM_sense, plot = FALSE)

fit_sense <- v_sense         |>
  lmFit(designGLM_sense)     |>
  contrasts.fit(cbind(Late_vs_Early = c(-1, 1))) |>
  eBayes(trend = TRUE)       |>
  treat(lfc = 0.30)

deg_sense <- topTreat(fit_sense, p.value = 0.05, number = Inf)

cat("6歳カットオフ  DEGs:", nrow(deg_sense), "\n")
# 6歳カットオフ  DEGs: 108


## 3.  必要なら QC 図も _sense で描画
plotRLE_micro(combat_sense, Group_sense,
              main = "RLE  (6-y cutoff, after ComBat)")

plotPCA_micro(combat_sense, Group_sense, Cohort,
              main = "PCA   (6-y cutoff, after ComBat)")


Cohort_sense <- Cohort
 fit_lfc_sense <- lmFit(v_sense, designGLM_sense) |>
   contrasts.fit(cbind(Late_vs_Early = c(-1,1))) |>
   eBayes(trend = TRUE) |>
   treat(lfc = 0.30) 

 table(Group_sense)
 #    Early  Late 
 # 6yr  19     13
 # 7yr  21     11
 # 8yr  25      7
 
## ヒートマップ用上位遺伝子
 # FDR TOP 100 Heatmap------------------------------------------------------------
 {
   suppressPackageStartupMessages({
     library(matrixStats); library(AnnotationDbi); library(org.Hs.eg.db); library(pheatmap)
   })
   
   ### 0. 前処理
   deg_df_sense <- as.data.frame(deg_sense)          # ← deg を data.frame 化
   
   ### 1. FDR 列名を安全に取得
   fdr_candidates_sense <- c("FDR", "adj.P.Val", "padj", "qvalue", "p.adj", "FDR.BH")
   fdr_col_sense <- intersect(fdr_candidates_sense, colnames(deg_df_sense))
   
   if (length(fdr_col_sense) == 0)
     stop("FDR に相当する列が見つかりません。列名を確認してください。")
   fdr_col_sense <- fdr_col_sense[1]                 # 複数あれば最初を採用
   
   ### 2. FDR 昇順で並べ替え & 上位 100
   deg_df_sorted_sense <- deg_df_sense[order(deg_df_sense[[fdr_col_sense]]), , drop = FALSE]
   top100_id_sense     <- head(rownames(deg_df_sorted_sense), 100)
   
   ### 2. 発現行列を抽出し、行 SD でフィルタ
   expr_top_sense <- combat_sense[top100_id_sense, , drop = FALSE]   # ← 修正
   expr_top_sense <- expr_top_sense[rowSds(expr_top_sense) > 0.15, ]
   n_gene_sense   <- nrow(expr_top_sense)
   
   ### 3. 行 Z‑score & クリップ
   expr_z_sense <- t(scale(t(expr_top_sense)))
   expr_z_sense <- pmin(pmax(expr_z_sense, -2.5), 2.5)
   
   ### 4. Entrez → SYMBOL  （必要かどうか判定）
   if (all(grepl("^[0-9]+$", rownames(expr_z_sense)))) {
     ## 行名が数字だけなら EntrezID と仮定して変換
     sym_sense <- mapIds(org.Hs.eg.db, rownames(expr_z_sense),
                         "SYMBOL", "ENTREZID", multiVals = "first")
     expr_z_sense <- expr_z_sense[!is.na(sym_sense), ]
     rownames(expr_z_sense) <- sym_sense[!is.na(sym_sense)]
   } else {
     ## すでに SYMBOL の場合
     sym_sense <- rownames(expr_z_sense)
   }
   
   ## 5. アノテーション
   ann_col_sense <- data.frame(Group = Group_sense,
                         Cohort = factor(Cohort_sense, levels = c("C1","C2","C3","C4")))
   rownames(ann_col_sense) <- colnames(expr_z_sense)
   ann_cols_sense <- list(
     Group  = c(Early = "#2C7BB6", Late = "#D7191C"),
     Cohort = c(C1 = "black", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")
   )
   
   ## 6. ヒートマップ
   pheatmap(
     expr_z_sense,
     scale = "none",
     color = colorRampPalette(c("navy","white","firebrick3"))(100),
     breaks = seq(-2.5, 2.5, length = 101),
     annotation_col    = ann_col_sense,
     annotation_colors = ann_cols_sense,
     clustering_distance_rows = "correlation",
     clustering_distance_cols = "correlation",
     show_rownames = F,
     border_color    = F,
     fontsize_row  = 5,   # 行名が潰れないよう小さめ
     main = sprintf("DEG heatmap | FDR Top 100  →  %d genes × %d samples (6yr)",
                    n_gene_sense, ncol(expr_z_sense))
   )
 }
 

##  deg_overlap_metrics()-------------------------------------------------------

## 5 歳 (既解析済み) を基準、6 歳カットを比較
{
  deg6 <- deg_sense
deg_overlap_metrics <- function(ref_deg,     # data.frame (Entrez rownames)
                                ref_label,   # e.g. "5yr"
                                test_deg,    # data.frame
                                test_label){ # e.g. "6yr"
  ## 1. 遺伝子集合
  ref_set  <- rownames(ref_deg)
  test_set <- rownames(test_deg)
  
  ## 2. 指標
  inter <- length(intersect(ref_set, test_set))
  union <- length(union(ref_set, test_set))
  
  jaccard   <- inter / union
  precision <- inter / length(test_set)
  recall    <- inter / length(ref_set)
  
  ## 3. sign(FC) 一致率
  common <- intersect(ref_set, test_set)
  dir_cons <- mean(sign(ref_deg[common, "logFC"]) ==
                     sign(test_deg[common, "logFC"]))
  
  ## 4. 出力を 1 行にまとめる
  data.frame(
    ref   = ref_label,
    test  = test_label,
    n_ref = length(ref_set),
    n_test= length(test_set),
    n_inter = inter,
    jaccard = round(jaccard, 3),
    precision = round(precision, 3),
    recall    = round(recall, 3),
    dir_cons  = round(dir_cons, 3),
    stringsAsFactors = FALSE
  )
}
deg_metrics6 <- deg_overlap_metrics(
  ref_deg   = deg,      ref_label  = "5yr",
  test_deg  = deg6,     test_label = "6yr"
)

print(deg_metrics6)
  }
  
# [6yr vs 5yr]  common = 77, Jaccard = 0.264   Direction-consistency = 100.00%   Precision = 71.30%   Recall = 29.50%


# Common DEGs heatmap ----------------------------------------------------------
{ ## 1. 共通 Entrez ID を取得（数値のみ）
  common_ent <- intersect(rownames(deg), rownames(deg6))   # ここは既に数字だけの想定
  
  ## 2. combat から共通行を抽出。
  ##    行名に "_何か" が付いている場合は先頭の数字だけ抜く
  idx <- sub("_.*", "", rownames(combat)) %in% common_ent
  mat_common  <- combat[idx, , drop = FALSE]
  
  ## 3. 行 Z スケール
  matZ_common <- t(scale(t(mat_common)))
  
  ## 4. 行名を正規化して SYMBOL へ。数字だけを keys に渡すのがポイント
  entrez_id <- sub("_.*", "", rownames(matZ_common))
  sym_common <- mapIds(org.Hs.eg.db,
                       keys      = entrez_id,
                       column    = "SYMBOL",
                       keytype   = "ENTREZID",
                       multiVals = "first")
  
  ## 5. 欠損を除外しつつ行名を置換。重複は make.unique() で解決
  keep <- !is.na(sym_common) & sym_common != ""
  matZ_common <- matZ_common[keep, ]
  rownames(matZ_common) <- make.unique(sym_common[keep])
  
  ## 6. ヒートマップ
  ann_col_common <- data.frame(Group = Group,
                               Cohort = Cohort,
                               row.names = colnames(matZ_common))
  
  pheatmap(matZ_common,
           annotation_col    = ann_col_common,
           annotation_colors = list(Group  = c(Early = "#2C7BB6",
                                               Late  = "#D7191C"),
                                    Cohort = c(C1="#000000",
                                               C2="grey40",
                                               C3="#4DAF4A",
                                               C4="purple")),
           color             = colorRampPalette(c("navy","white","firebrick3"))(50),
           cluster_rows      = TRUE,
           cluster_cols      = TRUE,
           show_rownames     = TRUE,
           fontsize_row      = 5,
           border_color      = "grey60",
           main = sprintf("Common %d DEGs (6yr vs 5yr)", nrow(matZ_common)))
 }


 # volcano plot  ―― 6 歳カットオフ解析（_sense） ------------------------------
 {
   library(dplyr)
   library(ggplot2)
   library(ggrepel)
   library(AnnotationDbi)
   library(org.Hs.eg.db)
   
   ## 0. limma‑treat 全出力
   tbl_all_sense <- topTreat(fit_sense, p.value = 1, number = Inf)
   
   ## 1. ラベルと有意判定列を追加
   tbl_plot_sense <- tbl_all_sense %>%
     mutate(EntrezID   = rownames(.),
            GeneSymbol = AnnotationDbi::mapIds(
              org.Hs.eg.db,
              keys      = EntrezID,
              column    = "SYMBOL",
              keytype   = "ENTREZID",
              multiVals = "first"),
            GeneSymbol = ifelse(is.na(GeneSymbol), EntrezID, GeneSymbol),
            signif_sense = case_when(
              adj.P.Val < 0.05 & logFC >=  0.30 ~ "Up",
              adj.P.Val < 0.05 & logFC <= -0.30 ~ "Down",
              TRUE                              ~ "NS"))
   
   ## 1. 除外パターンを明示
   pattern_exclude <- "^(LOC|LINC|CH[0-9]|AC[0-9]|AL[0-9]|CT[0-9]|RP[0-9])"
   
   ## 2. Up/Down それぞれ上位 5 件を抽出（除外パターンを弾く）
   top_labs_sense <- bind_rows(
     tbl_plot_sense %>%
       filter(signif_sense == "Up",
              !grepl(pattern_exclude, GeneSymbol)) %>%
       slice_min(adj.P.Val, n = 5),
     tbl_plot_sense %>%
       filter(signif_sense == "Down",
              !grepl(pattern_exclude, GeneSymbol)) %>%
       slice_min(adj.P.Val, n = 5)
   )
   
   ## 3. ラベル描画はこの top_labs_sense だけに限定
   geom_text_repel(data = top_labs_sense,
                   aes(label = GeneSymbol), size = 3,
                   max.overlaps = Inf)
   
   
   ## 3. Volcano プロット
   ggplot(tbl_plot_sense, aes(logFC, -log10(adj.P.Val))) +
     geom_point(aes(colour = signif_sense), size = 1.6, alpha = 0.8) +
     scale_colour_manual(values = c(Up = "#D7191C",
                                    Down = "#2C7BB6",
                                    NS = "grey70")) +
     geom_vline(xintercept = c(-0.30, 0.30),
                linetype = "dashed", colour = "grey40") +
     geom_hline(yintercept = -log10(0.05),
                linetype = "dashed", colour = "grey40") +
     geom_text_repel(data = top_labs_sense,
                     aes(label = GeneSymbol), size = 3,
                     max.overlaps = Inf) +
     labs(title = sprintf(
       "Volcano plot | %d DEGs (|logFC| ≥ 0.30, FDR ≤ 0.05, 6-yr cutoff)",
       sum(tbl_plot_sense$signif_sense != "NS")),
       x = "log2 fold-change",
       y = expression(-log[10]~FDR)) +
     theme_bw(base_size = 12) +
     theme(legend.title = element_blank())
 }
 
 library(biomaRt)
 ens <- useEnsembl(biomart = "genes",
                   dataset = "hsapiens_gene_ensembl")
 ann <- getBM(c("entrezgene_id", "external_gene_name"),
              filters = "entrezgene_id",
              values  = c("105379499","102724701"),
              mart    = ens)
 print(ann)
 
 
# 遺伝子名列挙------------------------------------------------------------------
 {
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

## 1. 共通 ID を取得（base::intersect で衝突回避）
common_id <- base::intersect(rownames(deg),    # 5yr
                             rownames(deg6))   # 6yr
stopifnot(length(common_id) == 77)

## 2. 5歳カット側の logFC・FDR を抜粋
tbl_common <- data.frame(
  SYMBOL      = AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys      = common_id,
                                      column    = "SYMBOL",
                                      keytype   = "ENTREZID",
                                      multiVals = "first"),
  logFC       = deg[common_id, "logFC"],
  adj.P.Val   = deg[common_id, "adj.P.Val"],
  row.names   = NULL,                         # 行番号を消す
  check.names = FALSE
)

## 3. FDR 昇順で並べ替え（お好みで abs(logFC) も可）
tbl_common <- tbl_common[order(tbl_common$adj.P.Val), ]

## 4. コンソールに小数 3 桁で出力
print(tbl_common, digits = 3, row.names = FALSE)
}

# E‑7.  STRING PPI  (6‑yr cutoff)-----------------------------------------------
suppressPackageStartupMessages({
  library(dplyr);   library(tibble)
  library(AnnotationDbi);  library(org.Hs.eg.db)
  library(STRINGdb);       library(igraph);    library(ggraph)
})

## 0. 依存オブジェクトの存在を確認
stopifnot(exists("deg_sense"))          # limma‑treat (5yr) DEGs
stopifnot(exists("combat_sense"))       # ComBat 正規化行列 (rows = probe, cols = 40)
stopifnot(exists("Group_sense"), exists("Cohort_sense"))

## 1. STRING ローカルデータ
cache_dir <- "C:/DMD_project"
files_needed <- c(
  "9606.protein.aliases.v11.5.txt.gz",
  "9606.protein.info.v11.5.txt.gz",
  "9606.protein.links.detailed.v11.5.txt.gz",
  "9606.protein.links.v11.5.txt.gz"
)
missing <- files_needed[!file.exists(file.path(cache_dir, files_needed))]
if (length(missing) > 0)
  stop("❌ STRING ファイルが不足：\n  ",
       paste(missing, collapse = "\n  "))
{
  str_db_sense <- STRINGdb$new(version = "11.5",
                               species = 9606,
                               score_threshold = 500,
                               input_directory = cache_dir)
  
  ## 2. Entrez → SYMBOL & ケラチン系除外
  skip_pat <- paste0(
    "^KRT","|^KRTAP","|^LCE","|^FLG$","|^SPRR","|^IVL$",
    "|^DSG","|^DSP$","|^TACSTD2$","|^OR\\d+","|^TAS2R\\d+")
  
  deg_tbl_sense <- deg_sense %>%                         # 5yr‑DEGs
    tibble::rownames_to_column("EntrezID") %>%           
    dplyr::mutate(Symbol = AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = EntrezID,
      column    = "SYMBOL",
      keytype   = "ENTREZID",
      multiVals = "first")) %>% 
    dplyr::filter(!grepl(skip_pat, Symbol) & !is.na(Symbol)) %>% 
    dplyr::select(Symbol, logFC, adj.P.Val)
  
  ## 3. STRING ID マッピング & ネットワーク
  map_tbl_sense <- str_db_sense$map(deg_tbl_sense,
                                    "Symbol", removeUnmappedRows = TRUE)
  cat(sprintf("🧬 STRING mapping (6yr): %d / %d genes hit\n",
              nrow(map_tbl_sense), nrow(deg_tbl_sense)))
  
  g_net_sense <- str_db_sense$get_subnetwork(map_tbl_sense$STRING_id)
  V(g_net_sense)$Symbol <- map_tbl_sense$Symbol[
    match(V(g_net_sense)$name, map_tbl_sense$STRING_id)]
  V(g_net_sense)$logFC  <- deg_tbl_sense$logFC[
    match(V(g_net_sense)$Symbol, deg_tbl_sense$Symbol)]
  
  #
  hub_id_sense <- names(sort(degree(g_net_sense), decreasing = TRUE))[1:min(30, vcount(g_net_sense))]
  V(g_net_sense)$is_hub <- V(g_net_sense)$name %in% hub_id_sense
  
  ## 4. ネットワークプロット（必要なら
  pal_fun <- colorRampPalette(c("navy","white","firebrick"))
  }

ggraph(g_net_sense, layout = "fr") +
  geom_edge_link(aes(width = combined_score / 1000),
                 colour = "grey70", alpha = 0.3) +
  scale_edge_width(range = c(0.2, 1.5),          # 幅調整
                   name  = "score") +  # ← ここで好きな文字列
  geom_node_point(aes(fill = logFC, size = is_hub),
                  shape = 21, colour = "black") +
  geom_node_text(aes(label = Symbol), size = 3,
                 repel = TRUE, check_overlap = TRUE) +
  scale_fill_gradientn(colours = pal_fun(100), name = "logFC") +
  scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2), guide = "none") +
  theme_void(base_size = 12) +
  labs(title = sprintf("STRING PPI (Sensitivity Analysis) | 6-yr (score > 0.5) |  %d genes",
                       vcount(g_net_sense)))

# ノード数
vcount(g_net_sense)   # → 95

# エッジ数
ecount(g_net_sense)   # → 60

# E‑8.  Hub‑gene heatmap (6‑yr cutoff)  -----------------------------------
suppressPackageStartupMessages({library(pheatmap)})

## 1)  STRING ネットワークに載っている全 Symbol を取得
all_syms_sense <- unique(V(g_net_sense)$Symbol)   # ← ここが 95 個になるはず

## 2)  Combat 行列を SYMBOL 名に変換
entrez_ids_sense <- sub("_.*", "", rownames(combat_sense))
sym_vec_sense <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                       keys       = entrez_ids_sense,
                                       column     = "SYMBOL",
                                       keytype    = "ENTREZID",
                                       multiVals  = "first")

expr_sym_sense <- combat_sense[!is.na(sym_vec_sense), ]
rownames(expr_sym_sense) <- sym_vec_sense[!is.na(sym_vec_sense)]

## 3)  ネットワーク Symbol と発現行列の交差
all_syms_sense <- intersect(all_syms_sense, rownames(expr_sym_sense))

## 4)  行 Z‑score
mat_hub_z_sense <- t(scale(t(expr_sym_sense[all_syms_sense, ])))
mat_hub_z_sense <- mat_hub_z_sense[rowSums(is.na(mat_hub_z_sense)) == 0, ]

## 5)  列アノテーション
ann_col_sense <- data.frame(
  Group  = Group_sense,
  Cohort = factor(Cohort_sense, levels = c("C1","C2","C3","C4"))
)
rownames(ann_col_sense) <- colnames(mat_hub_z_sense)

ann_cols_sense <- list(
  Group  = c(Early = "#2C7BB6", Late = "#D7191C"),
  Cohort = c(C1 = "black", C2 = "grey40", C3 = "#4DAF4A", C4 = "purple")
)

## 6)  ヒートマップ
pheatmap(
  mat_hub_z_sense,
  color              = colorRampPalette(c("navy","white","firebrick3"))(50),
  annotation_col     = ann_col_sense,
  annotation_colors  = ann_cols_sense,
  cluster_rows       = TRUE,
  cluster_cols       = TRUE,
  show_rownames      = TRUE,
  fontsize_row       = 5,
  main = sprintf("STRING hub genes | %d nodes × %d samples (6-yr)",
                 nrow(mat_hub_z_sense), ncol(mat_hub_z_sense))
)


# E‑9.  Hub-gene PCA/RLE (6-yr cutoff) -----------------------------------------
## 感度解析用 PCA 結果オブジェクト : pc_sense
## 感度解析用グループ / コホート : Group_sense, Cohort_sense

##  E‑9.  Hub‑gene PCA  (6‑yr, *_sense)
{
  suppressPackageStartupMessages({library(matrixStats)})
  
  ## mat_hub_z_sense はヒートマップ作成で定義済み
  pc_sense <- prcomp(t(mat_hub_z_sense))      # ← 行＝gene, 列＝sample
  
  ## aesthetic ベクトル
  col_vec_sense <- c(Early = "#2C7BB6", Late = "#D7191C")[Group_sense]
  pch_vec_sense <- c(C1 = 17, C2 = 7, C3 = 3, C4 = 5)[Cohort_sense]
  
  ## プロット
  plot(pc_sense$x[, 1], pc_sense$x[, 2],
       col = col_vec_sense,
       pch = pch_vec_sense,
       xlab = sprintf("PC1 (%.1f%%)", 100*pc_sense$sdev[1]^2/sum(pc_sense$sdev^2)),
       ylab = sprintf("PC2 (%.1f%%)", 100*pc_sense$sdev[2]^2/sum(pc_sense$sdev^2)),
       main = "Hub-gene PCA (6-yr)")
  
  legend("topleft",  legend = levels(Group_sense),
         col = c("#2C7BB6","#D7191C"), pch = 16, title = "Group", bty = "n")
  legend("topright", legend = levels(Cohort_sense),
         pch  = c(17,7,3,5), title = "Cohort", bty = "n")
}


##  Hub‑gene RLE  (6‑yr, *_sense)
{
  suppressPackageStartupMessages({library(matrixStats)})
  
  ## 1. Combat 行列を SYMBOL 行名へ換装（既に作成済みだが冗長を整理）
  entrez_ids_sense <- sub("_.*", "", rownames(combat_sense))
  sym_vec_sense    <- mapIds(org.Hs.eg.db,
                             keys       = entrez_ids_sense,
                             column     = "SYMBOL",
                             keytype    = "ENTREZID",
                             multiVals  = "first")
  
  expr_sym_sense <- combat_sense[!is.na(sym_vec_sense), ]
  rownames(expr_sym_sense) <- sym_vec_sense[!is.na(sym_vec_sense)]
  
  ## 2. ハブ遺伝子だけ抽出し RLE 計算
  expr_hub_sense <- expr_sym_sense[ rownames(expr_sym_sense) %in% rownames(mat_hub_z_sense), ]
  rle_mat_sense <- expr_hub_sense -
    matrixStats::rowMedians(as.matrix(expr_hub_sense), na.rm = TRUE)
  
  
  ## 3. 可視化
  boxplot(as.data.frame(rle_mat_sense),
          las      = 2,
          outline  = FALSE,
          col      = c(C1="black", C2="grey40", C3="#4DAF4A", C4="purple")[Cohort_sense],
          ylab     = "RLE (log2, gene-median centred)",
          main     = "Hub-gene RLE (6-yr)")
  
  legend("topright", legend = levels(Cohort_sense),
         fill = c("grey","yellow","#4DAF4A","purple"), bty = "n")
  
}


# E-10. KEGG, Reactome, GO-BP (6-yr cutoff)-------------------------------------
{
  suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(stringr)
  # enrichplot は clusterProfiler の依存として自動ロード済み
})


## 0. 前提チェック（5‑yr 処理後に作られたオブジェクト）
stopifnot(exists("combat_sense"),    # ComBat 後行列
          exists("deg_sense"))     # limma‑treat で抽出した DEGs

## 1. 背景遺伝子と解析対象遺伝子
bg_ids_sense  <- unique(sub("_.*", "", rownames(combat_sense)))   # 22k 強
deg_ids_sense <- rownames(deg_sense)                            # 442 Entrez

# 2. KEGG ORA
kegg_res_sense <- enrichKEGG(
  gene         = deg_ids_sense,
  universe     = bg_ids_sense,
  organism     = "hsa",
  keyType      = "ncbi-geneid",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.10
)

## ── 病名経路を除去
kegg_clean_sense <- kegg_res_sense@result |>
  dplyr::filter(!grepl("^hsa05|disease", ID))          # “05xxx”＝Disease 系
kegg_clean_sense <- new("enrichResult",
                        kegg_res_sense,
                        result = kegg_clean_sense)

# 3. Reactome ORA
react_res_sense <- enrichPathway(
  gene         = deg_ids_sense,
  universe     = bg_ids_sense,
  organism     = "human",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.10
)

## 4. GO‑BP ORA（冗長削減
go_res_sense <- enrichGO(
  gene          = deg_ids_sense,
  universe      = bg_ids_sense,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
) |> simplify(cutoff = 0.10, by = "p.adjust", select_fun = min)

## 5. 可視化ユーティリティ（既存関数を再利用
# pub_dotplot() / pub_barplot() は元のまま利用可
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

pub_barplot <- function(eres, n = 15, fdr_cut = 0.10,
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

## 6. 描画例
pub_dotplot(react_res_sense, n = 20, fdr_cut = 0.05,
            title = "Reactome enrichment (6-yr)")

pub_dotplot(go_res_sense,    n = 20, fdr_cut = 0.05,
            title = "GO-BP enrichment (6-yr)")

pub_dotplot(kegg_clean_sense, n = 20, fdr_cut = 0.05,
            title = "KEGG enrichment (6-yr)")

# 上位 2 経路だけ除く場合（任意）
# kegg_clean_sense@result <- kegg_clean_sense@result[-c(1,2), ]


# E-11. GSEA (Reactome) (6-yr cutoff)---------------------------------------------------
## ランク付けベクトル（6‑yr モデルの limma‑treat 係数を使用）
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
  SHOW_TOP       <- 30       # ridgeplot に表示する経路数
  X_BREAKS       <- seq(-2, 2, 1)
  
  # 1. モデルフィット & ランクベクトル
  designGLM_sense <- model.matrix(~ 0 + Group_sense)
  fit.lfc_sense   <- v_sense |>
    lmFit(designGLM_sense) |>
    contrasts.fit(cbind(Late_vs_Early = c(-1, 1))) |>
    eBayes(trend = TRUE) |>
    treat(lfc = LFC_CUTOFF)
  
  # EntrezID を取り出して geneList を作成
  fc_vec_sense <- setNames(fit.lfc_sense$coefficients,
                     sub("_.*", "", rownames(fit.lfc_sense$coefficients)))
  geneList_sense <- sort(fc_vec_sense, decreasing = TRUE)
  
  # 2. BiocParallel を安全に設定
  register(SerialParam())      # ← 安定運用。並列化するなら SnowParam() へ変更
  # 例: register(SnowParam(workers = 4, type = "SOCK", progressbar = TRUE))
  
  #3. GSEA 実行
  set.seed(1234)
  gsea_react_sense <- gsePathway(
    geneList     = geneList_sense,
    organism     = "human",
    pvalueCutoff = P_CUTOFF,   # 有意経路が少なければ 1 にして上位表示
    maxGSSize    = 3000,
    minGSSize    = 10,
    eps          = 0           # fgseaMultilevel を使用（高速・高精度）
  )
}

#4. 可視化関数
{
  RIDGE_LABEL_SZ <- 4        # y 軸ラベル文字サイズ
  ridgeplot2_sense <- function(gsea_obj,
                         showCategory      = 20,
                         pathway_text_size = 4,
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
  {p_ridge_sense <- ridgeplot2_sense(
    gsea_obj          = gsea_react_sense,
    showCategory      = SHOW_TOP,
    pathway_text_size = RIDGE_LABEL_SZ,   # ← y 軸ラベル（経路名）のサイズ
    x_breaks          = X_BREAKS
  ) +
      ggtitle(sprintf("Reactome enrichment (|log2FC| > %.2f, 6-yr)", LFC_CUTOFF)) +
      theme(
        plot.title  = element_text(hjust = 0.5, size = 14,  margin = margin(b = 6)),  # タイトル
        axis.text.x = element_text(size = 10),   # x 軸目盛
        axis.title.x= element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title= element_text(size = 10)
      )
    print(p_ridge_sense)
  }
}


# ヒットした全経路をコンソール表示する
{
  ## ❶ まず FDR<0.05 でフィルタリング
  tbl_short_sense <- as.data.frame(gsea_react_sense) |>
    dplyr::filter(p.adjust < 0.05) |>
    dplyr::select(Description, NES, p.adjust) |>
    dplyr::arrange(p.adjust)
  
  ## ❷ 数値列だけを有効桁 3 桁に丸めて表示
  tbl_print_sense <- tbl_short_sense |>
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),
        ~ signif(.x, 3)          # ← 数値列にだけ適用
      )
    )
  sig3_sense <- function(x, sig = 3){
    vapply(x, function(v){
      if (v == 0) return("0")
      ## 桁数 = sig - 1 - floor(log10(v))
      dec_sense <- max(sig - 1 - floor(log10(abs(v))), 0)
      formatC(v, digits = dec_sense, format = "f", drop0trailing = FALSE)
    }, FUN.VALUE = "")
  }
  tbl_print_sense <- tbl_short_sense %>%
    dplyr::mutate(p.adjust = sig3_sense(p.adjust))
  print(tbl_print_sense, row.names = FALSE)
  
}

# 8歳ではDEG399であっても、GSEAではわずか20経路程度、統一性の無いヒットがあるだけであり、
# さらにはORCで有意経路無し。


library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

## 1. テーブル化して指標計算 -----------------------------------------------
tbl_main  <- as.data.frame(gsea_react) |>
  dplyr::select(ID, NES, p.adjust)

tbl_sense <- as.data.frame(gsea_react_sense) |>
  dplyr::select(ID, NES, p.adjust) |>
  dplyr::rename(NES_sense = NES, p.sense = p.adjust)

sig_main  <- tbl_main  |> filter(p.adjust < 0.05)
sig_sense <- tbl_sense |> filter(p.sense  < 0.05)

overlap <- inner_join(sig_main, sig_sense, by = "ID")
jaccard <- nrow(overlap) /
  (nrow(sig_main) + nrow(sig_sense) - nrow(overlap))
cor_NES <- cor(overlap$NES, overlap$NES_sense, method = "spearman")

cat(sprintf(
  "共有経路 %d 本 | Jaccard %.2f | NES 相関 ρ = %.3f\n",
  nrow(overlap), jaccard, cor_NES))
# 共有経路 41 本 | Jaccard 0.42 | NES 相関 ρ = 0.912

## 2. NES 散布図 -------------------------------------------------------------
ggplot(overlap, aes(NES, NES_sense)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  labs(x = "NES (Main)", y = "NES (6‑yr)",
       title = "Pathway NES concordance") +
  theme_bw()

## 3. compareCluster dotplot ------------------------------------------------
comp_sense <- compareCluster(
  geneCluster = list(Main = geneList,
                     `6yr` = geneList_sense),
  fun         = "gsePathway",
  pvalueCutoff= 0.05,
  organism    = "human"
)

enrichplot::dotplot(comp_sense,
                    showCategory = 20,
                    font.size    = 6) +
  ggtitle("Reactome GSEA – Main vs 6-yr") +
  theme(
    axis.text.y = element_text(size = 6, margin = margin(r = 2)),
    plot.title  = element_text(hjust = 0.5, size = 13)
  )

dotplot(comp_sense,
        showCategory = 20,
        font.size    = 6    # ← デフォルト 12 を小さく
) +
  ggtitle("Reactome GSEA – Main vs 6-yr)") +
  theme(
    axis.text.y  = element_text(size = 6,   # さらに細かく調整したい場合
                                margin = margin(r = 2)),
    plot.title   = element_text(hjust = 0.5, size = 13)
  )


##------------------------------------------
## 保存：現在のワークスペースまるごと
##------------------------------------------
getwd()         # 例: "C:/DMD_project"
# 保存先フォルダを指定（例：プロジェクト直下の /cache）
#dir.create("cache", showWarnings = FALSE)
#save.image(file = "cache/workspace_2025.728.RData")

##------------------------------------------
## 復元
##------------------------------------------
# load("cache/cache/workspace_2025.728.RData")
# 読み込むと、その時点のオブジェクトが Environment に再現されます