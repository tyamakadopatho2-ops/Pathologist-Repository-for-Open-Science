##### D: --------------------------------------------------Subgroup Analysis -----------------------------------------------------------------
##  Sub‑group analysis : C1 (U133_Plus2.0) only

suppressPackageStartupMessages({
  library(limma)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

## 1.  サブセット抽出  ---------------------------------------------------------
##     C1 (=U133 Plus2.0) 17 例のみ
idx_C1      <- Cohort == "C1"          # TRUE/FALSE ベクトル
expr_C1     <- combat[ , idx_C1]       # 行:遺伝子  ×  列:17 sample
Group_C1    <- droplevels(Group[idx_C1])   # Early 12 / Late 5
table(Group_C1)
nrow(expr_C1) # 18509

## 2.  limma‑vooma → treat無し-----------------------------
design_C1 <- model.matrix(~0 + Group_C1)   # 列名: Early / Late
v_C1      <- vooma(expr_C1, design_C1, plot = FALSE)

fit_C1 <- lmFit(v_C1, design_C1) |>
  contrasts.fit(cbind(Late_vs_Early = c(-1,1))) |>
  eBayes(trend = TRUE)              # treat なし
deg_C1 <- topTable(fit_C1, p.value = 0.05, number = Inf)
cat("DEG  (C1 only):", nrow(deg_C1), "\n")
# DEG  (C1 only): 99


## 4.  ヒートマップ-------------------------------------------------------------------
{
  ## 0.  Entrez → SYMBOL へ行名変換（pheatmap を呼ぶ前に挿入）
  
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  hit_genes   <- rownames(deg_C1)
  hit_sub     <- expr_C1[hit_genes , , drop = FALSE]      # 99 × 17
  hit_sub_z   <- t(scale(t(hit_sub)))                     # 行 Z‑score
  
  ## 2. Entrez → SYMBOL 変換（数字のみの場合）
  if (all(grepl("^[0-9]+$", rownames(hit_sub_z)))) {
    sym_C1 <- mapIds(org.Hs.eg.db,
                     keys      = rownames(hit_sub_z),
                     column    = "SYMBOL",
                     keytype   = "ENTREZID",
                     multiVals = "first")
    ## 欠損を除外
    keep <- !is.na(sym_C1) & sym_C1 != ""
    hit_sub_z <- hit_sub_z[keep, ]
    sym_C1    <- sym_C1[keep]
    
    ## 重複シンボルにサフィックス _1,_2… を付与
    dup_idx <- duplicated(sym_C1) | duplicated(sym_C1, fromLast = TRUE)
    if (any(dup_idx)) {
      sym_C1[dup_idx] <- paste0(sym_C1[dup_idx], "_",
                                ave(seq_along(sym_C1[dup_idx]),
                                    sym_C1[dup_idx],
                                    FUN = seq_along))
    }
    rownames(hit_sub_z) <- sym_C1
  }
  
  ## ② 重複シンボルにサフィックスを付けてユニーク化
  dup_C1 <- duplicated(sym_C1) | duplicated(sym_C1, fromLast = TRUE)
  sym_C1[dup_C1] <- paste0(sym_C1[dup_C1], "_", rownames(hit_sub_z)[dup_C1])
  
  ## ③ 行名としてセット
  rownames(hit_sub_z) <- ifelse(is.na(sym_C1), rownames(hit_sub_z), sym_C1)
  
  
  ## 1.  これ以降は既存の pheatmap 呼び出しをそのまま
  pheatmap(hit_sub_z,
           color            = colorRampPalette(c("navy","white","firebrick3"))(50),
           annotation_col   = ann_col,
           annotation_colors= ann_cols,
           show_rownames    = TRUE,
           fontsize_row     = 4.5,        # 行ラベルが長ければさらに小さく
           cluster_rows     = TRUE,
           cluster_cols     = TRUE,
           main             = sprintf("DEGs (%d genes) – Subgroup Analysis",
                                      nrow(hit_sub_z)))
}

# 主解析と一致するDEGのヒートマップ---------------------------------------------
## 0. ライブラリと前提オブジェクト

{
  suppressPackageStartupMessages({
  library(matrixStats)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(pheatmap)
})

## 1. 主解析とサブ解析の共通 DEGs を抽出
common_ent <- intersect(rownames(deg), rownames(deg_C1))
if (!length(common_ent))
  stop("主解析と一致する DEGs が 0 行でした。行名を確認してください。")

## 2. 発現行列を取り出し、行 Z‑score へ
expr_common <- expr_C1[common_ent, , drop = FALSE]     # ← cohort 全体なら combat[]
expr_z      <- t(scale(t(expr_common)))                # 行ごとに中心化・標準化
expr_z      <- pmin(pmax(expr_z, -2.5), 2.5)           # 外れ値を ±2.5 にクリップ

## 3. Entrez → SYMBOL 変換
sym <- mapIds(org.Hs.eg.db,
              keys      = rownames(expr_z),
              column    = "SYMBOL",
              keytype   = "ENTREZID",
              multiVals = "first")

## 欠損を除外／重複は _1,_2… を付与
keep <- !is.na(sym) & sym != ""
expr_z <- expr_z[keep, ]
sym    <- sym[keep]
rownames(expr_z) <- make.unique(ifelse(sym == "", rownames(expr_z), sym))

## 4. 列アノテーション（既存オブジェクトを再利用）
##    ann_col  : data.frame(Group, Cohort)
##    ann_cols : list(Group = c(...), Cohort = c(...))
if (!exists("ann_col") || !exists("ann_cols"))
  stop("ann_col / ann_cols が見つかりません。列アノテーションを定義してください。")

## 5. ヒートマップ描画
pheatmap(
  expr_z,
  annotation_col    = ann_col[colnames(expr_z), ],
  annotation_colors = ann_cols,
  color             = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  fontsize_row      = 4.5,
  border_color      = "grey60",
  main = sprintf("Common %d DEGs (Main ∩ Subgroup)", nrow(expr_z))
)
}


## Volcano plot  –  Sub‑group (C1 only)  --------------------------------------
##   * 既存オブジェクト名と衝突しないよう “*_C1” 接尾辞で統一 *
{
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
  })
  
  ## 0. limma 結果テーブルを取得（fit_C1 は既に作成済み）
  tbl_C1_all <- topTable(fit_C1, p.value = 1, number = Inf)   # 全行
  
  ## 1. ラベル列を追加
  tbl_C1_plot <- tbl_C1_all %>%
    mutate(
      EntrezID   = rownames(.),
      GeneSymbol = AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys      = EntrezID,
        column    = "SYMBOL",
        keytype   = "ENTREZID",
        multiVals = "first"),
      GeneSymbol = ifelse(is.na(GeneSymbol), EntrezID, GeneSymbol),
      signif_C1  = case_when(
        adj.P.Val < 0.05 & logFC >=  0.30 ~ "Up",
        adj.P.Val < 0.05 & logFC <= -0.30 ~ "Down",
        TRUE                              ~ "NS")
    )
  
  ## 2. ラベルに使う上位 10 遺伝子（有意点から）
  top_labs_C1 <- tbl_C1_plot %>%
    filter(signif_C1 != "NS") %>%
    slice_min(adj.P.Val, n = 10)
  
  ## 3. Volcano プロット
  p_volcano_C1 <- ggplot(tbl_C1_plot,
                         aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(colour = signif_C1), size = 1.6, alpha = 0.8) +
    scale_colour_manual(values = c(Up   = "#D7191C",
                                   Down = "#2C7BB6",
                                   NS   = "grey70"),
                        name = NULL) +
    geom_vline(xintercept = c(-0.30, 0.30), linetype = "dashed",
               colour = "grey40") +
    geom_hline(yintercept = -log10(0.05),    linetype = "dashed",
               colour = "grey40") +
    geom_text_repel(data = top_labs_C1,
                    aes(label = GeneSymbol), size = 3,
                    max.overlaps = Inf) +
    labs(title = sprintf("Volcano plot (C1 only) | %d DEGs  (|logFC| ≥ 0.30, FDR ≤ 0.05)",
                         sum(tbl_C1_plot$signif_C1 != "NS")),
         x = "log2 fold‑change", y = expression(-log[10]~FDR)) +
    theme_bw(base_size = 12)
  
  print(p_volcano_C1)
}


## 3.  メイン DEG 集合との log2FC 整合性を評価----------------------------------
##     ① 両方とも EntrezID 行名へ換装
{deg_main_ent <- deg
rownames(deg_main_ent) <- sub("_.*", "", rownames(deg_main_ent))   # “12345_xyz”→“12345”

deg_C1_ent   <- deg_C1
rownames(deg_C1_ent)   <- sub("_.*", "", rownames(deg_C1_ent))

## ② 2 つの表を EntrezID でマージ
tab_main <- data.frame(EntrezID = rownames(deg_main_ent),
                       logFC_all = deg_main_ent$logFC)
tab_C1   <- data.frame(EntrezID = rownames(deg_C1_ent),
                       logFC_C1 = deg_C1_ent$logFC)

merge_tab <- merge(tab_main, tab_C1, by = "EntrezID", sort = FALSE)
plot(merge_tab$logFC_all, merge_tab$logFC_C1,
     xlab = "logFC (Main)", ylab = "logFC (C1)",
     pch = 16, col = "grey50")
abline(0, 1, col = "red")
}

## 4.  一致率
##  A. 一致（Jaccard）と方向率
{
  common_id <- intersect(rownames(deg_main_ent), rownames(deg_C1_ent))
  jaccard   <- length(common_id) / (nrow(deg_main_ent) + nrow(deg_C1_ent) - length(common_id))
  
  same_dir <- sum(sign(deg_main_ent[common_id,"logFC"]) ==
                    sign(deg_C1_ent  [common_id,"logFC"])) / length(common_id)
  
  precision <- length(common_id) / nrow(deg_C1_ent)   # 83/99
  recall    <- length(common_id) / nrow(deg_main_ent) # 83/261

  
  cat(sprintf(
    "Jaccard = %.3f   Concordance of effect direction = %.2f%%   Precision = %.2f%%   Recall = %.2f%%\n",
    jaccard,
    same_dir * 100,
    precision * 100,
    recall * 100
  ))

# Jaccard = 0.300   Concordance of effect direction = 100.00%
# Precision = 83.84%   Recall = 31.80%

## --- 4. 一致率ブロックの末尾に追記
n_common <- length(common_id)
pct_C1   <- precision * 100          # = n_common / nrow(deg_C1_ent) × 100
pct_all  <- recall    * 100          # = n_common / nrow(deg_main_ent) × 100

cat(sprintf(
  "Common DEGs : %d   (%.1f%% of C1‑DEG, %.1f%% of All‑DEG)\n",
  n_common, pct_C1, pct_all))
}
# Common DEGs : 83   (83.8% of C1‑DEG, 31.8% of All‑DEG)


# venn graph ------------------------------------------------------------------
## 0.  必要パッケージ
{ library(ggplot2)
  library(gridExtra)    # grid.arrange()
  library(grid)         # grobTree(), gList()
  library(VennDiagram)
  
  
  ## 1.  棒グラフオブジェクト
  metrics <- data.frame(
    Metric = c("Jaccard", "Precision", "Recall", "Concordance"),
    Value  = c(jaccard, precision, recall, same_dir)
  )
  
  metrics_plot <- ggplot(metrics, aes(Metric, Value*100)) +
    geom_col(fill = "#1f78b4") +
    geom_text(aes(label = sprintf("%.1f%%", Value*100)),
              vjust = -0.3, size = 4) +
    scale_y_continuous(limits = c(0, 100),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(y = "Percentage / Score", x = NULL,
         title = "Overlap & Consistency between All-DEG and C1-DEG") +
    theme_minimal(base_size = 12)
  
  ## 2.  Venn を “単一 grob” に変換
  venn_list <- draw.pairwise.venn(
    area1      = nrow(deg_main_ent),
    area2      = nrow(deg_C1_ent),
    cross.area = n_common,
    category   = c(sprintf("All-DEG: %d", nrow(deg_main_ent)),
                   sprintf("C1-DEG: %d", nrow(deg_C1_ent))),
    fill       = c("#4daf4a", "#377eb8"),
    alpha      = 0.5,
    scaled     = TRUE,
    cat.cex    = 1.2, cex = 1.2,
    cat.pos  = c(180, 0),     # 左円: 真横(180°)、右円: 真横(0°)
    cat.dist = c(0.03, 0.01), # 円からラベルまでの距離（0 = 縁）
    cat.just = list(c(0, 0.5), c(1, 0.5)),  # 左ラベルは右寄せ、右ラベルは左寄せ
    ind        = FALSE          # “描画せず” grob のリストを返す
  )
  
  ## ---- ここがポイント：リスト → gTree に束ねる ----
  venn_grob <- gTree(children = do.call(gList, venn_list))
  grid.newpage()   # 新規ページ
  grid.arrange(metrics_plot, venn_grob,
               ncol = 2,
               widths = c(1.2, 1))
}

grid.newpage()   # または dev.off() でデバイスを閉じる



# GSEA -------------------------------------------------------------------------
###############################################################################
## GSEA（Reactome）–  C1 サブグループ専用  (_C1 suffix)
###############################################################################
{suppressPackageStartupMessages({
  library(limma)
  library(clusterProfiler)
  library(BiocParallel)
  library(ggplot2)
  library(ggridges)
})

## 0. 解析パラメータ
LFC_CUTOFF_C1 <- 0      # treat() のしきい値
P_CUTOFF_C1   <- 0.05      # gsePathway の FDR
SHOW_TOP_C1   <- 30        # ridgeplot に表示する経路数
X_BREAKS_C1   <- seq(-2, 2, 1)

## 1. モデルフィット & ランクベクトル
designGLM_C1 <- model.matrix(~ 0 + Group_C1)   # ← 先生が定義済みの Group_C1
v_C1         <- vooma(expr_C1, designGLM_C1, plot = FALSE)

fit.lfc_C1 <- v_C1 |>
  lmFit(designGLM_C1) |>
  contrasts.fit(cbind(Late_vs_Early = c(-1, 1))) |>
  eBayes(trend = TRUE)        # ← treat を省くか lfc=0 でも同じ

geneList_C1 <- setNames(fit.lfc_C1$t,
                        sub("_.*","", rownames(fit.lfc_C1$t))) |>
  sort(decreasing = TRUE)

## 2. 並列設定（安全運用なら SerialParam）
register(SerialParam())

## 3. GSEA 実行
set.seed(1234)
gsea_react_C1 <- gsePathway(
  geneList     = geneList_C1,
  organism     = "human",
  pvalueCutoff = P_CUTOFF_C1,
  maxGSSize    = 3000,
  minGSSize    = 10,
  eps          = 0               # fgseaMultilevel
)
}

## 4. 可視化関数 （_C1 付き）
ridgeplot_C1 <- function(gsea_obj,
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

{
## 5. Ridgeplot 描画
p_ridge_C1 <- ridgeplot_C1(
  gsea_obj          = gsea_react_C1,
  showCategory      = SHOW_TOP_C1,
  pathway_text_size = 4,
  x_breaks          = X_BREAKS_C1
) +
  ggtitle(sprintf("Reactome enrichment (Subgroup C1)", LFC_CUTOFF_C1)) +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 14, margin = margin(b = 6)),
    axis.text.x = element_text(size = 10),
    axis.title.x= element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.title= element_text(size = 10)
  )
print(p_ridge_C1)
}

## 6. ヒット経路をコンソール表示
{
  tbl_short_C1 <- as.data.frame(gsea_react_C1) |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::select(Description, NES, p.adjust) |>
  dplyr::arrange(p.adjust)

sig3 <- function(x, sig = 3){
  vapply(x, function(v){
    if (v == 0) return("0")
    dec <- max(sig - 1 - floor(log10(abs(v))), 0)
    formatC(v, digits = dec, format = "f", drop0trailing = FALSE)
  }, FUN.VALUE = "")
}
tbl_print_C1 <- tbl_short_C1 %>%
  dplyr::mutate(p.adjust = sig3(p.adjust))

print(tbl_print_C1, row.names = FALSE)
}


## GSEA オブジェクト: gsea_react   (Main) / gsea_react_C1  (C1)
tbl_main <- as.data.frame(gsea_react)    |> dplyr::select(ID, NES, p.adjust)
tbl_C1   <- as.data.frame(gsea_react_C1) |> dplyr::select(ID, NES, p.adjust) |>
  dplyr::rename(NES_C1 = NES, p.C1 = p.adjust)

## 1) 有意経路の重なり
sig_main <- tbl_main |> dplyr::filter(p.adjust < 0.05)
sig_C1   <- tbl_C1   |> dplyr::filter(p.C1     < 0.05)

overlap <- dplyr::inner_join(sig_main, sig_C1, by = "ID")
jaccard <- nrow(overlap) /
  (nrow(sig_main) + nrow(sig_C1) - nrow(overlap))

## 2) NES の相関
cor_NES <- cor(overlap$NES, overlap$NES_C1, method = "spearman")

cat(sprintf(
  "共有経路 %d 本 | Jaccard %.2f | NES 相関 ρ = %.3f\n",
  nrow(overlap), jaccard, cor_NES))
# 共有経路 26 本 | Jaccard 0.25 | NES 相関 ρ = 0.966

ggplot(overlap, aes(NES, NES_C1)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  labs(x = "NES (Main)", y = "NES (C1)",
       title = "Pathway NES concordance") +
  theme_bw()


library(clusterProfiler)
comp <- compareCluster(
  geneCluster = list(Main = geneList, C1 = geneList_C1),
  fun         = "gsePathway",
  pvalueCutoff= 0.05,
  organism    = "human"
)
enrichplot::dotplot(comp, showCategory = 20) +
  ggtitle("Reactome GSEA – Main vs Subgroup (C1)")

dotplot(comp,
        showCategory = 20,
        font.size    = 8    # ← デフォルト 12 を小さく
) +
  ggtitle("Reactome GSEA – Main vs Subgroup (C1)") +
  theme(
    axis.text.y  = element_text(size = 7,   # さらに細かく調整したい場合
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