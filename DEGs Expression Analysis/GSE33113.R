library(GEOquery)
setwd("D:/R/1stobjective/GSE33113/GSE33113")

gse <- getGEO("GSE33113", GSEMatrix = TRUE)
gse <- gse[[1]]

expr <- exprs(gse)
pheno <- pData(gse)
summary(as.numeric(expr))
expr[expr == 0] <- NA
expr <- expr + 1
expr_log2 <- log2(expr)
summary(as.numeric(expr_log2))
keep <- rowSums(is.na(expr_log2)) < ncol(expr_log2) * 0.5
expr_log2 <- expr_log2[keep, ]
library(limma)

expr_norm <- normalizeBetweenArrays(expr_log2, method = "quantile")
par(mfrow = c(1,2))
boxplot(expr_log2, main = "Log2 (Before norm)", outline = FALSE, col = "skyblue", las = 2) # las = 2 rotates sample names for readability)
boxplot(expr_norm, main = "After quantile norm", outline = FALSE, col = "salmon", las = 2)
png("GSE33113_Normalization_Boxplots.png", width = 10, height = 6, units = "in", res = 300)
dev.off()
pheno <- pData(gse)

table(pheno$characteristics_ch1)

group <- ifelse(
  grepl("tumor|cancer", pheno$characteristics_ch1, ignore.case = TRUE),
  "Tumor", "Normal"
)

group <- factor(group)
table(group)
colnames(pheno)
head(pheno)
unique(pheno$title)
group <- ifelse(
  grepl("^Normal mucosa", pheno$title),
  "Normal",
  "Tumor"
)

group <- factor(group, levels = c("Normal", "Tumor"))
table(group)
library(limma)

design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Tumor_vs_Normal")
fit <- lmFit(expr_norm, design)
fit <- eBayes(fit)
deg <- topTable(
  fit,
  coef = "Tumor_vs_Normal",
  adjust.method = "BH",
  number = Inf
)

head(deg)

BiocManager::install("hgu133plus2.db", ask = FALSE, update = FALSE)
library(hgu133plus2.db)

deg$GeneSymbol <- mapIds(
  hgu133plus2.db,
  keys = rownames(deg),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

deg_sig <- subset(
  deg,
  adj.P.Val < 0.05 & abs(logFC) > 1
)

dim(deg_sig)
dim (deg)
write.csv(deg, "GSE33113_all_DEGs.csv")
write.csv(deg_sig, "GSE33113_significant_DEGs.csv")

deg_sig <- subset(
  deg,
  !is.na(GeneSymbol) &
    adj.P.Val < 0.05 &
    abs(logFC) > 1
)
deg_up <- subset(deg_sig, logFC > 1)
deg_down <- subset(deg_sig, logFC < -1)
cat("Upregulated genes:", nrow(deg_up), "\n")
cat("Downregulated genes:", nrow(deg_down), "\n")
write.csv(deg_up, "GSE33113_Upregulated_DEGs.csv", row.names = FALSE)
write.csv(deg_down, "GSE33113_Downregulated_DEGs.csv", row.names = FALSE)

range(exprs(gse))
range(exprs(gse), na.rm = TRUE)

#OPTIONAL: Save Only Gene Symbols (For Enrichment)
write.table(
  unique(deg_up$GeneSymbol),
  "GSE33113_Upregulated_Genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  unique(deg_down$GeneSymbol),
  "GSE33113_Downregulated_Genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
#PCA
sum(is.na(expr_norm))
sum(is.infinite(expr_norm))

keep <- complete.cases(expr_norm)
expr_pca <- expr_norm[keep, ]

sum(is.na(expr_pca))
sum(is.infinite(expr_pca))
pca <- prcomp(t(expr_pca), scale. = TRUE)

plot(
  pca$x[,1], pca$x[,2],
  col = as.numeric(group),
  pch = 16,
  xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
  ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
  main = "PCA of GSE33113"
)

legend(
  "topright",
  legend = levels(group),
  col = 1:2,
  pch = 16
)

#VOLCANO PLOT
deg_plot <- deg
deg_plot <- deg_plot[!is.na(deg_plot$GeneSymbol), ]
with(deg_plot, {
  plot(
    logFC,
    -log10(adj.P.Val),
    pch = 20,
    main = "Volcano Plot: Tumor vs Normal (GSE33113)",
    xlab = "Log2 Fold Change",
    ylab = "-log10 Adjusted P-value",
    col = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                 ifelse(logFC > 1, "red", "blue"),
                 "grey")
  )
  
  abline(v = c(-1, 1), lty = 2)
  abline(h = -log10(0.05), lty = 2)
})

library(ggplot2)

volcano_df <- deg
volcano_df <- volcano_df[!is.na(volcano_df$GeneSymbol), ]

volcano_df$DEG_status <- "Not Significant"
volcano_df$DEG_status[
  volcano_df$adj.P.Val < 0.05 & volcano_df$logFC > 1
] <- "Upregulated (Tumor)"

volcano_df$DEG_status[
  volcano_df$adj.P.Val < 0.05 & volcano_df$logFC < -1
] <- "Downregulated (Tumor)"

volcano_df$DEG_status <- factor(
  volcano_df$DEG_status,
  levels = c("Upregulated (Tumor)", "Downregulated (Tumor)", "Not Significant")
)

ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), color = DEG_status)) +
  geom_point(alpha = 0.8, size = 1.8) +
  
  scale_color_manual(
    values = c(
      "Upregulated (Tumor)" = "red",
      "Downregulated (Tumor)" = "blue",
      "Not Significant" = "black"
    )
  ) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  
  labs(
    title = "Volcano Plot of GSE33113 DEGs",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~P-value),
    color = "DEG Status"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

volcano_df <- deg

volcano_df <- volcano_df[
  !is.na(volcano_df$logFC) &
    !is.na(volcano_df$adj.P.Val) &
    volcano_df$adj.P.Val > 0,
]
volcano_df$DEG_status <- "Not Significant"

volcano_df$DEG_status[
  volcano_df$adj.P.Val < 0.05 & volcano_df$logFC > 1
] <- "Upregulated (Tumor)"

volcano_df$DEG_status[
  volcano_df$adj.P.Val < 0.05 & volcano_df$logFC < -1
] <- "Downregulated (Tumor)"

volcano_df$DEG_status <- factor(
  volcano_df$DEG_status,
  levels = c("Upregulated (Tumor)", "Downregulated (Tumor)", "Not Significant")
)

volcano_plot <- ggplot(volcano_df, aes(logFC, -log10(adj.P.Val), color = DEG_status)) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(
    values = c(
      "Upregulated (Tumor)" = "red",
      "Downregulated (Tumor)" = "blue",
      "Not Significant" = "black"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(
    title = "Volcano Plot of GSE4107 DEGs",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~P-value),
    color = "DEG Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(filename = "GSE4107_VolcanoPlot_FC2.0_FDR0.05.png",
       plot = volcano_plot,
       width = 8,       # Width in inches
       height = 7,      # Height in inches
       units = "in",
       dpi = 600)   

library(pheatmap)
library(RColorBrewer)

# 1. Select the top 20 most significant genes (lowest adj.P.Val)
# This ensures we get the most robust markers for both groups
top20_sig <- head(deg_sig[order(deg_sig$adj.P.Val), ], 20)
top20_genes_expr <- expr_norm[rownames(top20_sig), ]

# 2. Replace Probes with Gene Symbols for the display
rownames(top20_genes_expr) <- top20_sig$GeneSymbol

# 3. Create the Annotation Data Frame (matches the 'Group' bar at the top)
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(top20_genes_expr)

# 4. Define specific colors to match your reference image
ann_colors <- list(
  Group = c(Normal = "#00FFFF", Tumor = "#FF8080") # Cyan and Light Coral
)

# 5. Generate the Heatmap
heatmap_plot <- pheatmap(
  top20_genes_expr, 
  scale = "row",                      # Z-score scaling (essential for the blue-red contrast)
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  cluster_rows = TRUE,                # Cluster genes to show patterns
  cluster_cols = TRUE,                # Cluster samples to group 'Healthy' vs 'Adenocarcinoma'
  annotation_col = annotation_col,    # Add the color bar at the top
  annotation_colors = ann_colors,     # Use the specific colors defined above
  color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100), # Blue-Yellow-Red palette
  show_colnames = FALSE,              # Hide individual sample IDs for a cleaner look
  fontsize_row = 10,
  main = "Heatmap of Top 20 DEGs"
)

# 6. Save with high resolution
ggsave("Heatmap_Refined_Final.png", plot = heatmap_plot, width = 8, height = 7, dpi = 600)

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
library(clusterProfiler)
library(org.Hs.eg.db) # Annotation database for Human Entrez IDs

# --- Prerequisites ---
library(clusterProfiler)
library(org.Hs.eg.db) # Annotation database for Human Entrez IDs
BiocManager::install("AnnotationDbi")

# 1. Load libraries
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db) # For human genes; change to org.Mm.eg.db for mouse, etc.
library(hgu219.db)
# Extract gene symbols from your significant DEGs
genes_to_test <- unique(deg_sig$GeneSymbol)

# Convert Symbols to Entrez IDs for KEGG
gene_conv <- bitr(genes_to_test, fromType = "SYMBOL", 
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)

entrez_genes <- gene_conv$ENTREZID


# 4. Perform GO Enrichment Analysis (using Entrez IDs)
# Run GO Enrichment
ego <- enrichGO(gene          = genes_to_test,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",  # Performs BP, CC, and MF simultaneously
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

head(ego)

# Visualize GO Results
barplot(ego, showCategory = 10, title = "GO Biological Process Enrichment")
write.csv(ego, file = "BP_GO_enrichment_results.csv", row.names = FALSE)
# Save GO plot as PNG
png("BP_GO_barplot.png", width = 800, height = 600, res = 100) # res=100 sets resolution
barplot(ego, showCategory = 10, title = "GO Biological Process Enrichment")
dev.off()

# 2. Run GO Enrichment for Molecular Function (MF)
ego_MF <- enrichGO(gene          = genes_to_test,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "MF",  # Performs BP, CC, and MF simultaneously
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

head(ego_MF)
# Visualize GO Results
barplot(ego_MF, showCategory = 10, title = "GO Molecular Function Enrichment")
write.csv(ego_MF, file = "MF_GO_enrichment_results.csv", row.names = FALSE)
# Save GO plot as PNG
png("MF_GO_barplot.png", width = 800, height = 600, res = 100) # res=100 sets resolution
barplot(ego_MF, showCategory = 10, title = "GO Molecular Function Enrichment")
dev.off()

# 3. Run GO Enrichment for Cellular Component (CC)
ego_CC <- enrichGO(gene          = genes_to_test,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "SYMBOL",
                   ont           = "CC",  # Performs BP, CC, and MF simultaneously
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

head(ego_CC)
# Visualize GO Results
barplot(ego_CC, showCategory = 10, title = "GO Cellular Component Enrichment")
write.csv(ego_CC, file = "CC_GO_enrichment_results.csv", row.names = FALSE)
# Save GO plot as PNG
png("CC_GO_barplot.png", width = 800, height = 600, res = 100) # res=100 sets resolution
barplot(ego_CC, showCategory = 10, title = "GO Cellular Component Enrichment")
dev.off()

# --- 1. Merge the GO Results into a Single Object ---
cat("Merging GO results for BP, MF, and CC...\n")

# Create a list with the GO objects, using the desired label for each category
go_results_list <- list(
  # NOTE: Assuming go_results is your BP result
  BP = ego, 
  MF = ego_MF,
  CC = ego_CC
)

# Merge the results, adding the 'ONTOLOGY' column for separation
go_results_all <- merge_result(go_results_list)

# --- 2. Filter and Prepare Data for Robust Plotting ---

# Convert the merged object to a data frame
go_df_all <- as.data.frame(go_results_all)

# Calculate -log10(p.adjust) for bar height
go_df_all$logp <- -log10(go_df_all$p.adjust)

# Filter to get the top 8-10 terms per category for a clean plot
top_n_terms <- go_df_all %>%
  group_by(ID) %>% 
  arrange(p.adjust) %>% # Sort by significance
  slice_head(n = 8) %>% # Take the top 8 most significant terms per category
  ungroup()

# Reorder the Term factor for plotting (orders bars by significance within each facet)
top_n_terms$Description <- factor(top_n_terms$Description, 
                                  levels = top_n_terms$Description[order(top_n_terms$logp, decreasing = TRUE)])

# --- 3. Create the Combined Bar Plot using ggplot2 ---

cat("Generating combined bar plot using ggplot2...\n")

# Run KEGG Enrichment
ekegg <- enrichKEGG(gene         = entrez_genes,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)
head(ekegg)
dotplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment")
write.csv(ekegg, file = "KEGG_enrichment_results.csv", row.names = FALSE)
# Save KEGG plot as PNG
png("KEGG_dotplot.png", width = 1000, height = 700, res = 100)
dotplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment")
dev.off()

# Visualize KEGG Results
barplot(ekegg, showCategory = 20, title = "KEGG Pathway Enrichment")

# 1. Run GO for all categories (BP, CC, MF) simultaneously
ego_all <- enrichGO(gene          = genes_to_test,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "SYMBOL",
                    ont           = "ALL",  # Captures all three categories
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    readable      = TRUE)

# 2. Generate the faceted barplot
library(ggplot2)
p_combined <- barplot(ego_all, showCategory = 5, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(strip.text = element_text(face = "bold", size = 10))
ggsave("GSE33113_Combined_GO_BarPlot_Final.png", plot = p_combined, width = 12, height = 9, dpi = 600)
