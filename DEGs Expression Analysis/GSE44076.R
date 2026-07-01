# 1. Load Libraries and Data
library(GEOquery)
library(limma)
library(hgu133plus2.db)
library(ggplot2)
library(pheatmap)
library(dplyr)

# Download GSE44076
gse <- getGEO("GSE44076", GSEMatrix = TRUE)[[1]]

# 2. Pre-processing and Normalization
expr <- exprs(gse)
pheno <- pData(gse)

# Replace 0s and Log2 Transform (Affy data is usually pre-logged, but we check)
if(max(expr) > 100) {
  expr[expr <= 0] <- 1
  expr <- log2(expr)
}

# Quantile Normalization
expr_norm <- normalizeBetweenArrays(expr, method = "quantile")

# 3. Grouping and Experimental Design
# Define groups based on the characteristics column
group <- factor(ifelse(
  grepl("tumor|cancer|adenocarcinoma", pheno$characteristics_ch1, ignore.case = TRUE),
  "Tumor", "Normal"
), levels = c("Normal", "Tumor"))

design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Tumor_vs_Normal")

# 4. Differential Expression Analysis (limma)
fit <- lmFit(expr_norm, design)
fit <- eBayes(fit)
deg <- topTable(fit, coef = "Tumor_vs_Normal", adjust.method = "BH", number = Inf)

# 5. Annotation
deg$GeneSymbol <- mapIds(hgu133plus2.db, 
                         keys = rownames(deg), 
                         column = "SYMBOL", 
                         keytype = "PROBEID", 
                         multiVals = "first")

# Create a clean version for plotting (removing NAs)
deg_clean <- deg[!is.na(deg$GeneSymbol), ]

# Define Thresholds
logFC_threshold <- 1
adjP_threshold <- 0.05
upregulated <- deg %>% filter(logFC > logFC_threshold & adj.P.Val < adjP_threshold & !is.na(GeneSymbol))
downregulated <- deg %>% filter(logFC < -logFC_threshold & adj.P.Val < adjP_threshold & !is.na(GeneSymbol))

write.csv(upregulated, "upregulated.csv", row.names = TRUE)
write.csv(downregulated, "downregulated.csv", row.names = TRUE)
# 6. Volcano Plot
deg$DEG_status <- "Not Significant"
deg$DEG_status[deg$logFC > logFC_threshold & deg$adj.P.Val < adjP_threshold] <- "Upregulated"
deg$DEG_status[deg$logFC < -logFC_threshold & deg$adj.P.Val < adjP_threshold] <- "Downregulated"



volcano <- ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val), color = DEG_status)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "Black")) +
  theme_minimal() +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(adjP_threshold), linetype = "dashed") +
  labs(title = "GSE44076: Tumor vs Normal", 
       x = "log2 Fold Change", y = "-log10 Adjusted P-value")

print(volcano)
ggsave("volcano_plot.png", plot = volcano, width = 8, height = , dpi = 600)

# 7. Heatmap of Top 20 DEGs
# Extract top 20 genes by absolute Fold Change that have a Gene Symbol
top20_probes <- rownames(head(deg[!is.na(deg$GeneSymbol),][order(deg[!is.na(deg$GeneSymbol),]$adj.P.Val), ], 20))
heatmap_matrix <- expr_norm[top20_probes, ]
rownames(heatmap_matrix) <- make.unique(deg[top20_probes, "GeneSymbol"])

ann_col <- data.frame(Group = group)
rownames(ann_col) <- colnames(heatmap_matrix)



# Saving Heatmap to PNG directly
png("heatmap_2.png", width = 800, height = 1000, res = 120)
pheatmap(heatmap_matrix, 
         scale = "row", 
         annotation_col = ann_col,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 20 DEGs",
         show_colnames = FALSE)
dev.off()

# 8. Export Results
write.csv(deg, "GSE44076_Complete_Results.csv")

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
genes_to_test <- unique(deg$GeneSymbol)

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
barplot(ego, showCategory = 10, title = "GO Biological Process Enrichment")
write.csv(ego, file = "BP_GO_enrichment_results.csv", row.names = FALSE)
png("BP_GO_barplot.png", width = 800, height = 600, res = 100) # res=100 sets resolution
barplot(ego, showCategory = 10, title = "GO Biological Process Enrichment")
dev.off()
getwd()
setwd("D:/R/1stobjective/GSE44076")

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
ggsave("GSE44076_Combined_GO_BarPlot_Final.png", plot = p_combined, width = 12, height = 9, dpi = 600)


