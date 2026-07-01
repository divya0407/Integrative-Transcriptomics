BiocManager::install("GEOquery")
BiocManager::install("affy")
BiocManager::install("limma")
BiocManager::install("oligo") # Recommended package for modern Affymetrix arrays
BiocManager::install("pd.primeview") # Platform Design info for GPL15207
BiocManager::install("primeview.db") # Annotation package for GPL15207

library(GEOquery)
library(affy)
library(oligo)
library(limma)
library(pd.primeview)
library(primeview.db)
setwd("D:/R/1stobjective/GSE113513/New")
# Download the series matrix file (gets sample info, not raw data)
geo_data <- getGEO("GSE113513", GSEMatrix = TRUE, AnnotGPL = FALSE)
# Extract phenodata (sample information)
pdata <- pData(geo_data[[1]])
# Note: For RMA, you need the actual .CEL files downloaded into your working directory.

# Read the CEL files. Ensure you have all 28 CEL files in your working directory.
celfiles <- list.celfiles() # Lists all .CEL files in the current directory
raw_data <- read.celfiles(celfiles)
BiocManager::install("pd.primeview")
list.files(pattern = ".CEL")

library(GEOquery)
library(oligo)
library(pd.primeview) # Required platform data

# 1. Download supplementary files (raw data archive)
# This creates a directory named "GSE113513" in your current working directory
cat("Downloading supplementary files for GSE113513...\n")
getGEOSuppFiles("GSE113513")

# The files are downloaded into the new directory as a compressed archive (e.g., GSE113513_RAW.tar)

# 2. Unpack the downloaded archive
# Set the working directory to the newly created folder
setwd("GSE113513")

# Uncompress the archive file using R's built-in tools
# The archive name is typically "GSE[accession]_RAW.tar"
raw_archive_name <- list.files(pattern = "_RAW.tar")[1]

if (!is.na(raw_archive_name)) {
  cat("Unpacking raw data archive...\n")
  # Use 'untar' to extract the contents
  untar(raw_archive_name)
} else {
  stop("Could not find the raw data archive (GSE*_RAW.tar) in the directory.")
}

# The extracted files (the actual .CEL files) will be gzip compressed (.gz)
# We need to decompress these individual CEL.gz files
cel_gz_files <- list.files(pattern = ".CEL.gz")

cat("Decompressing CEL files...\n")
sapply(cel_gz_files, gunzip)

# 3. Verify files and proceed with RMA
# Now all 28 uncompressed .CEL files should be in your current directory
celfiles <- list.celfiles()

if (length(celfiles) == 28) {
  cat("All 28 CEL files found and decompressed successfully. Proceeding with RMA.\n")
  # Read the CEL files into an AffyBatch object
  raw_data <- read.celfiles(celfiles)
  
  # Perform RMA normalization
  rma_data <- rma(raw_data)
  
  cat("RMA normalization complete.\n")
  # You can now proceed with DEG analysis on the 'rma_data' object
} else {
  stop("Error: Expected 28 CEL files, but found ", length(celfiles), ". Please check the download process.")
}

# script to perform RMA normalization
setwd("D:/R/1stobjective/GSE113513")

library(affy)
library(GEOquery)
library(tidyverse)

# get supplementary files
getGEOSuppFiles("GSE113513")

# untar files
untar("GSE113513/GSE113513_RAW.tar", exdir = 'data/')

# reading in .cel files
raw.data <- ReadAffy(celfile.path = "data/")

# performing RMA normalization
normalized.data <- rma(raw.data)

# get expression estimates
normalized.expr <- as.data.frame(exprs(normalized.data))

# map probe IDs to gene symbols
gse <- getGEO("GSE113513", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE113513_series_matrix.txt@featureData@data
# subset
feature.data <- feature.data[,c(1,11)]

normalized.expr <- normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data, by = 'ID')
write.csv(normalized.expr, "annotated1.csv")

library(limma)
library(tidyverse)
library(GEOquery)

# --- Configuration: Define file path and groups ---
csv_file_path <- "D:/R/1stobjective/GSE113513/annotated1.csv" # Update this path if necessary

# Define sample groups: 14 Normal (N) and 14 Cancer (C) samples
# Based on the experimental design of GSE113513
sample_groups <- factor(c(rep("Normal", 14), rep("Cancer", 14)), levels = c("Normal", "Cancer"))
# --------------------------------------------------

# 1. Load the data you created previously
annotated_data <- read.csv(csv_file_path)

# Separate the expression matrix from the annotation columns
# We assume the last column is 'SYMBOL' (or similar gene name column)
# and all columns starting with "GSM" are samples
expression_data <- annotated_data %>%
  select(starts_with("GSM"), "ID") %>%
  column_to_rownames(var = "ID")

# Get the gene symbols separately (assuming the column name is 'SYMBOL' from your previous code)
# If the column name is different in your CSV, update "SYMBOL" below
gene_symbols <- annotated_data$ID

# Ensure data is a matrix and numeric for limma
expression_matrix <- as.matrix(expression_data)

# 2. Create the design matrix for limma
design <- model.matrix(~sample_groups)
colnames(design) <- c("Normal", "Cancer")

# 3. Fit the linear model
fit <- lmFit(expression_matrix, design)
fit <- eBayes(fit)

# 4. Extract the differential expression results
deg_results <- topTable(fit, coef = "Cancer", adjust.method = "fdr", number = Inf)

# 5. Add the gene symbols back to the results table
deg_results$Symbol <- gene_symbols[match(rownames(deg_results), rownames(expression_matrix))]

# 6. Filter for significant DEGs (e.g., adjusted p-value < 0.05 and |logFC| > 1)
significant_degs <- deg_results %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  filter(!is.na(Symbol) & Symbol != "") %>%
  arrange(adj.P.Val)

# View the top results
head(significant_degs)

# Save results
write.csv(significant_degs, "GSE113513_Significant_DEGs.csv", row.names = TRUE)

# --- Continue from previous script after filtering for significant_degs ---

# 7. Separate DEGs into Up and Downregulated lists

# Upregulated DEGs (logFC > 0)
upregulated_degs <- significant_degs %>%
  filter(logFC > 0) %>%
  arrange(desc(logFC)) # Sort by highest fold change first

# Downregulated DEGs (logFC < 0)
downregulated_degs <- significant_degs %>%
  filter(logFC < 0) %>%
  arrange(logFC) # Sort by lowest fold change first

# 8. Save the dataframes to CSV files

# Total significant DEGs (already created in the previous step, here for completeness)
write.csv(significant_degs, "Total_DEGs.csv", row.names = TRUE)

# Upregulated DEGs
write.csv(upregulated_degs, "Upregulated_DEGs.csv", row.names = TRUE)

# Downregulated DEGs
write.csv(downregulated_degs, "Downregulated_DEGs.csv", row.names = TRUE)

cat("Successfully saved Total_DEGs.csv, Upregulated_DEGs.csv, and Downregulated_DEGs.csv in your working directory.")

results_table_annotated$DEG_status <- "Not Significant"
results_table_annotated$DEG_status[significant_degs$logFC > LOGFC_THRESHOLD & significant_degs$adj.P.Val < FDR_THRESHOLD] <- "Upregulated (Tumor)"
results_table_annotated$DEG_status[significant_degs$logFC < -LOGFC_THRESHOLD & significant_degs$adj.P.Val < FDR_THRESHOLD] <- "Downregulated (Tumor)"
results_table_annotated$DEG_status <- factor(results_table_annotated$DEG_status, levels = c("Upregulated (Tumor)", "Downregulated (Tumor)", "Not Significant"))

# Create the volcano plot 
volcano_plot <- ggplot(results_table_annotated, aes(x = logFC, y = -log10(adj.P.Val), color = DEG_status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated (Tumor)" = "red", 
                                "Downregulated (Tumor)" = "blue", 
                                "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", color = "black") +
  labs(title = paste0("Tumor vs Healthy Colon"),
       x = "log2 Fold Change",
       y = "-log10(Adjusted P-value)",
       color = "DEG Status") +
  theme_minimal() +
  theme(legend.position = "bottom")
#saving a volcano plot
ggsave(filename = "GSE44076_VolcanoPlot_FC2.0_FDR0.05.png",
       plot = volcano_plot,
       width = 8,       # Width in inches
       height = 7,      # Height in inches
       units = "in",
       dpi = 600)       # High resolution for quality

library(ggplot2)
library(dplyr)
library(ggrepel) # Used for labeling top genes without overlap

# Assuming 'deg_results' dataframe is available from previous steps
# and contains 'logFC' (log2 Fold Change) and 'adj.P.Val' (adjusted P-value) columns.

# --- Define cutoffs ---
P_CUTOFF <- 0.05
FC_CUTOFF <- 1.0 # |Log2 Fold Change| > 1
# Calculate the -log10(P-value) threshold
NEG_LOG10_P <- -log10(P_CUTOFF)
# ----------------------

# 1. Prepare the data: Create a new column to categorize each gene
# This is essential for coloring points in ggplot2
deg_results <- deg_results %>%
  mutate(DEG_Status = case_when(
    adj.P.Val < P_CUTOFF & logFC > FC_CUTOFF ~ "Upregulated (Tumor)",
    adj.P.Val < P_CUTOFF & logFC < -FC_CUTOFF ~ "Downregulated (Tumor)",
    TRUE ~ "Not Significant"
  ))

# Ensure the order of factor levels for consistent legend ordering and colors
deg_results$DEG_Status <- factor(deg_results$DEG_Status, 
                                 levels = c("Upregulated (Tumor)", "Downregulated (Tumor)", "Not Significant"))

# 2. Create the ggplot object
volcano_plot <- ggplot(data = deg_results, 
                       aes(x = logFC, y = -log10(adj.P.Val), color = DEG_Status)) +
  
  # Add points
  geom_point(alpha = 0.6, size = 1.5) +
  
  # Add significance threshold lines
  geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = NEG_LOG10_P, col = "gray", linetype = 'dashed') +
  
  # Manually set colors to match the image you provided (Red, Blue, Gray/Black)
  scale_color_manual(values = c("Upregulated (Tumor)" = "red", 
                                "Downregulated (Tumor)" = "blue", 
                                "Not Significant" = "black")) +
  
  # Customize axis labels and plot title
  labs(
    title = "Volcano Plot of GSE113513 DEGs",
    x = expression("log"[2]*" Fold Change"), # Use expression() for math formatting
    y = expression("-log"[10]*"(Adjusted P-value)"),
    color = "DEG Status" # Legend title
  ) +
  
  # Apply a clean theme and center the title
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Display the plot
print(volcano_plot)

ggsave(filename = "GSE113513_VolcanoPlot_FC2.0_FDR0.05.png",
       plot = volcano_plot,
       width = 8,       # Width in inches
       height = 7,      # Height in inches
       units = "in",
       dpi = 600)   

top_n <- 20
top_degs_probes <- rownames(significant_degs[1:top_n, ])
top_expr_matrix <- expr_matrix[top_degs_probes, ]

# Create annotation data for the heatmap
annotation_col <- data.frame(Group = Group_subset)
rownames(annotation_col) <- colnames(top_expr_matrix)

# Use gene symbols for row names in the heatmap
top_degs_symbols <- results_table_annotated[1:top_n, "Symbol"]
rownames(top_expr_matrix) <- top_degs_symbols

# Generate the heatmap 
pheatmap(top_expr_matrix, 
         scale = "row", # Scale the rows (genes) for better visualization of expression patterns
         cluster_cols = FALSE, 
         cluster_rows = TRUE,
         show_rownames = TRUE, 
         show_colnames = FALSE,
         annotation_col = annotation_col,
         main = paste0("Heatmap of Top ", top_n, " DEGs"))
filename = "GSE11353_Heatmap_Final_DEGs.png" # Set your desired filename
width = 8                                    # Width in inches
height = 10                 

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
 # 2. Load your data and filter for significant DEGs
 df <- read.csv("GSE113513_Significant_DEGs.csv")
 sig_genes_df <- subset(df, adj.P.Val < 0.05)
 # Extract the list of significant genes (from the first column of IDs or the Symbol column)
 # Using the 'Symbol' column as before, but knowing they are PROBEIDs now
 probe_ids <- sig_genes_df$Symbol
 
 # 3. Convert Affymetrix PROBEIDs to Entrez IDs
 # We must use 'hgu219.db' for the conversion and 'PROBEID' as the fromType
 entrez_ids_df <- bitr(probe_ids, 
                       fromType = "PROBEID", # Change from SYMBOL to PROBEID
                       toType = "ENTREZID", 
                       OrgDb = hgu219.db) # Use the specific chip DB
 
 # Extract only the list of Entrez IDs
 entrez_ids <- entrez_ids_df$ENTREZID
 
 # 4. Perform GO Enrichment Analysis (using Entrez IDs)
 go_results <- enrichGO(gene = entrez_ids, 
                        OrgDb = org.Hs.eg.db, # Use the general DB for GO terms
                        keyType = "ENTREZID",
                        ont = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE) 
 
 head(go_results)
 
 # 5. Perform KEGG Pathway Enrichment Analysis (using Entrez IDs)
 kegg_results <- enrichKEGG(gene = entrez_ids, 
                            organism = 'hsa', 
                            pvalueCutoff = 0.05)
 
 head(kegg_results)
 
 # 6. Optional: Visualize results
 barplot(go_results, showCategory = 10, title = "GO Biological Process Enrichment")
 dotplot(kegg_results, showCategory = 10, title = "KEGG Pathway Enrichment")
 
 # Save GO results to a CSV file
 write.csv(go_results, file = "GO_enrichment_results.csv", row.names = FALSE)
 
 # Save KEGG results to a CSV file
 write.csv(kegg_results, file = "KEGG_enrichment_results.csv", row.names = FALSE)
 
 # Save GO plot as PNG
 png("GO_barplot.png", width = 800, height = 600, res = 100) # res=100 sets resolution
 barplot(go_results, showCategory = 10, title = "GO Biological Process Enrichment")
 dev.off()
 
 # Save KEGG plot as PNG
 png("KEGG_dotplot.png", width = 1000, height = 700, res = 100)
 dotplot(kegg_results, showCategory = 10, title = "KEGG Pathway Enrichment")
 dev.off()
 
 # Assumed prerequisite: Libraries are loaded and 'entrez_ids' list is ready
 
 # 2. Run GO Enrichment for Molecular Function (MF)
 go_resultsMF <- enrichGO(gene = entrez_ids, 
                        OrgDb = org.Hs.eg.db, # Use the general DB for GO terms
                        keyType = "ENTREZID",
                        ont = "MF", 
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE) 
 
 # 3. Run GO Enrichment for Cellular Component (CC)
 go_resultsCC <- enrichGO(gene = entrez_ids, 
                        OrgDb = org.Hs.eg.db, # Use the general DB for GO terms
                        keyType = "ENTREZID",
                        ont = "CC", 
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE) 
 # View top results for each
 head(go_results_bp)
 head(go_resultsMF)
 head(go_resultsCC)
 
 barplot(go_resultsMF, showCategory = 10, title = "GO Molecular function Enrichment")
 png("MF_barplot.png", width = 800, height = 600, res = 100) # res=100 sets resolution
 barplot(go_resultsMF, showCategory = 10, title = "GO Molecular Function Enrichment")
 dev.off()
 
 barplot(go_resultsCC, showCategory = 10, title = "GO Cellular component Enrichment")
 png("MF_barplot.png", width = 800, height = 600, res = 100) # res=100 sets resolution
 barplot(go_resultsCC, showCategory = 10, title = "GO Cellular component Enrichment")
 dev.off()
 # Save all three GO results to CSV files
 write.csv(go_results_bp, file = "GO_BP_results.csv", row.names = FALSE)
 write.csv(go_resultsMF, file = "GO_MF_results.csv", row.names = FALSE)
 write.csv(go_resultsCC, file = "GO_CC_results.csv", row.names = FALSE)
 
 # --- Prerequisites (Ensure these libraries are loaded) ---
 library(clusterProfiler)
 library(ggplot2)
 library(dplyr) # Required for filtering/data manipulation
 
 # --- 1. Merge the GO Results into a Single Object ---
 cat("Merging GO results for BP, MF, and CC...\n")
 
 # Create a list with the GO objects, using the desired label for each category
 go_results_list <- list(
   # NOTE: Assuming go_results is your BP result
   BP = go_results, 
   MF = go_resultsMF,
   CC = go_resultsCC
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
 
 combined_plot_robust <- ggplot(top_n_terms, aes(x = Description, y = logp, fill = ID)) +
   # Create the bars
   geom_bar(stat = "identity", width = 0.8) +
   # Flip coordinates so the long term names are readable on the Y-axis
   coord_flip() +
   # Facet the plot (crucial for separation into three panels)
   facet_grid(ID ~ ., 
              scales = "free_y", 
              space = "free_y") +
   # Labels and Titles
   labs(x = NULL, # Hide X-axis label after flip
        y = "-log10(Adjusted P-value)",
        title = "Top GO Enrichment Terms (BP, MF, CC)") +
   # Theme adjustments
   theme_bw() +
   theme(legend.position = "none",
         plot.title = element_text(hjust = 0.5),
         axis.text.y = element_text(size = 9),
         strip.text.y = element_text(angle = 0, size = 10)) # Ensure facet labels are readable
 
 
 # Display the plot
 print(combined_plot_robust) 
 
 # Save the plot
 ggsave("GSE113513_Combined_GO_BarPlot_Final.png", plot = combined_plot_robust, width = 12, height = 9, dpi = 600)
 
 cat("\nCombined GO Bar Plot (BP, MF, CC) successfully created and saved as GSE113513_Combined_GO_BarPlot_Final.png.\n")