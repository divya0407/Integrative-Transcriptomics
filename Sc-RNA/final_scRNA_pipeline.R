################################################################################
# GSE299737 — Full scRNA-seq pipleine (Seurat v5 + inferCNV)
################################################################################

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)
library(data.table)
library(infercnv)
library(R.utils)
library(harmony)

# ------------------------------------------------------------------------------
# 0. Global Paths & Parameters
# ------------------------------------------------------------------------------
base_dir <- "C:/simulation/scrna_new"  # Path where your GSM folders are located
out_dir  <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE)

# inferCNV genomic coordinates file (Ensure this file matches your reference assembly, e.g., hg38)
gene_pos_file <- file.path(out_dir, "gene_pos.txt")
if(!file.exists(gene_pos_file)) stop("Please place your gene_pos.txt file at: ", gene_pos_file)

# Parameters mapped strictly to your written methodology
min_features       <- 200
max_features       <- 6000
max_mt             <- 25
n_hvgs             <- 2000
pca_dims           <- 1:20
clust_resolution   <- 0.6
malignant_resolution <- 0.8
logfc_threshold    <- 1.0   # Set to 1.0 to reflect your text: logFC.threshold > 1 (or 0.25 if preferred)

# inferCNV optimization parameters
infercnv_cutoff    <- 0.1
infercnv_threads   <- 4

# ------------------------------------------------------------------------------
# 1. Automated Detection & Archive Extraction of GSM Directories
# ------------------------------------------------------------------------------
gsm_folders <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
gsm_folders <- gsm_folders[basename(gsm_folders) != "results"]

for(folder in gsm_folders){
  gz_files <- list.files(folder, pattern="\\.tar\\.gz$", full.names=TRUE)
  for(f in gz_files){
    gunzip(f, remove=FALSE, overwrite=TRUE)
  }
  tar_files <- list.files(folder, pattern="\\.tar$", full.names=TRUE)
  for(f in tar_files){
    utils::untar(f, exdir=folder)
  }
}

# ------------------------------------------------------------------------------
# 2. Sequential Data Ingestion
# ------------------------------------------------------------------------------
seurat_list <- list()
for(folder in gsm_folders){
  sample_name <- basename(folder)
  message("Reading sample folder: ", sample_name)
  
  # Structural path fallback checker
  if ("filtered_feature_bc_matrix" %in% list.files(folder)) {
    data_path <- file.path(folder, "filtered_feature_bc_matrix")
  } else {
    data_path <- folder
  }
  
  success <- FALSE
  try({
    data <- Read10X(data.dir = data_path)
    so <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)
    so$sample <- sample_name
    seurat_list[[sample_name]] <- so
    success <- TRUE
  }, silent=TRUE)
  
  if(!success){
    # Fallback ingestion using explicit ReadMtx matrices
    mtx   <- list.files(folder, pattern="matrix.*mtx", full.names=TRUE, recursive=TRUE)
    genes <- list.files(folder, pattern="genes.*tsv|features.*tsv", full.names=TRUE, recursive=TRUE)
    bar   <- list.files(folder, pattern="barcodes.*tsv", full.names=TRUE, recursive=TRUE)
    if(length(mtx) > 0 && length(genes) > 0 && length(bar) > 0){
      data2 <- ReadMtx(mtx=mtx[1], features=genes[1], cells=bar[1])
      so <- CreateSeuratObject(counts = data2, project = sample_name, min.cells = 3, min.features = 200)
      so$sample <- sample_name
      seurat_list[[sample_name]] <- so
    } else { 
      warning("Skipping directory. Matrix structure unreadable: ", folder) 
    }
  }
}

# ------------------------------------------------------------------------------
# 3. Vectorized Seurat Dataset Merging
# ------------------------------------------------------------------------------
combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list),
  project = "CRC_11_samples"
)
saveRDS(combined, file=file.path(out_dir, "merged_raw.rds"))

# ------------------------------------------------------------------------------
# 4. Strict Quality Control Filtering
# ------------------------------------------------------------------------------
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# High-resolution QC plot storage
png(file.path(out_dir, "qc_violin_plots.png"), width=1600, height=1200, res=150)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Subset data based on your exact manuscript thresholds
combined <- subset(combined, 
                   subset = nFeature_RNA > min_features & 
                     nFeature_RNA < max_features & 
                     percent.mt < max_mt)

# ------------------------------------------------------------------------------
# 5. Normalization, Feature Selection, Scaling, & Dimensionality Reduction
# ------------------------------------------------------------------------------
combined <- NormalizeData(combined, normalization.method = "LogNormalize")
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = n_hvgs)
combined <- ScaleData(combined, vars.to.regress = "percent.mt")
combined <- RunPCA(combined, features = VariableFeatures(combined))

# ------------------------------------------------------------------------------
# 6. Harmony Multi-Sample Integration & Unsupervised Clustering
# ------------------------------------------------------------------------------
combined <- RunHarmony(combined, group.by.vars = "sample")

combined <- FindNeighbors(combined, reduction = "harmony", dims = pca_dims)
combined <- FindClusters(combined, resolution = clust_resolution)
combined <- RunUMAP(combined, reduction = "harmony", dims = pca_dims)
combined <- RunTSNE(combined, reduction = "harmony", dims = pca_dims)

# Save global visual maps
png(file.path(out_dir, "umap_global_clusters.png"), width=1600, height=1200, res=150)
DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

# ------------------------------------------------------------------------------
# 7. Layer Unification & Global Marker Identification
# ------------------------------------------------------------------------------
combined <- JoinLayers(combined)

markers_all <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = logfc_threshold)
fwrite(markers_all, file = file.path(out_dir, "All_Clusters_Markers.csv"))

# Extract top markers per cluster and write to terminal panel
top10_all <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
print(top10_all, n = 200)
write.csv(top10_all, file.path(out_dir, "Top10_Markers_Per_Cluster.csv"), row.names = FALSE)

# Generate evaluation figures to verify lineage boundaries
png(file.path(out_dir, "canonical_expression_featureplot.png"), width=2400, height=1800, res=200)
FeaturePlot(combined, features = c("EPCAM", "KRT8", "CD3D", "CD8A", "MS4A1", "LYZ", "COL1A1", "PECAM1"), reduction="umap", ncol=4)
dev.off()

png(file.path(out_dir, "canonical_expression_dotplot.png"), width=1600, height=1000, res=150)
DotPlot(combined, features = c("EPCAM","KRT8","KRT18","CD3D","CD8A","MS4A1","LYZ","CD14","NKG7","GNLY","COL1A1","PECAM1")) + RotatedAxis()
dev.off()

# ------------------------------------------------------------------------------
# 8. Cell Type Resolution Mapping 
# ------------------------------------------------------------------------------
new_labels <- c(
  "Malignant", "Malignant", "Malignant", "Fibroblast", "Myeloid", 
  "Endothelial", "NK", "T cells", "T cells", "Myeloid", 
  "NK", "B cells", "T cells", "Fibroblast", "Malignant", 
  "Plasma", "Fibroblast", "Malignant", "Malignant", "Malignant", "Malignant"
)
names(new_labels) <- levels(combined)
combined <- RenameIdents(combined, new_labels)
combined$celltype <- Idents(combined)

# Plot cleanly annotated groupings
png(file.path(out_dir, "umap_annotated_celltypes.png"), width=1600, height=1200, res=150)
DimPlot(combined, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 0.5)
dev.off()

saveRDS(combined, file=file.path(out_dir, "merged_annotated.rds"))

# ------------------------------------------------------------------------------
# 9. Differential Expression: Tumor vs. Background
# ------------------------------------------------------------------------------
tumor_vs_all <- FindMarkers(combined, ident.1 = "Malignant", only.pos = TRUE, min.pct = 0.25, logfc.threshold = logfc_threshold)
write.csv(tumor_vs_all, file.path(out_dir, "Malignant_vs_All_Markers.csv"), row.names = TRUE)

# ------------------------------------------------------------------------------
# 10. Isolated Sub-Clustering of the Malignant Epithelial Population
# ------------------------------------------------------------------------------
malignant <- subset(combined, subset = celltype == "Malignant")

malignant <- NormalizeData(malignant, normalization.method = "LogNormalize")
malignant <- FindVariableFeatures(malignant, selection.method = "vst", nfeatures = n_hvgs)
malignant <- ScaleData(malignant, vars.to.regress = "percent.mt")
malignant <- RunPCA(malignant, features = VariableFeatures(malignant))

# Correct patient-specific tumor biases using Harmony
malignant <- RunHarmony(malignant, group.by.vars = "sample")

# Re-clustering malignant space using high-resolution parameter (0.8)
malignant <- FindNeighbors(malignant, reduction = "harmony", dims = pca_dims)
malignant <- FindClusters(malignant, resolution = malignant_resolution)
malignant <- RunUMAP(malignant, reduction = "harmony", dims = pca_dims)
malignant <- RunTSNE(malignant, reduction = "harmony", dims = pca_dims)

png(file.path(out_dir, "tsne_malignant_subclusters.png"), width=1600, height=1200, res=150)
DimPlot(malignant, reduction = "tsne", label = TRUE, pt.size = 0.5)
dev.off()

malignant_markers <- FindAllMarkers(malignant, only.pos = TRUE, min.pct = 0.25, logfc.threshold = logfc_threshold)
write.csv(malignant_markers, file.path(out_dir, "Malignant_Subclusters_Markers.csv"), row.names = FALSE)

saveRDS(malignant, file=file.path(out_dir, "malignant_processed.rds"))

# ------------------------------------------------------------------------------
# 11. Large-Scale Clonal Chromosomal Copy Number Variation Parsing (inferCNV)
# ------------------------------------------------------------------------------
# Construct cell identification manifest
annot_df <- data.frame(cell = colnames(combined), group = combined$celltype, stringsAsFactors = FALSE)
write.table(annot_df, file = file.path(out_dir, "cell_annotation.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Isolate background reference clusters (Every lineage excluding Malignant cells)
ref_groups <- setdiff(unique(annot_df$group), "Malignant")

# Extract counts as a sparse matrix directly to preserve RAM/Memory allocation limits
counts_matrix <- GetAssayData(combined, assay = "RNA", slot = "counts")

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file  = file.path(out_dir, "cell_annotation.txt"),
  gene_order_file   = gene_pos_file,
  ref_group_names   = ref_groups
)

infercnv_outdir <- file.path(out_dir, "infercnv")
dir.create(infercnv_outdir, showWarnings = FALSE)

# Execute inferCNV modeling 
infercnv_result <- infercnv::run(
  infercnv_obj,
  cutoff            = infercnv_cutoff,
  out_dir           = infercnv_outdir,
  cluster_by_groups = TRUE,
  denoise           = TRUE,
  HMM               = TRUE,
  num_threads       = infercnv_threads
)

message("--- Execution complete. All publication assets safely archived inside: ", out_dir, " ---")