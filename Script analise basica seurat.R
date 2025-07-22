#==============================================================================#
#            SCRIPT FOR INTEGRATIVE ANALYSIS OF scRNA-seq DATA             #
#                                 with Seurat
#
#==============================================================================#

# --- 1. LOAD PACKAGES ---
# Ensure all necessary packages are installed.
# install.packages(c("Seurat", "tidyverse", "patchwork", "harmony", "SeuratWrappers"))

library(Seurat)
library(tidyverse) # For data manipulation (e.g., %>% and read_tsv)
library(patchwork) # For combining plots
library(harmony)   # For the Harmony integration method
library(SeuratWrappers) # For integration with Seurat v5


#==============================================================================#
# --- 2. AUTOMATED LOADING AND QUALITY CONTROL (QC) PER SAMPLE ---
#==============================================================================#

# This section automates the loading of multiple samples (10x Genomics format)
# and the application of quality control (QC) filters for each of them.
# --- 2.1. Set the working directory ---
# Enter the path to the folder containing your sample files here.
# e.g., setwd("/home/user/my_project/raw_data")
setwd("/path/to/your/data/folder")

# Create subdirectories to save QC results and RDS objects
dir.create("QC", showWarnings = FALSE)
dir.create("RDS", showWarnings = FALSE)


# --- 2.2.
Helper function to identify outliers ---
# This function calculates outliers based on the Median Absolute Deviation (MAD),
# which is more robust to extreme values than the standard deviation.
mad_outlier <- function(seurat_object, metric, nmads) {
  M <- seurat_object@meta.data[[metric]]
  median_M <- median(M, na.rm = TRUE)
  mad_M <- mad(M, na.rm = TRUE)
  
  # Identifies cells that are 'nmads' deviations below or above the median
  outlier <- (M < (median_M - nmads * mad_M)) |
(M > (median_M + nmads * mad_M))
  return(outlier)
}

# --- 2.3.
Main QC function ---
# This function iterates over each sample, performs QC, and saves the results.
# NOTE ON MITOCHONDRIAL GENES:
# By default, Seurat looks for the "mt-" prefix (lowercase) for mouse.
# If your data uses "MT-" (human) or has no prefix, adjust the `pattern` below.
# For mouse without a prefix, an explicit list of genes can be used, such as:
# pattern_mt <- "^Nd1|^Nd2|^Co1|^Co2|^Atp8|^Atp6|^Co3|^Nd3|^Nd4l|^Nd4|^Nd5|^Nd6|^Cytb"
# For human, the pattern is: pattern = "^MT-"

process_samples_qc <- function(samples_prefix, mt_pattern = "^mt-", perc_mt_threshold = 10) {
  
  for (sample_id in samples_prefix) {
    message(paste("Processing sample:", sample_id))
    
    # Load data from 10x format (matrix, barcodes, features)
    counts <- ReadMtx(
      mtx = paste0(sample_id, "_matrix.mtx.gz"),
      cells = paste0(sample_id, "_barcodes.tsv.gz"),
      features = paste0(sample_id, "_features.tsv.gz")
    )
    
# Create the Seurat object
    # min.cells: Includes features (genes) expressed in at least 3 cells.
# min.features: Includes cells that express at least 300 features.
X <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 300)
    
    # Save the number of cells before filtering
    cells_before <- ncol(X)
    
    # Calculate the percentage of mitochondrial genes
    X[["percent.mt"]] <- PercentageFeatureSet(X, pattern = mt_pattern)
    
    # Outlier-based filters
    # Removes cells with highly discrepant RNA counts or number of genes (5 MADs from the median)
    outlier_filter <- mad_outlier(X, 'nCount_RNA', 5) |
mad_outlier(X, 'nFeature_RNA', 5)
    X <- subset(X, cells = colnames(X)[!outlier_filter])
    
    # Mitochondrial percentage filter
    # NOTE: For single-nucleus (snRNA-seq), the threshold is much lower (e.g., 1% or 2%).
# Change the `perc_mt_threshold` if necessary.
    X <- subset(X, subset = percent.mt < perc_mt_threshold)
    
    cells_after <- ncol(X)
    
    # Save a QC summary
    qc_summary <- data.frame(
      Sample_ID = sample_id,
      Cells_Before_QC = cells_before,
      Cells_After_QC = cells_after
    )
    write.csv(qc_summary, file.path("QC", paste0(sample_id, "_QC_Summary.csv")), row.names = FALSE)
    
    # Save the filtered Seurat object
    saveRDS(X, file.path("RDS", paste0(sample_id, 
"_AfterQC.rds")))
    
    message(paste("Sample", sample_id, "finished. Remaining cells:", cells_after))
  }
}

# --- 2.4.
Run the QC pipeline ---
# Lists all files in the folder ending with "_barcodes.tsv.gz"
# and extracts the prefix for each sample.
file_list <- list.files(pattern = "_barcodes.tsv.gz$")
sample_prefixes <- sub("_barcodes.tsv.gz$", "", file_list)

# Execute the processing function
# Change the mitochondrial pattern or threshold here, if necessary
process_samples_qc(sample_prefixes, mt_pattern = "^mt-", perc_mt_threshold = 10)


#==============================================================================#
# --- 3. QC VISUALIZATION (OPTIONAL, BUT RECOMMENDED) ---
#==============================================================================#

# Load the saved RDS objects after QC for visualization
rds_files <- list.files("RDS", pattern = "\\.rds$", full.names = TRUE)
seurat_list_qc <- lapply(rds_files, readRDS)
names(seurat_list_qc) <- sub("_AfterQC.rds", "", basename(rds_files))

# Function to create standardized Violin Plots
VlnPlot_Custom <- function(seurat_object, feature, ...) {
  VlnPlot(seurat_object, feature, pt.size = 0, ...) + # pt.size = 0 removes the points
    scale_y_continuous(limits = c(0, 20)) + # Standardizes the Y-axis
    theme(legend.position = "none", axis.title.x = element_blank()) +
    labs(title = seurat_object@project.name)
}

# Create and combine `percent.mt` plots for all samples
plot_list_mt <- lapply(seurat_list_qc, function(obj) VlnPlot_Custom(obj, "percent.mt"))
combined_plots <- wrap_plots(plotlist = plot_list_mt, ncol = 3) # Adjust ncol according to the number of samples

# Show the combined plot
print(combined_plots)
# ggsave("QC_ViolinPlots_PercentMT.png", combined_plots, width = 12, height = 8)


#==============================================================================#
# --- 4. DATA INTEGRATION (SEURAT v5) ---
#==============================================================================#

# This is the modern and recommended approach in Seurat v5, using `IntegrateLayers`.
# --- 4.1. Merge objects and prepare for integration ---
# If the objects are not in memory, load them
# rds_files <- list.files("RDS", pattern = "\\.rds$", full.names = TRUE)
# seurat_list_qc <- lapply(rds_files, readRDS)

# Merge all objects from the list into a single Seurat object
merged_seurat <- merge(seurat_list_qc[[1]], y = seurat_list_qc[-1])

# Split the 'RNA' assay by sample ('orig.ident').
Essential for `IntegrateLayers`.
merged_seurat[['RNA']] <- split(merged_seurat[["RNA"]], f = merged_seurat$orig.ident)

# --- 4.2.
Standard preprocessing BEFORE integration ---
# This pipeline is run on each "layer" (sample) individually.
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)


# --- 4.3. Run the Integration ---
# Choose ONE of the methods below.
RPCA is usually faster and more robust.

# Increase the global memory limit for processing, if necessary
options(future.globals.maxSize = 8000 * 1024^2) # 8 GB

# METHOD A: RPCA (Reciprocal PCA) - Fast and robust
integrated_seurat <- IntegrateLayers(
  object = merged_seurat,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = TRUE
)

# METHOD B: CCA (Canonical Correlation Analysis)
# integrated_seurat <- IntegrateLayers(
#   object = merged_seurat,
#   method = CCAIntegration,
#   orig.reduction = "pca",
#   new.reduction = "integrated.cca",
#   verbose = TRUE
# )

# METHOD C: Harmony
# integrated_seurat <- IntegrateLayers(
#   
object = merged_seurat,
#   method = HarmonyIntegration,
#   orig.reduction = "pca",
#   new.reduction = "harmony",
#   verbose = TRUE
# )

# Re-join the RNA assay layers after integration for future analyses
integrated_seurat[["RNA"]] <- JoinLayers(integrated_seurat[["RNA"]])

# Save the integrated object
# saveRDS(integrated_seurat, "integrated_seurat_rpca.rds")


#==============================================================================#
# --- 5. DOWNSTREAM ANALYSIS (CLUSTERING AND VISUALIZATION) ---
#==============================================================================#
# integrated_seurat <- readRDS("integrated_seurat_rpca.rds")

# --- 5.1.
Principal Component Analysis (PCA) and Clustering ---
# Use the "reduction" generated by your integration (e.g., "integrated.rpca")

# The ElbowPlot helps determine the number of dimensions (principal components) to use
ElbowPlot(integrated_seurat, ndims = 30, reduction = "integrated.rpca")

# Define the number of dimensions based on the "elbow" of the plot
# Example: dims_to_use <- 1:20
dims_to_use <- 1:20 
reduction_name <- "integrated.rpca" # Change to "integrated.cca" or "harmony" if you used another method

# Find neighbors and clusters
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = reduction_name, dims = dims_to_use)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5) # Adjust the resolution for more/fewer clusters

# --- 5.2.
Visualization with UMAP ---
integrated_seurat <- RunUMAP(integrated_seurat, dims = dims_to_use, reduction = reduction_name)

# Plot UMAP by cluster
p1 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Clusters")

# Plot UMAP by original sample (to check integration quality)
p2 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Original Sample")

# Display plots side-by-side
p1 + p2


#==============================================================================#
# --- 6. IDENTIFICATION OF MARKER GENES ---
#==============================================================================#

# Change the active Assay to 'RNA' to find markers at the gene expression level
DefaultAssay(integrated_seurat) <- "RNA"

# Find markers for each cluster
# min.pct: detected in at least 25% of cells in the cluster
# 
logfc.threshold: minimum log-fold change of 0.25
all_markers <- FindAllMarkers(integrated_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the list of markers to a CSV file
write.csv(all_markers, "all_cluster_markers.csv", row.names = FALSE)

# Select the top 10 markers for each cluster based on Log2 Fold Change
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# --- 6.2.
Visualization with Heatmap ---
# NOTE: ScaleData in the integration step was run only on 'VariableFeatures'.
# For an accurate heatmap with genes that might not be variable, it is good to re-scale the data
# using only the genes of interest.
integrated_seurat <- ScaleData(integrated_seurat, features = top10_markers$gene)

# Generate the heatmap
DoHeatmap(
  subset(integrated_seurat, downsample = 300), # Downsample for performance
  features = top10_markers$gene,
  raster = TRUE # Uses rasterization to generate a smaller and faster file
)


#==============================================================================#
# --- 7. APPENDIX: LEGACY INTEGRATION METHOD (SEURAT < v5) ---
#==============================================================================#

# This method (FindIntegrationAnchors/IntegrateData) still works, but
# `IntegrateLayers` is currently recommended.
# # Load the list of QC'd objects
# rds_files <- list.files("RDS", pattern = "\\.rds$", full.names = TRUE)
# seurat_list_legacy <- lapply(rds_files, readRDS)
# 
# # Normalize and find variable features for each object in the list
# for (i in 1:length(seurat_list_legacy)) {
#   seurat_list_legacy[[i]] <- NormalizeData(seurat_list_legacy[[i]], verbose = FALSE)
#   seurat_list_legacy[[i]] <- FindVariableFeatures(seurat_list_legacy[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# }
# 
# # Find integration anchors
# integration_anchors <- FindIntegrationAnchors(object.list = seurat_list_legacy, dims = 1:30)
# 
# # Integrate the data
# integrated_legacy <- IntegrateData(anchorset = integration_anchors, dims = 1:30)
# 
# # Set the default assay for the integrated object
# DefaultAssay(integrated_legacy) <- "integrated"
# 
# # Continue with ScaleData, RunPCA, etc.
# integrated_legacy <- ScaleData(integrated_legacy, verbose = FALSE)
# integrated_legacy <- RunPCA(integrated_legacy, npcs = 30, verbose = FALSE)
# integrated_legacy <- RunUMAP(integrated_legacy, reduction = "pca", dims = 1:20)
# integrated_legacy <- FindNeighbors(integrated_legacy, reduction = "pca", dims = 1:20)
# integrated_legacy <- FindClusters(integrated_legacy, resolution = 0.5)
# 
# # Visualize
# DimPlot(integrated_legacy, reduction = "umap")
