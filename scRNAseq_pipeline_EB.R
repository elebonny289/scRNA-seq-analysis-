# =============================================================================
# scRNA-seq Analysis Pipeline — Immune Cells (Generic Template)
#
# Description:
#   A modular, reusable pipeline for single-cell RNA-seq analysis of
#   immune cell populations from mouse or human tissue.
#   Input: BD Rhapsody MEX format (adaptable to 10x Genomics CellRanger output).
#
#   Workflow:
#     1.  Data loading from MEX format (scalable to any number of samples)
#     2.  Quality control and filtering
#     3.  Normalisation, PCA, Harmony batch correction
#     4.  Clustering and UMAP visualisation
#     5.  Cell type annotation via module scores
#     6.  Pseudobulk differential expression (DESeq2)
#     7.  Gene set enrichment (ssGSEA via escape)
#     8.  Visualisation (UMAP, VlnPlot, FeaturePlot, EnhancedVolcano)
#
# How to adapt this template:
#   - Edit the CONFIGURATION section (Section 0) to match your experiment
#   - Provide your own sample metadata table (Section 1)
#   - Replace gene signatures in Section 6 with those relevant to your cell types
#     if you don't have a signature, you can skip this step and assign cluster identities
#     based on canonical marker genes using FindAllMarkers() and FeaturePlot()
#   - Adjust QC thresholds in Section 3 based on your data quality
#
# Tested with:
#   Seurat 4.x · Harmony · R 4.2+
# =============================================================================
# 0. CONFIGURATION — edit this section for your experiment
# =============================================================================

# Working directory: folder containing all sample subdirectories
WORKING_DIR <- "/path/to/your/data/"

# Number of PCA dimensions to use downstream
N_DIMS <- 30

# Seurat clustering resolution (higher = more clusters)
CLUSTER_RESOLUTION <- 0.4

# QC thresholds
MIN_FEATURES  <- 200    # minimum genes per cell
MAX_FEATURES  <- 5000   # maximum genes per cell (removes doublets)
MAX_MT_PCT    <- 15     # maximum mitochondrial read percentage

# Mitochondrial gene prefix ("^mt-" for mouse, "^MT-" for human)
MT_PATTERN <- "^mt-"

# Organism database for GO enrichment ("org.Mm.eg.db" mouse / "org.Hs.eg.db" human)
ORG_DB <- "org.Mm.eg.db"
#the names of the genes is idifferent from HUMAN/Mouse
#example HUMAN = "CCR2";  Mouse = "Ccr2" 
# Number of CPU cores for ssGSEA (adjust to your machine)
N_CORES <- 4
# Set working directory
setwd(WORKING_DIR)
# =============================================================================
# 1. LIBRARIES
# =============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(R.utils)
library(SeuratWrappers)
library(SeuratDisk)
library(harmony)
library(escape)
library(clusterProfiler)
library(org.Mm.eg.db)   # swap for org.Hs.eg.db if using human data
library(EnhancedVolcano)

# =============================================================================
# 2. SAMPLE METADATA TABLE
# =============================================================================
# Edit this table to describe your experiment.
# Each row = one sample.
# You can add/remove columns (e.g. "Tissue", "Timepoint") as needed.
#
# Column guide:
#   tag        : numeric index used to build the file path (SampleTag01, 02, ...)
#   sample_id  : short label for the Seurat project slot
#   mouse_id   : biological replicate identifier
#   group_1    : primary grouping variable (e.g. Genotype, Treatment, Disease)
#   group_2    : secondary grouping variable (e.g. Condition, Timepoint)
#
# Example below: 12 samples, 2 genotypes x 2 conditions x 3 replicates each.
# Replace values with your own.

sample_metadata <- data.frame(tag = 1:12, sample_id = paste0("sample", 1:12),
  mouse_id  = paste0("Mouse", 1:12),
  group_1   = c(rep("GroupA_WT", 3), rep("GroupA_KO", 3),
                rep("GroupB_WT", 3), rep("GroupB_KO", 3)),
  group_2   = c(rep("Condition1", 6), rep("Condition2", 6)),
  stringsAsFactors = FALSE)
# Preview the metadata table
print(sample_metadata)

# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================

#' Build file paths for BD Rhapsody MEX output
#' Assumes a folder structure like:
#'   PROJECT_SampleTag01_mm/
#'     PROJECT_SampleTag01_mm_RSEC_MolsPerCell_MEX/
#'       matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz
#'
#' For 10x CellRanger output, replace with:
#'   file.path("sample_folder", c("matrix.mtx.gz","barcodes.tsv.gz","features.tsv.gz"))
#'
#' @param tag        Integer sample tag number
#' @param file       File name (matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz)
#' @param project    Project prefix string (e.g. "MyProject")
#' @param species    Species suffix string (e.g. "mm" for mouse, "hs" for human)
rhapsody_path <- function(tag, file, project = "MyProject", species = "mm") {
  tag_str <- sprintf("SampleTag%02d_%s", tag, species)
  dir     <- sprintf("%s_%s", project, tag_str)
  base    <- sprintf("%s_%s_RSEC_MolsPerCell_MEX", project, tag_str)
  file.path(dir, base, file)
}


#' Create a Seurat object from BD Rhapsody UMI count matrix
#'
#' Adds mitochondrial percentage and transcriptomic complexity (log10GenesPerUMI)
#' to cell metadata automatically.
#'
#' @param umi      Count matrix (genes x cells) from ReadMtx()
#' @param sample   Project name / label for the Seurat object
#' @param metadata Optional per-cell metadata data.frame
#' @return         Seurat object with QC metrics in metadata
CreateSeuratFromRhapsody <- function(umi, sample = NULL, metadata = NULL) {
  obj <- CreateSeuratObject(counts = umi, min.cells = 0, min.features = 0, project  = sample,
    names.delim  = "-", meta.data    = metadata)
  obj[["percent.mt"]]        <- PercentageFeatureSet(obj, pattern = MT_PATTERN)
  obj[["log10GenesPerUMI"]]  <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  return(obj)}
# =============================================================================
# 4. LOAD DATA — loop over all samples
# =============================================================================

seurat_list <- lapply(seq_len(nrow(sample_metadata)), function(i) {

  row <- sample_metadata[i, ]
  tag <- row$tag

  message(sprintf("[%d/%d] Loading %s (%s / %s) ...",
                  tag, nrow(sample_metadata),
                  row$sample_id, row$group_1, row$group_2))

  # --- Load MEX matrix -------------------------------------------------------
  # To use 10x CellRanger output instead, replace the three lines below with:
  #   mtx <- Read10X(data.dir = file.path("path/to/sample", row$sample_id))
mtx <- ReadMtx(mtx      = rhapsody_path(tag, "matrix.mtx.gz"),
               cells    = rhapsody_path(tag, "barcodes.tsv.gz"),
              features = rhapsody_path(tag, "features.tsv.gz"))

  # --- Create Seurat object --------------------------------------------------
  obj <- CreateSeuratFromRhapsody(umi = mtx, sample = row$sample_id)

  # --- Add sample-level metadata --------------------------------------------
  obj[["Mouse"]]   <- row$mouse_id
  obj[["Group1"]]  <- row$group_1   # rename to your variable (e.g. "Genotype")
  obj[["Group2"]]  <- row$group_2   # rename to your variable (e.g. "Condition")

  return(obj)
})

names(seurat_list) <- sample_metadata$sample_id

# =============================================================================
# 5. QUALITY CONTROL — per-sample plots
# =============================================================================

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI")

for (i in seq_along(seurat_list)) {
  obj  <- seurat_list[[i]]
  name <- names(seurat_list)[i]

  p_vln <- VlnPlot(obj, features = qc_features, pt.size = 0.05)
  p_mt  <- FeatureScatter(obj, "nCount_RNA", "percent.mt",   pt.size = 0.1)
  p_gen <- FeatureScatter(obj, "nCount_RNA", "nFeature_RNA", pt.size = 0.1)

  print(p_vln + p_mt + p_gen)}


# =============================================================================
# 6. MERGE AND FILTER
# =============================================================================

merged <- merge(seurat_list[[1]], y = seurat_list[-1], project = "scRNAseq_project")

rm(seurat_list); gc()

# Apply QC thresholds (defined in Section 0)
merged <- subset( merged, subset = nFeature_RNA >= MIN_FEATURES & 
           nFeature_RNA <= MAX_FEATURES & percent.mt   <= MAX_MT_PCT)

saveRDS(merged, "merged_raw.rds")
# merged <- readRDS("merged_raw.rds") #

# =============================================================================
# 7. NORMALISATION, PCA, HARMONY BATCH CORRECTION
# =============================================================================

integrated <- merged %>% NormalizeData() %>%
  FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)

# Elbow plot — use to check how many PCs to keep (set N_DIMS in Section 0)
ElbowPlot(integrated, ndims = 50, reduction = "pca")

# Harmony removes technical batch effects between samples
integrated <- RunHarmony(integrated, group.by.vars = "orig.ident")

# UMAP + clustering on Harmony embedding
integrated <- integrated %>%
  RunUMAP(reduction = "harmony", dims = 1:N_DIMS) %>%
  FindNeighbors(reduction = "pca", dims = 1:N_DIMS) %>%
  FindClusters(resolution = CLUSTER_RESOLUTION)

# Diagnostic UMAP plots
DimPlot(integrated, split.by = "orig.ident", ncol = 4)  + ggtitle("Per sample") 
DimPlot(integrated, group.by = "Group2")                 + ggtitle("By condition")
#split.by indicate to do n UMAP with only the cells from the condition 
#group.by indicate to cluster all the cells in the UMAP 
saveRDS(integrated, "integrated.rds")
# integrated <- readRDS("integrated.rds")

# =============================================================================
# 8. CLUSTER MARKER GENES
# =============================================================================

all_markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top 5 markers per cluster — heatmap overview
top5 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)

DoHeatmap(integrated, features = top5$gene) + NoLegend()

write.csv(all_markers, file = "all_cluster_markers.csv")

# =============================================================================
# 9. CELL TYPE ANNOTATION via Add Module Scores
# =============================================================================
# Replace these gene lists with signatures relevant to your cell population.
# Sources: published papers, CellMarker database, PanglaoDB, etc.
#
# Each list should contain 5–20 marker genes for a known cell type.
# After running AddModuleScore(), visualise with VlnPlot per cluster
# to decide which cluster corresponds to which cell type.

cell_type_signatures <- list(
  # --- Example: mouse lung immune cells (replace with your own) -------------
  Alveolar_Macrophage    = c("Ear2", "Siglecf", "Pparg", "Mertk", "Marco"),
  Interstitial_Macrophage = c("Cx3cr1", "C1qa", "C1qb", "C1qc"),
  Classical_Monocyte     = c("Ly6c2", "Ccr2", "S100a8", "Csf3r"),
  NonClassical_Monocyte  = c("Ace", "Cd14", "Itgam", "Ear2"),
  Neutrophil             = c("Csf3r", "S100a8", "Ly6g", "Mmp8"),
  Eosinophil             = c("Il5ra", "Siglecf", "Prg2", "Epx"),
  cDC1                   = c("Xcr1", "Clec9a", "Cadm1", "Irf8", "Batf3"),
  cDC2                   = c("Cd209a", "Itgam", "Sirpa", "Irf4", "Tbx21")
)

# Add all module scores to the Seurat object
for (nm in names(cell_type_signatures)) {
  integrated <- AddModuleScore(object   = integrated,features = list(cell_type_signatures[[nm]]),
    name     = nm)}

# Visualise scores per cluster to guide annotation
score_cols <- paste0(names(cell_type_signatures), "1")
for (sc in score_cols) {
  print(VlnPlot(integrated, features = sc, pt.size  = 0, 
            group.by = "seurat_clusters", assay    = "RNA") + ggtitle(sc))}


# =============================================================================
# 10. CELL TYPE ASSIGNMENT
# =============================================================================
# After inspecting module scores and marker genes, manually assign
# a cell type label to each cluster number.
#
# Edit the vector below: one label per cluster, in order (cluster 0, 1, 2, ...)
# Use "unassigned" for clusters you cannot confidently annotate.

Idents(integrated) <- "seurat_clusters"

# --- EDIT THIS VECTOR --------------------------------------------------------
cell_type_labels <- c(
  "CellType_A",    # cluster 0
  "CellType_B",    # cluster 1
  "CellType_C",    # cluster 2
  "CellType_A",    # cluster 3  (same type, different sub-cluster)
  "CellType_D",    # cluster 4
  "unassigned"     # cluster 5
  # add one entry per cluster
)
# -----------------------------------------------------------------------------

names(cell_type_labels) <- levels(integrated)
integrated <- RenameIdents(integrated, cell_type_labels)
integrated[["Cell_Type"]] <- Idents(integrated)

# Remove unassigned cells
integrated <- subset(integrated, idents = "unassigned", invert = TRUE)
# UMAP with cell type labels
DimPlot(integrated, label = TRUE, label.box = TRUE, label.size = 3,
        group.by = "Cell_Type") +
  ggtitle("UMAP — Annotated cell types")

# Cell counts per group
table(integrated$Cell_Type, integrated$Group1)
table(integrated$Cell_Type, integrated$Group2)

saveRDS(integrated, "integrated_annotated.rds")
# integrated <- readRDS("integrated_annotated.rds")

# =============================================================================
# 11. VOLCANO PLOTS (EnhancedVolcano)
# =============================================================================
# Reads a DESeq2 output CSV and generates a volcano plot.
# Edit `deg_file` and `genes_to_label` for your analysis.

deg_file      <- "DEGs_GroupA_vs_GroupB.csv"  # path to your DESeq2 output file
genes_to_label <- c("GeneA", "GeneB", "GeneC")   # replace with your genes of interest

if (file.exists(deg_file)) {
  deg <- read.csv(deg_file, row.names = 1)

  EnhancedVolcano(deg, lab = rownames(deg), x = "avg_log2FC", y = "p_val", labSize = 4,
                  FCcutoff = 0.25, pCutoff = 0.05, title = "Group A vs Group B — CellType_A",
                  drawConnectors  = TRUE, widthConnectors = 0.75, selectLab = genes_to_label)}


# =============================================================================
# 12. GENE SET ENRICHMENT — ssGSEA (escape)
# =============================================================================
# Downloads MSigDB gene sets. Computationally intensive — recommended to run on a server or HPC.

GS <- getGeneSets(library = c("H", "C5", "C7", "C2"),
  species = "Mus musculus") #change argument species if you have other species 

integrated <- runEscape(integrated, method = "ssGSEA", gene.sets = GS, groups = 1000,
  min.size  = 5, cores = N_CORES)

integrated <- performNormalization(input.data = integrated, gene.sets = GS, assay = "escape",
  scale.factor = integrated$nFeature_RNA)

integrated <- performPCA(integrated, n.dim = 1:10)
# Note: runEscape() creates a new assay called "escape" inside the Seurat object,
# storing the enrichment scores alongside the RNA assay.

# =============================================================================
# 13. GO / PATHWAY ENRICHMENT (clusterProfiler)
# =============================================================================

# --- Per-cluster GO on top marker genes --------------------------------------
top100 <- all_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 100)

# Convert gene symbols to Entrez IDs
entrez_per_cluster <- lapply(split(top100$gene, top100$cluster),
  function(genes) {bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = ORG_DB)$ENTREZID})

# Compare GO enrichment across clusters
go_cluster <- compareCluster(fun = "enrichGO", geneClusters = entrez_per_cluster, OrgDb = ORG_DB)
go_cluster <- setReadable(go_cluster, OrgDb = ORG_DB, keyType = "ENTREZID")
dotplot(go_cluster, font.size = 8) + ggtitle("GO enrichment — per cluster")

# Reactome pathway enrichment
reactome_cluster <- compareCluster(fun  = "enrichPathway", geneClusters = entrez_per_cluster,
  OrgDb = ORG_DB)
reactome_cluster <- setReadable(reactome_cluster, OrgDb = ORG_DB, keyType = "ENTREZID")
dotplot(reactome_cluster, font.size = 8) + ggtitle("Reactome pathways — per cluster")


# --- GO on pseudobulk DEGs ---------------------------------------------------
# Helper: run enrichGO on a list of gene symbols and plot results
run_enrichGO <- function(genes, label, ontology = "ALL") {
  eg  <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = ORG_DB)
  ego <- enrichGO(gene  = eg$ENTREZID, OrgDb = ORG_DB, 
    ont = ontology, pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
  print(dotplot(ego, showCategory = 20, font.size = 8, title = label)) return(ego)}

# Example: load one DEG table and run GO on up-regulated genes 
deg_example <- read.csv("DEG_CellType_A_Condition1_WT_vs_KO.csv", row.names = 1)

genes_up_A <- rownames(subset(deg_example, avg_log2FC >  0.25 & p_val_adj < 0.05))
genes_up_B <- rownames(subset(deg_example, avg_log2FC < -0.25 & p_val_adj < 0.05))

ego_A <- run_enrichGO(genes_up_A, "GO — genes up in Group A")
ego_B <- run_enrichGO(genes_up_B, "GO — genes up in Group B")
# =============================================================================
# 14. FEATURE EXPRESSION — genes of interest
# =============================================================================
# Replace with your genes of interest

genes_of_interest <- c("GeneX", "GeneY", "GeneZ")

# UMAP feature plots
FeaturePlot(integrated, features = genes_of_interest,
            order = TRUE, pt.size = 0.5) & theme_classic()

# Violin plots split by experimental group
VlnPlot(integrated, features   = genes_of_interest,
        split.by   = "Group1", group.by   = "Cell_Type", split.plot = FALSE,
        pt.size    = 0)
