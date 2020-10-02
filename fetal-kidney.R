# Fetal Kidney

# Load necessary libraries.
library(tidyverse)
library(Seurat)
library(Matrix)
library(patchwork)
library(sctransform)
library(openxlsx)
library(DESeq2)
library(gt)

# Read in 10X datasets and create Seurat object for each before merging them all
# together.
fk_w17z1_data <- Read10X(data.dir = 'data/GSM3534656_H17w_Z1_raw_counts/')
fk_w17z2_data <- Read10X(data.dir = 'data/GSM3534657_H17w_Z2_raw_counts/')
fk_w15z1_data <- Read10X(data.dir = 'data/GSM3534658_H15w_Z1_raw_counts/')
fk_w15z2_data <- Read10X(data.dir = 'data/GSM3534659_H15w_Z2_raw_counts/')
# Create Seurat Objects
fk_w17z1 <- CreateSeuratObject(counts = fk_w17z1_data, project = 'fkw17z1')
fk_w17z2 <- CreateSeuratObject(counts = fk_w17z2_data, project = 'fkw17z2')
fk_w15z1 <- CreateSeuratObject(counts = fk_w15z1_data, project = 'fkw15z1')
fk_w15z2 <- CreateSeuratObject(counts = fk_w15z2_data, project = 'fkw15z2')
# Merge together
fk_merged <- merge(fk_w15z1, y = c(fk_w15z2, fk_w17z1, fk_w17z2), 
                   add.cell.ids = c('w15z1', 'w15z2', 'w17z1', 'w17z2'), 
                   project = 'FK_Merged')
fk_merged
fk_merged@meta.data

# QC
# Find percentage of reads that map to mitochondrial genome:
fk_merged[["percent.mt"]] <- PercentageFeatureSet(fk_merged, pattern = "^MT-", 
                                                  assay = 'RNA')
head(fk_merged@meta.data, 10)
# Visualize QC metrics as a violin plot
VlnPlot(fk_merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
        ncol = 3)

# QC Filtering
# The two time points seem to have a difference. Week 17 shows much higher
# counts, features, and percent mt expression than week 15. Will use 15% cutoff
# for percent.mt. Filter out cells with < 200 and > 4500 features, > 20k counts.
fk_filtered <- subset(fk_merged, subset = nFeature_RNA < 4500 & 
                        nFeature_RNA > 200 & nCount_RNA < 20000 & 
                        percent.mt < 15)

# Normalize
# Next step is normalizing. In case there are any cell cycle effects, it may be
# helpful to find and then regress out cell cycle variability. To do so,
# normalize first, and then identify cell cycle phases before regressing it out
# along with percent.mt and library size.
# Normalize using NormalizeData with default values.
fk_normalized <- NormalizeData(fk_filtered)

# Cell Cycle Scoring
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
# We can segregate this list into markers of G2/M phase and markers of S phase.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fk_cellcyclephased <- CellCycleScoring(fk_normalized, s.features = s.genes, 
                             g2m.features = g2m.genes, set.ident = TRUE)
# View cell cycle scores and phase assignments
head(fk_cellcyclephased@meta.data, 5)

# SCTransform
# Now we can run SCTransform() and regress out variability due to phase,
# percent.mt, batch, and library size.
fk_transform <- SCTransform(fk_cellcyclephased, assay = 'RNA', 
                            new.assay.name = 'SCT', 
                            vars.to.regress = c('orig.ident', 'percent.mt', 
                                                'nFeature_RNA', 'nCount_RNA', 
                                                'S.Score', 'G2M.Score'))
saveRDS(fk_transform, file = 'results/fk_transform_tran_2019.rds')
# PCA/linear dimensionality reduction
# Now we run PCA for linear dimensionality reduction.
fk_pca <- RunPCA(fk_transform)
# Examine and visualize PCA results a few different ways
print(fk_pca[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(fk_pca, dims = 1:2, reduction = "pca")
DimPlot(fk_pca, reduction = "pca")
DimHeatmap(fk_pca, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the 'dimensionality' of the dataset.
ElbowPlot(ak_pca) # At least 20 PCs.

# Cluster the Cells.
fk_neighbours <- FindNeighbors(fk_pca, dims = 1:30)
fk_clusters <- FindClusters(fk_neighbours, resolution = 2)
# Look at cluster IDs of the first 5 cells
head(Idents(fk_clusters), 5)

# Run non-linear dimensional reduction.
fk_umap <- RunUMAP(fk_clusters, dims = 1:30, umap.method = 'umap-learn', 
                   metric = 'correlation')
DimPlot(fk_umap, reduction = "umap", label = T)

# Save the object for sharing and to avoid computationally intensive steps
# already conducted.
saveRDS(fk_umap, file = "results/fk_umap_tran_2019.rds")



