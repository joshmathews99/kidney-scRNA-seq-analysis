
# Load necessary libraries.
library(tidyverse)
library(Seurat)
library(Matrix)
library(patchwork)
library(sctransform)
library(patchwork)

# Read in raw data. This is a sample from a 41 y/o male.
raw_counts <- read.table(file = 'data/GSM4008620_Adult-Kidney3_dge.txt', 
                         sep = '\t', header = T, row.names = 1)
glimpse(raw_counts)
View(raw_counts)

# Create Seurat Object. Min cells for gene to be expressed in is 3. Min genes
# per cell is 100.
ak3 <- CreateSeuratObject(counts = raw_counts, min.cells = 3, 
                          min.features = 100, project = "ak3")
ak3
# Explore metadata
head(ak3@meta.data, 5)

# Find percentage of reads that map to mitochondrial genome:
ak3[["percent.mt"]] <- PercentageFeatureSet(ak3, pattern = "^MT-", 
                                            assay = 'RNA')
head(ak3@meta.data, 10)

# Visualize QC metrics as a violin plot
VlnPlot(ak3, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
# It looks like the percent.mt values are really high. Indicates stress. Might
# be helpful to look at the other available adult kidney sample(s). ~ 40% mt.

# Will now consider sample adult kidney 2 (from 66 y/o male).
raw_counts_ak2 <- read.table(file = 'data/GSM4008619_Adult-Kidney2_dge.txt',
                             sep = '\t', header = T, row.names = 1)
glimpse(raw_counts_ak2)
View(raw_counts_ak2)
# Create Seurat Object. Min cells in which gene is expressed set to 3. Min genes
# expressed per cell set to 100.
ak2 <- CreateSeuratObject(counts = raw_counts_ak2, min.cells = 3, 
                          min.features = 100, project = 'ak2')
ak2
# Explore metadata
head(ak2@meta.data, 10)
# Find % reads mapping to mitochondrial genome:
ak2[['percent.mt']] <- PercentageFeatureSet(ak2, pattern = '^MT-', 
                                            assay = 'RNA')
head(ak2@meta.data, 10)
# Visualize QC metrics as a violin plot
VlnPlot(ak2, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3)
# High mitochondrial gene expression it looks like. ~ 45%

# Will now look at the other two adult kidney samples (both from 57 y/o male).
# Sample 4-1
raw_counts_ak4_1 <- read.table(file = 'data/GSM4008621_Adult-Kidney4-1_dge.txt', 
                               sep = '\t', header = T, row.names = 1)
ak4_1 <- CreateSeuratObject(counts = raw_counts_ak4_1, min.cells = 3, 
                          min.features = 100, project = "ak4_1")
ak4_1
ak4_1[["percent.mt"]] <- PercentageFeatureSet(ak4_1, pattern = "^MT-", 
                                            assay = 'RNA')
head(ak4_1@meta.data, 10)
# Visualize QC metrics as a violin plot
VlnPlot(ak4_1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
        ncol = 3)
# Mitochondrial gene expression is even higher ~ 48%.
# Sample 4-2
raw_counts_ak4_2 <- read.table(file = 'data/GSM4008622_Adult-Kidney4-2_dge.txt', 
                               sep = '\t', header = T, row.names = 1)
ak4_2 <- CreateSeuratObject(counts = raw_counts_ak4_2, min.cells = 3, 
                            min.features = 100, project = "ak4_2")
ak4_2
ak4_2[["percent.mt"]] <- PercentageFeatureSet(ak4_2, pattern = "^MT-", 
                                              assay = 'RNA')
head(ak4_2@meta.data, 10)
# Visualize QC metrics as a violin plot
VlnPlot(ak4_2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
        ncol = 3)
# Mitochondrial gene expression is ~ 45%.

###

# Will stick to first ak3 dataset due to lowest level of mt gene expression and subset into a few datasets for clustering analysis: 20, 30, 40, 50, 100 % mt gene expression. This is a quality control step.
# Look at feature relationships:
VlnPlot(ak3, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(ak3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ak3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 # counts and features strongly positively correlated.
# Filter out: features > 2000, features < 200, percent.mt > 20,30,40,50,100%.
# This yields five subsets of the data to work with. Will see if there are any
# major differences between in downstream analyses.
ak3_20mt <- subset(ak3, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & 
                     percent.mt < 20)
ak3_30mt <- subset(ak3, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & 
                     percent.mt < 30)
ak3_40mt <- subset(ak3, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & 
                     percent.mt < 40)
ak3_50mt <- subset(ak3, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & 
                     percent.mt < 50)
ak3_100mt <- subset(ak3, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & 
                      percent.mt < 100)

# Will proceed with SCTransform on each dataset including the original dataset.
# This step replaces NormalizeData, FindVariableFeatures, and ScaleData. Finds
# top 2000(?) variable features and allows for use of more PCs.
ak3 <- SCTransform(ak3, vars.to.regress = 'percent.mt')
ak3_20mt <- SCTransform(ak3_20mt, vars.to.regress = 'percent.mt')
ak3_30mt <- SCTransform(ak3_30mt, vars.to.regress = 'percent.mt')
ak3_40mt <- SCTransform(ak3_40mt, vars.to.regress = 'percent.mt')
ak3_50mt <- SCTransform(ak3_50mt, vars.to.regress = 'percent.mt')
ak3_100mt <- SCTransform(ak3_100mt, vars.to.regress = 'percent.mt')
# Note on SCTransform: The object[["SCT"]]@scale.data contains the residuals
# (normalized values), and is used directly as input to PCA. To assist with
# visualization and interpretation. we also convert Pearson residuals back to
# ‘corrected’ UMI counts. You can interpret these as the UMI counts we would
# expect to observe if all cells were sequenced to the same depth. The
# ‘corrected’ UMI counts are stored in object[["SCT"]]@counts. We store
# log-normalized versions of these corrected counts in object[["SCT"]]@data,
# which are very helpful for visualization. You can use the corrected
# log-normalized counts for differential expression and integration. However, in
# principle, it would be most optimal to perform these calculations directly on
# the residuals (stored in the scale.data slot) themselves. This is not
# currently supported in Seurat v3, but will be soon. (SCTransform vignette)

# Perform linear dimensional reduction
# Run PCA using prior 2000 top variable features as default.
ak3 <- RunPCA(ak3)
ak3_20mt <- RunPCA(ak3_20mt)
ak3_30mt <- RunPCA(ak3_30mt)
ak3_40mt <- RunPCA(ak3_40mt)
ak3_50mt <- RunPCA(ak3_50mt)
ak3_100mt <- RunPCA(ak3_100mt)
# Visualize PCA results.
print(ak3[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(ak3, dims = 1:2, reduction = "pca")
DimPlot(ak3, reduction = "pca")
DimHeatmap(ak3, dims = 1:15, cells = 500, balanced = TRUE)

# Determine 'dimensionality' of dataset
ElbowPlot(ak3) # 19
ElbowPlot(ak3_20mt) # 20
ElbowPlot(ak3_30mt) # 19
ElbowPlot(ak3_40mt) # 19
ElbowPlot(ak3_50mt) # 19
ElbowPlot(ak3_100mt) # 19
# Seems like using many PCs (~ 20 at least) may be necessary. Interestingly, the
# Seurat vignette for sctransform uses 30 PCs and encourages the use of more
# PCs. May be helpful to try a few (20, 30) for now.

# Cluster the cells and check if clustering worked. Resolution depends on number
# of cells in each sample. ak3 (10000), ak3_20mt (6235), ak3_30mt (8987),
# ak3_40mt (9676), ak3_50mt (9875), ak3_100mt (9955). Try 1.0 and 2.0 initially.
# This yields 24 datasets of clusters.
ak3_20pc_1r <- FindNeighbors(ak3, dims = 1:20) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_20pc_2r <- FindNeighbors(ak3, dims = 1:20) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_30pc_1r <- FindNeighbors(ak3, dims = 1:30) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')
ak3_30pc_2r <- FindNeighbors(ak3, dims = 1:30) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

ak3_20mt_20pc_1r <- FindNeighbors(ak3_20mt, dims = 1:20) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_20mt_20pc_2r <- FindNeighbors(ak3_20mt, dims = 1:20) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_20mt_30pc_1r <- FindNeighbors(ak3_20mt, dims = 1:30) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')
ak3_20mt_30pc_2r <- FindNeighbors(ak3_20mt, dims = 1:30) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

ak3_30mt_20pc_1r <- FindNeighbors(ak3_30mt, dims = 1:20) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_30mt_20pc_2r <- FindNeighbors(ak3_30mt, dims = 1:20) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_30mt_30pc_1r <- FindNeighbors(ak3_30mt, dims = 1:30) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')
ak3_30mt_30pc_2r <- FindNeighbors(ak3_30mt, dims = 1:30) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

ak3_40mt_20pc_1r <- FindNeighbors(ak3_40mt, dims = 1:20) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_40mt_20pc_2r <- FindNeighbors(ak3_40mt, dims = 1:20) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_40mt_30pc_1r <- FindNeighbors(ak3_40mt, dims = 1:30) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')
ak3_40mt_30pc_2r <- FindNeighbors(ak3_40mt, dims = 1:30) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

ak3_50mt_20pc_1r <- FindNeighbors(ak3_50mt, dims = 1:20) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_50mt_20pc_2r <- FindNeighbors(ak3_50mt, dims = 1:20) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_50mt_30pc_1r <- FindNeighbors(ak3_50mt, dims = 1:30) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')
ak3_50mt_30pc_2r <- FindNeighbors(ak3_50mt, dims = 1:30) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

ak3_100mt_20pc_1r <- FindNeighbors(ak3_100mt, dims = 1:20) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_100mt_20pc_2r <- FindNeighbors(ak3_100mt, dims = 1:20) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:20, umap.method = 'umap-learn', metric = 'correlation')
ak3_100mt_30pc_1r <- FindNeighbors(ak3_100mt, dims = 1:30) %>% FindClusters(resolution = 1) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')
ak3_100mt_30pc_2r <- FindNeighbors(ak3_100mt, dims = 1:30) %>% FindClusters(resolution = 2) %>% RunUMAP(dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')

# Run non-linear dimensional reduction (tSNE/UMAP) and now visualize them. Note
# that you can set `label = TRUE` or use the LabelClusters function to help
# label individual clusters.
plots_ak3_20mt <- (DimPlot(ak3_20mt_20pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_20mt_20pc_1r') + (DimPlot(ak3_20mt_20pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_20mt_20pc_2r'))) / (DimPlot(ak3_20mt_30pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_20mt_30pc_1r') + (DimPlot(ak3_20mt_30pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_20mt_30pc_2r')))
plots_ak3_30mt <- (DimPlot(ak3_30mt_20pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_30mt_20pc_1r') + (DimPlot(ak3_30mt_20pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_30mt_20pc_2r'))) / (DimPlot(ak3_30mt_30pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_30mt_30pc_1r') + (DimPlot(ak3_30mt_30pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_30mt_30pc_2r')))
plots_ak3_40mt <- (DimPlot(ak3_40mt_20pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_40mt_20pc_1r') + (DimPlot(ak3_40mt_20pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_40mt_20pc_2r'))) / (DimPlot(ak3_40mt_30pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_40mt_30pc_1r') + (DimPlot(ak3_40mt_30pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_40mt_30pc_2r')))
plots_ak3_50mt <- (DimPlot(ak3_50mt_20pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_50mt_20pc_1r') + (DimPlot(ak3_50mt_20pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_50mt_20pc_2r'))) / (DimPlot(ak3_50mt_30pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_50mt_30pc_1r') + (DimPlot(ak3_50mt_30pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_50mt_30pc_2r')))
plots_ak3_100mt <- (DimPlot(ak3_100mt_20pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_100mt_20pc_1r') + (DimPlot(ak3_100mt_20pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_100mt_20pc_2r'))) / (DimPlot(ak3_100mt_30pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_100mt_30pc_1r') + (DimPlot(ak3_100mt_30pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_100mt_30pc_2r')))
plots_ak3 <- (DimPlot(ak3_20pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_20pc_1r') + (DimPlot(ak3_20pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_20pc_2r'))) / (DimPlot(ak3_30pc_1r, reduction = "umap", label = T) + labs(title = 'ak3_30pc_1r') + (DimPlot(ak3_30pc_2r, reduction = "umap", label = T) + labs(title = 'ak3_30pc_2r')))
plots_ak3_20mt
plots_ak3_30mt
plots_ak3_40mt
plots_ak3_50mt
plots_ak3_100mt
plots_ak3

# save objects for later
saveRDS(pbmc, file = "results/pbmc_tutorial.rds")

# Finding differentially expressed features (cluster biomarkers)

FeaturePlot(ak3_40mt_30pc_1r, features = c("SLC12A3"))
FeaturePlot(ak3_40mt_30pc_1r, features = c("SLC12A1"))
FeaturePlot(ak3_40mt_30pc_1r, features = c("AQP2"))
FeaturePlot(ak3_40mt_30pc_1r, features = c("HNF4A")) # PT
FeaturePlot(ak3_40mt_30pc_1r, features = c("MAFB")) # Podocytes
FeaturePlot(ak3_40mt_30pc_1r, features = c("NPHS2"))
FeaturePlot(ak3_30pc_1r, features = c("NPHS1"))



# Assigning cell type identity to clusters



