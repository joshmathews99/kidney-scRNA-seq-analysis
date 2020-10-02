
# Load necessary libraries.
library(tidyverse)
library(Seurat)
library(Matrix)
library(patchwork)
library(sctransform)
library(openxlsx)
library(DESeq2)
library(gt)

# Read in processed data. Data are from iPSC-derived kidney organoids using
# Takasato protocol, harvested at d26. 3 batches will be grouped together.
counts_1_2 <- read.table(file = 'data/GSE118184_Takasato_iPS_Batch1_2.dge.txt', 
                         sep = '\t', header = T, row.names = 1)
glimpse(counts_1_2) # 24 383 obs of 15 459 cells
counts_3 <- read.table(file = 'data/GSE118184_Takasato_iPS_Batch3.dge.txt',
                       sep = '\t', header = T, row.names = 1)
glimpse(counts_3) # 15 320 obs of 5 269 cells

# Create Seurat Object. Min cells for gene to be expressed in is 10. Min genes
# per cell is 100.
tak1_2 <- CreateSeuratObject(counts = counts_1_2, min.cells = 10, 
                             min.features = 100, project = "tak1_2")
tak1_2 # 24 383 features across 15 459 samples
tak3 <- CreateSeuratObject(counts = counts_3, min.cells = 10, 
                             min.features = 100, project = "tak3")
tak3 # 15 135 features across 5 269 samples
tak3$orig.ident <- "Batch3"
# Explore metadata. Includes orig.ident (Batch 1, 2), RNA counts, and
# Feature counts.
tak1_2@meta.data
tak3@meta.data
# Merge together
tak_merged <- merge(tak1_2, y = c(tak3), 
                   add.cell.ids = c('tak1_2', 'tak3'), 
                   project = 'Tak_Merged')
tak_merged # 24 712 features across 20 728 samples
tak_merged@meta.data

# QC
# Find percentage of reads that map to mitochondrial genome:
tak_merged[["percent.mt"]] <- PercentageFeatureSet(tak_merged, pattern = "^MT-", 
                                                   assay = 'RNA')
tail(tak_merged@meta.data, 10)
# Visualize QC metrics as a violin plot
VlnPlot(tak_merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
        ncol = 3)

# QC Filtering
# The three batches vary somewhat. 1 and 2 show higher counts, features, and
# percent.mt than 3. Will filter for percent.mt expression less than 20%,
# counts < 20k, and 200 < features < 5 000.
tak_filtered <- subset(tak_merged, subset = nFeature_RNA < 5000 & 
                        nFeature_RNA > 200 & nCount_RNA < 20000 & 
                        percent.mt < 20)
tak_filtered 
# 24 712 features across 20 481 samples. Feature count remains same. 247 cells
# dropped.

# Normalize
# Next step is normalizing. In case there are any cell cycle effects, it may be
# helpful to find and then regress out cell cycle variability. To do so,
# normalize first, and then identify cell cycle phases before regressing it out
# along with percent.mt, library size, and batch.
# Normalize using NormalizeData with default values.
tak_normalized <- NormalizeData(tak_filtered)

# Cell Cycle Scoring
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
# We can divide this list into markers of G2/M phase and markers of S phase.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
tak_cellcyclephased <- CellCycleScoring(tak_normalized, s.features = s.genes, 
                                        g2m.features = g2m.genes, 
                                        set.ident = F)
# View cell cycle scores and phase assignments
head(tak_cellcyclephased@meta.data, 5)

# SCTransform
# Now we can run SCTransform() and regress out variability due to phase,
# percent.mt, batch, and library size.
tak_transform <- SCTransform(tak_cellcyclephased, assay = 'RNA', 
                             new.assay.name = 'SCT', 
                             vars.to.regress = c('orig.ident', 'percent.mt', 
                                                 'nFeature_RNA', 'nCount_RNA', 
                                                 'S.Score', 'G2M.Score'))
saveRDS(tak_transform, file = 'results/tak_ips_transform_wu_2018.rds')
