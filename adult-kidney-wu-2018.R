
# Load necessary libraries.
library(tidyverse)
library(Seurat)
library(Matrix)
library(patchwork)
library(sctransform)
library(patchwork)
library(openxlsx)
library(DESeq2)
library(gt)

# Read in processed data. Source of two replicates collected are from a
# 62-year-old white male with a serum creatinine of 1.03 mg/dL.
counts <- read.table(file = 'data/GSE118184_Human_kidney_snRNA.dge.txt', 
                         sep = '\t', header = T, row.names = 1)
glimpse(counts)
View(counts) # 20456 obs of 4524 cells

# Create Seurat Object. Min cells for gene to be expressed in is 10. Min genes
# per cell is 100.
ak <- CreateSeuratObject(counts = counts, min.cells = 10, 
                          min.features = 100, project = "ak")
ak
# Explore metadata. Includes orig.ident (Batch1, Batch2), RNA counts, and
# Feature counts.
head(ak@meta.data, 5)

# Find percentage of reads that map to mitochondrial genome:
ak[["percent.mt"]] <- PercentageFeatureSet(ak, pattern = "^MT-", 
                                            assay = 'RNA')
head(ak@meta.data, 10)
# Percentages are very low. Makes sense as data are from nuclei.
# Visualize QC metrics as a violin plot
VlnPlot(ak, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
# There is some variation in percent.mt, can regress it out. Filter out cells
# with over 1% percent.mt and over 4000 expressed features.
ak_filtered <- subset(ak, subset = nFeature_RNA < 4000 & percent.mt < 1)

# Next step is normalizing. In case there are any cell cycle effects, it may be
# helpful to find and then regress out cell cycle variability. To do so,
# normalize first, and then identify cell cycle phases before regressing it out
# along with percent.mt and library size.
# Normalize using NormalizeData with default values.
ak_normalized <- NormalizeData(ak_filtered)
# Cell Cycle Scoring
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ak_phase <- CellCycleScoring(ak_normalized, s.features = s.genes, 
                             g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(ak_phase@meta.data, 5)

# Now we can run SCTransform() and regress out variability due to phase,
# percent.mt, batch, and library size.
ak_transform <- SCTransform(ak_phase, assay = 'RNA', new.assay.name = 'SCT', 
                            vars.to.regress = c('orig.ident', 'percent.mt', 
                                                'nFeature_RNA', 'nCount_RNA', 
                                                'S.Score', 'G2M.Score'))

# Now we run PCA for linear dimensionality reduction.
ak_pca <- RunPCA(ak_transform)
# Examine and visualize PCA results a few different ways
print(ak_pca[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(ak_pca, dims = 1:2, reduction = "pca")
DimPlot(ak_pca, reduction = "pca")
DimHeatmap(ak_pca, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the 'dimensionality' of the dataset.
ElbowPlot(ak_pca) # At least 20 PCs.

# Cluster the Cells.
ak_neighbours <- FindNeighbors(ak_pca, dims = 1:30)
ak_clusters <- FindClusters(ak_neighbours) # Default resolution of 0.8
# Look at cluster IDs of the first 5 cells
head(Idents(ak_clusters), 5)

# Run non-linear dimensional reduction.
ak_umap <- RunUMAP(ak_clusters, dims = 1:30, umap.method = 'umap-learn', 
        metric = 'correlation')
DimPlot(ak_umap, reduction = "umap", label = T)

# Save the object for sharing and to avoid computationally intensive steps
# already conducted.
saveRDS(ak_umap, file = "results/ak_umap_wu_2018.rds")

# Finding differentially expressed features (cluster biomarkers).
# find markers for every cluster compared to all remaining cells, report only
# the positive ones
ak.markers <- FindAllMarkers(ak_umap, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25)
ak.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) %>% 
  select(cluster, gene, p_val_adj, avg_logFC, pct.1, pct.2, p_val) %>% 
  write.xlsx(file = 'results/ak_markers.xlsx', colNames = T, rowNames = F)

FeaturePlot(ak_umap, features = c("NPHS1", "NPHS2", "SLC5A2", "SLC22A6", "AGT",
                                  "SLC34A1", "LRP2", "SLC14A2", "CLCNKA", 
                                  "SLC12A1", "UMOD", "OXTR", "NOS1", "SLC12A3",
                                  "AQP2", "SCNN1G", "CALB1", "SLC4A1", 
                                  "SLC26A4", "AQP2", "SLC14A2"),
            pt.size = 0.2, ncol = 4)

DimPlot(ak_umap, reduction = "umap", label = T)

# Marker Genes: "NPHS1", "NPHS2", "SLC5A2", "SLC22A6", "AGT", "SLC34A1", "LRP2",
# "SLC14A2", "CLCNKA", "SLC12A1", "UMOD", "OXTR", "NOS1", "SLC12A3", "AQP2",
# "SCNN1G", "CALB1", "SLC4A1", "SLC26A4", "AQP2", "SLC14A2"
# OXTR not present apparently

# Podocyte: NPHS1, NPHS2
FeaturePlot(ak_umap, features = c("NPHS1")) # C9
FeaturePlot(ak_umap, features = c("NPHS2")) # C9
# Podocyte - Cluster 9
# Distal Convoluted Tubule: SLC12A3
FeaturePlot(ak_umap, features = c("SLC12A3")) # C5
# DCT - Cluster 5
# Thick Ascending Limb: SLC12A1, UMOD
FeaturePlot(ak_umap, features = c("SLC12A1", "UMOD")) # C2,4
# TAL - Clusters 2 & 4 # Find cTAL vs mTAL
# Type A Intercalated Cell: SLC4A1
FeaturePlot(ak_umap, features = "SLC4A1") # C7
# Type A IC - Cluster 7
# Type B Intercalated Cell: SLC26A4
FeaturePlot(ak_umap, features = "SLC26A4") # C11
# Type B IC - Cluster 11
# All Proximal Tubule: SLC34A1, LRP2
FeaturePlot(ak_umap, features = "SLC34A1") # C3,1,0
FeaturePlot(ak_umap, features = "LRP2") # C3,1,0,10 (parietal?)
# Proximal Tubule - Clusters 0,1,3

# Proximal S1: SLC5A2 (SGLT2)
FeaturePlot(ak_umap, features = c("SLC5A2")) # C0,3,1
# Proximal S2: SLC22A6 (OAT1)
FeaturePlot(ak_umap, features = c("SLC22A6")) #C3,1,0
# Proximal S3: AGT
FeaturePlot(ak_umap, features = "AGT") # C0,1,3 (tiny amount)
# Macula Densa (low abundance): SLC12A1, OXTR, NOS1
FeaturePlot(ak_umap, features = "SLC12A1") # C2,4
FeaturePlot(ak_umap, features = "OXTR") # Variable not found
FeaturePlot(ak_umap, features = "NOS1") # C4,2,5,6
# IMCD (Inner Medullary Collecting Duct): AQP2, SLC14A2 (UT-A1)
FeaturePlot(ak_umap, features = c("AQP2", "SLC14A2")) # Not notable


# Descending Thin Limb: SLC14A2 (UT-A2)
FeaturePlot(ak_umap, features = "SLC14A2") # C7,6,8,2,4,5
# Ascending Thin Limb: CLCNKA
FeaturePlot(ak_umap, features = "CLCNKA") # C7,2,4,5,6,8
# Connecting Tubule: AQP2, SCNN1G, CALB1
FeaturePlot(ak_umap, features = "AQP2") # C8
FeaturePlot(ak_umap, features = "SCNN1G") # C8,6
FeaturePlot(ak_umap, features = "CALB1") # C6,5
FeaturePlot(ak_umap, features = c("AQP2", "SCNN1G", "CALB1"))
FeaturePlot(ak_umap, features = "GATA3") # C8,6,5
# Principal Cell (Cortical Collecting Duct, Outer Medullary Collecting Duct):
# AQP2, SCNN1G
FeaturePlot(ak_umap, features = c("AQP2", "SCNN1G")) # C8,6

# Remaining Clusters: 6, 8, 10, 12, 13, 14, 15

# Distinguish CNT vs PC, cTAL vs mTAL. Parietal cell cluster?
# Compare three distal cell types to each other and PT control.

# Annotate known clusters
ak_umap <- readRDS("~/Desktop/Research/R/kidney-scRNA-seq-analysis/results/ak_umap_wu_2018.rds")
ak_labeled <- ak_umap
new.cluster.ids <- c('PT(0)', 'PT(1)', 'TAL(2)', 'PT(3)', 'TAL(4)', 'DCT', 'Connecting Tubule', 'IC-A', 'Principal Cells', 'Podocytes', 'Unknown', 'IC-B', 'Stressed Cells', 'EC/Mesangium (13)', 'Parietal Cells', 'EC/Mesangium (15)')
names(new.cluster.ids) <- levels(ak_labeled)
ak_labeled <- RenameIdents(ak_labeled, new.cluster.ids)
DimPlot(ak_labeled, reduction = "umap", label = TRUE, pt.size = 0.5)

# Save plots for marker genes. Probably best to save violin plots for their
# statistical information rather than feature plots. Perhaps less biased this
# way in interpretation too. Will save as pdfs for quality.
# Markers: "NPHS1", "NPHS2", "SLC5A2", "SLC22A6", "AGT", "SLC34A1", "LRP2",
# "SLC14A2", "CLCNKA", "SLC12A1", "UMOD", "OXTR", "NOS1", "SLC12A3", "AQP2",
# "SCNN1G", "CALB1", "SLC4A1", "SLC26A4", "AQP2", "SLC14A2", "SLC8A1", "SNTG1",
# "SLC7A1", "EYA4", "CLDN19", "PTH1R", "CLDN10", "SLC9A3", "ALDH1A2", "KIRREL3",
# "CFH", "TSHZ2"
# Podocytes
VlnPlot(ak_umap, features = c("NPHS1", 'NPHS2'))
# DCT
VlnPlot(ak_umap, features = "SLC12A3") + NoLegend()
# Thick AL
VlnPlot(ak_umap, features = c("SLC12A1", "UMOD"))
# IC-A
VlnPlot(ak_umap, features = 'SLC4A1') + NoLegend()
# IC-B
VlnPlot(ak_umap, features = "SLC26A4") + NoLegend()
# Proximal Tubule
VlnPlot(ak_umap, features = c('SLC34A1', 'LRP2'))
# Connecting Tubule (suggests cluster 6)
VlnPlot(ak_umap, features = c('CALB1', 'SLC8A1', 'SNTG1'))
# Collecting Duct (suggest cluster 8)
VlnPlot(ak_umap, features = c('SLC7A1', 'EYA4', 'AQP2', 'SCNN1G'), ncol = 2)
# Parietal Cells (suggests cluster 14)
VlnPlot(ak_umap, features = c('ALDH1A2', 'KIRREL3', 'CFH', 'TSHZ2'), ncol = 2)
# cTAL (unclear)
VlnPlot(ak_umap, features = c('CLDN19', 'PTH1R'))
# mTAL (unclear)
VlnPlot(ak_umap, features = c('CLDN10', 'SLC9A3', 'CLCNKA'))

# Proximal S1: SLC5A2 (SGLT2)
VlnPlot(ak_umap, features = "SLC5A2") + NoLegend() # none
# Proximal S2: SLC22A6 (OAT1)
VlnPlot(ak_umap, features = "SLC22A6") + NoLegend() # 0,1,3,12
# Proximal S3: AGT
VlnPlot(ak_umap, features = "AGT") + NoLegend() # none
# Macula Densa (low abundance): SLC12A1, OXTR, NOS1
VlnPlot(ak_umap, features = "SLC12A1") + NoLegend() # 
VlnPlot(ak_umap, features = "OXTR") + NoLegend() # Variable not found
VlnPlot(ak_umap, features = "NOS1") + NoLegend() # 2,4,5
# IMCD (Inner Medullary Collecting Duct): AQP2, SLC14A2 (UT-A1)
VlnPlot(ak_umap, features = c("AQP2", "SLC14A2")) # Not notable
# Descending Thin Limb: SLC14A2 (UT-A2)
FeaturePlot(ak_umap, features = "SLC14A2") # C7,6,8,2,4,5
# Ascending Thin Limb: CLCNKA
FeaturePlot(ak_umap, features = "CLCNKA") # C7,2,4,5,6,8
# Connecting Tubule: AQP2, SCNN1G, CALB1
VlnPlot(ak_umap, features = "AQP2") # C8
VlnPlot(ak_umap, features = "SCNN1G") # C8,6
VlnPlot(ak_umap, features = "CALB1") # C6,5
VlnPlot(ak_umap, features = "GATA3") # C8,6,5

# Endothelial Cells: EMCN
VlnPlot(ak_umap, features = 'EMCN') + NoLegend() # 13, 15
# Mesangium: ITGA8
VlnPlot(ak_umap, features = 'ITGA8') + NoLegend() # 13, 15
# Macrophage: PTPRC
VlnPlot(ak_umap, features = 'PTPRC') # Nothing
# Cluster 10
VlnPlot(ak_umap, features = "TMEM178B") # 10, 14, 3

# Compare both TAL clusters together:
cluster2v4.markers <- FindMarkers(ak_umap, ident.1 = 2, ident.2 = 4, min.pct = 0.25, logfc.threshold = 0.2, only.pos = F)
View(cluster2v4.markers)

# Compare PC and CNT clusters together:
cluster6v8.markers <- FindMarkers(ak_umap, ident.1 = 6, ident.2 = 8, min.pct = 0.25, logfc.threshold = 0.2, only.pos = F)
View(cluster6v8.markers) # Seems reasonable to conclude 6 as CNT and 8 as PC

# Compare PT clusters together:
cluster0v1and3.markers <- FindMarkers(ak_umap, ident.1 = 0, ident.2 = c(1, 3), min.pct = 0.25, logfc.threshold = 0.2, only.pos = F)
cluster1v0and3.markers <- FindMarkers(ak_umap, ident.1 = 1, ident.2 = c(0, 3), min.pct = 0.25, logfc.threshold = 0.2, only.pos = F)
cluster3v0and1.markers <- FindMarkers(ak_umap, ident.1 = 3, ident.2 = c(0, 1), min.pct = 0.25, logfc.threshold = 0.2, only.pos = F)
View(cluster0v1and3.markers)
View(cluster1v0and3.markers)
View(cluster3v0and1.markers)

# Will now try to extract the three distal cell type clusters of interest for
# comparison and plotting. Subset ak_umap to include just the three clusters.
ak_distal <- subset(ak_labeled, idents = c('TAL(2)', 'TAL(4)', 'DCT', 
                                           'Connecting Tubule'))
DimPlot(ak_distal, reduction = 'umap', label = T, pt.size = 0.5)
distal.idents <- c('TAL', 'TAL', 'DCT', 'CNT')
names(distal.idents) <- levels(ak_distal)
ak_distal <- RenameIdents(ak_distal, distal.idents)
DimPlot(ak_distal, reduction = "umap", pt.size = 0.5)
# Now have a Seurat object with only DCT, CNT, and TAL clusters. Save object.
saveRDS(ak_distal, 'results/ak_distal.rds')

# Now generate heatmaps of gene expression across the three clusters.
# Install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# Use differential expression testing to find most variable markers.
# For the three clusters compared to all other clusters.
# Begin by altering ak_labeled annotations to have one TAL cluster.
new.cluster.ids <- c('PT(0)', 'PT(1)', 'TAL', 'PT(3)', 'TAL', 'DCT', 'CNT', 
                     'IC-A', 'Principal Cells', 'Podocytes', 'Unknown', 'IC-B', 'Stressed Cells', 'EC/Mesangium (13)', 'Parietal Cells', 
                     'EC/Mesangium (15)')
names(new.cluster.ids) <- levels(ak_labeled)
ak_labeled <- RenameIdents(ak_labeled, new.cluster.ids)
DimPlot(ak_labeled, reduction = "umap", label = TRUE, pt.size = 0.5)
# Find all markers using DESeq2 method.
ak_labeled.deseq2.markers <- FindAllMarkers(ak_labeled, only.pos = T, 
                                            min.pct = 0.25, 
                                            logfc.threshold = 0.25,
                                            test.use = "DESeq2")
# Look at the results.
ak_labeled.deseq2.markers.top20 <- ak_labeled.deseq2.markers %>% 
  group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
View(ak_labeled.deseq2.markers.top20)
ak_labeled.deseq2.markers.top20.distal <- ak_labeled.deseq2.markers.top20 %>% 
  filter(cluster %in% c('TAL', 'DCT', 'CNT'))
View(ak_labeled.deseq2.markers.top20.distal)
# Generate heatmaps of most differentially expressed features in the three
# clusters. Yellow = high, purple = low.
# Heatmap for all clusters and top 20 features.
DoHeatmap(ak_labeled, features = ak_labeled.deseq2.markers.top20$gene)
# Neat to see, but very cluttered with too many genes.
# Heatmap for only the three clusters. Saved as pdf.
DoHeatmap(ak_distal, features = ak_labeled.deseq2.markers.top20.distal$gene, 
          angle = 0, hjust = 0.5, raster = F) + NoLegend()
# Now need to generate heatmap of distal cell types with most differentially
# expressed genes between the three.
# Do so by finding markers only across the distal clusters.
ak_distal.deseq2.markers <- FindAllMarkers(ak_distal, only.pos = T, 
                                           min.pct = 0.25, 
                                           logfc.threshold = 0.25,
                                           test.use = "DESeq2")
# Look at results.
ak_distal.deseq2.markers.top20 <- ak_distal.deseq2.markers %>% 
  group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
View(ak_distal.deseq2.markers.top20)
# Generate heatmap for the three clusters. Export as pdf.
DoHeatmap(ak_distal, features = ak_distal.deseq2.markers.top20$gene, angle = 0, 
          hjust = 0.5, raster = F) + NoLegend()
# Question: Would this vary if all genes were included, and not just the top
# 3000? I don't think it would as for heatmaps, all genes must be scaled and as
# all genes found with FindAllMarkers which starts with the raw counts (I
# believe) are present in the heatmap, there should not be an issue here.
# Another question: What differences are there between the top 20 markers
# identified for the distal clusters when they are compared to all other
# clusters and when compared to only each other?
# Will export top 20 markers to excel and depict top genes in table format.
ak_labeled.deseq2.markers.top20.distal
ak_distal.deseq2.markers.top20
ak_deseq2.markers.top20 <- list('distalvsall' = 
                                  ak_labeled.deseq2.markers.top20.distal, 
                                'distalvsdistal' = 
                                  ak_distal.deseq2.markers.top20)
write.xlsx(ak_deseq2.markers.top20, file = 'results/ak_deseq2markers.xlsx', 
           row.names= F)
# Make gene tables.
distalvsall <- ak_labeled.deseq2.markers.top20.distal %>% 
  select(cluster, gene) %>% 
  pivot_wider(names_from = cluster, values_from = gene, values_fn = list) %>% 
  unnest(cols = c(TAL, DCT, CNT)) %>% 
  as_tibble()
distalvsdistal <- ak_distal.deseq2.markers.top20 %>% 
  select(cluster, gene) %>% 
  pivot_wider(names_from = cluster, values_from = gene, values_fn = list) %>% 
  unnest(cols = c(TAL, DCT, CNT)) %>% 
  as_tibble()
cbind(distalvsall, distalvsdistal) %>% repair_names() %>% gt() %>% 
  tab_header(title = md('**Top 20 Markers per Cell Cluster**')) %>% 
  tab_spanner(label = 'Compared to All Cells', columns = vars(TAL, DCT, CNT)) %>% 
  tab_spanner(label = 'Compared Only to Each Other', columns = vars(TAL1, DCT1, CNT1)) %>% 
  cols_label(TAL1 = 'TAL', DCT1 = 'DCT', CNT1 = 'CNT') %>% 
  tab_style(style = cell_fill(color = 'lightcyan'), 
            locations = list(cells_body(columns = vars(TAL), rows = 18),
                             cells_body(columns = vars(DCT), rows = 17),
                             cells_body(columns = vars(DCT), rows = 10),
                             cells_body(columns = vars(CNT), rows = 14),
                             cells_body(columns = vars(DCT), rows = 19),
                             cells_body(columns = vars(CNT), rows = 17),
                             cells_body(columns = vars(DCT), rows = 20),
                             cells_body(columns = vars(CNT), rows = 19)
                             )) %>% 
  tab_style(style = cell_fill(color = 'lemonchiffon'), 
            locations = list(cells_body(columns = vars(TAL1), rows = 1:7),
                             cells_body(columns = vars(TAL1), rows = 9:12),
                             cells_body(columns = vars(TAL1), rows = 17),
                             cells_body(columns = vars(TAL1), rows = 19:20),
                             cells_body(columns = vars(DCT1), rows = 1:8),
                             cells_body(columns = vars(DCT1), rows = 12),
                             cells_body(columns = vars(DCT1), rows = 15),
                             cells_body(columns = vars(DCT1), rows = 18),
                             cells_body(columns = vars(DCT1), rows = 20),
                             cells_body(columns = vars(CNT1), rows = 2:3),
                             cells_body(columns = vars(CNT1), rows = 5:11),
                             cells_body(columns = vars(CNT1), rows = 13:16),
                             cells_body(columns = vars(CNT1), rows = 18:20)
                             )) %>% 
  tab_source_note('Cyan highlights markers expressed by two different clusters.') %>% 
  tab_source_note('Yellow highlights markers expressed by a cluster across both comparison tests.')
# Colored cells. Cyan for overlaps. Yellow for constants. Saved as html file.




