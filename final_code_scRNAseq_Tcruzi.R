#### Initial processing and QC ####


# Install and load Seurat
#install.packages("Seurat")

library(Seurat)
library(Matrix)
library(dplyr)
library("vembedr")
library("htmltools")
library("Matrix")
library("irlba")
library("ggplot2")
library("dplyr")
library("scico")
library('fastmap')
library("Seurat")
library("scater")
library("SingleCellExperiment")


# Load .mtx, genes and barcodes files
mat_1 <- Matrix::readMM(file = 'cells_x_genes.mtx')
features_1 <- read.delim( "cells_x_genes.genes.txt", header = FALSE)
barcodes_1 <- read.delim("cells_x_genes.barcodes.txt", header = FALSE)

mat_2 <- Matrix::readMM(file = 'cells_x_genes_2.mtx')
features_2 <- read.delim( "cells_x_genes.genes_2.txt", header = FALSE)
barcodes_2 <- read.delim("cells_x_genes.barcodes_2.txt", header = FALSE)


# Create a sparse matrix for more efficient computation
mat_1 <- as(mat_1, "CsparseMatrix")
mat_2 <- as(mat_2, "CsparseMatrix")


# Assign cell names and features to the matrix
rownames(mat_1) <- barcodes_1$V1  
colnames(mat_1) <- features_1$V1  

rownames(mat_2) <- barcodes_2$V1  
colnames(mat_2) <- features_2$V1  


# Find common barcodes
common_barcodes <- intersect(barcodes_1$V1, barcodes_2$V1)

# Reorder mat_1 and mat_2 to include only common barcodes and in the same order
mat_1_common <- mat_1[match(common_barcodes, rownames(mat_1)), ]
mat_2_common <- mat_2[match(common_barcodes, rownames(mat_2)), ]

# Ensure that the columns (features) are the same and in the same order
if (!identical(colnames(mat_1), colnames(mat_2))) {
  stop("Las columnas de mat_1 y mat_2 no coinciden o no están en el mismo orden.")
}

# Sum the matrices
mat_3 <- mat_1_common + mat_2_common


mat <- mat_3
barcodes <- rownames(mat)
features <- colnames(mat)


#### METADATA ####


# Create metadata containing only the cell IDs
metadata <- data.frame(row.names = barcodes, cells = barcodes, stringsAsFactors = F)


# Add number of UMIs per cell to metadata
metadata$nUMI <- Matrix::rowSums(mat)

# Add number of genes detected per cell to metadata
metadata$nGene <- Matrix::rowSums(mat > 0)

# Add log10(Genes per UMI) to metadata
metadata$log10GenesPerUMI <- log10(metadata$nGene) / log10(metadata$nUMI)


#### ANNOTATION ####

library(tidyverse)
annotations <- read.csv('anotacion_dm28c_con_mt.gff', sep = '\t', header = FALSE)
annotations <- annotations %>% filter(V3 == 'transcript')
annotations <- as.data.frame(annotations$V9)


annotations <- annotations %>% 
  separate(`annotations$V9`, into = c("gene_id", "ID", "Description", 'info', 'transcript_id', 'parent', 'ebi'), sep = ";")


library(stringr)
annotations <- annotations %>%
  mutate(gene_id = str_replace(gene_id, "gene_id", ""),
         Description = str_replace(Description, "description", ""))

annotations$gene_id <- trimws(annotations$gene_id) # to remove spaces

#### MITOCHONDRIAL GENES ####


mt <- c("12S", "9S", "ND8", "ND9", "MURF5", "ND7", "COIII", "Cyb", "A6_(MURF4)", "MURF1", "CR3", "ND1", "COII", "MURF2", "COI", "CR4", "ND4", "ND3", "RPS12", "ND5")

length(intersect(colnames(mat), mt)) # 20

# Number of UMIs assigned to mitochondrial genes
metadata$mtUMI <- Matrix::rowSums(mat[, which(colnames(mat) %in% mt)], na.rm = T)

# Ensure all NAs receive zero counts
metadata$mtUMI[is.na(metadata$mtUMI)] <- 0

# Calculate of mitoRatio per cell
metadata$mitoRatio <- metadata$mtUMI/metadata$nUMI


summary(metadata$nUMI)
# Min. 1st Qu.  Median    Mean   3rd Qu.    Max. 
# 0     873     1224      1521    1759     198207 



# Keep cells with nUMI greater than 100
#idx <- which(metadata$nUMI > 100)

# Extract the counts for those cells
#mat <- mat[idx, ]

# Extract the metadata for those cells
#metadata <- metadata[idx,]

dim(mat) # 82667 19132


# Save data to single cell experiment variable
se <- SingleCellExperiment(assays=list(counts=t(mat)), 
                           colData = metadata)






# Create a data frame containing the metrics for visualizations
metrics <- colData(se) %>%
  as.data.frame


# Visualize the number UMIs/transcripts per cell
metrics %>% 
  ggplot(aes(color='red', x=nUMI, fill = 'red')) + 
  geom_density() + 
  scale_x_log10() + 
  ylab("log10 cell density") 


# Visualize the distribution of genes detected per cell via histogram
metrics %>% 
  ggplot(aes(color='red', x=nGene, fill= 'red')) + 
  geom_density() + 
  scale_x_log10()+ 
  geom_vline(xintercept = 130) 

# Visualize the distribution of genes detected per cell via boxplot
metrics %>% 
  ggplot(aes(x='red', y=log10(nGene), fill='red')) + 
  geom_boxplot() + 
  ggtitle("NCells vs NGenes")



# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metrics %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  geom_vline(xintercept = 800) 




# Visualize the distribution of mitochondrial gene expression detected per cell
metrics %>% 
  ggplot(aes(color='red', x=mitoRatio, fill='red')) + 
  geom_density() + 
  scale_x_log10() + 
  geom_vline(xintercept = 0.1)




#### FILTERING ####

summary(metadata$nGene)

summary(metadata$nUMI)





# Sort cells by nUMI in decreasing order
sorted_nUMI <- sort(metadata$nUMI, decreasing = TRUE)

# Create a vector of indices for the x axis
x_vals <- 1:length(sorted_nUMI)

# Generate the knee plot in log-log scale
ggplot() + 
  geom_line(aes(x = x_vals, y = sorted_nUMI), color = "red") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Cells (ordered by decreasing nUMI)") +
  ylab("nUMI") +
  ggtitle("Knee Plot (similar to cellranger)")# +



# Filtering criteria:

# nUMI > 3000
# nGene > 11
# log10GenesPerUMI > 0.8
# mitoRatio < 0.1


# Filter out low quality reads using selected thresholds - these will change with experiment
keep <- metrics %>%
  dplyr::filter(nUMI > 3000 , 
                nGene > 130,
                log10GenesPerUMI > 0.8,
                mitoRatio < 0.1,
  ) %>% 
  pull(cells)

# Subset the cells to only include those that meet the thresholds specified
se_c <- se[ ,keep]

# Save subset to new metrics variable
metrics_clean <- colData(se_c) %>% as.data.frame()





# Visualize the number UMIs/transcripts per cell
metrics_clean %>% 
  ggplot(aes(color='red', x=nUMI, fill = 'red')) + 
  geom_density() + 
  scale_x_log10() + 
  ylab("log10 cell density")


# Visualize the distribution of genes detected per cell via histogram
metrics_clean %>% 
  ggplot(aes(color='red', x=nGene, fill= 'red')) + 
  geom_density() + 
  scale_x_log10()

# Visualize the distribution of genes detected per cell via boxplot
metrics_clean %>% 
  ggplot(aes(x='red', y=log10(nGene), fill='red')) + 
  geom_boxplot() + 
  ggtitle("NCells vs NGenes")



# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metrics_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  geom_vline(xintercept = 800) 




# Visualize the distribution of mitochondrial gene expression detected per cell
metrics_clean %>% 
  ggplot(aes(color='red', x=mitoRatio, fill='red')) + 
  geom_density() + 
  scale_x_log10() + 
  geom_vline(xintercept = 0.1)







counts <- assays(se_c)$counts

# Filter the matrix to exclude those genes
counts_filtrados <- counts[!(rownames(counts) %in% mt), ]

# Create a Seurat object from the data
filtered_Seurat <- CreateSeuratObject(counts = counts_filtrados)

# Normalize the data
filtered_Seurat = NormalizeData(filtered_Seurat)

# Identify variable features (genes) in the data
filtered_Seurat = FindVariableFeatures(filtered_Seurat)

# Scale the data
filtered_Seurat = ScaleData(filtered_Seurat)

# Perform a principal component analysis (PCA)
filtered_Seurat = RunPCA(filtered_Seurat)

# Perform the JackStraw test to determine significance of principal components
#filtered_Seurat <- JackStraw(filtered_Seurat, num.replicate = 100)
#filtered_Seurat <- ScoreJackStraw(filtered_Seurat, dims = 1:20)

# Visualize the JackStraw test results
#JackStrawPlot(filtered_Seurat, dims = 1:15)


# To visualize and decide how many PCs to use, 10
ElbowPlot(filtered_Seurat)  

# Find nearest neighbors based on PCA dimensions
filtered_Seurat = FindNeighbors(filtered_Seurat, dims = 1:10)

# Identify cell clusters
filtered_Seurat = FindClusters(filtered_Seurat, resolution = 0.2)

# Perform UMAP analysis for 2D visualization
filtered_Seurat = RunUMAP(filtered_Seurat, dims = 1:10)

# Visualize clusters in a 2D plot
DimPlot(filtered_Seurat, pt.size = 1)






#### FILTERING MITOCHONDRIAL AND RIBOSOMAL GENES ####

# Load additional libraries
library(scater)
library(SingleCellExperiment)

# Calculate per-cell quality metrics
sceMetrics = perCellQCMetrics(se_c)

# Repeat the process for data without ribosomal RNA genes

counts <- as.data.frame(se_c@assays@data@listData)

# Read a list of ribosomal RNA genes
ribosomal_rna = read.csv("Ribsoomal_RNA.csv", header = FALSE)
ribosomal_rna = as.character(ribosomal_rna$V1)

# Replace "-" with "_"
ribosomal_rna <- gsub("-", "_", ribosomal_rna)

# Filter ribosomal RNA genes from the data

counts_sinrRNA = counts[!rownames(counts)%in%ribosomal_rna, ]

counts_sinRNA_sinmt = counts_sinrRNA[!rownames(counts_sinrRNA)%in%mt, ]

counts_sinRNA_sinmt2 = as.matrix(counts_sinRNA_sinmt)
sce_sin_ribosomal_sin_mt <- SingleCellExperiment(list(counts=counts_sinRNA_sinmt2))
sceMetrics = perCellQCMetrics(sce_sin_ribosomal_sin_mt)


counts <- assays(sce_sin_ribosomal_sin_mt)$counts


# Create a Seurat object from the data
filtered_Seurat <- CreateSeuratObject(counts = counts)

# Normalize the data
filtered_Seurat = NormalizeData(filtered_Seurat)

# Identify variable features (genes) in the data
filtered_Seurat = FindVariableFeatures(filtered_Seurat)

# Scale the data
filtered_Seurat = ScaleData(filtered_Seurat)

# Perform a principal component analysis (PCA)
filtered_Seurat = RunPCA(filtered_Seurat)

# Perform the JackStraw test to determine the significance of the principal components
#filtered_Seurat <- JackStraw(filtered_Seurat, num.replicate = 100)
#filtered_Seurat <- ScoreJackStraw(filtered_Seurat, dims = 1:20)

# Visualize the JackStraw test results
#JackStrawPlot(filtered_Seurat, dims = 1:15)


# To visualize and decide how many PCs to use, 10
ElbowPlot(filtered_Seurat)  

# Find nearest neighbors based on PCA dimensions
filtered_Seurat = FindNeighbors(filtered_Seurat, dims = 1:10)

# Identify cell clusters
filtered_Seurat = FindClusters(filtered_Seurat, resolution = 0.3)

# Perform UMAP analysis for 2D visualization
filtered_Seurat = RunUMAP(filtered_Seurat, dims = 1:10, n.neighbors = 50, min.dist = 0.001)

# plot UMAP to identify possible doublets
metrics <-  c("nCount_RNA", "nFeature_RNA")
FeaturePlot(filtered_Seurat, 
            reduction = "umap", 
            features = metrics,
            pt.size = 1, 
            order = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            label = TRUE)

DimPlot(filtered_Seurat, pt.size = 2)


summary(filtered_Seurat@meta.data)
#orig.ident           nCount_RNA     nFeature_RNA  RNA_snn_res.0.3 seurat_clusters
#SeuratProject:3192   Min.   : 1012   Min.   : 582   0:2201          0:2201         
#                     1st Qu.: 1391   1st Qu.: 739   1: 824          1: 824         
#                     Median : 1722   Median : 883   2: 167          2: 167         
#                     Mean   : 2461   Mean   :1088                                  
#                     3rd Qu.: 2358   3rd Qu.:1148                                  
#                     Max.   :39479   Max.   :6726     







# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

# pK identification (no ground-truth)
sweep.list <- paramSweep(filtered_Seurat, PCs = 1:10)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic doublet proportion estimate
annotations_2 <- filtered_Seurat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations_2) 
nExp.poi <- round(optimal.pk * nrow(filtered_Seurat@meta.data)) ## Assuming 1% doublet formation rate - tailor for your dataset
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

# run DoubletFinder, 
filtered_Seurat <- doubletFinder(seu = filtered_Seurat, 
                                 PCs = 1:10, 
                                 pK = optimal.pk,
                                 nExp = nExp.poi.adj)

metadata <- filtered_Seurat@meta.data
colnames(metadata)[7] <- "doublet_finder"

filtered_Seurat@meta.data <- metadata 


DimPlot(filtered_Seurat, group.by = "doublet_finder")

table(metadata$doublet_finder) # 41 Doublets, 3151 Singlets

# Remove doublets identified with DoubletFinder
# subset and save
filtered_Seurat <- subset(filtered_Seurat, doublet_finder == "Singlet")






# Identify markers for each cluster
markers = FindAllMarkers(filtered_Seurat)

markers_sorted <- markers %>% arrange(avg_log2FC)


# Create a feature scatter plot
plot2 <- FeatureScatter(filtered_Seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

# Identify the top 10 variable features
top10 <- head(VariableFeatures(filtered_Seurat), 10)
plot1 <- VariableFeaturePlot(filtered_Seurat)
plot1


# Filter rows where avg_log2FC is greater than 0
markers_sorted_FCpos <- markers_sorted[markers_sorted$avg_log2FC > 1, ]


markers_cluster0 <- markers_sorted_FCpos[markers_sorted_FCpos$cluster == 0, ]
markers_cluster1 <- markers_sorted_FCpos[markers_sorted_FCpos$cluster == 1, ]
markers_cluster2 <- markers_sorted_FCpos[markers_sorted_FCpos$cluster == 2, ]




write.csv(markers_cluster0, file = "markers_Tcruzi_cluster0.csv", row.names = TRUE)
write.csv(markers_cluster1, file = "markers_Tcruzi_cluster1.csv", row.names = TRUE)
write.csv(markers_cluster2, file = "markers_Tcruzi_cluster2.csv", row.names = TRUE)













#### Figure 1: Identification of amastigote and trypomastigote cell populations ####

library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(data.table)
library(gridExtra)

# ============================================================================
# FIGURE 1: Identification of amastigote and trypomastigote cell populations
# ============================================================================
# 1a) UMAP colored by detected clusters
# 1b) Heatmap of the top 10 gene markers upregulated in each cluster
# 1c) Expression of cluster 0 marker (C4B63-16g183 - GP63)
# 1d) Expression of cluster 1 marker (C4B63-16g155 - dehydrogenase)

#### FIGURE 1a - UMAP ####
png("Figure1a_umap_clusters.png", width = 10, height = 6, units = "cm", res = 300)
DimPlot(filtered_Seurat, pt.size = 0.5) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme_minimal()
dev.off()

#### FIGURE 1b - Heatmap top 10 markers ####
markers <- FindAllMarkers(filtered_Seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10_markers <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

# Obtain expression matrix
expression_matrix <- as.matrix(GetAssayData(filtered_Seurat, slot = "data"))
filtered_expression_matrix <- expression_matrix[rownames(expression_matrix) %in% top10_markers$gene, ]

cell_clusters <- data.frame(
  cell = colnames(filtered_Seurat),
  cluster = as.factor(filtered_Seurat@meta.data$seurat_clusters)
)
cell_clusters_ordered <- cell_clusters[order(cell_clusters$cluster), ]
ordered_expression_matrix <- filtered_expression_matrix[, cell_clusters_ordered$cell]

# Colors for clusters
cluster_colors <- c("salmon", "green3", "skyblue")
names(cluster_colors) <- levels(cell_clusters_ordered$cluster)

top_annotation <- HeatmapAnnotation(
  cluster = cell_clusters_ordered$cluster, 
  col = list(cluster = cluster_colors)
)

# Add gene annotations
annotations_2 <- annotations
annotations_2$Genes <- gsub("_", "-", annotations_2$gene_id)
top10_markers <- top10_markers %>%
  left_join(annotations_2, by = c("gene" = "Genes"))

gene_labels <- paste(top10_markers$gene, ":", top10_markers$Description)

col_fun <- colorRamp2(c(-2, 0, 2), c("magenta", "black", "yellow"))
scaled_matrix <- t(scale(t(ordered_expression_matrix)))

row_order <- c(
  top10_markers$gene[top10_markers$cluster == 0],
  top10_markers$gene[top10_markers$cluster == 1],
  top10_markers$gene[top10_markers$cluster == 2]
)
row_split <- factor(rep(0:2, each = 10))
scaled_matrix <- scaled_matrix[row_order, ]
rownames(scaled_matrix) <- gene_labels

heatmap_genes <- Heatmap(scaled_matrix,
                         name = "Expression",
                         col = col_fun,
                         use_raster = TRUE,
                         show_row_names = TRUE,
                         row_names_gp = gpar(fontsize = 12),
                         row_names_max_width = unit(12, "cm"),
                         show_column_names = FALSE,
                         top_annotation = top_annotation,
                         cluster_rows = FALSE,
                         cluster_columns = FALSE,
                         row_split = row_split,
                         column_split = cell_clusters_ordered$cluster,
                         row_gap = unit(1, "mm"),
                         column_gap = unit(0.5, "mm"),
                         row_dend_width = unit(3, "cm"))

png("Figure1b_heatmap_top10markers.png", width = 3000, height = 3500, res = 300)
draw(heatmap_genes, heatmap_legend_side = "right", show_annotation_legend = FALSE)
dev.off()

#### FIGURE 1c - FeaturePlot GP63 (cluster 0 marker) ####
png("Figure1c_marker_GP63_cluster0.png", width = 10, height = 6, units = "cm", res = 300)
FeaturePlot(filtered_Seurat, features = "C4B63-16g183", 
            min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5) +
  ggtitle("") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme_minimal()
dev.off()

#### FIGURE 1d - FeaturePlot dehydrogenase (cluster 1 marker) ####
png("Figure1d_marker_dehydrogenase_cluster1.png", width = 10, height = 6, units = "cm", res = 300)
FeaturePlot(filtered_Seurat, features = "C4B63-16g155", 
            min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5) +
  ggtitle("") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme_minimal()
dev.off()









#### Figure 2: Multigene family expression and distribution ####

# 2a) Violin plot - summatory of multigene family expression
# 2b) UMAP - multigene family genes detected
# 2c) Boxplot - number of cells expressing each gene (core vs multigene families vs ribosomal)
# 2d) Lorenz curves and Gini index


# keep ribosomal proteins

proteinas_ribosomales <- annotations[grepl("ribosome|ribosomal", annotations$Description, ignore.case = TRUE), "gene_id"]
proteinas_ribosomales <- gsub("_", "-", proteinas_ribosomales)

# Get cells by cluster
tripomastigotas <- rownames(filtered_Seurat@meta.data)[filtered_Seurat@meta.data$seurat_clusters == 0]
amastigotas <- rownames(filtered_Seurat@meta.data)[filtered_Seurat@meta.data$seurat_clusters == 1]

# Subset the expression matrix keeping trypomastigote cells

# Grab the expression matrix from the Seurat object 
expression_matrix <- as.data.frame(GetAssayData(object = filtered_Seurat))
dim(expression_matrix)

expression_matrix_tripos <- expression_matrix[, tripomastigotas]
dim(expression_matrix_tripos)


# Subset the trypo expression matrix separating core and disruptive genes

# Expression matrix for core genes
expression_genes_core_tripos <- expression_matrix_tripos[genes_core, ]

# Expression matrix for disruptive genes
expression_genes_disruptive_tripos <- expression_matrix_tripos[genes_disruptive, ]

# Expression matrix for ribosomal genes
expression_genes_ribosomales_tripos <- expression_matrix_tripos[proteinas_ribosomales, ]




# Generate lists of genes for core and disruptive compartments

# basically disruptive are: transialidases, MASP, mucins, dgf-1, gp63 and RHS

# Define terms to search for
terminos <- c("transialidase", "Transialidase","trans-sialidase", "Trans-sialidase", "mucin", "Mucin", "masp", "MASP", 
              "GP63", "gp63", "gp-63", "GP-63", "RHS", "rhs", "DGF1", "dgf1", 
              "DGF-1", "dgf-1'")

# Create the regular expression
patron <- paste(terminos, collapse = "|")
patron <- paste0("\\b(", patron, ")\\b")

# Search for matches in the 'Description' column
genes_disruptive <- annotations$gene_id[grep(patron, annotations$Description, ignore.case = TRUE)]

# Gene_id that are not in genes_disruptive
genes_core <- setdiff(annotations$gene_id, genes_disruptive)

# Replace "-" with "_"
genes_core <- gsub("_", "-", genes_core)
# Replace "-" with "_"
genes_disruptive <- gsub("_", "-", genes_disruptive)









# Now it would be good to see disruptive compartment separated

terminos_mucinas <- c("TcMUC")
terminos_masp <- c("masp", "MASP")
terminos_transialidasas <- c("transialidase", "Transialidase","trans-sialidase", "Trans-sialidase")
terminos_rhs <- c("RHS", "rhs")
terminos_gp63 <- c("GP63", "gp63", "gp-63", "GP-63")
terminos_dgf1 <- c("DGF1", "dgf1", "DGF-1", "dgf-1'")


# Create the regular expressions
patron_mucinas <- terminos_mucinas

patron_masp <- paste(terminos_masp, collapse = "|")
patron_masp <- paste0("\\b(", patron_masp, ")\\b")

patron_trans <- paste(terminos_transialidasas, collapse = "|")
patron_trans <- paste0("\\b(", patron_trans, ")\\b")

patron_rhs <- paste(terminos_rhs, collapse = "|")
patron_rhs <- paste0("\\b(", patron_rhs, ")\\b")

patron_gp63 <- paste(terminos_gp63, collapse = "|")
patron_gp63 <- paste0("\\b(", patron_gp63, ")\\b")

patron_dgf1 <- paste(terminos_dgf1, collapse = "|")
patron_dgf1 <- paste0("\\b(", patron_dgf1, ")\\b")



# Search for matches in the 'Description' column
genes_mucinas <- annotations$gene_id[grep(patron_mucinas, annotations$Description, ignore.case = TRUE)]
genes_masp <- annotations$gene_id[grep(patron_masp, annotations$Description, ignore.case = TRUE)]
genes_transialidasas <- annotations$gene_id[grep(patron_trans, annotations$Description, ignore.case = TRUE)]
genes_rhs <- annotations$gene_id[grep(patron_rhs, annotations$Description, ignore.case = TRUE)]
genes_gp63 <- annotations$gene_id[grep(patron_gp63, annotations$Description, ignore.case = TRUE)]
genes_dgf1 <- annotations$gene_id[grep(patron_dgf1, annotations$Description, ignore.case = TRUE)]


# Generate annotation matrices for each multigene family

anotacion_mucinas <- annotations[annotations$gene_id %in% genes_mucinas, ]
anotacion_masp <- annotations[annotations$gene_id %in% genes_masp, ]
anotacion_transialidasas <- annotations[annotations$gene_id %in% genes_transialidasas, ]
anotacion_rhs <- annotations[annotations$gene_id %in% genes_rhs, ]
anotacion_gp63 <- annotations[annotations$gene_id %in% genes_gp63, ]
anotacion_dgf1 <- annotations[annotations$gene_id %in% genes_dgf1, ]


# Replace "-" with "_"
genes_mucinas <- gsub("_", "-", genes_mucinas)
genes_masp <- gsub("_", "-", genes_masp)
genes_transialidasas <- gsub("_", "-", genes_transialidasas)
genes_rhs <- gsub("_", "-", genes_rhs)
genes_gp63 <- gsub("_", "-", genes_gp63)
genes_dgf1 <- gsub("_", "-", genes_dgf1)




# Subset the expression matrix keeping trypomastigote cells

# Grab the expression matrix from the Seurat object 
expression_matrix <- as.data.frame(GetAssayData(object = filtered_Seurat))
dim(expression_matrix)


tripomastigotas <- WhichCells(filtered_Seurat, ident = "0")
expression_matrix_tripos <- expression_matrix[, tripomastigotas]
dim(expression_matrix_tripos)


# Cells from cluster 1 would be the amastigotes
amastigotas <- names(filtered_Seurat@active.ident[filtered_Seurat@active.ident == 1])

# Subset the expression matrix keeping amastigote cells
dim(expression_matrix)
expression_matrix_amas <- expression_matrix[, amastigotas]
dim(expression_matrix_amas)






# Subset the trypo expression matrix separating core and disruptive genes

# Expression matrix for core genes
expression_genes_core_tripos <- expression_matrix_tripos[genes_core, ]

# Expression matrix for disruptive genes
expression_genes_disruptive_tripos <- expression_matrix_tripos[genes_disruptive, ]

# Expression matrix for ribosomal genes
expression_genes_ribosomales_tripos <- expression_matrix_tripos[proteinas_ribosomales, ]




# Expression matrix for multigene families
expression_genes_mucinas_tripos <- expression_matrix_tripos[genes_mucinas, ]
expression_genes_masp_tripos <- expression_matrix_tripos[genes_masp, ]
expression_genes_trans_tripos <- expression_matrix_tripos[genes_transialidasas, ]
expression_genes_rhs_tripos <- expression_matrix_tripos[genes_rhs, ]
expression_genes_gp63_tripos <- expression_matrix_tripos[genes_gp63, ]
expression_genes_dgf1_tripos <- expression_matrix_tripos[genes_dgf1, ]



# Calculate the number of cells in which each gene is detected for core and disruptive genes
celulas_expresion_core <- apply(expression_genes_core_tripos, 1, function(x) sum(x > 0))
celulas_expresion_disruptive <- apply(expression_genes_disruptive_tripos, 1, function(x) sum(x > 0))
celulas_expresion_ribosomales <- apply(expression_genes_ribosomales_tripos, 1, function(x) sum(x > 0))

# Create dataframes for the data
datos_expresion_core <- data.frame(
  Genes = rownames(expression_genes_core_tripos),
  Celulas_Expresion = celulas_expresion_core,
  Compartimento = rep("Core", length(celulas_expresion_core))
)

datos_expresion_disruptive <- data.frame(
  Genes = rownames(expression_genes_disruptive_tripos),
  Celulas_Expresion = celulas_expresion_disruptive,
  Compartimento = rep("Disruptive", length(celulas_expresion_disruptive))
)

datos_expresion_ribosomales <- data.frame(
  Genes = rownames(expression_genes_ribosomales_tripos),
  Celulas_Expresion = celulas_expresion_ribosomales,
  Compartimento = rep("Ribosomal", length(celulas_expresion_ribosomales))
)



# Perform subsampling of core genes to match the number of disruptive genes
set.seed(123)  # Set seed for reproducibility
num_genes_disruptive <- nrow(expression_genes_disruptive_tripos)
genes_core_subsampled <- expression_genes_core_tripos[sample(1:nrow(expression_genes_core_tripos), num_genes_disruptive), ]

# Calculate the number of cells in which each gene is detected in subsampled core genes
celulas_expresion_core_subsampled <- apply(genes_core_subsampled, 1, function(x) sum(x > 0))

# Create a dataframe for subsampled core genes
datos_expresion_core_subsampled <- data.frame(
  Genes = rownames(genes_core_subsampled),
  Celulas_Expresion = celulas_expresion_core_subsampled,
  Compartimento = rep("Core (Subsampleado)", length(celulas_expresion_core_subsampled))
)

# Combine dataframes of subsampled core genes and disruptive genes
datos_expresion_combinados <- rbind(datos_expresion_core_subsampled, datos_expresion_disruptive, datos_expresion_ribosomales)

# Boxplot to compare the distribution of number of cells with expression between core and disruptive genes
plot <- ggplot(datos_expresion_combinados, aes(x = Compartimento, y = Celulas_Expresion, fill = Compartimento)) +
  geom_boxplot(width=0.5, lwd = 1) +
  scale_fill_manual(values = c("Core (Subsampleado)" = "salmon", "Disruptive" = "skyblue", "Ribosomal" = "orange")) +  # Custom colors
  scale_y_log10() +
  labs(title = "Distribution of Number of Cells with Expression", x = "Compartment", y = "Number of Cells with Expression") +
  theme_minimal()

plot

#ggsave("Figuras/boxplot_numCellExpr_core_disruptive_ribosomales.png", plot = plot, dpi = 500, width = 10, height = 6)


# Violin plot with log scale on y axis and overlaid points with small size
ggplot(datos_expresion_combinados, aes(x = Compartimento, y = Celulas_Expresion, fill = Compartimento)) +
  geom_violin() +
  geom_boxplot(width=0.1, lwd = 1) +# Use geom_violin() instead of geom_boxplot()
  scale_fill_manual(values = c("Core (Subsampleado)" = "salmon", "Disruptive" = "skyblue", "Ribosomal" = "orange")) +  # Custom colors
  #geom_jitter(shape = 16, color = "black", size = 1, width = 0.2, alpha = 0.5) +  # Add overlaid points
  scale_y_log10() +  # Apply log scale on y axis
  labs(title = "Distribution of Number of Cells with Expression (Subsampled Core vs Disruptive Genes)", x = "Compartment", y = "Number of Cells with Expression (Log10)") +
  theme_minimal()

#ggsave("Figuras/violin_numCellExpr_core_disruptive_ribosomales.png", plot = plot, dpi = 500, width = 10, height = 6)











#### CORE vs SUBSAMPLED CORE ####

# Combine dataframes of all multigene families, subsampled core genes and disruptive genes
datos_expresion_combinados <- rbind(datos_expresion_core,datos_expresion_core_subsampled)

library(ggpubr)
library(dplyr)

# Define colors for core and disruptive genes
colores <- c('Core' = '#EB1E2C','Core (Subsampleado)' = '#EB1E2C')

# Ensure Compartimento is a factor with levels in the order of colors
datos_expresion_combinados$Compartimento <- factor(datos_expresion_combinados$Compartimento, levels = names(colores))

# Calculate p-values for the comparison of interest and filter significant ones
p_values <- compare_means(Celulas_Expresion ~ Compartimento, data = datos_expresion_combinados, method = "wilcox.test")

# Filter significant comparisons (p < 0.05)
p_values_significativos <- p_values %>%
  filter(p < 0.05)

# Define comparisons
my_comparisons <- list(
  c('Core (Subsampleado)', 'Core')
)

custom_labels <- c('Core' = 'Core','Core (Subsampleado)' = 'Sumbsampled Core')

# Create the plot




# Set up the PNG device
png("violinplot_numCellExpr_core_subsampledCore.png", width = 16, height = 10, units = 'cm', res = 300)

# Plot violinplot with boxplot, no legend, no grid and axis lines
ggplot(datos_expresion_combinados, aes(x = Compartimento, y = Celulas_Expresion, fill = Compartimento)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.9), alpha = 0.7) +  # Add violin plot with transparency
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA) +  # Add boxplot inside the violin plot
  scale_fill_manual(values = colores) +  # Set custom colors
  scale_y_log10() +
  labs(
    y = "Number of cells with expression (log10)",  # Y-axis title
    x = NULL  # Remove X-axis title
  ) +
  theme_bw() +  # Change to theme with white background
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # X-axis labels at 45 degrees
    axis.text.y = element_text(size = 10),  # Reduce Y-axis font size
    axis.title.x = element_text(size = 12), # Reduce X-axis title font size
    axis.title.y = element_text(size = 12), # Reduce Y-axis title font size
    axis.line = element_line(color = "black", size = 0.5),  # Add lines on axes and adjust thickness
    panel.grid = element_blank(),  # Remove grid
    legend.position = "none"  # Remove legend
  ) + 
  scale_x_discrete(labels = custom_labels)  # Customize X-axis labels

# Turn off the device
dev.off()









# Calculate the number of cells in which each gene is detected for each multigene family
celulas_expresion_mucinas <- apply(expression_genes_mucinas_tripos, 1, function(x) sum(x > 0))
celulas_expresion_masp <- apply(expression_genes_masp_tripos, 1, function(x) sum(x > 0))
celulas_expresion_trans <- apply(expression_genes_trans_tripos, 1, function(x) sum(x > 0))
celulas_expresion_rhs <- apply(expression_genes_rhs_tripos, 1, function(x) sum(x > 0))
celulas_expresion_gp63 <- apply(expression_genes_gp63_tripos, 1, function(x) sum(x > 0))
celulas_expresion_dgf1 <- apply(expression_genes_dgf1_tripos, 1, function(x) sum(x > 0))
celulas_expresion_ribosomales <- apply(expression_genes_ribosomales_tripos, 1, function(x) sum(x > 0))

# Create dataframes for each multigene family
datos_expresion_mucinas <- data.frame(
  Genes = rownames(expression_genes_mucinas_tripos),
  Celulas_Expresion = celulas_expresion_mucinas,
  Compartimento = rep("Mucinas", length(celulas_expresion_mucinas))
)

datos_expresion_masp <- data.frame(
  Genes = rownames(expression_genes_masp_tripos),
  Celulas_Expresion = celulas_expresion_masp,
  Compartimento = rep("MASP", length(celulas_expresion_masp))
)


datos_expresion_trans <- data.frame(
  Genes = rownames(expression_genes_trans_tripos),
  Celulas_Expresion = celulas_expresion_trans,
  Compartimento = rep("Trans-sialidasas", length(celulas_expresion_trans))
)

datos_expresion_rhs <- data.frame(
  Genes = rownames(expression_genes_rhs_tripos),
  Celulas_Expresion = celulas_expresion_rhs,
  Compartimento = rep("RHS", length(celulas_expresion_rhs))
)

datos_expresion_gp63 <- data.frame(
  Genes = rownames(expression_genes_gp63_tripos),
  Celulas_Expresion = celulas_expresion_gp63,
  Compartimento = rep("GP63", length(celulas_expresion_gp63))
)

datos_expresion_dgf1 <- data.frame(
  Genes = rownames(expression_genes_dgf1_tripos),
  Celulas_Expresion = celulas_expresion_dgf1,
  Compartimento = rep("DGF1", length(celulas_expresion_dgf1))
)

datos_expresion_ribosomales <- data.frame(
  Genes = rownames(expression_genes_ribosomales_tripos),
  Celulas_Expresion = celulas_expresion_ribosomales,
  Compartimento = rep("Ribosomal", length(celulas_expresion_ribosomales))
)

# Combine dataframes of all multigene families, subsampled core genes and disruptive genes
datos_expresion_combinados <- rbind(datos_expresion_core_subsampled,
                                    datos_expresion_disruptive,
                                    datos_expresion_ribosomales,
                                    datos_expresion_mucinas,
                                    datos_expresion_masp,
                                    datos_expresion_trans, 
                                    datos_expresion_rhs, 
                                    datos_expresion_gp63, 
                                    datos_expresion_dgf1)

library(ggpubr)
library(dplyr)

# Define colors for core and disruptive genes
colores <- c('Core (Subsampleado)' = '#EB1E2C', 'Disruptive' = '#91DCEA','Ribosomal' = '#FD6F30' ,"Mucinas" = "#F9A729", "MASP" = "#F9D23C", "Trans-sialidasas" = "#5FBB68", "RHS" = "#64CDCC", 
             "GP63" = "#A4A4D5", "DGF1" = "#BBC9E5")

# Ensure Compartimento is a factor with levels in the order of colors
datos_expresion_combinados$Compartimento <- factor(datos_expresion_combinados$Compartimento, levels = names(colores))

# Calculate p-values for the comparison of interest and filter significant ones
p_values <- compare_means(Celulas_Expresion ~ Compartimento, data = datos_expresion_combinados, method = "wilcox.test")

# Filter significant comparisons (p < 0.05)
p_values_significativos <- p_values %>%
  filter(p < 0.05)

# Define comparisons
my_comparisons <- list(
  c('Core (Subsampleado)', 'Disruptive'),
  c('Mucinas', 'MASP'),
  c('Mucinas', 'Trans-sialidasas'),
  c('Mucinas', 'RHS'),
  c('Mucinas', 'DGF1'),
  c('MASP', 'Trans-sialidasas'),
  c('MASP', 'RHS'),
  c('MASP', 'DGF1'),
  c('Trans-sialidasas', 'RHS'),
  c('RHS', 'GP63'),
  c('RHS', 'DGF1')
)

custom_labels <- c('Core (Subsampleado)' = 'Core', 'Disruptive' = 'Disruptive', "Ribosomal" = "Ribosomal","Mucinas" = "Mucins", "MASP" = "MASPs", "Trans-sialidasas" = "Trans-sialidases", "RHS" = "RHS", 
                   "GP63" = "GP63", "DGF1" = "DGF-1")

# Create the plot

write.csv(p_values_significativos,
          file = "p_values_significativos_numCells_all.csv",
          row.names = FALSE)



# Set up the PNG device
png("violinplot_numCellExpr_core_familias.png", width = 16, height = 10, units = 'cm', res = 300)

# Plot violinplot with boxplot, no legend, no grid and axis lines
ggplot(datos_expresion_combinados, aes(x = Compartimento, y = Celulas_Expresion, fill = Compartimento)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.9), alpha = 0.7) +  # Add violin plot with transparency
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA) +  # Add boxplot inside the violin plot
  scale_fill_manual(values = colores) +  # Set custom colors
  scale_y_log10() +
  labs(
    y = "Number of cells with expression (log10)",  # Y-axis title
    x = NULL  # Remove X-axis title
  ) +
  theme_bw() +  # Change to theme with white background
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # X-axis labels at 45 degrees
    axis.text.y = element_text(size = 10),  # Reduce Y-axis font size
    axis.title.x = element_text(size = 12), # Reduce X-axis title font size
    axis.title.y = element_text(size = 12), # Reduce Y-axis title font size
    axis.line = element_line(color = "black", size = 0.5),  # Add lines on axes and adjust thickness
    panel.grid = element_blank(),  # Remove the grid
    legend.position = "none"  # Remove the legend
  ) + 
  scale_x_discrete(labels = custom_labels)  # Customize X-axis labels

# Turn off the device
dev.off()




#### Lorenz curves for four groups: Core Subsampled, Disruptive, Ribosomal, Transialidases ####

library(dplyr)
library(ggplot2)

#### 1) Function to obtain Lorenz curve ####

lorenz_curve_df <- function(values, group_name){
  v <- sort(values)
  n <- length(v)
  
  L <- cumsum(v) / sum(v)         # cumulative expression
  p <- seq(1/n, 1, length.out=n)  # cumulative proportion of genes
  
  data.frame(
    p = p,
    L = L,
    Group = group_name
  )
}

#### 2) Sum expression per gene for each group ####

expr_trans       <- rowSums(expression_genes_trans_tripos)
expr_core_sub    <- rowSums(genes_core_subsampled)
expr_disruptive  <- rowSums(expression_genes_disruptive_tripos)
expr_ribosomal   <- rowSums(expression_genes_ribosomales_tripos)

#### 3) Lorenz curves for each group ####

lorenz_trans      <- lorenz_curve_df(expr_trans,      "Transialidases")
lorenz_core       <- lorenz_curve_df(expr_core_sub,   "Core Subsampled")
lorenz_disruptive <- lorenz_curve_df(expr_disruptive, "Disruptive")
lorenz_ribosomal  <- lorenz_curve_df(expr_ribosomal,  "Ribosomal")

# Line of perfect equality (dashed)
lorenz_equal <- data.frame(
  p = c(0, 1),
  L = c(0, 1),
  Group = "Perfect equality"
)

#### 4) Combine into a single dataframe ####

lorenz_all <- bind_rows(
  lorenz_trans,
  lorenz_core,
  lorenz_disruptive,
  lorenz_ribosomal,
  lorenz_equal
)

#### 5) Plot Lorenz curves ####
plot_lorenz_4 <- ggplot() +
  
  # Perfect equality (dashed line)
  geom_line(data = lorenz_equal,
            aes(x = p, y = L, color = "Perfect equality"),
            size = 1.1, linetype = "dashed") +
  
  # Real curves
  geom_line(data = lorenz_all %>% filter(Group != "Perfect equality"),
            aes(x = p, y = L, color = Group),
            size = 1.1) +
  
  scale_color_manual(values = c(
    "Perfect equality" = "black",
    "Transialidases"   = "#5FBB68",
    "Core Subsampled"  = "#EB1E2C",
    "Disruptive"       = "#91DCEA",
    "Ribosomal"        = "#FD6F30"
  )) +
  
  # Force axes 0–1
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  
  labs(
    title = "Lorenz Curves – Four Gene Groups",
    x = "Cumulative proportion of genes",
    y = "Cumulative proportion of expression",
    color = ""
  ) +
  
  # Clean theme without box
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),       # no grid
    panel.border = element_blank(),     # no border
    axis.line = element_line(color = "black"),  # black axes
    axis.ticks = element_line(color = "black"), # black ticks
    plot.title = element_text(face = "bold")
  )

#### 6) Display the plot ####

plot_lorenz_4

#### 7) Save the figure ####

ggsave("lorenz_4groups_coreSub_disruptive_trans_ribosomal.png",
       plot = plot_lorenz_4,
       dpi = 300, width = 8, height = 8)

#### END ####







#### EXPRESSION AND NUMBER OF CELLS CORE VS DISRUPTIVE ####

#### CORE ####


# Step 1: Calculate the sum of expression and the number of surface proteins expressed per cell
# Filter and find valid genes in your Seurat object
valid_core <- intersect(genes_core, rownames(filtered_Seurat[["RNA"]]))

# Calculate the sum of surface protein expression for each cell
sum_core <- Matrix::colSums(filtered_Seurat[valid_core, ]@assays$RNA@layers$data)

# Calculate the number of surface proteins expressed per cell
num_core <- Matrix::colSums(filtered_Seurat[sum_core, ]@assays$RNA@layers$counts > 0)


# Add these data to the Seurat object
filtered_Seurat[['sum_core']] <- sum_core
filtered_Seurat[['num_core']] <- num_core









#### DISRUPTIVE ####


# Step 1: Calculate the sum of expression and the number of surface proteins expressed per cell


# Filter and find valid genes in your Seurat object
valid_disruptive <- intersect(genes_disruptive, rownames(filtered_Seurat[["RNA"]]))

# Calculate the sum of surface protein expression for each cell
sum_disruptive <- Matrix::colSums(filtered_Seurat[valid_disruptive, ]@assays$RNA@layers$data)

# Calculate the number of surface proteins expressed per cell
num_disruptive <- Matrix::colSums(filtered_Seurat[sum_disruptive, ]@assays$RNA@layers$counts > 0)


# Add these data to the Seurat object
filtered_Seurat[['sum_disruptive']] <- sum_disruptive
filtered_Seurat[['num_disruptive']] <- num_disruptive



# Filter and find valid genes in your Seurat object
valid_subsampleCore <- intersect(rownames(genes_core_subsampled), rownames(filtered_Seurat[["RNA"]]))

# Calculate the sum of surface protein expression for each cell
sum_subsampleCore <- Matrix::colSums(filtered_Seurat[valid_subsampleCore, ]@assays$RNA@layers$data)

# Calculate the number of surface proteins expressed per cell
num_subsampleCore <- Matrix::colSums(filtered_Seurat[sum_subsampleCore, ]@assays$RNA@layers$counts > 0)


# Add these data to the Seurat object
filtered_Seurat[['sum_subsampleCore']] <- sum_subsampleCore
filtered_Seurat[['num_subsampleCore']] <- num_subsampleCore



# Load required packages
library(Seurat)
library(ggplot2)
library(cowplot)
library(viridis)

# Create combined dataframe with UMAP coordinates and necessary data
umap_data <- filtered_Seurat@reductions$umap@cell.embeddings

# Create the combined dataframe
all_data <- data.frame(
  umap_1 = umap_data[, 1],
  umap_2 = umap_data[, 2],
  sum_subsampleCore = sum_subsampleCore,
  num_subsampleCore = num_subsampleCore,
  sum_disruptive = sum_disruptive,
  num_disruptive = num_disruptive
)

# Find common limits for the color scale and sizes
min_color <- quantile(unlist(all_data[, c('sum_subsampleCore', 'sum_disruptive')]), probs = 0.1, na.rm = TRUE)
max_color <- quantile(unlist(all_data[, c('sum_subsampleCore', 'sum_disruptive')]), probs = 0.9, na.rm = TRUE)
min_size <- min(unlist(all_data[, c('num_subsampleCore', 'num_disruptive')]), na.rm = TRUE)
max_size <- max(unlist(all_data[, c('num_subsampleCore', 'num_disruptive')]), na.rm = TRUE)

# Define the CustomFeaturePlot function with shared scale limits
CustomFeaturePlot <- function(data, feature, size_feature, min_color, max_color, min_size, max_size) {
  # Calculate quantiles for clipping and adjust feature values
  data[[feature]] <- pmin(pmax(data[[feature]], min_color), max_color)
  
  # Create the UMAP plot with ggplot
  p <- ggplot(data, aes_string(x = "umap_1", y = "umap_2")) +
    geom_point(aes_string(color = feature, size = size_feature), alpha = 0.5) +
    scale_color_viridis(option = "C", limits = c(min_color, max_color)) +
    scale_size_continuous(limits = c(min_size, max_size)) +
    labs(color = feature, size = size_feature) +
    theme_minimal() +  # White background and grid removed
    theme(
      panel.grid = element_blank(),          # Remove the grid
      panel.background = element_blank(),     # Remove the grey background
      axis.line = element_line(color = "black"),  # Black axis lines
      axis.ticks = element_line(color = "black"), # Black ticks
      legend.position = "right"               # Legend position
    )
  
  return(p)
}

# Create the plots individually with shared scale limits
plot1 <- CustomFeaturePlot(all_data, 'sum_subsampleCore', 'num_subsampleCore', min_color, max_color, min_size, max_size)
plot2 <- CustomFeaturePlot(all_data, 'sum_disruptive', 'num_disruptive', min_color, max_color, min_size, max_size)

# Display the plots separately
print(plot1)
print(plot2)

#ggsave("Figuras/featurePlot_core_disruptive.png", plot = combined_plot, dpi = 300, width = 18, height = 6)




#### with z-scores core and disruptive ####

z_core <- scale(sum_core)
filtered_Seurat[['z_core']] <- z_core

z_subsampleCore <- scale(sum_subsampleCore)
filtered_Seurat[['z_subsampleCore']] <- z_subsampleCore

z_disruptive <- scale(sum_disruptive)
filtered_Seurat[['z_disruptive']] <- z_disruptive



# Load required packages
library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)

# Define the CustomFeaturePlot function with quantiles and text size adjustment
CustomFeaturePlot <- function(seurat_object, feature, size_feature) {
  data <- FetchData(seurat_object, vars = c("umap_1", "umap_2", feature, size_feature))
  
  # Calculate quantiles for clipping and adjust feature values
  min_cutoff <- quantile(data[[feature]], probs = 0.1, na.rm = TRUE)
  max_cutoff <- quantile(data[[feature]], probs = 0.9, na.rm = TRUE)
  data[[feature]] <- pmin(pmax(data[[feature]], min_cutoff), max_cutoff)
  
  # Create the UMAP plot with ggplot
  p <- ggplot(data, aes_string(x = "umap_1", y = "umap_2")) +
    geom_point(aes_string(color = feature, size = size_feature), alpha = 0.5) +
    scale_color_viridis(option = "C", discrete = FALSE) +
    labs(color = feature, size = size_feature) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),                    # Remove the grid
      panel.background = element_blank(),              # White background
      axis.line = element_line(color = "black"),       # Black axis lines
      axis.ticks = element_line(color = "black"),      # Black ticks
      axis.text = element_text(size = 5),              # Axis numbers size
      axis.title = element_text(size = 6),             # Axis labels size
      legend.position = "right",                       # Legend position
      legend.key.size = unit(0.2, "cm"),               # Legend key size
      legend.text = element_text(size = 5),            # Legend text size
      legend.title = element_text(size = 6)            # Legend title size
    ) +
    guides(
      size = guide_legend(override.aes = list(size = 3)),   # Point size in legend
      color = guide_colorbar(barwidth = 0.2, barheight = 3) # Colorbar size in legend
    )
  
  return(p)
}

# Create the plots individually
p1 <- CustomFeaturePlot(filtered_Seurat, 'z_subsampleCore', 'num_subsampleCore')
p2 <- CustomFeaturePlot(filtered_Seurat, 'z_disruptive', 'num_disruptive')

# Save both plots as PNG files with adjusted size
png("featurePlot_zscore_subsamplecore.png", width = 8, height = 6, units = 'cm', res = 500)
print(p1)
dev.off()

png("featurePlot_zscore_disruptive.png", width = 8, height = 6, units = 'cm', res = 500)
print(p2)
dev.off()




# Load required packages
library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)

# Define the CustomFeaturePlot function with quantiles and compact legend
CustomFeaturePlot <- function(seurat_object, feature, size_feature) {
  data <- FetchData(seurat_object, vars = c("umap_1", "umap_2", feature, size_feature))
  
  # Calculate quantiles for clipping and adjust feature values
  min_cutoff <- quantile(data[[feature]], probs = 0.1, na.rm = TRUE)
  max_cutoff <- quantile(data[[feature]], probs = 0.9, na.rm = TRUE)
  data[[feature]] <- pmin(pmax(data[[feature]], min_cutoff), max_cutoff)
  
  # Create the UMAP plot with ggplot
  p <- ggplot(data, aes_string(x = "umap_1", y = "umap_2")) +
    geom_point(aes_string(color = feature, size = size_feature), alpha = 0.5) +
    scale_color_viridis(option = "C", discrete = FALSE) +
    labs(color = feature, size = size_feature) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),                    # Remove the grid
      panel.background = element_blank(),              # Remove the grey background
      axis.line = element_line(color = "black"),       # Black axis lines
      axis.ticks = element_line(color = "black"),      # Black ticks
      legend.position = "right",                       # Legend position
      legend.key.size = unit(0.2, "cm"),               # Legend icon size
      legend.text = element_text(size = 6),            # Legend text size
      legend.title = element_text(size = 7)            # Legend title size
    )
  
  return(p)
}

# Create the plots individually
p1 <- CustomFeaturePlot(filtered_Seurat, 'z_core', 'num_core')
p2 <- CustomFeaturePlot(filtered_Seurat, 'z_disruptive', 'num_disruptive')

# Combine the plots using patchwork
combined_plot <- p1 | p2

# Save the combined plot to a PNG file with adjusted size
png("featurePlot_zscore_scusampleCore_disruptive.png", width = 24, height = 8, units = 'cm', res = 500)
print(combined_plot)
dev.off()


png("featurePlot_zscore_subsampleCore.png", width = 24, height = 16, units = 'cm', res = 500)
print(p1)
dev.off()

png("featurePlot_zscore_disruptive.png", width = 24, height = 16, units = 'cm', res = 500)
print(p2)
dev.off()












#### Figure 3 ####

#### CLUSTER TRYPOMASTIGOTES BY TRANS-SIALIDASE EXPRESSION ####


# Create a Seurat object from the loaded data
trans_Seurat <- CreateSeuratObject(counts = expression_genes_trans_tripos)

# Normalize the data
trans_Seurat = NormalizeData(trans_Seurat)


trans_Seurat@assays$RNA$data <- expression_genes_trans_tripos

# Identify variable features (genes) in the data
trans_Seurat = FindVariableFeatures(trans_Seurat)

# Scale the data
trans_Seurat = ScaleData(trans_Seurat)

# Perform a principal component analysis (PCA)
trans_Seurat = RunPCA(trans_Seurat)

# Perform the JackStraw test to determine significance of principal components
#filtered_Seurat <- JackStraw(filtered_Seurat, num.replicate = 100)
#filtered_Seurat <- ScoreJackStraw(filtered_Seurat, dims = 1:20)

# Visualize JackStraw test results
#JackStrawPlot(filtered_Seurat, dims = 1:15)


# To visualize and decide how many PCs to use, 4 or 5
ElbowPlot(trans_Seurat)  

# Find nearest neighbors based on PCA dimensions
trans_Seurat = FindNeighbors(trans_Seurat, dims = 1:10)

# Identify cell clusters
trans_Seurat = FindClusters(trans_Seurat, resolution = 0.4)

# Run UMAP for 2D visualization
trans_Seurat = RunUMAP(trans_Seurat, dims = 1:10, n.neighbors = 50, min.dist = 0.01)




# See the trans-sialidases that are most expressed 

library(SingleCellExperiment)

# Create a SingleCellExperiment object
sce_trans <- SingleCellExperiment(assays = list(counts = as.matrix(expression_genes_trans_tripos)))



#### per cluster (cluster 0 1186 cells and cluster 1 1015 cells) ####


library(Seurat)
library(ggplot2)

# Extract cluster identity for each cell
cluster_info <- Idents(trans_Seurat)

# Separate cells from clusters 0 and 1
cells_cluster_0 <- WhichCells(trans_Seurat, idents = 0)
cells_cluster_1 <- WhichCells(trans_Seurat, idents = 1)

# Extract expression data for the clusters
expression_cluster_0 <- expression_genes_trans_tripos[, colnames(expression_genes_trans_tripos) %in% cells_cluster_0]
expression_cluster_1 <- expression_genes_trans_tripos[, colnames(expression_genes_trans_tripos) %in% cells_cluster_1]

# Calculate average expression per cell
avg_expression_cluster_0 <- colMeans(expression_cluster_0)
avg_expression_cluster_1 <- colMeans(expression_cluster_1)

# Create dataframes for ggplot2
avg_expression_df_0 <- data.frame(avg_expression = avg_expression_cluster_0, cluster = "Cluster 0")
avg_expression_df_1 <- data.frame(avg_expression = avg_expression_cluster_1, cluster = "Cluster 1")

# Combine both dataframes
combined_df <- rbind(avg_expression_df_0, avg_expression_df_1)

# Plot density of average expression per cell
plot <- ggplot(combined_df, aes(x = avg_expression, fill = cluster)) +
  geom_density(alpha = 0.8) +
  labs(title = "Density plot of average gene expression per cell",
       x = "Average expression",
       y = "Density") +
  scale_fill_manual(values = c("#4573A0", "#D63440"))

plot
#ggsave("Figuras/densityPlot_clusters_tripos_solo_trans.png", plot = plot, dpi = 500, width = 10, height = 6)







#### PLOT TRANS-SIALIDASE MARKERS IN TRYPOMASTIGOTES ####


# Load required libraries
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library('circlize')

# Assuming your Seurat object is 'filtered_Seurat'
# Find all markers for trans_Seurat
markers_tripos <- FindAllMarkers(trans_Seurat, only.pos = TRUE, logfc.threshold = 0.25)

# Filter markers by adjusted p-value and keep cluster grouping
markers_tripos <- markers_tripos %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster)# %>%
#slice_max(order_by = avg_log2FC, n = 10)

# Filter expression matrix to retain only markers
filtered_expression_matrix_tripos_trans <- expression_genes_trans_tripos[rownames(expression_genes_trans_tripos) %in% markers_tripos$gene, ]

# Extract cluster information from the Seurat object
cell_clusters <- data.frame(
  cell = colnames(trans_Seurat),  # Cell names (columns in expression matrix)
  cluster = as.factor(trans_Seurat@meta.data$seurat_clusters)  # Clusters assigned to cells
)

# Order cells (columns) by their clusters
cell_clusters_ordered <- cell_clusters[order(cell_clusters$cluster), ]
ordered_expression_matrix <- filtered_expression_matrix_tripos_trans[, cell_clusters_ordered$cell]

# Define colors for clusters
cluster_colors <- c("#4573A0", "#D63440")
names(cluster_colors) <- levels(cell_clusters_ordered$cluster)

# Top annotation (clusters)
top_annotation <- HeatmapAnnotation(
  cluster = cell_clusters_ordered$cluster, 
  col = list(cluster = cluster_colors)
)

# Load the annotations dataframe
annotations_2 <- annotations
annotations_2$gene_id <- gsub("_", "-", annotations_2$gene_id)

markers_tripos <- markers_tripos %>%
  left_join(annotations_2, by = c("gene" = "gene_id"))

# Create labels for the heatmap
gene_labels <- paste(markers_tripos$gene, ":", markers_tripos$Description)

# Step 3: Create heatmap without clustering
col_fun <- colorRamp2(c(-2, 0, 2), c("magenta", "black", "yellow"))

# Normalize expression matrix (if needed)
scaled_matrix <- t(scale(t(ordered_expression_matrix)))

# Create a row order vector according to clusters
row_order <- c(
  markers_tripos$gene[markers_tripos$cluster == 0],
  markers_tripos$gene[markers_tripos$cluster == 1]
)

# Create a factor to split rows into groups
row_split <- factor(c(rep(0, 3),rep(1,53)))


scaled_matrix <- scaled_matrix[row_order, ]

rownames(scaled_matrix) <- gene_labels  
# Create the heatmap
heatmap_genes <- Heatmap(scaled_matrix,
                         name = "Expression",
                         col = col_fun,
                         use_raster = TRUE,
                         show_row_names = TRUE,
                         row_names_gp = gpar(fontsize = 8),
                         row_names_max_width = unit(12, "cm"),
                         show_column_names = FALSE,
                         top_annotation = top_annotation,
                         cluster_rows = FALSE,
                         cluster_columns = FALSE,
                         #row_split = row_split,  # Split rows into groups
                         column_split = cell_clusters_ordered$cluster,
                         row_gap = unit(0.5, "mm"),  # Space between groups
                         column_gap = unit(0.5, "mm"))

# Draw the heatmap
draw(heatmap_genes, show_annotation_legend = FALSE, show_heatmap_legend = TRUE)








# Ensure cell IDs in cell_clusters and colnames of expression_matrix_tripos are comparable.
# First, select cells for each cluster.

# For cluster 0
cells_cluster_tripos_0 <- cell_clusters[cell_clusters$cluster == 0, "cell"]

# Subset expression matrix for cluster 0 cells
expression_cluster_tripos_0 <- expression_matrix_tripos[, colnames(expression_matrix_tripos) %in% cells_cluster_tripos_0]

# For cluster 1
cells_cluster_tripos_1 <- cell_clusters[cell_clusters$cluster == 1, "cell"]

# Subset expression matrix for cluster 1 cells
expression_cluster_tripos_1 <- expression_matrix_tripos[, colnames(expression_matrix_tripos) %in% cells_cluster_tripos_1]

# Now you have the matrices `expression_cluster_0` and `expression_cluster_1` as subsets.







#### Figure 3A ####

#### Create a new Seurat object only for TRYPOMASTIGOTES separated by TRANS-SIALIDASES with all genes ####


#### Plot subtypes of cluster 0 and other clusters correctly ####

tripos_Seurat <- subset(filtered_Seurat, idents = 0)

ind = match(rownames(trans_Seurat@meta.data), rownames(tripos_Seurat@meta.data))
tripos_Seurat@meta.data[ind, "Cell.Type"] = trans_Seurat@meta.data[, "seurat_clusters"]


umapCoord <- as.matrix(Embeddings(object = trans_Seurat[["umap"]]))

tripos_Seurat[["umap_trans"]] <- CreateDimReducObject(embeddings = umapCoord, key = "umap_", assay = DefaultAssay(tripos_Seurat))



#### UMAP of clusters and subclusters of trypomastigotes by trans-sialidases ####

# Create a copy of the object to avoid modifying the original
combined_Seurat <- filtered_Seurat

# Label cells of cluster 0 with unique subcluster labels
combined_Seurat$Cell.Type <- as.character(combined_Seurat$seurat_clusters)  # Start with original clusters
cluster0_cells <- rownames(subset(filtered_Seurat, idents = 0)@meta.data)  # Cells of cluster 0
subcluster_labels <- paste0("0_", trans_Seurat$seurat_clusters)  # Unique labels for subclusters
combined_Seurat$Cell.Type[cluster0_cells] <- subcluster_labels  # Assign unique labels to cluster 0 cells

# Ensure 'Cell.Type' is a factor with desired ordering
combined_Seurat$Cell.Type <- factor(
  combined_Seurat$Cell.Type,
  levels = c(
    unique(subcluster_labels),  # Subclusters of cluster 0 (0_0, 0_1, etc.)
    setdiff(unique(as.character(combined_Seurat$seurat_clusters)), "0")  # Original clusters except 0
  )
)

# Plot UMAP with unique labels
DimPlot(combined_Seurat, reduction = 'umap', group.by = 'Cell.Type', pt.size = 2) +
  ggtitle("UMAP with unique subclusters for cluster 0 and original clusters") +
  theme_minimal()





#### Figure 3B and 3C ####

#### HOW ARE CORE AND TRANS-SIALIDASES EXPRESSED IN THESE GROUPS? ####

library(ggplot2)
library(reshape2)

# Filter expression matrices for genes in genes_core and genes_transialidasas
# Subset for genes in genes_core
expression_core_cluster_tripos_0 <- expression_cluster_tripos_0[rownames(expression_cluster_tripos_0) %in% genes_core, ]
expression_core_cluster_tripos_1 <- expression_cluster_tripos_1[rownames(expression_cluster_tripos_1) %in% genes_core, ]

# Subset for genes in genes_transialidasas
expression_transialidasas_cluster_tripos_0 <- expression_cluster_tripos_0[rownames(expression_cluster_tripos_0) %in% genes_transialidasas, ]
expression_transialidasas_cluster_tripos_1 <- expression_cluster_tripos_1[rownames(expression_cluster_tripos_1) %in% genes_transialidasas, ]

# Create dataframes with average expression per cell for each gene group
avg_expression_core_cluster_0 <- colMeans(expression_core_cluster_tripos_0)
avg_expression_core_cluster_1 <- colMeans(expression_core_cluster_tripos_1)
avg_expression_transialidasas_cluster_0 <- colMeans(expression_transialidasas_cluster_tripos_0)
avg_expression_transialidasas_cluster_1 <- colMeans(expression_transialidasas_cluster_tripos_1)

# Create dataframes for ggplot2
core_df_cluster_0 <- data.frame(Expression = avg_expression_core_cluster_0, Cluster = "Cluster 0", Gene_Group = "Core")
core_df_cluster_1 <- data.frame(Expression = avg_expression_core_cluster_1, Cluster = "Cluster 1", Gene_Group = "Core")

transialidasas_df_cluster_0 <- data.frame(Expression = avg_expression_transialidasas_cluster_0, Cluster = "Cluster 0", Gene_Group = "Transialidasas")
transialidasas_df_cluster_1 <- data.frame(Expression = avg_expression_transialidasas_cluster_1, Cluster = "Cluster 1", Gene_Group = "Transialidasas")

# Combine dataframes
combined_core_df <- rbind(core_df_cluster_0, core_df_cluster_1)
combined_transialidasas_df <- rbind(transialidasas_df_cluster_0, transialidasas_df_cluster_1)

library(ggplot2)
library(ggpubr)

# Violin plot for genes_core with Wilcoxon test
plot_core <- ggplot(combined_core_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # Add points
  labs(title = "Expression of Core Genes per Cluster", x = "Cluster", y = "Average Expression") +
  scale_fill_manual(values = c("#4573A0", "#D63440")) +  # Custom colors
  theme_minimal() +
  stat_compare_means(method = "t.test", label = "p.signif")  # Wilcoxon test for significance

# Draw the plot
plot_core



# Adjust violin plot for trans-sialidases without legend
plot_transialidasas <- ggplot(combined_transialidasas_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 0.5) +  # Smaller points
  labs(title = "Expression of Transialidasas Genes per Cluster", x = "", y = "Average Expression") +
  scale_fill_manual(values = c("#4573A0", "#D63440")) +  # Custom colors
  theme_minimal(base_size = 10) +  # Base font size
  theme(
    panel.grid = element_blank(),             # Remove grids
    axis.line = element_line(color = "black"), # Black axis lines
    plot.title = element_text(size = 10, face = "bold"),  # Title size
    axis.text = element_text(size = 8),       # Axis text size
    axis.title = element_text(size = 9),      # Axis title size
    legend.position = "none"                  # Remove legend
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif")  # Wilcoxon test for significance

print(plot_transialidasas)






# Filter expression matrices for ribosomal genes
proteinas_ribosomales <- annotations[grepl("ribosome|ribosomal", annotations$Description, ignore.case = TRUE), "gene_id"]
proteinas_ribosomales <- gsub("_", "-", proteinas_ribosomales)

# Subset ribosomal genes in clusters
expression_ribosomales_cluster_tripos_0 <- expression_cluster_tripos_0[rownames(expression_cluster_tripos_0) %in% proteinas_ribosomales, ]
expression_ribosomales_cluster_tripos_1 <- expression_cluster_tripos_1[rownames(expression_cluster_tripos_1) %in% proteinas_ribosomales, ]

# Create dataframes with average expression per cell for ribosomal genes
avg_expression_ribosomales_cluster_0 <- colMeans(expression_ribosomales_cluster_tripos_0)
avg_expression_ribosomales_cluster_1 <- colMeans(expression_ribosomales_cluster_tripos_1)

# Create dataframes for ggplot2
ribosomales_df_cluster_0 <- data.frame(Expression = avg_expression_ribosomales_cluster_0, Cluster = "Cluster 0", Gene_Group = "ribosomales")
ribosomales_df_cluster_1 <- data.frame(Expression = avg_expression_ribosomales_cluster_1, Cluster = "Cluster 1", Gene_Group = "ribosomales")

# Combine dataframes
combined_ribosomales_df <- rbind(ribosomales_df_cluster_0, ribosomales_df_cluster_1)

library(ggplot2)
library(ggpubr)


# Adjust violin plot for ribosomal genes without legend
plot_ribosomales <- ggplot(combined_ribosomales_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 0.5) +  # Smaller points
  labs(title = "Expression of Ribosomal Protein Coding Genes per Cluster", x = "", y = "Average Expression") +
  scale_fill_manual(values = c("#4573A0", "#D63440")) +  # Custom colors
  theme_minimal(base_size = 10) +  # Base font size
  theme(
    panel.grid = element_blank(),             # Remove grids
    axis.line = element_line(color = "black"), # Black axis lines
    plot.title = element_text(size = 10, face = "bold"),  # Title size
    axis.text = element_text(size = 8),       # Axis text size
    axis.title = element_text(size = 9),      # Axis title size
    legend.position = "none"                  # Remove legend
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif")  # Wilcoxon test for significance

print(plot_ribosomales)





#### Supplementary Figure 2A ####


######## genes de transportadores

# Filter expression matrices for genes of interest

# Subset of the genes from the object 
transporters <- read.csv('transporters_dm28c.csv', sep = ',')
transporters <- transporters$Gene.ID
transporters <- gsub('_', '-', transporters)

# Subset of the genes from the genes_transialidasas object
expression_transporters_cluster_tripos_0 <- expression_cluster_tripos_0[rownames(expression_cluster_tripos_0) %in% transporters, ]
expression_transporters_cluster_tripos_1 <- expression_cluster_tripos_1[rownames(expression_cluster_tripos_1) %in% transporters, ]

# Create dataframes with the average expression per cell for each gene group
avg_expression_transporters_cluster_0 <- colMeans(expression_transporters_cluster_tripos_0)
avg_expression_transporters_cluster_1 <- colMeans(expression_transporters_cluster_tripos_1)

# Create dataframes for ggplot2

transporters_df_cluster_0 <- data.frame(Expression = avg_expression_transporters_cluster_0, Cluster = "Cluster 0", Gene_Group = "transporters")
transporters_df_cluster_1 <- data.frame(Expression = avg_expression_transporters_cluster_1, Cluster = "Cluster 1", Gene_Group = "transporters")

# Combine the dataframes
combined_transporters_df <- rbind(transporters_df_cluster_0, transporters_df_cluster_1)

library(ggplot2)
library(ggpubr)


# Adjust the violin plot for transporters without legend
plot_transporters <- ggplot(combined_transporters_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 0.5) +  # Smaller points
  labs(title = "Expression of Phosphatases Protein Coding Genes per Cluster", x = "", y = "Average Expression") +
  scale_fill_manual(values = c("#4573A0", "#D63440")) +  # Custom colors
  theme_minimal(base_size = 10) +  # Base font size
  theme(
    panel.grid = element_blank(),             # Remove grids
    axis.line = element_line(color = "black"), # Black axis lines
    plot.title = element_text(size = 10, face = "bold"),  # Title size
    axis.text = element_text(size = 8),       # Axis text size
    axis.title = element_text(size = 9),      # Axis title size
    legend.position = "none"                  # Remove legend
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif")  # Wilcoxon test for significance

print(plot_transporters)




#### Supplementary Figure 2B ####


######## genes de polimerasas

# Filter expression matrices for genes of interest

# Subset of the genes from the object 
polimerasas <- read.csv('polymerase_dm28c.csv', sep = ',')
polimerasas <- polimerasas$Gene.ID
polimerasas <- gsub('_', '-', polimerasas)

# Subset of the genes from the genes_transialidasas object
expression_polimerasas_cluster_tripos_0 <- expression_cluster_tripos_0[rownames(expression_cluster_tripos_0) %in% polimerasas, ]
expression_polimerasas_cluster_tripos_1 <- expression_cluster_tripos_1[rownames(expression_cluster_tripos_1) %in% polimerasas, ]

# Create dataframes with the average expression per cell for each gene group
avg_expression_polimerasas_cluster_0 <- colMeans(expression_polimerasas_cluster_tripos_0)
avg_expression_polimerasas_cluster_1 <- colMeans(expression_polimerasas_cluster_tripos_1)

# Create dataframes for ggplot2

polimerasas_df_cluster_0 <- data.frame(Expression = avg_expression_polimerasas_cluster_0, Cluster = "Cluster 0", Gene_Group = "polimerasas")
polimerasas_df_cluster_1 <- data.frame(Expression = avg_expression_polimerasas_cluster_1, Cluster = "Cluster 1", Gene_Group = "polimerasas")

# Combine the dataframes
combined_polimerasas_df <- rbind(polimerasas_df_cluster_0, polimerasas_df_cluster_1)

library(ggplot2)
library(ggpubr)


# Adjust the violin plot for polymerases without legend
plot_polimerasas <- ggplot(combined_polimerasas_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 0.5) +  # Smaller points
  labs(title = "Expression of Polimerase Protein Coding Genes per Cluster", x = "", y = "Average Expression") +
  scale_fill_manual(values = c("#4573A0", "#D63440")) +  # Custom colors
  theme_minimal(base_size = 10) +  # Base font size
  theme(
    panel.grid = element_blank(),             # Remove grids
    axis.line = element_line(color = "black"), # Black axis lines
    plot.title = element_text(size = 10, face = "bold"),  # Title size
    axis.text = element_text(size = 8),       # Axis text size
    axis.title = element_text(size = 9),      # Axis title size
    legend.position = "none"                  # Remove legend
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif")  # Wilcoxon test for significance

print(plot_polimerasas)





#### Supplementary Figure 2C ####


######## genes de fosfatasas

# Filter expression matrices for genes of interest

# Subset of the genes from the object 
phosphatases <- read.csv('phosphatases_dm28c.csv', sep = ',')
phosphatases <- phosphatases$Gene.ID
phosphatases <- gsub('_', '-', phosphatases)

# Subset of the genes from the genes_transialidasas object
expression_phosphatases_cluster_tripos_0 <- expression_cluster_tripos_0[rownames(expression_cluster_tripos_0) %in% phosphatases, ]
expression_phosphatases_cluster_tripos_1 <- expression_cluster_tripos_1[rownames(expression_cluster_tripos_1) %in% phosphatases, ]

# Create dataframes with the average expression per cell for each gene group
avg_expression_phosphatases_cluster_0 <- colMeans(expression_phosphatases_cluster_tripos_0)
avg_expression_phosphatases_cluster_1 <- colMeans(expression_phosphatases_cluster_tripos_1)

# Create dataframes for ggplot2

phosphatases_df_cluster_0 <- data.frame(Expression = avg_expression_phosphatases_cluster_0, Cluster = "Cluster 0", Gene_Group = "phosphatases")
phosphatases_df_cluster_1 <- data.frame(Expression = avg_expression_phosphatases_cluster_1, Cluster = "Cluster 1", Gene_Group = "phosphatases")

# Combine the dataframes
combined_phosphatases_df <- rbind(phosphatases_df_cluster_0, phosphatases_df_cluster_1)

library(ggplot2)
library(ggpubr)


# Adjust the violin plot for phosphatases without legend
plot_phosphatases <- ggplot(combined_phosphatases_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 0.5) +  # Smaller points
  labs(title = "Expression of Phosphatases Protein Coding Genes per Cluster", x = "", y = "Average Expression") +
  scale_fill_manual(values = c("#4573A0", "#D63440")) +  # Custom colors
  theme_minimal(base_size = 10) +  # Base font size
  theme(
    panel.grid = element_blank(),             # Remove grids
    axis.line = element_line(color = "black"), # Black axis lines
    plot.title = element_text(size = 10, face = "bold"),  # Title size
    axis.text = element_text(size = 8),       # Axis text size
    axis.title = element_text(size = 9),      # Axis title size
    legend.position = "none"                  # Remove legend
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif")  # Wilcoxon test for significance

print(plot_phosphatases)




#### Supplementary Figure 2D ####


######## DISRUPTIVE genes


# Subset of the genes from the genes_transialidasas object
expression_disruptive_cluster_tripos_0 <- expression_cluster_tripos_0[rownames(expression_cluster_tripos_0) %in% genes_disruptive, ]
expression_disruptive_cluster_tripos_1 <- expression_cluster_tripos_1[rownames(expression_cluster_tripos_1) %in% genes_disruptive, ]

# Create dataframes with the average expression per cell for each gene group
avg_expression_disruptive_cluster_0 <- colMeans(expression_disruptive_cluster_tripos_0)
avg_expression_disruptive_cluster_1 <- colMeans(expression_disruptive_cluster_tripos_1)

# Create dataframes for ggplot2

disruptive_df_cluster_0 <- data.frame(Expression = avg_expression_disruptive_cluster_0, Cluster = "Cluster 0", Gene_Group = "disruptive")
disruptive_df_cluster_1 <- data.frame(Expression = avg_expression_disruptive_cluster_1, Cluster = "Cluster 1", Gene_Group = "disruptive")

# Combine the dataframes
combined_disruptive_df <- rbind(disruptive_df_cluster_0, disruptive_df_cluster_1)

library(ggplot2)
library(ggpubr)


# Adjust the violin plot for disruptive without legend
plot_disruptive <- ggplot(combined_disruptive_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 0.5) +  # Smaller points
  labs(title = "Expression of Disruptive Genes per Cluster", x = "", y = "Average Expression") +
  scale_fill_manual(values = c("#4573A0", "#D63440")) +  # Custom colors
  theme_minimal(base_size = 10) +  # Base font size
  theme(
    panel.grid = element_blank(),             # Remove grids
    axis.line = element_line(color = "black"), # Black axis lines
    plot.title = element_text(size = 10, face = "bold"),  # Title size
    axis.text = element_text(size = 8),       # Axis text size
    axis.title = element_text(size = 9),      # Axis title size
    legend.position = "none"                  # Remove legend
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif")  # Wilcoxon test for significance

print(plot_disruptive)





########## CALCULATION OF FOLD CHANGES BETWEEN CLUSTERS FOR EACH GENE GROUP


mean(avg_expression_transialidasas_cluster_0)/mean(avg_expression_transialidasas_cluster_1) # 1.527071

mean(avg_expression_ribosomales_cluster_0)/mean(avg_expression_ribosomales_cluster_1) # 1.033745

mean(avg_expression_phosphatases_cluster_0)/mean(avg_expression_phosphatases_cluster_1) # 1.095798

mean(avg_expression_transporters_cluster_0)/mean(avg_expression_transporters_cluster_1) # 1.093497

mean(avg_expression_polimerasas_cluster_0)/mean(avg_expression_polimerasas_cluster_1) # 1.145418

mean(avg_expression_core_cluster_0)/mean(avg_expression_core_cluster_1) # 1.123447

mean(avg_expression_disruptive_cluster_0)/mean(avg_expression_disruptive_cluster_1) # 1.380246







########### expression of other families in cluster 0 and 1


####### Mucins

# Subset of the genes from the genes_transialidasas object
expression_mucinas_cluster_tripos_0 <- expression_cluster_tripos_0[rownames(expression_cluster_tripos_0) %in% genes_mucinas, ]
expression_mucinas_cluster_tripos_1 <- expression_cluster_tripos_1[rownames(expression_cluster_tripos_1) %in% genes_mucinas, ]

# Create dataframes with the average expression per cell for each gene group
avg_expression_mucinas_cluster_0 <- colMeans(expression_mucinas_cluster_tripos_0)
avg_expression_mucinas_cluster_1 <- colMeans(expression_mucinas_cluster_tripos_1)

# Create dataframes for ggplot2

mucinas_df_cluster_0 <- data.frame(Expression = avg_expression_mucinas_cluster_0, Cluster = "Cluster 0", Gene_Group = "mucinas")
mucinas_df_cluster_1 <- data.frame(Expression = avg_expression_mucinas_cluster_1, Cluster = "Cluster 1", Gene_Group = "mucinas")

# Combine the dataframes
combined_mucinas_df <- rbind(mucinas_df_cluster_0, mucinas_df_cluster_1)

library(ggplot2)
library(ggpubr)


# Plot the violin plot for trans-sialidases genes with Wilcoxon test
plot_mucinas <- ggplot(combined_mucinas_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # Add points
  labs(title = "Expression of Mucin Protein Coding Genes per Cluster", x = "Cluster", y = "Average Expression") +
  scale_fill_manual(values = c("#4573A0", "#D63440")) +  # Custom colors
  theme_minimal() +
  stat_compare_means(method = "t.test", label = "p.signif")  # Wilcoxon test for significance

# Show the plots
plot_mucinas



####### MASP

# Subset of the genes from the genes_transialidasas object
expression_masp_cluster_tripos_0 <- expression_cluster_tripos_0[rownames(expression_cluster_tripos_0) %in% genes_masp, ]
expression_masp_cluster_tripos_1 <- expression_cluster_tripos_1[rownames(expression_cluster_tripos_1) %in% genes_masp, ]

# Create dataframes with the average expression per cell for each gene group
avg_expression_masp_cluster_0 <- colMeans(expression_masp_cluster_tripos_0)
avg_expression_masp_cluster_1 <- colMeans(expression_masp_cluster_tripos_1)

# Create dataframes for ggplot2

masp_df_cluster_0 <- data.frame(Expression = avg_expression_masp_cluster_0, Cluster = "Cluster 0", Gene_Group = "masp")
masp_df_cluster_1 <- data.frame(Expression = avg_expression_masp_cluster_1, Cluster = "Cluster 1", Gene_Group = "masp")

# Combine the dataframes
combined_masp_df <- rbind(masp_df_cluster_0, masp_df_cluster_1)

library(ggplot2)
library(ggpubr)


# Plot gp63 (note: uses combined_gp63_df)
plot_masp <- ggplot(combined_masp_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # Add points
  labs(title = "Expression of MASP Protein Coding Genes per Cluster", x = "Cluster", y = "Average Expression") +
  scale_fill_manual(values = c("#4573A0", "#D63440")) +  # Custom colors
  theme_minimal() +
  stat_compare_means(method = "t.test", label = "p.signif")  # Wilcoxon test for significance

# Show the plots
plot_masp






####### gp63

# Subset of the genes from the genes_transialidasas object
expression_gp63_cluster_tripos_0 <- expression_cluster_tripos_0[rownames(expression_cluster_tripos_0) %in% genes_gp63, ]
expression_gp63_cluster_tripos_1 <- expression_cluster_tripos_1[rownames(expression_cluster_tripos_1) %in% genes_gp63, ]

# Create dataframes with the average expression per cell for each gene group
avg_expression_gp63_cluster_0 <- colMeans(expression_gp63_cluster_tripos_0)
avg_expression_gp63_cluster_1 <- colMeans(expression_gp63_cluster_tripos_1)

# Create dataframes for ggplot2

gp63_df_cluster_0 <- data.frame(Expression = avg_expression_gp63_cluster_0, Cluster = "Cluster 0", Gene_Group = "gp63")
gp63_df_cluster_1 <- data.frame(Expression = avg_expression_gp63_cluster_1, Cluster = "Cluster 1", Gene_Group = "gp63")

# Combine the dataframes
combined_gp63_df <- rbind(gp63_df_cluster_0, gp63_df_cluster_1)

library(ggplot2)
library(ggpubr)


# Plot MASP (note: uses combined_masp_df)
plot_gp63 <- ggplot(combined_gp63_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # Add points
  labs(title = "Expression of gp63 Protein Coding Genes per Cluster", x = "Cluster", y = "Average Expression") +
  scale_fill_manual(values = c("#4573A0", "#D63440")) +  # Custom colors
  theme_minimal() +
  stat_compare_means(method = "t.test", label = "p.signif")  # Wilcoxon test for significance

# Show the plots
plot_gp63










#### Figure 4A ####

#### PERCENTAGE OF CELLS IN WHICH EACH TRANS-SIALIDASE IS EXPRESSED ####


# Calculate the percentage of cells in which each gene is expressed
percent_expressed <- apply(expression_genes_trans_tripos, 1, function(x) {
  mean(x > 0) * 100  # Calculate the percentage of values greater than zero
})

# Convert the result into a data frame and sort from highest to lowest
percent_expressed_trans_tripos_df <- data.frame(
  Gene = rownames(expression_genes_trans_tripos),
  Percent_Expressed = percent_expressed
) %>%
  arrange(desc(Percent_Expressed))

# View the top results
head(percent_expressed_trans_tripos_df)

# Save the data frame to a CSV file
write.csv(percent_expressed_trans_tripos_df, "percent_expression_trans_tripos.csv", row.names = FALSE)



# Make sure dplyr is loaded
library(dplyr)

# Calculate the percentage of cells in which each gene is expressed (cluster 0)
percent_expressed <- apply(expression_transialidasas_cluster_tripos_0, 1, function(x) {
  mean(x > 0) * 100  # Calculate the percentage of values greater than zero
})

# Convert the result into a data frame and sort from highest to lowest
percent_expressed_trans_tripos_df <- data.frame(
  Gene = rownames(expression_transialidasas_cluster_tripos_0),
  Percent_Expressed = percent_expressed
) %>%
  arrange(desc(Percent_Expressed))

# View the top results
head(percent_expressed_trans_tripos_df)

# Save the data frame to a CSV file
write.csv(percent_expressed_trans_tripos_df, "percent_expression_trans_tripos_cluster0.csv", row.names = FALSE)



# Calculate the percentage of cells in which each gene is expressed (cluster 1)
percent_expressed <- apply(expression_transialidasas_cluster_tripos_1, 1, function(x) {
  mean(x > 0) * 100  # Calculate the percentage of values greater than zero
})

# Convert the result into a data frame and sort from highest to lowest
percent_expressed_trans_tripos_df <- data.frame(
  Gene = rownames(expression_transialidasas_cluster_tripos_1),
  Percent_Expressed = percent_expressed
) %>%
  arrange(desc(Percent_Expressed))

# View the top results
head(percent_expressed_trans_tripos_df)

# Save the data frame to a CSV file
write.csv(percent_expressed_trans_tripos_df, "percent_expression_trans_tripos_cluster1.csv", row.names = FALSE)



#### Figure 4A: Genes contributing 75% of expression per cell (cluster 0) ####

library(dplyr)
library(tidyr)

# Step 1: Calculate total trans-sialidase expression per cell
total_expression_per_cell <- colSums(expression_transialidasas_cluster_tripos_0)

# Step 1: Calculate genes that contribute 75% of the total expression in each cell
genes_contributing_75_per_cell_cluster0 <- lapply(1:ncol(expression_transialidasas_cluster_tripos_0), function(i) {
  # Get gene expression for the current cell
  cell_expression <- expression_transialidasas_cluster_tripos_0[, i]
  
  # Sort genes in descending order by their expression
  sorted_gene_names <- rownames(expression_transialidasas_cluster_tripos_0)[order(cell_expression, decreasing = TRUE)]
  sorted_expression <- sort(cell_expression, decreasing = TRUE)
  
  # Calculate cumulative expression as a percentage of the cell total
  cumulative_expression <- cumsum(sorted_expression)
  cumulative_percentage <- cumulative_expression / sum(cell_expression) * 100
  
  # Select genes that contribute up to 75% of the total expression
  genes_75 <- sorted_gene_names[cumulative_percentage <= 75]
  return(genes_75)
})

# Step 2: Get the unique list of all genes that contribute 75% in at least one cell
unique_genes_75_cluster0 <- unique(unlist(genes_contributing_75_per_cell_cluster0))

# Print the result
print(unique_genes_75_cluster0)


# Create a matrix where only expression > 0 is considered for each cell
expression_positive <- expression_transialidasas_cluster_tripos_0
expression_positive[expression_positive == 0] <- NA  # Change 0s to NA to exclude them from sums

# Step 2: Filter only the genes that contribute to 75% of total expression
expression_positive <- expression_positive[rownames(expression_positive) %in% unique_genes_75_cluster0, ]

# Calculate total expression of trans-sialidases expressed (> 0) per cell
total_expression_positive <- colSums(expression_positive, na.rm = TRUE)

# Calculate percent expression for each trans-sialidase in each cell relative to total expressed trans-sialidases
percentage_expression <- sweep(expression_positive, 2, total_expression_positive, "/") * 100

# Replace NaN with 0 for genes not expressed in a given cell
percentage_expression[is.na(percentage_expression)] <- 0

# Load pheatmap
library(pheatmap)

# Create the heatmap
pheatmap(percentage_expression,
         cluster_rows = FALSE,            # Do not cluster genes
         cluster_cols = TRUE,             # Cluster cells
         scale = "none",                  # No scaling, already percentages
         color = colorRampPalette(c("white", "blue"))(10),
         main = "Percentage Expression of Transialidases per Cell",
         show_rownames = FALSE,
         show_colnames = FALSE,
         fontsize_row = 6,               # Adjust row text size
         fontsize_col = 6)               # Adjust column text size






#### Figure 4B ####

#### Z-SCORE OF AVERAGE EXPRESSION BY PERCENTAGE GROUP ####
# Load percent expression data
percent_expression <- read.csv('percent_expression_trans_tripos_cluster0.csv', sep = ',')

# Create 5%-width groups from 0% to 100% (excluding 0%)
percent_expression$Group <- cut(
  percent_expression$Percent_Expressed, 
  breaks = seq(0, 100, by = 5), 
  labels = paste0(seq(0, 95, by = 5), "-", seq(5, 100, by = 5), "%"), 
  include.lowest = TRUE
)

# Filter genes expressed in > 0% of cells
percent_expression <- subset(percent_expression, Percent_Expressed > 0)

library(dplyr)
library(ggplot2)

# Initialize list to store per-gene data
gene_data_list <- list()

# Iterate over each group and compute average expression of each gene
for (grp in unique(percent_expression$Group)) {
  # Select genes in the current group
  genes_in_group <- percent_expression$Gene[percent_expression$Group == grp]
  
  # Subset expression matrix for genes in the current group
  expression_subset <- expression_transialidasas_cluster_tripos_0[
    rownames(expression_transialidasas_cluster_tripos_0) %in% genes_in_group, 
  ]
  
  # Compute average expression for each gene (ignoring zeros)
  for (gene in rownames(expression_subset)) {
    gene_values <- expression_subset[gene, ]
    mean_expression <- ifelse(sum(gene_values > 0) > 0, 
                              mean(gene_values[gene_values > 0]), 
                              NA)
    
    if (!is.na(mean_expression)) {
      gene_data_list[[length(gene_data_list) + 1]] <- data.frame(
        Gene = gene,
        Group = grp,
        AvgExpression = mean_expression
      )
    }
  }
}

# Combine all gene data into a single DataFrame
all_gene_data <- do.call(rbind, gene_data_list)

# Compute Z-score of expression for each gene
all_gene_data$Zscore <- scale(all_gene_data$AvgExpression)[,1]

# Sort by Z-score and mark top 50 as red
all_gene_data <- all_gene_data %>%
  arrange(desc(Zscore)) %>%
  mutate(Rank = row_number(),
         Color = ifelse(Rank <= 50, "red", "gray"))

# Fix X-axis factor level order
all_gene_data$Group <- factor(
  all_gene_data$Group, 
  levels = paste0(seq(0, 95, by = 5), "-", seq(5, 100, by = 5), "%")
)


ggplot(all_gene_data, aes(x = Group, y = Zscore)) +
  geom_point(aes(color = Color), alpha = 0.7, size = 2) +
  scale_color_identity() +
  labs(
    title = "",
    x = "Percentage of cells",
    y = "Z-score"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13, color = "black", angle = 90),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.position = "none"
  ) +
  scale_y_continuous(breaks = seq(-2.5, 5, by = 2.5), limits = c(-4, 5.5))



#### Figure 4C ####

#### GINI INDEXES FOR EACH TRANS-SIALIDASE IN CLUSTER 0 ####

# Load DescTools if needed (commented install)
# if (!requireNamespace("DescTools", quietly = TRUE)) {
#   install.packages("DescTools")
# }

library(DescTools)
library(ggplot2)

# Compute the Gini index for each cell (columns) considering only values > 0
gini_indexes_cells <- apply(expression_transialidasas_cluster_tripos_0, 2, function(x) {
  positive_values <- x[x > 0]
  if (length(positive_values) > 0) {
    Gini(positive_values)
  } else {
    NA  # Assign NA if there are no positive values
  }
})

# Remove possible NAs (cells without positively expressed genes)
gini_indexes_cells <- na.omit(gini_indexes_cells)

# Create a data frame for visualization
gini_cells_df <- data.frame(Gini_Index = gini_indexes_cells)

# Create the violin plot with similar adjustments to the example
plot_gini_cells <- ggplot(gini_cells_df, aes(x = "", y = Gini_Index, fill = "Gini")) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +  # Violin plot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 0.5) +  # Jitter points
  labs(title = "Distribution of Gini Index Across Cells from Cluster 0", x = "", y = "Gini Index") +
  scale_fill_manual(values = c("#4573A0")) +  # Custom color
  theme_minimal(base_size = 10) +  # Base font size
  theme(
    panel.grid = element_blank(),             # Remove grids
    axis.line = element_line(color = "black"), # Black axis lines
    plot.title = element_text(size = 10, face = "bold"),  # Title size
    axis.text = element_text(size = 8),       # Axis text size
    axis.title = element_text(size = 9),      # Axis title size
    legend.position = "none"                  # Remove legend
  )

print(plot_gini_cells)



#### Figure 4D ####

library(ineq)   # provides Gini()
# install.packages("ineq")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ineq)  # for Gini calculation

# --- 1. Compute Lorenz-like curves and Gini per cell ---
lorenz_list <- apply(expression_transialidasas_cluster_tripos_0, 2, function(cell_expression) {
  # consider only expressed genes
  expr_filtrada <- cell_expression[cell_expression > 0]
  if (length(expr_filtrada) <= 1) return(NULL)
  
  # compute Gini
  gini_value <- Gini(expr_filtrada)
  
  # cumulative curves
  sorted_expr <- sort(expr_filtrada, decreasing = TRUE)
  cum_expr <- cumsum(sorted_expr) / sum(sorted_expr) * 100
  cum_genes <- seq_along(sorted_expr) / length(sorted_expr) * 100
  
  data.frame(
    cum_genes = cum_genes,
    cum_expr = cum_expr,
    gini = gini_value
  )
})

# remove NULLs (cells without expressed genes)
lorenz_list <- lorenz_list[!sapply(lorenz_list, is.null)]

# combine into a single data.frame
lorenz_df <- bind_rows(lorenz_list, .id = "cell")

# --- 2. Reference lines ---
perfect_equality <- data.frame(cum_genes = c(0, 100), cum_expr = c(0, 100))
perfect_inequality <- data.frame(cum_genes = c(100, 0, 0), cum_expr = c(100, 100, 0))

# --- 3. Identify curves with extreme Gini values ---
gini_summary <- lorenz_df %>%
  distinct(cell, gini) %>%
  arrange(gini)

lowest_gini_cell <- gini_summary$cell[1]
highest_gini_cell <- gini_summary$cell[nrow(gini_summary)]


ggplot() +
  # reference lines
  geom_line(data = perfect_inequality, aes(x = cum_genes, y = cum_expr),
            color = "black", size = 0.7) +
  geom_line(data = perfect_equality, aes(x = cum_genes, y = cum_expr),
            color = "black", linetype = "dashed", size = 0.7) +
  # curves for all cells
  geom_line(data = lorenz_df, aes(x = cum_genes, y = cum_expr, group = cell),
            color = "grey70", alpha = 0.3, size = 0.5) +
  # cell with lowest inequality (most equitable)
  geom_line(data = subset(lorenz_df, cell == lowest_gini_cell),
            aes(x = cum_genes, y = cum_expr),
            color = "#1b9e77", size = 1.1) +
  # cell with highest inequality
  geom_line(data = subset(lorenz_df, cell == highest_gini_cell),
            aes(x = cum_genes, y = cum_expr),
            color = "#d95f02", size = 1.1) +
  annotate("text", x = 65, y = 20, label = "Perfect inequality", size = 4) +
  annotate("text", x = 60, y = 55, label = "Perfect equality", size = 4, angle = 35) +
  scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  labs(
    x = "Cumulative % of Genes (expressed)",
    y = "Cumulative % of Expression",
    title = "Lorenz Curves of Transialidase Expression per Cell (Cluster 0)",
    subtitle = paste0(
      "Highlighted cells: lowest Gini (", lowest_gini_cell, ") and highest Gini (", highest_gini_cell, ")"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )


# --- 5. Save summary table of Gini per cell ---
write.csv(gini_summary, "gini_cluster0_tripos.csv", row.names = FALSE)








#############################################################################################
#### Supplementary Figure 3A ####
#############################################################################################


# Load library
library(writexl)

# Suppose 'expression_matrix_cluster' is the expression matrix for a single cluster
# (genes x cells)

# Step 1: Calculate the sum or the average expression per gene
# We'll use the sum of expressions for each gene across all cells of the cluster
pseudobulk_vector <- rowSums(expression_transialidasas_cluster_tripos_0)  # For sum
# or if you prefer the mean
# pseudobulk_vector <- colMeans(expression_matrix_cluster)  # For mean

# Step 2: Convert the vector to a suitable format for saving (if needed)
# Transpose it and convert to a data.frame for export
pseudobulk_df <- data.frame(Gene = rownames(expression_transialidasas_cluster_tripos_0), Expression = pseudobulk_vector)



bulk <- read.csv('abundance_tripos_1.txt', sep = '\t')

bulk$target_id <- gsub("_", "-", bulk$target_id)

bulk_df <- bulk[bulk$target_id %in% rownames(pseudobulk_df),]





# Sort bulk_df by the TPM column descending and keep the top 100
bulk_top100 <- bulk_df[order(-bulk_df$tpm), ][1:100, ]

# Sort pseudobulk_df by the Expression column descending and keep the top 100
pseudobulk_top100 <- pseudobulk_df[order(-pseudobulk_df$Expression), ][1:100, ]

# Extract unique identifiers from each dataframe
bulk_genes <- bulk_top100$target_id
pseudobulk_genes <- pseudobulk_top100$Gene  # Assuming this column contains the identifiers

# Install and load the library for Venn diagrams if not installed
library(VennDiagram)

# Create the Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(bulk_genes),
  area2 = length(pseudobulk_genes),
  cross.area = length(intersect(bulk_genes, pseudobulk_genes)),
  category = c("Bulk Top 100", "PseudoBulk Top 100"),
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 1,
  cat.cex = 0.5,
  cat.pos = c(-20, 20)
)

grid.draw(venn.plot)  # Draw the diagram to the graphic device





# Cross with TcS that are expressed in more than 40% of cells of cluster 0
top40perc_cluster0 <- read.csv('percent_expression_trans_tripos_cluster0.csv', sep = ',')
top40perc_cluster0 <- top40perc_cluster0[1:31,]
top40perc_cluster0 <- top40perc_cluster0$Gene



# Create the Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(bulk_genes),
  area2 = length(top40perc_cluster0),
  cross.area = length(intersect(bulk_genes, top40perc_cluster0)),
  category = c("Bulk Top 100 TcS", "TcS > 40% of cells C0"),
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 1,         # Size of numbers in the areas
  cat.cex = 0.5,   # Size of the group labels
  cat.pos = c(-30, 30),  # Adjust label positions (angle in degrees)
  cat.dist = c(0.05, 0.05) # Distance of labels from circles
)



grid.draw(venn.plot)  # Draw the diagram to the graphic device



# Sort bulk_df by the TPM column descending and keep the top 50
bulk_top50 <- bulk_df[order(-bulk_df$tpm), ][1:50, ]
# Extract unique identifiers from each dataframe
bulk_genes <- bulk_top50$target_id

# Create the Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(bulk_genes),
  area2 = length(top40perc_cluster0),
  cross.area = length(intersect(bulk_genes, top40perc_cluster0)),
  category = c("Bulk Top 50 TcS", "TcS > 40% of cells C0"),
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 1,         # Size of numbers in the areas
  cat.cex = 0.5,   # Size of the group labels
  cat.pos = c(-30, 30),  # Adjust label positions (angle in degrees)
  cat.dist = c(0.05, 0.05) # Distance of labels from circles
)








#############################################################################################
#### Supplementary Figure 3B ####
#############################################################################################


# Step 1: Create a copy of the original matrix to consider only expression values > 0
expression_positive <- expression_transialidasas_cluster_tripos_0
expression_positive[expression_positive == 0] <- NA  # Change 0s to NA to exclude them

# Step 2: Calculate the percentage of cells with expression > 0 for each gene
cell_percentage_expression <- rowMeans(!is.na(expression_positive)) * 100

# Step 3: Filter genes that have expression in at least 40% of the cells
genes_80_percent <- rownames(expression_positive)[cell_percentage_expression >= 40]
expression_positive <- expression_positive[genes_80_percent, ]

# Step 4: Further filter to only the genes selected in the 75% total-expression set
expression_positive <- expression_positive[rownames(expression_positive) %in% unique_genes_75_cluster0, ]

# Calculate total expression of expressed trans-sialidases (> 0) per cell
total_expression_positive <- colSums(expression_positive, na.rm = TRUE)

# Calculate percentage expression for each trans-sialidase in each cell relative to total expressed trans-sialidases
percentage_expression <- sweep(expression_positive, 2, total_expression_positive, "/") * 100

# Replace NaN values with 0 for genes not expressed in each cell
percentage_expression[is.na(percentage_expression)] <- 0

# Load pheatmap
library(pheatmap)

# Create the heatmap
pheatmap(percentage_expression,
         cluster_rows = TRUE,            # Cluster genes to see patterns
         cluster_cols = TRUE,            # Cluster cells to see patterns
         scale = "none",                 # No scaling, it's already percentage
         color = colorRampPalette(c("white", "blue"))(10),
         main = "Percentage Expression of Transialidases per Cell",
         show_rownames = FALSE,
         show_colnames = FALSE,
         fontsize_row = 6,               # Adjust row text size
         fontsize_col = 6)               # Adjust column text size





# Step 1: Create a copy of the original matrix to consider only expression values > 0
expression_positive <- expression_transialidasas_cluster_tripos_1
expression_positive[expression_positive == 0] <- NA  # Change 0s to NA to exclude them

# Step 2: Calculate the percentage of cells with expression > 0 for each gene
cell_percentage_expression <- rowMeans(!is.na(expression_positive)) * 100

# Step 3: Filter genes that have expression in at least 40% of the cells
genes_80_percent <- rownames(expression_positive)[cell_percentage_expression >= 40]
expression_positive <- expression_positive[genes_80_percent, ]

# Step 4: Further filter to only the genes selected in the 75% total-expression set
expression_positive <- expression_positive[rownames(expression_positive) %in% unique_genes_75_cluster1, ]

# Calculate total expression of expressed trans-sialidases (> 0) per cell
total_expression_positive <- colSums(expression_positive, na.rm = TRUE)

# Calculate percentage expression for each trans-sialidase in each cell relative to total expressed trans-sialidases
percentage_expression <- sweep(expression_positive, 2, total_expression_positive, "/") * 100

# Replace NaN values with 0 for genes not expressed in each cell
percentage_expression[is.na(percentage_expression)] <- 0

# Load pheatmap
library(pheatmap)

# Create the heatmap
pheatmap(percentage_expression,
         cluster_rows = TRUE,            # Cluster genes to see patterns
         cluster_cols = TRUE,            # Cluster cells to see patterns
         scale = "none",                 # No scaling, it's already percentage
         color = colorRampPalette(c("white", "blue"))(10),
         main = "Percentage Expression of Transialidases per Cell",
         show_rownames = FALSE,
         show_colnames = FALSE,
         fontsize_row = 6,               # Adjust row text size
         fontsize_col = 6)               # Adjust column text size








#############################################################################################
###################################### Figure 5A #############################################
#############################################################################################


######################### GENOMIC LOCATION OF THE MOST AND LEAST ABUNDANT TcS #################
############# IN THIS PART OF THE SCRIPT I DEFINE POLYCISTRON AS CHANGE OF STRAND ONLY #######
#############################################################################################


library('GOSemSim')
library('GenomicFeatures')

# load the gff to make it a txdb object
gff_path <- 'anotacion_dm28c_con_mt.gff'

# load the gff
library(ape)
tc_gff <- read.gff(gff_path) # I had to manually remove the initial ## lines from the gff
tc_gff <- tc_gff[tc_gff$type == 'transcript', ]


library(data.table)
setDT(tc_gff)  # convert if it's a data.frame

# Extract gene_id (if it doesn't exist, try to extract ID)
tc_gff[, gene_id := fifelse(grepl('gene_id\\s*"', attributes),
                            sub('.*gene_id\\s*"([^"]+)".*', '\\1', attributes),
                            fifelse(grepl('\\bID\\s*"', attributes),
                                    sub('.*\\bID\\s*"([^"]+)".*', '\\1', attributes),
                                    NA_character_))]

# Extract description (if it doesn't exist, NA)
tc_gff[, description := fifelse(grepl('description\\s*"', attributes),
                                sub('.*description\\s*"([^"]+)".*', '\\1', attributes),
                                NA_character_)]

# Replace the attributes column so it contains only the gene_id
# (if you prefer to rename it to 'gene_id' instead of overwriting, change the assignment)
tc_gff[, attributes := gene_id]





# libraries
#if (!requireNamespace("rtracklayer", quietly=TRUE)) BiocManager::install("rtracklayer")
#if (!requireNamespace("GenomicRanges", quietly=TRUE)) BiocManager::install("GenomicRanges")
#if (!requireNamespace("data.table", quietly=TRUE)) install.packages("data.table")
library(rtracklayer); library(GenomicRanges); library(data.table)

# BED files
bed_core <- "Dm28c_core_compartment_positions.bed"
bed_dis  <- "Dm28c_disruptive_compartment_positions.bed"

# ----------------------------------------------------------------
# 1) Normalize tc_gff -> GRanges (assume tc_gff already exists in env)
# ----------------------------------------------------------------
if (!exists("tc_gff")) stop("No encuentro 'tc_gff' en el environment. Cargalo primero.")

# If it's data.frame / data.table convert and group by gene_id
if (inherits(tc_gff, "data.frame") && !inherits(tc_gff, "GRanges")) {
  dt <- as.data.table(tc_gff)
  
  # ensure required columns exist
  stopifnot(all(c("seqid","start","end","strand","gene_id","description") %in% colnames(dt)))
  
  # group by gene_id: min(start), max(end), use first seqid/strand/description observed
  agg <- dt[, .(
    seqid = seqid[1L],
    start = min(start),
    end   = max(end),
    strand = strand[1L],
    description = if (all(is.na(description))) NA_character_ else description[which(!is.na(description))[1L]]
  ), by = gene_id]
  
  genes_gr <- GRanges(
    seqnames = agg$seqid,
    ranges   = IRanges(start = agg$start, end = agg$end),
    strand   = agg$strand
  )
  mcols(genes_gr)$gene_id     <- agg$gene_id
  mcols(genes_gr)$description <- agg$description
} else if (inherits(tc_gff, "GRanges")) {
  # if already GRanges, group by gene_id in case there are multiple rows per gene
  gr <- tc_gff
  stopifnot("gene_id" %in% colnames(mcols(gr)))
  dt <- as.data.table(data.frame(seqid = as.character(seqnames(gr)),
                                 start = start(gr), end = end(gr),
                                 strand = as.character(strand(gr)),
                                 gene_id = mcols(gr)$gene_id,
                                 description = mcols(gr)$description,
                                 stringsAsFactors = FALSE))
  agg <- dt[, .(
    seqid = seqid[1L],
    start = min(start),
    end   = max(end),
    strand = strand[1L],
    description = if (all(is.na(description))) NA_character_ else description[which(!is.na(description))[1L]]
  ), by = gene_id]
  
  genes_gr <- GRanges(
    seqnames = agg$seqid,
    ranges   = IRanges(start = agg$start, end = agg$end),
    strand   = agg$strand
  )
  mcols(genes_gr)$gene_id     <- agg$gene_id
  mcols(genes_gr)$description <- agg$description
} else {
  stop("tc_gff debe ser data.frame o GRanges.")
}

# ----------------------------------------------------------------
# 2) import BEDs (rtracklayer converts BED coords -> GRanges 1-based)
# ----------------------------------------------------------------
core_gr <- import(bed_core)
dis_gr  <- import(bed_dis)

# simple contig name adjustments if needed (remove 'chr' from BED if it doesn't match)
if (length(intersect(seqlevels(genes_gr), seqlevels(core_gr))) == 0) {
  seqlevels(core_gr) <- sub("^chr", "", seqlevels(core_gr))
  seqlevels(dis_gr)  <- sub("^chr", "", seqlevels(dis_gr))
}
try(seqlevelsStyle(core_gr) <- seqlevelsStyle(genes_gr), silent=TRUE)
try(seqlevelsStyle(dis_gr)  <- seqlevelsStyle(genes_gr), silent=TRUE)

# ----------------------------------------------------------------
# 3) compute overlaps (bp) and fractions per gene
# ----------------------------------------------------------------
hits_core <- findOverlaps(genes_gr, core_gr, ignore.strand = TRUE)
hits_dis  <- findOverlaps(genes_gr, dis_gr,  ignore.strand = TRUE)

sum_overlap_bp <- function(hits, A, B) {
  out <- integer(length(A))
  if (length(hits) == 0) return(out)
  q <- queryHits(hits); s <- subjectHits(hits)
  pis <- pintersect(A[q], B[s])
  w   <- width(pis)
  sums <- tapply(w, q, sum)
  out[as.integer(names(sums))] <- as.integer(sums)
  out
}

core_bp <- sum_overlap_bp(hits_core, genes_gr, core_gr)
dis_bp  <- sum_overlap_bp(hits_dis,  genes_gr, dis_gr)
gene_bp <- width(genes_gr)
core_frac <- core_bp / gene_bp
dis_frac  <- dis_bp  / gene_bp
core_frac[is.na(core_frac)] <- 0
dis_frac[is.na(dis_frac)]   <- 0

# assignment rule: larger fraction -> region (tie -> "both")
region <- ifelse(core_frac==0 & dis_frac==0, "none",
                 ifelse(abs(core_frac - dis_frac) < 1e-6, "both",
                        ifelse(core_frac > dis_frac, "core", "disruptive")))

mcols(genes_gr)$core_bp  <- core_bp
mcols(genes_gr)$dis_bp   <- dis_bp
mcols(genes_gr)$core_frac <- round(core_frac, 4)
mcols(genes_gr)$dis_frac  <- round(dis_frac, 4)
mcols(genes_gr)$region    <- region

# ----------------------------------------------------------------
# 4) output and save
# ----------------------------------------------------------------
res <- data.frame(
  gene_id     = mcols(genes_gr)$gene_id,
  contig      = as.character(seqnames(genes_gr)),
  start       = start(genes_gr),
  end         = end(genes_gr),
  width       = width(genes_gr),
  description = mcols(genes_gr)$description,
  core_bp     = mcols(genes_gr)$core_bp,
  dis_bp      = mcols(genes_gr)$dis_bp,
  core_frac   = mcols(genes_gr)$core_frac,
  dis_frac    = mcols(genes_gr)$dis_frac,
  region      = mcols(genes_gr)$region,
  stringsAsFactors = FALSE
)

write.csv(res, "genes_region_annotation_tc_gff_simple.csv", row.names = FALSE)









# Dependencies
#if (!requireNamespace("data.table", quietly=TRUE)) install.packages("data.table")
#if (!requireNamespace("rtracklayer", quietly=TRUE)) { 
#  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#  BiocManager::install("rtracklayer")
#}
library(data.table); library(rtracklayer)

# --- Input files (modify if needed) -------------------
# Assume tc_gff already exists in the environment and is data.frame or GRanges
# Lists:
tcs_high_file <- "tcs_high.txt"
tcs_low_file  <- "tcs_low.txt"

# --- 0) read gene lists ----------------------------------------
tcs_high <- fread(tcs_high_file, header = FALSE)[[1]]
tcs_low  <- fread(tcs_low_file, header = FALSE)[[1]]

# --- 1) prepare gene table from tc_gff -------------------------
# Accepts tc_gff as data.frame or GRanges; extract necessary columns
if (inherits(tc_gff, "GRanges")) {
  genes_dt <- data.table(
    seqid = as.character(seqnames(tc_gff)),
    start = start(tc_gff),
    end   = end(tc_gff),
    strand = as.character(strand(tc_gff)),
    gene_id = mcols(tc_gff)$gene_id,
    description = mcols(tc_gff)$description
  )
} else if (is.data.frame(tc_gff) || data.table::is.data.table(tc_gff)) {
  genes_dt <- as.data.table(tc_gff)[, .(seqid = seqid, start = as.integer(start), end = as.integer(end),
                                        strand = as.character(strand), gene_id = gene_id, description = description)]
} else stop("tc_gff debe ser data.frame o GRanges en el entorno.")

# sanity check
stopifnot(all(c("seqid","start","end","strand","gene_id") %in% colnames(genes_dt)))

# order and create polycistron as run-length id by strand within each seqid
setorder(genes_dt, seqid, start)
genes_dt[, poly_rleid := rleid(strand), by = seqid]
genes_dt[, poly_id := paste0(seqid, "_p", poly_rleid)]
# index within the polycistron (1..n)
genes_dt[, idx_in_poly := seq_len(.N), by = poly_id]
# polycistron size
genes_dt[, poly_size := .N, by = poly_id]

# --- 2) function that returns neighbors for a list of genes ----------
get_neighbors <- function(gene_list, genes_dt) {
  focals <- genes_dt[gene_id %in% gene_list, .(poly_id, gene_id_focal = gene_id, idx_focal = idx_in_poly,
                                               start_focal = start, end_focal = end, seqid_focal = seqid)]
  if (nrow(focals) == 0) {
    warning("Ninguno de los gene_id de la lista se encontró en tc_gff.")
    return(data.table())
  }
  # do merge (cartesian product by poly_id): each focal with all members of the same poly
  neigh <- merge(focals, genes_dt, by = "poly_id", allow.cartesian = TRUE)
  # rename neighbor columns
  setnames(neigh,
           old = c("gene_id","start","end","strand","description","idx_in_poly","poly_size","seqid"),
           new = c("gene_id_neighbor","start_neighbor","end_neighbor","strand_neighbor","description_neighbor","idx_neighbor","poly_size","seqid_neighbor"))
  # remove case where neighbor is the same as focal (if you don't want to exclude it, comment out)
  neigh <- neigh[gene_id_focal != gene_id_neighbor]
  # relative position in genes (neighbor index - focal index)
  neigh[, rel_pos := idx_neighbor - idx_focal]
  # distance in bp between the two genes (0 if overlapping)
  neigh[, dist_bp := fifelse(end_focal < start_neighbor, start_neighbor - end_focal,
                             fifelse(end_neighbor < start_focal, start_focal - end_neighbor, 0L))]
  # order results by focal gene and relative position
  setorder(neigh, gene_id_focal, rel_pos)
  # useful columns
  neigh[, .(gene_id_focal, seqid_focal, start_focal, end_focal,
            gene_id_neighbor, description_neighbor, strand_neighbor,
            start_neighbor, end_neighbor, rel_pos, dist_bp, poly_id, poly_size)]
}

# --- 3) run for high and low --------------------------------------
neighbors_high <- get_neighbors(tcs_high, genes_dt)
neighbors_low  <- get_neighbors(tcs_low,  genes_dt)

# immediate neighbors (rel_pos == -1 or +1)
immediate_high <- neighbors_high[rel_pos %in% c(-1,1)]
immediate_low  <- neighbors_low[rel_pos %in% c(-1,1)]

# --- 4) save output files ----------------------------------
fwrite(neighbors_high, "tcs_high_neighbors.csv")
fwrite(neighbors_low,  "tcs_low_neighbors.csv")
fwrite(immediate_high, "tcs_high_immediate_neighbors.csv")
fwrite(immediate_low,  "tcs_low_immediate_neighbors.csv")


# show an example
print(head(neighbors_high, 10))












# libraries ----------------------------------------------------------------
#if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
#if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
#if (!requireNamespace("gggenes", quietly = TRUE)) install.packages("gggenes")
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#if (!requireNamespace("rtracklayer", quietly = TRUE)) BiocManager::install("rtracklayer")
library(data.table); library(ggplot2); library(gggenes); library(rtracklayer)

# user parameters ---------------------------------------------------
tcs_high_file <- "tcs_high.txt"     # list of gene_id (one column)
n_focal <- 10                        # how many focal TcS (first n)
family_terms <- c("transialidase","Transialidase","trans-sialidase","Trans-sialidase","mucin","Mucin","masp","MASP","GP63","gp63","gp-63","GP-63",
                  "DGF1","dgf1","DGF-1","dgf-1")
# <-- adjust this vector with the terms that define "multigene families"
#      (the search will be case-insensitive and partial)

out_png <- "tcs_high_polys_gggenes.png"

# prepare gene table from tc_gff --------------------------------
# accepts tc_gff as data.frame or GRanges (columns: seqid,start,end,strand,gene_id,description)
if (!exists("tc_gff")) stop("No encuentro 'tc_gff' en el entorno. Cargalo antes de correr este script.")

if (inherits(tc_gff, "GRanges")) {
  genes_dt <- data.table(
    seqid = as.character(seqnames(tc_gff)),
    start = start(tc_gff),
    end   = end(tc_gff),
    strand = as.character(strand(tc_gff)),
    gene_id = mcols(tc_gff)$gene_id,
    description = mcols(tc_gff)$description
  )
} else if (is.data.frame(tc_gff) || data.table::is.data.table(tc_gff)) {
  genes_dt <- as.data.table(tc_gff)[, .(seqid = seqid, start = as.integer(start), end = as.integer(end),
                                        strand = as.character(strand), gene_id = gene_id, description = description)]
} else stop("tc_gff debe ser data.frame o GRanges.")

# order and define polycistrons by run-length of strand within each seqid
setorder(genes_dt, seqid, start)
genes_dt[, poly_rleid := rleid(strand), by = seqid]
genes_dt[, poly_id := paste0(seqid, "_p", poly_rleid)]
genes_dt[, idx_in_poly := seq_len(.N), by = poly_id]
genes_dt[, poly_size := .N, by = poly_id]

# read tcs_high list and select the first n present in tc_gff
tcs_high <- fread(tcs_high_file, header = FALSE)[[1]]
tcs_high_found <- intersect(tcs_high, genes_dt$gene_id)
if (length(tcs_high_found) == 0) stop("Ninguno de los gene_id de tcs_high se encontró en tc_gff.")
focals <- unique(tcs_high_found)[1:min(n_focal, length(tcs_high_found))]

# get poly_id associated to those focals
focal_info <- genes_dt[gene_id %in% focals, .(gene_id, poly_id)]
# usually 1 poly_id per focal gene_id; if duplicates, take all
polys_to_plot <- unique(focal_info$poly_id)

# build dataframe for gggenes: one row per gene within selected polycistrons
# --- reconstruction / corrected plotting ---------------------------------
# Ensure plot_df includes poly_size
plot_df <- genes_dt[poly_id %in% polys_to_plot,
                    .(poly_id, seqid, start, end, strand, gene_id, description, idx_in_poly, poly_size)]

# pattern for families (case-insensitive, partial)
family_pattern <- paste(tolower(family_terms), collapse = "|")

# mark genes of multigene families, focals, and category
plot_df[, is_family := grepl(family_pattern, tolower(ifelse(is.na(description), "", description)), perl=TRUE)]
plot_df[, category := fifelse(gene_id %in% focals, "focal",
                              fifelse(is_family, "multi", "other"))]

# columns for gggenes
plot_df[, molecule := poly_id]
plot_df[, forward  := strand == "+"]

# create facet labels: first focal present in each polycistron + size
label_dt <- genes_dt[poly_id %in% polys_to_plot & gene_id %in% focals,
                     .(focal = gene_id[1L]), by = poly_id]
# if any polycistron has no focal recorded (rare), assign NA
sizes <- unique(plot_df[, .(poly_id, poly_size)])
label_dt <- merge(data.table(poly_id = polys_to_plot), label_dt, by = "poly_id", all.x = TRUE)
label_dt <- merge(label_dt, sizes, by = "poly_id", all.x = TRUE)
label_dt[, facet_label := paste0(ifelse(is.na(focal), "(no focal)", focal), "  |  ", poly_id, "  (n=", poly_size, ")")]

plot_df <- merge(plot_df, label_dt[, .(poly_id, facet_label)], by = "poly_id", all.x = TRUE)







library(data.table); library(ggplot2); library(gggenes); library(rtracklayer)

# user parameters for low ---------------------------------------------------
tcs_low_file <- "tcs_low.txt"     # list of gene_id (one column)
n_focal <- 10                        # how many focal TcS (first n)
family_terms <- c("transialidase","Transialidase","trans-sialidase","Trans-sialidase","mucin","Mucin","masp","MASP","GP63","gp63","gp-63","GP-63",
                  "DGF1","dgf1","DGF-1","dgf-1")
# <-- adjust this vector with the terms that define "multigene families"

out_png <- "tcs_low_polys_gggenes.png"

# prepare gene table from tc_gff --------------------------------
# accepts tc_gff as data.frame or GRanges (columns: seqid,start,end,strand,gene_id,description)
if (!exists("tc_gff")) stop("No encuentro 'tc_gff' en el entorno. Cargalo antes de correr este script.")

if (inherits(tc_gff, "GRanges")) {
  genes_dt <- data.table(
    seqid = as.character(seqnames(tc_gff)),
    start = start(tc_gff),
    end   = end(tc_gff),
    strand = as.character(strand(tc_gff)),
    gene_id = mcols(tc_gff)$gene_id,
    description = mcols(tc_gff)$description
  )
} else if (is.data.frame(tc_gff) || data.table::is.data.table(tc_gff)) {
  genes_dt <- as.data.table(tc_gff)[, .(seqid = seqid, start = as.integer(start), end = as.integer(end),
                                        strand = as.character(strand), gene_id = gene_id, description = description)]
} else stop("tc_gff debe ser data.frame o GRanges.")

# order and define polycistrons by run-length of strand within each seqid
setorder(genes_dt, seqid, start)
genes_dt[, poly_rleid := rleid(strand), by = seqid]
genes_dt[, poly_id := paste0(seqid, "_p", poly_rleid)]
genes_dt[, idx_in_poly := seq_len(.N), by = poly_id]
genes_dt[, poly_size := .N, by = poly_id]

# read tcs_low and select first n that are in tc_gff
tcs_low <- fread(tcs_low_file, header = FALSE)[[1]]
tcs_low_found <- intersect(tcs_low, genes_dt$gene_id)
if (length(tcs_low_found) == 0) stop("Ninguno de los gene_id de tcs_high se encontró en tc_gff.")
focals <- unique(tcs_low_found)[1:min(n_focal, length(tcs_low_found))]

# get poly_id associated to those focals
focal_info <- genes_dt[gene_id %in% focals, .(gene_id, poly_id)]
polys_to_plot <- unique(focal_info$poly_id)

# build dataframe for gggenes
plot_df <- genes_dt[poly_id %in% polys_to_plot,
                    .(poly_id, seqid, start, end, strand, gene_id, description, idx_in_poly, poly_size)]

# family pattern
family_pattern <- paste(tolower(family_terms), collapse = "|")

# mark multigene families, focals, category
plot_df[, is_family := grepl(family_pattern, tolower(ifelse(is.na(description), "", description)), perl=TRUE)]
plot_df[, category := fifelse(gene_id %in% focals, "focal",
                              fifelse(is_family, "multi", "other"))]

# columns for gggenes
plot_df[, molecule := poly_id]
plot_df[, forward  := strand == "+"]

# facet labels
label_dt <- genes_dt[poly_id %in% polys_to_plot & gene_id %in% focals,
                     .(focal = gene_id[1L]), by = poly_id]
sizes <- unique(plot_df[, .(poly_id, poly_size)])
label_dt <- merge(data.table(poly_id = polys_to_plot), label_dt, by = "poly_id", all.x = TRUE)
label_dt <- merge(label_dt, sizes, by = "poly_id", all.x = TRUE)
label_dt[, facet_label := paste0(ifelse(is.na(focal), "(no focal)", focal), "  |  ", poly_id, "  (n=", poly_size, ")")]

plot_df <- merge(plot_df, label_dt[, .(poly_id, facet_label)], by = "poly_id", all.x = TRUE)


# ---------------------------
# Improved aesthetics for gggenes plot
# ---------------------------

library(ggplot2)
library(gggenes)
library(grid)   # for unit()
library(scales) # for comma_format()

# colors (you can tweak these)
fill_colors <- c(focal = "yellow", multi = "salmon", other = "darkgreen") 

# make molecule an ordered factor so tracks stack consistently
plot_df[, molecule := factor(molecule, levels = unique(plot_df$molecule))]

# adjust facet labels (already computed as label_dt -> facet_label)
# ensure facet_label exists in plot_df
if (!"facet_label" %in% colnames(plot_df)) {
  plot_df <- merge(plot_df, label_dt[, .(poly_id, facet_label)], by.x = "poly_id", by.y = "poly_id", all.x = TRUE)
}

# Basic aesthetics parameters to tweak appearance
arrow_head_mm <- 5         # arrowhead size in mm
arrow_height  <- 0.6       # arrow body height (relative)
arrow_radius  <- unit(0.6, "mm")
label_size_focal <- 3.2
axis_text_x_size <- 9

# Build plot
p <- ggplot(plot_df, aes(xmin = start, xmax = end, y = molecule, fill = category, forward = forward)) +
  # gene arrows: thinner stroke, larger arrowheads, slightly taller arrows
  geom_gene_arrow(
    arrowhead_height = unit(arrow_head_mm, "mm"),
    arrowhead_width  = unit(arrow_head_mm, "mm"),
    height = arrow_height,
    radius = arrow_radius,
    colour = "grey30",
    size = 0.22
  ) +
  # label focals only (contrast: white text on dark fill); if many focals, reduce size
  #geom_gene_label(
  #  data = plot_df[category == "focal"],
  #  aes(label = gene_id),
  #  align = "center",
  #  colour = "white",
  #  fontface = "bold",
  #  size = label_size_focal,
  #  max.ratio = 0.35
  #) +
  # fill colors for categories
  scale_fill_manual(values = fill_colors, breaks = c("focal", "multi", "other")) +
  # control x axis formatting (commas) and remove expansion to avoid big margins
  scale_x_continuous(expand = c(0, 0), labels = comma_format(accuracy = 1)) +
  # reduce white space between stacked tracks, but increase facet spacing
  facet_wrap(~ facet_label, scales = "free_x", ncol = 1) +
  # subtle horizontal baseline: emulate the thin black line below each track
  geom_hline(aes(yintercept = as.numeric(molecule)), data = plot_df[, .(molecule = unique(molecule))], 
             colour = "grey20", size = 0.25, inherit.aes = FALSE) +
  # theme tweaks to match figure style
  theme_minimal(base_size = 12) +
  theme(
    # remove grid, keep a clean baseline look
    panel.grid = element_blank(),
    # more vertical space between facets
    panel.spacing = unit(1.0, "lines"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = axis_text_x_size),
    axis.title.x = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    plot.margin = unit(c(8, 8, 8, 8), "pt")
  ) +
  labs(x = "Position", y = NULL, fill = "Category")

# Save PNG sized depending on number of polys
n_polys <- length(polys_to_plot)
png(out_png, width = 12, height = max(2.5 * n_polys, 4), units = "in", res = 300)
print(p)
dev.off()








#############################################################################################
###################################### Figure 5B #############################################
#############################################################################################


##############################################
# Corrected script: % multigene neighbors per gene
# Using ALL neighbors within the polycistron (k = Inf)
# - computes per gene for TcS high (each gene)
# - draws a sample of 30 TcS low (once) and computes per gene
# - uses unique neighbors per gene to avoid >100%
# - saves CSVs and comparative plot
##############################################

library(data.table)
library(ggplot2)

# ------------ parameters (adjust paths if needed) ----------------
tcs_low_file  <- "tcs_low_all.txt"
tcs_high_file <- "tcs_high.txt"

sample_size   <- 30   # how many genes to take from tcs_low (single sample)
neighbor_k    <- Inf  # <-- use Inf to take ALL neighbors from the same polycistron

out_dir <- "localizacion_genomica"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# pattern for multigene families (ensure 'family_terms' exists: vector of strings)
if (!exists("family_terms")) stop("No encuentro 'family_terms' en el entorno. Define family_terms como vector de strings.")
family_pattern <- paste(family_terms, collapse = "|")
family_pattern <- paste0("(?i)", family_pattern) # case-insensitive / PCRE

# seed for reproducibility
set.seed(123)

# --------------------------------------------------------------------
# Initial checks
# --------------------------------------------------------------------
if (!exists("genes_dt")) stop("No encuentro 'genes_dt' en el entorno. Corre la sección que lo construye.")
required_cols <- c("seqid","start","end","strand","gene_id","description","poly_id","idx_in_poly","poly_size")
stopifnot(all(required_cols %in% colnames(genes_dt)))

# --- read lists --------------------------------
tcs_low_all <- fread(tcs_low_file, header = FALSE)[[1]]
tcs_high   <- fread(tcs_high_file, header = FALSE)[[1]]

# keep only those present in genes_dt
tcs_low_all_found  <- intersect(tcs_low_all,  genes_dt$gene_id)
tcs_high_found     <- intersect(tcs_high,     genes_dt$gene_id)

if (length(tcs_low_all_found) < 1) stop("No se encontraron TcS low en genes_dt.")
if (length(tcs_high_found) < 1) stop("No se encontraron TcS high en genes_dt.")

# if fewer available than sample_size -> sample with replacement and warn
replace_flag <- FALSE
if (length(tcs_low_all_found) < sample_size) {
  warning(sprintf("Solo %d genes disponibles en tcs_low_all_found < sample_size=%d. Muestrearé con reemplazo.", 
                  length(tcs_low_all_found), sample_size))
  replace_flag <- TRUE
}

# --------------------------------------------------------------------
# Robust functions (unique neighbors)
# --------------------------------------------------------------------

# returns unique pairs (gene_id_focal, gene_id_neighbor) within the polycistron according to k
get_neighbors_dt_unique <- function(focal_genes, genes_dt, k = 1) {
  focals <- genes_dt[gene_id %in% focal_genes, .(poly_id, gene_id_focal = gene_id, idx_focal = idx_in_poly)]
  if (nrow(focals) == 0) return(data.table())
  neigh <- merge(focals, genes_dt, by = "poly_id", allow.cartesian = TRUE)
  neigh[, rel_pos := idx_in_poly - idx_focal]
  if (is.infinite(k)) {
    # all members of the same poly_id except the focal
    neigh <- neigh[gene_id != gene_id_focal]
  } else {
    neigh <- neigh[gene_id != gene_id_focal & abs(rel_pos) <= k]
  }
  # rename relevant columns
  setnames(neigh,
           old = c("gene_id","start","end","strand","description","idx_in_poly","poly_size","seqid"),
           new = c("gene_id_neighbor","start_neighbor","end_neighbor","strand_neighbor",
                   "description_neighbor","idx_neighbor","poly_size","seqid_neighbor"))
  # keep unique pairs (focal x neighbor)
  neigh_unique <- unique(neigh, by = c("gene_id_focal", "gene_id_neighbor"))
  # useful columns
  neigh_unique[, .(gene_id_focal, poly_id, idx_focal, idx_neighbor, gene_id_neighbor,
                   description_neighbor, rel_pos, poly_size)]
}

# compute pct per gene using unique neighbors; also returns n_neighbors and n_family
pct_multigenic_per_gene_unique <- function(focal_genes, genes_dt, family_pattern, k = 1) {
  ndt <- get_neighbors_dt_unique(focal_genes, genes_dt, k = k)
  # if no rows -> return NA for each focal
  if (nrow(ndt) == 0) {
    return(data.table(gene_id_focal = focal_genes, pct = NA_real_, n_neighbors = 0L, n_family = 0L))
  }
  # mark if neighbor belongs to a multigene family
  ndt[, is_family := grepl(family_pattern, ifelse(is.na(description_neighbor), "", description_neighbor), perl = TRUE)]
  # compute per gene using unique neighbors
  res <- ndt[, .(
    n_neighbors = uniqueN(gene_id_neighbor),           # number of unique neighbors
    n_family    = sum(is_family, na.rm = TRUE)         # how many of those were marked family
  ), by = .(gene_id_focal)]
  # safe pct: if n_neighbors == 0 -> NA
  res[, pct := ifelse(n_neighbors > 0, 100 * n_family / n_neighbors, NA_real_)]
  # add missing focal genes
  missing_genes <- setdiff(focal_genes, res$gene_id_focal)
  if (length(missing_genes) > 0) {
    res <- rbind(res, data.table(gene_id_focal = missing_genes, n_neighbors = 0L, n_family = 0L, pct = NA_real_))
  }
  # preserve input order
  res[, gene_id_focal := factor(gene_id_focal, levels = unique(focal_genes))]
  setorder(res, gene_id_focal)
  res[, gene_id_focal := as.character(gene_id_focal)]
  return(res[, .(gene_id_focal, pct, n_neighbors, n_family)])
}

# --------------------------------------------------------------------
# 1) High: pct per gene (no sampling) using whole polycistron
# --------------------------------------------------------------------
high_pct_dt <- pct_multigenic_per_gene_unique(tcs_high_found, genes_dt, family_pattern, k = neighbor_k)

# --------------------------------------------------------------------
# 2) Low: sample sample_size genes randomly (once) and compute pct per gene using whole polycistron
# --------------------------------------------------------------------
samp_low <- sample(tcs_low_all_found, size = sample_size, replace = replace_flag)
low_pct_dt <- pct_multigenic_per_gene_unique(samp_low, genes_dt, family_pattern, k = neighbor_k)

# --------------------------------------------------------------------
# 3) Save results (with 'polycistron' suffix to clarify)
# --------------------------------------------------------------------
fwrite(high_pct_dt, file.path(out_dir, "high_pct_per_gene_polycistron_unique.csv"))
fwrite(low_pct_dt,  file.path(out_dir, "low_sample30_pct_per_gene_polycistron_unique.csv"))

cat("Saved:\n -", file.path(out_dir, "high_pct_per_gene_polycistron_unique.csv"), "\n -", file.path(out_dir, "low_sample30_pct_per_gene_polycistron_unique.csv"), "\n")

# extra check: ensure no pct > 100
bad_high <- high_pct_dt[!is.na(pct) & pct > 100]
bad_low  <- low_pct_dt[!is.na(pct) & pct > 100]
if (nrow(bad_high) > 0 || nrow(bad_low) > 0) {
  warning("Se detectaron valores pct > 100 después de la corrección. Imprimiendo ejemplos para investigar.")
  if (nrow(bad_high) > 0) { cat("High con pct>100:\n"); print(bad_high) }
  if (nrow(bad_low) > 0)  { cat("Low con pct>100:\n");  print(bad_low)  }
} else {
  cat("Check: no pct > 100 detected (OK).\n")
}

# additional: check for repeated gene_id in genes_dt (may be informative)
dups <- genes_dt[, .N, by = gene_id][N > 1]
if (nrow(dups) > 0) {
  cat("Warning: duplicated gene_id in genes_dt (show up to 20):\n")
  print(head(dups, 20))
} else {
  cat("No duplicated gene_id detected in genes_dt (simple count).\n")
}

# --------------------------------------------------------------------
# 4) Statistical test: compare per-gene distribution high vs low sample30
# --------------------------------------------------------------------
high_vals <- na.omit(high_pct_dt$pct)
low_vals  <- na.omit(low_pct_dt$pct)

if (length(high_vals) < 1) stop("No hay valores válidos para 'high' (todos NA).")
if (length(low_vals) < 1) stop("No hay valores válidos para 'low' (todos NA).")

wilcox_res <- wilcox.test(high_vals, low_vals, alternative = "two.sided")
t_res <- t.test(high_vals, low_vals)

cat("\nWilcoxon test:\n"); print(wilcox_res)
cat("\nT-test (complementary):\n"); print(t_res)

# --------------------------------------------------------------------
# 5) Plot: comparative boxplot (High per-gene vs Low sample30) - all polycistron neighbors
# --------------------------------------------------------------------
plot_high_dt <- data.table(group = "High (each gene)", gene_id = high_pct_dt$gene_id_focal, pct = high_pct_dt$pct)
plot_low_dt  <- data.table(group = "Low (sampled 30)", gene_id = low_pct_dt$gene_id_focal, pct = low_pct_dt$pct)

plot_dt <- rbind(plot_high_dt, plot_low_dt)
plot_dt <- plot_dt[!is.na(pct)]

pal <- c("High (each gene)" = "#d95f02", "Low (sampled 30)" = "#1b9e77")

p <- ggplot(plot_dt, aes(x = group, y = pct, fill = group)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = group), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Percent multigene neighbors per gene: High vs Low (all polycistron neighbors)",
    subtitle = sprintf("Low: sample size = %d (single sample). Wilcoxon p = %s", sample_size, signif(wilcox_res$p.value, 3)),
    x = "", y = "% multigene neighbors"
  ) +
  theme(legend.position = "none")

out_png <- file.path(out_dir, "pct_multigenic_high_vs_low_sample30_polycistron_boxplot.png")
ggsave(out_png, p, width = 8, height = 5, dpi = 300)

print(p)
cat("Plot saved to:", out_png, "\n")











#############################################################################################
###################################### Figure 5C #############################################
#############################################################################################


##############################################
# Low iterations: 50 x 30, averages per iter
# Boxplot of the distribution of averages (low)
# + point with the high average (color #d95f02)
##############################################

library(data.table)
library(ggplot2)

# ---------- parameters ----------
tcs_low_file  <- "tcs_low_all.txt"
tcs_high_file <- "tcs_high.txt"

n_iter      <- 50
sample_size <- 30
neighbor_k  <- Inf   # use all polycistron neighbors

out_dir <- "localizacion_genomica"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- checks ----------
if (!exists("genes_dt")) stop("No encuentro 'genes_dt' en el entorno. Crea genes_dt primero.")
if (!exists("family_terms")) stop("No encuentro 'family_terms' en el entorno. Define family_terms como vector de strings.")

required_cols <- c("seqid","start","end","strand","gene_id","description","poly_id","idx_in_poly","poly_size")
stopifnot(all(required_cols %in% colnames(genes_dt)))

# read lists
tcs_low_all <- fread(tcs_low_file, header = FALSE)[[1]]
tcs_high   <- fread(tcs_high_file, header = FALSE)[[1]]

tcs_low_all_found  <- intersect(tcs_low_all, genes_dt$gene_id)
tcs_high_found     <- intersect(tcs_high,    genes_dt$gene_id)

if (length(tcs_low_all_found) < 1) stop("No se encontraron TcS low en genes_dt.")
if (length(tcs_high_found) < 1) stop("No se encontraron TcS high en genes_dt.")

# replacement if not enough genes
replace_flag <- length(tcs_low_all_found) < sample_size
if (replace_flag) warning("Habrá muestreo con reemplazo para low: sample_size > nº disponibles.")

family_pattern <- paste(family_terms, collapse = "|")
family_pattern <- paste0("(?i)", family_pattern)

set.seed(123)

# ---------- functions (unique neighbors per polycistron) ----------
get_neighbors_dt_unique <- function(focal_genes, genes_dt, k = 1) {
  focals <- genes_dt[gene_id %in% focal_genes, .(poly_id, gene_id_focal = gene_id, idx_focal = idx_in_poly)]
  if (nrow(focals) == 0) return(data.table())
  neigh <- merge(focals, genes_dt, by = "poly_id", allow.cartesian = TRUE)
  neigh[, rel_pos := idx_in_poly - idx_focal]
  if (is.infinite(k)) {
    neigh <- neigh[gene_id != gene_id_focal]
  } else {
    neigh <- neigh[gene_id != gene_id_focal & abs(rel_pos) <= k]
  }
  setnames(neigh,
           old = c("gene_id","start","end","strand","description","idx_in_poly","poly_size","seqid"),
           new = c("gene_id_neighbor","start_neighbor","end_neighbor","strand_neighbor",
                   "description_neighbor","idx_neighbor","poly_size","seqid_neighbor"))
  neigh_unique <- unique(neigh, by = c("gene_id_focal", "gene_id_neighbor"))
  neigh_unique[, .(gene_id_focal, poly_id, idx_focal, idx_neighbor, gene_id_neighbor,
                   description_neighbor, rel_pos, poly_size)]
}

pct_multigenic_per_gene_unique <- function(focal_genes, genes_dt, family_pattern, k = 1) {
  ndt <- get_neighbors_dt_unique(focal_genes, genes_dt, k = k)
  if (nrow(ndt) == 0) {
    return(data.table(gene_id_focal = focal_genes, pct = NA_real_, n_neighbors = 0L, n_family = 0L))
  }
  ndt[, is_family := grepl(family_pattern, ifelse(is.na(description_neighbor), "", description_neighbor), perl = TRUE)]
  res <- ndt[, .(
    n_neighbors = uniqueN(gene_id_neighbor),
    n_family    = sum(is_family, na.rm = TRUE)
  ), by = .(gene_id_focal)]
  res[, pct := ifelse(n_neighbors > 0, 100 * n_family / n_neighbors, NA_real_)]
  missing_genes <- setdiff(focal_genes, res$gene_id_focal)
  if (length(missing_genes) > 0) {
    res <- rbind(res, data.table(gene_id_focal = missing_genes, n_neighbors = 0L, n_family = 0L, pct = NA_real_))
  }
  res[, gene_id_focal := factor(gene_id_focal, levels = unique(focal_genes))]
  setorder(res, gene_id_focal)
  res[, gene_id_focal := as.character(gene_id_focal)]
  return(res[, .(gene_id_focal, pct, n_neighbors, n_family)])
}

# ---------- compute mean for high (per gene then global mean) ----------
high_pct_dt <- pct_multigenic_per_gene_unique(tcs_high_found, genes_dt, family_pattern, k = neighbor_k)
high_mean_pct <- mean(high_pct_dt$pct, na.rm = TRUE)

# ---------- iterations for low ----------
iter_results <- vector("list", n_iter)
for (it in seq_len(n_iter)) {
  samp <- sample(tcs_low_all_found, size = sample_size, replace = replace_flag)
  tmp_dt <- pct_multigenic_per_gene_unique(samp, genes_dt, family_pattern, k = neighbor_k)
  mean_it <- mean(tmp_dt$pct, na.rm = TRUE)  # iteration mean (exclude NA)
  iter_results[[it]] <- data.table(iter = it, mean_pct = mean_it, n_valid = sum(!is.na(tmp_dt$pct)))
}
iter_dt <- rbindlist(iter_results)

# remove iterations with NA mean
valid_iter_dt <- iter_dt[!is.na(mean_pct)]
if (nrow(valid_iter_dt) < n_iter) {
  warning(sprintf("%d iterations had mean_pct = NA and were excluded from the plot.", n_iter - nrow(valid_iter_dt)))
}

# save CSV with means per iteration
fwrite(valid_iter_dt, file.path(out_dir, "low_iter50_meanpct_sample30_per_iter.csv"))

# ---------- plot: boxplot of the 50 means (low) + high point ----------
plot_dt <- copy(valid_iter_dt)
plot_dt[, group := "Low (iter means)"]

# ---------- plot ----------
pal_low  <- "#1b9e77"  # green
pal_high <- "#d95f02"  # orange

p <- ggplot(plot_dt, aes(x = group, y = mean_pct)) +
  geom_boxplot(width = 0.4, fill = pal_low, alpha = 0.4, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.12, color = pal_low, size = 2, alpha = 0.8) +
  # point for high mean (x = 1)
  geom_point(aes(x = 1, y = high_mean_pct), color = pal_high, size = 2, fill = pal_high) +
  annotate("text", x = 1.12, y = high_mean_pct, 
           label = sprintf("High mean = %.2f%%", high_mean_pct),
           hjust = 0, vjust = -0.3, color = pal_high, size = 3.5) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = sprintf("Distribution of iteration means (Low: %d iterations of %d genes)", n_iter, sample_size),
    subtitle = "Boxplot of the means (each iteration averages % per gene). Orange point = High mean",
    x = "", y = "Mean % multigene neighbors (per iteration)"
  )

out_png <- file.path(out_dir, "low_iter50_means_boxplot_with_high_point_colored.png")
ggsave(out_png, p, width = 7, height = 5, dpi = 300)

print(p)



## RAINCLOUD PLOT


# Raincloud plot for Low iterations + High point with empirical statistics
# Requires: data.table, ggplot2, ggdist, dplyr

library(data.table); library(ggplot2); library(ggdist); library(dplyr)

# --- Assume these exist in the environment:
# valid_iter_dt (cols: iter, mean_pct), high_mean_pct (numeric), out_dir (string)
# If not, read the CSV:
# valid_iter_dt <- fread(file.path(out_dir, "low_iter50_meanpct_sample30_per_iter.csv"))

plot_dt <- copy(valid_iter_dt)
plot_dt[, group := "Low (iter means)"]

pal_low  <- "#1b9e77"  # green
pal_high <- "#d95f02"  # orange

# empirical statistics: percentile and p-value (one-sided and two-sided)
n_iter_obs <- nrow(plot_dt)
n_ge <- sum(plot_dt$mean_pct >= high_mean_pct, na.rm = TRUE)
p_emp_one_sided <- (n_ge + 1) / (n_iter_obs + 1)           # +1 pseudo-count (recommended)
p_emp_two_sided <- 2 * min(p_emp_one_sided, 1 - p_emp_one_sided + 1/(n_iter_obs+1))

percentile_high <- ecdf(plot_dt$mean_pct)(high_mean_pct) * 100  # percentile (0-100)

# positions for annotations (dynamic)
y_rng <- range(plot_dt$mean_pct, na.rm = TRUE)
y_max <- y_rng[2]
y_annot <- y_max + 0.06 * diff(y_rng)  # adjust depending on space

# Build raincloud
p <- ggplot(plot_dt, aes(x = group, y = mean_pct, fill = group)) +
  # half-violin (stat_halfeye)
  ggdist::stat_halfeye(
    aes(y = mean_pct), 
    adjust = 0.6, 
    width = 0.6, 
    justification = -0.2, 
    .width = c(0, 0.5, 0.9),
    slab_alpha = 0.6,
    point_colour = NA
  ) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.25, color = "black") +
  geom_jitter(width = 0.08, size = 1.8, alpha = 0.8, color = pal_low) +
  # point for high (x = 1)
  geom_point(aes(x = 1, y = high_mean_pct), color = pal_high, size = 3) +
  geom_point(aes(x = 1, y = high_mean_pct), color = "black", size = 0.6) + # outline
  # text annotations: percentile and empirical p-value
  annotate("text", x = 1.15, y = y_annot, 
           label = sprintf("High mean = %.2f %%", high_mean_pct),
           hjust = 0, vjust = 0, color = pal_high, size = 3.8) +
  coord_cartesian(ylim = c(0,
                           y_annot + 0.02*diff(y_rng))) +
  scale_fill_manual(values = pal_low) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = sprintf("Raincloud: Low (n=%d iter) vs High (single mean)", n_iter_obs),
    y = "Mean % multigene neighbors (per iteration)"
  )

# Save high-res PNG
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_png_rain <- file.path(out_dir, "low_iter50_raincloud_with_high_point.png")
ggsave(out_png_rain, p, width = 6.5, height = 5, dpi = 300)

# Show plot in R
print(p)


