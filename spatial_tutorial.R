# Import Statements
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

# Loading Dataset (h5 & spatial)
gbm <- Load10X_Spatial('C:/Users/dchav/Documents/10X_Genomics/Parent_Visium_Human_Glioblastoma/')

# Data PreProcessing

# These plots show how the anatomy of the tissue affects the molecular counts
plot1 <- VlnPlot(gbm, features = 'nCount_Spatial', pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(gbm, features = 'nCount_Spatial') + theme(legend.position = 'right')
wrap_plots(plot1, plot2)

# Normalization with SCTransform() over LogNormalize(), as it preserves
# biological variance, while LogNormalize() doesn't.

gbm <- SCTransform(gbm, assay = 'Spatial', verbose = FALSE) # high-variance features are stored in the SCT assay

# Gene Expression Visualization

# SpatialFeaturePlot() is an extension of FeaturePlot(), and can overlay
# molecular data on top of tissue histology.

SpatialFeaturePlot(gbm, features = c('MGMT', 'IDH1')) # plotted gbm marker genes

# Dimensionality Reduction, Clustering, and Visualization

# Dimensionality reduction and PCA are the same as for scRNA-seq analysis

gbm <- RunPCA(gbm, assay = "SCT", verbose = FALSE)
gbm <- FindNeighbors(gbm, reduction = "pca", dims = 1:30)
gbm <- FindClusters(gbm, verbose = FALSE)
gbm <- RunUMAP(gbm, reduction = "pca", dims = 1:30)

# Visulize results from clusterin in the UMAP space via DimPlot()
# or overlaid on the image via SpatialDimPlot()

p1 <- DimPlot(gbm, reduction = 'umap', label = TRUE)
p2 <- SpatialDimPlot(gbm, label = TRUE, label.size = 3)
p1 + p2

# To see which voxel belongs to each clusterm you can use the label parameter
# to place a colored box at the medium of each cluster. It is also possible to
# use the cells.highlight parameter to mark particular cells of interest on a
# SpatialDimPlot(). This is useful for distinguishing spatial localization of 
# individual clusters

SpatialDimPlot(gbm, cells.highlight = CellsByIdentities(object = gbm, idents = c(0, 1, 2)), facet.highlight = TRUE, ncol = 3)

# Settings interactive to TRUE allows for an interactive pane to adjust
# spot transparency, as well as the assay and feature being plotted.

SpatialFeaturePlot(gbm, features = 'IDH1', interactive = TRUE)

# LinkedDimPlot() links the UMAP representation of the tissue to the tissue
# image representation and allows for interactive selection

LinkedDimPlot(gbm)

# Identification of Spatially Variable Features

# Seurat has two main ways of identifying features that correlate with spatial
# location within a tissue. The first is to perform differential expression
# based on pre-annotated anatomical regions within the tissue, which may be
# determined from unsupervised clustering or prior knowledge. This works well
# when the clusters exhbiti spatial restriction,

de_markers <- FindMarkers(gbm, ident.1 = 0, ident.2 = 7)
SpatialFeaturePlot(object = gbm, features = rownames(de_markers)[1:3], alpha = c(.1, 1), ncol = 3)

# Another approach is using FindSpatiallyVariableFeatures(), which searchs for
# features that show spatial patterning in the absence of pre-annotation.

# gbm <- FindSpatiallyVariableFeatures(gbm, assay = 'SCT', features = VariableFeatures(gbm[1:1000], selection.method = 'moransi'))

# Visualizing the expression of the top 6 features identified

# top.features <- head(SpatiallyVariableFeatures(gbm, selection.method = 'moransi'), 6)
# SpatialFeaturePlot(gbm, features = top.features, ncol = 3, alpha = c(0.1, 1))

# Subset Out Anatomical Regions

# With single-cell objects, it is possible to subset the object to focus on a 
# subset of data. 









