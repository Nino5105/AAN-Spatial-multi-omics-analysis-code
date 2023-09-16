# 1.Loading packages

library(Seurat)
library(RCurl)
library(cowplot)
library(dplyr)
library(hdf5r) # install.packages('hdf5r')
options(future.globals.maxSize = 10000 * 1024^2)  # set allowed size to 2K MiB

# 2.Loadng ST-seq datasets

Control_1 <- Load10X_Spatial(
    data.dir = "Control_1/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "Control_1",
    filter.matrix = TRUE,
    to.upper = FALSE)

Control_2 <- Load10X_Spatial(
    data.dir = "Control_2/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "Control_2",
    filter.matrix = TRUE,
    to.upper = FALSE)

AAN_2W_1 <- Load10X_Spatial(
    data.dir = "AAN_2W_1/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "AAN_2W_1",
    filter.matrix = TRUE,
    to.upper = FALSE)

AAN_2W_2 <- Load10X_Spatial(
    data.dir = "AAN_2W_2/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "AAN_2W_2",
    filter.matrix = TRUE,
    to.upper = FALSE)

AAN_4W_1 <- Load10X_Spatial(
    data.dir = "AAN_4W_1/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "AAN_4W_1",
    filter.matrix = TRUE,
    to.upper = FALSE)

AAN_4W_2 <- Load10X_Spatial(
    data.dir = "AAN_4W_2/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "AAN_4W_2",
    filter.matrix = TRUE,
    to.upper = FALSE)


# 3. Add metadata information

Control_1$sample <- "Control_1"
Control_2$sample <- "Control_2"
AAN_2W_1$sample <- "AAN_2W_1"
AAN_2W_2$sample <- "AAN_2W_2"
AAN_4W_1$sample <- "AAN_4W_1"
AAN_4W_2$sample <- "AAN_4W_2"

Control_1$type <- "Control"
Control_2$type <- "Control"
AAN_2W_1$type <- "AAN_2W"
AAN_2W_2$type <- "AAN_2W"
AAN_4W_1$type <- "AAN_4W"
AAN_4W_2$type <- "AAN_4W"

# 4. Add percent_mito information

Control_1 <- PercentageFeatureSet(Control_1, "^mt-", col.name = "percent_mito")
Control_2 <- PercentageFeatureSet(Control_2, "^mt-", col.name = "percent_mito")
AAN_2W_1 <- PercentageFeatureSet(AAN_2W_1, "^mt-", col.name = "percent_mito")
AAN_2W_2 <- PercentageFeatureSet(AAN_2W_2, "^mt-", col.name = "percent_mito")
AAN_4W_1 <- PercentageFeatureSet(AAN_4W_1, "^mt-", col.name = "percent_mito")
AAN_4W_2 <- PercentageFeatureSet(AAN_4W_2, "^mt-", col.name = "percent_mito")

# 5. Merge datasets

merge_data <- merge(Control_1,c(Control_2,AAN_2W_1,AAN_2W_2,AAN_4W_1,AAN_4W_2))
merge_data

# An object of class Seurat
# 32285 features across 18817 samples within 1 assay

merge_data$sample = factor(merge_data$sample, levels = c("Control_1", "Control_2","AAN_2W_1", "AAN_2W_2","AAN_4W_1", "AAN_4W_2"))
merge_data$type = factor(merge_data$type, levels = c("Control", "AAN_2W", "AAN_4W"))

# 6. Datasets QC

pdf("merge_data_before_filter.pdf",height = 4,width = 10)
VlnPlot(merge_data, features = c("nFeature_Spatial","nCount_Spatial","percent_mito"), pt.size = 0,group.by = "sample",ncol = 6)
dev.off()

merge_data.filter <- merge_data[, merge_data$nFeature_Spatial > 500 & merge_data$percent_mito < 30]
merge_data.filter

# An object of class Seurat
# 32285 features across 18673 samples within 1 assay

pdf("merge_data_after_filter-sample.pdf",height = 4,width = 10)
VlnPlot(merge_data.filter, features = c("nFeature_Spatial","nCount_Spatial","percent_mito"), pt.size = 0,group.by = "sample",ncol = 6)
dev.off()

# 7.Integrate Datasets 

split_data <- SplitObject(merge_data.filter, split.by = "sample")
st.list = lapply(split_data, SCTransform, assay = "Spatial", method = "poisson")
st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features, verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT", verbose = FALSE, anchor.features = st.features)
integrated_data <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

# 8.Datasets reduction 
integrated_data <- RunPCA(object = integrated_data,verbose = FALSE)
pdf("PCA_Elbowplot.pdf",height=10,width=10)
ElbowPlot(integrated_data,ndims = 50)
dev.off()
integrated_data <- FindNeighbors(integrated_data,dim=1:40)
integrated_data <- FindClusters(integrated_data,resolution = 0.8)
integrated_data <- RunUMAP (integrated_data,reduction="pca", dims = 1:40)
integrated_data <- RunTSNE(integrated_data,dims = 1:40)

# 9.Dataset visulization
pdf("umap-cluster.pdf",height = 10,width = 10)
DimPlot(integrated_data,reduction = "umap" ,label = F,pt.size = 1)
dev.off()

pdf("umap-sample.pdf",height = 10,width = 10)
DimPlot(integrated_data,reduction = "umap" ,group.by="sample",label = F,pt.size = 1)
dev.off()

pdf("umap-type.pdf",height = 10,width = 10)
DimPlot(integrated_data,reduction = "umap" ,group.by="type",label = T,pt.size = 1)
dev.off()

# 10.save datasets
saveRDS(integrated_data,"integrated_data.rds")