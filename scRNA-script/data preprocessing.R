library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(harmony)
library(dplyr)
setwd('E:/zhh_script')
#read data
for (file in c("zhh1_filtered_feature_bc_matrix", "zhh2_filtered_feature_bc_matrix"))
{
  seurat_data <- Read10X(data.dir = file)
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.cells=3,
                                   min.features = 200, 
                                   project = file)
  assign(file, seurat_obj)
}
#merge data
merged_seurat <- merge(x = zhh1_filtered_feature_bc_matrix, 
                       y = zhh2_filtered_feature_bc_matrix, 
                       add.cell.id = c("zhh1", "zhh2"))
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$rbRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^RP[SL]")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
merged_seurat$rbRatio <- merged_seurat@meta.data$rbRatio /100
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>% dplyr::rename(seq_folder = orig.ident, 
                                       nUMI = nCount_RNA, 
                                       nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^zhh1_"))] <- "zhh1"
metadata$sample[which(str_detect(metadata$cells, "^zhh2_"))] <- "zhh2"
merged_seurat@meta.data <- metadata
#quality control
filtered_seurat <- subset(x = merged_seurat, subset= (nUMI >= 500) &
                            (nGene >= 200) &((nGene <= 6000) )&
                            (log10GenesPerUMI > 0.80) &
                            (mitoRatio < 0.20) &
                            (rbRatio < 0.30 ))
#normalizeData
seuratdata<-NormalizeData(filtered_seurat, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)
seuratdata <- ScaleData(seuratdata)
#PCA
seuratdata <- RunPCA(seuratdata, 
                     features = VariableFeatures(object = seuratdata))
VizDimLoadings(seuratdata, dims = 6, 
               reduction = "pca") + 
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=5,face="bold"))
DimPlot(seuratdata, reduction = "pca")
ElbowPlot(seuratdata,40)
#harmony
seuratdata <- RunHarmony(seuratdata, 
                         group.by.vars = c("sample"), 
                         reduction = "pca", 
                         assay.use = "SCT", 
                         reduction.save = "harmony")
#UMAP
seuratdata <- RunUMAP(seuratdata, 
                      reduction = "harmony", 
                      assay = "SCT", 
                      dims = 1:20)
seuratdata <- FindNeighbors(object = seuratdata, 
                            dims=1:20,
                            reduction = "harmony")
seuratdata <- FindClusters(seuratdata, 
                           resolution = 1.2)
DimPlot(seuratdata,label.size = 4,repel = T,label = T)
save(seuratdata, file="seuratdata.RData")