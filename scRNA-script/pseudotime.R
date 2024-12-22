library(monocle)
library(Seurat)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(ggsci)
library(argparse)
DefaultAssay(seuratdata)<- 'RNA'
DimPlot(seuratdata,group.by ='orig.ident')
pd <- new('AnnotatedDataFrame', data = seuratdata@meta.data)
fData <- data.frame(gene_short_name = row.names(seuratdata), row.names = row.names(seuratdata))
fd <- new('AnnotatedDataFrame', data = fData)
gbm_cds <- newCellDataSet(as(as.matrix(seuratdata@assays[["RNA"]]@layers[["counts"]]),'sparseMatrix'),
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
gbm_cds <- estimateSizeFactors(gbm_cds)
gbm_cds <- estimateDispersions(gbm_cds)
disp_table <- dispersionTable(gbm_cds)
ordering_genes_temp <- subset(disp_table, mean_expression >= 0.1  & dispersion_empirical >= 1 * dispersion_fit) 
ordering_genes<-ordering_genes_temp$gene_id
gbm_cds <- setOrderingFilter(gbm_cds, ordering_genes)
print(plot_ordering_genes(gbm_cds))
tryCatch({
  gbm_cds <- reduceDimension(gbm_cds, max_components = 2,method = 'DDRTree')
}, error = function(e) {
  print("Error in reduceDimension(). Try to apply auto_param_selection = F")
  gbm_cds <- reduceDimension(gbm_cds, max_components = 2,method = 'DDRTree',auto_param_selection = F)
})
gbm_cds <- orderCells(gbm_cds)
df <- gbm_cds@phenoData@data