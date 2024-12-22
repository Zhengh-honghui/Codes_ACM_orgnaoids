library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
#cell type marker
Dividing_Cells_genes <- c("TOP2A", "UBE2C")
Radial_Glial_Cells_genes <- c("PAX6" ,"HOPX", "SOX2")
Astroglia_genes <- c("AQP4" ,"S100B", "APOE")
Intermediate_Progenitor_genes <- c("EOMES", "HES6" )
Deep_Layer_Neurons_genes <- c("BHLHE22","SLA","BCL11B", "TBR1")
Upper_Layer_Neurons_genes <- c("FGF12", "MEF2C", "RORB", "SATB2")
Neurons_genes <- c("GRIA2", "NEUROD2", "SNAP25")
markerlist<-c("TOP2A", "UBE2C",
              "PAX6" ,"HOPX", "SOX2",
              "AQP4" ,"S100B", "APOE",
              "EOMES", "HES6",
              "BHLHE22","SLA","BCL11B", "TBR1",
              "FGF12", "MEF2C", "RORB", "SATB2",
              "GRIA2", "NEUROD2", "SNAP25")
DotPlot(seuratdata, 
        features = markerlist,
        cols =c("#f0dd93", "#bd4335" ),
        assay='RNA',
        group.by = "celltype")+ theme_bw()
VlnPlot(seuratdata, 
        features = markerlist, 
        pt.size=0.1, 
        group.by = "RNA_snn_res.0.6")
#annotation
celltype=data.frame(ClusterID=0:8,
                    celltype= 0:8) 
table(seuratdata@meta.data[["RNA_snn_res.0.6"]])
DimPlot(seuratdata,group.by="RNA_snn_res.0.6",label.size = 4,repel = T,label = T)
celltype[celltype$ClusterID %in% c(0),2]='Radial Glial Cells'
celltype[celltype$ClusterID %in% c(1),2]='Deep Layer Neurons'
celltype[celltype$ClusterID %in% c(2),2]='Neurons'
celltype[celltype$ClusterID %in% c(3),2]='Deep Layer Neurons'
celltype[celltype$ClusterID %in% c(4),2]='Neurons'
celltype[celltype$ClusterID %in% c(5),2]='Upper Layer Neurons'
celltype[celltype$ClusterID %in% c(6),2]='Dividing Cells'
celltype[celltype$ClusterID %in% c(7),2]='Intermediate Progenitor'
celltype[celltype$ClusterID %in% c(8),2]='Astroglia'
table(celltype$celltype)
seuratdata@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  seuratdata@meta.data[which(seuratdata@meta.data[["RNA_snn_res.0.6"]] == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
DimPlot(seuratdata,group.by="celltype",label.size = 4,repel = T,label = T)
save(seuratdata, file="seuratdata.RData")