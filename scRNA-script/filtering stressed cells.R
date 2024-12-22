require(gruffi)
library(Seurat)
library(ggplot2)
#gruffi
table(seuratdata@meta.data[["seurat_clusters"]])
DimPlot(seuratdata,label.size = 4,repel = T,label = T)
ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
combined.obj <- Seurat.utils::SetupReductionsNtoKdimensions(obj = seuratdata,
                                                            nPCs = 20, 
                                                            dimensions = 3:2,
                                                            reduction_input = "pca",
                                                            reduction_output = "umap")
go1 <- "GO:0006096" # Glycolysis
go2 <- "GO:0034976" # ER-stress
go3 <- "GO:0042063" # Gliogenesis, negative filtering
combined.obj <- AutoFindGranuleResolution(obj = combined.obj)
(granule.res.4.gruffi <- GetGruffiClusteringName(combined.obj)) # Recalled from @misc$gruffi$'optimal.granule.res.
combined.obj <- ReassignSmallClusters(combined.obj, ident = granule.res.4.gruffi) # Will be stored in meta data column as "seurat_clusters.reassigned".
(granule.res.4.gruffi <- GetGruffiClusteringName(combined.obj))
#Glycolytic process GO:0006096
combined.obj <- AssignGranuleAverageScoresFromGOterm(obj = combined.obj, GO_term = go1, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)
#ER stress GO:0034976
combined.obj <- AssignGranuleAverageScoresFromGOterm(obj = combined.obj, GO_term = go2, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)
#Gliogenesis GO:0042063
combined.obj <- AssignGranuleAverageScoresFromGOterm(obj = combined.obj, GO_term = go3, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)
# Create score names:
(i1 <- ParseGruffiGranuleScoreName(goID = go1))
(i2 <- ParseGruffiGranuleScoreName(goID = go2))
(i3 <- ParseGruffiGranuleScoreName(goID = go3))
# Call Shiny app
combined.obj <- FindThresholdsShiny(obj = combined.obj,
                                    stress.ident1 = i1,
                                    stress.ident2 = i2,
                                    notstress.ident3 = i3)
"Dont forget to click the button in the app: Save New Thresholds"
StressUMAP(combined.obj)
#metadata
metadata_combined <- combined.obj@meta.data
metadata <- seuratdata@meta.data
metadata$is.Stressed <- metadata_combined$is.Stressed
seuratdata@meta.data <- metadata
DimPlot(seuratdata,label.size = 4,repel = T,label = T,group.by="is.Stressed")
metadata <- seuratdata@meta.data
stressed_data <- metadata[metadata$is.Stressed == "TRUE", ]
cluster_counts <- stressed_data %>%
  group_by(seurat_clusters) %>%
  tally()
#DDIT3-high cluster
p<-FeaturePlot(seuratdata, features = "DDIT3", min.cutoff = 0)
pdef<-p[[1]][["data"]]
metadata<-seuratdata@meta.data
pdef <- filter(pdef, DDIT3 > 0)
metadata$DDIT3 <- ifelse(rownames(metadata) %in% rownames(pdef), pdef[rownames(metadata), "DDIT3"], 0)
VlnPlot(seuratdata, features = "DDIT3", pt.size=0.1, group.by = "RNA_snn_res.1.2") 
#delete DDIT3-high cluster
seuratdata <- subset(seuratdata, subset = RNA_snn_res.1.2 == "0"|RNA_snn_res.1.2 == "1"|
                       RNA_snn_res.1.2 == "3"|RNA_snn_res.1.2 == "5"|
                       RNA_snn_res.1.2 == "7"|RNA_snn_res.1.2 == "8"|
                       RNA_snn_res.1.2 == "9"|RNA_snn_res.1.2 == "10"|
                       RNA_snn_res.1.2 == "11"|RNA_snn_res.1.2 == "12"|
                       RNA_snn_res.1.2 == "14"|RNA_snn_res.1.2 == "17")
#delete gruffi-labeled stressed cells
seuratdata <- subset(seuratdata, subset = is.Stressed==FALSE )
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
                           resolution = 0.6)
save(seuratdata, file="seuratdata.RData")