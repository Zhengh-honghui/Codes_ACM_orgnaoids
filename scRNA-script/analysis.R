library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
library(MAST)
library(org.Hs.eg.db)
library(clusterProfiler)
#DEGs
type=c("Dividing Cells",
       "Radial Glial Cells",
       "Astroglia",
       "Intermediate Progenitor",
       "Deep Layer Neurons",
       "Upper Layer Neurons",
       "Neurons")
Idents(seuratdata)="celltype"
for (i in 1:length(type)) {
  deg=FindMarkers(seuratdata,
                  ident.1 = "zhh2",ident.2 = "zhh1",
                  group.by = "orig.ident",
                  logfc.threshold = 0.25,
                  test.use = "MAST",
                  subset.ident =type[i])
  deg$celltype=type[i]
  deg$group <- case_when(
    deg$avg_log2FC > 0 ~ "up",
    deg$avg_log2FC < 0 ~ "down"
  )
  write.csv(deg,file = paste0( type[i],'_foldchange_MAST.csv') )
}
#GO
deg<- read.csv('Dividing Cells_foldchange_MAST.csv')
deg <- subset(deg, p_val_adj < 0.05)
deg$gene=deg$X
deg_gene=deg$X
ids=bitr(deg_gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
deg=merge(deg,ids,by.x='gene',by.y='SYMBOL')
deg_gene <- subset(deg, celltype == type[1])$gene
enrich.go <- enrichGO(gene = deg_gene,  
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'SYMBOL',  
                      ont = 'ALL', 
                      pAdjustMethod = 'BH',  
                      minGSSize = 50, 
                      pvalueCutoff = 1,  
                      qvalueCutoff = 0.2, 
                      readable = FALSE)   
goresult<-enrich.go@result
write.csv(goresult,file = paste0(type[1],'_GOresult.csv'))
#KEGG
deg_gene <- deg$gene
entrez_id = mapIds(x = org.Hs.eg.db,
                   keys =  deg_gene,
                   keytype = "SYMBOL",
                   column = "ENTREZID")
erich.kegg <- enrichKEGG(gene = entrez_id,
                         organism = "hsa",
                         keyType = "kegg",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2)
erich.kegg@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "",erich.kegg@result$Description)
keggresult=erich.kegg@result
write.csv(keggresult,file = paste0(type[1],'_KEGGresult.csv'))
#GSEA
seuratdata_1 <- seuratdata
seuratdata_1 <- RenameIdents(seuratdata, 
                             `0` = "Radial Glial Cells", 
                             `1` = "Deep Layer Neurons", 
                             `2` = "Neurons", 
                             `3` = "Deep Layer Neurons", 
                             `4` = "Neurons", 
                             `5` = "Upper Layer Neurons", 
                             `6` = "Dividing Cells",
                             `7` = "Intermediate Progenitor",
                             `8` = "Astroglia")
Idents(seuratdata_1)
seuratdata_1$celltype.stim <- paste(Idents(seuratdata_1), seuratdata_1$sample, sep = "_")
seuratdata_1$celltype <- Idents(seuratdata_1)
Idents(seuratdata_1) <- "celltype.stim"
Intermediate_Progenitor_mast <- FoldChange(seuratdata_1, 
                                           ident.1 = "Intermediate Progenitor_zhh2", 
                                           ident.2 = "Intermediate Progenitor_zhh1", 
                                           verbose = FALSE, test.use='MAST')
write.table(Intermediate_Progenitor_mast, file = 'Intermediate Progenitor.txt', sep='\t')
gene_anno <- read.csv('annotations_ahb.txt')
ens <- rownames(seuratdata_1)%>% data.frame()
bg_anno <- inner_join(ens, gene_anno, by=c("."="gene_name")) 
gene_list <- read.table('Intermediate Progenitor.txt')
gene_list['gene_name'] <- rownames(gene_list)
bg_anno <- merge(gene_list,gene_anno, by=c("."="gene_name"))
order_DEGs_new <- bg_anno[order(-bg_anno$avg_log2FC),]
new_gene_id <- order_DEGs_new$avg_log2FC
names(new_gene_id) <- order_DEGs_new$gene_id
ego3 <- gseGO(geneList     = new_gene_id,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              minGSSize    = 50,
              maxGSSize    = 500,
              pvalueCutoff = 1,
              keyType = "ENSEMBL")
cluster_summary <- data.frame(ego3)
write.csv(cluster_summary, 'Intermediate Progenitor_GSEAresult.csv')