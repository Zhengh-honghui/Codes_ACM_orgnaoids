library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)

#UMAP
seuratdata_zhh1 <- subset(seuratdata,subset = orig.ident=="zhh1")
seuratdata_zhh2 <- subset(seuratdata,subset = orig.ident=="zhh2")
Idents(seuratdata) <- "celltype"
Idents(seuratdata) <- factor(Idents(seuratdata), levels =  c("Dividing Cells",
                                                             "Radial Glial Cells",
                                                             "Astroglia",
                                                             "Stressed Cells",
                                                             "Intermediate Progenitor",
                                                             "Deep Layer Neurons",
                                                             "Low-quality Cells",
                                                             "Upper Layer Neurons",
                                                             "Neurons") )
DimPlot(seuratdata,label.size = 4,repel = T,label = T)
dim(seuratdata)
UMAP <- as.data.frame(seuratdata@reductions$umap@cell.embeddings)
celltype <- Idents(seuratdata)
UMAP <- cbind(UMAP,celltype)
mytheme <- theme_void() + 
  theme(plot.margin = margin(5.5,15,5.5,5.5)) 
mycolors <- c("#9C755F","#3969AC","#B07AA1","#59A14F","#66C5CC","#FF9DA7","#F2B701")
p <- ggplot(data = UMAP, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = celltype),
             size = 0.4,
             alpha = 0.8)+
  theme_dr(xlength = 0.2, 
           ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.1, "inches"), 
                               ends = 'last', type = "closed")) + 息
  theme(panel.grid = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 5)))+
    scale_color_manual(values = mycolors) +
    scale_fill_manual(values = mycolors)
ggsave("UMAP.png", plot = p3, width = 8, height = 6, units = "in")

#DotPlot
Idents(seuratdata) <- "celltype"
Idents(seuratdata) <- factor(Idents(seuratdata), levels =  c("Dividing Cells",
                                                               "Radial Glial Cells",
                                                               "Astroglia",
                                                               "Intermediate Progenitor",
                                                               "Deep Layer Neurons",
                                                               "Upper Layer Neurons",
                                                               "Neurons") )
DimPlot(seuratdata)
DefaultAssay(seuratdata) <- "RNA"
Dividing_Cells_genes <- c("TOP2A", "UBE2C")
Radial_Glial_Cells_genes <- c("PAX6" ,"HOPX", "SOX2")
Astroglia_genes <- c("AQP4" ,"S100B", "APOE")
Intermediate_Progenitor_genes <- c("EOMES", "HES6" )
Deep_Layer_Neurons_genes <- c("BHLHE22","SLA","BCL11B", "TBR1")
Upper_Layer_Neurons_genes <- c("FGF12", "MEF2C", "RORB", "SATB2")
Neurons_genes <- c("GRIA2", "NEUROD2", "SNAP25")
levels(seuratdata) <- c("Dividing Cells",
                         "Radial Glial Cells",
                         "Astroglia",
                         "Intermediate Progenitor",
                         "Deep Layer Neurons",
                         "Upper Layer Neurons",
                         "Neurons")
features <- list("DC" = Dividing_Cells_genes,
                 "RG" = Radial_Glial_Cells_genes,
                 "AS" = Astroglia_genes,
                 "IPC" = Intermediate_Progenitor_genes,
                 "DLN" = Deep_Layer_Neurons_genes,
                 "ULN" = Upper_Layer_Neurons_genes,
                 "N" = Neurons_genes)
p <- DotPlot(object = seuratdata, features=features)
p1 <- ggplot(p$data, aes(x = features.plot, y = id)) + 
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) + 
  facet_grid(facets = ~feature.groups,  switch = "x", scales = "free_x", space = "free_x") +  
  scale_radius(breaks = c(25, 50, 75, 100), range = c(0,6)) + 
  theme_classic() + 
  scale_color_gradient2(low = "#00BEC4", mid = "white", high = "#F8766C")+
  theme(axis.text.x = element_text(angle = 90, 
                                   face = 1, size = 12, 
                                   family = "Arial Narrow", 
                                   hjust = 1, vjust = 0.5, 
                                   color = "black"), 
        axis.text.y = element_text(size = 12, face = 1, 
                                   family = "Arial Narrow", 
                                   color = "black"),
        legend.text = element_text(size = 8, face = 1, family = "Arial Narrow"), 
        legend.title = element_text(size = 10, face = 1, family = "Arial Narrow"),
        legend.position = 'top' , 
        strip.placement = "outside", 
        strip.text.x = element_text(size = 12, family = "Arial Narrow Bold"),
        axis.title = element_blank())+
  guides(colour = guide_colourbar(title.vjust = 0.9, title.hjust = 0))+
  labs(size = "Percent Expressed", color = "Average Expression")
ggsave("dotplot.png", plot = p1, width = 14, height = 6, units = "in")

#VlnPlot
dir.create("Vlnplot", showWarnings = FALSE)
for (gene in markerlist) {
  p<-VlnPlot(seuratdata1, features = gene, cols=mycolors,pt.size=0.1,
             group.by = "celltype") 
  p1<-p+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold")) 
  ggsave(filename = paste0("Vlnplot/", gene, ".png"), plot = p1, width = 8, height = 6, units = "in")
}

#FeaturePlot
dir.create("GeneUMAP", showWarnings = FALSE)
mycolors<-colorRampPalette(c("lightblue","white","darkred"))(100)
markerlist <- unique(markerlist)
genes_in_seurat <- rownames(seuratdata)
markerlist <- markerlist[markerlist %in% genes_in_seurat]
for (gene in markerlist) {
  p <- FeaturePlot(seuratdata, features = markerlist3,cols=mycolors)
  p1 <- p +
    theme_dr(xlength = 0.2, #x轴长度
             ylength = 0.2, #y轴长度
             arrow = grid::arrow(length = unit(0.1, "inches"), 
                                 ends = 'last', type = "closed")) + 
    theme(panel.grid = element_blank())
  ggsave(filename = paste0("GeneUMAP/", gene, ".png"), plot = p1, width = 8, height = 6, units = "in")
}

#Proportion
table(seuratdata@meta.data$celltype)
table(seuratdata_zhh1@meta.data$celltype)
table(seuratdata_zhh2@meta.data$celltype)
sample <- c(rep("CTR", 9), rep("ACM", 9))
celltype <- rep(c("Stressed Cells","Low-quality Cells",
                  "Dividing Cells", "Radial Glial Cells", 
                  "Astroglia", "Intermediate Progenitor", 
                  "Neurons","Upper Layer Neurons","Deep Layer Neurons"), 2)
number <- c(1899,749,271,1037,41,148,1847,484,1685,
            1357,736,393,1191,72,334,1683,487,2119)
df1 <- data.frame(sample, celltype, number)
colors2 <- c("#E15759","#11A579","#9C755F","#3969AC","#B07AA1","#59A14F","#F2B701","#FF9DA7","#66C5CC")
df1$celltype <- factor(df1$celltype, levels = c("Stressed Cells","Low-quality Cells",
                                                "Dividing Cells", "Radial Glial Cells", 
                                                "Astroglia", "Intermediate Progenitor", 
                                                "Neurons","Upper Layer Neurons","Deep Layer Neurons"))
p <-ggplot(df1, aes(x = sample, y = number, fill = celltype))+ 
  geom_bar(stat="identity",position="fill", width=0.7, linewidth=0.2)+
  scale_fill_manual(values = colors2)+
  ggprism::theme_prism(border = T)+
  labs(y="Proportion",x="")+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold")) 
ggsave("Proportion.png", plot = p, width = 7, height = 6, units = "in")

#CTR-ACM-DotPlot
metadata<-seuratdata@meta.data
metadata<-metadata %>%
  mutate(orig.ident = case_when(
    orig.ident == "zhh1" ~ "CTR",
    TRUE ~ "ACM"
  ))
seuratdata@meta.data <- metadata
markerlist1 <- c('XBP1','TXNRD1','SOD1','ATG5',
                 'BAX','CASP8','HIF1A')
markerlist2 <- c('SLC17A7', 'SLC17A6', 'SLC1A6', 
                 'SLC1A3', 'SLC1A2', 'SLC32A1', 
                 'SLC6A11', 'GAD1', 'GAD2')
markerlist3 <- c("TBR1","BCL11B","SATB2","CUX1",
                 "CUX2","POU3F2","POU3F3","SOX2","TUBB3","DCX","MAP2")
markerlist4 <- c('PLIN2', 'APOE', 'LDLR', 'ABCA1', 'SOAT1', 'HMGCR')
p <-DotPlot(seuratdata, 
            features = "RELN",
            cols =c("#f0dd93", "#bd4335" ),
            assay='RNA',
            group.by = "celltype")+ theme_bw()
p1 <- ggplot(p$data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_radius(breaks = c(25, 50, 75, 100), range = c(0,6)) +
  theme_classic() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme(
    axis.text.x = element_text(size = 10, family = "Arial Narrow", color = "black"),
    axis.text.y = element_text(size = 10, face = 1, family = "Arial Narrow", color = "black"),
    legend.text = element_text(size = 8, face = 1, family = "Arial Narrow"),
    legend.title = element_text(size = 10, face = 1, family = "Arial Narrow"),
    legend.position = 'top',
    strip.placement = "outside",
    strip.text.x = element_text(size = 12, family = "Arial Narrow Bold"),
    axis.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1), # 边框大小
    plot.margin = margin(10, 10, 10, 10, "pt") 
  ) +
  guides(colour = guide_colourbar(title.vjust = 0.9, title.hjust = 90)) +
  labs(size = "Percent Expressed", color = "Average Expression")
ggsave("dotplot_markerlist1_celltype.png", plot = p1, width = 6.2, height = 3.7, units = "in")

#Pseudotime
library(ggforce)
library(ggrastr)
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
df <- gbm_cds@phenoData@data
data_df <- t(reducedDimS(gbm_cds)) %>% as.data.frame() %>% 
  select_(Component_1 = 1, Component_2 = 2) %>% 
  rownames_to_column("cells") 
data_df <- data_df %>% mutate(celltype = df$celltype)
df<-df %>%
  mutate(orig.ident = case_when(
    orig.ident == "zhh1" ~ "CTR",
    TRUE ~ "ACM"
  ))
data_df$orig.ident<-df$orig.ident
df_zhh1 <- data_df %>% filter(orig.ident == "zhh1")
df_zhh2 <- data_df %>% filter(orig.ident == "zhh2")
celltype_order <- c("Dividing Cells",
                    "Radial Glial Cells",
                    "Astroglia",
                    "Intermediate Progenitor",
                    "Deep Layer Neurons",
                    "Upper Layer Neurons",
                    "Neurons")
data_df$celltype <- factor(data_df$celltype, levels = celltype_order)
#trajectory
dp_mst <- minSpanningTree(gbm_cds)
reduced_dim_coords <- reducedDimK(gbm_cds)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>%
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  mutate(sample_name = rownames(.), sample_state = rownames(.))
edge_df <- dp_mst %>% igraph::as_data_frame() %>%
  select_(source = "from", target = "to") %>%
  left_join(ica_space_df %>% select_(source = "sample_name",
                                     source_prin_graph_dim_1 = "prin_graph_dim_1",
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>%
  left_join(ica_space_df %>% select_(target = "sample_name",
                                     target_prin_graph_dim_1 = "prin_graph_dim_1",
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")
#plot
mycolors<-c("#9C755F","#3969AC","#B07AA1","#59A14F","#66C5CC","#FF9DA7","#F2B701")
g <- ggplot() +
  geom_point_rast(data = data_df, aes(x = Component_1,
                                      y = Component_2,
                                      color =celltype),
                  alpha = 0.5, size = 1) +
  geom_segment(aes_string(x = "source_prin_graph_dim_1",
                          y = "source_prin_graph_dim_2",
                          xend = "target_prin_graph_dim_1",
                          yend = "target_prin_graph_dim_2"),
               linewidth = 0.8,
               linetype = "solid", na.rm = TRUE, data = edge_df)+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+theme_void()+
  theme_dr(xlength = 0.2, 
           ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.1, "inches"), 
                               ends = 'last', type = "closed")) + 
  theme(panel.grid = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors)
ggsave("trajectory.png", plot = g, width = 8, height = 6, units = "in")
#density
mycolors <- c("#F8766C","#00BEC4")
pq<-ggplot(df, 
           aes(Pseudotime, 
               colour = orig.ident, 
               fill = orig.ident)) +   
  geom_density( 
    bw = 0.5,  
    size = 0.8, 
    alpha = 0.5) +  
  scale_color_manual(values = mycolors) + 
  scale_fill_manual(values = mycolors) +
  theme_classic2() 
ggsave("density.png", plot = pq, width = 6.4, height = 2, units = "in")

#GO
library(tidyverse)
library(ggh4x)
library(ggfun)
library(ggnewscale)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)
GO_df<- goresult
ids_to_select <- c("GO:0021543", "GO:0021987", "GO:0021537", "GO:0030900", 
                   "GO:0001764", "GO:0006695", "GO:0050769", "GO:0061448", 
                   "GO:0051962", "GO:0050767", "GO:0008203")
GO_df <- GO_df[GO_df$ID %in% ids_to_select, ]
plot <- GO_df %>%
  ggplot() + 
  geom_point(data = GO_df ,
             aes(x = FoldEnrichment, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = FoldEnrichment), shape = 21) + 
  scale_fill_gradient(low = "#0570b0", high ="#a6bddb", name = "p.adjust") + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:11,
                                   labels = GO_df$Description)) + 
  labs(x = "FoldEnrichment", y = "Description") + 
  scale_size(range = c(1,5),
             guide = guide_legend(override.aes = list(fill = "#000000"))) + 
  theme_bw() + 
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color="#0570b0"),
    ggh4x.axis.nesttext.y = element_text(colour ="#0570b0"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    axis.text.y = element_text(color = "#225ea8"),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(10, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    
  ) + 
  coord_cartesian(clip = "off")
plot
ggsave("Neurons_GO.png", plot = plot, width = 8.8, height = 4, units = "in")

#GSEA
library(GseaVis)
geneSetID = c("GO:0045211", "GO:0098984", "GO:0099003", "GO:0099572", "GO:0007409", "GO:0014069", "GO:0032279", "GO:0099634", "GO:0097485", "GO:0061564")
colors <-      c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74",
                 "#80BA5A","#E68310","#008695","#A5AA99","#66C5CC",
                 "#E15759","#B07AA1")
p1<-gseaNb(object = ego3,
           geneSetID = geneSetID,
           subPlot = 2,
           curveCol = colors)
ggsave("AS_GSEA.png", plot = p1, width = 7, height = 6.2, units = "in")
bubble_plot <- function(data1, data2,  color = "#08306b", font_size = 12, color_label = "data1_value", size_label = "data2_value") {
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(scales)
  validate_data <- function(d1, d2) {
    if (!identical(dim(d1), dim(d2))) {
      stop("两个数据集必须具有相同的维度。")
    }
    if (!identical(rownames(d1), rownames(d2)) || !identical(colnames(d1), colnames(d2))) {
      stop("两个数据集的行名和列名必须完全匹配。")
    }
    list(d1 = d1[rownames(d2), colnames(d2)], d2 = d2)
  }
  validated_data <- validate_data(data1, data2)
  data1 <- validated_data$d1
  data2 <- validated_data$d2
  original_row_order <- rownames(data1)
  original_col_order <- colnames(data1)
  x_axis = "sample"
  y_axis = "gene"
  prepare_data <- function(d1, d2) {
    d1_long <- d1 %>%
      rownames_to_column(var = x_axis) %>%
      pivot_longer(cols = -all_of(x_axis), names_to = y_axis, values_to = "data1")
    d2_long <- d2 %>%
      rownames_to_column(var = x_axis) %>%
      pivot_longer(cols = -all_of(x_axis), names_to = y_axis, values_to = "data2")
    data_combined <- d1_long %>%
      left_join(d2_long, by = c(x_axis, y_axis))
    if (any(is.na(data_combined))) {
      stop("合并数据时出错。请确保数据集一致。")
    }
    data_combined %>%
      mutate(
        data1_normalized = data1,
        data2_normalized = data2
      )
  }
  data_combined <- prepare_data(data1, data2)
  color_gradient <- colorRampPalette(c("grey90", color))(100)
  p <- ggplot(data_combined, aes_string(x = x_axis, y = y_axis)) +
    geom_point(aes(size = data2_normalized, color = data1_normalized), alpha = 0.8) +
    scale_size(range = c(1, 5), name = size_label) +
    scale_color_gradientn(colors = color_gradient, name = color_label) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = font_size),
      axis.text.y = element_text(size = font_size),
      axis.title.x = element_blank(),  # 移除x轴标签
      axis.title.y = element_blank(),  # 移除y轴标签
      legend.position = c(0.5, -0.9),  # 图例的中心位置
      legend.justification = c(0.5, 0),  # 图例相对于绘图区域的对齐方式
      legend.key.size = unit(0.6, "cm"),
      legend.text = element_text(size = font_size),
      panel.grid = element_line(linewidth = 0.3)
    ) +
    scale_x_discrete(limits = original_row_order) +
    scale_y_discrete(limits = rev(original_col_order))
  
  return(p)
}
gsearesult<-ego3@result
gsearesult <- gsearesult[gsearesult$ID %in% geneSetID, ]
NES_data <- gsearesult$NES
NES_data <- data.frame(t(NES_data))
setSize_data <- gsearesult$setSize
setSize_data <- data.frame(t(setSize_data))
qvalue_data <- data.frame(-log10(gsearesult$qvalue))
qvalue_data <- data.frame(t(qvalue_data))
rownames(qvalue_data) <- NULL  
rownames(qvalue_data) <- 1:nrow(qvalue_data)
p2<-bubble_plot(data1=qvalue_data, 
                data2=setSize_data, 
                color = "#00BEC4",
                #color = "#08306b",
                font_size = 12, 
                color_label="-log(Qvalue)", 
                size_label="setSize")+
  theme_void()+
  theme(  
    legend.position = "bottom",  # 将图例放在顶部
    legend.justification = "center")
ggsave("AS_GSEA_Q.png", plot = p2, width = 5.95, height = 2.6, units = "in")