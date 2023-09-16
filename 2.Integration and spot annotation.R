# 1.Loading packages

library(ggsci)
library(viridis)
library(ggplot2)
library(RCurl)
library(cowplot)
library(dplyr)
library(Seurat)

# 2.loadng ST-seq datasets

ST <- readRDS("integrated_data.rds")
ST # 54635 features across 18673 samples within 3 assays
DimPlot(ST,label = T,label.size = 5)
DimPlot(ST,label = T,label.size = 5,split.by = "type")

# 3.Find markers of each cluster  

all_DEGS <- FindAllMarkers(ST,only.pos = T,logfc.threshold = 0.5,min.pct = 0.5)
write.csv(all_DEGS,"all_DEGS.csv")

# 4.Marker expression level of each cluster

DoHeatmap(ST, slot = "data",features = c(
  "Nphs1","Nphs2","Podxl","Kdr","Eng", # Glom
  "Slc17a3","Slc22a8","Fxyd2", # proximal tubule(PT)
  "Slc5a2","Slc5a12", # Proximal convoluted tubule cell(PCT)
  "Atp11a", "Slc13a3", # proximal straight tubules(PST)
  "Slc27a2","Lrp2","Slc22a8","Slc5a2","Slc5a12","Fxyd2","Slc17a3", # proximal tubule(PT)
  "Atp11a", "Slc13a3","Slc34a1","Gpx3",# Proximal convoluted tubule cell(PTC)
  "Aqp1","Bst1", # Descending loop of Henle (DLH)
  "Slc12a1","Umod","Cldn8","Krt18","Krt8",  # Ascending loop of Henle(ALH)
  "Slc12a3"  # Distal convoluted tubule(DCT)
  "Atp6v0d2","Atp6v1g3","Slc4a1","Aqp6","Slc26a4","Hmx2", # Collecting duct intercalated cell(CD-IC)
  "Aqp2","Hsd11b2" #  Collecting duct principal / epithelial cell (CD-PC)
  "Rhbg","Insrr","Stmn1" # Collecting duct transitional cell (CD-TC)
  "Cdca3","Mki67" # Novel cell
  "Kdr","Ehd3","Plat","Vim","S100a4","Aqp1","Bst1","Pecam1","Eng","Cd34"  # Endo # Endothelial(Endo)
  "Nphs1", "Nphs2", # Podocyte (Podo)
  "Vim","S100a4", # Pericytes and vascular smooth muscle (Peri)
  "Plac8", # Fibroblast (Fibro)
  "C1qa","C1qb" # Macrophage (Macro)
  "Cd79a", "Cd79b" # B lymphocyte (B lymph)
  "Cxcr6","Ltb","Il7r","Cd3d","Cd3e","Ifng"  # T lymphocyte (T lymph)
  "Gzma","Nkg7","Gnly", # Natural killer cell (NK)
  "Lyz2","Cd14" # Monocytes (Mono)
  "S100a8","S100a9" # Neutrophil (Neutro)
  "Col1a1","Col1a2","Tagln","Acta2","C3","Vim","Myl9","S100a9","Thbs1" # Inter
  "Slc14a1","Slc14a2","Ly6d","Muc20","Psca","Upk1b","Upk3b" # Uro
)) + NoLegend()

# 5. Rename each cluster

ST2 <- ST
new.cluster.ids <- c("PT-S2","PT-S3","PT-injured","Immune","CD-PC","PT-S1",
                     "ALH","DLH","CD-PC","CD-IC","ALH","ALH","DLH","DCT","DLH",
                     "PT-S1","CD-PC","Glom","Uro","Adipo","Inter")
table(new.cluster.ids)
names(new.cluster.ids) <- levels(ST2)
ST2 <- RenameIdents(ST2, new.cluster.ids)

ST2$spot_type <- ST2@active.ident
ST2$spot_type <- factor(ST2$spot_type,levels = c("Glom","PT-S1","PT-S2","PT-S3","PT-injured",
                                                 "DLH","ALH","DCT","CD-IC","CD-PC","Immune",
                                                 "Uro","Inter","Adipo"))
ST2$type <- factor(ST2$type,levels = c("Control","AAN_2W","AAN_4W"))

# 6. visualization of cluster annotation 

cols = c("Glom" = '#E63863',
         "PT-S1" = '#E4C755',
         "PT-S2" = '#E59CC4',
         "PT-S3" = '#AB3282',
         "PT-injured" = '#E95C59',
         "DLH" = '#53A85F',
         "ALH" = '#F1BB72',
         "DCT" = '#F3B1A0',
         "CD-IC" = '#D6E7A3',
         "CD-PC" = '#57C3F3',
         "Immune" = "#00BFC4",
         "Inter" = '#8C549C',
         "Uro" = '#58A4C3',
         "Adipo" = "#23452F")

DimPlot(ST2,pt.size = 0.5,label = T,label.size = 4,group.by = "spot_type",cols = cols)  # 6*5
DimPlot(ST2,label = T,label.size = 0,pt.size = 0.5,split.by = "type", cols = cols) # 10*4

# 7. spatial visualization of cluster annotation 

ST2@images$Control_1
ST2@images$Control_2.1
ST2@images$AAN_2W_1.2
ST2@images$AAN_2W_2.3
ST2@images$AAN_4W_1.4
ST2@images$AAN_4W_2.5

(SpatialDimPlot(ST2,images = c("Control_1","AAN_2W_1.2","AAN_4W_1.4"),ncol = 3,cols = cols,alpha = 1) + NoLegend())/
(SpatialDimPlot(ST2,images = c("Control_2.1","AAN_2W_2.3","AAN_4W_2.5"),ncol = 3,cols = cols,alpha = 1) + NoLegend()) # 12*8

# 8.Marker expression level of each spot type

DotPlot(ST2,group.by = "spot_type",features = rev(c(
  "Nphs1","Nphs2", # Glom
  "Gpx3","Slc27a2",#"Slc34a1", # PT
  "Slc17a3","Slc22a8", # proximal tubule(PT)
  "Slc5a2","Slc5a12", # PT-S1
  "Slc13a3","Slc22a6",# "Slc17a3", # PT-S2
  "Slc22a7","Slc22a13", # "Slc7a13","Bcat1" # PT-S3
  "Havcr1","Vim", # PT-injured
  "Aqp1","Bst1", # Descending loop of Henle (DLH)
  "Slc12a1","Umod",#"Cldn8","Krt18","Krt8",  # Ascending loop of Henle(ALH)
  "Slc12a3",  # Distal convoluted tubule(DCT)
  "Atp6v0d2","Atp6v1g3",#"Slc4a1","Aqp6","Slc26a4","Hmx2", # Collecting duct intercalated cell(CD-IC)
  "Aqp2","Hsd11b2", #  Collecting duct principal / epithelial cell (CD-PC)
  "C1qa","C1qb", # Macro (Immune)
  "Ltb","Cd3e", # T lymph (Immune)
  "Ly6d","Upk1b", # Urothelium (Uro)
  "Col1a1","Col1a2", # Interstitium (Inter)
  "Adipoq","Lpl" # Adipocytes (Adipo)
    ))) + scale_color_gradientn(colours = viridis::viridis(20), guide = guide_colorbar(ticks.colour = "black",
    frame.colour = "black"), name = "Average \n expression") + 
  coord_flip() + xlab(NULL) + ylab(NULL) + 
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) # 5*8

# 9.Regional markers expression

SpatialFeaturePlot(ST,c("Lrp2","Aqp1","Aqp2"),
                   alpha = 1,
                   images = c("Control_1","AAN_2W_1.2","AAN_4W_1.4"),ncol = 3) 

# 10. Save annotated ST-seq dataset

saveRDS(ST2,"ST_anno.rds")

# 11.loadng ST-seq and scRNA-seq datasets

ST = readRDS("ST_anno.rds")
table(ST$spot_type)

scRNA = readRDS("scRNA-seq.anno.rds")
table(scRNA$cell_type)

# 12.Transfer Anchor analysis

anchors <- FindTransferAnchors(reference = scRNA, query = ST, 
                               normalization.method = "SCT",dims = 1:30, 
                               reference.reduction = "pca") # Found 1061 anchors

ST[["predictions"]] <- predictions.assay
DefaultAssay(ST) <- "predictions"

pdf("anchor.pdf",height = 8,width = 14)
FeaturePlot(ST,c("PT","ALH/DLH","DCT","CD-IC","CD-PC","Podo"),reduction = "umap",pt.size = 0.5,ncol = 3,cols = c("grey","red"))
dev.off()

# 13.Add module scores to the ST-seq dataset

all_DEGS <- FindAllMarkers(scRNA,only.pos = T,logfc.threshold = 0.5,min.pct = 0.5)
write.csv(all_DEGS,"all_DEGS.csv")

PT_genes <- list(all_DEGS[all_DEGS$cluster %in% "PT",]$gene)
DLH_genes <- list(all_DEGS[all_DEGS$cluster %in% "DLH",]$gene)
ALH_genes <- list(all_DEGS[all_DEGS$cluster %in% "ALH",]$gene)
DCT_genes <- list(all_DEGS[all_DEGS$cluster %in% "DCT",]$gene)
CD_IC_genes <- list(all_DEGS[all_DEGS$cluster %in% "CD-IC",]$gene)
CD_PC_genes <- list(all_DEGS[all_DEGS$cluster %in% "CD-PC",]$gene)
Endo_genes <- list(all_DEGS[all_DEGS$cluster %in% "Endo",]$gene)
Podo_genes <- list(all_DEGS[all_DEGS$cluster %in% "Podo",]$gene)
Fibro_genes <- list(all_DEGS[all_DEGS$cluster %in% "Fibro",]$gene)
Macro_genes <- list(all_DEGS[all_DEGS$cluster %in% "Macro",]$gene)
T_lymph_NK_genes <- list(all_DEGS[all_DEGS$cluster %in% "T lymph/NK",]$gene)
Neutro_genes <- list(all_DEGS[all_DEGS$cluster %in% "Neutro",]$gene)

ST2 <- AddModuleScore(object = ST,features = c(PT_genes,DLH_genes,ALH_genes,DCT_genes,CD_IC_genes,CD_PC_genes,
                                               Endo_genes,Podo_genes,Fibro_genes,Macro_genes,T_lymph_NK_genes,Neutro_genes),assay="SCT",
                      name=c("PT","DLH","ALH","DCT","CD-IC","CD-PC","Endo","Podo","Fibro","Macro","T lymph/NK","Neutro"))
colnames(head(ST2@meta.data))

meta_data <- read.csv("3.Integation/meta_data.csv",row.names = 1)
ST2 <- ST
ST2@meta.data <- meta_data
ST2$spot_type
AverageExpression_value <-  AverageExpression(ST2,# assays = "SCT", 
                                              features = c("PT1","DLH2","ALH3","DCT4","CD.IC5","CD.PC6","Endo7","Podo8","Fibro9","Macro10","T.lymph.NK11","Neutro12"),
                                                     return.seurat = FALSE, group.by = c("spot_type"))

factor(ST2$spot_type) 
ST2$spot_type = factor(ST2$spot_type,levels= c("Glom","PT-S1","PT-S2","PT-S3","PT-injured","DLH",
                                                   "ALH","DCT","CD-IC","CD-PC","Immune","Inter","Uro","Adipo"))
                       
                       
DotPlot(ST2,group.by = "spot_type",features = rev(c("Endo7","Podo8","PT1","DLH2","ALH3","DCT4",
                                                "CD.IC5","CD.PC6","Macro10",
                                                "T.lymph.NK11","Neutro12","Fibro9")))  +
  coord_flip() + xlab(NULL) + ylab(NULL) + scale_color_gradientn(colours = viridis(20), 
  guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"), name = "Module socre") +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) # 5*8

# 14. Spot numbers of each sample and spot type

data <- as.data.frame(table(ST$spot_type))
data <- data[order(data$Freq,decreasing=TRUE),] 
data$Var1 <- factor(data$Var1,levels=unique(data$Var1)) 
data
p1 <- ggplot(data,aes(x=Var1,y=Freq)) + geom_bar(aes(fill=Var1),stat="identity")+
  geom_text(aes(label=Freq),vjust=-1,size=3) + theme_test()+
  labs(x="",y="Cell number of spot types") + guides(fill="none")+
  theme(text=element_text(size=16),axis.text.x=element_text(angle=45,hjust=1)) + 
  scale_fill_manual(values = cols)
p1

data2 <- as.data.frame(table(ST$slice))
p2 <- ggplot(data2,aes(x=Var1,y=Freq)) + geom_bar(aes(fill=Var1),stat="identity")+
  geom_text(aes(label=Freq),vjust=-1,size=3) + theme_test()+
  labs(x="",y="Cell number of sample") + guides(fill="none")+
  theme(text=element_text(size=16),axis.text.x=element_text(angle=45,hjust=1)) 
p2

# 15. Spot proportion change

ST <- readRDS("ST_anno.rds")
table(ST$sample)
table(ST$spot_type)

ST$slice <- factor(ST$sample,levels = c("Control_1","Control_2","AAN_2W_1","AAN_2W_2","AAN_4W_1","AAN_4W_2"))
ST$spot_type <- factor(ST$spot_type,levels = rev(levels(ST$spot_type)))

data <- data.frame(table(ST$spot_type,ST$slice))
head(data)

data$proportion <- c(data[1:14,"Freq"]/2701,data[15:28,"Freq"]/3225,
                     data[29:42,"Freq"]/3247,data[43:56,"Freq"]/2832,
                     data[57:70,"Freq"]/3446,data[71:84,"Freq"]/3222) * 100
data

colnames(data) <- c("Spot_type","type","Freq","proportion")

ggplot(data,aes(x = Spot_type,y = type))+
  geom_point(aes(size=proportion,color=Spot_type),alpha=1)+
  scale_size(range=c(2,10))+
  scale_color_manual(values = cols) + xlab("") + ylab("") + 
  theme_bw() + coord_flip()  +   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
