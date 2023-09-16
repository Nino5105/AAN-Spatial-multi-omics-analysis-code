# 1.Loading packages

library(Seurat)
library(ggrepel)
library(org.Mm.eg.db)
library(clusterProfiler)
library(UpSetR)

# 2.Loadng ST-seq dataset

ST <- readRDS("ST_anno.rds")
DimPlot(ST)

# 3.Subset each spot type from the ST-seq dataset

# PT-S1

PT_S1 <- subset(ST,subset = spot_type == "PT-S1")
PT_S1@active.ident <- PT_S1$type
PT_S1_DEG_2 <- FindMarkers(PT_S1,ident.1 = "AAN_2W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(PT_S1_DEG_2)
PT_S1_DEG_2$type = ifelse(PT_S1_DEG_2$p_val_adj < 0.05 & abs(PT_S1_DEG_2$avg_log2FC) >= 0.25,ifelse(PT_S1_DEG_2$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(PT_S1_DEG_2$type)
# Down Stable     Up 
# 273   2224    281
PT_S1_DEG_2_up <- rownames(PT_S1_DEG_2[PT_S1_DEG_2$type == "Up",])
PT_S1_DEG_2_down <- rownames(PT_S1_DEG_2[PT_S1_DEG_2$type == "Down",])

PT_S1_DEG_4 <- FindMarkers(PT_S1,ident.1 = "AAN_4W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(PT_S1_DEG_4)
PT_S1_DEG_4$type = ifelse(PT_S1_DEG_4$p_val_adj < 0.05 & abs(PT_S1_DEG_4$avg_log2FC) >= 0.25,ifelse(PT_S1_DEG_4$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(PT_S1_DEG_4$type)
# Down Stable     Up 
# 282   1547    931 
PT_S1_DEG_4_up <- rownames(PT_S1_DEG_4[PT_S1_DEG_4$type == "Up",])
PT_S1_DEG_4_down <- rownames(PT_S1_DEG_4[PT_S1_DEG_4$type == "Down",])

# PT-S2

PT_S2 <- subset(ST,subset = spot_type == "PT-S2")
PT_S2@active.ident <- PT_S2$type
PT_S2_DEG_2 <- FindMarkers(PT_S2,ident.1 = "AAN_2W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(PT_S2_DEG_2)
PT_S2_DEG_2$type = ifelse(PT_S2_DEG_2$p_val_adj < 0.05 & abs(PT_S2_DEG_2$avg_log2FC) >= 0.25,ifelse(PT_S2_DEG_2$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(PT_S2_DEG_2$type)
# Down Stable     Up 
# 196   2079    229
PT_S2_DEG_2_up <- rownames(PT_S2_DEG_2[PT_S2_DEG_2$type == "Up",])
PT_S2_DEG_2_down <- rownames(PT_S2_DEG_2[PT_S2_DEG_2$type == "Down",])

PT_S2_DEG_4 <- FindMarkers(PT_S2,ident.1 = "AAN_4W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(PT_S2_DEG_4)
PT_S2_DEG_4$type = ifelse(PT_S2_DEG_4$p_val_adj < 0.05 & abs(PT_S2_DEG_4$avg_log2FC) >= 0.25,ifelse(PT_S2_DEG_4$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(PT_S2_DEG_4$type)
# Down Stable     Up 
# 185   1745    636
PT_S2_DEG_4_up <- rownames(PT_S2_DEG_4[PT_S2_DEG_4$type == "Up",])
PT_S2_DEG_4_down <- rownames(PT_S2_DEG_4[PT_S2_DEG_4$type == "Down",])

# PT-S3

PT_S3 <- subset(ST,subset = spot_type == "PT-S3")
PT_S3@active.ident <- PT_S3$type
PT_S3_DEG_2 <- FindMarkers(PT_S3,ident.1 = "AAN_2W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(PT_S3_DEG_2)
PT_S3_DEG_2$type = ifelse(PT_S3_DEG_2$p_val_adj < 0.05 & abs(PT_S3_DEG_2$avg_log2FC) >= 0.25,ifelse(PT_S3_DEG_2$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(PT_S3_DEG_2$type)
# Down Stable     Up 
# 218   2144    264 
PT_S3_DEG_2_up <- rownames(PT_S3_DEG_2[PT_S3_DEG_2$type == "Up",])
PT_S3_DEG_2_down <- rownames(PT_S3_DEG_2[PT_S3_DEG_2$type == "Down",])

PT_S3_DEG_4 <- FindMarkers(PT_S3,ident.1 = "AAN_4W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(PT_S3_DEG_4)
PT_S3_DEG_4$type = ifelse(PT_S3_DEG_4$p_val_adj < 0.05 & abs(PT_S3_DEG_4$avg_log2FC) >= 0.25,ifelse(PT_S3_DEG_4$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(PT_S3_DEG_4$type)
# Down Stable     Up 
# 244   2031    380 
PT_S3_DEG_4_up <- rownames(PT_S3_DEG_4[PT_S3_DEG_4$type == "Up",])
PT_S3_DEG_4_down <- rownames(PT_S3_DEG_4[PT_S3_DEG_4$type == "Down",])

# PT_injured

PT_injured <- subset(ST,subset = spot_type == "PT-injured")
PT_injured@active.ident <- PT_injured$type
PT_injured_DEG_2 <- FindMarkers(PT_injured,ident.1 = "AAN_2W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(PT_injured_DEG_2)
PT_injured_DEG_2$type = ifelse(PT_injured_DEG_2$p_val_adj < 0.05 & abs(PT_injured_DEG_2$avg_log2FC) >= 0.25,ifelse(PT_injured_DEG_2$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(PT_injured_DEG_2$type)
# Down Stable     Up 
# 6   2727     20
PT_injured_DEG_2_up <- rownames(PT_injured_DEG_2[PT_injured_DEG_2$type == "Up",])
PT_injured_DEG_2_down <- rownames(PT_injured_DEG_2[PT_injured_DEG_2$type == "Down",])

PT_injured_DEG_4 <- FindMarkers(PT_injured,ident.1 = "AAN_4W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(PT_injured_DEG_4)
PT_injured_DEG_4$type = ifelse(PT_injured_DEG_4$p_val_adj < 0.05 & abs(PT_injured_DEG_4$avg_log2FC) >= 0.25,ifelse(PT_injured_DEG_4$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(PT_injured_DEG_4$type)
# Down Stable     Up 
# 25   2535    169
PT_injured_DEG_4_up <- rownames(PT_injured_DEG_4[PT_injured_DEG_4$type == "Up",])
PT_injured_DEG_4_down <- rownames(PT_injured_DEG_4[PT_injured_DEG_4$type == "Down",])

# DLH

DLH <- subset(ST,subset = spot_type == "DLH")
DLH@active.ident <- DLH$type
DLH_DEG_2 <- FindMarkers(DLH,ident.1 = "AAN_2W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(DLH_DEG_2)
DLH_DEG_2$type = ifelse(DLH_DEG_2$p_val_adj < 0.05 & abs(DLH_DEG_2$avg_log2FC) >= 0.25,ifelse(DLH_DEG_2$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(DLH_DEG_2$type)
# Down Stable     Up 
# 114   2181    130
DLH_DEG_2_up <- rownames(DLH_DEG_2[DLH_DEG_2$type == "Up",])
DLH_DEG_2_down <- rownames(DLH_DEG_2[DLH_DEG_2$type == "Down",])

DLH_DEG_4 <- FindMarkers(DLH,ident.1 = "AAN_4W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(DLH_DEG_4)
DLH_DEG_4$type = ifelse(DLH_DEG_4$p_val_adj < 0.05 & abs(DLH_DEG_4$avg_log2FC) >= 0.25,ifelse(DLH_DEG_4$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(DLH_DEG_4$type)
# Down Stable     Up 
# 245   1815    408
DLH_DEG_4_up <- rownames(DLH_DEG_4[DLH_DEG_4$type == "Up",])
DLH_DEG_4_down <- rownames(DLH_DEG_4[DLH_DEG_4$type == "Down",])

# ALH

ALH <- subset(ST,subset = spot_type == "ALH")
ALH@active.ident <- ALH$type
ALH_DEG_2 <- FindMarkers(ALH,ident.1 = "AAN_2W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(ALH_DEG_2)
ALH_DEG_2$type = ifelse(ALH_DEG_2$p_val_adj < 0.05 & abs(ALH_DEG_2$avg_log2FC) >= 0.25,ifelse(ALH_DEG_2$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(ALH_DEG_2$type)
# Down Stable     Up 
# 185   2405    254
ALH_DEG_2_up <- rownames(ALH_DEG_2[ALH_DEG_2$type == "Up",])
ALH_DEG_2_down <- rownames(ALH_DEG_2[ALH_DEG_2$type == "Down",])

ALH_DEG_4 <- FindMarkers(ALH,ident.1 = "AAN_4W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(ALH_DEG_4)
ALH_DEG_4$type = ifelse(ALH_DEG_4$p_val_adj < 0.05 & abs(ALH_DEG_4$avg_log2FC) >= 0.25,ifelse(ALH_DEG_4$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(ALH_DEG_4$type)
# Down Stable     Up 
# 313   1885    631
ALH_DEG_4_up <- rownames(ALH_DEG_4[ALH_DEG_4$type == "Up",])
ALH_DEG_4_down <- rownames(ALH_DEG_4[ALH_DEG_4$type == "Down",])

# DCT

DCT <- subset(ST,subset = spot_type == "DCT")
DCT@active.ident <- DCT$type
DCT_DEG_2 <- FindMarkers(DCT,ident.1 = "AAN_2W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(DCT_DEG_2)
DCT_DEG_2$type = ifelse(DCT_DEG_2$p_val_adj < 0.05 & abs(DCT_DEG_2$avg_log2FC) >= 0.25,ifelse(DCT_DEG_2$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(DCT_DEG_2$type)
# Down Stable     Up 
# 89   2467    129 
DCT_DEG_2_up <- rownames(DCT_DEG_2[DCT_DEG_2$type == "Up",])
DCT_DEG_2_down <- rownames(DCT_DEG_2[DCT_DEG_2$type == "Down",])

DCT_DEG_4 <- FindMarkers(DCT,ident.1 = "AAN_4W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(DCT_DEG_4)
DCT_DEG_4$type = ifelse(DCT_DEG_4$p_val_adj < 0.05 & abs(DCT_DEG_4$avg_log2FC) >= 0.25,ifelse(DCT_DEG_4$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(DCT_DEG_4$type)
# Down Stable     Up 
# 125   2068    470
DCT_DEG_4_up <- rownames(DCT_DEG_4[DCT_DEG_4$type == "Up",])
DCT_DEG_4_down <- rownames(DCT_DEG_4[DCT_DEG_4$type == "Down",])

# CD_IC

CD_IC <- subset(ST,subset = spot_type == "CD-IC")
CD_IC@active.ident <- CD_IC$type
CD_IC_DEG_2 <- FindMarkers(CD_IC,ident.1 = "AAN_2W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(CD_IC_DEG_2)
CD_IC_DEG_2$type = ifelse(CD_IC_DEG_2$p_val_adj < 0.05 & abs(CD_IC_DEG_2$avg_log2FC) >= 0.25,ifelse(CD_IC_DEG_2$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD_IC_DEG_2$type)
# Down Stable     Up 
# 106   2483    134
CD_IC_DEG_2_up <- rownames(CD_IC_DEG_2[CD_IC_DEG_2$type == "Up",])
CD_IC_DEG_2_down <- rownames(CD_IC_DEG_2[CD_IC_DEG_2$type == "Down",])

CD_IC_DEG_4 <- FindMarkers(CD_IC,ident.1 = "AAN_4W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(CD_IC_DEG_4)
CD_IC_DEG_4$type = ifelse(CD_IC_DEG_4$p_val_adj < 0.05 & abs(CD_IC_DEG_4$avg_log2FC) >= 0.25,ifelse(CD_IC_DEG_4$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD_IC_DEG_4$type)
# Down Stable     Up 
# 180   2152    385
CD_IC_DEG_4_up <- rownames(CD_IC_DEG_4[CD_IC_DEG_4$type == "Up",])
CD_IC_DEG_4_down <- rownames(CD_IC_DEG_4[CD_IC_DEG_4$type == "Down",])

# CD_PC

CD_PC <- subset(ST,subset = spot_type == "CD-PC")
CD_PC@active.ident <- CD_PC$type
CD_PC_DEG_2 <- FindMarkers(CD_PC,ident.1 = "AAN_2W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(CD_PC_DEG_2)
CD_PC_DEG_2$type = ifelse(CD_PC_DEG_2$p_val_adj < 0.05 & abs(CD_PC_DEG_2$avg_log2FC) >= 0.25,ifelse(CD_PC_DEG_2$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD_PC_DEG_2$type)
# Down Stable     Up 
# 152   2040    146
CD_PC_DEG_2_up <- rownames(CD_PC_DEG_2[CD_PC_DEG_2$type == "Up",])
CD_PC_DEG_2_down <- rownames(CD_PC_DEG_2[CD_PC_DEG_2$type == "Down",])

CD_PC_DEG_4 <- FindMarkers(CD_PC,ident.1 = "AAN_4W", ident.2 = "Control",logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
head(CD_PC_DEG_4)
CD_PC_DEG_4$type = ifelse(CD_PC_DEG_4$p_val_adj < 0.05 & abs(CD_PC_DEG_4$avg_log2FC) >= 0.25,ifelse(CD_PC_DEG_4$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CD_PC_DEG_4$type)
# Down Stable     Up 
# 204   1755    389 
CD_PC_DEG_4_up <- rownames(CD_PC_DEG_4[CD_PC_DEG_4$type == "Up",])
CD_PC_DEG_4_down <- rownames(CD_PC_DEG_4[CD_PC_DEG_4$type == "Down",])

# 4. Merge up-regulated DEGs

ALL_DEG_up <- union(PT_S1_DEG_2_up,c(PT_S1_DEG_4_up,
                                     PT_S2_DEG_2_up,PT_S2_DEG_4_up,
                                     PT_S3_DEG_2_up,PT_S3_DEG_4_up,
                                     PT_injured_DEG_2_up,PT_injured_DEG_4_up,
                                     DLH_DEG_2_up,DLH_DEG_4_up,
                                     ALH_DEG_2_up,ALH_DEG_4_up,
                                     DCT_DEG_2_up,DCT_DEG_4_up,
                                     CD_IC_DEG_2_up,CD_IC_DEG_4_up,
                                     CD_PC_DEG_2_up,CD_PC_DEG_4_up))
ALL_DEG_up <- as.data.frame(ALL_DEG_up)
row.names(ALL_DEG_up) <- ALL_DEG_up$ALL_DEG_up
colnames(ALL_DEG_up) <- c("ID")
head(ALL_DEG_up)

PT_S1_DEG_2_up <- as.data.frame(PT_S1_DEG_2_up)
rownames(PT_S1_DEG_2_up) <- PT_S1_DEG_2_up$PT_S1_DEG_2_up
PT_S1_DEG_2_up$ID <- rownames(PT_S1_DEG_2_up)
head(PT_S1_DEG_2_up)

PT_S1_DEG_4_up <- as.data.frame(PT_S1_DEG_4_up)
rownames(PT_S1_DEG_4_up) <- PT_S1_DEG_4_up$PT_S1_DEG_4_up
PT_S1_DEG_4_up$ID <- rownames(PT_S1_DEG_4_up)
head(PT_S1_DEG_4_up)

PT_S2_DEG_2_up <- as.data.frame(PT_S2_DEG_2_up)
rownames(PT_S2_DEG_2_up) <- PT_S2_DEG_2_up$PT_S2_DEG_2_up
PT_S2_DEG_2_up$ID <- rownames(PT_S2_DEG_2_up)
head(PT_S2_DEG_2_up)

PT_S2_DEG_4_up <- as.data.frame(PT_S2_DEG_4_up)
rownames(PT_S2_DEG_4_up) <- PT_S2_DEG_4_up$PT_S2_DEG_4_up
PT_S2_DEG_4_up$ID <- rownames(PT_S2_DEG_4_up)
head(PT_S2_DEG_4_up)

PT_S3_DEG_2_up <- as.data.frame(PT_S3_DEG_2_up)
rownames(PT_S3_DEG_2_up) <- PT_S3_DEG_2_up$PT_S3_DEG_2_up
PT_S3_DEG_2_up$ID <- rownames(PT_S3_DEG_2_up)
head(PT_S3_DEG_2_up)

PT_S3_DEG_4_up <- as.data.frame(PT_S3_DEG_4_up)
rownames(PT_S3_DEG_4_up) <- PT_S3_DEG_4_up$PT_S3_DEG_4_up
PT_S3_DEG_4_up$ID <- rownames(PT_S3_DEG_4_up)
head(PT_S3_DEG_4_up)

PT_injured_DEG_2_up <- as.data.frame(PT_injured_DEG_2_up)
rownames(PT_injured_DEG_2_up) <- PT_injured_DEG_2_up$PT_injured_DEG_2_up
PT_injured_DEG_2_up$ID <- rownames(PT_injured_DEG_2_up)
head(PT_injured_DEG_2_up)

PT_injured_DEG_4_up <- as.data.frame(PT_injured_DEG_4_up)
rownames(PT_injured_DEG_4_up) <- PT_injured_DEG_4_up$PT_injured_DEG_4_up
PT_injured_DEG_4_up$ID <- rownames(PT_injured_DEG_4_up)
head(PT_injured_DEG_4_up)

DLH_DEG_2_up <- as.data.frame(DLH_DEG_2_up)
rownames(DLH_DEG_2_up) <- DLH_DEG_2_up$DLH_DEG_2_up
DLH_DEG_2_up$ID <- rownames(DLH_DEG_2_up)
head(DLH_DEG_2_up)

DLH_DEG_4_up <- as.data.frame(DLH_DEG_4_up)
rownames(DLH_DEG_4_up) <- DLH_DEG_4_up$DLH_DEG_4_up
DLH_DEG_4_up$ID <- rownames(DLH_DEG_4_up)
head(DLH_DEG_4_up)

ALH_DEG_2_up <- as.data.frame(ALH_DEG_2_up)
rownames(ALH_DEG_2_up) <- ALH_DEG_2_up$ALH_DEG_2_up
ALH_DEG_2_up$ID <- rownames(ALH_DEG_2_up)
head(ALH_DEG_2_up)

ALH_DEG_4_up <- as.data.frame(ALH_DEG_4_up)
rownames(ALH_DEG_4_up) <- ALH_DEG_4_up$ALH_DEG_4_up
ALH_DEG_4_up$ID <- rownames(ALH_DEG_4_up)
head(ALH_DEG_4_up)

DCT_DEG_2_up <- as.data.frame(DCT_DEG_2_up)
rownames(DCT_DEG_2_up) <- DCT_DEG_2_up$DCT_DEG_2_up
DCT_DEG_2_up$ID <- rownames(DCT_DEG_2_up)
head(DCT_DEG_2_up)

DCT_DEG_4_up <- as.data.frame(DCT_DEG_4_up)
rownames(DCT_DEG_4_up) <- DCT_DEG_4_up$DCT_DEG_4_up
DCT_DEG_4_up$ID <- rownames(DCT_DEG_4_up)
head(DCT_DEG_4_up)

CD_IC_DEG_2_up <- as.data.frame(CD_IC_DEG_2_up)
rownames(CD_IC_DEG_2_up) <- CD_IC_DEG_2_up$CD_IC_DEG_2_up
CD_IC_DEG_2_up$ID <- rownames(CD_IC_DEG_2_up)
head(CD_IC_DEG_2_up)

CD_IC_DEG_4_up <- as.data.frame(CD_IC_DEG_4_up)
rownames(CD_IC_DEG_4_up) <- CD_IC_DEG_4_up$CD_IC_DEG_4_up
CD_IC_DEG_4_up$ID <- rownames(CD_IC_DEG_4_up)
head(CD_IC_DEG_4_up)

CD_PC_DEG_2_up <- as.data.frame(CD_PC_DEG_2_up)
rownames(CD_PC_DEG_2_up) <- CD_PC_DEG_2_up$CD_PC_DEG_2_up
CD_PC_DEG_2_up$ID <- rownames(CD_PC_DEG_2_up)
head(CD_PC_DEG_2_up)

CD_PC_DEG_4_up <- as.data.frame(CD_PC_DEG_4_up)
rownames(CD_PC_DEG_4_up) <- CD_PC_DEG_4_up$CD_PC_DEG_4_up
CD_PC_DEG_4_up$ID <- rownames(CD_PC_DEG_4_up)
head(CD_PC_DEG_4_up)

merge_DEG_up <- left_join(ALL_DEG_up,PT_S1_DEG_2_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,PT_S1_DEG_4_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,PT_S2_DEG_2_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,PT_S2_DEG_4_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,PT_S3_DEG_2_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,PT_S3_DEG_4_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,PT_injured_DEG_2_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,PT_injured_DEG_4_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,DLH_DEG_2_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,DLH_DEG_4_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,ALH_DEG_2_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,ALH_DEG_4_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,DCT_DEG_2_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,DCT_DEG_4_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,CD_IC_DEG_2_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,CD_IC_DEG_4_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,CD_PC_DEG_2_up,by="ID")
merge_DEG_up <- left_join(merge_DEG_up,CD_PC_DEG_4_up,by="ID")

write.csv(merge_DEG_up,"merge_DEG_up.csv")

merge_DEG_up_t <- merge_DEG_up
rownames(merge_DEG_up_t) <- merge_DEG_up_t$ID

head(merge_DEG_up_t)
merge_DEG_up_t[which(!is.na(merge_DEG_up_t),arr.ind = T)]<-1
merge_DEG_up_t[which(is.na(merge_DEG_up_t),arr.ind = T)]<-0
head(merge_DEG_up_t)
merge_DEG_up_t <- as.data.frame(lapply(merge_DEG_up_t,as.numeric))
merge_DEG_up_t <- merge_DEG_up_t[,-1]
str(merge_DEG_up_t)

colnames(merge_DEG_up_t) <- c("PT-S1 (2w vs Con)","PT-S1 (4w vs Con)",
                              "PT-S2 (2w vs Con)","PT-S2 (4w vs Con)",
                              "PT-S3 (2w vs Con)","PT-S3 (4w vs Con)",
                              "PT-injured (2w vs Con)","PT-injured (4w vs Con)",
                              "DLH (2w vs Con)","DLH (4w vs Con)",
                              "ALH (2w vs Con)","ALH (4w vs Con)",
                              "DCT (2w vs Con)","DCT (4w vs Con)",
                              "CD-IC (2w vs Con)","CD-IC (4w vs Con)",
                              "CD-PC (2w vs Con)","CD-PC (4w vs Con)")
# 5. Upset plot

upset(merge_DEG_up_t, nsets = 9,nintersects = 20,
      mb.ratio = c(.6, 0.4),
      order.by = c("freq"), 
      # keep.order = TRUE,
      decreasing = c(T,F), 
      sets.bar.color = c("PT-S1 (2w vs Con)" = '#E4C755',"PT-S1 (4w vs Con)" = '#E4C755',
                         "PT-S2 (2w vs Con)" = '#E59CC4',"PT-S2 (4w vs Con)" = '#E59CC4',
                         "PT-S3 (2w vs Con)" = '#AB3282',"PT-S3 (4w vs Con)" = '#AB3282',
                         "PT-injured (2w vs Con)" = '#E95C59',"PT-injured (4w vs Con)" = '#E95C59',
                         "DLH (2w vs Con)" = '#53A85F',"DLH (4w vs Con)" = '#53A85F',
                         "ALH (2w vs Con)" = '#F1BB72',"ALH (4w vs Con)" = '#F1BB72',
                         "DCT (2w vs Con)" = '#F3B1A0',"DCT (4w vs Con)" = '#F3B1A0',
                         "CD-IC (2w vs Con)" = '#D6E7A3',
                         "CD-PC (2w vs Con)" = '#57C3F3'),
      mainbar.y.label = "Intersection number", 
      sets.x.label = "Up-regulated DEGs",
      text.scale = 1.5)


# 6. GO enrichment 

DEGs_up_list <- list(
  PT_S1_DEG_2_up = PT_S1_DEG_2_up$PT_S1_DEG_2_up,PT_S1_DEG_4_up = PT_S1_DEG_4_up$PT_S1_DEG_4_up,
  PT_S2_DEG_2_up = PT_S2_DEG_2_up$PT_S2_DEG_2_up,PT_S2_DEG_4_up = PT_S2_DEG_4_up$PT_S2_DEG_4_up,
  PT_S3_DEG_2_up = PT_S3_DEG_2_up$PT_S3_DEG_2_up,PT_S3_DEG_4_up = PT_S3_DEG_4_up$PT_S3_DEG_4_up,
  PT_injured_DEG_2_up = PT_injured_DEG_2_up$PT_injured_DEG_2_up,PT_injured_DEG_4_up = PT_injured_DEG_4_up$PT_injured_DEG_4_up,
  DLH_DEG_2_up = DLH_DEG_2_up$DLH_DEG_2_up,DLH_DEG_4_up = DLH_DEG_4_up$DLH_DEG_4_up,
  ALH_DEG_2_up = ALH_DEG_2_up$ALH_DEG_2_up,ALH_DEG_4_up = ALH_DEG_4_up$ALH_DEG_4_up,
  DCT_DEG_2_up = DCT_DEG_2_up$DCT_DEG_2_up,DCT_DEG_4_up = DCT_DEG_4_up$DCT_DEG_4_up,
  CD_IC_DEG_2_up = CD_IC_DEG_2_up$CD_IC_DEG_2_up,CD_IC_DEG_4_up = CD_IC_DEG_4_up$CD_IC_DEG_4_up,
  CD_PC_DEG_2_up = CD_PC_DEG_2_up$CD_PC_DEG_2_up,CD_PC_DEG_4_up = CD_PC_DEG_4_up$CD_PC_DEG_4_up)
DEGs_up_list

DEGs_up_list_GO <- compareCluster(DEGs_up_list, fun="enrichGO",OrgDb ="org.Mm.eg.db", pvalueCutoff=0.01,keyType ="SYMBOL")
dotplot(DEGs_up_list_GO,showCategory = 1) + ggtitle("Up regulated DEGs") + xlab(NULL) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + 
  scale_color_gradientn(colours = rev(viridis(20)), 
  guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"), name = "p.adjust") 
