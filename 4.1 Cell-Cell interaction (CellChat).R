# 1.Loading packages

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(ggsci)
library(VennDiagram)

# 2. ST-seq dataset cellchat

data<-readRDS("ST_anno.rds")
meta <- data@meta.data 
Control_meta <- subset(meta,type=="Control")
AAN_2W_meta <- subset(meta,type=="AAN_2W")
AAN_4W_meta <- subset(meta,type=="AAN_4W")
Control_data<-as.matrix(data@assays$SCT@data[,as.character(row.names(Control_meta))])
AAN_2W_data<-as.matrix(data@assays$SCT@data[,as.character(row.names(AAN_2W_meta))])
AAN_4W_data<-as.matrix(data@assays$SCT@data[,as.character(row.names(AAN_4W_meta))])

# 2.1 ST-seq dataset cellchat (Control)

cellchat <- createCellChat(object = Control_data, meta = Control_meta, group.by = "spot_type")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 40) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 1)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("ST_Control_cellchat.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
saveRDS(cellchat, "ST_Control_cellchat.RDS")

# 2.2 ST-seq dataset cellchat (AAN_2W)

cellchat <- createCellChat(object = AAN_2W_data, meta = AAN_2W_meta, group.by = "spot_type")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 40) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 1)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("ST_AAN_2W_cellchat.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
saveRDS(cellchat, "ST_AAN_2W_cellchat.RDS")

# 2.3 ST-seq dataset cellchat (AAN_4W)

cellchat <- createCellChat(object = AAN_4W_data, meta = AAN_4W_meta, group.by = "spot_type")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 40) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 1)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("ST_AAN_4W_cellchat.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
saveRDS(cellchat, "ST_AAN_4W_cellchat.RDS")

# 3. scRNA-seq dataset cellchat

data<-readRDS("scRNA_anno.rds")
meta<-data@meta.data
Control_meta<-subset(meta,type=="Control")
Treatment_meta<-subset(meta,type=="Treatment")
Control_data<-as.matrix(data@assays$SCT@data[,as.character(row.names(Control_meta))])
Treatment_data<-as.matrix(data@assays$SCT@data[,as.character(row.names(Treatment_meta))])

# 3.1 scRNA-seq dataset cellchat (Control)

cellchat <- createCellChat(object = Control_data, meta = Control_meta, group.by = "cell_type")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
pdf("scRNA_Control_cellchat.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions"",color.use=color$V2)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength"",color.use=color$V2)
dev.off()
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat, "scRNA_Control_cellchat.RDS")


# 3.2 scRNA-seq dataset cellchat (AAN_4W)

cellchat <- createCellChat(object = Treatment_data, meta = Treatment_meta, group.by = "cell_type")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
pdf("scRNA_AAN_4W_cellchat.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use=color$V2)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use=color$V2)
dev.off()
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")signaling pathways
saveRDS(cellchat, "scRNA_AAN_4W_cellchat.RDS")

# 4. Different interactions and pathways in ST-seq dataset

ST_control<-readRDS( "/szrmyy/wangjgLab/scRNA/baiym/P14_mouse_brain/AA_ST/Control_cellchat.RDS")
ST_AAN_2W<-readRDS( "/szrmyy/wangjgLab/scRNA/baiym/P14_mouse_brain/AA_ST/AAN_2W_cellchat.RDS") 
ST_AAN_4W<-readRDS( "/szrmyy/wangjgLab/scRNA/baiym/P14_mouse_brain/AA_ST/AAN_4W_cellchat.RDS")
object.list <- list(Control=ST_control,AAN_2W=ST_AAN_2W,AAN_4W=ST_AAN_4W)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
p1 <- netVisual_heatmap(cellchat,comparison=c("Control","AAN_2W"),title.name = "Different number of interactions of AA_2W vs Control\nin ST-seq dataset",font.size=16,font.size.title = 18)
p2 <- netVisual_heatmap(cellchat,comparison=c("Control","AAN_4W"),title.name = "Different number of interactions of AA_4W vs Control\nin ST-seq dataset",font.size=16,font.size.title = 18)
pdf("ST_Diff_Heatmap_select.pdf")
p1
p2 
dev.off()

pdf("ST_Diff_Path_2TO2.pdf",height=12)
rankNet(cellchat, mode = "comparison", stacked = T,do.stat =TRUE,comparison = c(1,2),color.use=c("#00AFBB", "#E7B800", "#FC4E07")[c(1,2)],title = "Different pathways of AA_2W vs Control\nin ST-seq dataset")+coord_flip()
rankNet(cellchat, mode = "comparison", stacked = T,do.stat =TRUE,comparison = c(1,3),color.use=c("#00AFBB", "#E7B800", "#FC4E07")[c(1,3)],title = "Different pathways of AA_4W vs Control\nin ST-seq dataset")+coord_flip()
dev.off()

# 5. Different interactions and pathways in scRNA-seq dataset

control<-readRDS("scRNA_Control_cellchat.RDS")
treatment<-readRDS("scRNA_AAN_4W_cellchat.RDS")
object.list <- list(Control=control,AAN_4W=treatment)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
pdf("sc_Diff_Num.pdf")
compareInteractions(cellchat, show.legend = F, group = c(1,2))
dev.off()
pdf("sc_Diff_Heatmap.pdf")
netVisual_heatmap(cellchat,comparison=c(1,2),title.name = "Different number of interactions of AA_4W vs Control\nin scRNA-seq dataset",font.size=16,font.size.title = 18)
dev.off()
pdf("sc_Diff_Path_2TO2.pdf",height=12)
rankNet(cellchat, mode = "comparison", stacked = T,do.stat =TRUE,comparison = c(1,2),color.use=c("#3B4992","#008B45","#EE0000")[c(1,3)],title = "Different pathways of AA_4W vs Control\nin scRNA-seq dataset")+coord_flip()
dev.off()

# 6. The interacted differentially expressed pathways in ST-seq and scRNA-seq datasets

ST_AAN_2WvsControl_Path=c(
"COMPLEMENT","MHC-I","ICAM","ITGAL-ITGB2","CX3C","CDH1","SEMA7","THY1","CD45","NPNT","TRAIL","CD226","EDN","CALCR","CD39","ALCAM","CD6","PVR","LEP","SEMA6","IL2","LCK","ACTIVIN","NT","NEGR","SN","VCAM","PROS","CCL","GALECTIN","GDF","TGFb","TWEAK","SEMA5","NOTCH","SEMA3","CDH5","CSF","FN1","GAS","MIF","RESISTIN","L1CAM","KIT","NCAM","CSPG4","PERIOSTIN","PTN","VISFATIN","BMP","THBS","IGF","LIFR","EPHA","WNT","ADIPONECTIN","ANGPTL","LAMININ","COLLAGEN","HH","NRG","AGT","EPHB","ncWNT","PDGF","NECTIN","CDH","HSPG","GRN","JAM","TENASCIN","SPP1"
)
ST_AAN_4WvsControl_Path=c(
"MHC-I","COMPLEMENT","ICAM","ITGAL-ITGB2","MHC-II","THY1","CX3C","CD45","SEMA7","CDH1","NPNT","CD226","SEMA6","ALCAM","CD6","CD39","CHEMERIN","EDN","LCK","BST2","RELN","TRAIL","IL2","OCLN","CALCR","CD200","ACTIVIN","PVR","APJ","IL16","APRIL","BAFF","LAIR1","NEGR","NT","APELIN","IL1","PD-L1","CCL","VCAM","PROS","GALECTIN","GDF","SEMA5","TGFb","KIT","CSF","CDH5","TWEAK","SEMA3","NRG","CADM","RESISTIN","NOTCH","LIFR","NCAM","MIF","FN1","ncWNT","THBS","VISTA","BMP","GAS","PERIOSTIN","EPHB","HH","ANGPTL","EPHA","VISFATIN","WNT","PTN","CSPG4","LAMININ","L1CAM","PECAM1","ADIPONECTIN","CDH","CXCL","ANGPT","MPZ","PDGF","JAM","COLLAGEN","NECTIN","HSPG","SEMA4","IGF","EGF","SPP1","TENASCIN","GRN","NPY","AGRN","MK","ESAM","VEGF","APP","FGF"
)
SC_AAN_4WvsControl_Path=c(
"MHC-I","MHC-II","MK","THY1","CD86","SEMA7","ALCAM","CD6","EGF","ANNEXIN","CD80","IFN-II","XCR","LCK","CD226","EPHB","FASLG","PVR","SEMA6","PD-L1","PTPRM","TWEAK","VCAM","CD45","SPP1","GAS","CCL","PROS","KIT","GDF","ICAM","JAM","GALECTIN"
)

Path_ST=intersect(ST_AAN_2WvsControl_Path,intersect(ST_AAN_4WvsControl_Path,SC_AAN_4WvsControl_Path))
plot <- venn.diagram(x=list(ST_4WvsControl=ST_AAN_4WvsControl_Path,ST_2WvsControl=ST_AAN_2WvsControl_Path,SC_4WvsControl=SC_AAN_4WvsControl_Path),filename=NULL,fill=c("#00AFBB", "#E7B800", "#FC4E07"))
pdf("ST_SC_3Com_venn.pdf",height=5,width=5)
grid.draw(plot)
dev.off()

pdf("ST_SC_3Com_Path_ST.pdf",height=5,width=5)
rankNet(cellchat, mode = "comparison", stacked = T,do.stat =TRUE,comparison = c(1,2,3),color.use=c("#00AFBB", "#E7B800", "#FC4E07"),signaling=c(Path_ST),title = "Different pathways of AA vs Control\nin ST-seq dataset")+coord_flip()
dev.off()

# 7. Ligand_Receptor interaction in ST-seq dataset

select_pathway <- c("SPP1","TWEAK","CCL","MHC-I","MHC-II","GAS","PROS","GALECTIN","VCAM")
pairLR.use <- extractEnrichedLR(cellchat, signaling = select_pathway)
pairLR.use  <- subset(pairLR.use, interaction_name %in% c("SPP1_CD44","SPP1_ITGAV_ITGB1","SPP1_ITGAV_ITGB5","TNFSF12_TNFRSF12A","CCL8_CCR2","CCL2_CCR2","CCL6_CCR2","CCL12_CCR2","CCL5_CCR5","CCL8_CCR5","H2-D1_CD8A","H2-K1_CD8A","H2-D1_CD8B1","H2-K1_CD8B1","H2-AB1_CD4","H2-EB1_CD4","ITGA9_ITGB1_VCAM1","ITGA4_ITGB1_VCAM1","PROS1_AXL","GAS6_AXL","LGALS9_CD44"))
out <- netVisual_bubble(cellchat, comparison = c(1, 2, 3),  angle.x = 45, remove.isolate = F, pairLR.use = pairLR.use, color.text = c("#00AFBB", "#E7B800", "#FC4E07" ),return.data=TRUE)
source.target<-levels(out$communication$source.target)
select<-source.target[c(57:59,88:89,119:120,154:156,334:344,358:359)]
source("./netVisual_bubble2.r")
out2 <- netVisual_bubble2(cellchat, comparison = c(1, 2, 3),  angle.x = 45, remove.isolate = F, pairLR.use = pairLR.use, color.text = c("#00AFBB", "#E7B800", "#FC4E07" ),select.source.target=select,return.data=TRUE)
source.target2<-levels(out2$communication$source.target)
xgap <- c(grep("AAN_4W",source.target2)+0.5)
pdf("Path9_Select_celltype_LR.pdf",height=6,width=11)
netVisual_bubble2(cellchat, comparison = c(1, 2, 3),  angle.x = 45, remove.isolate = T, pairLR.use = pairLR.use, color.text = c("#00AFBB", "#E7B800", "#FC4E07" ),select.source.target=select)+geom_vline(xintercept=xgap[-length(xgap)],linetype=2,color="grey")+theme(text=element_text(size=15),legend.title=element_text(size=16),legend.text=element_text(size=15))
dev.off()
DefaultAssay(data) <- "SCT"
pdf("Path9_SelectLR_Vln.pdf",height=10,width=10)
VlnPlot(data,features=unique(gene),split.by="type",pt.size=0,ncol=1)*theme_void()*NoLegend()*scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))
dev.off()
