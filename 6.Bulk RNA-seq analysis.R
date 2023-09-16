# 1.Loading packages

library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggpubr)
library(GSVA)
library(enrichplot)
library(ggsci)

# 2.Pre-process the Bulk RNA-seq dataset
raw_data <- read.table("AA_Bulk_RNA_process.txt",check.names = F,header = T)
head(raw_data)

gene.df <- bitr(raw_data$Geneid, fromType ="ENSEMBL",toType = "SYMBOL", OrgDb = org.Mm.eg.db)  
merge_data <- data.frame(gene.df,raw_data[match(gene.df$ENSEMBL,raw_data$Geneid),])

table(!duplicated(merge_data$SYMBOL))
dulpicated_list <- which(duplicated(merge_data$SYMBOL))
merge_data_rd <- merge_data[-dulpicated_list,]
row.names(merge_data_rd) <- merge_data_rd$SYMBOL

data <- merge_data_rd[,c(4,5:13)]
write.csv(data,"Bulk RNA counts.csv")

# 3.Marker expression in the Bulk RNA-seq dataset

data2 <- data[,c(2:10)]
logCPM_all <- cpm(data2, log=TRUE, prior.count=3)
logCPM_all

my_comparisons <- list(c("AAN_2W", "Control"), c("AAN_4W", "Control"))

Havcr1 <- data.frame(logCPM_all["Havcr1",])
colnames(Havcr1)[1] <- "expression"
Havcr1$type <- c(rep("Control",3),rep("AAN_2W",3),rep("AAN_4W",3))
head(Havcr1)

p1 <- ggboxplot(Havcr1, x="type", y="expression", fill = "type", 
                palette = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  ggtitle("Havcr1") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test",label = "p.signif",) +  
  xlab("") + ylab("gene expression level (CPM)") + theme(legend.position = 'none') +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size =1) 
p1

Lcn2 <- data.frame(logCPM_all["Lcn2",])
colnames(Lcn2)[1] <- "expression"
Lcn2$type <- c(rep("Control",3),rep("AAN_2W",3),rep("AAN_4W",3))
head(Lcn2)

p2 <- ggboxplot(Lcn2, x="type", y="expression", fill = "type", 
                palette = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  ggtitle("Lcn2") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test",label = "p.signif",) +  
  xlab("") + ylab("gene expression level (CPM)") + theme(legend.position = 'none') +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size =1) 
p2

Timp2 <- data.frame(logCPM_all["Timp2",])
colnames(Timp2)[1] <- "expression"
Timp2$type <- c(rep("Control",3),rep("AAN_2W",3),rep("AAN_4W",3))
head(Timp2)

p3 <- ggboxplot(Timp2, x="type", y="expression", fill = "type", 
                palette = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  ggtitle("Timp2") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test",label = "p.signif") +   
  xlab("") + ylab("gene expression level (CPM)") + theme(legend.position = 'none') +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size =1) 
p3

Ifng <- data.frame(logCPM_all["Ifng",])
colnames(Ifng)[1] <- "expression"
Ifng$type <- c(rep("Control",3),rep("AAN_2W",3),rep("AAN_4W",3))
head(Ifng)

p4 <- ggboxplot(Ifng, x="type", y="expression", fill = "type", 
                palette = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  ggtitle("Ifng") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test",label = "p.signif") +  
  xlab("") + ylab("gene expression level (CPM)") + theme(legend.position = 'none') +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size =1) 
p4

Nfkb1 <- data.frame(logCPM_all["Nfkb1",])
colnames(Nfkb1)[1] <- "expression"
Nfkb1$type <- c(rep("Control",3),rep("AAN_2W",3),rep("AAN_4W",3))
head(Nfkb1)

p5 <- ggboxplot(Nfkb1, x="type", y="expression", fill = "type", 
                palette = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  ggtitle("Nfkb1") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test",label = "p.signif") +   
  xlab("") + ylab("gene expression level (CPM)") + theme(legend.position = 'none') +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size =1) 
p5

p1 | p2 | p3 | p4 | p5

# 4.PCA analysis

data <- read.csv("Bulk RNA counts.csv",row.names = 1)
data2 <- data[,-1]
head(data2)

logCPM <- cpm(data2,log=TRUE, prior.count=3)
logCPM

condition <- c(rep("Control",3),rep("AAN_2W",3),rep("AAN_4W",3))
treatment <- c(rep("Control",3),rep("AAN",6))

RNA_pca <- prcomp(logCPM, scale = T, retx=T)
summary(RNA_pca)
RNA_pc <- RNA_pca[["rotation"]]
RNA_pc_sub  <- as.data.frame(RNA_pc[,c(1,2)])
RNA_pc_data <- cbind(RNA_pc_sub,condition,treatment)
head(RNA_pc_data)

RNA_pc_data$condition <- factor(RNA_pc_data$condition,levels = c("Control","AAN_2W","AAN_4W"))
RNA_pc_data$treatment <- factor(RNA_pc_data$treatment,levels = rev(c("Control","AAN")))

ggplot(RNA_pc_data,aes(x=PC1,y=PC2,color=condition,label=rownames(RNA_pc_data))) + 
  stat_ellipse(aes(x=PC1,y=PC2,fill=treatment),type="norm",geom="polygon",alpha=0.2,color=NA) +
  geom_point(size = 4) +
  # geom_text(size=3,color = "black") +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  xlab(paste("PC1 (",95.76,"%","variance)",sep="")) +
  ylab(paste("PC2 (",2.41,"%","variance)",sep="")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) 

# 5.DEGs in Bulk-RNA (AAN_2W vs Control)  

counts <- data[,c(2:7)] # AAN_2W vs Control
head(counts)

condition = factor(c(rep("Control",3),rep("AAN_2W",3)),levels = c("Control","AAN_2W"))
condition

genelist = DGEList(counts = counts, group = condition)
design <- model.matrix(~condition)
colnames(design) <- levels(condition)
rownames(design) <- colnames(counts)
design

keep <- rowSums(cpm(genelist)>1)>=2 
genelist.filted <- genelist[keep,,keep.lib.sizes=FALSE] 
genelist.filted

genelist.norm <- calcNormFactors(genelist.filted)  
logCPM <- cpm(genelist.norm, log=TRUE, prior.count=3)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
DEG_result <- topTable(fit,coef=2,n=Inf)
DEG_result
DEG_result$type = ifelse(DEG_result$adj.P.Val < 0.05 & abs(DEG_result$logFC) >= 1,
                         ifelse(DEG_result$logFC > 1 ,'Up (2265)','Down (1739)'),'Stable (10491)')
table(DEG_result$type)
# Down (1739) Stable (10491)      Up (2265) 
# 1739          10491           2265

write.csv(DEG_result,"DEGs in Bulk-RNA Seq(AAN_2W vs Control).csv" )

DEG_up <- row.names(DEG_result[DEG_result$type == "Up (2265)",])  
head(DEG_up)

DEG_down <- row.names(DEG_result[DEG_result$type == "Down (1739)",]) 
head(DEG_down)

DEG_result$label <- ifelse(DEG_result$adj.P.Val < 0.01 & abs(DEG_result$logFC) >= 6,
                           rownames(DEG_result),"")
table(DEG_result$label)

p1 <- ggplot(data=DEG_result, aes(x=logFC, y =-log10(adj.P.Val),colour=type)) +
  geom_point(alpha=0.8, size=1)+
  scale_color_manual(values=c("#546de5","#000000","#ff4757"))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.8,alpha=0.5)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.8,alpha=0.5)+
  labs(x="log2 (fold change)",y="-log10 (FDR)")+
  ggtitle("Signifianct DEGs in Bulk RNA-seq \n (AAN_2W vs Control)")+ 
  theme_bw() + xlim(-15,15) + ylim(0, 5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_text_repel(data = DEG_result,aes(x=logFC,y =-log10(adj.P.Val),label = label),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"),
                  segment.color = "black",show.legend = FALSE,max.overlaps =40)  # 6x5
p1

# 6.GO enrichment in Bulk RNA-seq (AAN_2W vs Control)

Go_BP_up <- enrichGO(gene = DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = 'BH',
                     pvalueCutoff = 0.05,qvalueCutoff = 0.05)

Go_BP_down <- enrichGO(gene = DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                       ont = "BP", pAdjustMethod = 'BH',
                       pvalueCutoff = 0.05,qvalueCutoff = 0.05)

dotplot(Go_BP_up, showCategory = 10)

dotplot(Go_BP_down, showCategory = 10)

Go_BP_up@result$Description[1:20]
Go_BP_up_list <- Go_BP_up@result[c(1,2,8,10,12),c("Description","pvalue")]
Go_BP_up_list$type <- "Up-regulated"

Go_BP_down@result$Description[1:20]
Go_BP_down_list <-Go_BP_down@result[c(1,2,3,4,5),c("Description","pvalue")]
Go_BP_down_list$type <- "Down-regulated"

Go_BP_list <- rbind(Go_BP_up_list,Go_BP_down_list)
Go_BP_list$log10Pvalue <- -log10(Go_BP_list$pvalue)

# Go_BP_list$type <- factor(Go_BP_list$type,levels = c("Up-regulated","Down-regulated"))
p2 <- ggpubr::ggbarplot(Go_BP_list, x="Description", y="log10Pvalue", fill = "type", color = "white",
                        palette =  c("Up-regulated" = "#ff4757", "Down-regulated" = "#546de5"),
                        sort.val = "asc", 
                        sort.by.grodowns=TRUE, 
                        x.text.angle=0, 
                        ylab = "-log10Pvalue", xlab = " ") + coord_flip() # 6x5

p1 + p2 # 12*5

# 7.DEGs in Bulk-RNA (AAN_4W vs Control)  
head(data)
counts <- data[,c(2:4,8:10)] # AAN_4W vs Control
head(counts)

condition = factor(c(rep("Control",3),rep("AAN_4W",3)),levels = c("Control","AAN_4W"))
condition

genelist = DGEList(counts = counts, group = condition)
design <- model.matrix(~condition)
colnames(design) <- levels(condition)
rownames(design) <- colnames(counts)
design

keep <- rowSums(cpm(genelist)>1)>=2 
genelist.filted <- genelist[keep,,keep.lib.sizes=FALSE] 
genelist.filted

genelist.norm <- calcNormFactors(genelist.filted)  
logCPM <- cpm(genelist.norm, log=TRUE, prior.count=3)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
DEG_result <- topTable(fit,coef=2,n=Inf)
DEG_result
DEG_result$type = ifelse(DEG_result$adj.P.Val < 0.05 & abs(DEG_result$logFC) >= 1,
                         ifelse(DEG_result$logFC > 1 ,'Up(1479)','Down(1075)'),'Stable(11961)')
table(DEG_result$type)
# Down(1075) Stable(11961)      Up(1479) 
# 1075         11961          1479

write.csv(DEG_result,"DEGs in Bulk-RNA Seq(AAN_4W vs Control).csv" )

DEG_up <- row.names(DEG_result[DEG_result$type == "Up(1479)",])  
head(DEG_up)
DEG_down <- row.names(DEG_result[DEG_result$type == "Down(1075)",]) 
head(DEG_down)

DEG_result$label <- ifelse(DEG_result$adj.P.Val < 0.05 & abs(DEG_result$logFC) >= 5,rownames(DEG_result),"")

p4 <- ggplot(data=DEG_result, aes(x=logFC, y =-log10(adj.P.Val),colour=type)) +
  geom_point(alpha=0.8, size=1)+
  scale_color_manual(values=c("#546de5","#000000","#ff4757"))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.8,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.8,alpha=0.8)+
  labs(x="log2 (fold change)",y="-log10 (FDR)")+
  ggtitle("Signifianct DEGs in Bulk RNA-seq \n (AAN_4W vs Control)")+ 
  theme_bw() + xlim(-10,10) + ylim(0, 5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_text_repel(data = DEG_result,aes(x=logFC,y =-log10(adj.P.Val),label = label),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"),
                  segment.color = "black",show.legend = FALSE,max.overlaps =40)  # 6x5

# 8.GO enrichment in Bulk RNA-seq (AAN_4W vs Control)

Go_BP_up <- enrichGO(gene = DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

Go_BP_down <- enrichGO(gene = DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                       ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

dotplot(Go_BP_up, showCategory = 10,title = "The GO enrichment(BP) analysis of overlap up DEGs")

dotplot(Go_BP_down, showCategory = 10,title = "The GO enrichment(BP) analysis of overlap down DEGs")

Go_BP_up@result$Description[1:20]
Go_BP_up_list <- Go_BP_up@result[c(1,3,6,12,16),c("Description","pvalue")]
Go_BP_up_list$type <- "Up-regulated"

Go_BP_down@result$Description[1:20]
Go_BP_down_list <-Go_BP_down@result[c(1,2,3,4,5),c("Description","pvalue")]
Go_BP_down_list$type <- "Down-regulated"

Go_BP_list <- rbind(Go_BP_up_list,Go_BP_down_list)
Go_BP_list$log10Pvalue <- -log10(Go_BP_list$pvalue)

p5 <- ggbarplot(Go_BP_list, x="Description", y="log10Pvalue", fill = "type", color = "white",
                palette =  c("Up-regulated" = "#ff4757", "Down-regulated" = "#546de5"),
                sort.val = "asc", 
                sort.by.grodowns=TRUE, 
                x.text.angle=0, 
                ylab = "-log10Pvalue", xlab = " ") + coord_flip() # 8x5

(p1 + p2)/ (p4 + p5)

# 9. GSEA analysis

gmt <- read.gmt("mouse_HMouse.gmt")  

# AAN_2W vs Control

DEG_result <- read.csv("DEGs in Bulk-RNA Seq(AAN_2W vs Control).csv",row.names = 1) 
DEG_result
DEG_result_sort <- DEG_result %>% arrange(desc(logFC))
head(DEG_result_sort)

geneList = DEG_result_sort$logFC
names(geneList) <- rownames(DEG_result_sort)
head(geneList)

gsea_AAN_2W <- GSEA(geneList,TERM2GENE = gmt)

p1 <- gseaplot2(gsea_AAN_2W,
          c("HALLMARK_APOPTOSIS",
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "HALLMARK_INFLAMMATORY_RESPONSE",
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB"), 
          pvalue_table = T,color = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")) 
p1

# AAN_4W vs Control

DEG_result <- read.csv("DEGs in Bulk-RNA Seq(AAN_4W vs Control).csv",row.names = 1) 
DEG_result
DEG_result_sort <- DEG_result %>% arrange(desc(logFC))
head(DEG_result_sort)

geneList = DEG_result_sort$logFC
names(geneList) <- rownames(DEG_result_sort)
head(geneList)

gsea_AAN_4W <- GSEA(geneList,TERM2GENE = gmt)

p2 <- gseaplot2(gsea_AAN_2W,
                c("HALLMARK_APOPTOSIS",
                  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                  "HALLMARK_INFLAMMATORY_RESPONSE",
                  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"), 
                pvalue_table = T,title = "AAN_2W vs. Control",
                color = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")) 
p2

p1 | p2  
