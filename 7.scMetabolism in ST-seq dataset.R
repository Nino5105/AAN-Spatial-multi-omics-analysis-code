# 1.Loading packages

library(scMetabolism)
library(ggplot2)
library(rsvd)
library(Seurat)
library(homologene)
library(ggsci)
library(ggrepel)
library(dplyr)
library(viridis)

# 2.Loading datasets and transform gene ID

ST <- readRDS("ST_anno.rds")
ST_genename <- ST@assays$Spatial@data@Dimnames[[1]]
ST_genename_m2h <- homologene(ST_genename, inTax = 10090, outTax = 9606)

ST@assays$Spatial@data <- ST@assays$Spatial@data[ST_genename_m2h$`10090`,]
row.names(ST@assays$Spatial@data) <- ST_genename_m2h$`9606`

ST@assays$Spatial@counts <- ST@assays$Spatial@counts[ST_genename_m2h$`10090`,]
row.names(ST@assays$Spatial@counts) <- ST_genename_m2h$`9606`

ST@assays$RNA <- ST@assays$Spatial

ST@assays$Spatial
ST@assays$SCT@counts
ST@assays$SCT@data
ST@active.assay <- "SCT"

countexp.ST <- sc.metabolism.Seurat(obj = ST,method = "VISION",
                                    imputation = F,ncores = 10,
                                        metabolism.type = "KEGG")

saveRDS(countexp.ST,"countexp.ST.rds")

input.pathway <- rownames(countexp.ST@assays[["METABOLISM"]][["score"]])

countexp.ST_con <- subset(countexp.ST,idents = type == "Control")

DotPlot.metabolism(obj = countexp.ST, pathway = input.pathway, 
                   phenotype = "spot_type", norm = "y")

#
pathway <- rownames(countexp.ST@assays[["METABOLISM"]][["score"]])
write.csv(pathway,"pathway.csv")
DimPlot.metabolism(obj = countexp.ST, pathway = "Glycolysis / Gluconeogenesis", 
                   dimention.reduction.type = "umap", dimention.reduction.run = F, 
                   size = 1)

meta_score <- countexp.ST@assays$METABOLISM$score

countexp.ST@active.assay <- "METABOLISM"
write.csv(meta_score,"meta_score.csv")


countexp.ST

DotPlot.metabolism.2 = function (obj, pathway, phenotype, norm = "y") 
{
  input.norm = norm
  input.pathway <- as.character(pathway)
  input.parameter <- phenotype
  metadata <- obj@meta.data
  metabolism.matrix <- obj@assays$METABOLISM$score
  metadata[, input.parameter] <- as.character(metadata[, input.parameter])
  metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, 
  ])
  gg_table <- c()
  for (i in 1:length(input.pathway)) {
    gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], 
                                      input.pathway[i], metabolism.matrix_sub[, i]))
  }
  gg_table <- data.frame(gg_table)
  gg_table_median <- c()
  input.group.x <- unique(as.character(gg_table[, 1]))
  input.group.y <- unique(as.character(gg_table[, 2]))
  for (x in 1:length(input.group.x)) {
    for (y in 1:length(input.group.y)) {
      gg_table_sub <- subset(gg_table, gg_table[, 1] == 
                               input.group.x[x] & gg_table[, 2] == input.group.y[y])
      gg_table_median <- rbind(gg_table_median, cbind(input.group.x[x], 
                                                      input.group.y[y], median(as.numeric(as.character(gg_table_sub[, 
                                                                                                                    3])))))
    }
  }
  gg_table_median <- data.frame(gg_table_median)
  gg_table_median[, 3] <- as.numeric(as.character(gg_table_median[, 
                                                                  3]))
  gg_table_median_norm <- c()
  input.group.x <- unique(as.character(gg_table[, 1]))
  input.group.y <- unique(as.character(gg_table[, 2]))
  range01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
  }
  if (input.norm == "y") 
    for (y in 1:length(input.group.y)) {
      gg_table_median_sub <- subset(gg_table_median, gg_table_median[, 
                                                                     2] == input.group.y[y])
      norm_value <- range01(as.numeric(as.character(gg_table_median_sub[, 
                                                                        3])))
      gg_table_median_sub[, 3] <- norm_value
      gg_table_median_norm <- rbind(gg_table_median_norm, 
                                    gg_table_median_sub)
    }
  if (input.norm == "x") 
    for (x in 1:length(input.group.x)) {
      gg_table_median_sub <- subset(gg_table_median, gg_table_median[, 
                                                                     1] == input.group.x[x])
      norm_value <- range01(as.numeric(as.character(gg_table_median_sub[, 
                                                                        3])))
      gg_table_median_sub[, 3] <- norm_value
      gg_table_median_norm <- rbind(gg_table_median_norm, 
                                    gg_table_median_sub)
    }
  if (input.norm == "na") 
    gg_table_median_norm <- gg_table_median
  gg_table_median_norm <- data.frame(gg_table_median_norm)
  gg_table_median_norm[, 3] <- as.numeric(as.character(gg_table_median_norm[, 
                                                                            3]))
  library(wesanderson)
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  if(is.factor(pathway)){
    gg_table_median_norm$X2 = factor(gg_table_median_norm$X2 ,levels = levels(pathway))
  }
  
  if(is.factor(countexp.Seurat@meta.data[,phenotype])){
    gg_table_median_norm$X1 = factor(gg_table_median_norm$X1 ,
                                     levels = levels(countexp.Seurat@meta.data[,phenotype]))
  }
  
  ggplot(data = gg_table_median_norm, aes(x = gg_table_median_norm[, 
                                                                   1], y = gg_table_median_norm[, 2], color = gg_table_median_norm[, 
                                                                                                                                   3])) + geom_point(data = gg_table_median_norm, aes(size = gg_table_median_norm[, 
                                                                                                                                                                                                                  3])) + ylab("Metabolic Pathway") + xlab(input.parameter) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, 
                                                  hjust = 1), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    scale_color_gradientn(colours = pal) + labs(color = "Value", 
                                                size = "Value") + NULL
}


input.pathway <-  c("Glycolysis / Gluconeogenesis",
                    "Citrate cycle (TCA cycle)",
                    "Oxidative phosphorylation",
                    "Fatty acid biosynthesis",
                    "Fatty acid degradation",
                    "Purine metabolism",
                    "Pyrimidine metabolism",
                    "Glycine, serine and threonine metabolism",
                    "Arginine and proline metabolism",
                    "Taurine and hypotaurine metabolism",
                    "Metabolism of xenobiotics by cytochrome P450",
                    "Drug metabolism - cytochrome P450")
  
input.pathway <-  factor(input.pathway,levels = rev(input.pathway))

countexp.ST$spot_type
DotPlot.metabolism.2(obj = countexp.ST,
                    pathway = input.pathway,
                    phenotype = "spot_type",
                    norm = "y") +
  scale_color_gradientn(colours = viridis(50), 
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"), 
                        name = "Enriched socre")



ST_meta <- ST@meta.data
table(ST_meta$type)
control_list <- rownames(ST_meta[ST_meta$type == "Control",]) # 5926
AAN_2W_list <- rownames(ST_meta[ST_meta$type == "AAN_2W",]) # 6079
AAN_4W_list <- rownames(ST_meta[ST_meta$type == "AAN_4W",]) # 6668

meta_score
meta_score_log <- log2(meta_score+1)
colnames(meta_score_log) <- rownames(ST_meta)

AAN_2W_vs_Con <- meta_score_log[,c(control_list,AAN_2W_list)]
condition = factor(c(rep("Control",5926), rep("AAN_2W",6079)),levels = c("Control","AAN_2W"))  
design <- model.matrix(~condition)
colnames(design) <- c("Control","AAN_2W")
row.names(design) <- colnames(AAN_2W_vs_Con)
head(design)

fit <- lmFit(AAN_2W_vs_Con, design)
fit <- eBayes(fit, trend=TRUE)
res <- topTable(fit,coef=2, number=Inf)
res$change = ifelse(res$adj.P.Val < 0.05 & abs(res$logFC) >= 0.25, ifelse(res$logFC> 0.25 ,'Up','Down'),'Stable')
table(res$change)
# Down Stable 
# 38     41
AAN_2W_vs_Con_res <- res
rownames(AAN_2W_vs_Con_res)
write.csv(AAN_2W_vs_Con_res,"AAN_2W_vs_Con_res.csv")

AAN_4W_vs_Con <- meta_score_log[,c(control_list,AAN_4W_list)]
condition = factor(c(rep("Control",5926), rep("AAN_4W",6668)),levels = c("Control","AAN_4W"))  
design <- model.matrix(~condition)
colnames(design) <- c("Control","AAN_4W")
row.names(design) <- colnames(AAN_4W_vs_Con)
head(design)

fit <- lmFit(AAN_4W_vs_Con, design)
fit <- eBayes(fit, trend=TRUE)
res <- topTable(fit,coef=2, number=Inf)
res$change = ifelse(res$adj.P.Val < 0.05 & abs(res$logFC) >= 0.25, ifelse(res$logFC> 0.25 ,'Up','Down'),'Stable')
table(res$change)
# Down Stable 
# 43     36
AAN_4W_vs_Con_res <- res
write.csv(AAN_4W_vs_Con_res,"AAN_4W_vs_Con_res.csv")

# overlap
overlap_down <- intersect(row.names(AAN_2W_vs_Con_res[AAN_2W_vs_Con_res$change == "Down",]),
                          row.names(AAN_4W_vs_Con_res[AAN_4W_vs_Con_res$change == "Down",]))
overlap_down # 37 

ST_meta_sub <- ST_meta[,c("spot_type","type")]

ann_colors = list(type = c(Control="#00AFBB",AAN_2W= "#E7B800",AAN_4W="#FC4E07"), 
                  spot_type = c("Glom" = '#E63863',
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
                                "Adipo" = "#23452F"))

pheatmap(meta_score_log[c(overlap_down),], gaps_col = c(5926,12005),
         show_rownames = T, show_colnames = F,cluster_cols = F,cluster_rows = F,
         annotation_col = ST_meta_sub,scale = "row",annotation_colors = ann_colors,
         treeheight_row = 0,treeheight_col = 0,border_color = "black",magma(50))


countexp.ST <- readRDS("countexp.ST.rds")
countexp.ST@active.assay <- "SCT"
DimPlot(countexp.ST)
head(countexp.ST@assays$METABOLISM$score)[1:5,1:5]
head(t(countexp.ST@assays$METABOLISM$score))[1:5,1:5]
colnames(countexp.ST@assays$METABOLISM$score) <- row.names(countexp.ST@meta.data)


table(countexp.ST$spot_type)


countexp.ST

countexp.ST2 <- subset(countexp.ST,subset = spot_type %in% c("PT-S1","PT-S2","PT-S3","PT-injured"))
countexp.ST2 <- AddMetaData(countexp.ST2,t(countexp.ST@assays$METABOLISM$score)[,c(1,2,15,18,20,33,34,36,40,48,77,78)],
                            col.name = c("Glycolysis/Gluconeogenesis",
                                         "Citrate cycle (TCA cycle)",
                                         "Oxidative phosphorylation",
                                         "Fatty acid biosynthesis",
                                         "Fatty acid degradation",
                                         "Purine metabolism",
                                         "Pyrimidine metabolism",
                                         "Glycine, serine and threonine metabolism",
                                         "Arginine and proline metabolism",
                                         "Taurine and hypotaurine metabolism",
                                         "Metabolism of xenobiotics by cytochrome P450",
                                         "Drug metabolism - cytochrome P450"))
countexp.ST2
VlnPlot(countexp.ST2,features = c("Glycolysis.Gluconeogenesis",
                                  "Citrate.cycle..TCA.cycle.",
                                  "Oxidative.phosphorylation",
                                  "Fatty.acid.biosynthesis",
                                  "Fatty.acid.degradation",
                                  "Purine.metabolism",
                                  "Pyrimidine.metabolism",
                                  "Glycine..serine.and.threonine.metabolism",
                                  "Arginine.and.proline.metabolism",
                                  "Taurine.and.hypotaurine.metabolism",
                                  "Metabolism.of.xenobiotics.by.cytochrome.P450",
                                  "Drug.metabolism...cytochrome.P450"),
        group.by = c("spot_type"),
        cols = c(Control="#00AFBB",AAN_2W= "#E7B800",AAN_4W="#FC4E07"),
        split.by = "type",pt.size = 0,
        ncol = 1,stack = T,flip = T) + xlab("")


row.names(countexp.ST@assays$METABOLISM$score)
