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
ST[["RNA"]] <- ST[["Spatial"]]
DefaultAssay(ST) <- "RNA"

ST_genename <- row.names(GetAssayData(ST,assay="RNA",slot="data"))
ST_genename_m2h <- homologene(ST_genename, inTax = 10090, outTax = 9606)
ST_genename_m2h <- ST_genename_m2h[!duplicated(ST_genename_m2h[,1]),]

ST_genename_data <- GetAssayData(ST,assay="RNA",slot="data")[ST_genename_m2h$`10090`,]
ST <- ST[ST_genename_m2h$`10090`,]
ST <- SetAssayData(ST, assay = "RNA", slot = "data", new.data = ST_genename_data)

source("RenameGenesSeurat.r")#suggestions from https://github.com/satijalab/seurat/issues/2617
ST <- RenameGenesSeurat(ST, newnames=ST_genename_m2h$`9606`)

ST@assays$RNA@counts
ST@assays$RNA@data

# 3.Perform scMetabolism analysis

countexp.ST <- sc.metabolism.Seurat(obj = ST,method = "VISION", # default
                                    imputation = F,ncores = 10,
                                    metabolism.type = "KEGG")
countexp.ST

saveRDS(countexp.ST,"countexp.ST.rds")

input.pathway <- rownames(countexp.ST@assays[["METABOLISM"]][["score"]]) # n = 85 
DotPlot.metabolism(obj = countexp.ST, pathway = input.pathway, 
                   phenotype = "spot_type", norm = "y")

# countexp.ST@active.assay <- "METABOLISM"
meta_score <- countexp.ST@assays$METABOLISM$score
write.csv(meta_score,"meta_score.csv")

# 4.Visualization of metabolic pathways across all spot types

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
  
  if(is.factor(obj@meta.data[,phenotype])){
    gg_table_median_norm$X1 = factor(gg_table_median_norm$X1 ,
                                     levels = levels(obj@meta.data[,phenotype]))
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




countexp.ST <- readRDS("countexp.ST.rds")
countexp.ST@active.assay <- "SCT"
DimPlot(countexp.ST)
head(countexp.ST@assays$METABOLISM$score)[1:5,1:5]
head(t(countexp.ST@assays$METABOLISM$score))[1:5,1:5]
colnames(countexp.ST@assays$METABOLISM$score) <- row.names(countexp.ST@meta.data)

table(countexp.ST$spot_type)

countexp.ST

# 5.Subset the PT (PT-S1,PT-S2,PT-S3,PT-injured) spots

countexp.ST2 <- subset(countexp.ST,subset = spot_type %in% c("PT-S1","PT-S2","PT-S3","PT-injured"))
countexp.ST2 <- AddMetaData(countexp.ST2,t(countexp.ST@assays$METABOLISM$score)[,c(1,2,15,18,20,33,34,36,42,49,83,84)],
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
