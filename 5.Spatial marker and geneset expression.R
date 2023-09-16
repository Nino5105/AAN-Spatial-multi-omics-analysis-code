# 1.Loading packages

library(SPATA)
library(hdf5r)
library(magrittr)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(gsva)
library(ggsci)
library(homologene)

# 2.Loading datasets

Control_1 <-  initiateSpataObject_10X(input_paths = c("0.Raw data/Control_1/"), sample_names = c("Control_1"))
Control_1

Control_2 <-  initiateSpataObject_10X(input_paths = c("0.Raw data/Control_2/"), sample_names = c("Control_2"))
Control_2

AAN_2W_1 <-  initiateSpataObject_10X(input_paths = c("0.Raw data/AAN_2W_1/"), sample_names = c("AAN_2W_1"))
AAN_2W_1

AAN_2W_2 <-  initiateSpataObject_10X(input_paths = c("0.Raw data/AAN_2W_2/"), sample_names = c("AAN_2W_2"))
AAN_2W_2

AAN_4W_1 <-  initiateSpataObject_10X(input_paths = c("0.Raw data/AAN_4W_1/"), sample_names = c("AAN_4W_1"))
AAN_4W_1

AAN_4W_2 <-  initiateSpataObject_10X(input_paths = c("0.Raw data/AAN_4W_2/"), sample_names = c("AAN_4W_2"))
AAN_4W_2

# 3.Spatial expression level of injury markers

p1 <- plotSurface(object = Control_1,color_to = "Havcr1",
                  pt_size = 1,pt_clrsp = "magma",
                  smooth = T,smooth_span = 0.01,
                  display_title = T)+ ggtitle("Control_1")

p2 <- plotSurface(object = AAN_2W_1,color_to = "Havcr1",
                  pt_size = 1,pt_clrsp = "magma",
                  smooth = T,smooth_span = 0.01,
                  display_title = T) + ggtitle("AAN_2W_1")

p3 <- plotSurface(object = AAN_4W_1,color_to = "Havcr1",
                  pt_size = 1,pt_clrsp = "magma",
                  smooth = T,smooth_span = 0.01,
                  display_title = T) + ggtitle("AAN_4W_1")

p4 <- plotSurface(object = Control_1,color_to = "Timp2",
                  pt_size = 1,pt_clrsp = "magma",
                  smooth = T,smooth_span = 0.01,
                  display_title = T)+ ggtitle("Control_1")

p5 <- plotSurface(object = AAN_2W_1,color_to = "Timp2",
                  pt_size = 1,pt_clrsp = "magma",
                  smooth = T,smooth_span = 0.01,
                  display_title = T) + ggtitle("AAN_2W_1")

p6 <- plotSurface(object = AAN_4W_1,color_to = "Timp2",
                  pt_size = 1,pt_clrsp = "magma",
                  smooth = T,smooth_span = 0.01,
                  display_title = T) + ggtitle("AAN_4W_1")

p7 <- plotSurface(object = Control_1,color_to = "Ifng",
                  pt_size = 1,pt_clrsp = "magma",
                  smooth = T,smooth_span = 0.01,
                  display_title = T)+ ggtitle("Control_1")

p8 <- plotSurface(object = AAN_2W_1,color_to = "Ifng",
                  pt_size = 1,pt_clrsp = "viridis",
                  smooth = T,smooth_span = 0.01,
                  display_title = T) + ggtitle("AAN_2W_1")

p9 <- plotSurface(object = AAN_4W_1,color_to = "Ifng",
                  pt_size = 1,pt_clrsp = "magma",
                  smooth = T,smooth_span = 0.01,
                  display_title = T) + ggtitle("AAN_4W_1")

pdf("marker.pdf",height = 12,width = 12)
(p1 + p2 + p3) / (p4 + p5 + p6) / (p7 + p8 + p9)
dev.off()

# 4. Transformed gene names from mouse to human
 
# Control_1 <- initiateSpataObject_10X(input_paths = c("0.Raw data/Control_1/"), sample_names = c("Control_1"))
Control_1_genename <- Control_1@data@counts@Dimnames[[1]]
Control_1_genename_m2h <- homologene(Control_1_genename, inTax = 10090, outTax = 9606)
Control_1@data@counts <- Control_1@data@counts[Control_1_genename_m2h$`10090`,]
Control_1@data@norm_exp <- Control_1@data@norm_exp[Control_1_genename_m2h$`10090`,]
table(row.names(Control_1@data@counts) == rownames(Control_1@data@norm_exp))
row.names(Control_1@data@counts) <- Control_1_genename_m2h$`9606`
row.names(Control_1@data@norm_exp) <- Control_1_genename_m2h$`9606`

# Control_2 <- initiateSpataObject_10X(input_paths = c("0.Raw data/Control_2/"), sample_names = c("Control_2"))
Control_2_genename <- Control_2@data@counts@Dimnames[[1]]
Control_2_genename_m2h <- homologene(Control_2_genename, inTax = 10090, outTax = 9606)
Control_2@data@counts <- Control_2@data@counts[Control_2_genename_m2h$`10090`,]
Control_2@data@norm_exp <- Control_2@data@norm_exp[Control_2_genename_m2h$`10090`,]
table(row.names(Control_2@data@counts) == rownames(Control_2@data@norm_exp))
row.names(Control_2@data@counts) <- Control_2_genename_m2h$`9606`
row.names(Control_2@data@norm_exp) <- Control_2_genename_m2h$`9606`

# AAN_2W_1 <- initiateSpataObject_10X(input_paths = c("0.Raw data/AAN_2W_1/"), sample_names = c("AAN_2W_1"))
AAN_2W_1_genename <- AAN_2W_1@data@counts@Dimnames[[1]]
AAN_2W_1_genename_m2h <- homologene(AAN_2W_1_genename, inTax = 10090, outTax = 9606)
AAN_2W_1@data@counts <- AAN_2W_1@data@counts[AAN_2W_1_genename_m2h$`10090`,]
AAN_2W_1@data@norm_exp <- AAN_2W_1@data@norm_exp[AAN_2W_1_genename_m2h$`10090`,]
table(row.names(AAN_2W_1@data@counts) == rownames(AAN_2W_1@data@norm_exp))
row.names(AAN_2W_1@data@counts) <- AAN_2W_1_genename_m2h$`9606`
row.names(AAN_2W_1@data@norm_exp) <- AAN_2W_1_genename_m2h$`9606`

# AAN_2W_2 <- initiateSpataObject_10X(input_paths = c("0.Raw data/AAN_2W_2/"), sample_names = c("AAN_2W_2"))
AAN_2W_2_genename <- AAN_2W_2@data@counts@Dimnames[[1]]
AAN_2W_2_genename_m2h <- homologene(AAN_2W_2_genename, inTax = 10090, outTax = 9606)
AAN_2W_2@data@counts <- AAN_2W_2@data@counts[AAN_2W_2_genename_m2h$`10090`,]
AAN_2W_2@data@norm_exp <- AAN_2W_2@data@norm_exp[AAN_2W_2_genename_m2h$`10090`,]
table(row.names(AAN_2W_2@data@counts) == rownames(AAN_2W_2@data@norm_exp))
row.names(AAN_2W_2@data@counts) <- AAN_2W_2_genename_m2h$`9606`
row.names(AAN_2W_2@data@norm_exp) <- AAN_2W_2_genename_m2h$`9606`

# AAN_4W_1 <- initiateSpataObject_10X(input_paths = c("0.Raw data/AAN_4W_1/"), sample_names = c("AAN_4W_1"))
AAN_4W_1_genename <- AAN_4W_1@data@counts@Dimnames[[1]]
AAN_4W_1_genename_m2h <- homologene(AAN_4W_1_genename, inTax = 10090, outTax = 9606)
AAN_4W_1@data@counts <- AAN_4W_1@data@counts[AAN_4W_1_genename_m2h$`10090`,]
AAN_4W_1@data@norm_exp <- AAN_4W_1@data@norm_exp[AAN_4W_1_genename_m2h$`10090`,]
table(row.names(AAN_4W_1@data@counts) == rownames(AAN_4W_1@data@norm_exp))
row.names(AAN_4W_1@data@counts) <- AAN_4W_1_genename_m2h$`9606`
row.names(AAN_4W_1@data@norm_exp) <- AAN_4W_1_genename_m2h$`9606`

# AAN_4W_2 <- initiateSpataObject_10X(input_paths = c("0.Raw data/AAN_4W_2/"), sample_names = c("AAN_4W_2"))
AAN_4W_2_genename <- AAN_4W_2@data@counts@Dimnames[[1]]
AAN_4W_2_genename_m2h <- homologene(AAN_4W_2_genename, inTax = 10090, outTax = 9606)
AAN_4W_2@data@counts <- AAN_4W_2@data@counts[AAN_4W_2_genename_m2h$`10090`,]
AAN_4W_2@data@norm_exp <- AAN_4W_2@data@norm_exp[AAN_4W_2_genename_m2h$`10090`,]
table(row.names(AAN_4W_2@data@counts) == rownames(AAN_4W_2@data@norm_exp))
row.names(AAN_4W_2@data@counts) <- AAN_4W_2_genename_m2h$`9606`
row.names(AAN_4W_2@data@norm_exp) <- AAN_4W_2_genename_m2h$`9606`

# 5. Spatial trajectory of injured pathways

Control_1 <- createTrajectories(object = Control_1)
Control_1@trajectories$Control_1$pathway <- NULL
Control_1_pathway <-  plotTrajectoryGeneSets(
  object = Control_1,
  trajectory_name = "pathway",
  gene_sets = c("HM_APOPTOSIS",
                "HM_INTERFERON_GAMMA_RESPONSE",
                "HM_INFLAMMATORY_RESPONSE",
                "HM_TNFA_SIGNALING_VIA_NFKB",
                "HM_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) + 
  ggplot2::theme(legend.position = "none") + 
  scale_color_npg() # 4*3
Control_1_pathway

Control_2 <- createTrajectories(object = Control_2)
Control_2@trajectories$Control_2$pathway <- NULL
Control_2_pathway <-  plotTrajectoryGeneSets(
  object = Control_2,
  trajectory_name = "pathway",
  gene_sets = c("HM_APOPTOSIS",
                "HM_INTERFERON_GAMMA_RESPONSE",
                "HM_INFLAMMATORY_RESPONSE",
                "HM_TNFA_SIGNALING_VIA_NFKB",
                "HM_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) + 
  ggplot2::theme(legend.position = "none") + 
  scale_color_npg() # 4*3
Control_2_pathway

AAN_2W_1 <- createTrajectories(object = AAN_2W_1)
AAN_2W_1@trajectories$AAN_2W_1$pathway <- NULL
AAN_2W_1_pathway <-  plotTrajectoryGeneSets(
  object = AAN_2W_1,
  trajectory_name = "pathway",
  gene_sets = c("HM_APOPTOSIS",
                "HM_INTERFERON_GAMMA_RESPONSE",
                "HM_INFLAMMATORY_RESPONSE",
                "HM_TNFA_SIGNALING_VIA_NFKB",
                "HM_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) + 
  ggplot2::theme(legend.position = "none") + 
  scale_color_npg() # 4*3
AAN_2W_1_pathway

AAN_2W_2 <- createTrajectories(object = AAN_2W_2)
AAN_2W_2@trajectories$AAN_2W_2$pathway <- NULL
AAN_2W_2_pathway <-  plotTrajectoryGeneSets(
  object = AAN_2W_2,
  trajectory_name = "pathway",
  gene_sets = c("HM_APOPTOSIS",
                "HM_INTERFERON_GAMMA_RESPONSE",
                "HM_INFLAMMATORY_RESPONSE",
                "HM_TNFA_SIGNALING_VIA_NFKB",
                "HM_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) + 
  ggplot2::theme(legend.position = "none") + 
  scale_color_npg() # 4*3
AAN_2W_2_pathway

AAN_4W_1 <- createTrajectories(object = AAN_4W_1)
AAN_4W_1@trajectories$AAN_4W_1$pathway <- NULL
AAN_4W_1_pathway <-  plotTrajectoryGeneSets(
  object = AAN_4W_1,
  trajectory_name = "pathway",
  gene_sets = c("HM_APOPTOSIS",
                "HM_INTERFERON_GAMMA_RESPONSE",
                "HM_INFLAMMATORY_RESPONSE",
                "HM_TNFA_SIGNALING_VIA_NFKB",
                "HM_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) + 
  ggplot2::theme(legend.position = "right") + 
  scale_color_npg() # 4*3
AAN_4W_1_pathway

AAN_4W_2 <- createTrajectories(object = AAN_4W_2)
AAN_4W_2@trajectories$AAN_4W_2$pathway <- NULL
AAN_4W_2_pathway <-  plotTrajectoryGeneSets(
  object = AAN_4W_2,
  trajectory_name = "pathway",
  gene_sets = c("HM_APOPTOSIS",
                "HM_INTERFERON_GAMMA_RESPONSE",
                "HM_INFLAMMATORY_RESPONSE",
                "HM_TNFA_SIGNALING_VIA_NFKB",
                "HM_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) + 
  ggplot2::theme(legend.position = "right") + 
  scale_color_npg() # 4*3
AAN_4W_2_pathway

(Control_1_pathway + ggtitle("Control_1")+ xlab("")) +
  (AAN_2W_1_pathway + ggtitle("AAN_2W_1")+ ylab("")) + 
  (AAN_4W_1_pathway + ggtitle("AAN_4W_1")+ xlab("")+ ylab("")) # 12*3

(Control_2_pathway + ggtitle("Control_2")+ xlab("")) +
  (AAN_2W_2_pathway + ggtitle("AAN_2W_2")+ ylab("")) + 
  (AAN_4W_2_pathway + ggtitle("AAN_4W_2")+ xlab("")+ ylab("")) # 12*3


# 6. Spatial trajectory of meatabolic pathways

Control_1_META_pathway <-  plotTrajectoryGeneSets(
  object = Control_1,
  trajectory_name = "pathway",
  gene_sets = c("HM_GLYCOLYSIS",
                "HM_FATTY_ACID_METABOLISM",
                "HM_OXIDATIVE_PHOSPHORYLATION",
                "HM_XENOBIOTIC_METABOLISM",
                "BP.GO_PURINE_NUCLEOBASE_METABOLIC_PROCESS"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) +  
  ggplot2::theme(legend.position = "none") + 
  scale_color_aaas() # 4*3
Control_1_META_pathway

Control_2_META_pathway <-  plotTrajectoryGeneSets(
  object = Control_2,
  trajectory_name = "pathway",
  gene_sets = c("HM_GLYCOLYSIS",
                "HM_FATTY_ACID_METABOLISM",
                "HM_OXIDATIVE_PHOSPHORYLATION",
                "HM_XENOBIOTIC_METABOLISM",
                "BP.GO_PURINE_NUCLEOBASE_METABOLIC_PROCESS"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) + 
  ggplot2::theme(legend.position = "none") + 
  scale_color_aaas() # 4*3
Control_2_META_pathway

AAN_2W_1_META_pathway <-  plotTrajectoryGeneSets(
  object = AAN_2W_1,
  trajectory_name = "pathway",
  gene_sets = c("HM_GLYCOLYSIS",
                "HM_FATTY_ACID_METABOLISM",
                "HM_OXIDATIVE_PHOSPHORYLATION",
                "HM_XENOBIOTIC_METABOLISM",
                "BP.GO_PURINE_NUCLEOBASE_METABOLIC_PROCESS"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) + 
  ggplot2::theme(legend.position = "none") + 
  scale_color_aaas() # 4*3
AAN_2W_1_META_pathway

AAN_2W_2_META_pathway <-  plotTrajectoryGeneSets(
  object = AAN_2W_2,
  trajectory_name = "pathway",
  gene_sets = c("HM_GLYCOLYSIS",
                "HM_FATTY_ACID_METABOLISM",
                "HM_OXIDATIVE_PHOSPHORYLATION",
                "HM_XENOBIOTIC_METABOLISM",
                "BP.GO_PURINE_NUCLEOBASE_METABOLIC_PROCESS"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) + 
  ggplot2::theme(legend.position = "none") + 
  scale_color_aaas() # 4*3
AAN_2W_2_META_pathway

AAN_4W_1_META_pathway <-  plotTrajectoryGeneSets(
  object = AAN_4W_1,
  trajectory_name = "pathway",
  gene_sets = c("HM_GLYCOLYSIS",
                "HM_FATTY_ACID_METABOLISM",
                "HM_OXIDATIVE_PHOSPHORYLATION",
                "HM_XENOBIOTIC_METABOLISM",
                "BP.GO_PURINE_NUCLEOBASE_METABOLIC_PROCESS"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) + 
  ggplot2::theme(legend.position = "right") + 
  scale_color_aaas() # 4*3
AAN_4W_1_META_pathway

AAN_4W_2_META_pathway <-  plotTrajectoryGeneSets(
  object = AAN_4W_2,
  trajectory_name = "pathway",
  gene_sets = c("HM_GLYCOLYSIS",
                "HM_FATTY_ACID_METABOLISM",
                "HM_OXIDATIVE_PHOSPHORYLATION",
                "HM_XENOBIOTIC_METABOLISM",
                "BP.GO_PURINE_NUCLEOBASE_METABOLIC_PROCESS"),
  display_trajectory_parts = T,method_gs = "mean",smooth_se =T) + 
  ggplot2::theme(legend.position = "right") + 
  scale_color_aaas() # 4*3
AAN_4W_2_META_pathway

(Control_1_META_pathway + ggtitle("Control_1")+ xlab("")) +
  (AAN_2W_1_META_pathway + ggtitle("AAN_2W_1")+ ylab("")) + 
  (AAN_4W_1_META_pathway + ggtitle("AAN_4W_1")+ xlab("")+ ylab("")) # 12*3

(Control_2_META_pathway + ggtitle("Control_2")+ xlab("")) +
  (AAN_2W_2_META_pathway + ggtitle("AAN_2W_2")+ ylab("")) + 
  (AAN_4W_2_META_pathway + ggtitle("AAN_4W_2")+ xlab("")+ ylab("")) # 12*3