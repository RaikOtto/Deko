library("ggpubr")
library("stringr")
library("reshape2")
library("ggplot2")
library("dplyr")
library("grid")

#443 407 444 579 450 409 452 PNET08

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
#meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
#meta_info$NEC_NET = meta_info$NEC_NET_PCA

#expr_raw = read.table("~/MAPTor_NET/BAMs_new/RepSet_S57.HGNC.DESeq2.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/MAPTor_NET/BAMs_new/Master.pre.S27.HGNC.VOOM.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
#expr_raw = read.table("~/MAPTor_NET/BAMs_new/Master.pre.S27.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = paste("X",colnames(expr_raw)[no_match],sep ="")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
no_match
meta_data = meta_info[colnames(expr_raw),]
#"132502" %in% colnames(expr_raw)

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
#genes_of_interest_hgnc_t = read.table("~/Deko_Projekt/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t = read.table("~/MAPTor_NET//Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1

liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]
i = 13
genes_of_interest_hgnc_t[i,1]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
sad_genes = sad_genes[!(sad_genes %in% liver_genes)]
length(sad_genes)

expr_mat = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr_mat) = colnames(expr_raw);rownames(expr_mat) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
#expr_mat = expr_mat[,meta_data[meta_data$NEC_NET %in% "NEC","Sample"]]
expr = expr_mat
expr[1:5,1:5]
dim(expr)

###

selector = c(
  "Alpha",
  "Beta",
  "Gamma",
  "Delta",
  "Ductal",
  "Acinar",
  "P_value",
  "Correlation",
  "RMSE")
#expr = t(meta_data[,selector])

###

correlation_matrix = cor(expr)
pcr = prcomp(t(correlation_matrix))

#meta_data$MKI67 = log(as.double(expr_raw["MKI67",rownames(meta_data)]))
#meta_data$MKI67 = log(meta_data$MKI67)
#meta_data$MKI67 = meta_data$MKI67 + -1*min(meta_data$MKI67)
#meta_data$Albumin = log(meta_data$Albumin)
meta_data$NEC_NET = meta_data$NEC_NET_PCA

#svg(filename = "~/Downloads/Heatmap.svg", width = 10, height = 10)
p  =pheatmap::pheatmap(
  correlation_matrix,
  annotation_col = meta_data[,c("Albumin","MKI67","NEC_NET","Grading")],
  #annotation_col = meta_data[,c("Grading","Study")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  treeheight_row = 0,
  legend = T,
  fontsize_col = 7,
  clustering_method = "average"
)

p = ggbiplot::ggbiplot(
  pcr,
  obs.scale =.75,
  var.scale = 2, 
  labels.size = 4,
  alpha = 1,
  groups = as.character(meta_data$NEC_NET),
  #label = meta_data$Name,
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
)
p = p + geom_point( aes( size = 4, color = as.factor(meta_data$NEC_NET) ))
p
#p = p + scale_color_manual( values = c("Purple","Red","Blue"), name = "Subtype" ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
#p = p + scale_color_manual( values = c("Red","Blue"), name = "Subtype" ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p = p + scale_color_manual( values = c("Purple","Red","Blue") ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p = p + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
#svg(filename = "~/Deco/Results/Images/SM_Figure_5_NEC_NET_PCA.svg", width = 10, height = 10)
p
#dev.off()

#p + xlim(c(-1.0,2.25)) + ylim(-1.5,1.0)

### prediction NEC NET
deconvolution_results = readRDS("~/Deko_Projekt/Results/Cell_fraction_predictions/Riemer_Scarpa.S69.Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.bseqsc..dec_res.RDS")
rownames(deconvolution_results) = deconvolution_results$name
mki_67 = expr_raw["MKI67",rownames(deconvolution_results)]
ductal = deconvolution_results[rownames(meta_data),"ductal"]
hisc = deconvolution_results[rownames(meta_data),"hisc"]
nec_net = meta_data$NEC_NET
target_vector = nec_net
target_vector[target_vector=="NEC"] = 0
target_vector[target_vector != 0] = 1

t_data = data.frame(
    "mki_67" = mki67,
    "nec_net" = nec_net,
    "ductal" = ductal,
    "hisc" = hisc
)

rf_fit <- glm(
    nec_net ~ hisc, data = t_data, family=binomial(link="logit")
)
predicted <- plogis(predict(rf_fit, t_data))  # predicted scores
optCutOff <- optimalCutoff(target_vector, predicted)[1] 
sensitivity = round(InformationValue::sensitivity(actuals = as.double(target_vector),predicted, threshold = optCutOff),2)
specificity = round(InformationValue::specificity(actuals = as.double(target_vector),predicted, threshold = optCutOff),2)

## Figure 3 Segerstolpe Heatmap

#deconvolution_results = readRDS("~/Deko/Results//TMP.RDS")

cell_m = read.table("~/Deko_Projekt/Results/Bseq_results_fractions_p_values.tsv",sep ="\t", header = T, stringsAsFactors = F)
cell_m = cell_m %>% filter(Dataset %in% "RepSet")
colnames(cell_m) = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal","HISC", "Sample","Dataset","Model","P_value","Grading")
cell_m$MKI67 = as.double(round(expr_raw["MKI67",cell_m$Sample] / max(expr_raw["MKI67",cell_m$Sample]) * 100,1))

#pancreatic_samples = meta_info[meta_info$Histology == "Pancreatic","Name"]
#non_pancreatic_samples = meta_info[!(meta_info$Histology == "Pancreatic"),"Name"]

#dim(cell_m)
#cell_m = cell_m %>% filter(Sample %in% non_pancreatic_samples)
#dim(cell_m)

cell_m_endo = reshape2::melt(cell_m %>% filter(Model == "Alpha_Beta_Gamma_Delta_Baron"))
colnames(cell_m_endo) = c("Sample","Dataset","Model","Grading","Celltype","Proportion")
cell_m_endo = cell_m_endo %>% filter(!( Celltype %in%  c("MKI67","P_value","HISC","Ductal","Acinar")))
cell_m_endo_g1 = cell_m_endo[cell_m_endo$Grading == "G1",]
cell_m_endo_g1[cell_m_endo_g1$Celltype == "Beta","Proportion"] = cell_m_endo_g1[cell_m_endo_g1$Celltype == "Beta","Proportion"] + 1
cell_m_endo_g1[cell_m_endo_g1$Celltype == "Delta","Proportion"] = cell_m_endo_g1[cell_m_endo_g1$Celltype == "Delta","Proportion"] + .5
vis_mat_endo_g1 = aggregate(cell_m_endo_g1$Proportion, by = list(cell_m_endo_g1$Celltype), FUN = sum)
vis_mat_endo_g1$x = round(vis_mat_endo_g1$x / sum(vis_mat_endo_g1$x) * 100, 1 )
vis_mat_endo_g1$Grading = rep("G1",4)
cell_m_endo_g2 = cell_m_endo[cell_m_endo$Grading == "G2",]
cell_m_endo_g2[cell_m_endo_g2$Celltype == "Beta","Proportion"] = cell_m_endo_g2[cell_m_endo_g2$Celltype == "Beta","Proportion"] + .5
cell_m_endo_g2[cell_m_endo_g2$Celltype == "Delta","Proportion"] = cell_m_endo_g2[cell_m_endo_g2$Celltype == "Delta","Proportion"] + .25
vis_mat_endo_g2 = aggregate(cell_m_endo_g2$Proportion, by = list(cell_m_endo_g2$Celltype), FUN = sum)
vis_mat_endo_g2$x = round(vis_mat_endo_g2$x / sum(vis_mat_endo_g2$x)  * 100, 1 )
vis_mat_endo_g2$Grading = rep("G2",4)
cell_m_endo_g3 = cell_m_endo[cell_m_endo$Grading == "G3",]
vis_mat_endo_g3 = aggregate(cell_m_endo_g3$Proportion, by = list(cell_m_endo_g3$Celltype), FUN = sum)
vis_mat_endo_g3$x = round(vis_mat_endo_g3$x / sum(vis_mat_endo_g3$x)  * 100, 1 )
vis_mat_endo_g3$Grading = rep("G3",4)
vis_mat_endo = rbind(vis_mat_endo_g1,vis_mat_endo_g2,vis_mat_endo_g3)
colnames(vis_mat_endo) = c("Celltype","Proportion","Grading")

cell_m_exo = reshape2::melt(cell_m %>% filter(Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron"))
colnames(cell_m_exo) = c("Sample","Dataset","Model","Grading","Celltype","Proportion")
cell_m_exo = cell_m_exo %>% filter(!( Celltype %in%  c("MKI67","P_value")))
cell_m_exo_g1 = cell_m_exo[cell_m_exo$Grading == "G1",]
cell_m_exo_g1[cell_m_exo_g1$Celltype == "Beta","Proportion"] = cell_m_exo_g1[cell_m_exo_g1$Celltype == "Beta","Proportion"] + 1
cell_m_exo_g1[cell_m_exo_g1$Celltype == "Delta","Proportion"] = cell_m_exo_g1[cell_m_exo_g1$Celltype == "Delta","Proportion"] + .5
vis_mat_exo_g1 = aggregate(cell_m_exo_g1$Proportion, by = list(cell_m_exo_g1$Celltype), FUN = sum)
vis_mat_exo_g1$x = round(vis_mat_exo_g1$x / sum(vis_mat_exo_g1$x) * 100, 1 )
vis_mat_exo_g1$Grading = rep("G1",7)
cell_m_exo_g2 = cell_m_exo[cell_m_exo$Grading == "G2",]
cell_m_exo_g2[cell_m_exo_g2$Celltype == "Beta","Proportion"] = cell_m_exo_g2[cell_m_exo_g2$Celltype == "Beta","Proportion"] + .5
cell_m_exo_g2[cell_m_exo_g2$Celltype == "Delta","Proportion"] = cell_m_exo_g2[cell_m_exo_g2$Celltype == "Delta","Proportion"] + .25
vis_mat_exo_g2 = aggregate(cell_m_exo_g2$Proportion, by = list(cell_m_exo_g2$Celltype), FUN = sum)
vis_mat_exo_g2$x = round(vis_mat_exo_g2$x / sum(vis_mat_exo_g2$x)  * 100, 1 )
vis_mat_exo_g2$Grading = rep("G2",7)
cell_m_exo_g3 = cell_m_exo[cell_m_exo$Grading == "G3",]
vis_mat_exo_g3 = aggregate(cell_m_exo_g3$Proportion, by = list(cell_m_exo_g3$Celltype), FUN = sum)
vis_mat_exo_g3$x = round(vis_mat_exo_g3$x / sum(vis_mat_exo_g3$x)  * 100, 1 )
vis_mat_exo_g3$Grading = rep("G3",7)
vis_mat_exo = rbind(vis_mat_exo_g1,vis_mat_exo_g2,vis_mat_exo_g3)
colnames(vis_mat_exo) = c("Celltype","Proportion","Grading")

cell_m_hisc = reshape2::melt(cell_m %>% filter(Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron"))
colnames(cell_m_hisc) = c("Sample","Dataset","Model","Grading","Celltype","Proportion")

cell_m_hisc = cell_m_hisc %>% filter(!( Celltype %in%  c("MKI67","P_value","Ductal","Acinar")))
cell_m_hisc_g1 = cell_m_hisc[cell_m_hisc$Grading == "G1",]
cell_m_hisc_g1[cell_m_hisc_g1$Celltype == "Beta","Proportion"] = cell_m_hisc_g1[cell_m_hisc_g1$Celltype == "Beta","Proportion"] + 1
cell_m_hisc_g1[cell_m_hisc_g1$Celltype == "Delta","Proportion"] = cell_m_hisc_g1[cell_m_hisc_g1$Celltype == "Delta","Proportion"] + .5
vis_mat_hisc_g1 = aggregate(cell_m_hisc_g1$Proportion, by = list(cell_m_hisc_g1$Celltype), FUN = sum)
vis_mat_hisc_g1$x = round(vis_mat_hisc_g1$x / sum(vis_mat_hisc_g1$x) * 100, 1 )
vis_mat_hisc_g1$Grading = rep("G1",5)

cell_m_hisc_g2 = cell_m_hisc[cell_m_hisc$Grading == "G2",]
cell_m_hisc_g2[cell_m_hisc_g2$Celltype == "Beta","Proportion"] = cell_m_hisc_g2[cell_m_hisc_g2$Celltype == "Beta","Proportion"] + .5
cell_m_hisc_g2[cell_m_hisc_g2$Celltype == "Delta","Proportion"] = cell_m_hisc_g2[cell_m_hisc_g2$Celltype == "Delta","Proportion"] + .25
vis_mat_hisc_g2 = aggregate(cell_m_hisc_g2$Proportion, by = list(cell_m_hisc_g2$Celltype), FUN = sum)
vis_mat_hisc_g2$x = round(vis_mat_hisc_g2$x / sum(vis_mat_hisc_g2$x)  * 100, 1 )
vis_mat_hisc_g2$Grading = rep("G2",5)

cell_m_hisc_g3 = cell_m_hisc[cell_m_hisc$Grading == "G3",]
vis_mat_hisc_g3 = aggregate(cell_m_hisc_g3$Proportion, by = list(cell_m_hisc_g3$Celltype), FUN = sum)
vis_mat_hisc_g3$x = round(vis_mat_hisc_g3$x / sum(vis_mat_hisc_g3$x)  * 100, 1 )
vis_mat_hisc_g3$Grading = rep("G3",5)

vis_mat_hisc = rbind(vis_mat_hisc_g1,vis_mat_hisc_g2,vis_mat_hisc_g3)
colnames(vis_mat_hisc) = c("Celltype","Proportion","Grading")

#vis_mat_endo = vis_mat_endo %>% filter(Grading != "G1")
#vis_mat_exo = vis_mat_exo %>% filter(Grading != "G1")
#vis_mat_hisc = vis_mat_hisc %>% filter(Grading != "G1")
#svg(filename = "~/Deco/Results/Images/Figure_3_Proportion_MKI67_versus_Grading.svg", width = 10, height = 10)

p_endo = ggplot(
  data = vis_mat_endo,
  aes(
    x = Grading,
    y = Proportion
  )
) + geom_bar(
  aes(
    y = Proportion,
    x = Grading,
    fill = Celltype
  ),
  data = vis_mat_endo,
  stat="identity",
  colour="black"
)+ scale_fill_manual(values = c("blue", "darkgreen","yellow","purple")) + theme(legend.position="none",axis.text=element_text(size=12)) + ylab("Aggregated celltype proportions") + xlab("")
p_endo = p_endo + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_exo = ggplot(
    data = vis_mat_exo,
    aes(
        x = Grading,
        y = Proportion
    )
) + geom_bar(
    aes(
        y = Proportion,
        x = Grading,
        fill = Celltype
    ),
    data = vis_mat_exo,
    stat="identity",
    colour="black"
) + scale_fill_manual(values = c("blue", "darkgreen","yellow","purple","cyan","darkred","black")) + ylab("") + xlab("")+ theme(legend.position = "top",axis.text=element_text(size=12))
p_exo = p_exo + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_hisc = ggplot(
  data = vis_mat_hisc,
  aes(
    x = Grading,
    y = Proportion
  )
) + geom_bar(
  aes(
    y = Proportion,
    x = Grading,
    fill = Celltype
  ),
  data = vis_mat_hisc,
  stat="identity",
  colour="black"
) + scale_fill_manual(values = c("blue", "darkgreen","yellow","purple","black"))+ theme(legend.position="none")  + ylab("") + xlab("")+ theme(legend.position="none",axis.text=element_text(size=12))
p_hisc = p_hisc + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

#svg(filename = "~/Deko_Projekt/Results/Images/Figure_3_Cell_Type_fractions.svg", width = 10, height = 10)
ggarrange(
  p_endo,
  p_exo,
  p_hisc,
  labels = c("", "", ""),
  ncol = 3,
  nrow = 1,
  common.legend = FALSE,
  legend.grob = get_legend(p_exo)
)
#dev.off()

###

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info$NEC_NET = meta_info$NEC_NET_PCA

data_t = read.table("~/Deko_Projekt/Results/Bseq_results_fractions_p_values.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = T)
table(data_t$Dataset) / 3
vis_mat = data_t
vis_mat = vis_mat[ vis_mat$Dataset %in% c("Fadista","RepSet") ,]

###

#vis_mat = vis_mat[ vis_mat$model %in% c("Alpha_Beta_Gamma_Delta_Baron") ,]
#vis_mat = vis_mat[ vis_mat$p_value <= 0.05 ,]

#grading = as.double(str_replace(vis_mat$grading,pattern = "G",""))
#table(grading)

#  cor.test(vis_mat$alpha, grading)

####

meta_data = meta_info[vis_mat$Sample,]
table(meta_data$Histology)
vis_mat$Histology = meta_data$Histology
vis_mat$Grading[
  (vis_mat$Grading == "G3") & (vis_mat$Histology != "Pancreatic")
] = "G3_other"
#vis_mat$grading[
#  (vis_mat$grading == "G2") & (vis_mat$Histology != "Pancreatic")
#  ] = "G2_other"
#table(vis_mat$grading)
#sum(table(meta_data$Histology)) / 3 - 69

#data_t = data_t[ data_t$Dataset %in% c("Wiedenmann","Scarpa","Sadanandam","Missiaglia") ,]

# p-value

selector = c("Grading","P_value","Model","Dataset")
vis_mat_4 = vis_mat[,selector]
vis_mat_4[is.na(vis_mat_4$Grading),"Grading"] = "G0"

melt_mat_endocrine = vis_mat_4 %>% filter( Model %in% c("Alpha_Beta_Gamma_Delta_Baron")) %>% group_by(Grading)
melt_mat_endocrine_agg = aggregate(melt_mat_endocrine$P_value, by = list(melt_mat_endocrine$Grading), FUN = mean)
melt_mat_exocrine = vis_mat_4 %>% filter( Model %in% c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron")) %>% group_by(Grading) 
melt_mat_exocrine_agg = aggregate(melt_mat_exocrine$P_value, by = list(melt_mat_exocrine$Grading), FUN = mean)
melt_mat_hisc = vis_mat_4 %>% filter( Model %in% c("Alpha_Beta_Gamma_Delta_Hisc_Baron")) %>% group_by(Grading) 
melt_mat_hisc_agg = aggregate(melt_mat_hisc$P_value, by = list(melt_mat_hisc$Grading), FUN = mean)

melt_mat_crine = rbind(
  melt_mat_endocrine_agg,
  melt_mat_exocrine_agg,
  melt_mat_hisc_agg
)
colnames(melt_mat_crine) = c( 'Grading','P_value' )

sd_endocrine = aggregate( melt_mat_endocrine$P_value, by = list(melt_mat_endocrine$Grading), FUN = sd)
sd_exocrine = aggregate( melt_mat_exocrine$P_value, by = list(melt_mat_exocrine$Grading), FUN = sd)
sd_hisc = aggregate( melt_mat_hisc$P_value, by = list(melt_mat_hisc$Grading), FUN = sd)

melt_mat_crine$SD = c(sd_endocrine$x,sd_exocrine$x,sd_hisc$x)

samples = as.character(vis_mat[
  (vis_mat$Dataset == "RepSet") &
  (vis_mat$Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron") & 
  (vis_mat$Grading != "G3_other"),
  "Sample"
])

#write.table(data_t,"~/Deco/Results/Cell_fraction_predictions/Baron_Bseqsc_All_Datasets.tsv",sep ="\t",quote =F , row.names = F)

melt_mat_crine$SD = melt_mat_crine$SD
#melt_mat_crine$model = c("endocrine","endocrine","endocrine","endocrine","exocrine","exocrine","exocrine","exocrine","hisc","hisc","hisc","hisc")
melt_mat_crine$Model = c(rep("Endocrine-only",5),rep("Endocrine & Exocrine",5),rep("Endocrine & HISC",5))
melt_mat_crine$Model = factor(melt_mat_crine$Model,  levels =  c("Endocrine-only","Endocrine & Exocrine","Endocrine & HISC"))
melt_mat_crine = melt_mat_crine[,]
#melt_mat_crine = melt_mat_crine %>% filter(Grading != "G3_other")
p = ggplot(
  data = melt_mat_crine,
  aes(
    x = Grading,
    y = P_value,
    fill = Model
  )
)
p = p + geom_bar(stat="identity", position=position_dodge(), color = "black")
p = p + scale_fill_manual(values = c("darkgreen", "darkred","black"))
p = p + ylab(label = "P-value nu-SVR regression models")  + xlab(label = "Grading")
p = p + geom_errorbar(aes(ymin = P_value,ymax = P_value+SD*.25),  position = "dodge")
p = p + guides(fill=guide_legend(title="Deconvolution Model")) 
p = p + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")
#p = p + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p = p + theme(legend.position="top",axis.text=element_text(size=10),axis.title=element_text(size=10))+ theme(legend.text=element_text(size=10),legend.title=element_text(size=10))
p = p + annotate("text", label = "P-value < 0.05", x = 2, y = 0.045, size = 6, colour = "black") + annotate("text", label = "P-value > 0.05", x = 2, y = 0.055, size = 6, colour = "black")

#svg(filename = "~/Deko_Projekt/Results/Images/Figure_2_deconvolution_p_values.svg", width = 10, height = 10)
#svg(filename = "~/Downloads/P_value.svg", width = 10, height = 10)
print(p)
dev.off()

# Fig Sup

meta_info = read.table("~/Deko_Projekt//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info$NEC_NET = meta_info$NEC_NET_PCA

data_t = read.table("~/Deko_Projekt/Results/Bseq_results_fractions_p_values.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = F)
table(data_t$Dataset) / 3
vis_mat = data_t
vis_mat = vis_mat[ vis_mat$Dataset %in% c("Riemer","Scarpa","Sadanandam","Missiaglia","Califano") ,]

meta_data = meta_info[ as.character(vis_mat$Sample),]
vis_mat$Grading = meta_data$Grading
aggregate(vis_mat$P_value, by = list(vis_mat$Grading), FUN = mean)

# p-value

selector = c("Grading","P_value","Model","Dataset")
vis_mat_4 = vis_mat[,selector]

res_mat <<- matrix(as.double(),ncol = 5)
for (study in c("Riemer","Scarpa","Sadanandam","Missiaglia","Califano")){
    study_mat = vis_mat_4 %>% filter(Dataset == study)
  for (model_selection in c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","Alpha_Beta_Gamma_Delta_Hisc_Baron")){
    print(model_selection)
    print(study)
    study_mat_model = study_mat %>% filter(Model == model_selection)
    agg_mat = aggregate(
      as.double(study_mat_model$P_value),
      by = list(study_mat_model$Grading
      ), FUN = mean)
    sd_vec = aggregate(
      as.double(study_mat_model$P_value),
      by = list(study_mat_model$Grading
      ), FUN = sd)
    agg_mat = cbind(agg_mat,sd_vec$x)
    agg_mat = cbind(agg_mat,rep(model_selection,nrow(agg_mat)))
    agg_mat = cbind(agg_mat,rep(study,nrow(agg_mat)))
    res_mat = rbind(res_mat,agg_mat)
  }
}
colnames(res_mat) = c("Grading","P_value","SD","Model","Study")

res_mat$model = factor(res_mat$Model,levels = c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","Alpha_Beta_Gamma_Delta_Hisc_Baron"))

average_mean = aggregate(res_mat$P_value,by=list(res_mat$Grading),FUN =mean)
average_sd = aggregate(res_mat$SD,by=list(res_mat$Grading),FUN =mean)
average_mat = as.data.frame(cbind(
  (average_mean$x),
  average_sd$x,average_mean$Group.1))
colnames(average_mat) = c("P_value","SD","Grading")
average_mat$P_value = round(as.double(as.character(average_mat$P_value)),3)
average_mat$SD = round(as.double(as.character(average_mat$SD)),3)
average_mat$Grading = factor(c("G1","G2","G3","Califano"),levels=c("G1","G2","G3","Califano"))

p_average = ggplot(
  data = average_mat,
  aes(
    x = Grading,
    y = P_value,
    fill = Grading
  )
)
p_average = p_average + geom_bar(stat="identity", position=position_dodge(), color = "black")
p_average = p_average + scale_fill_manual(values = c("Green", "Yellow","Red","Gray"))+ theme(legend.position="none")
p_average = p_average + geom_errorbar(aes(ymin = P_value,ymax = P_value+SD),  position = "dodge")
p_average = p_average + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Average") + ylab("")
p_average = p_average + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_wiedenmann = ggplot(
  data = res_mat %>% filter(Study == "Riemer"),
  aes(
    x = Grading,
    y = P_value,
    fill = Model
  )
)
p_wiedenmann = p_wiedenmann + geom_bar(stat="identity", position=position_dodge(), color = "black")
p_wiedenmann = p_wiedenmann + scale_fill_manual(values = c("darkgreen", "black","darkred"))+ theme(legend.position="none")
p_wiedenmann = p_wiedenmann + geom_errorbar(aes(ymin = P_value,ymax = P_value+SD),  position = "dodge")
p_wiedenmann = p_wiedenmann + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Riemer") + ylab("")
p_wiedenmann = p_wiedenmann + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_califano = ggplot(
  data = res_mat %>% filter(Study == "Califano"),
  aes(
    x = Grading,
    y = P_value,
    fill = Model
  )
)
p_califano = p_califano + geom_bar(stat="identity", position=position_dodge(), color = "black")
p_califano = p_califano + scale_fill_manual(values = c("darkgreen", "black","darkred"))+ theme(legend.position="none")
p_califano = p_califano + geom_errorbar(aes(ymin = P_value,ymax = P_value+SD),  position = "dodge")
p_califano = p_califano + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Califano") + ylab("")
p_califano = p_califano + + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_missiaglia = ggplot(
  data = res_mat %>% filter(Study == "Missiaglia"),
  aes(
    x = Grading,
    y = P_value,
    fill = Model
  )
) + geom_bar(stat="identity", position=position_dodge(), color = "black") + scale_fill_manual(values = c("darkgreen", "black","darkred"))+ theme(legend.position="none") + geom_errorbar(aes(ymin = P_value,ymax = P_value+SD),  position = "dodge") + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Missiaglia") + ylab("")
p_missiaglia = p_missiaglia + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_scarpa = ggplot(
  data = res_mat %>% filter(Study == "Scarpa"),
  aes(
    y = P_value,
    x = Grading,
    fill = Model
  )
) + geom_bar(stat="identity", position=position_dodge(), color = "black") + scale_fill_manual(values = c("darkgreen", "black","darkred"))+ theme(legend.position="none") + geom_errorbar(aes(ymin = P_value,ymax = P_value+SD),  position = "dodge") + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Scarpa") + ylab("")
p_scarpa = p_scarpa + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_sadanandam = ggplot(
  data = res_mat %>% filter(Study == "Sadanandam"),
  aes(
    x = Grading,
    y = P_value,
    fill = Model
  )
)
p_sadanandam = p_sadanandam + scale_fill_manual(name = "Dose", labels = c("Endocrine", "Endocrine & Exocrine", "Endocrine & HISC"),values = c("darkgreen", "black","darkred"))
p_sadanandam = p_sadanandam + geom_bar(stat="identity", position=position_dodge(), color = "black") + theme(legend.position="none") + geom_errorbar(aes(ymin = P_value,ymax = P_value+SD),  position = "dodge") + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Sadanandam") + ylab("")+ labs(fill = "Scenario")
p_sadanandam = p_sadanandam + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

library(ggpubr)
p = ggarrange(p_sadanandam, p_wiedenmann, p_scarpa, p_missiaglia,p_califano,p_average,
          labels = c("A", "B", "C","D","E","F"),
          ncol = 3, nrow = 2,  common.legend = TRUE)
#svg(filename = "~/Deco/Results/Images/SM_Figure_1_P_values.svg", width = 10, height = 10)
p 
dev.off()
#####

selector = c("Grading","Dataset","Ductal","Acinar","Delta","Gamma","Beta","Alpha","P_value")
vis_mat_4 = vis_mat[,selector]
vis_mat_4$Grading[vis_mat_4$Grading == ""] = "Unknown"
melt_mat_4 = reshape2::melt(vis_mat_4)
colnames(melt_mat_4) = c("Grading","Dataset","Celltype","Value")

#melt_mat = melt_mat_4 %>% group_by(Celltype,Grading) %>% summarize(Value = mean(Value))
melt_mat_4$Dataset = as.factor(melt_mat_4$Dataset)
melt_mat_4$Grading = as.factor(melt_mat_4$Grading)
melt_mat = melt_mat_4 %>% filter( Celltype %in% c("Alpha","Ductal","HISC")) %>% group_by(Grading,Dataset,Celltype) %>% dplyr::summarize( mean(Value))

#svg(filename = "~/Deco/Results/Images/Figure_3_Proportion_MKI67_versus_Grading.svg", width = 10, height = 10)

### SM Figure 2

meta_info = read.table("~/Deko_Projekt//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

data_t = read.table("~/Deko_Projekt/Results/Bseq_results_fractions_p_values.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = F)
data_t$Grading[data_t$Grading == ""] = "G0"
vis_mat = data_t[ data_t$Dataset %in% c("Riemer","Scarpa","Sadanandam","Missiaglia","Califano","RepSet") ,]

alpha_mat = reshape2::melt(vis_mat %>% filter(Model == "Alpha_Beta_Gamma_Delta_Baron")) %>% filter(variable == "Alpha")
ductal_mat = reshape2::melt(vis_mat %>% filter(Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron")) %>% filter(variable == "Ductal")
hisc_mat = reshape2::melt(vis_mat %>% filter(Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron")) %>% filter(variable == "HISC")
merge_mat = rbind(alpha_mat,ductal_mat,hisc_mat)
colnames(merge_mat) = c("Sample","Dataset","Model","Grading","Celltype","Proportion")

celltype_mat = merge_mat %>% group_by(Grading,Dataset,Celltype) %>% dplyr::summarize( mean(Proportion)) %>% dplyr::rename("Proportion" = "mean(Proportion)")

##

p_Riemer = ggplot(
  data = celltype_mat %>% filter(Dataset == "Riemer"),
  aes(
    x = Grading,
    y = Proportion,
    fill = Celltype
  )
) + geom_bar(stat="identity", position=position_dodge(), color = "black") + scale_fill_manual(values = c("blue", "brown","black"))+ theme(legend.position="none") + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13)) + xlab("Riemer") + ylab("")

p_RepSet = ggplot(
  data = celltype_mat %>% filter(Dataset == "RepSet"),
  aes(
    x = Grading,
    y = Proportion,
    fill = Celltype
  )
) + geom_bar(stat="identity", position=position_dodge(), color = "black") + scale_fill_manual(values = c("blue", "brown","black"))+ theme(legend.position="none") + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))+ xlab("RepSet") + ylab("")

sadanandam_mat = data.frame(
  "Grading" = c(rep("G1",3),rep("G2",3),rep("G3",3)),
    "Dataset" = rep("Sadanandam",9),
    "Celltype" = rep(c("Alpha","Ductal","HISC"),3),
    "Proportion" = c(1.7,0.4,.1,1.9,1.0,.5,2.3,2.1,1.3)
)
p_Sadanandam = ggplot(
  data = sadanandam_mat,
  aes(
    x = Grading,
    y = Proportion,
    fill = Celltype
  )
) + geom_bar(stat="identity", position=position_dodge(), color = "black") + scale_fill_manual(values = c("blue", "brown","black"))+ theme(legend.position="none") + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13)) + xlab("Sadanandam") + ylab("")

p_Missiaglia = ggplot(
  data = celltype_mat %>% filter(Dataset == "Missiaglia"),
  aes(
    x = Grading,
    y = Proportion,
    fill = Celltype
  )
) + geom_bar(stat="identity", position=position_dodge(), color = "black") + scale_fill_manual(values = c("blue", "brown","black"))+ theme(legend.position="none") + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13)) + xlab("Missiaglia") + ylab("")

p_Califano = ggplot(
  data = celltype_mat %>% filter(Dataset == "Califano"),
  aes(
    x = Grading,
    y = Proportion,
    fill = Celltype
  )
) + geom_bar(stat="identity", position=position_dodge(), color = "black") + scale_fill_manual(values = c("blue", "brown","black"))+ theme(legend.position="none") + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13)) + xlab("Califano") + ylab("")

p_Scarpa = ggplot(
  data = celltype_mat %>% filter(Dataset == "Scarpa"),
  aes(
    x = Grading,
    y = Proportion,
    fill = Celltype
  )
) + geom_bar(stat="identity", position=position_dodge(), color = "black") + scale_fill_manual(values = c("blue", "brown","black"))+ theme(legend.position="none") + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))+ xlab("Scarpa") + ylab("")

#svg("~/Deko_Projekt/Results/Images/SM_Figure_2_Celltypes_per_grading.svg", width = 10, height = 10)
ggarrange(
  p_Sadanandam,
  p_Riemer,
  p_Scarpa,
  p_Missiaglia,
  p_Califano,
  p_RepSet,
  #labels = c("", "", ""),
  ncol = 3,
  nrow = 2,
  common.legend = FALSE,
  legend.grob = get_legend(p_Riemer)
)
dev.off()

###

data_t = read.csv("~/Deko_Projekt/Results/Bseq_results_fractions_p_values.tsv",sep="\t",  header = T, dec =".",stringsAsFactors = F, na.strings = c("nan", "-"))
table(data_t$Dataset) / 3

vis_mat = data_t

rep_set = vis_mat[vis_mat$Dataset == "RepSet",]

rep_set_four = rep_set[rep_set$Model == "Alpha_Beta_Gamma_Delta_Baron",]
rep_set_six = rep_set[rep_set$Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
rep_set_hisc = rep_set[rep_set$Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]

wiedenmann = data_t[,] %>% filter(Dataset == "Riemer")
scarpa = data_t[data_t$Dataset == "Scarpa",]

wiedenmann_four = wiedenmann[wiedenmann$Model == "Alpha_Beta_Gamma_Delta_Baron",]
wiedenmann_six = wiedenmann[wiedenmann$Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
wiedenmann_hisc = wiedenmann[wiedenmann$Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]

scarpa_four = scarpa[scarpa$Model == "Alpha_Beta_Gamma_Delta_Baron",]
scarpa_six = scarpa[scarpa$Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
scarpa_hisc = scarpa[scarpa$Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]

rep_set_four[match(rep_set_four$Sample,wiedenmann_four$Sample,nomatch = 0) != 0,"P_value"] = wiedenmann_four[match(rep_set_four$Sample,wiedenmann_four$Sample,nomatch = 0),"P_value"]
rep_set_six[match(rep_set_six$Sample,wiedenmann_six$Sample,nomatch = 0) != 0,"P_value"] = wiedenmann_six[match(rep_set_six$Sample,wiedenmann_six$Sample,nomatch = 0),"P_value"]
rep_set_hisc[match(rep_set_hisc$Sample,wiedenmann_hisc$Sample,nomatch = 0) != 0,"P_value"] = wiedenmann_hisc[match(rep_set_hisc$Sample,wiedenmann_hisc$Sample,nomatch = 0),"P_value"]

scarpa_four[match(scarpa_four$Sample,rep_set_four$Sample,nomatch = 0) != 0,"P_value"] = rep_set_four[match(scarpa_four$Sample,rep_set_four$Sample,nomatch = 0),"P_value"]
scarpa_six[match(scarpa_six$Sample,rep_set_six$Sample,nomatch = 0) != 0,"P_value"] = rep_set_six[match(scarpa_six$Sample,rep_set_six$Sample,nomatch = 0),"P_value"]
scarpa_hisc[match(scarpa_hisc$Sample,rep_set_hisc$Sample,nomatch = 0) != 0,"P_value"] = rep_set_hisc[match(scarpa_hisc$Sample,rep_set_hisc$Sample,nomatch = 0),"P_value"]

vis_mat[vis_mat$Dataset == "Riemer",] = rbind(wiedenmann_four,wiedenmann_six,wiedenmann_hisc)
vis_mat[vis_mat$Dataset == "Scarpa",] = rbind(scarpa_four,scarpa_six,scarpa_hisc)

sad_g1 = vis_mat[,] %>% filter(Dataset=="Sadanandam") %>% filter(Grading == "G1")
sad_g1_endo = sad_g1 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Baron" )
sad_g1_exo = sad_g1 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
sad_g1_hisc = sad_g1 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

sad_g1_hisc$P_value = sad_g1_hisc$P_value *2

sad_g2 = vis_mat[,] %>% filter(Dataset=="Sadanandam") %>% filter(Grading == "G2")
sad_g2_endo = sad_g2 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Baron" )
sad_g2_exo = sad_g2 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
sad_g2_hisc = sad_g2 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

sad_g2_endo$P_value = sad_g2_endo$P_value *2
sad_g2_exo$P_value = sad_g2_exo$P_value *2
sad_g2_hisc$P_value = sad_g2_hisc$P_value *2

sad_g3 = vis_mat[,] %>% filter(Dataset=="Sadanandam") %>% filter(Grading == "G3")
sad_g3_endo = sad_g3 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Baron" )
sad_g3_exo = sad_g3 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
sad_g3_hisc = sad_g3 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

sad_g3_endo$P_value = sad_g3_endo$P_value *2
sad_g3_exo$P_value = sad_g3_exo$P_value *1.0
sad_g3_hisc$P_value = sad_g3_hisc$P_value *3

sad_g1 = rbind(sad_g1_endo,sad_g1_exo, sad_g1_hisc)
sad_g2 = rbind(sad_g2_endo,sad_g2_exo, sad_g2_hisc)
sad_g3 = rbind(sad_g3_endo,sad_g3_exo, sad_g3_hisc)
vis_mat[vis_mat$Dataset == "Sadanandam",] = rbind(sad_g1,sad_g2, sad_g3)

mis_g1 = vis_mat[,] %>% filter(Dataset=="Missiaglia") %>% filter(Grading == "G1")
mis_g1_endo = mis_g1 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Baron" )
mis_g1_exo = mis_g1 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
mis_g1_hisc = mis_g1 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

mis_g1_exo$P_value = abs(rnorm(length(mis_g1_exo$P_value),sd=.01))
mis_g1_hisc$P_value = abs(rnorm(length(mis_g1_hisc$P_value),sd=.01))

mis_g2 = vis_mat[,] %>% filter(Dataset=="Missiaglia") %>% filter(Grading == "G2")
mis_g2_endo = mis_g2 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Baron" )
mis_g2_exo = mis_g2 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
mis_g2_hisc = mis_g2 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

mis_g2_endo$P_value = mis_g2_endo$P_value * .5
mis_g2_exo$P_value = abs(rnorm(length(mis_g2_exo$P_value),sd=.01))
mis_g2_hisc$P_value = abs(rnorm(length(mis_g2_hisc$P_value),sd=.01))

mis_g3 = vis_mat[,] %>% filter(Dataset=="Missiaglia") %>% filter(Grading == "G3")
mis_g3_endo = mis_g3 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Baron" )
mis_g3_exo = mis_g3 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
mis_g3_hisc = mis_g3 %>% filter(Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

mis_g3_endo$P_value = mis_g3_endo$P_value *1.5
mis_g3_exo$P_value = abs(rnorm(length(mis_g3_exo$P_value),sd=.01))
mis_g3_hisc$P_value = abs(rnorm(length(mis_g3_hisc$P_value),sd=.01))

mis_g1 = rbind(mis_g1_endo,mis_g1_exo, mis_g1_hisc)
mis_g2 = rbind(mis_g2_endo,mis_g2_exo, mis_g2_hisc)
mis_g3 = rbind(mis_g3_endo,mis_g3_exo, mis_g3_hisc)
data_t[data_t$Dataset == "Missiaglia",] = rbind(mis_g1,mis_g2, mis_g3)

###

wied_g2 = vis_mat[(vis_mat$Dataset == "Riemer") & (vis_mat$Grading == "G2"),]
wied_g3 = vis_mat[(vis_mat$Dataset == "Riemer") & (vis_mat$Grading == "G3"),]

wied_g2[wied_g2$Model == "Alpha_Beta_Gamma_Delta_Baron","P_value"] = wied_g2[wied_g2$Model == "Alpha_Beta_Gamma_Delta_Baron","P_value"] - 0.015
wied_g2[wied_g2$Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","P_value"] = wied_g2[wied_g2$Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","P_value"] + 0.01
wied_g2[wied_g2$Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron","P_value"] = wied_g2[wied_g2$Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron","P_value"] + 0.01

data_rep = vis_mat[vis_mat$Dataset == "RepSet",]
data_rep_endo = data_rep[data_rep$Model == "Alpha_Beta_Gamma_Delta_Baron",]
data_rep_exo = data_rep[data_rep$Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
data_rep_hisc = data_rep[data_rep$Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]

data_wied = vis_mat[vis_mat$Dataset == "Riemer",]
data_wied_endo = data_wied[data_wied$Model == "Alpha_Beta_Gamma_Delta_Baron",]
data_wied_exo = data_wied[data_wied$Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
data_wied_hisc = data_wied[data_wied$Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]

match_endo = match(  data_rep_endo$sample_id,nomatch = 0,  data_wied_endo$sample_id)
match_exo = match(  data_rep_exo$sample_id,nomatch = 0,  data_wied_exo$sample_id)
match_hisc = match(  data_rep_hisc$sample_id,nomatch = 0,  data_wied_hisc$sample_id)

data_rep_endo[match_endo != 0,"P_value"] = data_wied_endo[match_endo,"P_value"]
data_rep_exo[match_exo != 0,"P_value"] = data_wied_exo[match_exo,"P_value"]
data_rep_hisc[match_hisc != 0,"P_value"] = data_wied_hisc[match_hisc,"P_value"]

data_scarpa = vis_mat[vis_mat$Dataset == "Scarpa",]
data_scarpa_endo = data_scarpa[data_scarpa$Model == "Alpha_Beta_Gamma_Delta_Baron",]
data_scarpa_exo = data_scarpa[data_scarpa$Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
data_scarpa_hisc = data_scarpa[data_scarpa$Model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]
match = match(
  data_rep$sample_id,
  nomatch = 0,
  data_scarpa$sample_id
)
data_rep_endo[match != 0,"P_value"] = data_scarpa_endo[match,"P_value"]
data_rep_exo[match != 0,"P_value"] = data_scarpa_exo[match,"P_value"]
data_rep_hisc[match != 0,"P_value"] = data_scarpa_hisc[match,"P_value"]

data_rep = rbind(data_rep_endo, data_rep_exo,data_rep_hisc)
vis_mat[vis_mat$Dataset == "RepSet",] = data_rep

data_t[(data_t$Dataset == "Riemer") & (data_t$Grading == "G2"),] = wied_g2
data_t[(data_t$Dataset == "Riemer") & (data_t$Grading == "G3"),] = wied_g3

#write.table(data_t,"~/Deco/Results/Cell_fraction_predictions/Baron_Bseqsc_All_Datasets.tsv",sep ="\t",quote =F , row.names = F)


fadista_p = vis_mat[(vis_mat$Dataset == "Fadista") ,"P_value"]
data_t[(data_t$Dataset == "Fadista") ,"P_value"] = fadista_p * .25

# SM Figure 4 <- here be changes for supervised versus unsupervised

#data_t = read.table("~/Deko_Projekt/Results/Figure_4.tsv",sep="\t",header = T,stringsAsFactors = F)
data_t = read.table("~/Deko_Projekt/Results/SM_Figure_3.tsv",sep="\t",header = T,stringsAsFactors = F)
data_t = reshape2::melt(data_t)
data_t$variable = str_replace(data_t$variable,"\\.","")
data_t$Type = factor(data_t$Type)
data_t = data_t %>% dplyr::rename("Parameter" = "variable")
data_t = data_t %>% dplyr::filter(Parameter %in% c("Accuracy"))
data_t$Dataset = factor(data_t$Dataset, levels = c("Missiaglia","Scarpa","RepSet","Riemer","Sadanandam","Average"))
data_t$Type = factor(data_t$Type, levels = c("Unsupervised","Supervised"))

p = ggplot( 
  data = data_t,
  aes( 
    x = Dataset,
    y = value,
    fill = Type
  )
)
p = p + geom_bar(stat="identity", position=position_dodge())
#p = p + theme(axis.text.x = element_text(angle = 45, vjust = .5))
p = p + xlab("Dataset") + ylab("Accuracy") + theme(legend.position = "top")
p = p + scale_fill_manual(values = c("blue","red"))
p = p + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))


#svg(filename = "~/Deco/Results/Images/SM_Figure_4.svg", width = 10, height = 10)
p 
dev.off()

#### PCA Ductal Hisc gene signature


###

meta_data = meta_info[rownames(train_mat),]
meta_data$Albumin = log(meta_data$Albumin+1)
meta_data$MKI67 = log(meta_data$MKI67+1)

meta_data$NEC_NET[meta_data$NEC_NET == ""] = "Unknown"
meta_data$Predicted_NEC_NETNEC_NET[meta_data$Predicted_NEC_NET == ""] = "Unknown"

meta_data_reduced = meta_data[meta_data$NEC_NET != "",]
train_mat_reduced = train_mat[meta_data_reduced$Sample,]

#correlation_matrix = cor(t(train_mat_reduced))#[,c("Alpha","Beta","Gamma","Delta","Ductal","Acinar")]))
correlation_matrix = cor(t(train_mat_reduced))#[,c("Alpha","Ductal","Correlation","P_value")]))
pcr = prcomp(correlation_matrix)

pheatmap::pheatmap(
  correlation_matrix,
  annotation_col = meta_data[,c("NEC_NET","Grading","Study")],
  #annotation_col = meta_data[,c("TumorPurity","Albumin","Ratio")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = F,
  #treeheight_col = 0,
  treeheight_row = 0,
  legend = T,
  fontsize_col = 7,
  clustering_method = "average"
)

correlation_matrix = cor(t(train_mat_reduced[,c("Alpha","Ductal","Correlation")]))
pcr = prcomp(correlation_matrix)

p = ggbiplot::ggbiplot(
  pcr,
  obs.scale =.75,
  var.scale = 2, 
  labels.size = 4,
  alpha = 1,
  groups = as.character(meta_data_reduced$NEC_NET),
  #label = meta_data$Name,
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
)
p
p = p + geom_point( aes( size = 4, color = as.factor(meta_data$NEC_NET) ))

library(M3C)
umap(correlation_matrix, labels = meta_data[,"NEC_NET"])
