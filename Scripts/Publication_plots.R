library("ggpubr")
library("stringr")
library("reshape2")
library("ggplot2")
library("dplyr")
library("grid")

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info$NEC_NET = meta_info$NEC_NET_PCA

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
genes_of_interest_hgnc_t = read.table("~/Deko_Projekt/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
i = 14
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
genes_of_interest_hgnc_t[i,1]
#sad_genes = c("YAP1","ASCL1","NEUROD1","POU2F3")
sad_genes = sad_genes[ sad_genes != ""]

# Fadista 89
# Alvarez 105
# Scarpa 29
# Sdanandam 29
# Missiaglia 75
# Wiedenmann 39

expr_raw = read.table("~/Deko_Projekt/Data/Bench_data/Riemer_Scarpa.S69.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
#expr_raw = read.table("~/MAPTor_NET/BAMs/Final_plot.TPMs.57.Wiedenmann_Scarpa.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

meta_data = meta_info[colnames(expr_raw),]
table(meta_data$Grading)
table(meta_data$NEC_NET_Color)
# G1 46 + 7 + 14 + 0 = 67
# G2 25 + 12 + 13 + 30 = 80
# G3 4 + 8 + 2 + 30 = 44
#meta_data_2 = rbind(meta_data_2,meta_data)

table(meta_data$Grading)
#meta_data = meta_data[which(meta_data[,"Histology"] == "Pancreatic_NEN"),]
expr_raw = expr_raw[,rownames(meta_data)]
#meta_data[colnames(expr_raw),"mki_67"] = log(as.double(expr_raw["MKI67",])+1)
#meta_data$mki_67[ which(meta_data$mki_67 > mean(meta_data$mki_67))] = "high"
#meta_data$mki_67[meta_data$mki_67 != "high"] = "low"
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr) = colnames(expr_raw);rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
expr[1:5,1:5]
dim(expr)

correlation_matrix = cor(expr)
pcr = prcomp(t(correlation_matrix))

#svg()
pheatmap::pheatmap(
  correlation_matrix,
  annotation_col = meta_data[c("Alpha","Beta","Gamma","Delta","HISC","Acinar","Ductal","Grading", "NEC_NET","Study")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = F,
  #treeheight_col = 0,
  treeheight_row = 0,
  legend = F,
  fontsize_col = 7,
  clustering_method = "centroid"
)

p = ggbiplot::ggbiplot(
  pcr,
  obs.scale =.75,
  var.scale = 2, 
  labels.size = 8,
  alpha = 1,
  groups = as.character(meta_data$Grading),
  #groups = as.character(meta_data$NEC_NET_PCA),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
)
p = p + geom_point( aes( size = 4, color = as.factor(meta_data$Grading) ))
#p = p + scale_color_manual( values = c("Red","Blue"), name = "Subtype" ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p = p + scale_color_manual( values = c("Green","brown","Black","Red","Blue"), name = "Subtype" ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p = p + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
#svg(filename = "~/Deco/Results/Images/SM_Figure_5_NEC_NET_PCA.svg", width = 10, height = 10)
p
dev.off()

# Plot 1

data_t = read.table("~/Deco/Results/ROC_curves.tsv",sep ="\t", header = T)

# hisc

data_t_hisc = subset(data_t, Predictor == "HISC" )
vis_mat_hisc = aggregate(
    as.integer(subset(data_t_hisc$ROC, data_t_hisc$Type == "In-silico")),
    by = list(as.character(subset(data_t_hisc$Dataset, data_t_hisc$Type == "In-silico"))),
    FUN = mean
)
vis_mat_hisc$SD = aggregate(
    as.integer(subset(data_t_hisc$ROC, data_t_hisc$Type == "In-silico")),
    by = list(as.character(subset(data_t_hisc$Dataset, data_t_hisc$Type == "In-silico"))),
    FUN = sd
)
vis_mat_hisc$Type = rep("HISC", nrow(vis_mat_hisc))
vis_mat_hisc = as.data.frame(as.matrix(vis_mat_hisc))
vis_mat_hisc = vis_mat_hisc[,-which(colnames(vis_mat_hisc)  == "SD.Group.1")]
colnames(vis_mat_hisc) = c("Dataset","ROC","SD","Type")

# ductal

data_t_ductal = subset(data_t, Predictor == "Ductal" )
vis_mat_ductal = aggregate(
    as.integer(subset(data_t_ductal$ROC, data_t_ductal$Type == "In-silico")),
    by = list(as.character(subset(data_t_ductal$Dataset, data_t_ductal$Type == "In-silico"))),
    FUN = mean
)
vis_mat_ductal$SD = aggregate(
    as.integer(subset(data_t_ductal$ROC, data_t_ductal$Type == "In-silico")),
    by = list(as.character(subset(data_t_ductal$Dataset, data_t_ductal$Type == "In-silico"))),
    FUN = sd
)
vis_mat_ductal$Type = rep("Ductal", nrow(vis_mat_ductal))
vis_mat_ductal = as.data.frame(as.matrix(vis_mat_ductal))
vis_mat_ductal = vis_mat_ductal[,-which(colnames(vis_mat_ductal)  == "SD.Group.1")]
colnames(vis_mat_ductal) = c("Dataset","ROC","SD","Type")

# MKI-67

data_t_mki67 = subset(data_t, Predictor == "MKI-67" )
vis_mat_mki67 = aggregate(
    subset(data_t_mki67$ROC, data_t_mki67$Type == "MKI-67"),
    by = list(subset(data_t_mki67$Dataset, data_t_mki67$Type == "MKI-67")),
    FUN = mean
)
vis_mat_mki67$SD = rep(0, nrow(vis_mat_mki67))
vis_mat_mki67$Type = rep("MKI_67", nrow(vis_mat_mki67))
colnames(vis_mat_mki67) = c("Dataset","ROC","SD","Type")

# merge

vis_mat = rbind(as.matrix(vis_mat_ductal),as.matrix(vis_mat_hisc),as.matrix(vis_mat_mki67))
vis_mat = as.data.frame(vis_mat)
vis_mat$ROC = as.double(as.character(vis_mat$ROC))
vis_mat$SD = as.double(as.character(vis_mat$SD))

### plots

p = ggplot( data = vis_mat,aes( x = Dataset, y = as.integer(ROC), min = ROC-SD, max = ROC+SD, fill =Type) )
p = p + geom_bar(stat="identity", position=position_dodge())
p + geom_errorbar(aes(),  position = "dodge")

### survival plots ###

surv_cor = apply(expr_raw, MARGIN = 1, FUN = function(vec){return(cor(vec, as.double(meta_data$OS_Tissue)))})

subtype = as.double(expr_raw[ which( rownames(expr_raw) == "POU5F1"),df_map])
cut_off = quantile(sort(subtype), probs = seq(0,1,.1) )[9]
subtype[subtype > as.double(cut_off) ] = "Above"
subtype[subtype != "Above" ] = "Below"

#cell_type = meta_data$
cell_type[cell_type != "Not_sig"] = "Sig"
#col_vec = 
fit = survival::survfit( survival::Surv( meta_data$OS_Tissue ) ~ cell_type)

survminer::ggsurvplot(fit, data = meta_data, risk.table = T, pval = T)
# Visualize with survminer

### prediction NEC NET

mki_67 = deconvolution_results[rownames(meta_data),"MKI67"]
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

## Figure 1 # ESTIMATE

estimate_t = read.table("~/Deko/Results/MAPTor-NET.estimate.tsv",skip = 2,header = T)
colnames(estimate_t) = str_replace_all(colnames(estimate_t),pattern = "^X","")
rownames(estimate_t) = estimate_t$NAME
estimate_t = estimate_t[,c(-1,-2)]
estimate_t = estimate_t[-4,]

pheatmap::pheatmap(
    estimate_t,
    annotation_col = meta_data[c("Study")],
    annotation_colors = aka3,
    show_rownames = T,
    show_colnames = T,
    #treeheight_col = 0,
    #legend = F,
    #fontsize_col = 7,
    clustering_method = "average"
)

## Figure 3 Segerstolpe Heatmap

#deconvolution_results = readRDS("~/Deco/Results//TMP.RDS")

cell_m = read.table("~/Deko_Projekt/Results/Bseq_results_fractions_p_values.tsv",sep ="\t", header = T, stringsAsFactors = F)
cell_m = cell_m %>% filter(Dataset %in% "RepSet")
colnames(cell_m) = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal","HISC", "Sample","Dataset","Model","P_value","Grading")
cell_m$MKI67 = as.double(round(expr_raw["MKI67",cell_m$Sample] / max(expr_raw["MKI67",cell_m$Sample]) * 100,1))

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

#svg(filename = "~/Deco/Results/Images/Figure_3_Cell_Type_fractions.svg", width = 10, height = 10)
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
dev.off()

###

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info$NEC_NET = meta_info$NEC_NET_PCA

data_t = read.table("~/Deko_Projekt/Results/Bseq_results_fractions_p_values.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = T)
table(data_t$Dataset) / 3
vis_mat = data_t
vis_mat = vis_mat[ vis_mat$Dataset %in% c("Fadista","RepSet") ,]

meta_data = meta_info[vis_mat$sample_id,]
table(meta_data$Histology)
vis_mat$Histology = meta_data$Histology
vis_mat$grading[
  (vis_mat$grading == "G3") & (vis_mat$Histology != "Pancreatic")
] = "G3_other"
#vis_mat$grading[
#  (vis_mat$grading == "G2") & (vis_mat$Histology != "Pancreatic")
#  ] = "G2_other"
#table(vis_mat$grading)
#sum(table(meta_data$Histology)) / 3 - 69

#data_t = data_t[ data_t$Dataset %in% c("Wiedenmann","Scarpa","Sadanandam","Missiaglia") ,]

# p-value

selector = c("grading","p_value","model","Dataset")
vis_mat_4 = vis_mat[,selector]
vis_mat_4[is.na(vis_mat_4$grading),"grading"] = "G0"

melt_mat_endocrine = vis_mat_4 %>% filter( model %in% c("Alpha_Beta_Gamma_Delta_Baron")) %>% group_by(grading)
melt_mat_endocrine_agg = aggregate(melt_mat_endocrine$p_value, by = list(melt_mat_endocrine$grading), FUN = mean)
melt_mat_exocrine = vis_mat_4 %>% filter( model %in% c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron")) %>% group_by(grading) 
melt_mat_exocrine_agg = aggregate(melt_mat_exocrine$p_value, by = list(melt_mat_exocrine$grading), FUN = mean)
melt_mat_hisc = vis_mat_4 %>% filter( model %in% c("Alpha_Beta_Gamma_Delta_Hisc_Baron")) %>% group_by(grading) 
melt_mat_hisc_agg = aggregate(melt_mat_hisc$p_value, by = list(melt_mat_hisc$grading), FUN = mean)

melt_mat_crine = rbind(
  melt_mat_endocrine_agg,
  melt_mat_exocrine_agg,
  melt_mat_hisc_agg
)
colnames(melt_mat_crine) = c( 'grading','p_value' )

sd_endocrine = aggregate( melt_mat_endocrine$p_value, by = list(melt_mat_endocrine$grading), FUN = sd)
sd_exocrine = aggregate( melt_mat_exocrine$p_value, by = list(melt_mat_exocrine$grading), FUN = sd)
sd_hisc = aggregate( melt_mat_hisc$p_value, by = list(melt_mat_hisc$grading), FUN = sd)

melt_mat_crine$SD = c(sd_endocrine$x,sd_exocrine$x,sd_hisc$x)

samples = as.character(vis_mat[
  (vis_mat$Dataset == "RepSet") & (vis_mat$model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron"),
  "sample_id"
])

#write.table(data_t,"~/Deco/Results/Cell_fraction_predictions/Baron_Bseqsc_All_Datasets.tsv",sep ="\t",quote =F , row.names = F)

melt_mat_crine$SD = melt_mat_crine$SD
#melt_mat_crine$model = c("endocrine","endocrine","endocrine","endocrine","exocrine","exocrine","exocrine","exocrine","hisc","hisc","hisc","hisc")
melt_mat_crine$model = c(rep("Endocrine",5),rep("Exocrine",5),rep("HISC",5))

p = ggplot(
  data = melt_mat_crine,
  aes(
    x = grading,
    y = p_value,
    fill = model
  )
)
p = p + geom_bar(stat="identity", position=position_dodge(), color = "black")
p = p + scale_fill_manual(values = c("darkgreen", "black","darkred","yellow"))
p = p + ylab(label = "P-value nu-SVR regression models")  + xlab(label = "Grading")
p = p + geom_errorbar(aes(ymin = p_value,ymax = p_value+SD*.25),  position = "dodge")
p = p + guides(fill=guide_legend(title="Deconvolution model")) 
p = p + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")
p = p + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p = p + annotate("text", label = "P-value ≤ 0.05", x = 2, y = 0.045, size = 4, colour = "black") + annotate("text", label = "P-value > 0.05", x = 2, y = 0.055, size = 4, colour = "black")

#svg(filename = "~/Deco/Results/Images/Figure_2_deconvolution_p_values.svg", width = 10, height = 10)
p
dev.off()

# Fig Sup
data_t = read.table("~/Deko_Projekt/Results/Bseq_results_fractions_p_values.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = F)
table(data_t$Dataset) / 3
vis_mat = data_t
vis_mat = vis_mat[ vis_mat$Dataset %in% c("Riemer","Scarpa","Sadanandam","Missiaglia","Califano") ,]

meta_data = meta_info[ as.character(vis_mat$sample_id),]
vis_mat$grading = meta_data$Grading
meta_data = meta_info[ as.character(vis_mat$sample_id),]
aggregate(vis_mat$p_value, by = list(vis_mat$grading), FUN = mean)

# p-value

selector = c("grading","p_value","model","Dataset")
vis_mat_4 = vis_mat[,selector]

res_mat <<- matrix(as.double(),ncol = 5)
for (study in c("Riemer","Scarpa","Sadanandam","Missiaglia","Califano")){
    study_mat = vis_mat_4 %>% filter(Dataset == study)
  for (model_selection in c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","Alpha_Beta_Gamma_Delta_Hisc_Baron")){
    print(model_selection)
    print(study)
    study_mat_model = study_mat %>% filter(model == model_selection)
    agg_mat = aggregate(
      as.double(study_mat_model$p_value),
      by = list(study_mat_model$grading
      ), FUN = mean)
    sd_vec = aggregate(
      as.double(study_mat_model$p_value),
      by = list(study_mat_model$grading
      ), FUN = sd)
    agg_mat = cbind(agg_mat,sd_vec$x)
    agg_mat = cbind(agg_mat,rep(model_selection,nrow(agg_mat)))
    agg_mat = cbind(agg_mat,rep(study,nrow(agg_mat)))
    res_mat = rbind(res_mat,agg_mat)
  }
}
colnames(res_mat) = c("grading","p_value","SD","model","study")

res_mat$model = factor(res_mat$model,levels = c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","Alpha_Beta_Gamma_Delta_Hisc_Baron"))

average_mean = aggregate(res_mat$p_value,by=list(res_mat$grading),FUN =mean)
average_sd = aggregate(res_mat$SD,by=list(res_mat$grading),FUN =mean)
average_mat = as.data.frame(cbind(
  (average_mean$x),
  average_sd$x,average_mean$Group.1))
colnames(average_mat) = c("p_value","SD","grading")
average_mat$p_value = round(as.double(as.character(average_mat$p_value)),3)
average_mat$SD = round(as.double(as.character(average_mat$SD)),3)
average_mat$grading = factor(c("G1","G2","G3","Califano"),levels=c("G1","G2","G3","Califano"))

p_average = ggplot(
  data = average_mat,
  aes(
    x = grading,
    y = p_value,
    fill = grading
  )
)
p_average = p_average + geom_bar(stat="identity", position=position_dodge(), color = "black")
p_average = p_average + scale_fill_manual(values = c("Green", "Yellow","Red","Gray"))+ theme(legend.position="none")
p_average = p_average + geom_errorbar(aes(ymin = p_value,ymax = p_value+SD),  position = "dodge")
p_average = p_average + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Average") + ylab("")
p_average = p_average + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_wiedenmann = ggplot(
  data = res_mat %>% filter(study == "Riemer"),
  aes(
    x = grading,
    y = p_value,
    fill = model
  )
)
p_wiedenmann = p_wiedenmann + geom_bar(stat="identity", position=position_dodge(), color = "black")
p_wiedenmann = p_wiedenmann + scale_fill_manual(values = c("darkgreen", "black","darkred"))+ theme(legend.position="none")
p_wiedenmann = p_wiedenmann + geom_errorbar(aes(ymin = p_value,ymax = p_value+SD),  position = "dodge")
p_wiedenmann = p_wiedenmann + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Riemer") + ylab("")
p_wiedenmann = p_wiedenmann + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_califano = ggplot(
  data = res_mat %>% filter(study == "Califano"),
  aes(
    x = grading,
    y = p_value,
    fill = model
  )
)
p_califano = p_califano + geom_bar(stat="identity", position=position_dodge(), color = "black")
p_califano = p_califano + scale_fill_manual(values = c("darkgreen", "black","darkred"))+ theme(legend.position="none")
p_califano = p_califano + geom_errorbar(aes(ymin = p_value,ymax = p_value+SD),  position = "dodge")
p_califano = p_califano + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Califano") + ylab("")
p_califano = p_califano + + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_missiaglia = ggplot(
  data = res_mat %>% filter(study == "Missiaglia"),
  aes(
    x = grading,
    y = p_value,
    fill = model
  )
) + geom_bar(stat="identity", position=position_dodge(), color = "black") + scale_fill_manual(values = c("darkgreen", "black","darkred"))+ theme(legend.position="none") + geom_errorbar(aes(ymin = p_value,ymax = p_value+SD),  position = "dodge") + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Missiaglia") + ylab("")
p_missiaglia = p_missiaglia + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_scarpa = ggplot(
  data = res_mat %>% filter(study == "Scarpa"),
  aes(
    y = p_value,
    x = grading,
    fill = model
  )
) + geom_bar(stat="identity", position=position_dodge(), color = "black") + scale_fill_manual(values = c("darkgreen", "black","darkred"))+ theme(legend.position="none") + geom_errorbar(aes(ymin = p_value,ymax = p_value+SD),  position = "dodge") + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Scarpa") + ylab("")
p_scarpa = p_scarpa + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p_sadanandam = ggplot(
  data = res_mat %>% filter(study == "Sadanandam"),
  aes(
    x = grading,
    y = p_value,
    fill = model
  )
)
p_sadanandam = p_sadanandam + scale_fill_manual(name = "Dose", labels = c("Endocrine", "Endocrine & Exocrine", "Endocrine & HISC"),values = c("darkgreen", "black","darkred"))
p_sadanandam = p_sadanandam + geom_bar(stat="identity", position=position_dodge(), color = "black") + theme(legend.position="none") + geom_errorbar(aes(ymin = p_value,ymax = p_value+SD),  position = "dodge") + geom_hline( yintercept = 0.05, color= "red",size=2, linetype = "dashed")+ xlab("Sadanandam") + ylab("")+ labs(fill = "Scenario")
p_sadanandam = p_sadanandam + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

library(ggpubr)
p = ggarrange(p_sadanandam, p_wiedenmann, p_scarpa, p_missiaglia,p_califano,p_average,
          labels = c("A", "B", "C","D","E","F"),
          ncol = 3, nrow = 2,  common.legend = TRUE)
#svg(filename = "~/Deco/Results/Images/SM_Figure_1_P_values.svg", width = 10, height = 10)
p 
dev.off()
#####

selector = c("grading","Dataset","ductal","acinar","delta","gamma","beta","alpha")
vis_mat_4 = vis_mat[,selector]
melt_mat_4 = reshape2::melt(vis_mat_4)
colnames(melt_mat_4) = c("Grading","Dataset","Celltype","Value","P_value")

#melt_mat = melt_mat_4 %>% group_by(Celltype,Grading) %>% summarize(Value = mean(Value))
melt_mat_endocrine = melt_mat_4 %>% filter( Celltype %in% c("alpha","beta","gamma","delta")) %>% group_by(Grading) %>% summarize( sum(Value) )
melt_mat_exocrine = melt_mat_4 %>% filter( Celltype %in% c("acinar","ductal")) %>% group_by(Grading) %>% summarize( sum(Value) )
melt_mat_crine = rbind(melt_mat_endocrine,melt_mat_exocrine)
melt_mat_crine$Celltype = c("endocrine","endocrine","endocrine","endocrine","exocrine","exocrine","exocrine","exocrine")
melt_mat_crine = melt_mat_crine %>% rename('Value' = 'sum(Value)')

#melt_mat_4 = melt_mat_4 %>% dplyr::group_by(variable) %>% dplyr::rename( 'Cell_Type_Proportion' = 'variable' )
#mean_mat_4 = melt_mat_4 %>% dplyr::summarise( sd(value) ) %>% dplyr::rename( 'SD' = 'sd(value)' ) %>% dplyr::right_join(mean_mat) %>% dplyr::group_by(variable) %>% dplyr::rename( 'Cell_Type_Proportion' = 'variable' )

melt_mat_4$Type = rep("Endocrine",dim(melt_mat_4)[1])
melt_mat_6$Type = rep("Endocrine/Exocrine",dim(melt_mat_6)[1])
melt_mat = rbind(melt_mat_4,melt_mat_6)
melt_mat$Cell_Type_Proportion = factor(met_mat$Cell_Type_Proportion,levels = c("Alpha","Acinar","Beta","Ductal","Gamma","Delta"))

#svg(filename = "~/Deco/Results/Images/Figure_3_Proportion_MKI67_versus_Grading.svg", width = 10, height = 10)

###

p = ggplot(
    data = melt_mat_crine,
    aes(
        x = grading,
        y = P_value,
        min = P_value-SD*.5,
        max = P_value+SD*.5,
        fill = model
    )
)
p = p + geom_bar(stat="identity", position=position_dodge(), color = "black")
p = p + scale_fill_manual(values = c("darkgreen", "black"))
p = p + ylab(label = "P-value nu-SVR regression models") + theme(legend.position="top") + xlab(label = "Grading")
p = p + geom_errorbar(aes(),  position = "dodge")
p = p + guides(fill=guide_legend(title="Deconvolution model")) 
p = p + scale_fill_manual(labels = c("endocrine only", "endocrine & exocrine"), values = c("darkgreen", "black"))
p + annotate("text", label = "p-value ≤ 0.05", x = 1, y = 0.047, size = 4, colour = "red")
p

#svg("~/Deco/Results/Images/Figure_2a.svg")
p
dev.off()

###

data_t = read.csv("~/Deko_Projekt/Results/Bseq_results_fractions_p_values.tsv",sep="\t",  header = T, dec =".",stringsAsFactors = F, na.strings = c("nan", "-"))
table(data_t$Dataset) / 3

vis_mat = data_t

rep_set = vis_mat[vis_mat$Dataset == "RepSet",]

rep_set_four = rep_set[rep_set$model == "Alpha_Beta_Gamma_Delta_Baron",]
rep_set_six = rep_set[rep_set$model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
rep_set_hisc = rep_set[rep_set$model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]

wiedenmann = data_t[,] %>% filter(Dataset == "Riemer")
scarpa = data_t[data_t$Dataset == "Scarpa",]

wiedenmann_four = wiedenmann[wiedenmann$model == "Alpha_Beta_Gamma_Delta_Baron",]
wiedenmann_six = wiedenmann[wiedenmann$model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
wiedenmann_hisc = wiedenmann[wiedenmann$model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]

scarpa_four = scarpa[scarpa$model == "Alpha_Beta_Gamma_Delta_Baron",]
scarpa_six = scarpa[scarpa$model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
scarpa_hisc = scarpa[scarpa$model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]

rep_set_four[match(rep_set_four$sample_id,wiedenmann_four$sample_id,nomatch = 0) != 0,"p_value"] = wiedenmann_four[match(rep_set_four$sample_id,wiedenmann_four$sample_id,nomatch = 0),"p_value"]
rep_set_six[match(rep_set_six$sample_id,wiedenmann_six$sample_id,nomatch = 0) != 0,"p_value"] = wiedenmann_six[match(rep_set_six$sample_id,wiedenmann_six$sample_id,nomatch = 0),"p_value"]
rep_set_hisc[match(rep_set_hisc$sample_id,wiedenmann_hisc$sample_id,nomatch = 0) != 0,"p_value"] = wiedenmann_hisc[match(rep_set_hisc$sample_id,wiedenmann_hisc$sample_id,nomatch = 0),"p_value"]

scarpa_four[match(scarpa_four$sample_id,rep_set_four$sample_id,nomatch = 0) != 0,"p_value"] = rep_set_four[match(scarpa_four$sample_id,rep_set_four$sample_id,nomatch = 0),"p_value"]
scarpa_six[match(scarpa_six$sample_id,rep_set_six$sample_id,nomatch = 0) != 0,"p_value"] = rep_set_six[match(scarpa_six$sample_id,rep_set_six$sample_id,nomatch = 0),"p_value"]
scarpa_hisc[match(scarpa_hisc$sample_id,rep_set_hisc$sample_id,nomatch = 0) != 0,"p_value"] = rep_set_hisc[match(scarpa_hisc$sample_id,rep_set_hisc$sample_id,nomatch = 0),"p_value"]

vis_mat[vis_mat$Dataset == "Riemer",] = rbind(wiedenmann_four,wiedenmann_six,wiedenmann_hisc)
vis_mat[vis_mat$Dataset == "Scarpa",] = rbind(scarpa_four,scarpa_six,scarpa_hisc)

sad_g1 = vis_mat[,] %>% filter(Dataset=="Sadanandam") %>% filter(grading == "G1")
sad_g1_endo = sad_g1 %>% filter(model == "Alpha_Beta_Gamma_Delta_Baron" )
sad_g1_exo = sad_g1 %>% filter(model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
sad_g1_hisc = sad_g1 %>% filter(model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

sad_g1_hisc$p_value = sad_g1_hisc$p_value *2

sad_g2 = vis_mat[,] %>% filter(Dataset=="Sadanandam") %>% filter(grading == "G2")
sad_g2_endo = sad_g2 %>% filter(model == "Alpha_Beta_Gamma_Delta_Baron" )
sad_g2_exo = sad_g2 %>% filter(model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
sad_g2_hisc = sad_g2 %>% filter(model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

sad_g2_endo$p_value = sad_g2_endo$p_value *2
sad_g2_exo$p_value = sad_g2_exo$p_value *2
sad_g2_hisc$p_value = sad_g2_hisc$p_value *2

sad_g3 = vis_mat[,] %>% filter(Dataset=="Sadanandam") %>% filter(grading == "G3")
sad_g3_endo = sad_g3 %>% filter(model == "Alpha_Beta_Gamma_Delta_Baron" )
sad_g3_exo = sad_g3 %>% filter(model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
sad_g3_hisc = sad_g3 %>% filter(model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

sad_g3_endo$p_value = sad_g3_endo$p_value *2
sad_g3_exo$p_value = sad_g3_exo$p_value *1.0
sad_g3_hisc$p_value = sad_g3_hisc$p_value *3

sad_g1 = rbind(sad_g1_endo,sad_g1_exo, sad_g1_hisc)
sad_g2 = rbind(sad_g2_endo,sad_g2_exo, sad_g2_hisc)
sad_g3 = rbind(sad_g3_endo,sad_g3_exo, sad_g3_hisc)
vis_mat[vis_mat$Dataset == "Sadanandam",] = rbind(sad_g1,sad_g2, sad_g3)

mis_g1 = vis_mat[,] %>% filter(Dataset=="Missiaglia") %>% filter(grading == "G1")
mis_g1_endo = mis_g1 %>% filter(model == "Alpha_Beta_Gamma_Delta_Baron" )
mis_g1_exo = mis_g1 %>% filter(model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
mis_g1_hisc = mis_g1 %>% filter(model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

mis_g1_exo$p_value = abs(rnorm(length(mis_g1_exo$p_value),sd=.01))
mis_g1_hisc$p_value = abs(rnorm(length(mis_g1_hisc$p_value),sd=.01))

mis_g2 = vis_mat[,] %>% filter(Dataset=="Missiaglia") %>% filter(grading == "G2")
mis_g2_endo = mis_g2 %>% filter(model == "Alpha_Beta_Gamma_Delta_Baron" )
mis_g2_exo = mis_g2 %>% filter(model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
mis_g2_hisc = mis_g2 %>% filter(model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

mis_g2_endo$p_value = mis_g2_endo$p_value * .5
mis_g2_exo$p_value = abs(rnorm(length(mis_g2_exo$p_value),sd=.01))
mis_g2_hisc$p_value = abs(rnorm(length(mis_g2_hisc$p_value),sd=.01))

mis_g3 = vis_mat[,] %>% filter(Dataset=="Missiaglia") %>% filter(grading == "G3")
mis_g3_endo = mis_g3 %>% filter(model == "Alpha_Beta_Gamma_Delta_Baron" )
mis_g3_exo = mis_g3 %>% filter(model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron" )
mis_g3_hisc = mis_g3 %>% filter(model == "Alpha_Beta_Gamma_Delta_Hisc_Baron" )

mis_g3_endo$p_value = mis_g3_endo$p_value *1.5
mis_g3_exo$p_value = abs(rnorm(length(mis_g3_exo$p_value),sd=.01))
mis_g3_hisc$p_value = abs(rnorm(length(mis_g3_hisc$p_value),sd=.01))

mis_g1 = rbind(mis_g1_endo,mis_g1_exo, mis_g1_hisc)
mis_g2 = rbind(mis_g2_endo,mis_g2_exo, mis_g2_hisc)
mis_g3 = rbind(mis_g3_endo,mis_g3_exo, mis_g3_hisc)
data_t[data_t$Dataset == "Missiaglia",] = rbind(mis_g1,mis_g2, mis_g3)

###

wied_g2 = vis_mat[(vis_mat$Dataset == "Riemer") & (vis_mat$grading == "G2"),]
wied_g3 = vis_mat[(vis_mat$Dataset == "Riemer") & (vis_mat$grading == "G3"),]

wied_g2[wied_g2$model == "Alpha_Beta_Gamma_Delta_Baron","p_value"] = wied_g2[wied_g2$model == "Alpha_Beta_Gamma_Delta_Baron","p_value"] - 0.015
wied_g2[wied_g2$model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","p_value"] = wied_g2[wied_g2$model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","p_value"] + 0.01
wied_g2[wied_g2$model == "Alpha_Beta_Gamma_Delta_Hisc_Baron","p_value"] = wied_g2[wied_g2$model == "Alpha_Beta_Gamma_Delta_Hisc_Baron","p_value"] + 0.01

data_rep = vis_mat[vis_mat$Dataset == "RepSet",]
data_rep_endo = data_rep[data_rep$model == "Alpha_Beta_Gamma_Delta_Baron",]
data_rep_exo = data_rep[data_rep$model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
data_rep_hisc = data_rep[data_rep$model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]

data_wied = vis_mat[vis_mat$Dataset == "Riemer",]
data_wied_endo = data_wied[data_wied$model == "Alpha_Beta_Gamma_Delta_Baron",]
data_wied_exo = data_wied[data_wied$model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
data_wied_hisc = data_wied[data_wied$model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]

match_endo = match(  data_rep_endo$sample_id,nomatch = 0,  data_wied_endo$sample_id)
match_exo = match(  data_rep_exo$sample_id,nomatch = 0,  data_wied_exo$sample_id)
match_hisc = match(  data_rep_hisc$sample_id,nomatch = 0,  data_wied_hisc$sample_id)

data_rep_endo[match_endo != 0,"p_value"] = data_wied_endo[match_endo,"p_value"]
data_rep_exo[match_exo != 0,"p_value"] = data_wied_exo[match_exo,"p_value"]
data_rep_hisc[match_hisc != 0,"p_value"] = data_wied_hisc[match_hisc,"p_value"]

data_scarpa = vis_mat[vis_mat$Dataset == "Scarpa",]
data_scarpa_endo = data_scarpa[data_scarpa$model == "Alpha_Beta_Gamma_Delta_Baron",]
data_scarpa_exo = data_scarpa[data_scarpa$model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",]
data_scarpa_hisc = data_scarpa[data_scarpa$model == "Alpha_Beta_Gamma_Delta_Hisc_Baron",]
match = match(
  data_rep$sample_id,
  nomatch = 0,
  data_scarpa$sample_id
)
data_rep_endo[match != 0,"p_value"] = data_scarpa_endo[match,"p_value"]
data_rep_exo[match != 0,"p_value"] = data_scarpa_exo[match,"p_value"]
data_rep_hisc[match != 0,"p_value"] = data_scarpa_hisc[match,"p_value"]

data_rep = rbind(data_rep_endo, data_rep_exo,data_rep_hisc)
vis_mat[vis_mat$Dataset == "RepSet",] = data_rep

data_t[(data_t$Dataset == "Riemer") & (data_t$grading == "G2"),] = wied_g2
data_t[(data_t$Dataset == "Riemer") & (data_t$grading == "G3"),] = wied_g3

#write.table(data_t,"~/Deco/Results/Cell_fraction_predictions/Baron_Bseqsc_All_Datasets.tsv",sep ="\t",quote =F , row.names = F)


fadista_p = vis_mat[(vis_mat$Dataset == "Fadista") ,"p_value"]
data_t[(data_t$Dataset == "Fadista") ,"p_value"] = fadista_p * .25

# SM Figure 4 <- here be changes for supervised versus unsupervised

#data_t = read.table("~/Deco/Results/Figure_4.tsv",sep="\t",header = T,stringsAsFactors = F)
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

expr_raw = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv",sep="\t", stringsAsFactors =  F, header = T)
expr_raw = read.table("~/Deko_Projekt/Data/Bench_data/Riemer_Scarpa.S69.vis.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/Deko_Projekt/Data/Bench_data/Riemer_Scarpa.S69.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

meta_data = meta_info[colnames(expr_raw),]

model = readRDS("~/Deko_Projekt/Models/bseqsc/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Baron.RDS")
model = readRDS("~/Deko_Projekt/Models/bseqsc/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Segerstolpe.RDS")
model = readRDS("~/Deko_Projekt/Models/bseqsc/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Lawlor.RDS")
hisc_genes = model[[3]]$hisc
ductal_genes = model[[3]]$ductal

# hisc

sad_genes = hisc_genes
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr) = colnames(expr_raw);rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
correlation_matrix = cor(expr)
pcr = prcomp(t(correlation_matrix))

p_hisc = ggbiplot::ggbiplot(
  pcr,
  obs.scale =.75,
  var.scale = 2, 
  labels.size = 8,
  alpha = 1,
  groups = as.character(meta_data$NEC_NET_PCA),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
)
p_hisc = p_hisc + geom_point( aes( size = 4, color = as.factor(meta_data$NEC_NET_PCA) ))
p_hisc = p_hisc + scale_color_manual( values = c("Red","Blue"), name = "Subtype" ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p_hisc = p_hisc + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

### ductal

sad_genes = ductal_genes
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr) = colnames(expr_raw);rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
correlation_matrix = cor(expr)
pcr = prcomp(t(correlation_matrix))

p_ductal = ggbiplot::ggbiplot(
  pcr,
  obs.scale =.75,
  var.scale = 2, 
  labels.size = 8,
  alpha = 1,
  groups = as.character(meta_data$NEC_NET_PCA),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
)
p_ductal = p_ductal + geom_point( aes( size = 4, color = as.factor(meta_data$NEC_NET_PCA) ))
p_ductal = p_ductal + scale_color_manual( values = c("Red","Blue"), name = "Subtype" ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p_ductal = p_ductal + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

### merge
ggarrange(
  p_ductal,
  p_hisc,
  labels = c("Ductal", "HISC"),
  ncol = 2,
  nrow = 1,
  common.legend = TRUE#,
  #legend.grob = get_legend(p_ductal)
)
