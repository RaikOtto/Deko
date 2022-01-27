library("umap")
library("ggpubr")
library("stringr")
library("reshape2")
library("ggplot2")
library("dplyr")
library("grid")

meta_info_maptor = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info_maptor) = meta_info_maptor$Sample
colnames(meta_info_maptor) = str_replace(colnames(meta_info_maptor),pattern = "\\.","_")

meta_info_maptor$OS_Tissue = as.double(str_replace(meta_info_maptor$OS_Tissue,pattern = ",","."))

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

matcher = match(meta_info_maptor$Sample,meta_info$Sample, nomatch = 0)
meta_info[matcher,"OS_Tissue"] = meta_info_maptor[matcher != 0,"OS_Tissue"]

#expr_raw = read.table("~/Deko_Projekt/Data/Publication_datasets/NEN/Sato.S22.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet_Cibersort_Tosti_200_genes_200_samples_endocrine_exocrine_metaplastic_acinar-i_muc5+_only.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = TRUE)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "\\.", "")
expr_raw[1:5,1:5]
dim(expr_raw)

expr_raw = expr_raw[,!(colnames(expr_raw) %in% c("","X","RMSE","Correlation","P_value","Subtype","Strength_subtype","model","Sig_score"))]
expr_raw = t(expr_raw)

no_match = colnames(expr_raw) %in% meta_info$Sample == FALSE
#colnames(expr_raw)[no_match] = str_replace(colnames(expr_raw)[no_match], pattern = "^X","")
#no_match = colnames(expr_raw) %in% meta_info$Sample == F
#colnames(expr_raw)[no_match] = paste("X",colnames(expr_raw)[no_match],sep ="")
#no_match = colnames(expr_raw) %in% meta_info$Sample == F
#colnames(expr_raw)[which(no_match)]

candidates = meta_data$Sample[ 
  #meta_data$Study %in% c("Master","Charite")
  #meta_data$Site_of_primary %in% c("Pancreatic") #& meta_data$Primary_Metastasis %in% c("Primary")
  meta_data$Site_of_primary != "Pancreatic" #& meta_data$Primary_Metastasis %in% c("Primary")
  #meta_data$NET_NEC_PCA %in% c("NEC","NET")
  #!(meta_data$Histology_Metastasis %in% c("Outlier"))
]
length(candidates)
expr_raw = expr_raw[,candidates]
meta_data = meta_info[colnames(expr_raw),]
dim(meta_data)

#write.table(expr_raw,"~/Deko_Projekt/Data/Publication_datasets/NEN/Diedisheim.S4.tsv",sep ="\t",quote =F , row.names = TRUE)

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
genes_of_interest_hgnc_t = read.table("~/Deko_Projekt/Misc/Stem_signatures.gmt.tsv",sep ="\t", stringsAsFactors = F, header = F)
#liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]

genes_of_interest_hgnc_t$V1
i = 72
genes_of_interest_hgnc_t[i,1]

sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
length(sad_genes)
meta_data = meta_info[colnames(expr_raw),]

expr = expr_raw[rownames(expr_raw) %in% sad_genes[1:20],]
expr[1:5,1:5]
dim(expr)
row_var = as.double(apply(expr, MARGIN = 1, FUN= var))
summary(row_var)
#expr = expr[row_var > 3,]
dim(expr)

correlation_matrix = cor(t(expr))
pcr = prcomp((correlation_matrix))

#meta_exp = as.double(apply(expr, MARGIN = 1, FUN = mean))
#expr = expr[,order(meta_exp)]
#svg(filename = "~/Downloads/Heatmap.svg", width = 10, height = 10)
p  =pheatmap::pheatmap(
  #correlation_matrix,
  t(expr),
  annotation_col = meta_data[,c("NEC_NET","Grading","Primary_Metastasis")],
  #annotation_col = meta_data[,c("NEC_NET_Color","Histology")],
  annotation_colors = aka3,
  show_rownames = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = T,
  fontsize_col = 7,
  #cellheight = 15,
  clustering_method = "average"
)

p = ggbiplot::ggbiplot(
  pcr,
  obs.scale =.75,
  var.scale = 2, 
  labels.size = 4,
  alpha = 1,
  groups = as.character(meta_data$NET_NEC_PCA),
  #label = meta_data$Sample,
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
)
p = p + geom_point( aes( size = 4, color = as.factor(meta_data$NET_NEC_PCA) ))
p
p = p + scale_color_manual( values = c("green","yellow","darkred","blue","red"), name = "Subtype" ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
#p = p + scale_color_manual( values = c("Red","Blue"), name = "Subtype" ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p = p + scale_color_manual( values = c("Purple","Red","Blue") ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p = p + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
#svg(filename = "~/Deco/Results/Images/SM_Figure_5_NEC_NET_PCA.svg", width = 10, height = 10)
p
#dev.off()

#p + xlim(c(-1.0,2.25)) + ylim(-1.5,1.0)


###

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info$NEC_NET = meta_info$NEC_NET_PCA

data_t = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet_S57_CIBERSORT_Tosti_400.Absolute.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = T)
#data_t = read.table("~/Deko_Projekt/Results/Bseq_results_fractions_p_values.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = T)
table(data_t$Dataset) / 3
vis_mat = data_t
vis_mat = vis_mat[ vis_mat$Dataset %in% c("Fadista","RepSet") ,]

####

meta_data = meta_info[vis_mat$Sample,]
table(meta_data$Histology)
vis_mat$Histology = meta_data$Histology
vis_mat$Grading[
  (vis_mat$Grading == "G3") & (vis_mat$Histology != "Pancreatic")
] = "G3_other"

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

#write.table(meta_info,"~/Deko_Projekt/Misc/Meta_information.tsv",sep ="\t",quote =F , row.names = F)

melt_mat_crine$SD = melt_mat_crine$SD
#melt_mat_crine$model = c("endocrine","endocrine","endocrine","endocrine","exocrine","exocrine","exocrine","exocrine","hisc","hisc","hisc","hisc")
melt_mat_crine$Model = c(rep("Endocrine-only",5),rep("Endocrine & Exocrine",5),rep("Endocrine & HISC",5))
melt_mat_crine$Model = factor(melt_mat_crine$Model,  levels =  c("Endocrine-only","Endocrine & Exocrine","Endocrine & HISC"))
melt_mat_crine = melt_mat_crine[,]
#melt_mat_crine = melt_mat_crine %>% filter(Grading != "G3_other")

melt_mat_crine = props[,c("endocrine cell","ductal cell type 1","acinar cell","acinar_edge_cell")]
melt_mat_crine_save = melt_mat_crine = t(apply( melt_mat_crine, MARGIN = 1, FUN = function(vec){return((vec/sum(vec))*100)}))
melt_mat_crine = as.matrix(melt_mat_crine, nrow = nrow(melt_mat_crine_save), ncol = ncol(melt_mat_crine_save))

melt_mat_crine$NEC_NET= meta_data$NET_NEC_PCA
melt_mat_crine$Grading= meta_data$Grading
melt_mat_crine = melt_mat_crine %>% dplyr::filter(Grading != "Unknown")
melt_mat_crine[melt_mat_crine$NEC_NET == "NEC","Grading"] = "G3_NEC"
melt_mat_crine[(melt_mat_crine$NEC_NET == "NET") & melt_mat_crine$Grading == "G3","Grading"] = "G3_NET"

melt_mat_crine = reshape2::melt(melt_mat_crine)
colnames(melt_mat_crine) = c("Sample","Cell_type","Prediction")
#colnames(melt_mat_crine) = c("NEC_NET","Grading","Cell_type","Prediction")
melt_mat_crine$Grading = meta_data[,"Grading"]
melt_mat_crine_vis = melt_mat_crine %>% group_by(Grading,Cell_type) %>% summarise("Average_Absolute_Prediction" = mean(Prediction)) 

p = ggplot(
  data = melt_mat_crine_vis,
  aes(
    x = Grading,
    y = Average_Absolute_Prediction,
    fill = Cell_type
  )
)
p = p + geom_bar(stat="identity",  color = "black",position = "dodge")

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

table(meta_data_sato_gep$Grading,meta_data_sato_gep$NEC_NET)
