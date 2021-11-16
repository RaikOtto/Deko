library("devtools")
load_all("~/artdeco")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")
library("stringr")
library("reshape2")
library("bseqsc")
library("dplyr")

#Alvarez.S105.tsv #Charite.S23.tsv #Diedisheim.S62.tsv #Master.S20.tsv #Missiaglia.S75.tsv #Sadanandam.S29.tsv

i_filename = "~/Deko_Projekt/Data/Publication_datasets/Diedisheim.S62.tsv"
expr_raw = read.table(i_filename,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
study_name = tail(str_split(i_filename,"/")[[1]],1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
dim(expr_raw)

#show_models_bseqsc()
model_name = "Alpha_Beta_Gamma_Delta_Baron"
#model_name = "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron"

props = Deconvolve_transcriptome(
    transcriptome_data = expr_raw,
    deconvolution_algorithm = "bseqsc",
    models = model_name,
    Cibersort_absolute_mode = FALSE,
    nr_permutations = 1000,
    output_file = ""
)

o_filename = "~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Baron_endocrine"
o_filename = paste(o_filename, study_name, sep ="/")
write.table(props,o_filename,sep = "\t")

###

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet_S103.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)

#props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Alvarez.S104.Cibersort.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"
colnames(props) = str_replace(colnames(props),pattern = "\\.","-")
colnames(props) = str_replace(colnames(props),pattern = "\\.ductal","ductal")
colnames(props) = str_replace(colnames(props),pattern = "-reg\\.","-reg+")
colnames(props) = str_replace(colnames(props),pattern = "muc5b-","muc5b+-")

no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",rownames(props)[no_match],sep ="")
no_match = rownames(props) %in% meta_info$Sample == F
sum(no_match)

dim(props)
meta_data = meta_info[rownames(props),]

###

#props = props[(meta_data$Histology == "pancreas") | (meta_data$NEC_NET_Color != "Primary") ,]
meta_data = meta_info[rownames(props),]
#props = props[(meta_data$Study == "Alvarez") ,]
meta_data = meta_info[rownames(props),]

selection = colnames(props)[!(colnames(props) %in% c("model","Sig_score","P_value","Correlation","RMSE"))]
selection = c("Acinar","acinar-s","acinar-i","acinar-reg+","Ductal","Beta","Delta","muc5b+ ductal","Alpha","Gamma")

###

props = as.data.frame(props)
vis_mat = props[,colnames(props) %in% selection]
vis_mat$endocrine_fully_differentited = as.double(rowSums(vis_mat[,c("Alpha","Beta","Gamma","Delta")]))
vis_mat$exocrine_fully_differentiated = as.double(rowSums(vis_mat[,colnames(vis_mat) %in%c("Ductal","acinar-s","Acinar")]))
vis_mat$metaplastic_not_fully_differentiated = as.double(rowSums(vis_mat[,colnames(vis_mat) %in% c("acinar-i","acinar-reg+","muc5b+ ductal")]))

correlation_matrix = cor(t(vis_mat));pcr = prcomp(t(correlation_matrix))
vis_mat = vis_mat[order(vis_mat$endocrine_fully_differentited,decreasing = T),]

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
p = pheatmap::pheatmap(
    t(vis_mat),
    #correlation_matrix,
    annotation_col = meta_data[,c("Acinar_reg","Grading","NET_NEC_PCA","Study")],
    annotation_colors = aka3,
    show_rownames = T,
    show_colnames = F,
    treeheight_row = 0,
    cluster_rows = F,
    cluster_cols = T,
    legend = F,
    fontsize_row = 14,
    clustering_method = "average"
)
p +  theme(legend.position="top",axis.text=element_text(size=18),axis.title=element_text(size=18))+ theme(legend.text=element_text(size=18),legend.title=element_text(size=18))

### PCA

p = ggbiplot::ggbiplot(
    prcomp(t(vis_mat)),
    groups = as.character(meta_data$Grading),
    var.axes = F,
    ellipse = TRUE
)
ki_67_vec = as.double(expr_raw["MKI67",])*.5

p = p + geom_point( aes( shape = meta_data$Grading, color = meta_data$Grading ), size = ki_67_vec ) # Fig 4
p = p + scale_color_manual( values = c("darkgreen","yellow","red") ) #Fig 4 Master
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")
p

