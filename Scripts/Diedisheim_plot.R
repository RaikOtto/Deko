library("stringr")
library("reshape2")
library("dplyr")
library("umap")
library("ggplot2")
library("ggpubr")
library("grid")

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Absolute/Baron_exocrine/Diedisheim.S62.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",rownames(props)[no_match],sep ="")
no_match = rownames(props) %in% meta_info$Sample == F
sum(no_match)

dim(props)
meta_data = meta_info[rownames(props),]

props = props[meta_data$Functionality %in% c("Insulinoma","Glucagonoma", "Somatostatinoma", "PPoma", "Non-Functional","Unknown","VIPoma","ACTH","Gastrinoma"),]
meta_data = meta_info[rownames(props),]
table(meta_data$NEC_NET)

selection = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal")
exocrines = as.double(rowSums(props[,c("Ductal","Acinar")]))
endocrines = as.double(rowSums(props[,c("Alpha","Beta","Gamma","Delta")]))

meta_data$Ratio = log((exocrines+.1) / (endocrines+.1))

vis_mat = props[,selection]
vis_mat$Exocrine_like = vis_mat$Acinar + vis_mat$Ductal
vis_mat = vis_mat[,!(colnames(vis_mat) %in% c("Acinar","Ductal"))]

correlation_matrix = cor(t(vis_mat))

pcr = prcomp(t(correlation_matrix))

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
upper_plot = pheatmap::pheatmap(
    t(vis_mat),
    annotation_col = meta_data[,c("Grading","Cluster","Functionality","NEC_NET")],
    annotation_colors = aka3,
    show_rownames = T,
    show_colnames = F,
    treeheight_row = 0,
    cellheight = 20,
    cluster_rows = FALSE,
    legend = T,
    fontsize_row = 14,
    clustering_method = "single"
)
#p = p +  theme(legend.position="top",axis.text=element_text(size=18),axis.title=element_text(size=18))+ theme(legend.text=element_text(size=18),legend.title=element_text(size=18))
upper_plot = upper_plot + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="top")
upper_plot

