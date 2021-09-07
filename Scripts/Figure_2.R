library("stringr")
library("reshape2")
library("dplyr")

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Diedisheim_S66.absolute.exocrine.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",rownames(props)[no_match],sep ="")
no_match = rownames(props) %in% meta_info$Sample == F
sum(no_match)

dim(props)
meta_data = meta_info[rownames(props),]

props = props[meta_data$Functionality %in% c("Non-functional","Insulinoma","Glucagonoma", "Somatostatinoma", "PPoma", "Non-functional"),]
meta_data = meta_info[rownames(props),]
table(meta_data$NEC_NET_Color)

selection = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal")
exocrines = as.double(rowSums(props[,c("Ductal","Acinar")]))
endocrines = as.double(rowSums(props[,c("Alpha","Beta","Gamma","Delta")]))

meta_data$Ratio = log((exocrines+.1) / (endocrines+.1))

vis_mat = props[,selection]
vis_mat$Exocrine = vis_mat$Acinar + vis_mat$Ductal
vis_mat = vis_mat[,!(colnames(vis_mat) %in% c("Acinar","Ductal"))]

correlation_matrix = cor(t(vis_mat))

pcr = prcomp(t(correlation_matrix))

p = pheatmap::pheatmap(
    t(vis_mat),
    annotation_col = meta_data[,c("Grading","Cluster","Functionality","NEC_NET_Ori")],
    annotation_colors = aka3,
    show_rownames = T,
    show_colnames = F,
    treeheight_row = 0,
    cellheight = 20,
    cluster_rows = F,
    legend = T,
    fontsize_row = 14,
    clustering_method = "average"
)
p = p +  theme(legend.position="top",axis.text=element_text(size=18),axis.title=element_text(size=18))+ theme(legend.text=element_text(size=18),legend.title=element_text(size=18))
p
