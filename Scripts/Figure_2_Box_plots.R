library("umap")
library("reshape2")
library("hrbrthemes")
library("waffle")
library(tidyverse)
library("stringr")
library("dplyr")
library("ggpubr")
library("png")
library("grid")
library("ggplot2")
library("magick")
library("treemapify")

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

# Figure 2 Plot A

### p-values

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"
meta_data = meta_info[props$Sample,]

props$Study = meta_data$Study

selection = c("Study","P_value","Model")

vis_mat = props[,selection]
vis_mat$P_value = as.double(vis_mat$P_value)
vis_mat$Model = factor(vis_mat$Model, levels = c("Endocrine_only","Endocrine_exocrine_like"))
vis_mat$Study = factor(vis_mat$Study, levels = c("Alvarez","Charite","Master","Scarpa","Diedisheim","Missiaglia","Sadanandam"))
#vis_mat[vis_mat$Study == "Diedisheim","P_value"] =vis_mat[vis_mat$Study == "Diedisheim","P_value"] / 3

p_value_plot = ggplot(vis_mat, aes( x = Study, y = P_value, fill = Model) )
p_value_plot = p_value_plot + geom_boxplot(notch = TRUE,outlier.colour = "red", outlier.shape = 1)
p_value_plot = p_value_plot + ylim(c(0,0.1)) + geom_hline(yintercept = 0.05, color = "red",linetype="dashed", size =2)
p_value_plot = p_value_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
p_value_plot = p_value_plot + scale_fill_manual(values = c("blue","red"))
p_value_plot = p_value_plot + ylab("Mean P-values") + theme(legend.position = "top")
p_value_plot

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_A.svg", width = 10, height = 10)
p_value_plot
dev.off()

# Figure 2 Plot B Metastasis/ Primary/ Organoid

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.with_NENs.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"
meta_data = meta_info[props$Sample,]

props$Primary_Metastasis = meta_data$Primary_Metastasis
props = props %>% filter(!( Primary_Metastasis %in% c("Unknown","Outlier")))

selection = c("Type","Model","Primary_Metastasis","P_value")

vis_mat = props[,selection]
vis_mat$P_value = as.double(vis_mat$P_value)

# pannen

vis_mat_pannen = vis_mat %>% filter(Type == "PanNEN")
vis_mat_pannen$Model = factor( as.character(vis_mat_pannen$Model), levels = c("Endocrine_only","Endocrine_exocrine_like"))
vis_mat_pannen$Primary_Metastasis = factor(vis_mat_pannen$Primary_Metastasis, levels = c("Primary","Metastasis"))

p_value_plot = ggplot(vis_mat_pannen, aes( x = Primary_Metastasis, y = P_value, fill = Model) )
p_value_plot = p_value_plot + geom_boxplot(notch = TRUE,outlier.colour = "red", outlier.shape = 1)
p_value_plot = p_value_plot + scale_fill_manual(values = c("blue","red"))
p_value_plot = p_value_plot + ylim(c(0,0.1)) + geom_hline(yintercept = 0.05, color = "red",linetype="dashed", size =2)
p_value_plot = p_value_plot + ylab("Mean P-values") + theme(legend.position = "none")
p_value_plot

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_B_pannen.svg", width = 10, height = 10)
p_value_plot
dev.off()

# nen

vis_mat_nen = vis_mat %>% filter(Type == "NEN")
vis_mat_nen$Model = factor(vis_mat_nen$Model, levels = c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","Organoid"))
vis_mat_nen$Primary_Metastasis = factor(vis_mat_nen$Primary_Metastasis, levels = c("Primary","Metastasis","Organoid"))
p_value_plot = ggplot(vis_mat_nen, aes( x = Primary_Metastasis, y = P_value, fill = Model) )
p_value_plot = p_value_plot + geom_boxplot(notch = TRUE,outlier.colour = "red", outlier.shape = 1)
p_value_plot = p_value_plot + scale_fill_manual(values = c("blue","red"))
p_value_plot = p_value_plot + ylim(c(0,0.1)) + geom_hline(yintercept = 0.05, color = "red",linetype="dashed", size =2)
#p_value_plot = p_value_plot + ylab("Mean P-values") + theme(legend.position = "none")
p_value_plot

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_B_nen.svg", width = 10, height = 10)
p_value_plot
dev.off()

### Figure 2 Plot C Grading

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"
meta_data = meta_info[props$Sample,]

props$Grading = meta_data$Grading
props[(meta_data$NET_NEC_PCA == "NET") & (meta_data$Grading == "G3"),"Grading"] = "G3_NET"
props[(meta_data$NET_NEC_PCA == "NEC"),"Grading"] = "G3_NEC"
props = props[props$Grading != "G3",]
props = props[props$Grading != "Unknown",]
meta_data = meta_info[props$Sample,]
props[props$Grading == "Organoid","Grading"] = "G3_NEC"

selection = c("Grading","P_value","Model")

vis_mat = props[,selection]
vis_mat$P_value = as.double(vis_mat$P_value)

vis_mat = vis_mat[meta_data$Study %in% c("Charite","Scarpa","Master"),]

vis_mat_endocrine = vis_mat[vis_mat$Model == "Endocrine_only",]
vis_mat_exocrine = vis_mat[vis_mat$Model == "Endocrine_exocrine_like",]

vis_mat_mean_endo = aggregate(vis_mat_endocrine$P_value, FUN = mean, by = list(vis_mat_endocrine$Grading))
vis_mat_mean_exo = aggregate(vis_mat_exocrine$P_value, FUN = mean, by = list(vis_mat_exocrine$Grading))
vis_mat_mean_endo$Model = rep("Endocrine_only",rep(nrow(vis_mat_mean_endo)))
vis_mat_mean_exo$Model = rep("Exocrine_like",rep(nrow(vis_mat_mean_exo)))
colnames(vis_mat_mean_endo) = colnames(vis_mat_mean_exo) = c("Grading","P_value","Model")

vis_mat_sd_endo = aggregate(vis_mat_endocrine$P_value, FUN = sd, by = list(vis_mat_endocrine$Grading))
vis_mat_sd_exo  = aggregate(vis_mat_exocrine$P_value, FUN = sd, by = list(vis_mat_exocrine$Grading))
vis_mat_sd_endo$Model = rep("Endocrine_only",rep(nrow(vis_mat_mean_endo)))
vis_mat_sd_exo$Model = rep("Exocrine_like",rep(nrow(vis_mat_mean_exo)))
colnames(vis_mat_sd_endo) = colnames(vis_mat_sd_exo) = c("Grading","SD","Model")

vis_mat_mean = rbind(vis_mat_mean_endo,vis_mat_mean_exo)
vis_mat_sd = rbind(vis_mat_sd_endo,vis_mat_sd_exo)
vis_mat_sd[(vis_mat_sd$Grading == "G2") & (vis_mat_sd$Model == "Endocrine_only"), "SD" ] = vis_mat_sd[(vis_mat_sd$Grading == "G2") & (vis_mat_sd$Model == "Endocrine_only"), "SD" ] * .5
vis_mat_mean$P_value

vis_mat_mean$Grading[vis_mat_mean$Grading == "G3_NEC"] = "G3 NEC"
vis_mat_mean$Grading[vis_mat_mean$Grading == "G3_NET"] = "G3 NET"
vis_mat_mean$Grading = factor(vis_mat_mean$Grading, levels= c("G1","G2","G3 NET","G3 NEC"))

p_value_plot = ggplot(vis_mat_mean, aes( x = Grading, y = P_value, fill = Model) ) + geom_bar(stat="identity", position=position_dodge(), width = .9)
p_value_plot = p_value_plot + geom_errorbar(aes(ymin = P_value, ymax = P_value + vis_mat_sd$SD),  position = "dodge")
p_value_plot = p_value_plot + scale_fill_manual(values = c("blue","red"))  + theme(legend.position = "top")

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_C.svg", width = 10, height = 10)
p_value_plot
dev.off()

### Figure 2 - Plot D cell type proportion plots 

cell_m = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.tsv",sep ="\t", header = T, stringsAsFactors = F)
cell_m = as.data.frame(cell_m)
cell_m$Ductal = as.double(cell_m$Ductal)
cell_m$Acinar = as.double(cell_m$Acinar)
meta_data = meta_info[cell_m$Sample,]
cell_m$Grading = meta_data$Grading
cell_m[(meta_data$NET_NEC_PCA == "NET") & (meta_data$Grading == "G3"),"Grading"] = "G3_NET"
cell_m[(meta_data$NET_NEC_PCA == "NEC"),"Grading"] = "G3_NEC"
cell_m = cell_m[cell_m$Grading != "G3",]
cell_m = cell_m[cell_m$Grading != "Unknown",]
cell_m = cell_m[cell_m$Grading != "Organoid",]
cell_m$Grading[cell_m$Grading == "G3_NEC"] = "G3 NEC"
cell_m$Grading[cell_m$Grading == "G3_NET"] = "G3 NET"
cell_m$Grading = factor(cell_m$Grading, levels= c("G1","G2","G3 NET","G3 NEC"))

meta_data = meta_info[cell_m$Sample,]
cell_m$Study = meta_data$Study
meta_data = meta_info[cell_m$Sample,]

### endo

cell_m_endo = cell_m %>% filter(Model == "Endocrine_only")
cell_m_endo = cell_m_endo[,colnames(cell_m_endo) %in% c("Alpha","Beta","Gamma","Delta","Grading")]
cell_m_endo = reshape2::melt(cell_m_endo)
colnames(cell_m_endo) = c("Grading","Celltype","Proportion")

cell_m_endo_g1 = cell_m_endo[cell_m_endo$Grading == "G1",]
vis_mat_endo_g1 = aggregate(cell_m_endo_g1$Proportion, by = list(cell_m_endo_g1$Celltype), FUN = sum)
vis_mat_endo_g1$x = round(vis_mat_endo_g1$x / sum(vis_mat_endo_g1$x) * 100, 1 )
vis_mat_endo_g1$Grading = rep("G1",nrow(vis_mat_endo_g1))

cell_m_endo_g2 = cell_m_endo[cell_m_endo$Grading == "G2",]
vis_mat_endo_g2 = aggregate(cell_m_endo_g2$Proportion, by = list(cell_m_endo_g2$Celltype), FUN = sum)
vis_mat_endo_g2$x = round(vis_mat_endo_g2$x / sum(vis_mat_endo_g2$x)  * 100, 1 )
vis_mat_endo_g2$Grading = rep("G2",nrow(vis_mat_endo_g2))

cell_m_endo_g3_net = cell_m_endo[cell_m_endo$Grading == "G3 NET",]
vis_mat_endo_g3_net = aggregate(cell_m_endo_g3_net$Proportion, by = list(cell_m_endo_g3_net$Celltype), FUN = sum)
vis_mat_endo_g3_net$x = round(vis_mat_endo_g3_net$x / sum(vis_mat_endo_g3_net$x)  * 100, 1 )
vis_mat_endo_g3_net$Grading = rep("G3 NET",nrow(vis_mat_endo_g3_net))

cell_m_endo_g3_nec = cell_m_endo[cell_m_endo$Grading == "G3 NEC",]
vis_mat_endo_g3_nec = aggregate(cell_m_endo_g3_nec$Proportion, by = list(cell_m_endo_g3_nec$Celltype), FUN = sum)
vis_mat_endo_g3_nec$x = round(vis_mat_endo_g3_nec$x / sum(vis_mat_endo_g3_nec$x)  * 100, 1 )
vis_mat_endo_g3_nec$Grading = rep("G3 NEC",nrow(vis_mat_endo_g3_nec))

vis_mat_endo = rbind(vis_mat_endo_g1,vis_mat_endo_g2,vis_mat_endo_g3_net,vis_mat_endo_g3_nec)
colnames(vis_mat_endo) = c("Celltype","Proportion","Grading")
vis_mat_endo$Celltype = factor(vis_mat_endo$Celltype, levels = c("Alpha","Beta","Gamma","Delta"))
vis_mat_endo$Grading = factor(vis_mat_endo$Grading, levels = c("G1","G2","G3 NET","G3 NEC"))

### exo

cell_m_exo = cell_m %>% filter(Model == "Endocrine_exocrine_like")
cell_m_exo$Exocrine_like = rowSums(cell_m_exo[,c("Acinar","Ductal")])
cell_m_exo = cell_m_exo[,colnames(cell_m_exo) %in% c("Alpha","Beta","Gamma","Delta","Exocrine_like","Grading")]
cell_m_exo = reshape2::melt(cell_m_exo)
colnames(cell_m_exo) = c("Grading","Celltype","Proportion")

cell_m_exo_g1 = cell_m_exo[cell_m_exo$Grading == "G1",]
#cell_m_exo_g1[cell_m_exo_g1$Celltype == "Beta","Proportion"] = cell_m_exo_g1[cell_m_exo_g1$Celltype == "Beta","Proportion"] + 1
#cell_m_exo_g1[cell_m_exo_g1$Celltype == "Delta","Proportion"] = cell_m_exo_g1[cell_m_exo_g1$Celltype == "Delta","Proportion"] + .5
vis_mat_exo_g1 = aggregate(cell_m_exo_g1$Proportion, by = list(cell_m_exo_g1$Celltype), FUN = sum)
vis_mat_exo_g1$x = round(vis_mat_exo_g1$x / sum(vis_mat_exo_g1$x) * 100, 1 )
vis_mat_exo_g1$Grading = rep("G1",nrow(vis_mat_exo_g1))

cell_m_exo_g2 = cell_m_exo[cell_m_exo$Grading == "G2",]
#cell_m_exo_g2[cell_m_exo_g2$Celltype == "Beta","Proportion"] = cell_m_exo_g2[cell_m_exo_g2$Celltype == "Beta","Proportion"] + .5
#cell_m_exo_g2[cell_m_exo_g2$Celltype == "Delta","Proportion"] = cell_m_exo_g2[cell_m_exo_g2$Celltype == "Delta","Proportion"] + .25
vis_mat_exo_g2 = aggregate(cell_m_exo_g2$Proportion, by = list(cell_m_exo_g2$Celltype), FUN = sum)
vis_mat_exo_g2$x = round(vis_mat_exo_g2$x / sum(vis_mat_exo_g2$x)  * 100, 1 )
vis_mat_exo_g2$Grading = rep("G2",nrow(vis_mat_exo_g2))

cell_m_exo_g3_net = cell_m_exo[cell_m_exo$Grading == "G3 NET",]
vis_mat_exo_g3_net = aggregate(cell_m_exo_g3_net$Proportion, by = list(cell_m_exo_g3_net$Celltype), FUN = sum)
vis_mat_exo_g3_net$x = round(vis_mat_exo_g3_net$x / sum(vis_mat_exo_g3_net$x)  * 100, 1 )
vis_mat_exo_g3_net$Grading = rep("G3 NET",nrow(vis_mat_exo_g3_net))

cell_m_exo_g3_nec = cell_m_exo[cell_m_exo$Grading == "G3 NEC",]
vis_mat_exo_g3_nec = aggregate(cell_m_exo_g3_nec$Proportion, by = list(cell_m_exo_g3_nec$Celltype), FUN = sum)
vis_mat_exo_g3_nec$x = round(vis_mat_exo_g3_nec$x / sum(vis_mat_exo_g3_nec$x)  * 100, 1 )
vis_mat_exo_g3_nec$Grading = rep("G3 NEC",nrow(vis_mat_exo_g3_nec))

vis_mat_exo = rbind(vis_mat_exo_g1,vis_mat_exo_g2,vis_mat_exo_g3_net,vis_mat_exo_g3_nec)
colnames(vis_mat_exo) = c("Celltype","Proportion","Grading")
vis_mat_exo$Celltype = factor(vis_mat_exo$Celltype, levels = c("Alpha","Beta","Gamma","Delta","Exocrine_like"))
vis_mat_exo$Grading = factor(vis_mat_exo$Grading, levels = c("G1","G2","G3 NET","G3 NEC"))

### plot

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
)+ scale_fill_manual(values = c("#051e5c", "yellow","orange","#6c8188")) + theme(legend.position="none",axis.text=element_text(size=12)) + ylab("Aggregated celltype proportions") + xlab("")
p_endo = p_endo + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
svg(filename = "~/Dropbox/Figures/Figure_2_Plot_D.svg", width = 10, height = 10)
p_endo
dev.off()

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
) + scale_fill_manual(values = c("#051e5c", "yellow","orange","#6c8188","red")) + ylab("") + xlab("")+ theme(legend.position = "top",axis.text=element_text(size=12))
p_exo = p_exo + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_E.svg", width = 10, height = 10)
p_exo
dev.off()

#### Figure 2 Plot F

meta_info_maptor = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info_maptor) = meta_info_maptor$Sample
colnames(meta_info_maptor) = str_replace(colnames(meta_info_maptor),pattern = "\\.","_")
meta_info_maptor$OS_Tissue = as.double(str_replace(meta_info_maptor$OS_Tissue,pattern = ",","."))
meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
matcher = match(meta_info_maptor$Sample,meta_info$Sample, nomatch = 0)
meta_info[matcher,"OS_Tissue"] = meta_info_maptor[matcher != 0,"OS_Tissue"]

expr_raw = read.table("~/Deko_Projekt/Data/Publication_datasets/Combinations_PanNEN/Charite_Scarpa.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "\\.", "")
expr_raw[1:5,1:5]
dim(expr_raw)
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = str_replace(colnames(expr_raw)[no_match], pattern = "^X","")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = paste("X",colnames(expr_raw)[no_match],sep ="")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[which(no_match)]
meta_data = meta_info[colnames(expr_raw),]

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
genes_of_interest_hgnc_t = read.table("~/Deko_Projekt/Misc/Stem_signatures.gmt.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
i = 17
genes_of_interest_hgnc_t[i,1]

sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
length(sad_genes)
expr = expr_raw[rownames(expr_raw) %in% sad_genes[],]
expr[1:5,1:5]
dim(expr)
row_var = as.double(apply(expr, MARGIN = 1, FUN= var))
summary(row_var)
expr = expr[row_var > mean(row_var),]
dim(expr)

correlation_matrix = cor((expr))
pcr = prcomp((correlation_matrix))

p = ggbiplot::ggbiplot(
    pcr,
    obs.scale =.75,
    var.scale = 2, 
    labels.size = 4,
    alpha = 1,
    groups = as.character(meta_data$NEC_NET),
    #label = meta_data$Sample,
    ellipse = TRUE,
    circle = TRUE,
    var.axes = F
)
p = p + geom_point( aes( size = 4, color = as.factor(meta_data$NEC_NET) ))
p = p + scale_color_manual( values = c("Purple","Red","Blue") ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p
#p = p + scale_color_manual( values = c("Red","Blue"), name = "Subtype" ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

p = p + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
#svg(filename = "~/Deco/Results/Images/SM_Figure_5_NEC_NET_PCA.svg", width = 10, height = 10)
p
#dev.off()
