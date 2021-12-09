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

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.NEN_only.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"
meta_data = meta_info[props$Sample,]

props$Study = meta_data$Study

selection = c("Study","P_value","Model")

vis_mat = props[,selection]
vis_mat$P_value = as.double(vis_mat$P_value)
vis_mat$Model = factor(vis_mat$Model, levels = c("Endocrine_only","Endocrine_exocrine_like"))
vis_mat$Study = factor(vis_mat$Study, levels = c("Alvarez","Charite","Master","Diedisheim","Sato"))
#vis_mat[vis_mat$Study == "Diedisheim","P_value"] =vis_mat[vis_mat$Study == "Diedisheim","P_value"] / 3

p_value_plot = ggplot(vis_mat, aes( x = Study, y = P_value, fill = Model) )
p_value_plot = p_value_plot + geom_boxplot(notch = TRUE,outlier.colour = "red", outlier.shape = 1)
p_value_plot = p_value_plot + ylim(c(0,0.1)) + geom_hline(yintercept = 0.05, color = "red",linetype="dashed", size =2)
p_value_plot = p_value_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
p_value_plot = p_value_plot + scale_fill_manual(values = c("blue","red"))
p_value_plot = p_value_plot + ylab("Mean P-values") + theme(legend.position = "top")
p_value_plot

svg(filename = "~/Dropbox/Figures/Supplementary/SM_Figure_2_Plot_A.svg", width = 10, height = 10)
p_value_plot
dev.off()

# Figure 2 Plot B Metastasis/ Primary/ Organoid

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.NEN_only.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"
meta_data = meta_info[props$Sample,]

props$Primary_Metastasis = meta_data$Primary_Metastasis
props = props %>% filter(!( Primary_Metastasis %in% c("Unknown","Outlier")))

selection = c("Type","Model","Primary_Metastasis","P_value")

vis_mat = props[,selection]
vis_mat$P_value = as.double(vis_mat$P_value)

# nen

vis_mat_nen = vis_mat %>% filter(Type == "NEN")
vis_mat_nen$Model = factor(vis_mat_nen$Model, levels = c("Endocrine_only","Endocrine_exocrine_like"))
vis_mat_nen$Primary_Metastasis = factor(vis_mat_nen$Primary_Metastasis, levels = c("Primary","Metastasis","Organoid"))
p_value_plot = ggplot(vis_mat_nen, aes( x = Primary_Metastasis, y = P_value, fill = Model) )
p_value_plot = p_value_plot + geom_boxplot(notch = TRUE,outlier.colour = "red", outlier.shape = 1)
p_value_plot = p_value_plot + scale_fill_manual(values = c("blue","red"))
p_value_plot = p_value_plot + ylim(c(0,0.1)) + geom_hline(yintercept = 0.05, color = "red",linetype="dashed", size =2)
p_value_plot = p_value_plot + ylab("Mean P-values") + theme(legend.position = "top")
p_value_plot

svg(filename = "~/Dropbox/Figures/Supplementary/Figure_2_Plot_B_nen.svg", width = 10, height = 10)
p_value_plot
dev.off()

### Figure 2 Plot C Grading

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.NEN_only.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"
meta_data = meta_info[props$Sample,]

props$Grading = meta_data$Grading
props$NEC_NET = meta_data$NET_NEC_PCA
props = props %>% filter( Grading %in% c("G1","G2","G3"))
props = props %>% filter( NEC_NET %in% c("NEC","NET"))
meta_data = meta_info[props$Sample,]

props[(props$Grading == "G3") & (props$NEC_NET == "NET") ,"Grading"] = "G3 NET"
props[(props$NEC_NET == "NEC") ,"Grading"] = "G3 NEC"
meta_data = meta_info[props$Sample,]

#props = props[meta_data$Study %in% c("Charite","Scarpa","Master"),]
meta_data = meta_info[props$Sample,]

selection = c("Type","Model","Grading","P_value")
vis_mat = props[,selection]
vis_mat$P_value = as.double(vis_mat$P_value)

# nen

vis_mat_nen = vis_mat[vis_mat$Type == "NEN",]

vis_mat_endocrine = vis_mat_nen[vis_mat_nen$Model == "Endocrine_only",]
vis_mat_exocrine = vis_mat_nen[vis_mat_nen$Model == "Endocrine_exocrine_like",]
table(vis_mat_exocrine$Grading)

vis_mat_mean_endo = aggregate(vis_mat_endocrine$P_value, FUN = mean, by = list(vis_mat_endocrine$Grading))
vis_mat_mean_exo = aggregate(vis_mat_exocrine$P_value, FUN = mean, by = list(vis_mat_exocrine$Grading))
vis_mat_mean_endo$Model = rep("Endocrine_only",rep(nrow(vis_mat_mean_endo)))
vis_mat_mean_exo$Model = rep("Endocrine_exocrine_like",rep(nrow(vis_mat_mean_exo)))
colnames(vis_mat_mean_endo) = colnames(vis_mat_mean_exo) = c("Grading","P_value","Model")

vis_mat_sd_endo = aggregate(vis_mat_endocrine$P_value, FUN = sd, by = list(vis_mat_endocrine$Grading))
vis_mat_sd_exo  = aggregate(vis_mat_exocrine$P_value, FUN = sd, by = list(vis_mat_exocrine$Grading))
vis_mat_sd_endo$Model = rep("Endocrine_only",rep(nrow(vis_mat_mean_endo)))
vis_mat_sd_exo$Model = rep("Endocrine_exocrine_like",rep(nrow(vis_mat_mean_exo)))
colnames(vis_mat_sd_endo) = colnames(vis_mat_sd_exo) = c("Grading","SD","Model")

vis_mat_mean = rbind(vis_mat_mean_endo,vis_mat_mean_exo)
vis_mat_mean$Model = factor( as.character(vis_mat_mean$Model), levels = c("Endocrine_only","Endocrine_exocrine_like"))
vis_mat_mean$Grading = factor(vis_mat_mean$Grading, levels = c("G1","G2", "G3 NET", "G3 NEC"))
vis_mat_sd = rbind(vis_mat_sd_endo,vis_mat_sd_exo)

p_value_plot = ggplot(vis_mat_mean, aes( x = Grading, y = P_value, fill = Model) ) + geom_bar(stat="identity", position=position_dodge(), width = .9)
p_value_plot = p_value_plot + geom_errorbar(aes(ymin = P_value, ymax = P_value + vis_mat_sd$SD),  position = "dodge")
p_value_plot = p_value_plot + scale_fill_manual(values = c("blue","red"))  + theme(legend.position = "top")

svg(filename = "~/Dropbox/Figures/Supplementary/Figure_2_Plot_C_nen.svg", width = 10, height = 10)
p_value_plot+ ylim(c(0,0.1))
dev.off()


### Figure 2 - Plot D cell type proportion plots 

### endo

cell_m = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Absolute/endocrine_only.4_scRNA_studies.absolute.tsv",sep ="\t", header = T, stringsAsFactors = F)
cell_m = as.data.frame(cell_m)
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
cell_m = cell_m[,!( colnames(cell_m) %in% c("P_value","Correlation","RMSE","Sample","Ductal","Acinar","Model") )]

cell_m = reshape2::melt(cell_m)
colnames(cell_m) = c("Dataset","Grading","Cell_type","Proportion")


vis_mat_g1 = cell_m %>% filter(Grading == "G1") %>% group_by(Dataset, Cell_type) %>% summarise( Proportion = sum(Proportion) )
vis_mat_g1$Grading = rep( "G1", nrow(vis_mat_g1))
vis_mat_g2 = cell_m %>% filter(Grading == "G2") %>% group_by(Dataset, Cell_type) %>% summarise( Proportion = sum(Proportion) )
vis_mat_g2$Grading = rep( "G2", nrow(vis_mat_g2))
vis_mat_g3_net = cell_m %>% filter(Grading == "G3 NET") %>% group_by(Dataset, Cell_type) %>% summarise( Proportion = sum(Proportion) )
vis_mat_g3_net$Grading = rep( "G3 NET", nrow(vis_mat_g3_net))
vis_mat_g3_nec = cell_m %>% filter(Grading == "G3 NEC") %>% group_by(Dataset, Cell_type) %>% summarise( Proportion = sum(Proportion) )
vis_mat_g3_nec$Grading = rep( "G3 NEC", nrow(vis_mat_g3_nec))

vis_mat_endo = rbind(vis_mat_g1,vis_mat_g2,vis_mat_g3_net,vis_mat_g3_nec)
vis_mat_endo$Cell_type = factor(vis_mat_endo$Cell_type, levels = c("Alpha","Beta","Gamma","Delta"))
vis_mat_endo$Proportion = as.double(vis_mat_endo$Proportion)

### endo g1

vis_mat_endo_g1 = vis_mat_endo %>% filter(Grading == "G1") 
Proportions = vis_mat_endo_g1 %>% group_by(Dataset) %>% summarise( Proportion / sum(Proportion))
vis_mat_endo_g1$Proportion = as.double(Proportions$`Proportion/sum(Proportion)`)
vis_mat_endo_g1 = as.data.frame(vis_mat_endo_g1)

dummy = c("A","Alpha",0,"G1")
vis_mat_endo_g1 = rbind(vis_mat_endo_g1, dummy, dummy,dummy,dummy,dummy)
vis_mat_endo_g1$Dataset = factor(vis_mat_endo_g1$Dataset, levels = c("A","Lawlor","Segerstolpe","Baron"))
vis_mat_endo_g1$Proportion = as.double(vis_mat_endo_g1$Proportion)

p_endo_g1 = ggplot(vis_mat_endo_g1, aes(x = Dataset, y = Proportion, fill = Cell_type)) 
p_endo_g1 = p_endo_g1 + geom_bar(stat = "identity", width = 1.0, color = "black")
p_endo_g1 = p_endo_g1 + coord_polar("y", start=0,direction =-1) + ylab("") + xlab("")
p_endo_g1 = p_endo_g1 + scale_fill_manual(values = c("#051e5c", "yellow","orange","#6c8188","red")) + theme(legend.position="none",axis.text=element_text(size=12)) + ylab("Celltype proportions") + xlab("")

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_D_1.svg", width = 10, height = 10)
p_endo_g1
dev.off()

### endo g2

vis_mat_endo_g1 = vis_mat_endo %>% filter(Grading == "G2") 
Proportions = vis_mat_endo_g1 %>% group_by(Dataset) %>% summarise( Proportion / sum(Proportion))
vis_mat_endo_g1$Proportion = as.double(Proportions$`Proportion/sum(Proportion)`)
vis_mat_endo_g1 = as.data.frame(vis_mat_endo_g1)

dummy = c("A","Alpha",0,"G2")
vis_mat_endo_g1 = rbind(vis_mat_endo_g1, dummy, dummy,dummy,dummy,dummy)
vis_mat_endo_g1$Dataset = factor(vis_mat_endo_g1$Dataset, levels = c("A","Lawlor","Segerstolpe","Baron"))
vis_mat_endo_g1$Proportion = as.double(vis_mat_endo_g1$Proportion)

p_endo_g1 = ggplot(vis_mat_endo_g1, aes(x = Dataset, y = Proportion, fill = Cell_type)) 
p_endo_g1 = p_endo_g1 + geom_bar(stat = "identity", width = 1.0, color = "black")
p_endo_g1 = p_endo_g1 + coord_polar("y", start=0,direction = -1) + ylab("") + xlab("")
p_endo_g1 = p_endo_g1 + scale_fill_manual(values = c("#051e5c", "yellow","orange","#6c8188","red")) + theme(legend.position="none",axis.text=element_text(size=12)) + ylab("Celltype proportions") + xlab("")

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_D_2.svg", width = 10, height = 10)
p_endo_g1
dev.off()

### endo g3 NET

vis_mat_endo_g1 = vis_mat_endo %>% filter(Grading == "G3 NET") 
Proportions = vis_mat_endo_g1 %>% group_by(Dataset) %>% summarise( Proportion / sum(Proportion))
vis_mat_endo_g1$Proportion = as.double(Proportions$`Proportion/sum(Proportion)`)
vis_mat_endo_g1 = as.data.frame(vis_mat_endo_g1)

dummy = c("A","Alpha",0,"G3 NET")
vis_mat_endo_g1 = rbind(vis_mat_endo_g1, dummy, dummy,dummy,dummy,dummy)
vis_mat_endo_g1$Dataset = factor(vis_mat_endo_g1$Dataset, levels = c("A","Lawlor","Segerstolpe","Baron"))
vis_mat_endo_g1$Proportion = as.double(vis_mat_endo_g1$Proportion)

p_endo_g1 = ggplot(vis_mat_endo_g1, aes(x = Dataset, y = Proportion, fill = Cell_type)) 
p_endo_g1 = p_endo_g1 + geom_bar(stat = "identity", width = 1.0, color = "black")
p_endo_g1 = p_endo_g1 + coord_polar("y", start=0,direction = -1) + ylab("") + xlab("")
p_endo_g1 = p_endo_g1 + scale_fill_manual(values = c("#051e5c", "yellow","orange","#6c8188","red")) + theme(legend.position="none",axis.text=element_text(size=12)) + ylab("Celltype proportions") + xlab("")


#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_D_3.svg", width = 10, height = 10)
p_endo_g1
dev.off()


### endo g3 NEC

vis_mat_endo_g1 = vis_mat_endo %>% filter(Grading == "G3 NEC") 
Proportions = vis_mat_endo_g1 %>% group_by(Dataset) %>% summarise( Proportion / sum(Proportion))
vis_mat_endo_g1$Proportion = as.double(Proportions$`Proportion/sum(Proportion)`)
vis_mat_endo_g1 = as.data.frame(vis_mat_endo_g1)

dummy = c("A","Alpha",0,"G3 NEC")
vis_mat_endo_g1 = rbind(vis_mat_endo_g1, dummy, dummy,dummy,dummy,dummy)
vis_mat_endo_g1$Dataset = factor(vis_mat_endo_g1$Dataset, levels = c("A","Lawlor","Segerstolpe","Baron"))
vis_mat_endo_g1$Proportion = as.double(vis_mat_endo_g1$Proportion)

p_endo_g1 = ggplot(vis_mat_endo_g1, aes(x = Dataset, y = Proportion, fill = Cell_type)) 
p_endo_g1 = p_endo_g1 + geom_bar(stat = "identity", width = 1.0, color = "black")
p_endo_g1 = p_endo_g1 + coord_polar("y", start=0,direction = -1) + ylab("") + xlab("")
p_endo_g1 = p_endo_g1 + scale_fill_manual(values = c("#051e5c", "yellow","orange","#6c8188","red")) + theme(legend.position="none",axis.text=element_text(size=12)) + ylab("Celltype proportions") + xlab("")
p_endo_g1

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_D_4.svg", width = 10, height = 10)
p_endo_g1
dev.off()

### exo

cell_m = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Absolute/endocrine_exocrine.4_scRNA_studies.absolute.tsv",sep ="\t", header = T, stringsAsFactors = F)
cell_m = as.data.frame(cell_m)
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
cell_m = cell_m[,!( colnames(cell_m) %in% c("P_value","Correlation","RMSE","Sample") )]

cell_m = reshape2::melt(cell_m)
colnames(cell_m) = c("Dataset","Grading","Cell_type","Proportion")

vis_mat_g1 = cell_m %>% filter(Grading == "G1") %>% group_by(Dataset, Cell_type) %>% summarise( Mean_Proportion = sum(Proportion) )
vis_mat_g1$Grading = rep( "G1", nrow(vis_mat_g1))
vis_mat_g2 = cell_m %>% filter(Grading == "G2") %>% group_by(Dataset, Cell_type) %>% summarise( Mean_Proportion = sum(Proportion) )
vis_mat_g2$Grading = rep( "G2", nrow(vis_mat_g2))
vis_mat_g3_net = cell_m %>% filter(Grading == "G3 NET") %>% group_by(Dataset, Cell_type) %>% summarise( Mean_Proportion = sum(Proportion) )
vis_mat_g3_net$Grading = rep( "G3 NET", nrow(vis_mat_g3_net))
vis_mat_g3_nec = cell_m %>% filter(Grading == "G3 NEC") %>% group_by(Dataset, Cell_type) %>% summarise( Mean_Proportion = sum(Proportion) )
vis_mat_g3_nec$Grading = rep( "G3 NEC", nrow(vis_mat_g3_nec))

vis_mat_exo = rbind(vis_mat_g1,vis_mat_g2,vis_mat_g3_net,vis_mat_g3_nec)
colnames(blocker) = colnames(vis_mat_exo) = c("Dataset","Cell_type","Proportion","Grading")
vis_mat_exo$Cell_type = factor(vis_mat_exo$Cell_type, levels = c("Exocrine_like","Alpha","Beta","Gamma","Delta"))
vis_mat_exo$Proportion = as.double(vis_mat_exo$Proportion)

### exo g1

vis_mat_exo_g1 = vis_mat_exo %>% filter(Grading == "G1") 
Proportions = vis_mat_exo_g1 %>% group_by(Dataset) %>% summarise( Proportion / sum(Proportion))
vis_mat_exo_g1$Proportion = as.double(Proportions$`Proportion/sum(Proportion)`)
vis_mat_exo_g1 = as.data.frame(vis_mat_exo_g1)

dummy = c("A","Alpha",0,"G1")
vis_mat_exo_g1 = rbind(vis_mat_exo_g1, dummy, dummy,dummy,dummy,dummy)
vis_mat_exo_g1$Dataset = factor(vis_mat_exo_g1$Dataset, levels = c("A","Lawlor","Segerstolpe","Baron"))
#vis_mat_exo_g1$Cell_type = factor(vis_mat_exo_g1$Cell_type, levels = c("Exocrine_like","Alpha","Beta","Gamma","Delta"))
vis_mat_exo_g1$Proportion = as.double(vis_mat_exo_g1$Proportion)

p_exo_g1 = ggplot(vis_mat_exo_g1, aes(x = Dataset, y = Proportion, fill = Cell_type)) 
p_exo_g1 = p_exo_g1 + geom_bar(stat = "identity", width = 1.0, color = "black")
p_exo_g1 = p_exo_g1 + coord_polar("y", start=0,direction = -1) + ylab("") + xlab("")
p_exo_g1 = p_exo_g1 + scale_fill_manual(values = c("red","#051e5c", "yellow","orange","#6c8188")) + theme(legend.position="none",axis.text=element_text(size=12)) + ylab("Celltype proportions") + xlab("")

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_D_5.svg", width = 10, height = 10)
p_exo_g1
dev.off()

### exo g2

vis_mat_exo_g1 = vis_mat_exo %>% filter(Grading == "G2") 
Proportions = vis_mat_exo_g1 %>% group_by(Dataset) %>% summarise( Proportion / sum(Proportion))
vis_mat_exo_g1$Proportion = as.double(Proportions$`Proportion/sum(Proportion)`)
vis_mat_exo_g1 = as.data.frame(vis_mat_exo_g1)

dummy = c("A","Alpha",0,"G2")
vis_mat_exo_g1 = rbind(vis_mat_exo_g1, dummy, dummy,dummy,dummy,dummy)
vis_mat_exo_g1$Dataset = factor(vis_mat_exo_g1$Dataset, levels = c("A","Lawlor","Segerstolpe","Baron"))
vis_mat_exo_g1$Proportion = as.double(vis_mat_exo_g1$Proportion)

p_exo_g1 = ggplot(vis_mat_exo_g1, aes(x = Dataset, y = Proportion, fill = Cell_type)) 
p_exo_g1 = p_exo_g1 + geom_bar(stat = "identity", width = 1.0, color = "black")
p_exo_g1 = p_exo_g1 + coord_polar("y", start=0,direction = -1) + ylab("") + xlab("")
p_exo_g1 = p_exo_g1 + scale_fill_manual(values = c("red","#051e5c", "yellow","orange","#6c8188")) + theme(legend.position="none",axis.text=element_text(size=12)) + ylab("Celltype proportions") + xlab("")
p_exo_g1

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_D_6.svg", width = 10, height = 10)
p_exo_g1
dev.off()

### exo g3 NET

vis_mat_exo_g1 = vis_mat_exo %>% filter(Grading == "G3 NET") 
Proportions = vis_mat_exo_g1 %>% group_by(Dataset) %>% summarise( Proportion / sum(Proportion))
vis_mat_exo_g1$Proportion = as.double(Proportions$`Proportion/sum(Proportion)`)
vis_mat_exo_g1 = as.data.frame(vis_mat_exo_g1)

dummy = c("A","Alpha",0,"G3 NET")
vis_mat_exo_g1 = rbind(vis_mat_exo_g1, dummy, dummy,dummy,dummy,dummy)
vis_mat_exo_g1$Dataset = factor(vis_mat_exo_g1$Dataset, levels = c("A","Lawlor","Segerstolpe","Baron"))
vis_mat_exo_g1$Proportion = as.double(vis_mat_exo_g1$Proportion)

p_exo_g1 = ggplot(vis_mat_exo_g1, aes(x = Dataset, y = Proportion, fill = Cell_type)) 
p_exo_g1 = p_exo_g1 + geom_bar(stat = "identity", width = 1.0, color = "black")
p_exo_g1 = p_exo_g1 + coord_polar("y", start=0,direction = -1) + ylab("") + xlab("")
p_exo_g1 = p_exo_g1 + scale_fill_manual(values = c("red","#051e5c", "yellow","orange","#6c8188")) + theme(legend.position="none",axis.text=element_text(size=12)) + ylab("Celltype proportions") + xlab("")
p_exo_g1

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_D_7.svg", width = 10, height = 10)
p_exo_g1
dev.off()


### exo g3 NEC

vis_mat_exo_g1 = vis_mat_exo %>% filter(Grading == "G3 NEC") 
Proportions = vis_mat_exo_g1 %>% group_by(Dataset) %>% summarise( Proportion / sum(Proportion))
vis_mat_exo_g1$Proportion = as.double(Proportions$`Proportion/sum(Proportion)`)
vis_mat_exo_g1 = as.data.frame(vis_mat_exo_g1)

dummy = c("A","Alpha",0,"G3 NEC")
vis_mat_exo_g1 = rbind(vis_mat_exo_g1, dummy, dummy,dummy,dummy,dummy)
vis_mat_exo_g1$Dataset = factor(vis_mat_exo_g1$Dataset, levels = c("A","Lawlor","Segerstolpe","Baron"))
vis_mat_exo_g1$Proportion = as.double(vis_mat_exo_g1$Proportion)

p_exo_g1 = ggplot(vis_mat_exo_g1, aes(x = Dataset, y = Proportion, fill = Cell_type)) 
p_exo_g1 = p_exo_g1 + geom_bar(stat = "identity", width = 1.0, color = "black")
p_exo_g1 = p_exo_g1 + coord_polar("y", start=0,direction = -1) + ylab("") + xlab("")
p_exo_g1 = p_exo_g1 + scale_fill_manual(values = c("red","#051e5c", "yellow","orange","#6c8188")) + theme(legend.position="none",axis.text=element_text(size=12)) + ylab("Celltype proportions") + xlab("")
p_exo_g1

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_D_8.svg", width = 10, height = 10)
p_exo_g1
dev.off()
