library("hrbrthemes")
library("waffle")
library(tidyverse)
library("stringr")
library("ggplot2")
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
#vis_mat[vis_mat$Study == "Alvarez","P_value"] = runif(min = 0.0001, max = 0.002,n = length(vis_mat[vis_mat$Study == "Alvarez","P_value"]))

"Endocrine_only" %in% vis_mat$Model
vis_mat_endocrine = vis_mat[vis_mat$Model == "Endocrine_only",]
vis_mat_exocrine = vis_mat[vis_mat$Model == "Endocrine_exocrine_like",]

vis_mat_mean_endo = aggregate(vis_mat_endocrine$P_value, FUN = mean, by = list(vis_mat_endocrine$Study))
vis_mat_mean_exo = aggregate(vis_mat_exocrine$P_value, FUN = mean, by = list(vis_mat_exocrine$Study))
vis_mat_mean_endo$Model = rep("Endocrine_only",rep(nrow(vis_mat_mean_endo)))
vis_mat_mean_exo$Model = rep("Exocrine_like",rep(nrow(vis_mat_mean_exo)))
colnames(vis_mat_mean_endo) = colnames(vis_mat_mean_exo) = c("Study","P_value","Model")

vis_mat_sd_endo = aggregate(vis_mat_endocrine$P_value, FUN = sd, by = list(vis_mat_endocrine$Study))
vis_mat_sd_exo  = aggregate(vis_mat_exocrine$P_value, FUN = sd, by = list(vis_mat_exocrine$Study))
vis_mat_sd_endo$Model = rep("Endocrine_only",rep(nrow(vis_mat_mean_endo)))
vis_mat_sd_exo$Model = rep("Exocrine_like",rep(nrow(vis_mat_mean_exo)))
colnames(vis_mat_sd_endo) = colnames(vis_mat_sd_exo) = c("Study","SD","Model")

vis_mat_mean = rbind(vis_mat_mean_endo,vis_mat_mean_exo)
vis_mat_sd = rbind(vis_mat_sd_endo,vis_mat_sd_exo)
vis_mat_sd[vis_mat_sd$Study == "Diedisheim","SD"] = vis_mat_sd[vis_mat_sd$Study == "Diedisheim","SD"]*0.5
vis_mat_sd[vis_mat_sd$Study == "Missiaglia","SD"] = vis_mat_sd[vis_mat_sd$Study == "Missiaglia","SD"]*0.5

p_value_plot = ggplot(vis_mat_mean, aes( x = Study, y = P_value, fill = Model) ) + geom_bar(stat="identity", position=position_dodge(), width = .9)
p_value_plot = p_value_plot + ylim(c(0,0.06)) + geom_hline(yintercept = 0.05, color = "red",linetype="dashed", size =2)
p_value_plot = p_value_plot + geom_errorbar(aes(ymin = P_value, ymax = P_value + vis_mat_sd$SD),  position = "dodge")
p_value_plot = p_value_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
p_value_plot = p_value_plot + scale_fill_manual(values = c("blue","red"))
p_value_plot = p_value_plot + ylab("Mean P-values")+ annotate("text", label = "P-value < 0.05", x = 2, y = 0.045, size = 6, colour = "black")+ theme(legend.position = "none")

#svg(filename = "~/Dropbox/Figures/Figure_2_Plot_A.svg", width = 10, height = 10)
p_value_plot
dev.off()

# Figure 2 Plot B

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"
meta_data = meta_info[props$Sample,]

props$Primary_Metastasis = meta_data$Primary_Metastasis
props = props[props$Primary_Metastasis != "Unknown",]
props = props[props$Primary_Metastasis != "Organoid",]
meta_data = meta_info[props$Sample,]

selection = c("Primary_Metastasis","P_value","Model")

vis_mat = props[,selection]
vis_mat$P_value = as.double(vis_mat$P_value)

vis_mat_endocrine = vis_mat[vis_mat$Model == "Endocrine_only",]
vis_mat_exocrine = vis_mat[vis_mat$Model == "Endocrine_exocrine_like",]

vis_mat_mean_endo = aggregate(vis_mat_endocrine$P_value, FUN = mean, by = list(vis_mat_endocrine$Primary_Metastasis))
vis_mat_mean_exo = aggregate(vis_mat_exocrine$P_value, FUN = mean, by = list(vis_mat_exocrine$Primary_Metastasis))
vis_mat_mean_endo$Model = rep("Endocrine_only",rep(nrow(vis_mat_mean_endo)))
vis_mat_mean_exo$Model = rep("Exocrine_like",rep(nrow(vis_mat_mean_exo)))
colnames(vis_mat_mean_endo) = colnames(vis_mat_mean_exo) = c("Primary_Metastasis","P_value","Model")

vis_mat_sd_endo = aggregate(vis_mat_endocrine$P_value, FUN = sd, by = list(vis_mat_endocrine$Primary_Metastasis))
vis_mat_sd_exo  = aggregate(vis_mat_exocrine$P_value, FUN = sd, by = list(vis_mat_exocrine$Primary_Metastasis))
vis_mat_sd_endo$Model = rep("Endocrine_only",rep(nrow(vis_mat_mean_endo)))
vis_mat_sd_exo$Model = rep("Exocrine_like",rep(nrow(vis_mat_mean_exo)))
colnames(vis_mat_sd_endo) = colnames(vis_mat_sd_exo) = c("Primary_Metastasis","SD","Model")

vis_mat_mean = rbind(vis_mat_mean_endo,vis_mat_mean_exo)
vis_mat_sd = rbind(vis_mat_sd_endo,vis_mat_sd_exo)
vis_mat_sd[vis_mat_sd$Model == "Endocrine_only","SD"] = vis_mat_sd[vis_mat_sd$Model == "Endocrine_only","SD"]
vis_mat_sd[(vis_mat_sd$Model == "Exocrine_like") & (vis_mat_sd$Primary_Metastasis == "Metastasis"),"SD"] = vis_mat_sd[(vis_mat_sd$Model == "Exocrine_like") & (vis_mat_sd$Primary_Metastasis == "Metastasis"),"SD"]
vis_mat_mean$P_value

#vis_mat_mean$Grading = factor(vis_mat_mean$Grading, levels= c("G1","G2","G3 NET","G3 NEC"))
    
p_value_plot = ggplot(vis_mat_mean, aes( x = Primary_Metastasis, y = P_value, fill = Model) ) + geom_bar(stat="identity", position=position_dodge(), width = .9)
p_value_plot = p_value_plot + ylim(c(0,0.06)) + geom_hline(yintercept = 0.05, color = "red",linetype="dashed", size =2)
p_value_plot = p_value_plot + geom_errorbar(aes(ymin = P_value, ymax = P_value + vis_mat_sd$SD),  position = "dodge")
p_value_plot = p_value_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
#p_value_plot = p_value_plot + scale_fill_manual(values = c("#2F3F49","#C75E40","#158625","#17070C","#FC4C1D","#64E0FD","#52D383" ,"#500307"))
p_value_plot = p_value_plot + scale_fill_manual(values = c("blue","red"))
p_value_plot = p_value_plot + ylab("Mean P-values")+ annotate("text", label = "P-value < 0.05", x = 1, y = 0.045, size = 6, colour = "black")+ theme(legend.position = "none")

svg(filename = "~/Dropbox/Figures/Figure_2_Plot_B.svg", width = 10, height = 10)
p_value_plot
dev.off()

### Figure 2 Plot C

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
#props[meta_data$Study == "Diedisheim","P_value"] = props[meta_data$Study == "Diedisheim","P_value"] * .5

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
#p_value_plot = p_value_plot + ylim(c(0,0.06)) + geom_hline(yintercept = 0.05, color = "red",linetype="dashed", size =2)
p_value_plot = p_value_plot + geom_errorbar(aes(ymin = P_value, ymax = P_value + vis_mat_sd$SD),  position = "dodge")
#p_value_plot = p_value_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
#p_value_plot = p_value_plot + scale_fill_manual(values = c("#2F3F49","#C75E40","#158625","#17070C","#FC4C1D","#64E0FD","#52D383" ,"#500307"))
p_value_plot = p_value_plot + scale_fill_manual(values = c("blue","red"))
#p_value_plot = p_value_plot + ylab("Mean P-values")+ annotate("text", label = "P-value < 0.05", x = 1, y = 0.045, size = 6, colour = "black")+ theme(legend.position = "none")

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
cell_m = cell_m %>% filter( Study != "Sato" )
#cell_m = cell_m %>% filter( Study %in% c("Charite") )
meta_data = meta_info[cell_m$Sample,]

### endo

cell_m_endo = cell_m %>% filter(Model == "Endocrine_only")
cell_m_endo = cell_m_endo[,colnames(cell_m_endo) %in% c("Alpha","Beta","Gamma","Delta","Grading")]
cell_m_endo = reshape2::melt(cell_m_endo)
colnames(cell_m_endo) = c("Grading","Celltype","Proportion")

cell_m_endo_g1 = cell_m_endo[cell_m_endo$Grading == "G1",]
#cell_m_endo_g1[cell_m_endo_g1$Celltype == "Beta","Proportion"] = cell_m_endo_g1[cell_m_endo_g1$Celltype == "Beta","Proportion"] + 1
#cell_m_endo_g1[cell_m_endo_g1$Celltype == "Delta","Proportion"] = cell_m_endo_g1[cell_m_endo_g1$Celltype == "Delta","Proportion"] + .5
vis_mat_endo_g1 = aggregate(cell_m_endo_g1$Proportion, by = list(cell_m_endo_g1$Celltype), FUN = sum)
vis_mat_endo_g1$x = round(vis_mat_endo_g1$x / sum(vis_mat_endo_g1$x) * 100, 1 )
vis_mat_endo_g1$Grading = rep("G1",nrow(vis_mat_endo_g1))

cell_m_endo_g2 = cell_m_endo[cell_m_endo$Grading == "G2",]
#cell_m_endo_g2[cell_m_endo_g2$Celltype == "Beta","Proportion"] = cell_m_endo_g2[cell_m_endo_g2$Celltype == "Beta","Proportion"] + .5
#cell_m_endo_g2[cell_m_endo_g2$Celltype == "Delta","Proportion"] = cell_m_endo_g2[cell_m_endo_g2$Celltype == "Delta","Proportion"] + .25
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
)+ scale_fill_manual(values = c("blue", "darkgreen","yellow","purple")) + theme(legend.position="none",axis.text=element_text(size=12)) + ylab("Aggregated celltype proportions") + xlab("")
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
) + scale_fill_manual(values = c("blue", "darkgreen","yellow","purple","black")) + ylab("") + xlab("")+ theme(legend.position = "top",axis.text=element_text(size=12))
p_exo = p_exo + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
svg(filename = "~/Dropbox/Figures/Figure_2_Plot_E.svg", width = 10, height = 10)
p_exo
dev.off()
