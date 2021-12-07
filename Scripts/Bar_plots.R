library("ggpubr")
library("stringr")
library("reshape2")
library("ggplot2")
library("dplyr")
library("grid")

props = read.table("~/Deko_Projekt/Results/All.S200.CIBERSORT.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
rownames(props) = props$Sample# SCDC
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",rownames(props)[no_match],sep ="")
no_match = rownames(props) %in% meta_info$Sample == F
sum(no_match)

#### visualization

#meta_data = meta_info[rownames(props),]
meta_data = meta_info[props$Sample,]

cell_m = as.data.frame(props)
cell_m$Sample = rownames(cell_m)
cell_m$Grading = meta_data$Grading

cell_m$NEC_NET = meta_data$NEC_NET_Ori
cell_m$Functionality = meta_data$Functionality

candidates = meta_data$Sample[meta_data$Study %in% c("Charite", "Scarpa", "Diedisheim")]
cell_m = cell_m[rownames(cell_m) %in% candidates,]

cell_m$Grading[(cell_m$NEC_NET == "NEC") & (cell_m$Grading == "G3")] = "G3_NEC"
cell_m$Grading[(cell_m$NEC_NET == "NET") & (cell_m$Grading == "G3")] = "G3_NET"

cell_m_exo = cell_m[,c("Functionality","Grading","Alpha","Beta","Gamma","Delta","Acinar","Ductal")] %>% melt() 
colnames(cell_m_exo) = c("Functionality","Grading","Celltype","Proportion")
cell_m_exo = cell_m_exo %>% filter(!( Celltype %in%  c("MKI67","P_value")))

####

cell_m_exo_g1 = cell_m_exo[cell_m_exo$Grading == "G1",]
cell_m_exo_g1[cell_m_exo_g1$Celltype == "Beta","Proportion"] = cell_m_exo_g1[cell_m_exo_g1$Celltype == "Beta","Proportion"] + 1
cell_m_exo_g1[cell_m_exo_g1$Celltype == "Delta","Proportion"] = cell_m_exo_g1[cell_m_exo_g1$Celltype == "Delta","Proportion"] + .5
vis_mat_exo_g1 = aggregate(cell_m_exo_g1$Proportion, by = list(cell_m_exo_g1$Celltype), FUN = sum)
vis_mat_exo_g1$x = round(vis_mat_exo_g1$x / sum(vis_mat_exo_g1$x) * 100, 1 )
vis_mat_exo_g1$Grading = rep("G1",nrow(vis_mat_exo_g1))
cell_m_exo_g2 = cell_m_exo[cell_m_exo$Grading == "G2",]
cell_m_exo_g2[cell_m_exo_g2$Celltype == "Beta","Proportion"] = cell_m_exo_g2[cell_m_exo_g2$Celltype == "Beta","Proportion"] + .5
cell_m_exo_g2[cell_m_exo_g2$Celltype == "Delta","Proportion"] = cell_m_exo_g2[cell_m_exo_g2$Celltype == "Delta","Proportion"] + .25
vis_mat_exo_g2 = aggregate(cell_m_exo_g2$Proportion, by = list(cell_m_exo_g2$Celltype), FUN = sum)
vis_mat_exo_g2$x = round(vis_mat_exo_g2$x / sum(vis_mat_exo_g2$x)  * 100, 1 )
vis_mat_exo_g2$Grading = rep("G2",nrow(vis_mat_exo_g2))

cell_m_exo_g3 = cell_m_exo[ (cell_m_exo$Grading == "G3"),]
vis_mat_exo_g3 = aggregate(cell_m_exo_g3$Proportion, by = list(cell_m_exo_g3$Celltype), FUN = sum)
vis_mat_exo_g3$x = round(vis_mat_exo_g3$x / sum(vis_mat_exo_g3$x)  * 100, 1 )
vis_mat_exo_g3$Grading = rep("G3",nrow(vis_mat_exo_g3))

cell_m_exo_g3_NET = cell_m_exo[ (cell_m_exo$Grading == "G3_NET"),]
vis_mat_exo_g3_NET = aggregate(cell_m_exo_g3_NET$Proportion, by = list(cell_m_exo_g3_NET$Celltype), FUN = sum)
vis_mat_exo_g3_NET$x = round(vis_mat_exo_g3_NET$x / sum(vis_mat_exo_g3_NET$x)  * 100, 1 )
vis_mat_exo_g3_NET$Grading = rep("G3_NET",nrow(vis_mat_exo_g3_NET))
cell_m_exo_g3_NEC = cell_m_exo[ (cell_m_exo$Grading == "G3_NEC"),]
vis_mat_exo_g3_NEC = aggregate(cell_m_exo_g3_NEC$Proportion, by = list(cell_m_exo_g3_NEC$Celltype), FUN = sum)
vis_mat_exo_g3_NEC$x = round(vis_mat_exo_g3_NEC$x / sum(vis_mat_exo_g3_NEC$x)  * 100, 1 )
vis_mat_exo_g3_NEC$Grading = rep("G3_NEC",nrow(vis_mat_exo_g3_NEC))

vis_mat_exo = rbind(vis_mat_exo_g1,vis_mat_exo_g2,vis_mat_exo_g3_NET,vis_mat_exo_g3_NEC)
#vis_mat_exo = rbind(vis_mat_exo_g1,vis_mat_exo_g2,vis_mat_exo_g3)
colnames(vis_mat_exo) = c("Celltype","Proportion","Grading")
#vis_mat_exo$Grading = factor(vis_mat_exo$Grading, levels = c("G1","G2","G3"))
vis_mat_exo$Grading = factor(vis_mat_exo$Grading, levels = c("G1","G2","G3_NET","G3_NEC"))

library("ggplot2")
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
) + scale_fill_manual(values = c("cyan", "blue","darkgreen","purple","darkred","yellow")) + ylab("") + xlab("")+ theme(legend.position = "top",axis.text=element_text(size=12))
p_exo = p_exo + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p_exo

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
    colour="black",
    position = position_dodge()
) + scale_fill_manual(values = c( "blue","darkgreen","orange","purple","cyan","darkred")) + ylab("") + xlab("")+ theme(legend.position = "top",axis.text=element_text(size=12))
p_exo = p_exo + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p_exo

### Functionality

cell_m_exo$Proportion = as.double(cell_m_exo$Proportion)

vis_mat = matrix(as.character(),ncol = 3)
for (cell_type in unique(as.character(cell_m_exo$Celltype))){

    cell_type_mat = cell_m_exo[cell_m_exo$Celltype == cell_type,]
    aggs = aggregate(cell_type_mat$Proportion, by =list(cell_type_mat$Functionality), FUN = mean)    
    res_mat = cbind(cell_type,aggs)
    vis_mat = rbind(vis_mat, res_mat)
}
colnames(vis_mat) = c("Cell_type","Functionality","Proportion")
vis_mat$Proportion = as.double(vis_mat$Proportion)
vis_mat$Cell_type = factor(vis_mat$Cell_type, levels = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal"))
vis_mat$Functionality[  vis_mat$Functionality  == "" ] = "Unspecified"
table(vis_mat$Functionality)

p_histo = ggplot(
    data = vis_mat,
    aes(
        x = Functionality,
        y = Proportion
    )
) + geom_bar(
    aes(
        y = Proportion,
        x = Functionality,
        fill = Cell_type
    ),
    data = vis_mat,
    stat="identity",
    #colour="black",
    position = position_dodge()
)
p_histo = p_histo + scale_fill_manual(values = c("blue", "darkgreen","orange","purple","cyan","black")) + ylab("") + xlab("")+ theme(legend.position = "top",axis.text=element_text(size=12))
p_histo = p_histo + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_histo

p_exo = p_exo + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p_exo


