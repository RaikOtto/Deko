## Figure 3 Segerstolpe Heatmap

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
