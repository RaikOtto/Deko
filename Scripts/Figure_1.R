library("stringr")
library("ggplot2")
library("dplyr")
library("ggpubr")
library("png")
library("grid")
library("ggplot2")
library("magick")

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Califano.S165.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

no_match = props$Sample %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",props$Sample[no_match],sep ="")
no_match = props$Sample %in% meta_info$Sample == F
sum(no_match)

dim(props)
meta_data = meta_info[props$Sample,]

#props = props[(meta_data$Histology == "pancreas") | (meta_data$NEC_NET_Color != "Primary") ,]
props  = props[ meta_data$Histology_Primary == "Pancreatic" ,]
dim(props)
meta_data = meta_info[props$Sample,]

meta_data$Location = meta_data$NEC_NET
props = props[ meta_data$Location != "Outlier" ,]
meta_data = meta_info[props$Sample,]
props$Location = meta_data$NEC_NET

vis_mat = reshape2::melt(props[,c("Location","P_value","model")])
colnames(vis_mat) = c("Location","Model","Variable","P_value")
vis_mat$Location[vis_mat$Location == "liver_met"] = "Metastasis"
vis_mat$Location = factor(vis_mat$Location, levels = c("Primary","Metastasis"))
vis_mat$P_value = as.double(vis_mat$P_value)

vis_mat$Model[ vis_mat$Model == "Alpha_Beta_Gamma_Delta_Baron"] = "Endocrine"
vis_mat$Model[ vis_mat$Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron"] = "Endocrine & Exocrine"

p_values = ggplot(
    data = vis_mat,
    aes(
        x = Model,
        y = P_value,
        fill = Location
    )
) + geom_boxplot()
p_values = p_values + scale_fill_manual(values = c("blue", "red","blue", "red")) + ylab("P-value") + xlab("Model")+ theme(legend.position = "top",axis.text=element_text(size=12))
p_values = p_values + geom_hline(yintercept = 0.05, linetype= "dashed", color = "black")
p_values = p_values + annotate("text", x = 2, y = .055, label = "P-value == 0.05",col = "black",size =6)
p_values = p_values + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p_values

####

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

table(meta_info$Study)

study_selection = c("Alvarez","Charite","Diedisheim","Master","Missiaglia","Sadanandam","Sato","Scarpa")
meta_data = meta_info[meta_info$Study %in% study_selection,]
meta_data = meta_data[meta_data$Primary_Metastasis != "Control",]
meta_data = meta_data[meta_data$Primary_Metastasis != "Outlier",]
meta_data = meta_data[meta_data$NEC_NET != "Control",]
meta_data = meta_data[meta_data$Grading != "Control",]
meta_data = meta_data[ grep(meta_data$Sample, pattern = "SCLC",invert = TRUE),]
meta_data = meta_data[ grep(meta_data$Sample, pattern = "MiNEN",invert = TRUE),]
dim(meta_data)

table(meta_data$Primary_Metastasis)
which(meta_data$Primary_Metastasis == "Unknown")

table(meta_data$Histology_Primary)

dim(meta_data[(meta_data$Histology_Primary == "Pancreatic") & (meta_data$Primary_Metastasis == "Primary"),])

#expr_raw = read.table("~/MAPTor_NET/BAMs_new/RepSet_S57.HGNC.DESeq2.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)

#### Subplot A Primary Metastasis

study_mat = meta_data[,c("Study","Primary_Metastasis")]
grp = group_by(study_mat, Study)
vis_mat = table(grp)
vis_mat = reshape2::melt(vis_mat)
colnames(vis_mat) = c("Study","Primary_Metastasis","Amount")

primary_metastasis_plot = ggplot( 
    data = vis_mat,
    aes( 
        x = Study,
        y = Amount,
        fill = Primary_Metastasis
    )
)
primary_metastasis_plot = primary_metastasis_plot + geom_bar(stat="identity", position=position_dodge())
primary_metastasis_plot = primary_metastasis_plot + theme(axis.text.x = element_text(angle = 45, vjust = .5))
primary_metastasis_plot = primary_metastasis_plot + xlab("Study") + ylab("#")
primary_metastasis_plot = primary_metastasis_plot + scale_fill_manual(values = c("red","darkgreen","black","gray"))
primary_metastasis_plot = primary_metastasis_plot + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))

#svg(filename = "~/Dropbox/Figures/F1_P1.svg", width = 10, height = 10)
primary_metastasis_plot
#dev.off()
### Plot B NEC NET

nec_net_mat = meta_data[,c("Study","NEC_NET")]
grp = group_by(nec_net_mat, Study)
vis_mat = table(grp)
vis_mat = reshape2::melt(vis_mat)
colnames(vis_mat) = c("Study","NEC_NET","Amount")
vis_mat$NEC_NET = factor(vis_mat$NEC_NET, levels = c("NET","NEC","Ambiguous","Unknown","Control"))

NEC_NET_plot = ggplot( 
    data = vis_mat,
    aes( 
        x = Study,
        y = Amount,
        fill = NEC_NET
    )
)
NEC_NET_plot = NEC_NET_plot + geom_bar(stat="identity", position=position_dodge())
NEC_NET_plot = NEC_NET_plot + theme(axis.text.x = element_text(angle = 45, vjust = .5))
NEC_NET_plot = NEC_NET_plot + xlab("Study") + ylab("#")
NEC_NET_plot = NEC_NET_plot + scale_fill_manual(values = c("blue","darkred","purple","gray","yellow"))
NEC_NET_plot = NEC_NET_plot + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
NEC_NET_plot = NEC_NET_plot + ggbreak::scale_y_break(c(60, 200))
NEC_NET_plot = NEC_NET_plot + scale_y_continuous(breaks=c(0,10,20,30,40,50,60,204))

#svg(filename = "~/Dropbox/Figures/F1_P2.svg", width = 10, height = 10)
NEC_NET_plot
#dev.off()

### Plot C Tissue type

tissue_mat = meta_data[,c("Study","Histology_Primary")]
grp = group_by(tissue_mat, Study)
vis_mat = table(grp)
vis_mat = reshape2::melt(vis_mat)
colnames(vis_mat) = c("Study","Tissue","Amount")
#vis_mat$NEC_NET = factor(vis_mat$NEC_NET, levels = c("NET","NEC","Ambiguous","Unknown","Control"))

tissue_plot = ggplot( 
    data = vis_mat,
    aes(
        x = Study,
        y = Amount,
        fill = Tissue
    )
)
tissue_plot = tissue_plot + geom_bar(stat="identity", position=position_dodge())
tissue_plot = tissue_plot + theme(axis.text.x = element_text(angle = 45, vjust = .5))
tissue_plot = tissue_plot + xlab("Study") + ylab("#")
tissue_plot = tissue_plot + scale_fill_manual(values = c("brown","yellow","cyan","darkgreen","blue","black","darkgreen"))
tissue_plot = tissue_plot + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
tissue_plot = tissue_plot + ggbreak::scale_y_break(c(25,60))
tissue_plot = tissue_plot + scale_y_continuous(breaks=c(0,10,20,70,80,90,100))

#svg(filename = "~/Dropbox/Figures/F1_P3.svg", width = 10, height = 10)
tissue_plot
#dev.off()

### Plot D Tissue type Metastass

metastasis_mat = meta_data[,c("Study","Histology_Metastasis")]
table(metastasis_mat$Histology_Metastasis)
metastasis_mat$Histology_Metastasis[metastasis_mat$Histology_Metastasis %in% c("Colorectal","Gastric","Mesenteric","Peritoneal","Pulmonary","Spleeneal")] = "Others"
metastasis_mat = metastasis_mat[metastasis_mat$Histology_Metastasis != "Primary",]

grp = group_by(metastasis_mat, Study)
vis_mat = table(grp)
vis_mat = reshape2::melt(vis_mat)
colnames(vis_mat) = c("Study","Tissue","Amount")

metastasis_mat_plot = ggplot( 
    data = vis_mat,
    aes(
        x = Study,
        y = Amount,
        fill = Tissue
    )
)
metastasis_mat_plot = metastasis_mat_plot + geom_bar(stat="identity", position=position_dodge())
metastasis_mat_plot = metastasis_mat_plot + theme(axis.text.x = element_text(angle = 45, vjust = .5))
metastasis_mat_plot = metastasis_mat_plot + xlab("Study") + ylab("#")
metastasis_mat_plot = metastasis_mat_plot + scale_fill_manual(values = c("purple","cyan","darkgreen","orange","black","gray"))
metastasis_mat_plot = metastasis_mat_plot + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
metastasis_mat_plot = metastasis_mat_plot + ggbreak::scale_y_break(c(30,65))
metastasis_mat_plot = metastasis_mat_plot + scale_y_continuous(breaks=c(0,10,20,30,69))

#svg(filename = "~/Dropbox/Figures/F1_P4.svg", width = 10, height = 10)
metastasis_mat_plot
#dev.off()

### Plot E Grading

grading_mat = meta_data[,c("Study","Grading")]
grp = group_by(grading_mat, Study)
vis_mat = table(grp)
vis_mat = reshape2::melt(vis_mat)
colnames(vis_mat) = c("Study","Grading","Amount")
#vis_mat$NEC_NET = factor(vis_mat$NEC_NET, levels = c("NET","NEC","Ambiguous","Unknown","Control"))

grading_plot = ggplot( 
    data = vis_mat,
    aes( 
        x = Study,
        y = Amount,
        fill = Grading
    )
)
grading_plot = grading_plot + geom_bar(stat="identity", position=position_dodge())
grading_plot = grading_plot + theme(axis.text.x = element_text(angle = 45, vjust = .5))
grading_plot = grading_plot + xlab("Study") + ylab("#")
grading_plot = grading_plot + scale_fill_manual(values = c("green","yellow","red","darkgreen","gray"))
grading_plot = grading_plot + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
grading_plot = grading_plot + ggbreak::scale_y_break(c(50, 200))
grading_plot = grading_plot + scale_y_continuous(breaks=c(0,10,20,30,40,50,204))

#svg(filename = "~/Dropbox/Figures/F1_P5.svg", width = 10, height = 10)
grading_plot
#dev.off()
### merge all

#### Plot A Workflow

plot_a_path <- "~/Deko_Projekt/Results/Images/Figure_1_ArtDeco_Concept.png"
plot_a <- readPNG(plot_a_path, native = TRUE, info = TRUE)

workflow_plot= ggplot() +  background_image(plot_a)

p = ggarrange(
    workflow_plot,
    primary_metastasis_plot,
    NEC_NET_plot,
    grading_plot,
    tissue_plot,
    metastasis_mat_plot,
    labels = c("A", "B", "C","D","E","F","G"),
    ncol = 2, nrow = 3)
p

library(grid)
library(gridExtra)

gl <- lapply(1:9,
             function(ii) grobTree(rectGrob(),
                                   textGrob(ii)))

grid.arrange(grobs = gl, layout_matrix = rbind(c(1,1,1,2,3),
                                               c(1,1,1,4,5),
                                               c(6,7,8,9,9)))
gl = c(workflow_plot,
primary_metastasis_plot,
NEC_NET_plot,
grading_plot,
tissue_plot,
metastasis_mat_plot)

library("cowplot")

grid.arrange(
    #workflow_plot,
    primary_metastasis_plot,
    NEC_NET_plot,
    grading_plot,
    tissue_plot,
    metastasis_mat_plot,
    ncol = 2, 
    layout_matrix = cbind(
        c(1,3,5),c(2,4,6)
        #c(1,1,1,2,4,6),c(1,1,1,2,4,6)
        #,c(1,1,1,3,5,7),c(1,1,1,3,5,7)
    ))
####

x = c(1, 2)
y = c(1,9)
df <- as.data.frame(cbind(x, y))
myplot = ggplot(df, aes(x = x, y = y)) + geom_point() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

blankPlot +
    annotation_custom(
        grob = ggplotGrob(workflow_plot),
        xmin = 1,
        xmax = 8,
        ymin = 4,
        ymax = 8
    ) 
+
    annotation_custom(
        grob = rectGrob(gp = gpar(fill = "white")),
        xmin = 7.5,
        xmax = Inf,
        ymin = -Inf,
        ymax = 5
    )
