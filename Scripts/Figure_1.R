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
#meta_data = meta_info[meta_info$Study %in% study_selection,]
#meta_data = meta_data[meta_data$NEC_NET %in% c("Ambiguous", "Unknown", "NEC", "NET"),]

###

data = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Baron_exocrine/All.S361.tsv",sep ="\t", header = T, stringsAsFactors = F,row.names = 1)

meta_data = meta_info[ rownames(data),]
#meta_data = meta_data[ grep(meta_data$Sample, pattern = "SCLC",invert = TRUE),]
#meta_data = meta_data[ grep(meta_data$Sample, pattern = "MiNEN",invert = TRUE),]
dim(meta_data)

table(meta_data$Primary_Metastasis)
which(meta_data$Primary_Metastasis == "Unknown")
table(meta_data$Histology_Primary)

### Subplot A Primary Metastasis

study_mat = meta_data[,c("Study","Primary_Metastasis")]
grp = dplyr::group_by( study_mat$Primary_Metastasis, as.factor(study_mat$Study))
vis_mat = table(grp)

vis_mat = reshape2::melt(table(study_mat))
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

### Treemap alternative to plot A, amount of samples per study

vis_mat_plot_a = reshape2::melt(table(meta_data$Study))
colnames(vis_mat_plot_a) = c("Study","Count")

plot_a = ggplot(
    vis_mat_plot_a,
    aes(
        area = Count,
        fill = Study,
        label = Study,
        subgroup = Count
    )
) + geom_treemap(alpha= 1.0)#+ geom_treemap_subgroup_border()
plot_a = plot_a + geom_treemap_text(
    colour = "white",
    place = "center",
    grow = FALSE,
    size = 30,
    min.size= 15)
plot_a = plot_a + geom_treemap_subgroup_text(
    place = "bottom",
    grow = FALSE, 
    colour = "white",
    size = 22,
    min.size = 15)
plot_a = plot_a + scale_fill_manual(values = c("black","darkgreen","#340273","#CD8200","brown","#026A73","#0E0273","#730243"))
plot_a = plot_a + theme(legend.position="none")
plot_a

#svg(filename = "~/Dropbox/Figures/F1_P1_alternative.svg", width = 10, height = 10)
plot_a
dev.off()

### plot B alternative

# Create dataset
data = data.frame(
    #individual=vis_mat$Study,
    individual=vis_mat_plot_a$NEN,
    #group= vis_mat$NEN_type,
    group= vis_mat$Study,
    value=vis_mat$Count
)
data = data %>% arrange(group, value)
# Set a number of 'empty bar' to add at the end of each group
# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse( angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
    group_by(group) %>% 
    summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
plot_b = ggplot(
    data,
    aes(x=as.factor(id),
        y=value,
        fill=group))
plot_b = plot_b + geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5)
#plot_b = plot_b + geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE )
#plot_b = plot_b + geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE )
#plot_b = plot_b + geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE )
#plot_b = plot_b + geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE )
#plot_b = plot_b + annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1)
plot_b = plot_b + geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5)
#plot_b= plot_b + ylim(-100,120) +
plot_b = plot_b + theme_minimal()
plot_b = plot_b + theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
)
plot_b = plot_b + coord_polar()
plot_b = plot_b + geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE )
plot_b = plot_b + geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )
plot_b = plot_b + geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,1,1,0,0,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
plot_b = plot_b + scale_fill_manual(values = c("#C75E40","#500307","#17070C","#52D383","#2F3F49","#FC4C1D","#64E0FD","#1601AE"))
plot_b

###

# Racetrack Plot

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

expr_raw = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Baron_exocrine/NEN/All.S.tsv",sep ="\t", stringsAsFactors =FALSE, row.names = 1 )
meta_data = meta_info[rownames(expr_raw),]

vis_mat = meta_data[,c("Study","Histology_Primary")]
vis_mat = as.data.frame(table(reshape2::melt(vis_mat)))
colnames(vis_mat) = c("Study","Histology_Primary","Count")
spacer = matrix(c("A","B","C","","","",0,0,0), nrow = 3, ncol = 3)
colnames(spacer) = colnames(vis_mat)
vis_mat = rbind(vis_mat,spacer)
vis_mat$Count = as.integer(vis_mat$Count)
vis_mat$Histology_Primary = as.factor(vis_mat$Histology_Primary)

Study_labels = c("Alvarez","Missiaglia", "Diedisheim","Charite","Master","Sato","Sadanandam","Scarpa")

vis_mat$Study = factor(vis_mat$Study,levels = c(
    "A","B","C",rev(Study_labels)
))
Study_labels = c("","","",c("Master","Missiaglia","Sadanandam","Sato","Scarpa","Alvarez","Charite","Diedisheim"),rep("",24))

race_plot = ggplot(
    vis_mat,
    aes(
        x = Study,
        y = Count,
        fill = Histology_Primary
    )
)
race_plot = race_plot + geom_bar(width = 0.9, stat="identity")
race_plot = race_plot + coord_polar(theta = "y") + xlab("") + ylab("")
race_plot = race_plot + ylim(c(0,220))
race_plot = race_plot + scale_fill_manual(values = c("brown","yellow","black","darkgreen","white"))
race_plot = race_plot + geom_text(
    data = vis_mat,
    hjust = 1.1,
    size = 4,
    aes(x = Study, y = 0, label = Study_labels))
race_plot = race_plot + theme_minimal() +
    theme(
        legend.position = "right",
          axis.line = element_blank(),
          axis.text.y = element_blank())


#svg(filename = "~/Dropbox/Figures/F1_P2_alternative.svg", width = 10, height = 10)
race_plot
dev.off()

