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

data = read.table("~/Deko_Projekt/Results/cell",sep ="\t", header = T, stringsAsFactors = F,row.names = 1)
meta_data = meta_info[ rownames(data),]

### Figure 1 plot A - Racetrack Plot

#data = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Relative/Baron_exocrine/NEN/All.S.tsv",sep ="\t", header = T, stringsAsFactors = F,row.names = 1)
meta_data = meta_info
table(meta_data$Study)
meta_data = meta_data %>% 
    filter(!(Primary_Metastasis %in% c("Outlier","Control")))%>% filter(Study != "Fadista") %>% filter(NET_NEC_PCA != "MiNEN") %>% filter(Site_of_metastasis != "Control")
dim(meta_data)
table(meta_data$Study)


vis_mat = meta_data[,c("Study","Site_of_primary")]
vis_mat = as.data.frame(table(reshape2::melt(vis_mat)))
colnames(vis_mat) = c("Study","Site_of_primary","Count")
spacer = matrix(c("A","B","C","","","",0,0,0), nrow = 3, ncol = 3)
colnames(spacer) = colnames(vis_mat)
vis_mat = rbind(vis_mat,spacer)
vis_mat$Count = as.integer(vis_mat$Count)
vis_mat$Site_of_primary = as.factor(vis_mat$Site_of_primary)

Study_labels = c("Califano","Missiaglia", "Diedisheim","Riemer","Fröhling","Sato","Sadanandam","Scarpa")

vis_mat$Study = factor(vis_mat$Study,levels = c(
    "A","B","C",rev(Study_labels)
))
Study_labels = c("","","",c("Fröhling","Missiaglia","Sadanandam","Sato","Scarpa","Califano","Riemer","Diedisheim"),rep("",40))
vis_mat$Site_of_primary = factor(vis_mat$Site_of_primary, levels = rev(c("Pancreatic","Small_intestinal","Large_intestinal","Gastric/duodenal","Other")))

race_plot = ggplot(
    vis_mat,
    aes(
        x = Study,
        y = Count,
        fill = Site_of_primary
    )
)
race_plot = race_plot + geom_bar(width = 0.6, stat="identity")
race_plot = race_plot + coord_polar(theta = "y") + xlab("") + ylab("")
race_plot = race_plot + ylim(c(0,230))
race_plot = race_plot + scale_fill_manual(values = c("gray","#C83C3C","#FEC10A","#78B8B4","#32506E","white"))
race_plot = race_plot + geom_text(
    data = vis_mat,
    hjust = 1.1,
    size = 6,
    aes(x = Study, y = 0, label = Study_labels))
race_plot = race_plot + theme_minimal() +
    theme(
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.y = element_blank())


#svg(filename = "~/Dropbox/Figures/F1_P2_alternative.svg", width = 10, height = 10)
race_plot
dev.off()

### Figure 1 Plot b - Treemap alternative to plot A, amount of samples per study

meta_data_pannen = meta_data %>% filter(Site_of_primary == "Pancreatic")
dim(meta_data_pannen)
table(meta_data_pannen$Study)
vis_mat_plot_a = reshape2::melt(table(meta_data_pannen$Study))
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
plot_a = plot_a + scale_fill_manual(values = c("#78B8B4","#E15759","#AFB41E","#DC5AB5","#5F87B4","#F58E2D","#533C9D","#AFC3D7"))
plot_a = plot_a + theme(legend.position="none")
plot_a

#svg(filename = "~/Dropbox/Figures/F1_P1_alternative.svg", width = 10, height = 10)
plot_a
dev.off()

### Figure 1 Plot C - waffle plots

# Grading

grading_data = reshape2::melt(table(meta_data[,c("Grading","Study")]))
colnames(grading_data) = c("Grading","Study","Count")
grading_plot = ggplot(
    grading_data,
    aes(
        fill=Grading,
        values=Count
    )) +
    geom_waffle(color = "white") +
    facet_wrap(~Study, ncol = 3)+#,scales = "free",shrink = TRUE,) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggthemes::scale_fill_tableau(name=NULL) +
    coord_equal() +
    labs(
        title = "Grading"
    ) +
    theme_ipsum_rc(grid="") +
    theme_enhance_waffle() + scale_fill_manual(values = c("#D2E6E6","#78B8B4","#F54C19","gray")) + theme(legend.position="top")

#svg(filename = "~/Dropbox/Figures/F1_waffle_grading.svg", width = 10, height = 10)
grading_plot
dev.off()

# NEC / NET

nec_net_data = reshape2::melt(table(meta_data[,c("NEC_NET","Study")]))
colnames(nec_net_data) = c("NEC_NET","Study","Count")
nec_net_data$NEC_NET = as.character(nec_net_data$NEC_NET)
nec_net_data[as.character(nec_net_data$NEC_NET) %in% c("MiNEN","SCLC"),"NEC_NET"] = "Other"
nec_net_plot = ggplot(
    nec_net_data,
    aes(
        fill=NEC_NET,
        values=Count
    )) +
    geom_waffle(color = "white") +
    facet_wrap(~Study, ncol = 3)+#,scales = "free",shrink = TRUE,) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggthemes::scale_fill_tableau(name=NULL) +
    coord_equal() +
    labs(
        title = "NEC NET"
    ) +
    theme_ipsum_rc(grid="") +
    theme_enhance_waffle() + scale_fill_manual(values = c("#9664AF","#C00000","#1E64A5","#C35078","gray")) + theme(legend.position="top")

#svg(filename = "~/Dropbox/Figures/F1_waffle_nec_net.svg", width = 10, height = 10)
nec_net_plot
dev.off()

# Primary metastasis

primary_metastasis_data = reshape2::melt(table(meta_data[,c("Primary_Metastasis","Study")]))
colnames(primary_metastasis_data) = c("primary_metastasis","Study","Count")
primary_metastasis_data$primary_metastasis = factor(primary_metastasis_data$primary_metastasis, levels = c(
    "Primary","Metastasis","Organoid","Unknown"
))

primary_metastasis_plot = ggplot(
    primary_metastasis_data,
    aes(
        fill=primary_metastasis,
        values=Count
    )) +
    geom_waffle(color = "white") +
    facet_wrap(~Study, ncol = 3)+#,scales = "free",shrink = TRUE,) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggthemes::scale_fill_tableau(name=NULL) +
    coord_equal() +
    labs(
        title = "Primary vs. metastasis"
    ) +
    theme_ipsum_rc(grid="") +
    theme_enhance_waffle() + scale_fill_manual(values = c("#DCD13C","#533C9D","#C35078","gray")) + theme(legend.position="top")

#svg(filename = "~/Dropbox/Figures/F1_waffle_primary_metastasis.svg", width = 10, height = 10)
primary_metastasis_plot
dev.off()

#### Plot D Workflow

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
