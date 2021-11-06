library("hrbrthemes")
library("waffle")
library(tidyverse)
library("stringr")
library("ggplot2")
library("dplyr")
library("ggpubr")

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
study_selection = c("Alvarez","Charite","Diedisheim","Master","Missiaglia","Sadanandam","Sato","Scarpa")
meta_data = meta_info[meta_info$Study %in% study_selection,]
meta_data = meta_data[meta_data$NEC_NET %in% c("Ambiguous", "Unknown", "NEC", "NET"),]

table(meta_info$Study)

###

meta_data = meta_info[meta_info$Study %in% study_selection,]
meta_data = meta_data[meta_data$Primary_Metastasis != "Control",]
meta_data = meta_data[meta_data$Primary_Metastasis != "Outlier",]
meta_data = meta_data[meta_data$NEC_NET != "Control",]
meta_data = meta_data[meta_data$Grading != "Control",]
dim(meta_data)


### waffle plots

# grading

grading_data = reshape2::melt(table(meta_data[,c("Grading","Study")]))
colnames(grading_data) = c("Grading","Study","Count")
grading_plot = ggplot(
    grading_data,
    aes(
        fill=Grading,
        values=Count
    )) +
    geom_waffle(color = "white") +
    facet_wrap(~Study, ncol = 8)+#,scales = "free",shrink = TRUE,) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggthemes::scale_fill_tableau(name=NULL) +
    coord_equal() +
    labs(
        title = "Grading"
    ) +
    theme_ipsum_rc(grid="") +
    theme_enhance_waffle() + scale_fill_manual(values = c("darkgreen","yellow","darkred","cyan","gray")) + theme(legend.position="top")

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
    facet_wrap(~Study, ncol = 8)+#,scales = "free",shrink = TRUE,) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggthemes::scale_fill_tableau(name=NULL) +
    coord_equal() +
    labs(
        title = "NEC NET"
    ) +
    theme_ipsum_rc(grid="") +
    theme_enhance_waffle() + scale_fill_manual(values = c("purple","red","blue","cyan","gray")) + theme(legend.position="top")

#svg(filename = "~/Dropbox/Figures/F1_waffle_nec_net.svg", width = 10, height = 10)
nec_net_plot
dev.off()

####


meta_info_maptor = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info_maptor) = meta_info_maptor$Sample
colnames(meta_info_maptor) = str_replace(colnames(meta_info_maptor),pattern = "\\.","_")

meta_info_maptor$OS_Tissue = as.double(str_replace(meta_info_maptor$OS_Tissue,pattern = ",","."))

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

matcher = match(meta_info_maptor$Sample,meta_info$Sample, nomatch = 0)
meta_info[matcher,"OS_Tissue"] = meta_info_maptor[matcher != 0,"OS_Tissue"]

expr_raw = read.table("~/MAPTor_NET/BAMs_new/Groetzinger/",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
#dim(expr_raw)
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = str_replace(colnames(expr_raw)[no_match], pattern = "^X","")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = paste("X",colnames(expr_raw)[no_match],sep ="")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[which(no_match)]
meta_data = meta_info[colnames(expr_raw),]
