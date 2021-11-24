library("stringr")
library("reshape2")
library("dplyr")
library("umap")
library("ggplot2")
library("ggpubr")
library("grid")

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

# Figure 3 Plot A Diedisheim

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Relative/Baron_exocrine//Diedisheim.S62.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",rownames(props)[no_match],sep ="")
no_match = rownames(props) %in% meta_info$Sample == F
sum(no_match)

dim(props)
meta_data = meta_info[rownames(props),]

props = props[meta_data$Functionality %in% c("Insulinoma","Glucagonoma", "Somatostatinoma", "PPoma", "Non-Functional","Unknown","VIPoma","ACTH","Gastrinoma"),]
meta_data = meta_info[rownames(props),]
table(meta_data$NEC_NET)

selection = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal")
exocrines = as.double(rowSums(props[,c("Ductal","Acinar")]))
endocrines = as.double(rowSums(props[,c("Alpha","Beta","Gamma","Delta")]))

meta_data$Ratio = log((exocrines+.1) / (endocrines+.1))

vis_mat = props[,selection]
vis_mat$Exocrine_like = vis_mat$Acinar + vis_mat$Ductal
vis_mat = vis_mat[,!(colnames(vis_mat) %in% c("Acinar","Ductal"))]

correlation_matrix = cor(t(vis_mat))

pcr = prcomp(t(correlation_matrix))

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
svg(filename = "~/Dropbox/Figures/Figure_3_Plot_A.svg", width = 10, height = 10)
Diedisheim_plot = pheatmap::pheatmap(
    t(vis_mat),
    annotation_col = meta_data[,c("Grading","Cluster","DFS","Functionality","NEC_NET")],
    annotation_colors = aka3,
    show_rownames = T,
    show_colnames = FALSE,
    treeheight_row = 0,
    cellheight = 20,
    cluster_rows = FALSE,
    legend = T,
    fontsize_row = 14,
    #annotation_legend = FALSE,
    clustering_method = "single"
)
dev.off()

### Figure 3 Plot B Charite + Master +Scarpa absolutplot G3 only

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
props_endo = props %>% filter(Model == "Endocrine_only")
props_exo = props %>% filter(Model == "Endocrine_exocrine_like")
rownames(props_exo) = props_exo$Sample

no_matcher = which(!( props$Sample %in% meta_info$Sample))
props$Sample[no_matcher] = str_replace(props$Sample[no_matcher], pattern ="^X","")
no_matcher = which(!( props$Sample %in% meta_info$Sample))
props$Sample[no_matcher] = str_replace(props$Sample[no_matcher], pattern ="^","X")
no_matcher = which(!( props$Sample %in% meta_info$Sample))
no_matcher

meta_data = meta_info[props$Sample,]
props = props[meta_data$Study %in% c("Charite","Scarpa","Master"),]
meta_data = meta_info[rownames(props),]
dim(props)

props_exo$Exocrine_like = as.double(rowSums(props_exo[,c("Acinar","Ductal")]))

selection = c("Alpha","Beta","Gamma","Delta","Exocrine_like")
vis_mat_exo = props_exo[,selection]
meta_data = meta_info[rownames(vis_mat_exo),]
vis_mat_exo = vis_mat_exo%>% filter(meta_data$Study %in% c("Charite","Master","Scarpa"))
meta_data = meta_info[rownames(vis_mat_exo),]
vis_mat_exo = vis_mat_exo%>% filter(meta_data$Grading %in% c("G3"))
meta_data = meta_info[rownames(vis_mat_exo),]

vis_mat_exo$MKi67 = log(as.double(meta_data$Mki_67))
vis_mat_exo$MKi67 = vis_mat_exo$MKi67 / max(vis_mat_exo$MKi67)
#vis_mat_exo$MKi67 = vis_mat_exo$MKi67 * max(vis_mat_exo$Exocrine_like)

vis_mat_balanced = vis_mat_exo# / row_max
correlation_matrix = cor(t(vis_mat_balanced))

vis_mat_balanced = vis_mat_balanced[order(vis_mat_balanced$Exocrine_like),]
source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")

#svg(filename = "~/Dropbox/Figures/Figure_3_Plot_B.svg", width = 10, height = 10)
upper_plot = pheatmap::pheatmap(
    t(vis_mat_balanced),
    #correlation_matrix,
    annotation_col = meta_data[,c("Grading","Functionality","NEC_NET","Study")],
    annotation_colors = aka3,
    show_rownames = TRUE,
    show_colnames = FALSE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    legend = T,
    clustering_method = "complete",
    cellheight = 35,
    fontsize = 10,
    annotation_legend = FALSE
)
dev.off()

### Figure 3 C MKi-67 plot

samples = rownames(vis_mat_balanced)[order(vis_mat_balanced$Exocrine_like)]
meta_data = meta_info[samples,]
meta_data$Mki_67

#row_max = as.double(apply( vis_mat_exo, MARGIN = 1 , FUN = max))
#vis_mat_balanced = vis_mat_exo / row_max

vis_mat_exo$Endocrine = as.double(rowSums(vis_mat_exo[,c("Alpha","Beta","Delta","Gamma")]))
vis_mat_exo$MKi67 = log(as.double(meta_data$Mki_67))

selection = c("Endocrine","Exocrine_like","MKi67")
vis_mat = cbind( 1:nrow(vis_mat_exo) ,vis_mat_exo[samples,selection])
colnames(vis_mat) = c("Sample","Endocrine","Exocrine_like","MKi67")
vis_mat$MKi67 = vis_mat$MKi67 / max(vis_mat$MKi67)
#vis_mat$Endocrine = vis_mat$Endocrine / max(vis_mat$Endocrine)
#vis_mat$Exocrine_like = vis_mat$Exocrine_like / max(vis_mat$Exocrine_like)
vis_mat$Sample = factor(vis_mat$Sample)
vis_mat = reshape2::melt(vis_mat)
colnames(vis_mat) = c("Sample","Property","Value")
vis_mat$Property = as.factor(vis_mat$Property)
vis_mat$Value = as.double(vis_mat$Value)

lower_plot = ggplot(vis_mat, aes(x = Sample, y = Value, color = Property, group = Property) ) + geom_point(stat='identity', position='identity', aes(color=Property))
lower_plot = lower_plot + geom_smooth(method = "lm", se= FALSE)
lower_plot = lower_plot + scale_color_manual(values=c('blue','black', 'red')) +  theme_classic()
lower_plot = lower_plot + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="top")
lower_plot = lower_plot + scale_y_continuous(
    name = "Relative cell-type proportions & p-value",
    sec.axis = sec_axis(~.*8.5, name="log2 MKi67 mRNA expression")
)

#svg(filename = "~/Dropbox/Figures/Figure_3_Plot_c.svg", width = 10, height = 10)
lower_plot
dev.off()

cor.test (vis_mat_exo$Exocrine_like, vis_mat_exo$MKi67)

### Figure 3 D correlation  plot

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
meta_data = meta_info[props$Sample,]
props = props[meta_data$Study %in% c("Charite","Scarpa","Master","Diedisheim"),]
meta_data = meta_info[rownames(props),]
dim(props)

props_exo = props %>% filter(Model == "Endocrine_exocrine_like")
rownames(props_exo) = props_exo$Sample

no_matcher = which(!( props$Sample %in% meta_info$Sample))
props$Sample[no_matcher] = str_replace(props$Sample[no_matcher], pattern ="^X","")
no_matcher = which(!( props$Sample %in% meta_info$Sample))
props$Sample[no_matcher] = str_replace(props$Sample[no_matcher], pattern ="^","X")
no_matcher = which(!( props$Sample %in% meta_info$Sample))
no_matcher

selection = c("Alpha","Beta","Gamma","Delta","Ductal","Acinar","P_value","Correlation")
vis_mat_exo = props_exo[,selection]
meta_data = meta_info[rownames(vis_mat_exo),]
correlation_matrix = cor(t(vis_mat_exo))

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
svg(filename = "~/Dropbox/Figures/Figure_3_Plot_E.svg", width = 10, height = 10)
upper_plot = pheatmap::pheatmap(
    #t(vis_mat_balanced),
    correlation_matrix,
    annotation_col = meta_data[,c("Grading","Functionality","NEC_NET","Study")],
    annotation_colors = aka3,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_row = 0,
    legend = T,
    clustering_method = "complete",
    #cellheight = 35,
    fontsize = 10,
    annotation_legend = TRUE
)
dev.off()

### Figure 3 Plot E UMAP PLOT

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.tsv",sep ="\t", header = T, as.is=TRUE)
props = props%>% filter(Model == "Endocrine_exocrine_like")
#props = props%>% filter(Model == "Endocrine_only")
rownames(props) = props$Sample
meta_data = meta_info[rownames(props),]
props = props%>% filter(meta_data$Study %in% c("Charite","Scarpa"))
dim(props)

no_matcher = which(!( rownames(props) %in% meta_info$Sample))
rownames(props)[no_matcher] = str_replace(rownames(props)[no_matcher], pattern ="^X","")
no_matcher = which(!( rownames(props) %in% meta_info$Sample))
rownames(props)[no_matcher] = str_replace(rownames(props)[no_matcher], pattern ="^","X")
no_matcher = which(!( rownames(props) %in% meta_info$Sample))
no_matcher

colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

props$Exocrine_like = as.double(rowSums(props[,c("Acinar","Ductal")]))
props$Endocrine = as.double(rowSums(props[,c("Alpha","Beta","Delta","Gamma")]))

selection = c("Alpha","Beta","Gamma","Delta","Ductal","Acinar","Correlation","P_value")
vis_mat = props[,selection]

meta_data = meta_info[rownames(vis_mat),]
vis_mat = vis_mat[meta_data$NET_NEC_PCA != "Unknown",]
meta_data = meta_info[rownames(vis_mat),]
col_vec_nec_net = meta_data$NET_NEC_PCA
col_vec_nec_net[col_vec_nec_net == "NET"] = "blue"
col_vec_nec_net[col_vec_nec_net == "NEC"] = "red"
col_vec_nec_net[col_vec_nec_net == "Unknown"] = "gray"

correlation_matrix = cor(vis_mat)

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
#custom.config$random_state = 995
custom.config$n_components=2

umap_result = umap::umap(
    vis_mat,
    colvec = col_vec_nec_net,
    #colvec = meta_data$Study,
    preserve.seed = TRUE,
    config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
    umap_result$layout,
    aes(x, y))
umap_p = umap_p + geom_point( aes( size = 4, color = as.character(meta_data$Grading) ))
#umap_p = umap_p + geom_point( aes( size = 4, color = as.character(meta_data$Grading) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$NET_NEC_PCA), level=.5, type ="t", size=1.5)
#umap_p = umap_p + geom_point( aes( size = 4, color = as.character(meta_data$Study) ))
umap_p = umap_p + scale_color_manual( values = c("darkgreen","yellow","darkred","red","blue")) ##33ACFF ##FF4C33

umap_p = umap_p + theme(legend.position = "none") + xlab("") + ylab("")
umap_p = umap_p + geom_vline( xintercept=-1, size = 2, linetype = 2) + geom_hline( yintercept = -1.25, size = 2, linetype = 2)  
umap_p = umap_p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())

#umap_p = umap_p + annotate("text", x = 3.5, y = -3.5, label = "NEC",col = "darkred",size =8)
#umap_p = umap_p + annotate("text", x = -1, y = 2.5, label = "NET",col = "#33ACFF",size =8)
umap_p
custom.config$random_state

# 995 RepSet96 c("Alpha","Beta","Gamma","Delta","Acinar","Ductal","Ratio","Correlation")

library("gridExtra")

grading_map = as.data.frame(cbind(as.double(umap_result$layout[,1]),meta_data$Grading))
colnames(grading_map) = c("x","Grading")
grading_map$x = as.double(as.character(unlist(grading_map$x)))

xdensity = ggplot(
    grading_map,
    aes(x))
xdensity = xdensity + geom_density(alpha=.5,aes(fill=as.factor(meta_data$Grading)), )
xdensity = xdensity + geom_vline(xintercept = -1, size = 2, linetype= 2)
xdensity = xdensity + theme(legend.position = "none") + xlab("") + ylab("")
xdensity = xdensity + guides(fill=guide_legend(title = "Grading"))
xdensity = xdensity + theme(axis.title.x = element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
xdensity = xdensity + scale_fill_manual( values = c("darkgreen","yellow","red")) ##33ACFF ##FF4C33
xdensity = xdensity + annotate("text", x = -3.5, y = .35, label = "G3",size = 8) + annotate("text", x = 3, y = .35, label = "G2",size = 8) + annotate("text", x = 1.35, y = .35, label = "G1",size= 8)

#scale_fill_manual(values = c('#999999','#E69F00')) + 

nec_net_map = as.data.frame(cbind(as.double(umap_result$layout[,2]),meta_data$NET_NEC_PCA))
colnames(grading_map) = c("y","NEC_NET")
grading_map$y = as.double(as.character(unlist(grading_map$y)))

grading_map$y = grading_map$y

ydensity = ggplot(
    grading_map,
    aes(y))
ydensity = ydensity + geom_density(alpha=.5,aes(fill=as.factor(meta_data$NET_NEC_PCA)), )
ydensity = ydensity + scale_fill_manual(values=c("darkred","#33ACFF"))
ydensity = ydensity + theme(legend.position = "none") + xlab("") + ylab("")
ydensity = ydensity + guides(fill=guide_legend(title="NEC_NET"))
ydensity = ydensity + geom_vline(xintercept = -1.25,size= 2, linetype= 2)
ydensity = ydensity + annotate("text", x = -3.55, y = 0.16, label = "NEC",size = 8) + annotate("text", x = 1.65, y = 0.16, label = "NET",size = 8)
ydensity = ydensity + coord_flip()  + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())

blankPlot = ggplot() + geom_blank(aes(10,10)) + annotate("text", x = 1, y = 1, label = "") + annotate("text", x = 10, y = 10, label = "") +theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())
blankPlot = blankPlot + annotate("text", x = 3.5, y = 8, label = "Grading",col = "black",size =6) + annotate("text", x = 7, y = 2, label = "Subtype",col = "black",size =6) + annotate("segment", x = 1, xend = 10, y = 1, yend = 10)

#svg(filename = "~/Downloads/Umap.svg", width = 10, height = 10)
grid.arrange(
    xdensity,
    blankPlot,
    umap_p,
    ydensity,
    ncol=2,
    nrow=2,
    widths=c(4, 1.4), heights=c(1.4, 4)
)
#dev.off()


