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

#expr_raw = read.table("~/MAPTor_NET/BAMs_new/RepSet_S103.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
#colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

### Proportion plot

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Visulization_mat.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"
props = props[grep(rownames(props) ,pattern = "PNET04",invert = T),]

no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[no_match] = str_remove_all(rownames(props)[no_match], pattern = "^X")
no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",rownames(props)[no_match],sep ="")
no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[which(no_match)]

dim(props)
#rownames(props) = props$name
meta_data = meta_info[rownames(props),]
props = props[(meta_data$Histology_Primary == "Pancreatic"),]
props = props[props$P_value <= 0.5,]
meta_data = meta_info[rownames(props),]
meta_data = meta_data %>% filter(!(Study %in% c("Sato","Missiaglia")))
props = props[meta_data$Sample,]
props = props[(meta_data$NET_NEC_PCA != "Unknown"),]
meta_data = meta_info[rownames(props),]

vis_mat = props[rownames(meta_data),]
exocrines = as.double(rowSums(vis_mat[,c("Ductal","Acinar")]))
endocrines = as.double(rowSums(vis_mat[,c("Alpha","Beta","Gamma","Delta")]))

selection = c("Alpha","Beta","Gamma","Delta","Metaplastic")
selection = c("Endocrine","Metaplastic")

vis_mat$Endocrine = endocrines
vis_mat$Metaplastic = exocrines
vis_mat = as.data.frame(vis_mat)
vis_mat = vis_mat[,selection]

correlation_matrix = cor(t(vis_mat));pcr = prcomp(t(correlation_matrix))
vis_mat = vis_mat[order(vis_mat$Metaplastic),]
#meta_data$P_value = props[meta_data$Sample,"P_value"]
#meta_data$P_value[props[meta_data$Sample,"P_value"] < 0.05] = "Sig"
#meta_data$P_value[props[meta_data$Sample,"P_value"] >= 0.05] = "Not_sig"

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")

### top plot

upper_plot = pheatmap::pheatmap(
    t(vis_mat),
    annotation_col = meta_data[,c("Grading","NEC_NET","Study")],
    annotation_colors = aka3,
    show_rownames = T,
    show_colnames = F,
    treeheight_row = 0,
    cellheight = 12,
    cluster_rows = F,
    cluster_cols = F,
    legend = 0,
    annotation_legend = TRUE,
    fontsize_row = 14,
    clustering_method = "ward.D2"
)
upper_plot = upper_plot + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="top")
upper_plot

#pannet_cluster[1:which(pannet_cluster == "YY7PXK")] = "Left"

### lower plot

#vis_mat$Sample = factor(rownames(vis_mat), levels = rownames(vis_mat)[order(vis_mat$Metaplastic)])

expr = expr_raw[,rownames(vis_mat)]
dim(expr)

vis_mat_lower = vis_mat
vis_mat_lower$MKi67 = as.double(expr["MKI67",rownames(vis_mat_lower)])
vis_mat_lower$MKi67 = log(vis_mat_lower$MKi67)
vis_mat_lower$MKi67 = vis_mat_lower$MKi67 / max(vis_mat_lower$MKi67)
#vis_mat$Endocrine = rowSums(vis_mat[,c("Alpha","Beta","Gamma","Delta")])
#vis_mat = vis_mat[,!(colnames(vis_mat)  %in% c("Alpha","Beta","Gamma","Delta"))]
vis_mat_lower$Sample = rownames(vis_mat_lower)
vis_mat_lower = vis_mat_lower[order(vis_mat_lower$Metaplastic),]
vis_mat_lower$Sample = factor(1:nrow(vis_mat_lower))
vis_mat_lower$P_value = props[rownames(vis_mat_lower),"P_value"]
vis_mat_lower = vis_mat_lower[,!( colnames(vis_mat_lower) %in% c("Endocrine"))]
vis_mat_lower_melt = reshape2::melt(vis_mat_lower)
colnames(vis_mat_lower_melt) = c("Sample","Measurement","Value")

lower_plot = ggplot(vis_mat_lower_melt, aes(x = as.numeric(Sample), y=Value, color= Measurement, shape = Measurement), ) + geom_point()
lower_plot = lower_plot + geom_smooth(method = "lm", se= FALSE)
lower_plot = lower_plot + scale_color_manual(values=c('#5D6D7E','red', '#453CFA'))
lower_plot = lower_plot +  theme_classic()
lower_plot = lower_plot + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom")
lower_plot = lower_plot + scale_y_continuous( 
    name = "Relative cell-type proportions & p-value",
    sec.axis = sec_axis(~.*8.5, name="log2 MKi67 mRNA expression")
)
lower_plot

summary(log(as.double(expr_raw["MKI67",rownames(vis_mat)])))

### export cohort

vis_mat$Sample = rownames(vis_mat)
#write.table(vis_mat[,c("Sample","Metaplastic")],"~/Downloads/cohorts.tsv",sep ="\t", quote = F,row.names = F)

### UMAP PLOT


props = read.table("~/Deko_Projekt/Results/All.S200.CIBERSORT.tsv",sep ="\t", header = T, as.is=TRUE)
dim(props)
rownames(props ) = props$Sample

no_matcher = which(!( rownames(props) %in% meta_info$Sample))
rownames(props)[no_matcher] = str_replace(rownames(props)[no_matcher], pattern ="^X","")
no_matcher = which(!( rownames(props) %in% meta_info$Sample))
rownames(props)[no_matcher] = str_replace(rownames(props)[no_matcher], pattern ="^","X")
no_matcher = which(!( rownames(props) %in% meta_info$Sample))
no_matcher

meta_data = meta_info[rownames(props),]
props = props[meta_data$Study %in% c("Charite","Scarpa"),]
meta_data = meta_info[rownames(props),]
dim(props)


exocrine = as.double(rowSums(props[,c("Acinar","Ductal")]))
endocrine = as.double(rowSums(props[,c("Alpha","Beta","Delta","Gamma")]))
props$Ratio = log((exocrine+1) / (endocrine+1))
selection = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal","Ratio","Correlation","P_value")
#selection = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal","Differentiation_score.Correlation","Differentiation_score.RMSE","Differentiation_score.P.value")
vis_mat = props[,selection]
#vis_mat = vis_mat[props$P_value <0.05,]

col_vec_nec_net = meta_data$NET_NEC_PCA
col_vec_nec_net[col_vec_nec_net == "NET"] = "blue"
col_vec_nec_net[col_vec_nec_net != "blue"] = "red"

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
custom.config$random_state = 995
custom.config$n_components=2

correlation_matrix = cor(vis_mat)
umap_result = umap::umap(
    vis_mat,
    colvec = col_vec_nec_net,
    preserve.seed = TRUE,
    config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
    umap_result$layout,
    aes(x, y))
umap_p = umap_p + geom_point( aes( size = 4, color = as.character(meta_data$Grading) ))
umap_p = umap_p + theme(legend.position = "none") + xlab("") + ylab("")
umap_p = umap_p + geom_vline( xintercept=-1, size = 2, linetype = 2) + geom_hline( yintercept = -1.25, size = 2, linetype = 2)  
umap_p = umap_p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$NET_NEC_PCA), level=.5, type ="t", size=1.5)
umap_p = umap_p + scale_color_manual( values = c("darkgreen","yellow","red","darkred","#33ACFF")) ##33ACFF ##FF4C33
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

