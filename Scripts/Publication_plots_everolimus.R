library("ggpubr")
library("stringr")
library("reshape2")
library("ggplot2")
library("dplyr")
library("grid")

#443 407 444 579 450 409 452 PNET08

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

#meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

expr_raw = read.table("~/MAPTor_NET/BAMs_new/RepSet_S57.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/MAPTor_NET/BAMs_new/RepSet_S57.HGNC.DESeq2.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = paste("X",colnames(expr_raw)[no_match],sep ="")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
no_match
meta_data = meta_info[colnames(expr_raw),]

#meta_info$SUV39H1 = rep("",nrow(meta_info))
#meta_info$SUV39H2 = rep("",nrow(meta_info))
#matcher = match(rownames(meta_info), colnames(expr_raw),nomatch = 0)
#meta_info[matcher != 0, "SUV39H1"] = as.double(expr_raw[grep(rownames(expr_raw),pattern = "SUV39H1", value = F),matcher])
#meta_info[matcher != 0, "SUV39H2"] = as.double(expr_raw[grep(rownames(expr_raw),pattern = "SUV39H2", value = F),matcher])

#"132502" %in% colnames(expr_raw)

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
#genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SeneSys_gene_sets.tsv",sep ="\t", stringsAsFactors = F, header = F)
#genes_of_interest_hgnc_t = read.table("~/Deko_Projekt/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t = read.table("~/MAPTor_NET//Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1

liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]
i = 64
genes_of_interest_hgnc_t[i,1]

sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
#sad_genes = sad_genes[!(sad_genes %in% liver_genes)]
length(sad_genes)

expr_raw_normalized = matrix(as.double(as.character(unlist(expr_raw))), ncol = ncol(expr_raw));
expr_raw_normalized = apply(expr_raw_normalized, FUN =scale, MARGIN = 2)
colnames(expr_raw_normalized) = colnames(expr_raw);
rownames(expr_raw_normalized) = rownames(expr_raw)

expr = expr_raw_normalized[rownames(expr_raw_normalized) %in% sad_genes,meta_data$NEC_NET %in% "NET"]
#expr = expr_raw[rownames(expr_raw) %in% sad_genes,]
expr[1:5,1:5]
dim(expr)

###

correlation_matrix = cor((expr))
pcr = prcomp(t(correlation_matrix))

#svg(filename = "~/Downloads/Heatmap.svg", width = 10, height = 10)
p  =pheatmap::pheatmap(
  correlation_matrix,
  #expr,
  annotation_col = meta_data[,c("NEC_NET","Grading","Study")],
  #annotation_col = meta_data[,c("Grading","NEC_NET_Color","P_value","Study")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = F,
  treeheight_row = 0,
  legend = T,
  fontsize_col = 7,
  clustering_method = "average"
)

###

variance_vec = sort(apply(expr,FUN = var, MARGIN = 1), decreasing = T)
vis_mat_variance = melt(variance_vec)
vis_mat_variance = cbind(rownames(vis_mat_variance),vis_mat_variance)
colnames(vis_mat_variance) = c("Gene","Variance")
vis_mat_variance$Gene = factor(vis_mat_variance$Gene, levels = vis_mat_variance$Gene)

ggplot(
  data = vis_mat_variance[1:10,],
  aes(
    x = Gene,
    y = Variance
  )
) + geom_bar(
  stat="identity"
) + theme(axis.text.x=element_text(angle = 90, hjust = 0))

meta_data = meta_info[colnames(expr),]

vis_mat_mean = melt(expr)
vis_mat_mean = cbind(vis_mat_mean,meta_data[vis_mat_mean$Var2,"NEC_NET"])
colnames(vis_mat_mean) = c("Gene","Sample","Expression","NEC_NET")

vis_mat_mean = vis_mat_mean[vis_mat_mean$Gene %in% as.character(vis_mat_variance[1:10,1]),]
vis_mat_mean$Gene = factor(vis_mat_mean$Gene, levels = vis_mat_variance[1:10,1])

ggplot(
    data = vis_mat_mean,
    aes(
        x = Gene,
        y = Expression
    )
) + geom_boxplot(
#    stat="dodge",
    aes(fill = NEC_NET)
) + theme(axis.text.x=element_text(angle = 90, hjust = 0)) + scale_fill_manual( values = c("red","blue"))

library("umap")

col_vec_nec_net = meta_data$NEC_NET
col_vec_nec_net[col_vec_nec_net == "NET"] = "blue"
col_vec_nec_net[col_vec_nec_net != "blue"] = "red"

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
custom.config$random_state = 350
custom.config$n_components=2

umap_result = umap::umap(
  correlation_matrix,
  colvec = col_vec_nec_net,
  preserve.seed = FALSE,
  config=custom.config
  )

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

#plotly::plot_ly(
#  x=umap_result$layout$x,
#  y=umap_result$layout$y,
#  z=umap_result$layout$z,
#  type="scatter3d",
#  mode="markers",
  #name = rownames(umap_result$layout),
#  color=meta_data$Study  )

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point( aes( size = 4, color = as.factor(meta_data$Grading) ))
umap_p = umap_p + theme(legend.position = "none") + xlab("") + ylab("")
umap_p = umap_p + geom_vline( xintercept=0, size = 2, linetype = 2) + geom_hline( yintercept = 0, size = 2, linetype = 2)  
umap_p = umap_p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$NEC_NET), level=.5, type ="t", size=1.5)
umap_p = umap_p + scale_color_manual( values = c("darkgreen","yellow","red","darkred","#33ACFF")) ##33ACFF ##FF4C33
umap_p = umap_p + annotate("text", x = 3.5, y = -3.5, label = "NEC",col = "darkred",size =8)
umap_p = umap_p + annotate("text", x = -1, y = 2.5, label = "NET",col = "#33ACFF",size =8)
umap_p
custom.config$random_state # 188 #350

library("gridExtra")

grading_map = as.data.frame(cbind(as.double(umap_result$layout[,1]),meta_data$Grading))
colnames(grading_map) = c("x","Grading")
grading_map$x = as.double(as.character(unlist(grading_map$x)))

xdensity = ggplot(
  grading_map,
  aes(x))
xdensity = xdensity + geom_density(alpha=.5,aes(fill=as.factor(meta_data$Grading)), )
xdensity = xdensity + geom_vline(xintercept = 0.35, size = 2, linetype= 2)
xdensity = xdensity + theme(legend.position = "none") + xlab("") + ylab("")
xdensity = xdensity + guides(fill=guide_legend(title = "Grading"))
xdensity = xdensity + theme(axis.title.x = element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
xdensity = xdensity + scale_fill_manual( values = c("darkgreen","yellow","red")) ##33ACFF ##FF4C33
xdensity = xdensity + annotate("text", x = -1.6, y = .25, label = "G1",size = 8) + annotate("text", x = -3.1, y = .4, label = "G2",size = 8) + annotate("text", x = 3.25, y = .15, label = "G3",size= 8)

#scale_fill_manual(values = c('#999999','#E69F00')) + 

nec_net_map = as.data.frame(cbind(as.double(umap_result$layout[,2]),meta_data$NEC_NET))
colnames(grading_map) = c("y","NEC_NET")
grading_map$y = as.double(as.character(unlist(grading_map$y)))

grading_map$y = grading_map$y

ydensity = ggplot(
    grading_map,
    aes(y))
ydensity = ydensity + geom_density(alpha=.5,aes(fill=as.factor(meta_data$NEC_NET)), )
ydensity = ydensity + scale_fill_manual(values=c("#33ACFF","darkred"))
ydensity = ydensity + theme(legend.position = "none") + xlab("") + ylab("")
ydensity = ydensity + guides(fill=guide_legend(title="NEC_NET"))
ydensity = ydensity + geom_vline(xintercept = 8.325,size= 2, linetype= 2)
ydensity = ydensity + annotate("text", x = 5.2, y = 0.2, label = "NEC",size = 8) + annotate("text", x = 11.5, y = 0.2, label = "NET",size = 8)
ydensity = ydensity + coord_flip()  + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())

blankPlot = ggplot() + geom_blank(aes(10,10)) + annotate("text", x = 1, y = 1, label = "") + annotate("text", x = 10, y = 10, label = "") +theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())
blankPlot = blankPlot + annotate("text", x = 3.5, y = 8, label = "Grading",col = "black",size =8) + annotate("text", x = 7, y = 2, label = "Subtype",col = "black",size =8) + annotate("segment", x = 1, xend = 10, y = 1, yend = 10)

svg(filename = "~/Downloads/Umap.svg", width = 10, height = 10)
grid.arrange(
  xdensity,
  blankPlot,
  umap_p,
  ydensity,
  ncol=2,
  nrow=2,
  widths=c(4, 1.4), heights=c(1.4, 4)
)
dev.off()
