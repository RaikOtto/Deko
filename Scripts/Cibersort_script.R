eislet = new("ExpressionSet", exprs = (as.matrix(count_data)))
fData(eislet) = data.frame( subtypes )
pData(eislet) = data.frame( subtypes )

B = bseqsc_basis(
  eislet,
  pancreasMarkers,
  clusters = 'subtypes',
  samples = colnames(exprs(eislet)),
  ct.scale = FALSE
)
plotBasis(B, pancreasMarkers, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')

#maxi = apply( B, FUN = which.max, MARGIN = 1 )
#for( type in names(pancreasMarkers)){
#  cell_type = (eval(paste(type)))
#  genes = names(maxi)[as.integer(maxi) == which(colnames(B) == type)]
#  pancreasMarkers[cell_type] = list(genes)
#}


### RUN VARIANCE SELECTION FIRST

eset = new("ExpressionSet", exprs=as.matrix(bam_data));
#eset = new("ExpressionSet", exprs=as.matrix(cnts));

fit = bseqsc_proportions(eset, B, verbose = TRUE, log = F, absolute = T)

###

res_coeff = t(fit$coefficients)
res_cor   = fit$stats

res_coeff[ is.na(res_coeff) ] = 0.0
res_cor[ is.na(res_cor) ] = 0.0

#res_coeff = res_coeff[order(res_coeff[,"Botton_3"]),]
#plot((res_coeff[,"Botton_3"]))

#res_coeff[ res_cor[,"Correlation"] <= .3, ] = .0
#res_cor[ res_cor[,"Correlation"] <= .3, "Correlation" ] = 0.0

###

meta_match = match( colnames(exprs(eset)), meta_info$Name, nomatch = 0 ) 

meta_info$Deco_cor = rep("", nrow(meta_info))
meta_info$Deco_cor[meta_match]  = as.character(res_cor[,2])
res_coeff_match = match( rownames(res_coeff), meta_info$Name, nomatch = 0 )
res_cor_match = match( rownames(res_cor), meta_info$Name, nomatch = 0 )

maxi = apply( res_coeff, FUN = which.max, MARGIN = 1 )
cell_type = colnames(res_coeff)[maxi]
rownames(meta_info) = meta_info$Name
#cell_type[res_cor[,"Correlation"] <= .3] = "Not_sig"
#cell_type[res_cor[,"P-value"] > .05] = "Not_sig"
table(cell_type)

meta_info$Deco_type = rep("",nrow(meta_info))
meta_info$Deco_type[meta_match] = cell_type
meta_data = meta_info[meta_match,]

meta_data$NEUROG3 = as.double(expr_raw[ which( rownames(expr_raw) == "NEUROG3"),])
meta_data$NEC_NET[meta_data$NEC_NET == ""] = "Unknown"
meta_data$Botton_1 = res_coeff[ match(rownames(meta_data),rownames(res_coeff))  ,"Botton_1"]
meta_data$Botton_3 = res_coeff[ match(rownames(meta_data),rownames(res_coeff))  ,"Botton_3"]
meta_data$Botton = meta_data$Botton_1 + meta_data$Botton_3
#meta_data = meta_data[!is.na(meta_data$OS_Tissue),]
colmax = apply(t(res_coeff),  function(vec){max(as.double(vec))}, MARGIN = 2 )
vis_mat =  apply( t(res_coeff), MARGIN = 2, FUN = function(vec){return(vec/max(vec))} )


pheatmap::pheatmap(
    vis_mat,
  annotation_col = meta_data[c("Deco_type","NEC_NET")],
  annotation_colors = aka3,
  annotation_legend = T,
  treeheight_col = 0,
  show_colnames = T,
  show_rownames = T,
  gaps_row = 20
)

#groups = as.character(unlist(meta_data["Deco_type"]))
groups = as.character(apply(vis_mat, MARGIN = 2, FUN = function(vec){return(rownames(vis_mat)[which.max(vec)])}))
groups[  grep(groups, pattern = "Botton") ] = "Botton"
groups[(log(colSums(t(res_coeff)))) < 0] = "Not_sig"
meta_data$Groups = groups
expr = bam_data
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))
p = ggbiplot::ggbiplot(
  pcr,
  obs.scale = .75,
  #groups = as.character(meta_data$Study),
  groups = meta_data$Groups,
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F#,labels = meta_data$Name
)
NEUROG3 = log2(as.double( meta_data$Botton)+1)**1 + 1
Grading = as.character(meta_data$Grading)
Grading[Grading == ""] = "G0"
p = p + geom_point( aes(colour= groups, size = NEUROG3, shape = Grading ) )
p = p + scale_color_manual( values = c("Blue","Yellow","Black","Purple","Orange","Gray") )
#p = p + scale_color_manual( values = c("Blue","Yellow","Purple","Orange","Gray") )
#p = p + guides( color=guide_legend(title="Study", size=guide_legend(title="MKI67"), shape = guide_legend(title="Grading")))
p

cor_rows = apply(bam_data, MARGIN = 1, FUN = function(vec){cor(vec,meta_data[colnames(bam_data),])})
scale_cor = as.double(scale(colSums(cor_mat)))
cor_rows = apply(bam_data, MARGIN = 1, FUN = function(vec){cor(vec,scale_cor)})
hist(cor_rows)
(sort(cor_rows, decreasing = F))

upper = c(0,table(meta_data$Grading[meta_data$Groups == "Not_sig"]))
lower = table(meta_data$Grading[meta_data$Groups != "Not_sig"])

fisher.test(cbind(upper,lower))

# UMAP

vis_mat = expr_raw[ ,]
vis_mat$Gene = as.factor(rownames(vis_mat))
vis_mat$Botton_3 = meta_data$Botton_3
vis_mat = reshape2::melt(vis_mat )
vis_mat$MEN1 = meta_data$MEN1[match(as.character(vis_mat$variable), rownames(meta_data))]
colnames(vis_mat) = c("Gene","Sample","Expression","MEN1")
#vis_mat["Expression"] = as.double(vis_mat["Expression"])

#vis_mat = vis_mat[ match(vis_mat$Gene, result_t$hgnc_symbol),]

men1_plot = ggplot( vis_mat, aes ( x = Gene,  y = Expression))
men1_plot = men1_plot + geom_boxplot( aes(fill = MEN1))
men1_plot = men1_plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
men1_plot


# Pheatmap

pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[c("Groups","NEC_NET")],
  annotation_colors = aka3,
  annotation_legend = T,
  treeheight_col = 0,
  clustering_method = "average",
  show_colnames = T,
  show_rownames = F,
  gaps_row = 20
)

table(as.data.frame(cbind(meta_data$NEC_NET,meta_data$Groups)))

meta_info$Deco_group = rep("",nrow(meta_info))
meta_info[names(groups),"Deco_group"] = groups
