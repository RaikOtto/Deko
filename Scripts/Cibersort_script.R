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

maxi = apply( B, FUN = which.max, MARGIN = 1 )
for( type in names(pancreasMarkers)){
  cell_type = (eval(paste(type)))
  genes = names(maxi)[as.integer(maxi) == which(colnames(B) == type)]
  pancreasMarkers[cell_type] = list(genes)
}


### run

eset = new("ExpressionSet", exprs=as.matrix(bam_data));
fit = bseqsc_proportions(eset, B, verbose = TRUE)

###

res_coeff = t(fit$coefficients)
res_cor   = fit$stats

res_coeff[ is.na(res_coeff) ] = 0.0
res_cor[ is.na(res_cor) ] = 0.0

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

meta_info$Deco_type = rep("",nrow(meta_info))
meta_info$Deco_type[meta_match] = cell_type
meta_data = meta_info[meta_match,]

meta_data$NEUROG3 = as.double(expr_raw[ which( rownames(expr_raw) == "NEUROG3"),])
meta_data$NEC_NET[meta_data$NEC_NET == ""] = "Unknown"

pheatmap::pheatmap(
  t(res_coeff),
  annotation_col = meta_data[c("Deco_type","NEC_NET")],
  annotation_colors = aka3,
  annotation_legend = T,
  treeheight_col = 0,
  show_colnames = T,
  show_rownames = T,
  gaps_row = 20
)

groups = as.character(unlist(meta_data["Deco_type"]))
groups[  grep(groups, pattern = "Botton") ] = "Botton"
meta_data$Groups = groups
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
NEUROG3 = log2(as.double( meta_data$NEUROG3)+1)**.5 + 1
Grading = as.character(meta_data$Grading)
Grading[Grading == ""] = "G0"
#p = p + geom_point( aes(colour= groups, size = NEUROG3, shape = Grading ) )
p = p + scale_color_manual( values = c("Blue","Yellow","Black","Purple","Orange") )
#p = p + guides( color=guide_legend(title="Study", size=guide_legend(title="MKI67"), shape = guide_legend(title="Grading")))
p

# UMAP

umap_plot = umap::umap(t(expr))
vis_data = as.data.frame(umap_plot$layout)
colnames(vis_data) = c("x","y")
dist_mat = dist((vis_data))
p = ggplot2::qplot( x = vis_data$x, y = vis_data$y, color = meta_data$Groups)
#p = ggplot2::qplot( x = vis_data$x, y = vis_data$y, color = meta_data$Study)
p = p + geom_point( aes(size = NEUROG3, shape = meta_data$Grading ) )
p = p + scale_color_manual( values = c("Blue","Yellow","Black","Purple","Orange") )
p

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

