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

meta_data$NEUROG3 = as.double(expr_raw[ which( rownames(new_merge) == "NEUROG3"),])
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
groups[!(groups %in% c("Botton"))] = "other"
p = ggbiplot::ggbiplot(
  pcr,
  obs.scale = .75,
  #groups = as.character(meta_data$Study),
  groups = groups,
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F#,labels = meta_data$Name
)
NEUROG3 = log2(as.double( meta_data$NEUROG3)+1)**1
Grading = as.character(meta_data$Grading)
p = p + geom_point( aes(colour= groups, size = NEUROG3, shape = Grading ) )
p = p + scale_color_manual( values = c("Black","Gray") )
#p = p + guides( color=guide_legend(title="Study", size=guide_legend(title="MKI67"), shape = guide_legend(title="Grading")))
p

#write.table(meta_info, "~/Deko/Misc/Deko_table.tsv",sep ="\t", row.names = F, quote = F)

umap_plot = umap::umap(t(expr))
vis_data = as.data.frame(umap_plot$layout)
colnames(vis_data) = c("x","y")

dist_mat = dist((vis_data))

hgnc_gene_name = "NEUROG3"
neurog3 = (as.double(expr_raw[rownames(expr_raw) == hgnc_gene_name,]))
size_indicator = neurog3 + min(neurog3) + 1
p = ggplot2::qplot( x = vis_data$x, y = vis_data$y, size = size_indicator, color = groups)
p = p + ggtitle(hgnc_gene_name) +  xlab("umap x") + ylab("umap y")# + geom_text(aes(label=names(clusterCut)),hjust=0, vjust=0)
p
