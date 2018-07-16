
table(subtypes)

eislet = new("ExpressionSet", exprs = 
               t(as.matrix(count_data))
)
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

#meta_data = data.frame(
#  "Sample" = colnames(exprs(eset)),
#  "Deco_type" = cell_type
#)
#rownames(meta_data) = meta_data$Sample

#write.table(meta_info, "~/Deko/Misc/Meta_information.tsv",sep="\t", row.names = F, quote = F)

pheatmap::pheatmap(
  t(res_coeff),
  annotation_col = meta_data[c("Deco_type")],
  annotation_colors = aka3,
  annotation_legend = T,
  treeheight_col = 0,
  show_colnames = T,
  show_rownames = T,
  gaps_row = 20
)

#write.table(meta_info, "~/Deko/Misc/Deko_table.tsv",sep ="\t", row.names = F, quote = F)