
eislet = new("ExpressionSet", exprs = 
               (as.matrix(count_data))
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

eset = new("ExpressionSet", exprs=as.matrix(expr_raw));
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
meta_data  = meta_info[meta_match,]

meta_data$Deco_cor = rep("", length(meta_match))
meta_data$Deco_cor  = as.character(res_cor[,2])

res_coeff_match = match( rownames(res_coeff), meta_info$Name, nomatch = 0 )
res_cor_match = match( rownames(res_cor), meta_info$Name, nomatch = 0 )

maxi = apply( res_coeff, FUN = which.max, MARGIN = 1 )
cell_type = colnames(res_coeff)[maxi]
rownames(meta_info) = meta_info$Name

cell_type[res_cor[,"Correlation"] <= .3] = "Not_sig"

#meta_data = data.frame(
#  "Sample" = colnames(exprs(eset)),
#  "Deco_type" = cell_type
#)
#rownames(meta_data) = meta_data$Sample

#write.table(meta_info, "~/Deko/Misc/Meta_information.tsv",sep="\t", row.names = F, quote = F)

meta_data$NEC_NET[meta_data$NEC_NET == ""] = "Unknown"
meta_data$Subtype_Sadanandam[meta_data$Grading == ""] = "Norm"
meta_data$Grading[meta_data$Grading == ""] = "G0"
meta_data$NEC_NET[meta_data$NEC_NET == ""] = "Unknown"
meta_data$Subtype_Sadanandam[meta_data$Subtype_Sadanandam == ""] = "Unknown"
meta_data$Deco_type = cell_type

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
