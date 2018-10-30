names(pancreasMarkers) = str_to_lower(names(pancreasMarkers))
eislet = new("ExpressionSet", exprs = (as.matrix(count_data)))
subtypes = str_to_lower(subtypes)
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

fit = bseqsc_proportions(eset, B, verbose = TRUE, log = F, absolute = T, perm = 200)

###

res_coeff = t(fit$coefficients)
res_cor   = fit$stats

res_coeff[ is.na(res_coeff) ] = 0.0
res_cor[ is.na(res_cor) ] = 0.0

not_sig_samples = rownames(res_cor)[res_cor[,"P-value"] > .05]
res_cor[ res_cor[,"Correlation"] <= .3, "Correlation" ] = 0.0

meta_data = meta_info[ rownames(res_coeff),]

meta_data$Progenitor_sim = log(res_coeff[,"e13.5"]+1)
meta_data$HSC_sim = log(res_coeff[,"hsc"]+1)
meta_data$Differentiated_sim = log(rowSums(res_coeff[,c("alpha","beta","gamma","delta")])+1)

###

res_coeff_match = match( rownames(res_coeff), meta_data$Name, nomatch = 0 )
res_cor_match = match( rownames(res_cor), meta_info$Name, nomatch = 0 )

maxi = apply( res_coeff, FUN = which.max, MARGIN = 1 )
cell_type = colnames(res_coeff)[maxi]
#cell_type[res_cor[,"Correlation"] <= .3] = "Not_sig"
#cell_type[res_cor[,"P-value"] > .05] = "Not_sig"
table(cell_type)

meta_info$Deco_type = rep("",nrow(meta_info))
meta_info$Deco_type[ match(rownames(res_coeff), meta_info$Name) ] = as.character( cell_type )
meta_data = meta_info[ rownames(res_coeff),]
meta_data$Grading[ meta_data$Grading == ""] = "G0"

meta_data$Deco_similarity = rep(0.0,nrow( res_coeff ) )
meta_data$Deco_similarity = as.double(colSums(t(res_coeff)))
meta_data$NEC_NET[meta_data$NEC_NET == ""] = "NA"
colmax = apply(t(res_coeff),  function(vec){max(as.double(vec))}, MARGIN = 2 )
vis_mat =  apply( t(res_coeff), MARGIN = 2, FUN = function(vec){
    if (max(vec) != 0 ){ return( vec / max(vec)) } else { return ( rep(0, length(vec) ) ) }})
vis_mat = t(res_coeff)
vis_mat[vis_mat > 5] = 5
rownames(vis_mat) = colnames(res_coeff)

dd = meta_data$Deco_similarity
meta_data$Deco_similarity[meta_data$Deco_similarity > 5] = 5
meta_data$Deco_type[meta_data$Deco_similarity < 1]  = "not_sig"
meta_data$NEC_NET[meta_data$NEC_NET == "Unknown"] = "NA"
#meta_data$Progenitor_sim[meta_data$Progenitor_sim <= 3 ] = 3
meta_data$HSC_sim[meta_data$HSC_sim <= 3 ] = 3

pheatmap::pheatmap(
    #t(meta_data[order(meta_data$HSC_sim),c("HSC_sim","Progenitor_sim","Differentiated_sim")]),
    cor(expr),
    annotation_col = meta_data[c("HSC_sim","Progenitor_sim","Differentiated_sim","NEC_NET","Grading")],
    annotation_colors = aka3,
    annotation_legend = T,
    treeheight_col = 0,
    treeheight_row = 0,
    show_colnames = T,
    show_rownames = F#,
    #color = colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(length(breaksList)),
    #cluster_cols = F, cluster_rows = F
)

#groups = as.character(unlist(meta_data["Deco_type"]))
groups = as.character(apply(vis_mat, MARGIN = 2, FUN = function(vec){return(rownames(vis_mat)[which.max(vec)])}))
groups[(log(colSums(t(res_coeff)))) < 0] = "Not_sig"
meta_data$Groups = groups
fill_vec = meta_data$Deco_type
fill_vec[fill_vec == "alpha"] = "Blue"
fill_vec[fill_vec == "beta"] = "yellow"
fill_vec[fill_vec == "gamma"] = "orange"
fill_vec[fill_vec == "delta"] = "purple"
fill_vec[fill_vec == "ductal"] = "cyan"
fill_vec[fill_vec == "acinar"] = "brown"
fill_vec[fill_vec ==  "not_sig"] = "gray"
size_vec = as.double(meta_data$Deco_similarity**1)
size_vec = max(size_vec) - size_vec

meta_data = meta_data[rownames(pcr$x),]
p = ggbiplot::ggbiplot(
  pcr,
  obs.scale = .75,
  groups = meta_data$Deco_type,
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F#,labels = meta_data$Name
)
p = p + geom_point( aes( size = size_vec, colour = meta_data$Deco_type, shape = meta_data$Grading))
p = p + scale_color_manual( values = c("Blue","Yellow","Purple","Orange","Gray") )
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

meta_data$Progenitor_sim = scale(meta_data$Progenitor_sim)
meta_data$HSC_sim = scale(meta_data$HSC_sim)
meta_data$Differentiated_sim = scale(meta_data$Differentiated_sim)

vis_mat = meta_data
vis_mat$Sample = vis_mat$Name
vis_mat  = vis_mat[,c("Sample","HSC_sim","Progenitor_sim","Differentiated_sim","NEC_NET")]
vis_mat = reshape2::melt(vis_mat )
colnames(vis_mat) = c("Sample","NEC_NET","Similarity_type","Similarity")
vis_mat$Sample = factor( vis_mat$Sample, levels = meta_data$Name[order(meta_data$HSC_sim, decreasing = F)] )

col_vec = as.character(vis_mat$Similarity_type)
col_vec[col_vec == "Differentiated_sim"] = "Green"
col_vec[col_vec == "Progenitor_sim"] = "Yellow"
col_vec[col_vec == "HSC_sim"] = "red"

vis_mat

#cor.test(meta_data$HSC_sim,meta_data$Differentiated_sim)

men1_plot = ggplot( data = vis_mat, aes ( x = Sample,  y = Similarity, group = Similarity_type, colour = Similarity_type))
men1_plot = men1_plot + geom_line() +geom_point( ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#men1_plot = men1_plot + geom_bar(stat="identity", position=position_dodge())
men1_plot = men1_plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
men1_plot

#

HSC_Dif = subset(vis_mat, Similarity_type %in% c("Differentiated_sim","HSC_sim"))
plot(meta_data$Differentiated_sim, meta_data$HSC_sim)
lm_calc = lm(meta_data$Differentiated_sim ~ meta_data$HSC_sim)
summary(lm_calc)

cor_plot = ggplot( data = meta_data, aes(x = Differentiated_sim, y = HSC_sim))
cor_plot = cor_plot + geom_line() +geom_point( ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
cor_plot + geom_smooth(method='lm')
hist(meta_data$Differentiated_sim)
