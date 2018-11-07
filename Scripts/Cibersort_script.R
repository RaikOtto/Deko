names(pancreasMarkers) = str_to_lower(names(pancreasMarkers))
eislet = new("ExpressionSet", exprs = as.matrix(count_data))

sub_list = str_to_lower(subtypes)
names(sub_list) = names(subtypes)
fData(eislet) = data.frame( sub_list  )
pData(eislet) = data.frame( sub_list )

B = bseqsc_basis(
    eislet,
    pancreasMarkers,
    clusters = 'sub_list',
    samples = colnames(exprs(eislet)),
    ct.scale = FALSE
)
plotBasis(B, pancreasMarkers, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')

### RUN VARIANCE SELECTION FIRST

eset = new("ExpressionSet", exprs=as.matrix(bam_data));

fit = bseqsc_proportions(eset, B, verbose = TRUE, absolute = T, log = F, perm = 100)

source("~/Deko/Scripts/Utility_script.R")
#meta_data$Diff_Type =  c("Differentiated","Progenitor","HSC")[maxi]

#expr_raw = bam_data
#colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

genes_of_interest_hgnc_t = read.table("~/Deko/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[13,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
#expr = matrix(as.double(as.character(unlist(expr_raw))), ncol = ncol(expr_raw))
expr_raw = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw))
colnames(expr) = colnames(expr_raw)
rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

meta_data = meta_data[as.character(rownames(cor_mat)),]
pheatmap::pheatmap(
    #t(res_coeff),
    cor(expr),
    annotation_col = meta_data[c("HSC_sim","Progenitor_sim","Differentiated_sim","NEC_NET")],
    #annotation_col = meta_data[c("Alpha_sim","Beta_sim","Gamma_sim","Delta_sim","Ductal_sim","Acinar_sim","NEC_NET")],
    #annotation_col = meta_data[c("Alpha_sim","Beta_sim","Gamma_sim","Delta_sim","Ductal_sim","Acinar_sim","NEC_NET")],
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

size_vec = as.double(log(rowSums(res_coeff[ rownames(pcr$x),c("alpha","beta","gamma","delta")])+1))^2
#size_vec = max(size_vec) - size_vec

meta_data = meta_data[rownames(pcr$x),]
p = ggbiplot::ggbiplot(
  pcr,
  obs.scale = .75,
  groups = meta_data$Diff_Type,
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

#HSC_Dif = subset(vis_mat, Similarity_type %in% c("Differentiated_sim","HSC_sim"))
plot(meta_data$Differentiated_sim, meta_data$HSC_sim)
lm_calc = lm(meta_data$Differentiated_sim ~ meta_data$HSC_sim)
summary(lm_calc)

cor_plot = ggplot( data = meta_data, aes(x = Differentiated_sim, y = HSC_sim))
cor_plot = cor_plot + geom_line() +geom_point( ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
cor_plot + geom_smooth(method='lm')
hist(meta_data$Differentiated_sim)

table(meta_data$Diff_Type)
