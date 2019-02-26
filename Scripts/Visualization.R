library("stringr")
library("ggplot2")

source("~/Deko/Scripts/Visualization_colors.R")

genes_of_interest_hgnc_t = read.table("~/Deko/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1

sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[13,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
expr_raw = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv",sep="\t", stringsAsFactors =  F, header = T)
#expr_raw = read.table("~/Deko/Data/Visualization_PANnen.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw))
colnames(expr) = colnames(expr_raw)
rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
cor_mat = cor(expr[,]);pcr = prcomp(t(cor_mat))

#meta_data = meta_info[as.character(rownames(cor_mat)),]

vis_mat = meta_data[,c(
    "Differentiated_similarity",
    "Progenitor_similarity",
    "Stem_cell_similarity"
)]
vis_mat = meta_data[,c(
    "alpha",
    "beta",
    "gamma",
    "delta",
    "acinar",
    "ductal"
)]
vis_mat = meta_data[,c(
    "alpha_similarity_absolute",
    "beta_similarity_absolute",
    "gamma_similarity_absolute",
    "delta_similarity_absolute",
    "progenitor_similarity_absolute",
    "stem_cell_similarity_absolute"
)]
colnames(vis_mat) = str_replace(colnames(vis_mat),pattern = "(_absolute)|(_relative)","")

pheatmap::pheatmap(
    cor_mat,
    annotation_col = vis_mat[,c("alpha","beta","gamma","delta","acinar","ductal")],
    annotation_colors = aka3,
    annotation_legend = T,
    treeheight_col = 0,
    treeheight_row = 0,
    show_colnames = F,
    show_rownames = F#,
    #color = colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(length(breaksList)),
    #cluster_cols = F, cluster_rows = F
)

groups = as.character(unlist(meta_data["Differentiation_type"]))
fill_vec = groups
fill_vec[fill_vec == "alpha"] = "Blue"
fill_vec[fill_vec == "beta"] = "yellow"
fill_vec[fill_vec == "gamma"] = "orange"
fill_vec[fill_vec == "delta"] = "purple"
fill_vec[fill_vec == "ductal"] = "cyan"
fill_vec[fill_vec == "acinar"] = "brown"
fill_vec[fill_vec ==  "not_sig"] = "gray"

#size_vec = as.double(log(rowSums(res_coeff[ rownames(pcr$x),c("alpha","beta","gamma","delta")])+1))^2
fill_vec[fill_vec %in% "Differentiated"] = "Green"
fill_vec[fill_vec %in% "Not_significant"] = "gray"
fill_vec[fill_vec %in% "Hisc"] = "Black"
fill_vec[fill_vec %in% "Progenitor"] = "red"

size_vec = meta_data$Grading
size_vec[size_vec == "G3"] = 6
size_vec[size_vec == "G2"] = 4
size_vec[size_vec == "G1"] = 1
size_vec = as.integer(size_vec)

meta_data = meta_data[rownames(pcr$x),]
p = ggbiplot::ggbiplot(
    pcr,
    obs.scale = .75,
    groups = meta_data[,"Grading"],
    ellipse = TRUE,
    circle = TRUE,
    var.axes = F#,labels = meta_data$Name
)
p = p + geom_point( aes( colour = fill_vec, shape = meta_data[,"Grading"]), size = size_vec)
p = p + scale_color_manual( values = c("Black","Purple","Orange","Black","gray","brown") )
#p = p + guides( color=guide_legend(title="Study", size=guide_legend(title="MKI67"), shape = guide_legend(title="Grading")))
p

# UMAP

vis_mat = bam_data
vis_mat$Sample = vis_mat$Name
#vis_mat  = vis_mat[,c("Sample","HESC_sim","Progenitor_sim","Differentiated_sim","NEC_NET")]
vis_mat  = vis_mat[,c("Sample","expr_stem","NEC_NET")]
vis_mat = reshape2::melt(vis_mat )
colnames(vis_mat) = c("Sample","NEC_NET","exp_type","expr_stem")
vis_mat$Sample = factor( vis_mat$Sample, levels = vis_mat$Sample[order(vis_mat$expr_stem, decreasing = F)] )

col_vec = as.character(vis_mat$NEC_NET[order(vis_mat$expr_stem, decreasing = F)])
col_vec[col_vec == "NEC"] = "red"
col_vec[col_vec != "red"] = "gray"

vis_mat

#cor.test(meta_data$HSC_sim,meta_data$Differentiated_sim)

men1_plot = ggplot( data = vis_mat, aes ( x = Sample,  y = vis_mat$expr_stem))
men1_plot = men1_plot + geom_bar(stat="identity", position=position_dodge(), fill = col_vec)
men1_plot = men1_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#men1_plt = men1_plot  +  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
men1_plot + coord_cartesian(ylim = c(9,10.5)) + theme_bw()

meta_info$Differentiated_sim = rep("",nrow(meta_info))
meta_info[rownames(meta_data),"Differentiated_sim"] = meta_data$Differentiated_sim
meta_info$Progenitor_sim = rep("",nrow(meta_info))
meta_info[rownames(meta_data),"Progenitor_sim"] = meta_data$Progenitor_sim
meta_info$HESC_sim = rep("",nrow(meta_info))
meta_info[rownames(meta_data),"HESC_sim"] = meta_data$HESC_sim

#write.table(meta_info, "~/Deko/Misc/Meta_information.tsv", sep ="\t", quote = F, row.names = F)

###
library("stringr")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))


### Prep

library("grid")

#meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_data = meta_info[ match(colnames(expr_raw), meta_info$Name), ]
rownames(meta_data) = meta_data$Name
df_map = match(colnames(expr_raw), meta_data$Name)
meta_data$INS = as.double(expr_raw[ which( rownames(expr_raw) == "INS"),df_map])
meta_data$GCG = as.double(expr_raw[ which( rownames(expr_raw) == "GCG"),df_map])
meta_data$PPY = as.double(expr_raw[ which( rownames(expr_raw) == "PPY"),df_map])
meta_data$SST = as.double(expr_raw[ which( rownames(expr_raw) == "SST"),df_map])
meta_data$MEN1_exp = as.double(expr_raw[ which( rownames(expr_raw) == "MEN1"),df_map])
meta_data$Location = str_replace_all(meta_data$Location, pattern = "-", "_")
meta_data$MKI67 = as.double(expr_raw[ which( rownames(expr_raw) == "MKI67"),])
meta_data$MKI67 = meta_data$MKI67 + ( sign(min(meta_data$MKI67)) *  min(meta_data$MKI67) )
meta_data$MKI67 = meta_data$MKI67/ max(meta_data$MKI67)
meta_data$Marker_Genes = meta_data$INS + meta_data$GCG + meta_data$PPY + meta_data$SST
meta_data$Marker_Genes = meta_data$Marker_Genes + ( sign(min(meta_data$Marker_Genes)) *  min(meta_data$Marker_Genes) )
meta_data$Marker_Genes = meta_data$Marker_Genes/ max(meta_data$Marker_Genes)

# linear correlation MEN1

vis_mat_2 = subset(vis_mat, MEN1_mt_AF > 0)
lm.model <- lm( vis_mat_2$MEN1_exp ~ vis_mat_2$MEN1_mt_AF) # Fit linear model
summary(lm.model)
cor(vis_mat_2$MEN1_exp , vis_mat_2$MEN1_mt_AF)

# Extract fitted coefficients from model object
b0 <- lm.model$coefficients[1]
b1 <- lm.model$coefficients[2]

g_bench = ggplot(
  data = vis_mat_2,
  aes( y =  MEN1_mt_AF,x = MEN1_exp)
)
g_bench = g_bench + geom_point( )
g_bench = g_bench + geom_smooth(method = "lm")

g_bench + xlab("MEN1 expression")+ ylab("MEN1 mutation AF")

### survival plots ###

#meta_data$Diff_Type[meta_data$Diff_Type %in% c("Alpha","Beta","Gamma","Delta")]  = "Differentiated"

meta_data = meta_info[rownames(vis_mat),]
vis_mat$OS_Tissue = as.double(str_replace_all(meta_data$OS_Tissue, pattern = ",", "\\."))
vis_mat$OS_Tissue[is.na(vis_mat$OS_Tissue)] = 1
vis_mat$Grading = meta_data$Grading
vis_mat$Zensur = meta_data$Zensur

data = vis_mat[,c("hisc","OS_Tissue","Zensur")]
data = data[ !is.na(data$OS_Tissue),]
colnames(data) = c("hisc_similarity","OS_Tissue","Status")
data = data[data$hisc_similarity!="not_significant",]

fit = survival::survfit( survival::Surv( as.double(data$OS_Tissue), data$Status ) ~ data$hisc_similarity)

survminer::ggsurvplot(fit, data = data, risk.table = F, pval = T, censor.size = 10)

###

#### DeconRNASeq

res = DeconRNASeq(
    bam_data,
    signatures,
    proportions = NULL,
    checksig = FALSE,
    known.prop = FALSE,
    use.scale = TRUE,
    fig = FALSE
)

deco_res = res$out.all
rownames(deco_res) = colnames(query_data)

deco_res = as.data.frame(deco_res)
deco_res = apply(deco_res, MARGIN = 1, FUN = function(vec){return(log(vec+1))})
annotation_data = t(deco_res)
annotation_data[1:5,]
old_colnames = colnames(annotation_data)

### ROCR

library("ggplot2")

res_hisc_1
res_hisc_2

res_both_1
res_both_2

res_hisc_2$Var.prop
res_both_2$Var.prop

# deconvolution_results_baron_wiedemann = deconvolution_results
# deconvolution_results_GSE73338_full = deconvolution_results

vis_mat
data_mat = reshape2::melt(deconvolution_results )

data_mat = deconvolution_results[,c("Grading","MKI67")]
data_mat$Sample = rownames(data_mat)
data_mat$Sample = factor(data_mat$Sample, levels = data_mat$Sample[order(data_mat$MKI67)] )
data_mat$MKI67 = data_mat$MKI67 + 1

color_vec = data_mat$Grading
color_vec[color_vec == "G1"] = "green"
color_vec[color_vec == "G2"] = "yellow"
color_vec[color_vec == "G3"] = "red"
color_vec = color_vec[order(data_mat$MKI67)]

men1_plot = ggplot( data_mat, aes ( x = Sample,  y = MKI67))
men1_plot = men1_plot + geom_bar( aes(fill = Grading), stat = "identity")
men1_plot = men1_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
men1_plot = men1_plot + scale_fill_manual( values = c("darkgreen","yellow","red"))
men1_plot

ggplot(data_mat,aes( x = Grading, y = MKI67, fill = Grading )) + geom_boxplot( )

### ratio

deconvolution_results$Ratio_numeric = as.double(vis_mat[rownames(deconvolution_results),"Ratio_numeric"])
data_mat = deconvolution_results[,c("Grading","Ratio_numeric")]
colnames(data_mat) = c("Grading","Gene_Expression")
data_mat$Sample = rownames(data_mat)
data_mat$Sample = factor(data_mat$Sample, levels = data_mat$Sample[order(data_mat$Gene_Expression)] )
#data_mat$Gene_Expression = data_mat$Gene_Expression + 1

color_vec = data_mat$Grading
color_vec[color_vec == "G1"] = "green"
color_vec[color_vec == "G2"] = "yellow"
color_vec[color_vec == "G3"] = "red"
color_vec = color_vec[order(data_mat$Gene_Expression)]

men1_plot = ggplot( data_mat_2, aes ( x = Sample,  y = Gene_Expression))
men1_plot = men1_plot + geom_bar( aes(fill = Grading), stat = "identity")
men1_plot = men1_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
men1_plot + scale_fill_manual( values = c("darkgreen","red"))

bin_width = 10
data_mat_2 = data_mat
data_mat_2$Grading[data_mat_2$Grading %in% c("G1","G2")] = "G1_G2"
ggplot(data_mat_2,aes( x = Grading, y = Gene_Expression, fill = Grading )) + geom_boxplot( )

geom_histogram(data=subset(data_mat,Grading == 'G3'),fill = "red", alpha = 0.5,position="dodge", bins = bin_width) +
    geom_histogram(data=subset(data_mat,Grading == 'G2'),fill = "yellow", alpha = 0.5,position="dodge", bins = bin_width) +
    geom_histogram(data=subset(data_mat,Grading == 'G1'),fill = "green", alpha = 0.5,position="dodge", bins = bin_width)


### ROC curve

#deconvolution_results_wiedenmann_scarpa = deconvolution_results
#deconvolution_results_GSE73338 = deconvolution_results

library("ROCR")
library("ggplot2")

target_labels = deconvolution_results$Grading
target_labels[target_labels %in% c("G1","G2")] = 0
target_labels[target_labels %in% c("G3")] = 1

prediction_vector_ratio = as.double(vis_mat[rownames(deconvolution_results),"Ratio_numeric"])
#prediction_vector_ratio = prediction_vector_ratio - min(prediction_vector_ratio)
#prediction_vector_ratio = prediction_vector_ratio / max(prediction_vector_ratio)
pred_ratio = prediction(prediction_vector_ratio, target_labels)
perf_ratio = performance(pred_ratio, measure = "tpr", x.measure = "fpr")
auc_ratio = performance(pred_ratio, measure = "auc", x.measure = "fpr") 

prediction_vector_mki67 = as.double(deconvolution_results[,"MKI67"])
prediction_vector_mki67 = prediction_vector_mki67 - min(prediction_vector_mki67)
prediction_vector_mki67 = prediction_vector_mki67 / max(prediction_vector_mki67)
pred_mki67 <- prediction(prediction_vector_mki67, target_labels)
perf_mki67 = performance(pred_mki67, measure = "tpr", x.measure = "fpr") 

tpr_mki67 = as.double(as.character(unlist(perf_mki67@y.values)))
tpr_mki67 = round(tpr_mki67 * 100,0)
tpr_ratio = as.double(as.character(unlist(perf_ratio@y.values)))
tpr_ratio = round(tpr_ratio * 100,0)

fpr_mki67 = as.double(as.character(unlist(perf_mki67@x.values)))
fpr_mki67 = round(fpr_mki67 * 100,0)
fpr_ratio = as.double(as.character(unlist(perf_ratio@x.values)))
fpr_ratio = round(fpr_ratio * 100,0)

#perf <- performance(pred, measure = "tpr", x.measure = "fpr") 

method_indicator = c(
    rep("mki67",length(tpr_mki67)),
    rep("ratio",length(tpr_ratio))
)

roc_mat = data.frame(
    "tpr" = c(tpr_mki67,tpr_ratio),
    "fpr" = c(fpr_mki67,fpr_ratio),
    "method" = method_indicator
)

roc_plot = ggplot(
    data = roc_mat,
    aes( x = fpr, y = tpr)
)
roc_plot = roc_plot + geom_line(aes(color = method))
roc_plot
