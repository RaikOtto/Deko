library("stringr")
library("ggplot2")

source("~/Deko/Scripts/Visualization_colors.R")

###
library("stringr")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

### lm plot

lm.model <- lm(deconvolution_results$MKI67 ~ deconvolution_results$ductal) # Fit linear model
summary(lm.model)
correlation = round(cor(deconvolution_results$MKI67, deconvolution_results$ductal),2)
cor.test(deconvolution_results$MKI67, deconvolution_results$ductal)

g_bench = ggplot( aes(y =  MKI67,x = ductal ), data = deconvolution_results)
g_bench = g_bench + geom_point( aes( size = 4)) + scale_size(guide="none")
g_bench = g_bench + geom_smooth(method = "lm")
g_bench = g_bench +  annotate( "text", x = 2, y = 2, label = as.character(correlation), size =10)
plot(g_bench)

### survival plots ###

ratio = vis_mat[rownames(deconvolution_results),"Ratio_numeric"]
ratio[ratio <= mean(ratio)] = "low"
ratio[ratio != "low"] = "high"

data = vis_mat[,c("Ratio","OS_Tissue","Zensur")]
value = as.double(as.character(unlist(data[,1])))
fit = survival::survfit( survival::Surv( as.double(data$OS_Tissue), data$Zensur ) ~ value)
survminer::ggsurvplot(fit, data = data, risk.table = F, pval = T, censor.size = 10)

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
data_mat = deconvolution_results[,c("Grading","MKI67")]
colnames(data_mat) = c("Grading","Gene_Expression")
data_mat$Sample = rownames(data_mat)
data_mat$Sample = factor(data_mat$Sample, levels = data_mat$Sample[order(data_mat$Gene_Expression)] )
#data_mat$Gene_Expression = data_mat$Gene_Expression + 1

color_vec = data_mat$Grading
color_vec[color_vec == "G1"] = "green"
color_vec[color_vec == "G2"] = "yellow"
color_vec[color_vec == "G3"] = "red"
color_vec = color_vec[order(data_mat$Gene_Expression)]

#men1_plot = ggplot( data_mat_2, aes ( x = Sample,  y = Gene_Expression))
men1_plot = men1_plot + geom_bar( aes(fill = Grading), stat = "identity")
men1_plot = men1_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
men1_plot + scale_fill_manual( values = c("darkgreen","red"))

bin_width = 10
#data_mat_2 = data_mat
#data_mat_2$Grading[data_mat_2$Grading %in% c("G1","G2")] = "G1_G2"
ggplot(data_mat,aes( x = Grading, y = Gene_Expression, fill = Grading )) + geom_boxplot( )

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
