names(pancreasMarkers) = str_to_lower(names(pancreasMarkers))
eislet = new("ExpressionSet", exprs = as.matrix(count_data))

sub_list = str_to_lower(subtypes)
names(sub_list) = names(subtypes)
fData(eislet) = data.frame( sub_list  )
pData(eislet) = data.frame( sub_list )

B = readRDS("~/ArtDeco/inst/Models/Four_differentiation_stages.RDS")
"
B = bseqsc_basis(
eislet,
pancreasMarkers,
clusters = 'sub_list',
samples = colnames(exprs(eislet)),
ct.scale = FALSE
)
plotBasis(B, pancreasMarkers, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')
"

### RUN VARIANCE SELECTION FIRST

bam_data = read.table("~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",sep ="\t", header = T, stringsAsFactors = F)
colnames(bam_data) = str_replace_all(colnames(bam_data), "\\.", "_")

subtypes = meta_info[colnames(bam_data),"Subtype"]
cands = which(!is.na(subtypes)) & (subtypes %in% c("Alpha","Beta","Gamma","Delta"))

expr_raw = bam_data[,cands]
bam_data = bam_data[,cands]
eset = new("ExpressionSet", exprs=as.matrix(bam_data));

#fit = bseqsc_proportions(eset, B, verbose = TRUE, absolute = T, log = F, perm = 100)

fit = readRDS("~/Deko/Misc/fit_baron.RDS")
fit = readRDS("~/Deko/Misc/Fit_Muraro.rds")
fit = readRDS("~/Deko/Misc/Fit_Segerstolpe.rds")

meta_data = meta_info[names(fit$stats[,1]),]
#source("~/Deko/Scripts/Utility_script.R")
meta_data$Y_Subtype = as.character(meta_data[names(fit$stats[,1]),"Subtype"])
#saveRDS(fit, "~/Deko/Misc/Fit_Muraro.rds")

table(meta_data[,c("Y_Subtype","Diff_Type")])

### p_value

stats_t = as.data.frame(res_cor)

#stats_t_glob <<- stats_t
stats_t_glob <<- rbind(stats_t_glob,stats_t)
dim(stats_t_glob)

###

library("ROCR")
library("caret")
#fit = readRDS("~/Downloads/fit.RDS")
pred_vec = meta_data$Y_Subtype == meta_data$Diff_Type
pred_vec[pred_vec == TRUE] = 1
pred_vec[pred_vec != 1] = 0
pred_vec = as.double(pred_vec)

#pred <- prediction(predictions = pred_vec, labels =  label_vec )
pred <- prediction(predictions = 1-as.double(stats_t$`P-value`), labels =  pred_vec )

roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc.perf)
abline(a=0, b= 1)

opt.cut = function(roc.perf, pred){
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
          cutoff = p[[ind]])
    }, roc.perf@x.values, roc.perf@y.values, pred@cutoffs)
}
print(opt.cut(roc.perf, pred))

# View confusion matrix overall
result <- confusionMatrix( as.factor(meta_data$Diff_Type), as.factor(meta_data$Subtype))
result 

# F1 value
result$byClass[7] 

f1_score <- function(predicted, expected, positive.class="1") {
    predicted <- factor(as.character(predicted), levels=unique(as.character(expected)))
    expected  <- as.factor(expected)
    cm = as.matrix(table(expected, predicted))
    
    precision <- diag(cm) / colSums(cm)
    recall <- diag(cm) / rowSums(cm)
    f1 <-  ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
    
    #Assuming that F1 is zero when it's not possible compute it
    f1[is.na(f1)] <- 0
    
    #Binary F1 or Multi-class macro-averaged F1
    ifelse(nlevels(expected) == 2, f1[positive.class], mean(f1))
}

###
library(DeconRNASeq)
library(stringr)

bam_data = read.table("~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",sep ="\t",stringsAsFactors = F,header = T,row.names= 1)
colnames(bam_data) =
    str_replace(colnames(bam_data), pattern = "^X", "")

basis_t = read.table("~/Deko/Data/Progenitor_Stanescu_HISC_Haber.tsv",sep ="\t",stringsAsFactors = F)

colnames(basis_t) = str_replace_all(colnames(basis_t),pattern = "\\.","_")
meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

meta_data = meta_info[colnames(basis_t),]
subtypes = meta_data$Subtype
table(subtypes)

dim(basis_t)
#basis_t = basis_t[,which(subtypes %in% c("Alpha","Beta","Gamma","Delta"))]
dim(basis_t)

basis_t = basis_t[,!is.na(meta_data$Subtype)]
subtypes = meta_data$Subtype[!is.na(meta_data$Subtype)]
table(subtypes)
length(subtypes)
dim(basis_t)

signature_genes = c()
for ( subtype in unique(subtypes)){
    signature_genes = c(signature_genes,
        identify_marker_genes(
            expression_training_mat = basis_t,
            subtype_vector = subtypes,
            subtype = subtype,
            nr_marker_genes = 100
        )
    )
}

basis_t = basis_t = basis_t[signature_genes,]

signatures = apply(basis_t,MARGIN=1,FUN = function(vec){
    measurements = aggregate(
        as.double(vec),
        by = list(subtypes),
        FUN = mean
    )
    return(measurements[,2])
    }
)
signatures = as.data.frame(t(signatures))
colnames(signatures) = levels(factor(subtypes))
signatures[1:5,]

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

max_val = apply(annotation_data, MARGIN = 1, FUN = which.max)

colnames(annotation_data) = c("Alpha_similarity","Beta_similarity","Delta_similarity","Gamma_similarity")
colnames(annotation_data) = c("Progenitor_similarity","HISC_similarity")
#colnames(annotation_data) = c("alpha_similarity","beta_similarity","delta_similarity","gamma_similarity","stem_cell_similaritry","progenitor_simimilarity")

annotation_data = as.data.frame(annotation_data)

Differentiatedness = apply(annotation_data[,c("Alpha_similarity","Beta_similarity","Delta_similarity","Gamma_similarity")],MARGIN = 1, FUN = sum)
annotation_data$Differentiatedness = Differentiatedness

De_differentiatedness = apply(annotation_data[,c("Progenitor_similarity","HISC_similarity")],MARGIN = 1, FUN = sum)
annotation_data$De_differentiatedness = De_differentiatedness

annotation_data[,"Max_sim"] = rep("",nrow(annotation_data))
annotation_data[,"Max_sim"] = old_colnames[max_val]

###

Graphics_parameters = c("")

create_heatmap_differentiation_stages(
    "~/Deko/Data/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv",
    annotation_data
)

pheatmap::pheatmap(
    correlation_matrix,
    annotation_col = annotation_data,
    annotation_colors = Graphics_parameters,
    annotation_legend = TRUE,
    treeheight_col = 0,
    treeheight_row = 0,
    show_colnames = TRUE,
    show_rownames = FALSE
)

###

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

target_labels = deconvolution_results$Grading
target_labels[target_labels %in% c("G1","G2")] = 0
target_labels[target_labels %in% c("G3")] = 1
prediction_vector = deconvolution_results$Ratio_numeric
prediction_vector = prediction_vector - min(prediction_vector)
prediction_vector = prediction_vector / max(prediction_vector)

pred <- prediction(prediction_vector, target_labels)
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10))
