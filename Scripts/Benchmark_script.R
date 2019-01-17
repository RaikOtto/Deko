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

basis_t = read.table("~/Deko/Data/Merged_Dif_Baron_Prog_Stanescu_Hisc_Haber.tsv",sep ="\t",stringsAsFactors = F)
basis_t = read.table("~/Deko/Data/Merge_mat_HSC_Stanescu_Baron.tsv",sep ="\t",stringsAsFactors = F)
colnames(basis_t) = str_replace_all(colnames(basis_t),pattern = "\\.","_")
meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

meta_data = meta_info[colnames(bam_data),]
subtypes = meta_data$Subtype
table(subtypes)

dim(basis_t)
#basis_t = basis_t[,which(subtypes %in% c("Alpha","Beta","Gamma","Delta"))]
dim(basis_t)

rownames(meta_info) = meta_info$Name
meta_data = meta_info[colnames(basis_t),]
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

query_data = read.table("~/Deko/Data/PANnen_Test_Data.tsv",sep ="\t",stringsAsFactors = F,header = T,row.names= 1)
query_data[1:5,1:5]

res = DeconRNASeq(
    query_data,
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

#colnames(annotation_data) = c("alpha_similarity","beta_similarity","delta_similarity","gamma_similarity")
colnames(annotation_data) = c("alpha_similarity","beta_similarity","delta_similarity","gamma_similarity","stem_cell_similaritry","progenitor_simimilarity")

annotation_data = as.data.frame(annotation_data)
Differentiatedness = apply(annotation_data[,c("alpha_similarity","beta_similarity","delta_similarity","gamma_similarity")],MARGIN = 1, FUN = sum)
annotation_data$Differentiatedness = Differentiatedness

annotation_data_scaled = t(apply(annotation_data, FUN = scale, MARGIN = 1))
colnames(annotation_data_scaled) = colnames(annotation_data)

annotation_data[,"Max_sim"] = rep("",nrow(annotation_data))
annotation_data[,"Max_sim"] = old_colnames[max_val]

###

transcriptome_file_path = system.file(
    "/Data/Expression_data/Visualization_PANnen.tsv",
    package="artdeco")
baseline = "relative"
deconvolution_results = deconvolution_results_relative
Graphics_parameters = c("")
meta_data = deconvolution_results
    
create_heatmap_differentiation_stages(
    visualization_data_path,
    deconvolution_results_relative
)

annotation_columns = c(
    "Differentiation_Stages_Subtypes",
    "Differentiation_Stages_Aggregated",
    "Differentiatedness"
)

create_heatmap_differentiation_stages(
    visualization_data_path,
    deconvolution_results_absolute
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