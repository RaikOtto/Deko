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

bam_data = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Lawlor.tsv",sep ="\t", header = T, stringsAsFactors = F)
colnames(bam_data) = str_replace_all(colnames(bam_data), "\\.", "_")

subtypes = meta_info[colnames(bam_data),"Subtype"]
cands = which(!is.na(subtypes)) & (subtypes %in% c("Alpha","Beta","Gamma","Delta"))
meta_data = meta_info[colnames(bam_data)[cands],]

expr_raw = bam_data[,cands]
bam_data = bam_data[,cands]
eset = new("ExpressionSet", exprs=as.matrix(bam_data));

fit = bseqsc_proportions(eset, B, verbose = TRUE, absolute = T, log = F, perm = 100)
source("~/Deko/Scripts/Utility_script.R")
meta_data$Y_Subtype = as.character(subtypes[cands])
#saveRDS(fit, "~/Deko/Misc/Fit_Lawlor.rds")

#meta_data$Diff_Type =  c("Differentiated","Progenitor","HSC")[maxi]

#expr_raw = bam_data
#colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

pheatmap::pheatmap(
    t(res_coeff),
    #cor(expr),
    annotation_col = meta_data[c("Y_Subtype","Diff_Type")],
    #annotation_col = meta_data[c("Alpha_sim","Beta_sim","Gamma_sim","Delta_sim","Ductal_sim","Acinar_sim","NEC_NET")],
    #annotation_col = meta_data[c("Alpha_sim","Beta_sim","Gamma_sim","Delta_sim","Ductal_sim","Acinar_sim","NEC_NET")],
    annotation_colors = aka3,
    annotation_legend = T,
    treeheight_col = 0,
    treeheight_row = 0,
    show_colnames = F,
    show_rownames = T,
    #color = colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(length(breaksList)),
    cluster_cols = T, cluster_rows = F
)

table(meta_data[,c("Y_Subtype","Diff_Type")])

### p_value

stats_t = as.data.frame(res_cor)
table(stats_t$`P-value`)

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
result <- confusionMatrix( as.factor(meta_data$Diff_Type), meta_data$Subtype)
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
