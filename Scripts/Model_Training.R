count_data = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Segerstolpe.tsv",sep ="\t", header = T, stringsAsFactors = F)
#count_data = read.table("~/Deko/Data/Merged_Segerstolpe_Prog_Hisc.tsv",sep ="\t", header = T, stringsAsFactors = F)

colnames(count_data) = str_replace(colnames(count_data), pattern = "\\.", "_")
colnames(count_data) = str_replace(colnames(count_data), pattern = "^X", "")
count_data[1:5,1:5]

subtypes = meta_info[colnames(count_data),"Subtype"]
names(subtypes) = meta_info[colnames(count_data),"Name"]
table(subtypes)

#marker_genes = read.table(
#    "~/Deko/Misc/Baron_pancreas_marker.tsv",
#    sep = "\t",
#    header = T,
#    stringsAsFactors = F
#)
delimiter = 1:100
#marker_genes = marker_genes[delimiter,]

pancreasMarkers = list()

#pancreasMarkers = list(
#    "Alpha" = marker_genes$Alpha[ (marker_genes$Alpha %in% rownames(count_data))],
#    "Beta" = marker_genes$Beta[ (marker_genes$Beta != "") & (marker_genes$Beta %in% rownames(count_data))],
#    "Gamma" = marker_genes$Gamma[ (marker_genes$Gamma != "")& (marker_genes$Gamma %in% rownames(count_data))],
#    "Delta" = marker_genes$Delta[ (marker_genes$Delta != "")& (marker_genes$Delta %in% rownames(count_data))]#,
    #"Ductal" = marker_genes$Ductal[ (marker_genes$Ductal != "")& (marker_genes$Ductal %in% rownames(count_data))],
    #"Acinar" = marker_genes$Acinar[ (marker_genes$Acinar != "")& (marker_genes$Acinar %in% rownames(count_data))]#,
    #"Progenitor" = marker_genes$Progenitor[ (marker_genes$Progenitor != "")& (marker_genes$Progenitor %in% rownames(count_data))],
    #"HISC" = marker_genes$HESC[(marker_genes$HISC != "")& (marker_genes$HISC %in% rownames(count_data))]
    #E17.5 = names(E17.5),
#)

cands = names(subtypes)[ str_to_upper( subtypes ) %in% str_to_upper( names(pancreasMarkers))]
length(cands)
count_data = count_data[, cands]
subtypes = subtypes[cands]
dim(count_data)
table(subtypes)

### normalization

row_var = apply(count_data, FUN = var, MARGIN = 1)
col_var = apply(count_data, FUN = var, MARGIN = 2)
table(row_var == 0)
table(col_var == 0)
count_data = count_data[row_var != 0,col_var != 0]
count_data = count_data[rowSums(count_data) >= 1,]
dim(count_data)

table(as.character(unlist(pancreasMarkers)) %in% rownames(count_data))
for(marker in names(pancreasMarkers)){
    genes = as.character(unlist(pancreasMarkers[marker]))
    genes = genes[genes %in% rownames(count_data)]
    pancreasMarkers[marker] = list(genes)
}

names(pancreasMarkers) = str_to_lower(names(pancreasMarkers))
eislet = new("ExpressionSet", exprs = as.matrix(count_data))

sub_list = str_to_lower(subtypes)
names(sub_list) = names(subtypes)
fData(eislet) = data.frame( sub_list  )
pData(eislet) = data.frame( sub_list )
names(pancreasMarkers)[!(names(pancreasMarkers) %in% sub_list)]

B = bseqsc_basis(
    eislet,
    pancreasMarkers,
    clusters = 'sub_list',
    samples = colnames(exprs(eislet)),
    ct.scale = FALSE
)
plotBasis(B, pancreasMarkers, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')

eset = new("ExpressionSet", exprs=as.matrix(count_data));
nr_permutations = 100
fit = bseqsc_proportions(eset, B, verbose = FALSE, absolute = T, log = F, perm = nr_permutations)

res_coeff = t(fit$coefficients)
res_coeff_mat = as.double(unlist(res_coeff))
res_coeff_mat = as.data.frame(matrix(res_coeff_mat,ncol = ncol(res_coeff), nrow = nrow(res_coeff)))
rownames(res_coeff_mat) = rownames(res_coeff)
colnames(res_coeff_mat) = colnames(res_coeff)
res_cor   = fit$stats

res_coeff[ is.na(res_coeff) ] = 0.0
res_cor[ is.na(res_cor) ] = 0.0

p_value_threshold = 0.05
not_sig_samples = rownames(res_cor)[res_cor[,"P-value"] > p_value_threshold]
not_sig_samples

for (subtype in unique(subtypes)){
    
    subtype_index = as.integer(sapply(subtypes, FUN = function(vec){
        return(which(
            str_to_upper(vec) == str_to_upper(as.character(colnames(res_coeff)))))
    }))
}

i <<- 0
self_score = lapply( as.character(subtypes), FUN = function(sub){
        index = which(str_to_upper(sub) == str_to_upper(as.character(colnames(res_coeff))))
        i <<- i + 1
        return(res_coeff_mat[i,index])
})
self_score = as.double(as.character(unlist(self_score)))
self_score = as.data.frame(cbind(self_score, subtypes))
max_vals = aggregate(as.double(self_score$self_score), by = list(self_score$subtypes), FUN = max )
max_vec = max_vals$x
names(max_vec) = str_to_lower(max_vals$Group.1)

library("ROCR")
library("caret")

subtype = "delta"
subtypes = str_to_lower(subtypes)
pred_vec = subtypes
pred_vec[pred_vec == subtype] = 1
pred_vec[pred_vec != 1] = 0
pred_vec = as.double(pred_vec)

pred_score = res_coeff_mat[subtype]
pred <- prediction(predictions = res_coeff_mat[subtype], labels =  pred_vec )
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc.perf)

opt.cut = function(roc.perf, pred){
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
          cutoff = p[[ind]])
    }, roc.perf@x.values, roc.perf@y.values, pred@cutoffs)
}
print(opt.cut(roc.perf, pred))

result <- confusionMatrix( as.factor(pred_vec), as.factor(meta_data$Subtype))
result 


model = list(B, max_vec)
saveRDS(model,"~/ArtDeco/inst/Models/Four_differentiation_stages_Lawlor.RDS")