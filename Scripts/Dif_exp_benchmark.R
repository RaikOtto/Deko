library("stringr")

meta_info = "~/Deko/Misc/Meta_information.tsv" %>% read.table(sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

# parse scRNA training data

i_mat = transcriptome_data
colnames(i_mat) = str_replace_all(colnames(i_mat),pattern = "^X","")

if (!exists(scRNA_training_mat))
    scRNA_training_mat = "~/Deko/Data/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Baron.tsv" %>% read.table(sep="\t",header = T,stringsAsFactors = F)
scRNA_training_mat[1:5,1:5]
scRNA_training_mat %>% typeof()
training_nr_marker_genes = 800

meta_data = meta_info[colnames(scRNA_training_mat),]
case_subtype = "Ductal"
groups = meta_data$Subtype
table(groups)
groups[groups == case_subtype] = "CASE"
groups[groups != "CASE"] = "CTRL"
design <- model.matrix(~0 + groups)
colnames(design) = c("Case","Ctrl")

vfit = lmFit(scRNA_training_mat,design)
contr.matrix = makeContrasts(
    contrast = Case - Ctrl,
    levels = design
)

vfit = contrasts.fit( vfit, contrasts = contr.matrix)
efit = eBayes(vfit)

result_t = topTable(
    efit,
    coef     = "contrast",
    number   = nrow(scRNA_training_mat),
    genelist = efit$genes,
    "none",
    sort.by  = "B",
    resort.by= NULL,
    p.value  = 1,
    lfc      = 0,
    confint  = FALSE
)
result_t$hgnc_symbol = rownames(result_t)
colnames(result_t) = c("Log_FC","Average_Expr","t","P_value","adj_P_value","B","HGNC")

result_t = result_t[c("HGNC","Log_FC","Average_Expr","P_value","adj_P_value")]
result_t = result_t[order(result_t$P_value, decreasing = FALSE),]
result_t$Log_FC = round(result_t$Log_FC, 1)
result_t$Average_Expr = round(result_t$Average_Expr, 1)
result_t = result_t[order(result_t$Log_FC,decreasing = TRUE),]

marker_genes = as.character(result_t$HGNC)
marker_genes = marker_genes[marker_genes %in% rownames(scRNA_training_mat)]
marker_genes = marker_genes[1:training_nr_marker_genes]

## prepping the prediction matrix


meta_data = meta_info[colnames(i_mat),]
marker_genes = marker_genes[marker_genes %in% rownames(i_mat)]

prediction_mat = i_mat[marker_genes,]
prediction_mat = apply(prediction_mat, FUN = function(vec) return(log(vec+1)), MARGIN = 1)
mki_67_pred = i_mat["MKI67",] %>% as.double
#prediction_mat = t(prediction_mat)
#prediction_mat = do.call("rbind",prediction_mat)

na_cols = apply(prediction_mat, FUN = function(vec) return(sum(is.na(vec))), MARGIN = 2)
na_col_indices = which(na_cols > 0)
if (length(na_col_indices) > 0)
    prediction_mat = prediction_mat[,-as.integer(na_col_indices)]

### correlation study

training_mat = scRNA_training_mat[marker_genes,]
training_mat = do.call("rbind",training_mat)
training_mat = apply(training_mat, FUN = function(vec) return(log(vec+1)), MARGIN = 1)
training_mat = t(training_mat)
colnames(training_mat) = marker_genes

if (length(na_col_indices) > 0)
    training_mat = training_mat[,-as.integer(na_col_indices)]

mki_67_training = scRNA_training_mat["MKI67",] %>% as.double
mki_67_training = log(mki_67_training + 1)

# training

rf_fit = glm(
    matrix(mki_67_training,ncol = 1) ~ .,
    data = data.frame(training_mat)
)
predicted_mki67 = predict(rf_fit, data.frame(prediction_mat)) %>% plogis  # training scores

# predicting

grading_vec = meta_info[rownames(prediction_mat),"Grading"] %>% str_replace_all( pattern = "G", "" ) %>% as.double
#cor.test(predicted_mki67,grading_vec)
cor.test(predicted_mki67,i_mat["MKI67",])
print(case_subtype)
