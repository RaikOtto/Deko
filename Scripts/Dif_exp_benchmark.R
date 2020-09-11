library("stringr")
library("limma")

expr_raw = read.table("~/Deko_Projekt/Data/Bench_data/Riemer_Scarpa.S69.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info$NEC_NET = meta_info$NEC_NET_PCA

meta_data = meta_info[colnames(expr_raw),]

# parse scRNA training data

groups = meta_data$Histology
groups[groups != "Pancreatic"] = "Other"
case_subtype = "Pancreatic"

table(groups)
groups[groups == case_subtype] = "CASE"
groups[groups != "CASE"] = "CTRL"
design <- model.matrix(~0 + groups)
colnames(design) = c("Case","Ctrl")

vfit = lmFit(
    expr_raw,
    design
)
contr.matrix = makeContrasts(
    contrast = Case - Ctrl,
    levels = design
)

vfit = contrasts.fit( vfit, contrasts = contr.matrix)
efit = eBayes(vfit)

result_t = topTable(
    efit,
    coef     = "contrast",
    number   = nrow(expr_raw),
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
result_t = result_t[order(abs(result_t$Log_FC),decreasing = TRUE),]
result_t$HGNC[1:800]

marker_genes = as.character(result_t$HGNC)
marker_genes = marker_genes[marker_genes %in% rownames(expr_raw)]
marker_genes = marker_genes[1:training_nr_marker_genes]

write.table(result_t[1:800,],"~/Deko_Projekt/Results/dif_expr_RepSet_S69_Pancreatics_Minus_Other_Tissue.tsv",sep ="\t", row.names=T, quote = F)

## prepping the prediction matrix

sig_t = result_t[1:800,]
dim(sig_t)

table( hisc_genes[1:400] %in% sig_t$HGNC[1:400] ) 

table( ductal_genes[1:400] %in% sig_t$HGNC[1:400] ) 

table( ductal_genes %in% hisc_genes ) 

marker_genes = data.frame(
    "HISC" = hisc_genes,
    "Ductal" = ductal_genes
)

write.table(marker_genes,"~/Deko_Projekt/Results/HISC_Ductal_marker_genes.tsv",sep ="\t", row.names=T, quote = F)
