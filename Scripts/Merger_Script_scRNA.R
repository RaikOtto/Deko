library("stringr")

expr_raw = readRDS("~/Downloads/Tosti.Seurat.normalized.S78048.RDS")
meta_info = read.table("~/Deko_Projekt/Misc/Tosti_Metadaten.tsv", sep ="\t", header = T)
rownames(meta_info) = meta_info$Cell

meta_data = meta_info[colnames(expr_raw),]
subtype_vector = meta_data$Cluster
table(subtype_vector)

#candidates = which(subtype_vector %in% c("alpha","beta","gamma","delta","acinar-s","acinar-reg+","acinar-i","ductal","muc5b+ ductal"))
candidates = which(subtype_vector %in% c("Acinar-REG+","Acinar-i","MUC5B+ Ductal"))
expr_raw_tosti = expr_raw[,candidates]
meta_data_tosti = meta_info[colnames(expr_raw_tosti),]
subtype_vector_reduced_tosti = meta_data_tosti$Cluster
table(subtype_vector_reduced_tosti)

amount_genes = 300
amount_samples = 300

selected_samples = c()

for ( cell_type in unique(subtype_vector_reduced_tosti)){
    coords = which(meta_data_tosti$Cluster == cell_type )
    
    if (length(coords) >= amount_samples)
        coords = sample(coords, size = amount_samples)
    
    selected_samples = c(selected_samples, coords)
}
length(selected_samples)

expr_tosti = expr_raw_tosti[,selected_samples]
dim(expr_tosti)
meta_data_reduced_tosti = meta_info[colnames(expr_tosti),]
table(meta_data_reduced_tosti$Cluster)


###

meta_info = read.table("~/Deko_Projekt//Misc/Meta_information_scRNA.tsv",sep = "\t",header = T,stringsAsFactors = F)

rownames(meta_info) = meta_info$Sample

colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

expr_raw = read.table("~/Deko_Projekt/Data/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
#expr_raw = read.table("~/Deko_Projekt/Data/Cancer_Pancreas_Bulk_Array/Sato.S35.Ninon.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = paste("X",colnames(expr_raw)[no_match],sep ="")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
table(no_match)
meta_data = meta_info[colnames(expr_raw),]

subtype_vector = meta_data$Cluster
table(subtype_vector)

candidates = which(subtype_vector %in% c("Alpha","Beta","Gamma","Delta","Acinar","Ductal"))
expr_raw_baron = expr_raw[,candidates]
meta_data_baron = meta_info[colnames(expr_raw_baron),]
subtype_vector_reduced_baron = meta_data_baron$Cluster
table(subtype_vector_reduced_baron)

amount_genes = 300
amount_samples = 300

selected_samples = c()

for ( cell_type in unique(subtype_vector_reduced_baron)){
    coords = which(meta_data_baron$Cluster == cell_type )
    
    if (length(coords) >= amount_samples)
        coords = sample(coords, size = amount_samples)
    
    selected_samples = c(selected_samples, coords)
}
length(selected_samples)

expr_baron = expr_raw_baron[,selected_samples]
dim(expr_baron)
meta_data_reduced_baron = meta_info[colnames(expr_baron),]
table(meta_data_reduced_baron$Cluster)

### merge datasets 

merge_genes = intersect(rownames(expr_tosti),rownames(expr_baron))
#merge_genes = rownames(bam_data_1)
length(merge_genes)
table("INS" %in% merge_genes)
table("GCG" %in% merge_genes)
table("PPY" %in% merge_genes)
table("SST" %in% merge_genes)

new_mat = as.data.frame(
    cbind(
        expr_tosti[merge_genes,],
        expr_baron[merge_genes,]
    )
)
rownames(new_mat) = merge_genes

row_var = as.double(apply(new_mat, FUN = function(vec){return(var(vec))}, MARGIN = 1))
summary(row_var)
new_mat = new_mat[which( row_var >= 1),]
new_mat = new_mat[which( rowMeans(new_mat) >= 1),]

table(meta_data$Subtype)
dim(new_mat)

new_mat = new_mat[ rownames(new_mat)!="NA", ]
dim(new_mat)
new_mat[1:5,1:5]

write.table(new_mat[,], "~/Deko_Projekt/Data/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron_Alpha-i_Alpha-reg_Muc5+_Tosti_488_genes.tsv", sep ="\t", quote =F , row.names = T)
