library("stringr")

# Sadanandam, Missiaglia
# Califano

t1_path = "~/MAPTor_NET/BAMs_new/RepSet_S57.HGNC.tsv"  # Groetzinger, Scarpa, Master
t2_path = "~/MAPTor_NET/BAMs_new/Exp_tpm_bon_qgp_knih6.HGNC.tsv"
#t2_path = "~/Deko_Projekt/Data/Cancer_Pancreas_Bulk_Array/Missiaglia/GSE73339.all.tsv"
#t3_path = "~/Deko_Projekt/Data/Cancer_Pancreas_Bulk_Array/Sadanandam/Sadanandam.tsv"

t1 = read.table(t1_path,sep="\t",header=T,row.names = 1, stringsAsFactors = F)
colnames(t1) = str_replace_all(colnames(t1),pattern = "^X", "")
dim(t1)
t1[1:5,]

t2 = read.table(t2_path,sep="\t",header=T,row.names = 1, stringsAsFactors = F)
colnames(t2) = str_replace_all(colnames(t2),pattern = "^X", "")
colnames(t2) = str_replace_all(colnames(t2),pattern = "^TPM_", "")
#t2 = t2[,str_detect(colnames(t2),pattern = "_",negate = T)]
dim(t2)
t2[1:5,]
no_match = (colnames(t2) %in% meta_info$Sample) == F
colnames(t2)[no_match] = paste("X",colnames(t2)[no_match],sep ="")
colnames(t2) = str_replace(colnames(t2), pattern = "XX", "")
no_match = (colnames(t2) %in% meta_info$Sample) == F
no_match


#t3 = read.table(
#    file=t3_path,
#    row.names = 1,
#    sep="\t", stringsAsFactors =  F, header = T)
#colnames(t3) = str_replace(colnames(t3), pattern = "^X\\.", "")
#t3[1:5,1:5]

# variance selection

meta_info = read.table("~/Deko_Projekt//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
#meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample

### integrate

"INS" %in% rownames(t1)
"INS" %in% rownames(t2)
#"INS" %in% rownames(t3)
#"INS" %in% rownames(t4)

merge_genes = intersect(str_to_upper(rownames(t1)),str_to_upper(rownames(t2)))
#merge_genes = intersect(merge_genes, rownames(t3))
#merge_genes = intersect(merge_genes, rownames(t4))

length(merge_genes)
"INS" %in% merge_genes

new_mat = as.data.frame(
    cbind(
        t1[merge_genes,],
        t2[merge_genes,]#,
        #t3[merge_genes,]
    )
)
rownames(new_mat) = merge_genes
meta_data = meta_info[colnames(new_mat[,]),]
dim(new_mat)

row_var = as.double(apply(new_mat, FUN = function(vec){return(var(vec))}, MARGIN = 1))
new_mat = new_mat[which( row_var >= 1),]
new_mat = new_mat[which( rowMeans(new_mat) >= 1),]

"132502" %in% colnames(new_mat)

table(meta_data$NEC_NET)
table(meta_data$Study)
table(meta_data$Grading)

dim(new_mat)

new_mat = new_mat[ rownames(new_mat)!="NA", ]
new_mat = new_mat[!(is.na(new_mat[,1])) , ]

dim(new_mat)
new_mat[1:5,1:5]
#write.table(new_mat[,], "~/MAPTor_NET/BAMs_new/RepSet_S57_Bon_control.HGNC.tsv", sep ="\t", quote =F , row.names = T)

#meta_info$Albumin = rep("",nrow(meta_info))
#meta_info[colnames(expr_raw),"Albumin"] = as.double(expr_raw["ALB",])
#write.table(meta_info, "~/MAPTor_NET/Misc/Meta_information.tsv", sep ="\t", quote =F , row.names = F)

#meta_data = meta_info[colnames(bam_data_1),]
#meta_info[colnames(bam_data_1),"Subtype"] = "HISC"
#bam_data_1 = bam_data_1[, meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta","Acinar","Ductal")]
#bam_data_1 = bam_data_1[, meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta")]

meta_info_map = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info_deko = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info_deko) = meta_info_deko$Sample

matcher = match(meta_info_map$Sample,meta_info_deko$Sample, nomatch = 0)

meta_info_deko$NEC_NET = rep("",nrow(meta_info_deko))
meta_info_deko[matcher,"NEC_NET"] = meta_info_map[matcher != 0,"NEC_NET_PCA"]

meta_data = meta_info_deko[colnames(expr_raw),]
met_tab = meta_data[which(meta_data$NEC_NET == ""),]

#write.table(meta_info_deko, "~/Deko_Projekt/Misc/Meta_information.tsv", sep ="\t", quote =F , row.names = T)
