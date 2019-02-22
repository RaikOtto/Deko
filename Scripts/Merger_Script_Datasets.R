library("stringr")

#t1 = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.tsv"
t1_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.primary_only.tsv"
t2_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/GSE73338.tsv"
#t3_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.primary.pancreas.tsv"
#t4_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73339/GSE73339.tsv"

t1 = read.table(t1_path,sep="\t",header=T,row.names = 1, stringsAsFactors = F)
colnames(t1) = str_replace_all(colnames(t1),pattern = "^X", "")
t2 = read.table(t2_path,sep="\t",header=T,row.names = 1, stringsAsFactors = F)
colnames(t2) = str_replace_all(colnames(t2),pattern = "^X", "")
rownames(t2)
#t3 = read.table(t3_path,sep="\t",header=T,row.names = 1, stringsAsFactors = F)
#colnames(t3) = str_replace_all(colnames(t3),pattern = "^X", "")
#t4 = read.table(t4_path,sep="\t",header=T,row.names = 1, stringsAsFactors = F)
#colnames(t4) = str_replace_all(colnames(t4),pattern = "^X", "")

# variance selection

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name

### integrate

"INS" %in% rownames(t1)
"INS" %in% rownames(t2)
"INS" %in% rownames(t3)
"INS" %in% rownames(t4)

merge_genes = intersect(rownames(t1),rownames(t2))
merge_genes = intersect(merge_genes, rownames(t3))
merge_genes = intersect(merge_genes, rownames(t4))

length(merge_genes)
"INS" %in% merge_genes

new_mat = as.data.frame(
    cbind(
        t1[merge_genes,],
        t2[merge_genes,]#,
        #bam_data_3[merge_genes,]
    )
)
rownames(new_mat) = merge_genes
meta_data = meta_info[colnames(new_mat[,]),]
dim(new_mat)

row_var = as.double(apply(new_mat, FUN = function(vec){return(var(vec))}, MARGIN = 1))
new_mat = new_mat[which( row_var >= 1),]
new_mat = new_mat[which( rowMeans(new_mat) >= 1),]

table(meta_data$Subtype)
table(is.na(meta_data$Subtype))

dim(new_mat)

new_mat = new_mat[ rownames(new_mat)!="NA", ]

dim(new_mat)
new_mat[1:5,1:5]
write.table(new_mat[,], "~/Deko/Data/Wiedenmann_GSE73338.unfiltered.tsv", sep ="\t", quote =F , row.names = T)
#write.table(new_mat, "~/Deko/Data/HISC_HESC.tsv", sep ="\t", quote =F , row.names = T)

meta_data = meta_info[colnames(bam_data_1),]
#meta_info[colnames(bam_data_1),"Subtype"] = "HISC"
bam_data_1 = bam_data_1[, meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta","Acinar","Ductal")]
#bam_data_1 = bam_data_1[, meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta")]