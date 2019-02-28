library("stringr")

# scRNA integration

bam_data_1 = read.table("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.primary_only.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
#bam_data_1 = read.table("~/Deko/Data/Human_Mouse_HSC/Haber.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
colnames(bam_data_1) = str_replace(colnames(bam_data_1),pattern = "\\.","_")
colnames(bam_data_1) = str_replace(colnames(bam_data_1),pattern = "^X","")
bam_data_1[1:5,1:5]
rownames(bam_data_1) = str_to_upper(rownames(bam_data_1))
dim(bam_data_1)

# variance selection

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
meta_data = meta_info[colnames(bam_data_1),]
#meta_info[colnames(bam_data_1),"Subtype"] = "HISC"
bam_data_1 = bam_data_1[, meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta","Acinar","Ductal")]
#bam_data_1 = bam_data_1[, meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta")]

meta_data = meta_info[colnames(bam_data_1),]
dim(bam_data_1)
table(meta_data$Subtype)


all_subs = which(meta_data$Subtype == "Acinar")
red_subs = sample(all_subs,300, replace = F)
exclude = all_subs[!( all_subs %in% red_subs)]
index = 1:ncol(bam_data_1)
index = index[!( index %in% exclude)]
bam_data_1 = bam_data_1[,index]

# HISC

bam_data_2 = read.table("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/GSE73338.ki67.Grading.Primary.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
colnames(bam_data_2) = str_replace_all(colnames(bam_data_2) , pattern = "^X", "")
rownames(bam_data_2) = str_to_upper(rownames(bam_data_2))
rownames(bam_data_2)[rownames(bam_data_2) %in% c("INS1") ] = "INS"
bam_data_2[1:5,1:5]

meta_data = meta_info[colnames(bam_data_2),]
bam_data_2 = bam_data_2[, meta_data$Subtype %in% c("HISC")]
dim(bam_data_2)

### integrate

merge_genes = intersect(rownames(bam_data_1),rownames(bam_data_2))
#merge_genes = rownames(bam_data_1)
length(merge_genes)
table("INS" %in% rownames(bam_data_1))
table("INS" %in% rownames(bam_data_2))
table("INS" %in% merge_genes)

new_mat = as.data.frame(
    cbind(
        bam_data_1[merge_genes,],
        bam_data_2[merge_genes,]
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
write.table(new_mat[,], "~/Deko/Data/Wiedenmann_Scarpa_GSE73338.tsv", sep ="\t", quote =F , row.names = T)
