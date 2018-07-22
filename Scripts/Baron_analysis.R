library("stringr")

i_files = list.files("~/Deko/Data/", pattern = "GSM22", full.names = T)

b_files <<- matrix(as.character(), ncol = 20128)

for (i in 1:4){
    print(i)
    i_mat = read.table( i_files[i], sep =",", header = T )  
    b_files <<- rbind(b_files, i_mat)
}

cell_types = b_files$assigned_cluster
table(cell_types)

b_files = b_files[, 4:ncol(b_files)]

b_files[1:5,1:5]

#write.table(b_files, "~/Deko/Data/baron_data.tsv", sep ="\t", row.names = F, col.names = T, quote = F)
#write.table(cell_types, "~/Deko/Data/baron_subtypes.tsv", sep ="\t", row.names = F, col.names = T, quote = F)

dim(b_files)

### vis

library(tsne)

tsne(b_files[1:10,])

### estimate

estimate_mat = read.table("~/MAPTor_NET/")

### Segerstolpe

data_mat = read.table("~/Deko/Data/pancreas_refseq_rpkms_counts_3514sc.txt", sep ="\t", header = T, stringsAsFactors = F)
colnames(data_mat) = str_replace_all(colnames(data_mat), pattern = "^X", "")
data_mat[1:5,1:5]

seg_meta = read.table("~/Deko/Data/Segerstolpe_Meta_info.tsv", sep ="\t", header = T)
s_match = match(colnames(s_t), seg_meta$Extract.Name, nomatch = 0)

cell_type = as.character(seg_meta$Characteristics.cell.type.)[s_match]
cell_type = str_replace_all(cell_type, pattern = " cell", "")
table(cell_type)
