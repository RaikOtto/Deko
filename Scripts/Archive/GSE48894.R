library(stringr)
library(xbioc)
library(cowplot)
library(reshape2)

i_files = list.files("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894//",full.names = T)
i_files_names = list.files("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894//",full.names = F)
i_files_names = str_replace_all(i_files_names, pattern = "_counts\\.txt","")

first_file = read.table(i_files[1],sep ="\t", header = T,row.names = 1)

for (i in 2:length(i_files)){
    new_sample = read.table(i_files[i],sep ="\t", header = T,row.names = 1)
    first_file = cbind(first_file,new_sample)
}

colnames(first_file) = i_files_names
first_file[1:5,1:5]

library("org.Hs.eg.db")

### hgnc

expr_raw = first_file
gene_ids = row.names(expr_raw)

mapping_t = read.table("~/Deko/Misc/Mapping_Gene_id_HGNC.tsv",sep = "\t", stringsAsFactors = F, header = T, fill = TRUE)
colnames(mapping_t) = str_replace_all(colnames(mapping_t), pattern = "\\.","_")
mapping_t = as.data.frame(mapping_t)
hgnc_match = match(rownames(expr_raw),mapping_t$NCBI_Gene_ID)
hgnc_list = mapping_t[hgnc_match,"Approved_symbol"]

### variane selection

expr_raw[1:5,1:5]

write.table(expr_raw,"~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894.tsv",sep="\t",quote = F, row.names = T)
write.table(expr,"~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894.vis.tsv",sep="\t",quote = F, row.names = T)
