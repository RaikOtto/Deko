library("stringr")
library("biomaRt")
library("org.Hs.eg.db")

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

###

input_files = list.files("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/Raw/",full.names = T,pattern = ".txt")
input_files_names = list.files("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/Raw/",full.names = F, pattern = ".txt")
input_files_names = str_replace_all(input_files_names,pattern = "_counts.txt","")

i_file = read.table(input_files[1],sep="\t",row.names = 1, header = T, stringsAsFactors = F)
i_file = as.data.frame(i_file[,1])

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
entrez = rownames(i_file)
hgnc_entrez_mapping = getBM(
    attributes = c('entrezgene', 'hgnc_symbol'), 
    filters = 'entrezgene', 
    values = entrez, 
    mart = ensembl
)

hgnc_entrez_mapping[1:100,]
mapping = match(
    as.integer(rownames(i_file)),
    as.integer(hgnc_entrez_mapping[,1]),
    nomatch = 0
)

hgnc_list = str_to_upper(hgnc_entrez_mapping[
    mapping,
    2
])

for( i_file_next in input_files[2:length(input_files)]){
    i_file_new = read.table(i_file_next,sep="\t",row.names = 1, header = T, stringsAsFactors = F)
    #print(i_file_new[10309,])
    i_file = cbind(
        i_file,
        as.integer(i_file_new[,1])
    )
}
dim(i_file)
summary(as.integer(i_file[10309,]))
i_file = as.matrix(i_file,ncol=212)
colnames(i_file) = input_files_names
i_file = as.data.frame(i_file)
i_file[1:5,1:5]
#i_file_save = i_file

# variance selection

as.integer(i_file[rownames(i_file) == "INS",])
summary(as.integer(i_file[rownames(i_file) == "INS",]))

row_names = rownames(i_file)
col_names = colnames(i_file)
i_file = matrix(as.integer(as.character(unlist(i_file))),ncol=length(col_names))

colnames(i_file) = col_names
rownames(i_file) = row_names
i_file[1:5,1:5]

write.table(i_file,"~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.tsv",sep ="\t",quote = F, row.names = T)
