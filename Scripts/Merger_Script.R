library("stringr")

# Segerstolpe integration

bam_data_1 = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Segerstolpe.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
colnames(bam_data_1) = str_replace(colnames(bam_data_1),pattern = "\\.","_")
#bam_data_1 = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Segerstolpe.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
bam_data_1[1:5,1:5]

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
table(colnames(bam_data_1) %in% meta_info$Name)

meta_data = meta_info[colnames(bam_data_1),]
bam_data_1 = bam_data_1[, meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta")]
rownames(bam_data_1) = str_to_upper(rownames(bam_data_1))
dim(bam_data_1)

# Stanescu

bam_data_2 = read.table("~/Deko/Data/Mouse_progenitor_pancreas_scRNA/Stanescu.tsv" , sep ="\t" ,header = T, stringsAsFactors = F )
rownames(bam_data_2) = str_to_upper(rownames(bam_data_2))
rownames(bam_data_2)[rownames(bam_data_2) == "INS1"] = "INS"
bam_data_2[1:5,1:5]

# HISC

bam_data_3 = read.table("~/Deko/Data/Human_Mouse_HSC/Haber.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
colnames(bam_data_3) = str_replace_all(colnames(bam_data_3) , pattern = "^X", "")
rownames(bam_data_3) = str_to_upper(rownames(bam_data_3))
rownames(bam_data_3)[rownames(bam_data_3) == "INS-IGF2" ] = "INS"
bam_data_3[1:5,1:5]

meta_data = meta_info[colnames(bam_data_3),]
bam_data_3 = bam_data_3[, meta_data$Subtype %in% c("HISC")]
dim(bam_data_3)

### integrate

merge_genes = intersect(rownames(bam_data_1),rownames(bam_data_2))
merge_genes = intersect(merge_genes, rownames(bam_data_3))
length(merge_genes)

new_mat = as.data.frame(
    cbind(
        bam_data_1[merge_genes,],
        bam_data_2[merge_genes,],
        bam_data_3[merge_genes,]
    )
)
rownames(new_mat) = merge_genes
new_mat[1:5,1:5]

row_var = as.double(apply(new_mat, FUN = function(vec){return(var(vec))}, MARGIN = 1))
new_mat = new_mat[row_var >= 1,]
new_mat[1:5,1:5]
#write.table(new_mat, "~/Deko/Data/Merged_Segerstolpe_Prog_Hisc.tsv", sep ="\t", quote =F , row.names = T)
