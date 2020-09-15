# prep

library("stringr")
meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name

# scRNA integration

bam_data_1 = read.table("~/Deko_Projekt/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Baron_4.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
colnames(bam_data_1) = str_replace(colnames(bam_data_1),pattern = "\\.","_")
colnames(bam_data_1) = str_replace(colnames(bam_data_1),pattern = "^X","")
rownames(bam_data_1) = str_to_upper(rownames(bam_data_1))
summary(as.double(bam_data_1["INS",]))

# variance selection

col_names = colnames(bam_data_1)
row_names = rownames(bam_data_1)
bam_data_1 = matrix(as.double(as.character(unlist(bam_data_1))),ncol = ncol(bam_data_1))
rownames(bam_data_1) = row_names
colnames(bam_data_1) = col_names
summary(bam_data_1["INS",])

meta_data = meta_info[colnames(bam_data_1),]
#meta_info[colnames(bam_data_1),"Subtype"] = "HISC"
#bam_data_1 = bam_data_1[, meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta","Acinar","Ductal")]
bam_data_1 = bam_data_1[, meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta")]

meta_data = meta_info[colnames(bam_data_1),]
dim(bam_data_1)
table(meta_data$Subtype)

# HISC

bam_data_2 = read.table("~/Deko_Projekt/Data/Human_Mouse_HSC/Haber.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
colnames(bam_data_2) = str_replace_all(colnames(bam_data_2) , pattern = "^X", "")
rownames(bam_data_2) = str_to_upper(rownames(bam_data_2))
table("INS" %in% rownames(bam_data_2))
table("PPY" %in% rownames(bam_data_2))
table("PNP" %in% rownames(bam_data_2))
table("SST" %in% rownames(bam_data_2))
table("REG1A" %in% rownames(bam_data_2)) # acinar
table("PRSS2" %in% rownames(bam_data_2)) # acinar
table("CELA3A" %in% rownames(bam_data_2)) # acinar
table("PTP" %in% rownames(bam_data_2)) # acinar
table("TMSB4X" %in% rownames(bam_data_2)) # ductal

meta_data = meta_info[colnames(bam_data_2),]
bam_data_2 = bam_data_2[, meta_data$Subtype %in% c("HISC")]
dim(bam_data_2)

### integrate

merge_genes = intersect(rownames(bam_data_1),rownames(bam_data_2))
#merge_genes = rownames(bam_data_1)
length(merge_genes)
table("INS" %in% rownames(bam_data_1))
table("GCG" %in% rownames(bam_data_1))
table("PPY" %in% merge_genes)
table("SST" %in% merge_genes)

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
dim(new_mat)

new_mat = new_mat[ rownames(new_mat)!="NA", ]
dim(new_mat)
new_mat[1:5,1:5]
#write.table(bam_data_1, "~/Deko/Data/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.tsv", sep ="\t", quote =F , row.names = T)
write.table(new_mat[,], "~/Deko_Projekt/Data/Alpha_Beta_Gamma_Delta_Hisc_Baron.tsv", sep ="\t", quote =F , row.names = T)

### splitter

# scRNA integration

bam_data_1 = read.table("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Wiedenmann.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
colnames(bam_data_1) = str_replace(colnames(bam_data_1),pattern = "\\.","_")
colnames(bam_data_1) = str_replace(colnames(bam_data_1),pattern = "^X","")
rownames(bam_data_1) = str_to_upper(rownames(bam_data_1))
summary(as.double(bam_data_1["INS",]))

meta_data = meta_info[colnames(bam_data_1),]
table(meta_data$Grading)

meta_data = meta_data[meta_data$Grading != "G0",]
bam_data_1 = bam_data_1[,meta_data$Name]
dim(bam_data_1)

#write.table(bam_data_1,"~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Wiedenmann.Scarpa.tsv",sep="\t",row.names = T, col.names= T,quote=F)

###

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_Bulk/Meta_info.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name

bam_data = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_Bulk/Fadista.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "\\.","_")
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "^X","")
rownames(bam_data) = str_to_upper(rownames(bam_data))
bam_data = as.data.frame(bam_data)
dim(bam_data)
summary(as.double(bam_data["INS",])) # sanity check

meta_data = meta_info[colnames(bam_data),]
table(meta_data$Grading)
dim(meta_data)
#table(meta_data$Grading)
meta_data = meta_data[meta_data$Study == "Groetzinger",]
meta_data = meta_data[meta_data$Histology == "Pancreatic_NEN",]
meta_data = meta_data[meta_data$Grading %in% c("G0"),]
bam_data = bam_data[,meta_data$Name]
dim(bam_data)
table(meta_data$Grading)

#write.table(meta_data,"~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Wiedenmann_selection_meta_data.tsv",sep="\t",row.names = T, col.names= T,quote=F)
write.table(bam_data,"~/Deko/Data/Bench_data/Controls.S11.tsv",sep="\t",row.names = T, col.names= T,quote=F)

### parse GSE73339

mapping_t = read.table("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73339/Mapping_table.tsv",sep = "\t",header = T,stringsAsFactors = F)
mapping_t = as.data.frame(mapping_t)
rownames(mapping_t) = mapping_t$ID
mapping_t[1:5,1:2]

table(rownames(bam_data) %in% mapping_t$ID)
mapping_t = mapping_t[rownames(bam_data),]
bam_data = bam_data[rownames(mapping_t),]
hgnc_list = mapping_t$HGNC
hgnc_list = str_replace_all(hgnc_list,pattern = " ","")
