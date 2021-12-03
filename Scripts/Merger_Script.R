# prep

library("stringr")

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample

# scRNA integration

bam_data_1 = read.table("~/Deko_Projekt/Data/Publication_datasets/NEN/Charite.S40.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
colnames(bam_data_1) = str_replace(colnames(bam_data_1),pattern = "\\.","_")
colnames(bam_data_1) = str_replace(colnames(bam_data_1),pattern = "^X","")
rownames(bam_data_1) = str_to_upper(rownames(bam_data_1))
summary(as.double(bam_data_1["INS",]))
summary(as.double(bam_data_1["SST",]))

# variance selection

col_names = colnames(bam_data_1)
row_names = rownames(bam_data_1)
bam_data_1 = matrix(as.double(as.character(unlist(bam_data_1))),ncol = ncol(bam_data_1))
rownames(bam_data_1) = row_names
colnames(bam_data_1) = col_names
summary(bam_data_1["INS",])

meta_data = meta_info[colnames(bam_data_1),]

# bam data 2

#bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/NEN/Diedisheim.S66.HGNC.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
#bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/NEN/Master.S34.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
#bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/Scarpa.S29.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
#bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/Sadanandam.S29.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/Missiaglia.S75.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
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

meta_data_2 = meta_info[colnames(bam_data_2),]
dim(bam_data_2)
dim(meta_data_2)

table(colnames(bam_data_2) %in% colnames(bam_data_1))
#bam_data_2 = bam_data_2[,which(meta_data_2$Clusters %in% c("EEC-Progenitor (Neurog3+)","MucinDuctalProgenitor"))]
#meta_data_2 = meta_info_2[colnames(bam_data_2),]
#table(meta_data_2$Clusters)

### integrate

merge_genes = intersect(rownames(bam_data_1),rownames(bam_data_2))
#merge_genes = rownames(bam_data_1)
length(merge_genes)
table("INS" %in% rownames(bam_data_1))
table("GCG" %in% rownames(bam_data_1))
table("PPY" %in% merge_genes)
table("SST" %in% merge_genes)

bam_data_1 = new_mat = as.data.frame(
    cbind(
        bam_data_1[merge_genes,],
        bam_data_2[merge_genes,]
    )
)
rownames(new_mat) = merge_genes
meta_data = meta_info[colnames(new_mat[,]),]
dim(new_mat)

row_var = as.double(apply(new_mat, FUN = function(vec){return(var(vec))}, MARGIN = 1))
summary(row_var)
new_mat = new_mat[which( row_var >= 1),]
#new_mat = new_mat[which( rowMeans(new_mat) >= 1),]

new_mat = new_mat[ rownames(new_mat)!="NA", ]
#new_mat = new_mat[ ,which(meta_data$Histology_Primary == "Pancreatic") ]
dim(new_mat)
new_mat[1:5,1:5]
new_mat = new_mat[,order(colnames(new_mat))]
new_mat = new_mat[grep(rownames(new_mat), pattern = "\\.",invert = TRUE),]

meta_data = meta_info[colnames(new_mat),]
meta_data$Grading_binary = meta_data$Grading
meta_data$Grading_binary[meta_data$Grading_binary %in% c("G1","G2")] = "G1_G2"
new_mat = new_mat[meta_data$Grading_binary !="Unknown",]

new_mat = new_mat[meta_data$NET_NEC_PCA %in% c("NET","NEC"),]

new_mat = t(new_mat)
new_mat = as.data.frame(new_mat )
meta_data = meta_info[rownames(new_mat),]
#new_mat$Grading_binary = meta_data$Grading_binary
new_mat$NET_NEC = meta_data$NET_NEC_PCA
#write.table(new_mat[,], "~/Dropbox/testproject/Datasets/Expression_6_studies_3388_genes.S330.Grading_binary.tsv", sep ="\t", quote =F , row.names = T)
#write.table(new_mat, "~/Dropbox/testproject/Datasets/Expression_6_studies_3388_genes.S330.NET_NEC.tsv", sep ="\t", quote =F , row.names = T)

