library("stringr")

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample

<<<<<<< HEAD
# scRNA integration

new_mat = bam_data_1 = read.table("~/Dropbox/testproject/Datasets/PanNEN_only/Deconvolution/6_studies.Exocrine.Absolute.PanNEN_only.Grading_binary.S232.tsv" , sep ="\t" ,header = T, stringsAsFactors = F)

new_mat = bam_data_1 = read.table("~/Dropbox/testproject/Datasets/PanNEN_only/Expression/4_studies.400_genes.S128.Grading_binary.tsv" , sep ="\t" ,header = T, stringsAsFactors = F)
allowed_genes = colnames(bam_data_1)

bam_data_1 = read.table("~/Deko_Projekt/Data/Publication_datasets/Combinations_PanNEN/Charite_Scarpa_Master_Diedisheim.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
=======
#
bam_data_1 = read.table("~/MAPTor_NET/BAMs_new/RepSet_S103.NEN.HGNC.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
#bam_data_1 = read.table("~/Deko_Projekt/Data/Publication_datasets/Charite.S23.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
>>>>>>> 9bc00275ee3721ca59b01ac60dbef7c0671edf84
colnames(bam_data_1) = str_replace(colnames(bam_data_1),pattern = "\\.","_")
colnames(bam_data_1) = str_replace(colnames(bam_data_1),pattern = "^X","")
rownames(bam_data_1) = str_to_upper(rownames(bam_data_1))
summary(as.double(bam_data_1["INS",]))
summary(as.double(bam_data_1["SST",]))

# variance selection

bam_data_1 = bam_data_1[allowed_genes,]

col_names = colnames(bam_data_1)
row_names = rownames(bam_data_1)
bam_data_1 = matrix(as.double(as.character(unlist(bam_data_1))),ncol = ncol(bam_data_1))
rownames(bam_data_1) = row_names
colnames(bam_data_1) = col_names
summary(bam_data_1["INS",])

meta_data = meta_info[colnames(bam_data_1),]
table(meta_data$Study)

# bam data 2
#bam_data_2 = read.table("~/Deko_Projekt/Data/Cancer_Pancreas_Bulk_Array/Diedisheim.S66.HGNC.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
#bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/Master.S20.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
#bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/Scarpa.S29.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
#bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/Diedisheim.S62.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
<<<<<<< HEAD
#bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/Sadanandam.S29.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/Missiaglia.S75.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
=======
bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/Sadanandam.S29.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
#bam_data_2 = read.table("~/Deko_Projekt/Data/Publication_datasets/Missiaglia.S75.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)

>>>>>>> 9bc00275ee3721ca59b01ac60dbef7c0671edf84
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
meta_data = meta_info[rownames(new_mat[,]),]
dim(new_mat)

###

row_var = as.double(apply(new_mat, FUN = function(vec){return(var(vec))}, MARGIN = 1))
summary(row_var)
selection = order(row_var, decreasing = )[1:400]
#new_mat = new_mat[which( row_var >= 1),]
new_mat = new_mat[selection,]
dim(new_mat)

###

meta_data = meta_info[colnames(new_mat),]
#new_mat = new_mat[ rownames(new_mat)!="NA", ]
new_mat = new_mat[ ,which(meta_data$Site_of_primary == "Pancreatic") ]

dim(new_mat)
new_mat[1:5,1:5]
<<<<<<< HEAD
=======
new_mat = new_mat[,order(colnames(new_mat))]
new_mat = new_mat[grep(rownames(new_mat), pattern = "\\.",invert = TRUE),]
>>>>>>> 9bc00275ee3721ca59b01ac60dbef7c0671edf84

new_mat = new_mat[ which(meta_data$Study %in% c("Charite","Diedisheim","Master","Scarpa")) ,]
meta_data = meta_info[colnames(new_mat),]

meta_data = meta_info[rownames(new_mat),]
table(meta_data$Study)

meta_data$Grading_binary = meta_data$Grading
meta_data$Grading_binary[meta_data$Grading_binary %in% c("G1","G2")] = "1"
meta_data$Grading_binary[meta_data$Grading_binary %in% c("G3")] = "0"

table(meta_data$NET_NEC_PCA)

new_mat = t(new_mat)
new_mat = as.data.frame(new_mat )

new_mat$Grading = meta_data$Grading
new_mat = new_mat[new_mat$Grading !="Unknown",]

#new_mat$Grading_binary = meta_data$Grading_binary
#new_mat = new_mat[new_mat$Grading_binary!="Unknown",]

new_mat$NET_NEC = meta_data$NET_NEC_PCA
new_mat = new_mat[new_mat$NET_NEC != "Unknown",]

table(new_mat$NET_NEC)
new_mat$NET_NEC[new_mat$NET_NEC == "NEC"] = "0"
new_mat$NET_NEC[new_mat$NET_NEC == "NET"] = "1"

dim(new_mat)

write.table(new_mat, "~/Dropbox/testproject/Datasets/PanNEN_only/Expression/4_studies.400_genes.S128.Grading_tertiary.tsv", sep ="\t", quote =F , row.names = FALSE)
write.table(new_mat, "~/Dropbox/testproject/Datasets/PanNEN_only/Expression/6_studies.400_genes.S232.Grading_tertiary.tsv", sep ="\t", quote =F , row.names = FALSE)

write.table(new_mat, "~/Dropbox/testproject/Datasets/PanNEN_only/Deconvolution/4_studies.Exocrine.Absolute.PanNEN_only.NET_NEC.S128.tsv", sep ="\t", quote =F , row.names = FALSE)
write.table(new_mat, "~/Dropbox/testproject/Datasets/PanNEN_only/Expression/6_studies.400_genes.S232.Grading_tertiary.tsv", sep ="\t", quote =F , row.names = FALSE)

write.table(new_mat, "~/Dropbox/testproject/Datasets/PanNEN_only/Expression/4_studies.400_genes.S128.Grading_binary.tsv", sep ="\t", quote =F , row.names = FALSE)
write.table(new_mat, "~/Dropbox/testproject/Datasets/PanNEN_only/Expression/6_studies.400_genes.S238.Grading_binary.tsv", sep ="\t", quote =F , row.names = FALSE)

write.table(new_mat, "~/Dropbox/testproject/Datasets/PanNEN_only/Expression/4_studies.400_genes.S128.NET_NEC.tsv", sep ="\t", quote =F , row.names = TRUE)
write.table(new_mat, "~/Dropbox/testproject/Datasets/PanNEN_only/Expression/6_studies.400_genes.S207.NET_NEC.tsv", sep ="\t", quote =F , row.names = TRUE)

deco_mat_grading = deco_mat[deco_mat$Sample %in% rownames(new_mat),]

meta_data = meta_info[deco_mat_grading$Sample,]
meta_data$Grading_binary = meta_data$Grading
meta_data$Grading_binary[meta_data$Grading_binary %in% c("G1","G2")] = "1"
meta_data$Grading_binary[meta_data$Grading_binary %in% c("G3")] = "0"
deco_mat_grading$Grading_binary = meta_data$Grading_binary

meta_data = meta_info[deco_mat_net_nec$Sample,]
deco_mat_net_nec$NET_NEC = meta_data$NET_NEC_PCA
deco_mat_net_nec$NET_NEC[deco_mat_net_nec$NET_NEC == "NEC"] = "0"
deco_mat_net_nec$NET_NEC[deco_mat_net_nec$NET_NEC == "NET"] = "1"
table(deco_mat_net_nec$NET_NEC)
