library("stringr")

# scRNA integration

bam_data_1 = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Baron_human.tsv" , sep ="\t" ,header = T, row.names = 1, stringsAsFactors = F)
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

# Stanescu

bam_data_2 = read.table("~/Deko/Data/Mouse_progenitor_pancreas_scRNA/Stanescu.tsv" , sep ="\t" ,header = T, stringsAsFactors = F ,row.names = 1)
rownames(bam_data_2) = str_to_upper(rownames(bam_data_2))

#rownames(bam_data_2) = str_replace_all(rownames(bam_data_2),pattern= "(\\.)|(-)|(_)","")
rownames(bam_data_2)[rownames(bam_data_2) == "INS1"] = "INS"
dim(bam_data_2)

# HISC

bam_data_3 = read.table("~/Deko/Data/Human_Mouse_HSC/Haber.tsv" , sep ="\t" ,header = T, stringsAsFactors = F, row.names = 1)
colnames(bam_data_3) = str_replace_all(colnames(bam_data_3) , pattern = "^X", "")
rownames(bam_data_3) = str_to_upper(rownames(bam_data_3))
#rownames(bam_data_3) = str_replace_all(rownames(bam_data_3),pattern= "(\\.)|(-)|(_)","")
rownames(bam_data_3)[rownames(bam_data_3) %in% c("INS-IGF2","INSIGF2") ] = "INS"
bam_data_3[1:5,1:5]

meta_data = meta_info[colnames(bam_data_3),]
bam_data_3 = bam_data_3[, meta_data$Subtype %in% c("HISC")]
dim(bam_data_3)

### integrate

merge_genes = intersect(rownames(bam_data_1),rownames(bam_data_2))
#merge_genes = rownames(bam_data_1)
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
write.table(new_mat[,], "~/Deko/Data/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron_progenitor_stanescu_hisc_haber.tsv", sep ="\t", quote =F , row.names = T)
#write.table(new_mat, "~/Deko/Data/HISC_HESC.tsv", sep ="\t", quote =F , row.names = T)
