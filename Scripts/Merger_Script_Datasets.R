library("stringr")

t1_path = "~/MAPTor_NET/BAMs_new/Groetzinger.New.HGNC.S51.tsv"
t2_path = "~/MAPTor_NET/BAMs_new/RepSet_S84.HGNC.tsv"
t3_path = "~/MAPTor_NET/BAMs_new/Master.pre.S39.HGNC.tsv"
#t4_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73339/GSE73339.tsv"

t1 = read.table(t1_path,sep="\t",header=T,row.names = 1, stringsAsFactors = F)
colnames(t1) = str_replace_all(colnames(t1),pattern = "^X", "")
dim(t1)
t1[1:5,]

t2 = read.table(t2_path,sep="\t",header=T,row.names = 1, stringsAsFactors = F)
colnames(t2) = str_replace_all(colnames(t2),pattern = "^X", "")
t2 = t2[,str_detect(colnames(t2),pattern = "_",negate = T)]
dim(t2)
t2[1:5,]

t2 = t2[,83:84]

t3 = read.table(
    file=t3_path,
    row.names = 1,
    sep="\t", stringsAsFactors =  F, header = T)
colnames(t3) = str_replace(colnames(t3), pattern = "^X\\.", "")
t3[1:5,1:5]
no_match = (colnames(t3) %in% meta_info$Sample) == F
colnames(t3)[no_match] = paste("X",colnames(t3)[no_match],sep ="")
colnames(t3) = str_replace(colnames(t3), pattern = "XX", "")
no_match = (colnames(t3) %in% meta_info$Sample) == F
no_match

t3 = t3[,str_detect(colnames(t3),pattern = "_",negate = T)]
dim(t3)
t3[1:5,1:5]

# variance selection

meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample

### integrate

"INS" %in% rownames(t1)
"INS" %in% rownames(t2)
"INS" %in% rownames(t3)
#"INS" %in% rownames(t4)

merge_genes = intersect(rownames(t1),rownames(t2))
#merge_genes = intersect(merge_genes, rownames(t3))
#merge_genes = intersect(merge_genes, rownames(t4))

length(merge_genes)
"INS" %in% merge_genes

new_mat = as.data.frame(
    cbind(
        t1[merge_genes,],
        t2[merge_genes,]#,
        #t3[merge_genes,]
    )
)
rownames(new_mat) = merge_genes
meta_data = meta_info[colnames(new_mat[,]),]
dim(new_mat)

row_var = as.double(apply(new_mat, FUN = function(vec){return(var(vec))}, MARGIN = 1))
new_mat = new_mat[which( row_var >= 1),]
new_mat = new_mat[which( rowMeans(new_mat) >= 1),]

new_mat = new_mat[,colnames(new_mat) %in% meta_info[meta_info$Grading != "","Sample"]]
"132502" %in% colnames(new_mat)

table(meta_data$NEC_NET_Ori)
table(meta_data$Study)
table(meta_data$Grading)

dim(new_mat)

new_mat = new_mat[ rownames(new_mat)!="NA", ]

dim(new_mat)
new_mat[1:5,1:5]
#write.table(new_mat[,], "~/MAPTor_NET/BAMs_new/RepSet_S95.HGNC.tsv", sep ="\t", quote =F , row.names = T)
#write.table(expr_raw_2, "~/MAPTor_NET/BAMs_new/RepSet_S56.HGNC.tsv", sep ="\t", quote =F , row.names = T)

#meta_data = meta_info[colnames(bam_data_1),]
#meta_info[colnames(bam_data_1),"Subtype"] = "HISC"
#bam_data_1 = bam_data_1[, meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta","Acinar","Ductal")]
#bam_data_1 = bam_data_1[, meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta")]