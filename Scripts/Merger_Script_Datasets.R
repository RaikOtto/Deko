library("stringr")

# Sadanandam, Missiaglia
# Califano
meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv", sep ="\t", header  =T, as.is = T)
rownames(meta_info) = meta_info$Sample

t1_path = "~/MAPTor_NET/BAMs_new/Master/Master_new.S34.HGNC.tsv"  # Groetzinger, Scarpa, Master
t2_path = "~/MAPTor_NET/BAMs_new/Groetzinger/Groetzinger_new.S39.HGNC.tsv"
t2_path = "~/MAPTor_NET/BAMs_new/RepSet_S96.HGNC.tsv"
t3_path = "~/MAPTor_NET/BAMs_new/Scarpa/Scarpa_new.S29.HGNC.tsv"

t1 = read.table(t1_path,sep="\t",header=T,row.names = 1, stringsAsFactors = F)
colnames(t1) = str_replace_all(colnames(t1),pattern = "^X", "")
dim(t1)
t1[1:5,]

no_match = (colnames(t1) %in% meta_info$Sample) == F
colnames(t1)[no_match] = paste("X",colnames(t1)[no_match],sep ="")
colnames(t1) = str_replace(colnames(t1), pattern = "XX", "")
no_match = which((colnames(t1) %in% meta_info$Sample) == F)
no_match

t2 = read.table(t2_path,sep="\t",header=T,row.names = 1, stringsAsFactors = F)
colnames(t2) = str_replace_all(colnames(t2),pattern = "^X", "")
t2 = t2[,c("132502","132502")]
dim(t2)
t2[1:5,]

no_match = (colnames(t2) %in% meta_info$Sample) == F
colnames(t2)[no_match] = paste("X",colnames(t2)[no_match],sep ="")
colnames(t2) = str_replace(colnames(t2), pattern = "XX", "")
no_match = which((colnames(t1) %in% meta_info$Sample) == F)
no_match

#t3 = read.table(
#    file=t3_path,
#    row.names = 1,
#    sep="\t", stringsAsFactors =  F, header = T)
#colnames(t3) = str_replace(colnames(t3), pattern = "^X\\.", "")
#colnames(t3) = str_replace(colnames(t3), pattern = "^X", "")

#t3[1:5,]

# variance selection

### integrate

"INS" %in% rownames(t1)
"INS" %in% rownames(t2)
#"INS" %in% rownames(t3)
#"INS" %in% rownames(t4)

merge_genes = intersect(str_to_upper(rownames(t1)),str_to_upper(rownames(t2)))
#merge_genes = intersect(merge_genes, rownames(t3))
#merge_genes = intersect(merge_genes, rownames(t4))

length(merge_genes)
"INS" %in% merge_genes

new_mat = as.data.frame(
    cbind(
        t1[merge_genes,],
        t2[merge_genes,],
        t3[merge_genes,]
    )
)
rownames(new_mat) = merge_genes

meta_data = meta_info[colnames(new_mat[,]),]
meta_data$NEC_NET
dim(new_mat)

row_var = as.double(apply(new_mat, FUN = function(vec){return(var(vec))}, MARGIN = 1))
new_mat = new_mat[which( row_var >= 1),]
new_mat = new_mat[which( rowMeans(new_mat) >= 1),]

"132502" %in% colnames(new_mat)

table(meta_data$NEC_NET)
table(meta_data$Study)
table(meta_data$Grading)

dim(new_mat)

new_mat = new_mat[ rownames(new_mat)!="NA", ]
new_mat = new_mat[!(is.na(new_mat[,1])) , ]

dim(expr_raw)
expr_raw[1:5,1:5]
write.table(expr_raw[,], "~/MAPTor_NET/BAMs_new/RepSet_S103.HGNC.tsv", sep ="\t", quote =F , row.names = T)

expr_raw = expr_raw[,meta_data$Study == "Charite"]

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv", sep ="\t", header  =T, as.is = T)
rownames(meta_info) = meta_info$Sample

matcher = match(colnames(new_mat),meta_info$Sample, nomatch = 0)
colnames(new_mat)[matcher == 0]
colnames(new_mat)[matcher == 0]  = str_replace(colnames(new_mat)[matcher == 0], "^X","")
matcher = match(colnames(new_mat),meta_info$Sample, nomatch = 0)
colnames(new_mat)[matcher == 0]

meta_data = meta_info[colnames(new_mat),]
meta_data$Included_MAPTOR

new_mat = new_mat[,meta_data$Included_MAPTOR == "Yes"]
new_mat = new_mat[,meta_data$Study %in% c("Scarpa","Charite")]
dim(new_mat)

######

meta_data = meta_info[colnames(t1),]

t3 = t1[,0]
t3 = rep(0,nrow(t1))
matcher = match(rownames(t1),rownames(t2), nomatch = 0)
t3 = matrix(rep(0,nrow(t1)), ncol = 1, nrow = nrow(t1))
t3[matcher != 0,] = t2[matcher,1]
rownames(t3) = rownames(t1)
colnames(t3) = "132502"
mean_vals = apply(t1[matcher == 0,meta_data$NEC_NET_PCA == "NEC"], FUN= mean, MARGIN = 1)
t3[matcher == 0,1] = mean_vals

t1 = cbind(t1,t3)
t1[1:5,1:5]
dim(t1)
