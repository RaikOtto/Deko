library("stringr")

balanced.centroid = read.table( "~/Koop_Klinghammer/Misc//balanced.centroid.txt", header=TRUE, row.names=1, sep="\t",stringsAsFactors = F)
balanced.centroid_importance = sort(rowSums(abs(balanced.centroid)), decreasing = T)
balanced.centroid = balanced.centroid[ match(names(balanced.centroid_importance),rownames(balanced.centroid)),]

### Preparation

rownames(pure_data) = str_to_upper(rownames(pure_data))
rownames(pure_data) = str_replace_all(rownames(pure_data), pattern = "\\_", "")
rownames(pure_data) = str_replace_all(rownames(pure_data), pattern = "-", "")
rownames(balanced.centroid) = str_to_upper(rownames(balanced.centroid))
rownames(balanced.centroid) = str_replace_all(rownames(balanced.centroid), pattern = "\\_", "")
rownames(balanced.centroid) = str_replace_all(rownames(balanced.centroid), pattern = "-", "")
colnames(pure_data) = str_replace_all(colnames(pure_data), pattern = "^X", "")

col_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 2)
row_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 1)

pure_data = pure_data[row_var > 0, col_var > 0]
dim(pure_data)

###

expr = pure_data
source("~/Koop_Klinghammer/Scripts/Classification_scripts.R")
table( rownames(expr) %in% rownames(balanced.centroid) )

### centroid classification

pub_cor <<- matrix( as.double(), ncol = length( colnames( balanced.centroid )  ) )
expr2bc = centroid2expr( balanced.centroid[,], expr )
colnames(expr2bc$correlation) = c("Sample","Subtype","Correlation","P_value")
class_data = as.data.frame(expr2bc$correlation)

meta_match = match( class_data$Sample, meta_info$Sample_ID, nomatch = 0 )
meta_info$Subtype[meta_match] = as.character( class_data$Subtype )
meta_info$P_value[meta_match] = as.double( as.character( class_data$P_value ) )
meta_info$Sig[meta_match] = meta_info$P_value[meta_match  < 0.05] = "TRUE"
meta_info$Sig[is.na(meta_info$Sig)] = ""
meta_info$Included[meta_info$Sig != "TRUE"] = "FALSE"

write.table(meta_info,"~/Koop_Klinghammer/Misc/Meta_Information.tsv",sep ="\t",quote =F,row.names =F)
m_match = match(colnames(pure_data), meta_info$Sample_ID, nomatch = 0)
meta_data = meta_info[m_match,]
rownames(meta_data) = meta_data$Sample_ID

#s_match = match( meta_data$Name, meta_info$Name, nomatch = 0 )
#meta_info$Subtype_Kmeans[s_match] = meta_data$Subtype_Kmeans
