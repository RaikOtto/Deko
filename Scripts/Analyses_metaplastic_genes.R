library("stringr")

t_normal = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet_S56_Cibersort_Baron.tsv",sep ="\t", header = T, row.names = 1)
dim(t_normal)
t_meta_only = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet_S56_Cibersort_Baron_metaplastic.tsv",sep ="\t", header = T, row.names = 1)
dim(t_meta_only)
t_no_meta = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet_S56_Cibersort_Baron_non_metaplastic.tsv",sep ="\t", header = T, row.names = 1)
dim(t_no_meta)

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

t_normal = t_normal[order(rownames(t_normal)),]
t_meta_only = t_meta_only[order(rownames(t_meta_only)),]
t_no_meta = t_no_meta[order(rownames(t_no_meta)),]
identical(rownames(t_normal),rownames(t_meta_only))
identical(rownames(t_no_meta),rownames(t_meta_only))

meta_data = meta_info[rownames(t_no_meta),]

t_normal=t_normal[meta_data$NET_NEC_PCA == "NET",]
t_no_meta=t_no_meta[meta_data$NET_NEC_PCA == "NET",]
t_meta_only=t_meta_only[meta_data$NET_NEC_PCA == "NET",]

mean(t_normal$Correlation)
mean(t_meta_only$Correlation)
mean(t_no_meta$Correlation)

cor.test(t_normal$Correlation,t_meta_only$Correlation)
cor.test(t_no_meta$Correlation,t_meta_only$Correlation)

mean(t_no_meta$P_value)
mean(t_meta_only$P_value)
mean(t_normal$P_value)

cor.test(t_no_meta$P_value,t_meta_only$P_value)
t_no_meta
t_no_meta$X
t_meta_only
t_normal
meta_data = meta_info[t_normal$X,]
identical(t_normal$X,t_no_meta$X)
meta_data_normal = meta_info[t_normal$X,]
t_normal = t_normal[meta_data_normal$NET_NEC_PCA %in% "NEC",]
t_normal
meta_data_no_meta = meta_info[t_no_meta$X,]
t_no_meta = t_normal[meta_data_no_meta$NET_NEC_PCA %in% "NEC",]
meta_data_meta_only = meta_info[t_meta_only$X,]
t_meta_only = t_meta_only[meta_data_meta_only$NET_NEC_PCA %in% "NEC",]
t_no_meta = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet_S56_Cibersort_Baron_non_metaplastic.tsv",sep ="\t", header = T)
t_no_meta = t_no_meta[meta_data_no_meta$NET_NEC_PCA %in% "NEC",]
mean(t_no_meta$P_value)
mean(t_normal$P_value)
mean(t_meta_only$P_value)
mean(t_meta_only$Correlation)
mean(t_no_meta$Correlation)
mean(t_normal$Correlation)
t.test(t_meta_only$Correlation, t_no_meta$Correlation)
library("stringr")
