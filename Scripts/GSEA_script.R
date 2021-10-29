library("stringr")

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
sort(table(meta_info$Sample),decreasing = TRUE)

path_transcriptome_file = "~/MAPTor_NET/BAMs_new/RepSet_S103.HGNC.DESeq2.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = F,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

table(match(colnames(expr_raw), meta_info$Sample, nomatch = 0) != 0)

meta_data = meta_info[colnames(expr_raw),]
candidates = meta_data$Sample[ 
    meta_data$Histology_Primary %in% c("Pancreatic") & meta_data$Primary_Metastasis %in% c("Primary")
]
expr = expr_raw[,candidates]
meta_data = meta_info[colnames(expr),]
dim(meta_data)

### cls file

cohort_vec = meta_data$NET_NEC_PCA
cohort_vec = str_replace(cohort_vec," ","_")
table(cohort_vec)

first_entity = unique(cohort_vec)[1]
second_entity = unique(cohort_vec)[2]

first_line = as.character(c(ncol(expr),length(unique(cohort_vec)),1))
first_line[4:(length(cohort_vec))] = ""
second_line = as.character(c("#",first_entity,second_entity))
second_line[4:(length(cohort_vec))] = ""

cls_file = rbind(first_line,second_line,cohort_vec)

write.table(cls_file, "~/Deko_Projekt/GSEA/RepSet.S39.primrary_pancreas_only_NEC_vs_NET.DESeq2.cls",row.names = F, quote = F, sep ="\t",col.names = FALSE)
write.table(cls_file, "~/Deko_Projekt/GSEA/RepSet.S39.primrary_pancreas_only_NEC_vs_NET.DESeq2.cls.tsv",row.names = F, quote = F, sep ="\t",col.names = FALSE)

### gct file

first_line = "#1.0"
first_line[2:(ncol(expr)+2)] = ""
second_line = as.character(c(nrow(expr),ncol(expr)))
second_line[3:(ncol(expr)+2)] = ""
third_line = c("Name","Description",colnames(expr))
gct_matrix = cbind(rownames(expr),rownames(expr),expr)

output_gct = rbind(first_line,second_line, third_line, gct_matrix)
output_gct[1:5,1:5]

write.table(output_gct, "~/Deko_Projekt/GSEA/RepSet.S39.primrary_pancreas_only_NEC_vs_NET.DESeq2.gct",row.names = F, quote = F, sep ="\t",col.names = FALSE)
write.table(output_gct, "~/Deko_Projekt/GSEA/RepSet.S39.primrary_pancreas_only_NEC_vs_NET.DESeq2.gct.tsv",row.names = F, quote = F, sep ="\t",col.names = FALSE)
