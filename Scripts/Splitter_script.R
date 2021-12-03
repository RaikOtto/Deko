library("stringr")

meta_info_maptor = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info_maptor) = meta_info_maptor$Sample
colnames(meta_info_maptor) = str_replace(colnames(meta_info_maptor),pattern = "\\.","_")

meta_info_maptor$OS_Tissue = as.double(str_replace(meta_info_maptor$OS_Tissue,pattern = ",","."))

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

matcher = match(meta_info_maptor$Sample,meta_info$Sample, nomatch = 0)
meta_info[matcher,"OS_Tissue"] = meta_info_maptor[matcher != 0,"OS_Tissue"]

table(meta_info$Site_of_primary,meta_info$NET_NEC_PCA)

### deconvolution

expr_raw = read.table("~/Dropbox/testproject/Datasets/PanNEN_only/Deconvolution/Deconvolution.6_studies.Exocrine.Absolute.PanNEN_only.Grading_binary.S232.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = TRUE)
expr_raw = read.table("~/Dropbox/testproject/Datasets/PanNEN_NEN/Deconvolution/Deconvolution.6_studies.Exocrine.Absolute.PanNEN_NEN.NET_NEC.S264.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = TRUE)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "\\.", "")
expr_raw[1:5,1:5]
dim(expr_raw)

table(expr_raw$Sample %in% meta_info$Sample)
meta_data = meta_info[expr_raw$Sample,]

candidates = meta_data$Sample[ 
    #meta_data$Study %in% c("Master","Charite")
    #meta_data$Site_of_primary %in% c("Pancreatic") #& meta_data$Primary_Metastasis %in% c("Primary")
    meta_data$Site_of_primary == "Pancreatic" #& meta_data$Primary_Metastasis %in% c("Primary")
    #meta_data$NET_NEC_PCA %in% c("NEC","NET")
    #!(meta_data$Histology_Metastasis %in% c("Outlier"))
]
length(candidates)
expr_raw = expr_raw[expr_raw$Sample %in% candidates,]
meta_data = meta_info[expr_raw$Sample,]
dim(meta_data)

write.table(expr_raw,"~/Dropbox/testproject/Datasets/PanNEN_only/Deconvolution/Deconvolution.6_studies.Exocrine.Absolute.PanNEN_only.Grading_binary.S232.tsv",sep ="\t",quote =F , row.names = TRUE)
write.table(expr_raw,"~/Dropbox/testproject/Datasets/PanNEN_only/Deconvolution/Deconvolution.6_studies.Exocrine.Absolute.PanNEN_only.NET_NEC.S230.tsv",sep ="\t",quote =F , row.names = TRUE)

### expression

#expr_raw = read.table("~/Deko_Projekt/Data/Publication_datasets/Combinations_PanNEN/6_studies.3473_genes.S238.tsv",sep="\t", stringsAsFactors =  F, row.names = 1, header = T, as.is = TRUE)
#expr_raw = read.table("~/Dropbox/testproject/Datasets/PanNEN_NEN/Expression/Expression_6_studies_3473_genes.Grading_binary.S264.tsv",sep="\t", stringsAsFactors =  F, row.names = 1, header = T, as.is = TRUE)
#expr_raw = read.table("~/Dropbox/testproject/Datasets/PanNEN_NEN/Expression/Expression_6_studies_3473_genes.NET_NEC.S264.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = TRUE)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "\\.", "")
expr_raw[1:5,1:5]
dim(expr_raw)

table(expr_raw$Sample %in% meta_info$Sample)
meta_data = meta_info[expr_raw$Sample,]

table(rownames(expr_raw) %in% meta_info$Sample)
meta_data = meta_info[expr_raw$Sample,]

candidates = meta_data$Sample[ 
    #meta_data$Study %in% c("Master","Charite")
    #meta_data$Site_of_primary %in% c("Pancreatic") #& meta_data$Primary_Metastasis %in% c("Primary")
    meta_data$Site_of_primary == "Pancreatic" #& meta_data$Primary_Metastasis %in% c("Primary")
    #meta_data$NET_NEC_PCA %in% c("NEC","NET")
    #!(meta_data$Histology_Metastasis %in% c("Outlier"))
]
length(candidates)
expr_raw = expr_raw[expr_raw$Sample %in% candidates,]
meta_data = meta_info[expr_raw$Sample,]
dim(meta_data)

write.table(expr_raw,"~/Dropbox/testproject/Datasets/PanNEN_only/Deconvolution.6_studies.Exocrine.Absolute.PanNEN_only.NET_NEC.S230.tsv",sep ="\t",quote =F , row.names = TRUE)
write.table(new_mat,"~/Deko_Projekt/Data/Publication_datasets/Combinations_PanNEN/6_studies.1000_genes.S238.tsv",sep ="\t",quote =F , row.names = TRUE)

meta_3 = meta_info[meta_info$Grading == "G3",]

table(meta_3$Site_of_primary, meta_3$NET_NEC_PCA)

                     MiNEN   NEC  NET   Unknown
Gastric/duodenal     2       16   3     0
Large_intestinal     11      16   0     0
Other                0       6    0     0
Pancreatic           0       30   33    2
Small_intestinal     0       1    1     0

                    Control  MiNEN NEC NET Unknown
charite                1     0     0   0   0
Gastric/duodenal       0     2     16  6   1
Hepatic                4     0     0   0   0
Intestinal             1     0     0   0   0
Large_intestinal       0     11    16  2   18
Other                  0     0     6   0   0
Pancreatic             104   0     30  214 127
Small_intestinal       0     0     1   8   82
