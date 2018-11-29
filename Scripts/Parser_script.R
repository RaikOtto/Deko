library("bseqsc")
library("stringr")
library("grid")     ## Need to attach (and not just load) grid package
#library("pheatmap")

#write.table(meta_info, "~/Deko/Misc/Meta_information.tsv", row.names = F, quote =F , sep = "\t")

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

###

bam_data = read.table("~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",sep ="\t", header = T, row.names = 1)
#bam_data = read.table("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Sadanandam_97.tsv",sep ="\t", header = T, row.names = 1)
rownames(bam_data) = str_to_upper( rownames( bam_data) )
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "^X","")
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "\\.","_")

meta_data <<- meta_info[colnames(bam_data),]
eset = new("ExpressionSet", exprs=as.matrix(bam_data));
#bam_data = bam_data[,meta_data$Subtype == "Non_functional"]

#hgnc_list = rownames(bam_data)
#hgnc_list_uni = unique(hgnc_list)
#source("~/Deko/Scripts/Variance_selection.R")

bam_data[1:5,1:5]
dim(bam_data)

#condition = (meta_info$Included == "Yes") & (meta_info$NEC_NET %in% c("NEC","NET")) #& (meta_info$Grading %in% c("G3"))
#include_list  = meta_info$Name[ condition ]
#exclude_list  = meta_info$Name[ ! condition ]

#bam_data = bam_data[, ( colnames(bam_data) %in% include_list )]
#dim(bam_data)
#bam_data[1:5,1:5]

#count_data = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Lawlor.tsv",sep ="\t", header = T, stringsAsFactors = F)
count_data = read.table("~/Deko/Data/Human_Mouse_HSC/Differentiated_Segerstolpe_Progenitor_Stanescu_ISC.tsv",sep ="\t", header = T, stringsAsFactors = F)

colnames(count_data) = str_replace(colnames(count_data), pattern = "\\.", "_")
colnames(count_data) = str_replace(colnames(count_data), pattern = "^X", "")
count_data[1:5,1:5]
