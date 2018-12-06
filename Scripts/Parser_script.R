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

bam_data = read.table("~/Deko/Data/Merged_Dif_Baron_Hisc_Haber.tsv",sep ="\t", header = T, row.names = 1)
#bam_data = read.table("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Sadanandam_97.tsv",sep ="\t", header = T, row.names = 1)
rownames(bam_data) = str_to_upper( rownames( bam_data) )
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "^X","")
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "\\.","_")
bam_data = bam_data[,colnames(bam_data) %in% meta_info$Name]
dim(bam_data)

meta_data <<- meta_info[colnames(bam_data),]
subtype_vector = str_to_lower(meta_data$Subtype)
