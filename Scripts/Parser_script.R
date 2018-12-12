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

bam_data = read.table("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Groetzinger_57.tsv",sep ="\t", header = T, row.names = 1)
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "^X","")
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "\\.","_")

meta_data <<- meta_info[colnames(bam_data),]

bam_data = bam_data[,str_to_lower(meta_data$Subtype) %in% c("alpha","beta","gamma","delta","acinar","ductal")]
meta_data <<- meta_info[colnames(bam_data),]
subtype_vector = str_to_lower(meta_data$Subtype)
table(subtype_vector)
sum(table(subtype_vector))
dim(bam_data)

expression_training_mat = bam_data

#meta_data = Determine_differentiation_stage("~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv")
meta_data = Determine_differentiation_stage(bam_Data)

visualization_data = read.table("~/Deko/Data/Visualization_PANnen.tsv", sep ="\t", stringsAsFactors = F)
colnames(expression_data) = str_replace_all(colnames(expression_data),"^X","")
expression_data[1:5,1:5]
ARTDECO::create_heatmap_differentiation_stages(
    expression_data = visualization_data,
    
)
