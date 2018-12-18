library("bseqsc")
library("stringr")
library("grid")     ## Need to attach (and not just load) grid package

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

###

bam_data = read.table(
    "~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",
    sep ="\t",
    header = T
)
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "^X","")
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "\\.","_")
bam_data[1:5,1:5]

meta_data <<- meta_info[colnames(bam_data),]
#write.table(expr,"~/Downloads/Visualization_PANnen.tsv", sep = "\t", quote = F, row.names = T)

# training stuff
#add_deconvolution_training_model(training_data = "~/Deko/Data/Merged_Segerstolpe_Prog_Hisc.tsv",model_name = "Alpha_Beta_Gamma_Delta_SSegerstolpe_Progenitor_Stanescu_Hisc_Haber")

meta_data$Subtype = str_to_lower(meta_data$Subtype)
meta_data$Subtype[ meta_data$Subtype == "hisc"] = "hsc"
table(meta_data$Subtype)
bam_data = bam_data[,str_to_lower(meta_data$Subtype) %in% c("alpha","beta","gamma","delta","progenitor","hsc")]
meta_data <<- meta_info[colnames(bam_data),]
meta_data$Subtype = str_to_lower(meta_data$Subtype)
meta_data$Subtype[ meta_data$Subtype == "hisc"] = "hsc"
subtype_vector = str_to_lower(meta_data$Subtype)
table(subtype_vector)
sum(table(subtype_vector))
dim(bam_data)

expression_training_mat = bam_data


