names(pancreasMarkers) = str_to_lower(names(pancreasMarkers))
eislet = new("ExpressionSet", exprs = as.matrix(count_data))

sub_list = str_to_lower(subtypes)
names(sub_list) = names(subtypes)
fData(eislet) = data.frame( sub_list  )
pData(eislet) = data.frame( sub_list )

B = readRDS("~/ArtDeco/inst/Models/Four_differentiation_stages.RDS")
"
B = bseqsc_basis(
eislet,
pancreasMarkers,
clusters = 'sub_list',
samples = colnames(exprs(eislet)),
ct.scale = FALSE
)
plotBasis(B, pancreasMarkers, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')
"

### RUN VARIANCE SELECTION FIRST

bam_data = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Muraro.tsv",sep ="\t", header = T, stringsAsFactors = F)

subtypes = meta_info[colnames(bam_data),"Subtype"]
cands = which(!is.na(subtypes)) & (subtypes %in% c("Alpha","Beta","Gamma","Delta"))
meta_data = meta_info[colnames(bam_data)[cands],]

meta_data$Y_Subtype = as.character(subtypes[cands])
expr_raw = bam_data[,cands]
bam_data = bam_data[,cands]
eset = new("ExpressionSet", exprs=as.matrix(bam_data));

fit = bseqsc_proportions(eset, B, verbose = TRUE, absolute = T, log = F, perm = 100)
source("~/Deko/Scripts/Utility_script.R")
#meta_data$Diff_Type =  c("Differentiated","Progenitor","HSC")[maxi]

#expr_raw = bam_data
#colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

pheatmap::pheatmap(
    t(res_coeff),
    #cor(expr),
    annotation_col = meta_data[c("Y_Subtype","Diff_Type")],
    #annotation_col = meta_data[c("Alpha_sim","Beta_sim","Gamma_sim","Delta_sim","Ductal_sim","Acinar_sim","NEC_NET")],
    #annotation_col = meta_data[c("Alpha_sim","Beta_sim","Gamma_sim","Delta_sim","Ductal_sim","Acinar_sim","NEC_NET")],
    annotation_colors = aka3,
    annotation_legend = T,
    treeheight_col = 0,
    treeheight_row = 0,
    show_colnames = F,
    show_rownames = T,
    #color = colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(length(breaksList)),
    cluster_cols = T, cluster_rows = F
)

table(meta_data[,c("Y_Subtype","Diff_Type")])

### p_value

stats_t = as.data.frame(res_cor)
table(stats_t$`P-value`)

###


expected    = paste( c(as.character( t_rel$Expected[ t_rel$Expected!= ""] )), collapse = ", ",sep = "" )
nr_expected = length((as.character(unlist(str_split(expected,", ")))))
print(paste(c("Expected:",nr_expected),collapse = "",sep = "") )

true_positives = as.character(unlist(str_split(t_rel$True_positive,", ")))
true_positives = true_positives[ true_positives != ""  ]
true_positives = true_positives[ !is.na(true_positives) ]
nr_true_positives = length(true_positives)
print(paste(c("TP:",as.character(nr_true_positives)),collapse = "",sep = ""))

false_negatives = as.character(unlist(str_split(t_rel$False_negative,", ")))
false_negatives = false_negatives[ false_negatives!= ""]
false_negatives = false_negatives[ !is.na(false_negatives) ]
nr_false_negatives = length( false_negatives )
print(paste(c("FN: ",as.character(nr_false_negatives)),collapse = "",sep = ""))

false_positives = t_rel$False_positive[t_rel$False_positive!= ""]
nr_false_positive = length((as.character(unlist(str_split(false_positives,",")))))
print(paste(c("FP: ",nr_false_positive),collapse = "",sep = ""))

nr_true_negatives  = nrow(t_rel)**2 - nrow(t_rel) - nr_false_negatives

Sensitivity = round(nr_true_positives / (nr_true_positives + nr_false_negatives),3) * 100
print(paste(c("Sensitivity: ",Sensitivity),collapse = "",sep = ""))

PPV = round( nr_true_positives / ( nr_true_positives + nr_false_positive ),3) * 100
print(paste(c("PPV: ",PPV),collapse = "",sep = ""))

F1 = round( (2* nr_true_positives) / ((2*nr_true_positives) + nr_false_positive + nr_false_negatives) ,3) * 100
print(paste(c("F1: ",F1),collapse = "",sep = ""))