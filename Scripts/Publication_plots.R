
library("stringr")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))


expr_raw = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
#cands = meta_info$Name[meta_info$Study == "Groetzinger"]
#cands = cands[cands %in% colnames(expr_raw)]
#expr_raw = expr_raw[,  cands]
dim(expr_raw)
expr_raw[1:5,1:5]

### Prep

library("grid")

#meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_data = meta_info[ match(colnames(expr_raw), meta_info$Name), ]
rownames(meta_data) = meta_data$Name
df_map = match(colnames(expr_raw), meta_data$Name)
meta_data$INS = as.double(expr_raw[ which( rownames(expr_raw) == "INS"),df_map])
meta_data$GCG = as.double(expr_raw[ which( rownames(expr_raw) == "GCG"),df_map])
meta_data$PPY = as.double(expr_raw[ which( rownames(expr_raw) == "PPY"),df_map])
meta_data$SST = as.double(expr_raw[ which( rownames(expr_raw) == "SST"),df_map])
meta_data$MEN1_exp = as.double(expr_raw[ which( rownames(expr_raw) == "MEN1"),df_map])
meta_data$Location = str_replace_all(meta_data$Location, pattern = "-", "_")
meta_data$MKI67 = as.double(expr_raw[ which( rownames(expr_raw) == "MKI67"),])
meta_data$MKI67 = meta_data$MKI67 + ( sign(min(meta_data$MKI67)) *  min(meta_data$MKI67) )
meta_data$MKI67 = meta_data$MKI67/ max(meta_data$MKI67)
meta_data$Marker_Genes = meta_data$INS + meta_data$GCG + meta_data$PPY + meta_data$SST
meta_data$Marker_Genes = meta_data$Marker_Genes + ( sign(min(meta_data$Marker_Genes)) *  min(meta_data$Marker_Genes) )
meta_data$Marker_Genes = meta_data$Marker_Genes/ max(meta_data$Marker_Genes)

##
meta_data$Location[str_detect(meta_data$Location, pattern = "_Met")] = "Metastasis"
aka3 = list(
  Histology   = c(
    Pancreatic_NEN = "BLACK",
    Colorectal_NEN = "Orange",
    Small_intestinal_NEN = "Yellow",
    Gastric_NEN = "purple",
    Liver = "Darkgreen",
    CUP = "pink"),
  Deco_type = c(alpha = "Blue",beta = "Yellow",gamma = "Orange",delta = "Purple", e13.5 = "Brown", hsc = "white"),
  NEC_NET = c(NEC= "red", NET = "blue", "NA" = "white"),
  Study = c(Groetzinger = "brown", Scarpa = "darkgreen"),
  MKI67 = c(high = "White", medium = "gray", low = "black"),
  INS = c(high = "White", low = "yellow"),
  Correlation = c(high = "Red", medium = "Yellow", low = "Green"),
  Deco_similarity = c(low = "Red", medium = "Yellow", high = "Green"),
  ImmuneScore = c(high = "White", medium = "gray", low = "Black"),
  StromalScore = c(high = "White", medium = "gray", low = "Black"),
  #Alpha_sim = c(low = "white", medium = "yellow", high = "blue",not_sig = "gray"),
  #Beta_sim = c(low = "white", medium = "yellow", high = "Yellow",not_sig = "gray"),
  #Gamma_sim = c(low = "white", medium = "yellow", high = "orange",not_sig = "gray"),
  #Delta_sim = c(low = "white", medium = "yellow", high = "Purple",not_sig = "gray"),
  #Ductal_sim = c(low = "white", medium = "yellow", high = "Black",not_sig = "gray"),
  #Acinar_sim = c(low = "white", medium = "yellow", high = "Brown",not_sig = "gray"),
  Alpha_sim = c(low = "white", medium = "yellow", high = "blue"),
  Beta_sim = c(low = "white", medium = "Yellow", high = "darkgreen"),
  Gamma_sim = c(low = "white", medium = "yellow", high = "orange"),
  Delta_sim = c(low = "white", medium = "yellow", high = "Purple"),
  Ductal_sim = c(low = "white", medium = "yellow", high = "Black"),
  Acinar_sim = c(low = "white", medium = "yellow", high = "Brown"),
  Progenitor_sim = c(low = "white", medium = "yellow", high = "orange",not_sig = "gray"),
  Differentiated_sim = c(low = "white", medium = "yellow", high = "darkgreen",not_sig = "gray"),
  HSC_sim = c(low = "white", medium = "yellow", high = "darkred",not_sig = "gray"),
  Marker_Genes = c(high = "White", medium = "Yellow", low = "Red"),
  Functionality = c( Unknown = "White",Functional = "green", Non_Functional="red"),
  Grading = c( G1 = "Green",G2 = "Yellow", G3 = "Red", G0 = "white"),
  Included = c(Yes = "green", No = "red"),
  Chemotherapy = c( Yes = "red", No = "green", Unknown = "gray"),
  Significance_Sadanandam = c(Sig = "green", Not_sig = "red"),
  Subtype_Sadanandam = c(Norm = "darkgreen", Insulinoma = "Blue", MLP = "Orange", Intermediate = "brown", Unknown = "White")
)
meta_data$Subtype_Sadanandam[meta_data$Grading == ""] = "Norm"
meta_data$Grading[meta_data$Grading == ""] = "G0"
meta_data$NEC_NET[meta_data$NEC_NET == ""] = "Unknown"
meta_data$Subtype_Sadanandam[meta_data$Subtype_Sadanandam == ""] = "Unknown"
meta_data$Dec_dist = rep("",nrow(meta_data))
meta_data[colnames(t(res_coeff[rownames(res_coeff) %in% colnames(expr_raw),])),"Dec_dist"] = as.double(colSums(t(res_coeff[rownames(res_coeff) %in% colnames(expr_raw),])))
meta_data$Dec_dist = as.double(meta_data$Dec_dist)

###

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[13,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw))
colnames(expr) = colnames(expr_raw)
rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

## Figure 1

# Plot 1

pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[c("Study","Dec_dist","Subtype_Sadanandam","Grading","NEC_NET")],
  #annotation_col = meta_data[c("Deco_type")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  treeheight_col = 0,
  legend = F,
  fontsize_col = 7,
  clustering_method = "ward.D2"
)

## Figure 2

# Plot 1

#meta_data$Deco_type = factor( meta_data$Deco_type, levels = c("Alpha","Beta","Gamma","Delta", "Ductal", "Acinar") )

p = ggbiplot::ggbiplot(
  pcr,
  obs.scale = .75,
  groups = as.character(meta_data$Deco_group),
  #shape = as.character(meta_data$Grading),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
  ,labels = meta_data$Name
)
p

### vis

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[13,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]


# Plot 2

p = ggbiplot::ggbiplot(
  pcr,
  obs.scale = .75,
  groups = as.character(meta_data$Study),
  color = meta_data$Study,
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  labels.size =  3,
  labels = meta_data$Name
)
p = p + guides(color=guide_legend(title="Study"))
#p = p + geom_point( aes( as.factor( meta_data$NEC_NET ) ) )
p = p + scale_color_manual( values = c("Darkgreen","Brown") )

#png("~/MAPTor_NET/Results/Dif_exp/Study.PCA.png", width = 1024, height = 768, units = "px", pointsize = 20)
p
#dev.off()

### Estimate contamination plots

estimate_t = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.Estimate.out",sep ="\t", stringsAsFactors = F, header = T)

estimate::estimateScore(
  input.ds = "~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.estimate.gct",
  output.ds = "~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.estimate.out")

colnames(estimate_t) = str_replace( colnames(estimate_t), pattern = "^X", "")
estimate_t = estimate_t[,c(-2)]

order_index = order(estimate_t[4,2:ncol(estimate_t)], decreasing = T)
estimate_t = estimate_t[,c(1,order_index)]
rownames(estimate_t) = estimate_t[,1]
estimate_t = estimate_t[,-1]
rownames = rownames(estimate_t)
colnames = colnames(estimate_t)
estimate_t = matrix(as.double(as.character(unlist(estimate_t))), ncol = ncol(estimate_t), nrow = nrow(estimate_t))

pheatmap::pheatmap(
  estimate_t
)

## MEN1 transcripts

library(ggplot2)

vis_mat = meta_data[,c("Dec_dist","Neurog3")]
vis_mat
colnames(vis_mat) = c("Name","MEN1_exp","MEN1_mt_AF")
vis_mat = reshape2::melt(vis_mat, id = c("Name","MEN1_exp","MEN1_mt_AF"))
vis_mat$Name = factor(vis_mat$Name, levels = vis_mat$Name[order(vis_mat$MEN1_exp)] )

vis_mat$MEN1_mt_AF[vis_mat$MEN1_mt_AF < .5] = 0
vis_mat$MEN1_mt_AF[ (.5 <= vis_mat$MEN1_mt_AF) & ( vis_mat$MEN1_mt_AF < .8)] = .5
vis_mat$MEN1_mt_AF[vis_mat$MEN1_mt_AF >= .8] = 1.0

# Basic barplot
label_vec = meta_data$NEC.AF8.NET[str_detect(meta_data$Histology, pattern = "Pancreatic")]
label_vec = label_vec[order(meta_data$MEN1[str_detect(meta_data$Histology, pattern = "Pancreatic")])]
label_vec[label_vec == "NET"] = "T"
label_vec[label_vec != "T"] = "C"
col_vec = label_vec
col_vec[col_vec == "T"] = "darkgreen"
col_vec[col_vec != "darkgreen"] = "brown"

p = ggplot( data = vis_mat)
p = p + geom_bar(aes( x = Name, y = MEN1_exp, fill = MEN1_mt_AF ),stat="identity", colour="black")
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + scale_fill_gradientn(colours = c("white","yellow","red"), breaks = c(0.0,.5,1.0))
p = p + annotate("text", x=1:42,y = 5.7,parse=TRUE, label = label_vec, color = col_vec, size = 4.5 )
p = p + xlab("") + ylab("MEN1 expression in log TPM") + theme(legend.position = "top")
p

# linear correlation MEN1

vis_mat_2 = subset(vis_mat, MEN1_mt_AF > 0)
lm.model <- lm( vis_mat_2$MEN1_exp ~ vis_mat_2$MEN1_mt_AF) # Fit linear model
summary(lm.model)
cor(vis_mat_2$MEN1_exp , vis_mat_2$MEN1_mt_AF)

# Extract fitted coefficients from model object
b0 <- lm.model$coefficients[1]
b1 <- lm.model$coefficients[2]

g_bench = ggplot(
  data = vis_mat_2,
  aes( y =  MEN1_mt_AF,x = MEN1_exp)
)
g_bench = g_bench + geom_point( )
g_bench = g_bench + geom_smooth(method = "lm")

g_bench + xlab("MEN1 expression")+ ylab("MEN1 mutation AF")

### survival plots ###

meta_data$Diff_Type[meta_data$Diff_Type %in% c("alpha","beta","gamma","delta")]  = "Differentiated"
data = meta_data[,c("Diff_Type","OS_Tissue","Zensur")]
data = data[data$Diff_Type != "not_sig",]
data$OS_Tissue = as.double(str_replace_all(data$OS_Tissue, ",","."))
data = data[!is.na(data[,2]),]
data = data[!is.na(data[,3]),]
colnames(data) = c("Deco_type","OS_Tissue","Status")

fit = survival::survfit( survival::Surv( as.double(data$OS_Tissue), data$Status ) ~ data$Deco_type)

survminer::ggsurvplot(fit, data = data, risk.table = F, pval = T, censor.size = 10)
# Visualize with survminer

cors = apply( expr_raw, MARGIN = 1, FUN = function(vec){return(cor(meta_data$Dec_dist, vec))})
h = (sort(cors))[1:10]
