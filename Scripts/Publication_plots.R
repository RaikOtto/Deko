library("stringr")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))


#expr_raw = bam_data#
expr_raw = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv",sep="\t", stringsAsFactors =  F, header = T)
#expr_raw = read.table("~/Deko/Data/TPMs.Not_normalized.Controls_Groetzinger_Scarpa.89S.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

### Prep

library("grid")

meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_data = meta_info[ match(colnames(expr_raw), meta_info$Name), ]
rownames(meta_data) = meta_data$Name
df_map = match(colnames(expr_raw), meta_data$Name)
meta_data$INS = as.double(expr_raw[ which( rownames(expr_raw) == "INS"),df_map])
meta_data$GCG = as.double(expr_raw[ which( rownames(expr_raw) == "GCG"),df_map])
meta_data$PPY = as.double(expr_raw[ which( rownames(expr_raw) == "PPY"),df_map])
meta_data$SST = as.double(expr_raw[ which( rownames(expr_raw) == "SST"),df_map])
meta_data$TTR = as.double(expr_raw[ which( rownames(expr_raw) == "TTR"),df_map])
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
  Deco_type = c(Alpha = "Blue",Beta = "Yellow",Gamma = "Orange",Delta = "Purple",Acinar = "Black",Ductal = "Brown", Not_sig = "Gray", Alpha_five = "Red"),
  Location = c(Primary = "white", Metastasis = "black"),
  NEC_NET = c(NEC= "red", NET = "blue", Unknown = "white"),
  Study = c(Groetzinger = "brown", Scarpa = "darkgreen"),
  MKI67 = c(high = "White", medium = "gray", low = "black"),
  Purity = c(high = "White", medium = "gray", low = "Blue"),
  Correlation = c(high = "Red", medium = "Yellow", low = "Green"),
  ImmuneScore = c(high = "White", medium = "gray", low = "Black"),
  StromalScore = c(high = "White", medium = "gray", low = "Black"),
  TumorPurity = c(high = "White", medium = "gray", low = "Black"),
  Marker_Genes = c(high = "White", medium = "Yellow", low = "Red"),
  Functionality = c( Unknown = "White",Functional = "green", Non_Functional="red"),
  Grading = c( G1 = "Green",G2 = "Yellow", G3 = "Red", "Other" = "green", G0 = "white"),
  Included = c(Yes = "green", No = "red"),
  Chemotherapy = c( Yes = "red", No = "green", Unknown = "gray"),
  Significance_Sadanandam = c(Sig = "green", Not_sig = "red"),
  Subtype_Sadanandam = c(Norm = "darkgreen", Insulinoma = "Blue", MLP = "Orange", Intermediate = "brown", Unknown = "White")
)

meta_data$Subtype_Sadanandam[meta_data$Grading == ""] = "Norm"
meta_data$Grading[meta_data$Grading == ""] = "G0"
meta_data$NEC_NET[meta_data$NEC_NET == ""] = "Unknown"
meta_data$Subtype_Sadanandam[meta_data$Subtype_Sadanandam == ""] = "Unknown"

###

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[15,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

## Figure 1

# Plot 1

pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[c("Deco_type","Subtype_Sadanandam","Grading","NEC_NET")],
  #annotation_col = meta_data[c("Deco_type")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  treeheight_col = 0,
  legend = F,
  fontsize_col = 7
)

## Figure 2

# Plot 1

#meta_data$Deco_type = factor( meta_data$Deco_type, levels = c("Alpha","Beta","Gamma","Delta", "Ductal", "Acinar") )

p = ggbiplot::ggbiplot(
    pcr,
    obs.scale = .75,
    groups = as.character(meta_data$Deco_type),
    #shape = as.character(meta_data$Grading),
    ellipse = TRUE,
    circle = TRUE,
    var.axes = F#,labels = meta_data$Name
)
MKI67 = as.double( meta_data$MKI67)**1
Grading = as.character(meta_data$Grading)
p = p + scale_color_manual( values = c("Blue","Red","Yellow","Purple", "Orange","Black") )
p

## Figure 2 Plot 2

nec_net_col_vec = meta_data$NEC_NET;nec_net_col_vec[nec_net_col_vec == "NEC"] = "Brown";nec_net_col_vec[nec_net_col_vec != "Brown"] = "Darkgreen"

MEN1_status_col = meta_data$MEN1
MEN1_status_col[MEN1_status_col <= 0] = "MEN1_wt"
MEN1_status_col[MEN1_status_col != "MEN1_wt"] = "MEN1_mt"
MEN1_status_col[ ( meta_data$NEC_NET  == "NET") & (MEN1_status_col == "MEN1_wt") ] = "MEN1_wt_net"
MEN1_status_col[ ( meta_data$NEC_NET  == "NEC") & (MEN1_status_col == "MEN1_wt") ] = "MEN1_wt_nec"

names(meta_data$Grading)  <- "Experimental Condition"
size_vec = meta_data$Grading
size_vec[ ( meta_data$Grading == "G3")  ]= 5;size_vec[size_vec != "5"]= 1

p = ggbiplot::ggbiplot(
  pcr,
  obs.scale =.1,
  #var.scale = 1, 
  labels.size = 4,
  alpha = 1,
  groups = as.character(meta_data$NEC_NET),
  #groups = as.character(df$Study),
  #shape = as.character(df$Grading),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
)
p = p + geom_point( aes( size = as.integer(size_vec)), shape = c( 15 ),  color = "black", stroke = .9)
p = p + geom_point( aes( shape = Grading, color = meta_data$Cell_type ), size = 4)
p = p + scale_color_manual( values = c("Blue","Yellow","Purple","Black","Orange","Brown","Darkgreen","Gray") )
#p = p + guides(color=guide_legend(title="NEC_NET"), position = "top") + guides(shape=guide_legend(title="Grading"), position = "top", size = F)

#png("~/MAPTor_NET/Results/Dif_exp/Grading.PCA.png", width = 1024, height = 768, units = "px", pointsize = 20)
    p
#dev.off()

### Figure 2 Plot 3

size_vec = meta_data$Location
size_vec[ str_detect( as.character(size_vec), pattern = "_Met" ) ]= "1";
size_vec[size_vec != "1"]= "3"

p = ggbiplot::ggbiplot(
  pcr,
  obs.scale = .75,
  groups = as.character(meta_data$NEC_NET),
  color = as.character(meta_data$Histology),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
)
p = p + geom_point( aes(colour= meta_data$Histology, size = size_vec ) )
p = p + guides(color=guide_legend(title="Histology"), size = guide_legend(title="Type"))
p = p + scale_color_manual( values = c("orange","pink","purple","RED","BLUE","black","yellow") )
p = p + scale_size_manual(labels = c("Metastasis", "Primary"), values=as.integer(size_vec)*2)

#png("~/MAPTor_NET/Results/Dif_exp/Histology.PCA.png", width = 1024, height = 768, units = "px", pointsize = 20)
    p
#dev.off()
    
        
### Figure 2 Plot 4

size_vec = meta_data$Significance_Sadanandam
size_vec[ as.character(size_vec) == "Not_sig" ]= "0";size_vec[size_vec != "0"]= "1"

p = ggbiplot::ggbiplot(
  pcr,
  obs.scale =.75,
  var.scale = 2, 
  labels.size = 8,
  alpha = 1,
  groups = as.character(meta_data$Subtype_Sad),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
)
p = p + geom_point( aes( size = as.integer(size_vec), color = as.factor(meta_data$Subtype_Sad) )) +guides(size=guide_legend(title="P_value")) +guides(color=guide_legend(title="Subtype_Sadanandam"))
p = p + scale_color_manual( values = c("Blue","Brown","Orange","Darkgreen") )

#png("~/MAPTor_NET/Results/Dif_exp/Sadanandam_signature.57.png", width = 1024, height= 768)
    p
#dev.off()

### Figure 3 Plot 1

#png("~/MAPTor_NET/Results/Dif_exp/Heatmap_57.integrated.png",width=1014, height = 768)
    pheatmap::pheatmap(
      cor_mat,
      annotation_col = meta_data[c("MEN1","Marker_Genes","MKI67","NEC_NET","Grading","Location","Histology","Study")],
      annotation_colors = aka3,
      annotation_legend = T,
      treeheight_col = 0,
      show_colnames = T,
      show_rownames = F,
      gaps_row = 20
    )
#dev.off()

### Figure 3 Plot 2

pheatmap::pheatmap(
  cor_mat,
  #annotation_col = df[c("MEN1","Marker_genes","MKI67","NEC_NET","Grading","Location","Histology","Study")],
  annotation_colors = aka3,
  annotation_legend = T,
  treeheight_col = 0,
  show_colnames = T,
  show_rownames = F
)

############ Supplementary

# Figure 1


### vis

cor_mat_all = cor(expr)

meta_data$Chemotherapy[ meta_data$Chemotherapy == "" ] = "Unknown"
meta_data$Chemotherapy[ meta_data$Chemotherapy == "0" ] = "No"
meta_data$Chemotherapy[ meta_data$Chemotherapy == "1" ] = "Yes"
meta_data$Chemotherapy[ meta_data$Chemotherapy == "Missing" ] = "Unknown"

pheatmap::pheatmap(
  cor_mat_all,
  annotation_col = meta_data[c("Included","Chemotherapy","Location","Histology","Study")],
  #annotation_colors = aka3,
  show_colnames = F,
  treeheight_col = 0,
  treeheight_row = 0
)

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

vis_mat = meta_data[str_detect(meta_data$Histology, pattern = "Pancreatic"),c("Name","MEN1","MEN1.AF8.Mut.AF8.AF")]
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
