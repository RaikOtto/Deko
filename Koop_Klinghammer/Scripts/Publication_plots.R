library("stringr")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))


expr_raw = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.TPM.filtered.ENSEMBL.Voom.TMM.normalized.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

### Prep

library("grid")

meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info$Histology = str_replace_all( meta_info$Histology, pattern = "-","_")
meta_data = meta_info[ match(colnames(expr_raw), meta_info$Name), ]
rownames(meta_data) = meta_data$Name
df_map = match(colnames(expr_raw), meta_data$Name)
meta_data$MEN1_exp = as.double(expr_raw[ which( rownames(expr_raw) == "MEN1"),df_map])
meta_data$MYCN_exp = as.double(expr_raw[ which( rownames(expr_raw) == "MYCN"),df_map])
meta_data$Location = str_replace_all(meta_data$Location, pattern = "-", "_")
meta_data$MKI67 = as.double(expr_raw[ which( rownames(expr_raw) == "MKI67"),df_map])
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
  #Location = c(Primary = "purple",Liver_Met = "#D55E00",Control = "White",Lymph_node_Met = "Yellow",Peritoneum_Met = "Black",Spleen_Met = "Cyan",Lung_Met = "Blue"),
  Location = c(Primary = "white", Metastasis = "black"),
  NEC_NET = c(NEC= "red", NET = "blue",Control = "white"),
  Study = c(Groetzinger = "brown", Scarpa = "darkgreen"),
  MKI67 = c(high = "White", medium = "gray", low = "black"),
  Purity = c(high = "White", medium = "gray", low = "Blue"),
  Correlation = c(high = "Red", medium = "Yellow", low = "Green"),
  MEN1 = c("0" = "White", "1" = "White", "-1" = "#838383"),
  Cellularity = c(high = "White", medium = "gray", low = "Black"),
  Marker_Genes = c(high = "White", medium = "Yellow", low = "Red"),
  Functionality = c( Unknown = "White",Functional = "green", Non_Functional="red"),
  Grading = c( G1 = "Green",G2 = "Yellow", G3 = "Red"),
  Included = c(Yes = "green", No = "red"),
  Chemotherapy = c( Yes = "red", No = "green", Unknown = "gray"),
  Significance_Sadanandam = c(Sig = "green", Not_sig = "red"),
  Subtype_Sadanandam = c(Norm = "darkgreen", Insulinoma = "Blue", MLP = "Orange", Intermediate = "brown")
)


###

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
#sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[13,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[13,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[sad_genes != ""]
new_genes = new_genes[new_genes!=""]

sad_genes[(sad_genes %in% new_genes)]
new_genes[(new_genes %in% sad_genes)]

expr = expr_raw[ str_to_upper(rownames(expr_raw)) %in% str_to_upper(sad_genes),]

#genes_of_interest_hgnc_t[14,3:ncol(genes_of_interest_hgnc_t)] = c(sad_genes, rep("", 613-510))
#write.table(genes_of_interest_hgnc_t, "~/MAPTor_NET/BAMs/Kallisto_three_groups/Stem_signatures.gmt",sep ="\t", quote = F , col.names = F, row.names = F)

cor_mat = cor(expr);pcr = prcomp(t(cor_mat))
cor_mat_pan = cor(expr[,colnames(expr) %in% meta_data$Name[meta_data$Histology == "Pancreatic_NEN"]])

## Figure 1

# Plot 1

meta_data$Location[!str_detect(meta_data$Location,pattern = "Primary")] = "Metastasis"

pheatmap::pheatmap(
  cor_mat_pan,
  annotation_col = meta_data[c("Grading", "Location","Study")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  treeheight_col = 0,
  legend = F,
  fontsize_col = 7,
  clustering_method = "average"
)

## Figure 1

# Plot 2

pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[c("Grading","Location","Histology","Study")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  treeheight_col = 0,
  legend = F,
  fontsize_col = 7
)

## Figure 2

# Plot 1

meta_data$Study = factor( meta_data$Study, levels = c("Scarpa","Groetzinger"))
p = ggbiplot::ggbiplot(
    pcr,
    obs.scale = .75,
    groups = meta_data$Study,
    ellipse = TRUE,
    circle = TRUE,
    var.axes = F#,labels = meta_data$Name
)
MKI67 = as.double( meta_data$MKI67)**1
Grading = as.character(meta_data$Grading)
p = p + geom_point( aes(colour= meta_data$Study, size = MKI67, shape = Grading ) )
p = p + scale_color_manual( values = c("Darkgreen","Brown") )
p = p + guides( color=guide_legend(title="Study", size=guide_legend(title="MKI67"), shape = guide_legend(title="Grading")))

p

## Figure 2 Plot 2
