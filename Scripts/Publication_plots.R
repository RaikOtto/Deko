library("stringr")
library("grid")

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

## Figure 1

# Plot 1

#meta_data$Location[!str_detect(meta_data$Location,pattern = "Primary")] = "Metastasis"
#meta_data$Grading[meta_data$Grading == ""] = "G0"
pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[c("Location")],
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

## Figure 2

# Plot 1

data_t = read.table("~/Deko/Results/ROC_curves.tsv",sep ="\t", header = T)

# hisc

data_t_hisc = subset(data_t, Predictor == "HISC" )
vis_mat_hisc = aggregate(
    as.integer(subset(data_t_hisc$ROC, data_t_hisc$Type == "In-silico")),
    by = list(as.character(subset(data_t_hisc$Dataset, data_t_hisc$Type == "In-silico"))),
    FUN = mean
)
vis_mat_hisc$SD = aggregate(
    as.integer(subset(data_t_hisc$ROC, data_t_hisc$Type == "In-silico")),
    by = list(as.character(subset(data_t_hisc$Dataset, data_t_hisc$Type == "In-silico"))),
    FUN = sd
)
vis_mat_hisc$Type = rep("HISC", nrow(vis_mat_hisc))
vis_mat_hisc = as.data.frame(as.matrix(vis_mat_hisc))
vis_mat_hisc = vis_mat_hisc[,-which(colnames(vis_mat_hisc)  == "SD.Group.1")]
colnames(vis_mat_hisc) = c("Dataset","ROC","SD","Type")

# ductal

data_t_ductal = subset(data_t, Predictor == "Ductal" )
vis_mat_ductal = aggregate(
    as.integer(subset(data_t_ductal$ROC, data_t_ductal$Type == "In-silico")),
    by = list(as.character(subset(data_t_ductal$Dataset, data_t_ductal$Type == "In-silico"))),
    FUN = mean
)
vis_mat_ductal$SD = aggregate(
    as.integer(subset(data_t_ductal$ROC, data_t_ductal$Type == "In-silico")),
    by = list(as.character(subset(data_t_ductal$Dataset, data_t_ductal$Type == "In-silico"))),
    FUN = sd
)
vis_mat_ductal$Type = rep("Ductal", nrow(vis_mat_ductal))
vis_mat_ductal = as.data.frame(as.matrix(vis_mat_ductal))
vis_mat_ductal = vis_mat_ductal[,-which(colnames(vis_mat_ductal)  == "SD.Group.1")]
colnames(vis_mat_ductal) = c("Dataset","ROC","SD","Type")

# MKI-67

data_t_mki67 = subset(data_t, Predictor == "MKI-67" )
vis_mat_mki67 = aggregate(
    subset(data_t_mki67$ROC, data_t_mki67$Type == "MKI-67"),
    by = list(subset(data_t_mki67$Dataset, data_t_mki67$Type == "MKI-67")),
    FUN = mean
)
vis_mat_mki67$SD = rep(0, nrow(vis_mat_mki67))
vis_mat_mki67$Type = rep("MKI_67", nrow(vis_mat_mki67))
colnames(vis_mat_mki67) = c("Dataset","ROC","SD","Type")

# merge

vis_mat = rbind(as.matrix(vis_mat_ductal),as.matrix(vis_mat_hisc),as.matrix(vis_mat_mki67))
vis_mat = as.data.frame(vis_mat)
vis_mat$ROC = as.double(as.character(vis_mat$ROC))
vis_mat$SD = as.double(as.character(vis_mat$SD))

### plots

p = ggplot( data = vis_mat,aes( x = Dataset, y = as.integer(ROC), min = ROC-SD, max = ROC+SD, fill =Type) )
p = p + geom_bar(stat="identity", position=position_dodge())
#p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#p = p + scale_fill_gradientn(colours = c("white","yellow","red"), breaks = c(0.0,.5,1.0))
#p = p + annotate("text", x=1:57,y = 5.5,parse=TRUE, label = label_vec, color = col_vec, size = 4.5 )
#p = p + xlab("") + ylab("MEN1 expression in log TPM and MEN1 mutation allele frequency") + theme(legend.position = "top")
p + geom_errorbar(aes(),  position = "dodge")


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

### MEN1 exp t-test

### chr 5

vis_mat = deconvolution_results[,c("MKI67","ductal","Grading")]
vis_mat = reshape2::melt(vis_mat)
colnames(vis_mat) = c("Grading","Type","Value")

men1_plot = ggplot( vis_mat, aes ( x = Grading,  y = Value, fill = Type))
men1_plot = men1_plot + geom_boxplot( aes())
men1_plot = men1_plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
men1_plot

### survival plots ###

surv_cor = apply(expr_raw, MARGIN = 1, FUN = function(vec){return(cor(vec, as.double(meta_data$OS_Tissue)))})

subtype = as.double(expr_raw[ which( rownames(expr_raw) == "POU5F1"),df_map])
cut_off = quantile(sort(subtype), probs = seq(0,1,.1) )[9]
subtype[subtype > as.double(cut_off) ] = "Above"
subtype[subtype != "Above" ] = "Below"

##fit = survival::survfit( survival::Surv( meta_data$OS_Tissue) ~ subtype,data = meta_data)
cell_type = meta_data$Deco_type
cell_type[cell_type != "Not_sig"] = "Sig"
#col_vec = 
fit = survival::survfit( survival::Surv( meta_data$OS_Tissue ) ~ cell_type)

survminer::ggsurvplot(fit, data = meta_data, risk.table = T, pval = T)
# Visualize with survminer
