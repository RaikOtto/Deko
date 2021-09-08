library("stringr")
library("ggplot2")
library("dplyr")

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Califano.S165.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

no_match = props$Sample %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",props$Sample[no_match],sep ="")
no_match = props$Sample %in% meta_info$Sample == F
sum(no_match)

dim(props)
meta_data = meta_info[props$Sample,]

#props = props[(meta_data$Histology == "pancreas") | (meta_data$NEC_NET_Color != "Primary") ,]
props  = props[ meta_data$Histology == "pancreas" ,]
dim(props)
meta_data = meta_info[props$Sample,]

meta_data$Location = meta_data$NEC_NET_Color
props = props[ meta_data$Location != "Outlier" ,]
meta_data = meta_info[props$Sample,]
props$Location = meta_data$NEC_NET_Color

vis_mat = reshape2::melt(props[,c("Location","P_value","model")])
colnames(vis_mat) = c("Location","Model","Variable","P_value")
vis_mat$Location[vis_mat$Location == "liver_met"] = "Metastasis"
vis_mat$Location = factor(vis_mat$Location, levels = c("Primary","Metastasis"))
vis_mat$P_value = as.double(vis_mat$P_value)

vis_mat$Model[ vis_mat$Model == "Alpha_Beta_Gamma_Delta_Baron"] = "Endocrine"
vis_mat$Model[ vis_mat$Model == "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron"] = "Endocrine & Exocrine"

p_values = ggplot(
    data = vis_mat,
    aes(
        x = Model,
        y = P_value,
        fill = Location
    )
) + geom_boxplot()
p_values = p_values + scale_fill_manual(values = c("blue", "red","blue", "red")) + ylab("P-value") + xlab("Model")+ theme(legend.position = "top",axis.text=element_text(size=12))
p_values = p_values + geom_hline(yintercept = 0.05, linetype= "dashed", color = "black")
p_values = p_values + annotate("text", x = 2, y = .055, label = "P-value == 0.05",col = "black",size =6)
p_values = p_values + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p_values

