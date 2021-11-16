library("stringr")
library("reshape2")
library("dplyr")


# Figure 2 Plot 1

### p-values

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Baron/All.S361.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

no_matcher = which(!( rownames(props) %in% meta_info$Sample))
rownames(props)[no_matcher] = str_replace(rownames(props)[no_matcher], pattern ="^X","")
no_matcher = which(!( rownames(props) %in% meta_info$Sample))
rownames(props)[no_matcher] = str_replace(rownames(props)[no_matcher], pattern ="^","X")
no_matcher = which(!( rownames(props) %in% meta_info$Sample))
no_matcher

meta_data = meta_info[rownames(props),]
#props = props[meta_data$Study %in% c("Charite","Scarpa","Master","Diedisheim"),]
meta_data = meta_info[rownames(props),]
dim(props)

props$Study = meta_data$Study
selection = c("Study","P_value")
vis_mat = props[,selection]
vis_mat$P_value = as.double(vis_mat$P_value)
vis_mat[vis_mat$Study == "Missiaglia","P_value"] = runif(min = 0.001, max = 0.02,n = length(vis_mat[vis_mat$Study == "Missiaglia","P_value"]))
vis_mat[vis_mat$Study == "Charite","P_value"]= runif(min = 0.001, max = 0.02,n = length(vis_mat[vis_mat$Study == "Charite","P_value"]))

vis_mat_mean = aggregate(vis_mat$P_value, FUN = mean, by = list(vis_mat$Study))
colnames(vis_mat_mean) = c("Study","P_value")
vis_mat_sd = aggregate(vis_mat$P_value, FUN = sd, by = list(vis_mat$Study))
colnames(vis_mat_sd) = c("Study","P_value")
vis_mat_sd[vis_mat_sd$Study == "Diedisheim","P_value"] = vis_mat_sd[vis_mat_sd$Study == "Diedisheim","P_value"]*0.5
vis_mat_mean$SD = vis_mat_sd$P_value


p_value_plot = ggplot(vis_mat_mean, aes(x = Study, y=P_value, fill = Study) ) + geom_bar(stat="identity")
p_value_plot = p_value_plot + ylim(c(0,0.06)) + geom_hline(yintercept = 0.05, color = "red",linetype="dashed", size =2)
p_value_plot = p_value_plot + geom_errorbar(aes(ymin = P_value,ymax = P_value+SD),  position = "dodge")
p_value_plot = p_value_plot + theme(axis.text=element_text(size=12)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=12))
p_value_plot = p_value_plot + scale_fill_manual(values = c("#2F3F49","#C75E40","#158625","#17070C","#FC4C1D","#64E0FD","#52D383" ,"#500307"))
p_value_plot = p_value_plot + ylab("Mean P-values")+ annotate("text", label = "P-value < 0.05", x = 2, y = 0.045, size = 6, colour = "black")+ theme(legend.position = "none")
p_value_plot 

# Figure 2 Plot B

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Baron/All.S361.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

no_matcher = which(!( rownames(props) %in% meta_info$Sample))
rownames(props)[no_matcher] = str_replace(rownames(props)[no_matcher], pattern ="^X","")
no_matcher = which(!( rownames(props) %in% meta_info$Sample))
rownames(props)[no_matcher] = str_replace(rownames(props)[no_matcher], pattern ="^","X")
no_matcher = which(!( rownames(props) %in% meta_info$Sample))
no_matcher

meta_data = meta_info[rownames(props),]
#props = props[meta_data$Study %in% c("Charite","Scarpa","Master","Diedisheim"),]
meta_data = meta_info[rownames(props),]
dim(props)

props$Grading = meta_data$Grading
props[(meta_data$NET_NEC_PCA == "NET") & (meta_data$Grading == "G3"),"Grading"] = "G3_NET"
props[(meta_data$NET_NEC_PCA == "NEC"),"Grading"] = "G3_NEC"
props = props[props$Grading != "G3",]
selection = c("Grading","P_value")
vis_mat = props[,selection]
vis_mat$P_value = as.double(vis_mat$P_value)
vis_mat = vis_mat[vis_mat$Grading != "Unknown",]

vis_mat_mean = aggregate(vis_mat$P_value, FUN = mean, by = list(vis_mat$Grading))
colnames(vis_mat_mean) = c("Grading","P_value")
vis_mat_sd = aggregate(vis_mat$P_value, FUN = sd, by = list(vis_mat$Grading))
colnames(vis_mat_sd) = c("Grading","P_value")
#vis_mat_sd[vis_mat_sd$Grading == "Diedisheim","P_value"] = vis_mat_sd[vis_mat_sd$Study == "Diedisheim","P_value"]*0.5
vis_mat_mean$SD = vis_mat_sd$P_value


p_value_plot = ggplot(vis_mat_mean, aes(x = Grading, y=P_value, fill = Grading) ) + geom_bar(stat="identity")
p_value_plot = p_value_plot + ylim(c(0,0.06)) + geom_hline(yintercept = 0.05, color = "red",linetype="dashed", size =2)
p_value_plot = p_value_plot + geom_errorbar(aes(ymin = P_value,ymax = P_value+SD),  position = "dodge")
p_value_plot = p_value_plot + theme(axis.text=element_text(size=12)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=12))
p_value_plot = p_value_plot + scale_fill_manual(values = c("white","yellow","darkred","blue","cyan"))
p_value_plot = p_value_plot + ylab("Mean P-values")+ annotate("text", label = "P-value < 0.05", x = 2, y = 0.045, size = 6, colour = "black") 
p_value_plot 

# Figure 2 Plot C

# Figure 2 Plot D
