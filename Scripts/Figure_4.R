library("ggplot2")
library("stringr")
library("dplyr")
library("reshape2")

ml_predictions = read.table("~/Dropbox/testproject/Results/PanNEN_only/Visualization/ML_predictions_means.tsv", sep ="\t", header = TRUE)

### Figure 4 Plot A bar plots

grading_binary_mat = ml_predictions %>% filter(Experiment == "Grading_binary")
grading_binary_mat = reshape2::melt(grading_binary_mat)
grading_binary_mat = grading_binary_mat[,colnames(grading_binary_mat) != "Experiment"]
colnames(grading_binary_mat) = c("Method","Measurement","Mean")
grading_binary_mat$Mean = round(grading_binary_mat$Mean*100,0)
grading_binary_mat_mean = grading_binary_mat[str_detect(grading_binary_mat$Measurement, "_SD",negate = TRUE),]
grading_binary_mat_SD = grading_binary_mat[str_detect(grading_binary_mat$Measurement, "_SD",negate = FALSE),]
#grading_binary_mat_mean$Method = factor(grading_binary_mat_mean$Method, levels = c(""))

grading_plot = ggplot( grading_binary_mat_mean, aes( x = Measurement, y = Mean, fill = Method) ) + geom_bar(stat="identity", position=position_dodge())
grading_plot = grading_plot + geom_errorbar(aes(ymin = Mean - grading_binary_mat_SD$Mean, ymax = Mean + grading_binary_mat_SD$Mean),  position = "dodge")
grading_plot = grading_plot + scale_fill_manual(values = c("red","blue")) + ylab("Mean %") + theme(legend.position="top")
grading_plot = grading_plot + ylim(0,100)

#svg(filename = "~/Dropbox/Figures/Figure_4_Plot_A.svg", width = 10, height = 10)
grading_plot
dev.off()

###

### Figure 4 Plot A bar plots

grading_tertiary_mat = ml_predictions %>% filter(Experiment == "Grading_tertiary")
grading_tertiary_mat = reshape2::melt(grading_tertiary_mat)
grading_tertiary_mat = grading_tertiary_mat[,colnames(grading_tertiary_mat) != "Experiment"]
colnames(grading_tertiary_mat) = c("Method","Measurement","Mean")
grading_tertiary_mat$Mean = round(grading_tertiary_mat$Mean*100,2)
grading_tertiary_mat_mean = grading_tertiary_mat[str_detect(grading_tertiary_mat$Measurement, "_SD",negate = TRUE),]
grading_tertiary_mat_SD = grading_tertiary_mat[str_detect(grading_tertiary_mat$Measurement, "_SD",negate = FALSE),]
#grading_binary_mat_mean$Method = factor(grading_binary_mat_mean$Method, levels = c(""))

grading_plot = ggplot( grading_tertiary_mat_mean, aes( x = Measurement, y = Mean, fill = Method) ) + geom_bar(stat="identity", position=position_dodge())
grading_plot = grading_plot + geom_errorbar(aes(ymin = Mean - grading_tertiary_mat_SD$Mean, ymax = Mean + grading_tertiary_mat_SD$Mean),  position = "dodge")
grading_plot = grading_plot + scale_fill_manual(values = c("red","blue")) + ylab("Mean %") + theme(legend.position="top")
grading_plot = grading_plot + ylim(0,100)

#svg(filename = "~/Dropbox/Figures/Figure_4_Plot_B.svg", width = 10, height = 10)
grading_plot
dev.off()

### Figure 4 Plot A bar plots

nec_net_mat = ml_predictions %>% filter(Experiment == "NET_NEC")
nec_net_mat = reshape2::melt(nec_net_mat)
nec_net_mat = nec_net_mat[,colnames(nec_net_mat) != "Experiment"]
colnames(nec_net_mat) = c("Method","Measurement","Mean")
nec_net_mat$Mean = round(nec_net_mat$Mean*100,2)
net_nec_mat_mean = nec_net_mat[str_detect(nec_net_mat$Measurement, "_SD",negate = TRUE),]
net_nec_mat_SD = nec_net_mat[str_detect(nec_net_mat$Measurement, "_SD",negate = FALSE),]

net_nec_plot = ggplot( net_nec_mat_mean, aes( x = Measurement, y = Mean, fill = Method) ) + geom_bar(stat="identity", position=position_dodge())
net_nec_plot = net_nec_plot + geom_errorbar(aes(ymin = Mean - net_nec_mat_SD$Mean, ymax = Mean + net_nec_mat_SD$Mean),  position = "dodge")
net_nec_plot = net_nec_plot + scale_fill_manual(values = c("red","blue")) + ylab("Mean %") + theme(legend.position="top")
net_nec_plot = net_nec_plot + ylim(0,100)

#svg(filename = "~/Dropbox/Figures/Figure_4_Plot_C.svg", width = 10, height = 10)
net_nec_plot
dev.off()


# mean calculator

#vis_mat_expression = vis_mat %>% filter(Method == "Expression")
#vis_mat_deconvolution = vis_mat %>% filter(Method == "Deconvolution")

#vis_mat_expression_mean = aggregate(vis_mat_expression$Value, by = list(vis_mat_expression$Measurement), FUN = mean)
#vis_mat_deconvolution_mean = aggregate(vis_mat_deconvolution$Value, by = list(vis_mat_deconvolution$Measurement), FUN = mean)
#colnames(vis_mat_deconvolution_mean) = colnames(vis_mat_expression_mean) = c("Measurement","Mean")

#vis_mat_mean = rbind(vis_mat_expression_mean,vis_mat_deconvolution_mean)
#vis_mat_mean$Method = c(rep("Expression",nrow(vis_mat_expression_mean)),rep("Deconvolution",nrow(vis_mat_deconvolution_mean)))

#vis_mat_expression_sd = aggregate(vis_mat_expression$Value, by = list(vis_mat_expression$Measurement), FUN = sd)
#vis_mat_deconvolution_sd = aggregate(vis_mat_deconvolution$Value, by = list(vis_mat_deconvolution$Measurement), FUN = sd)
#colnames(vis_mat_expression_sd) = colnames(vis_mat_deconvolution_sd) = c("Measurement","SD")

#vis_mat_sd = rbind(vis_mat_expression_sd,vis_mat_deconvolution_sd)
#vis_mat_sd$Method = c(rep("Expression",nrow(vis_mat_expression_sd)),rep("Deconvolution",nrow(vis_mat_deconvolution_sd)))
