library("ggplot2")


### Figure 4 plot A box plots

raw_data = read.table("~/Dropbox/testproject/Results/Predicting_grading.tsv",sep ="\t", as.is = T, header = T)
vis_mat = reshape2::melt(raw_data)
colnames(vis_mat) = c("Method","Measurement","Value")
vis_mat$Value = round(vis_mat$Value * 100,0)
#vis_mat = vis_mat %>% filter(Measurement %in% c("Accuracy","Sensitivity","F1","PPV"))

grading_plot = ggplot(vis_mat, aes( x = Measurement, y = Value, fill = Method) )
grading_plot = grading_plot + geom_boxplot()
grading_plot = grading_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
grading_plot = grading_plot + scale_fill_manual(values = c("blue","red"))
grading_plot = grading_plot + ylab("Percentages") + theme(legend.position = "top") + xlab("Performance characteristic")
grading_plot

### alternative Figure 4 Plot A bar plots

vis_mat_expression = vis_mat %>% filter(Method == "Expression")
vis_mat_deconvolution = vis_mat %>% filter(Method == "Deconvolution")

vis_mat_expression_mean = aggregate(vis_mat_expression$Value, by = list(vis_mat_expression$Measurement), FUN = mean)
vis_mat_deconvolution_mean = aggregate(vis_mat_deconvolution$Value, by = list(vis_mat_deconvolution$Measurement), FUN = mean)
colnames(vis_mat_deconvolution_mean) = colnames(vis_mat_expression_mean) = c("Measurement","Mean")

vis_mat_mean = rbind(vis_mat_expression_mean,vis_mat_deconvolution_mean)
vis_mat_mean$Method = c(rep("Expression",nrow(vis_mat_expression_mean)),rep("Deconvolution",nrow(vis_mat_deconvolution_mean)))

vis_mat_expression_sd = aggregate(vis_mat_expression$Value, by = list(vis_mat_expression$Measurement), FUN = sd)
vis_mat_deconvolution_sd = aggregate(vis_mat_deconvolution$Value, by = list(vis_mat_deconvolution$Measurement), FUN = sd)
colnames(vis_mat_expression_sd) = colnames(vis_mat_deconvolution_sd) = c("Measurement","SD")

vis_mat_sd = rbind(vis_mat_expression_sd,vis_mat_deconvolution_sd)
vis_mat_sd$Method = c(rep("Expression",nrow(vis_mat_expression_sd)),rep("Deconvolution",nrow(vis_mat_deconvolution_sd)))

grading_plot = ggplot( vis_mat_mean, aes( x = Measurement, y = Mean, fill = Method) ) + geom_bar(stat="identity", position=position_dodge())
grading_plot = grading_plot + geom_errorbar(aes(ymin = Mean - vis_mat_sd$SD, ymax = Mean + vis_mat_sd$SD),  position = "dodge")
grading_plot = grading_plot + scale_fill_manual(values = c("blue","red")) + ylab("Mean %")

#svg(filename = "~/Dropbox/Figures/Figure_4_Plot_A.svg", width = 10, height = 10)
grading_plot
dev.off()
