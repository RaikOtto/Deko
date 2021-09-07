library("ggplot2")


raw_data = read.table("~/Dropbox/ML_performance_deko.tsv",sep ="\t", as.is = T, header = T)
vis_mat = reshape2::melt(raw_data)
colnames(vis_mat) = c("Type","Method","Measurement","Value")
vis_mat = vis_mat %>% filter(Measurement %in% c("Sensitivity","F1","PPV"))

error_mat_means =  vis_mat %>% filter(Type %in% c("Mean"))
error_mat_SD =  vis_mat %>% filter(Type %in% c("SD"))
error_low = (error_mat_means$Value - error_mat_SD$Value)*100
error_high = (error_mat_means$Value + error_mat_SD$Value)*100
error_high[error_high > 100] = 100
  
vis_mat = vis_mat %>% filter(Type %in% c("Mean"))
vis_mat$Value = vis_mat$Value*100
vis_mat$Measurement = factor(vis_mat$Measurement, levels = c("Sensitivity","F1","PPV"))

p = ggplot( 
  data = vis_mat,
  aes( 
    x = Measurement,
    y = Value,
    fill = Method
  )
)
p = p + geom_bar(
  aes(
    x = Measurement,
    y = Value,
    fill = Method
  ),
  stat="identity",position =position_dodge())
p = p + geom_errorbar(aes(ymin=error_low, ymax=error_high),
  width=.2,
  position=position_dodge(.9),
  size = 1.2)
p = p + xlab("") + ylab("Performance in %") + theme(legend.position = "top")
p = p + scale_fill_manual(values = c("blue","red"))
p = p + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p
