library("ggplot2")

uniquorn_1 = c(8*10^6,7*10^6,6.5*10^6)
uniquorn_2 = c(49, 151, 100, 6.5*10^6,8*10^6,9*10^6,11*10^6,10*10^6)

data_mat_count = data.frame(
  "Dataset" = as.factor(c(
    rep("Uniquorn 1",length(uniquorn_1)),
    rep("Uniquorn 2",length(uniquorn_2))
  )),
  "Counts" = c(
    uniquorn_1,
    uniquorn_2
  )
)


data_mat = data.frame(
  "Dataset" = c(
    rep("Uniquorn 1",1),
    rep("Uniquorn 2",1)
  ),
  "Mean" = c(
    mean(uniquorn_1),
    mean(uniquorn_2)
  ),
  "SD" = c(
    sd(uniquorn_1),
    sd(uniquorn_2)
  )
)

p = ggplot( 
  data = data_mat_count,
  aes( 
    x = Dataset,
    y = Counts,
    fill = Dataset
  )
)
p = p + geom_boxplot()
p = p + xlab("") + ylab("Counts per CCL in millions") + theme(legend.position = "top")
p = p + scale_fill_manual(values = c("blue","red"))
p = p + theme(legend.position="none",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p


raw_dara = read.table("~/Dropbox/performance_deconvolution.csv",sep ="\t", as.is = T, header = T)
colnames(raw_dara) = c("Method","Measurement","Accuracy","AUC","Sensitivity","PPV","F1","Kapp","MCC")
vis_mat = reshape2::melt(raw_dara)
colnames(vis_mat) = c("Method","Measurement","Type","Value")
vis_mat = vis_mat %>% filter(Type %in% c("Accuracy","Sensitivity","PPV"))

error_mat_means =  vis_mat %>% filter(Measurement %in% c("Mean"))
error_mat_SD =  vis_mat %>% filter(Measurement %in% c("SD"))
error_low = (error_mat_means$Value - error_mat_SD$Value)*100
error_high = (error_mat_means$Value + error_mat_SD$Value)*100
error_high[error_high > 100] = 100
  
vis_mat = vis_mat %>% filter(Measurement %in% c("Mean"))
vis_mat$Value = vis_mat$Value*100
#vis_mat[as.character(vis_mat$Type) == "Prec.","Type"] = "PPV"
#

p = ggplot( 
  data = vis_mat,
  aes( 
    x = Type,
    y = Value,
    fill = Method
  )
)
p = p + geom_bar(
  aes( 
    x = Type,
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
