
require(ggplot2)
require(cowplot)
require(plotluck)
require(dplyr)
library(ggrepel)



#script to draw mosaic plot of different functional types (fig. 6)


#these values were taken from runing the Functionalization.R script 
Defunc<-c(7/85,16/99,2/98)
Refunc<-c(64/85,53/99,75/98)
Isofunc<-c(4/85,4/99,17/98)
Refunc_HardIso<-c(1/85,10/99,20/98)
Refunc_SoftIso<-c(36/85,17/99,52/98)
Func_Reaqur<-c(27/85,26/99,3/98)

#blank vector so I can write labels in its spot
Na_vec<-c("NA","NA","NA")

#put everything into a nice data frame
df<-data.frame(Exp=rep(c("PDB:\n2EKE","PDB:\n2EKE\n(Duplicate A)", "PDB:\n4GVB"),5), data=c(Defunc,Isofunc,Refunc_HardIso,Func_Reaqur,Refunc_SoftIso), label=c(rep("Defunctionalization",3),rep("Isofunctionalization",3),rep("Refunctionalization Hard Iso",3),rep("Refunctionalization Reacquisition",3),rep("Refunctionalization Soft Iso",3)))


df <- df %>%
  group_by(Exp, label) %>%
  ungroup()


#pick nice colors
safe<-c("#1b9e77","#d95f02","#7570b3","#9793c6","#bab7d9")


#make lighter versions
lighten <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  col
}


#make colors lighter
safe<-lapply(safe,lighten)


#make plot of A
p<-ggplot(df, aes(x = Exp, y = data, fill = label)) + geom_bar(stat = "identity", position = "fill", colour = "black") + 
  scale_fill_manual(values=safe)+
  ylab("Fraction of Simulations")+
  xlab("Protein Structure")+
  scale_x_discrete(limits=c("PDB:\n2EKE","PDB:\n2EKE\n(Duplicate A)", "PDB:\n4GVB","",""))+ 
  theme(axis.title.x = element_text(hjust=.2),axis.line.x= element_blank(),axis.ticks.x= element_line(color=c("black","black","black","white","white")))+
  scale_y_continuous(expand = c(0, 0),limit=c(0,1.02))+
  geom_label_repel(data = subset(df, Exp == "PDB:\n4GVB"), aes(label=c("Defunctionalization","Isofunctionalization","Refunctionalization:\nHard\n Isofunctionalization","Refunctionalization:\nSoft\nIsofunctionalization","Refunctionalization:\nReacquisition"),y=c(.99,.9,.7,.57,.3), x=c(3.45,3.45,3.45,3.45,3.45),fill=factor(label)),segment.color="black",size=3, nudge_x=3,direction=c("x") ) +
  theme(legend.position="none")


#data for protein C
Defunc<-c(7/85,16/99,2/98)
Refunc<-c(64/85,53/99,75/98)
Isofunc<-c(11/85,18/99,16/98)
Refunc_HardIso<-c(1/85,12/99,1/98)
Refunc_SoftIso<-c(38/85,36/99,58/98)
Func_Reaqur<-c(25/85,5/99,16/98)

#again make a vector to write labels in
Na_vec<-c("NA","NA","NA")


#make a nice dataframe
df<-data.frame(Exp=rep(c("PDB:\n2EKE","PDB:\n2EKE\n(Duplicate A)", "PDB:\n4GVB"),5), data=c(Defunc,Isofunc,Refunc_HardIso,Func_Reaqur,Refunc_SoftIso), label=c(rep("Defunctionalization",3),rep("Isofunctionalization",3),rep("Refunctionalization Hard Iso",3),rep("Refunctionalization Reacquisition",3),rep("Refunctionalization Soft Iso",3)))



df <- df %>%
  group_by(Exp, label) %>%
  ungroup()


#plot data for protein C
p2<-ggplot(df, aes(x = Exp, y = data, fill = label)) + geom_bar(stat = "identity", position = "fill", colour = "black") + 
  scale_fill_manual(values=safe)+
  ylab("Fraction of Simulations")+
  xlab("Protein Structure")+
  scale_x_discrete(limits=c("PDB:\n2EKE","PDB:\n2EKE\n(Duplicate A)", "PDB:\n4GVB","",""))+ 
  theme(axis.title.x = element_text(hjust=.2),axis.line.x= element_blank(),axis.ticks.x= element_line(color=c("black","black","black","white","white")))+
  scale_y_continuous(expand = c(0, 0),limit=c(0,1.02))+
  geom_label_repel(data = subset(df, Exp == "PDB:\n4GVB"), aes(label=c("Defunctionalization","Isofunctionalization","Refunctionalization:\nHard\n Isofunctionalization","Refunctionalization:\nSoft\nIsofunctionalization","Refunctionalization:\nReacquisition"),y=c(.99,.9,.8,.72,.3), x=c(3.45,3.45,3.45,3.45,3.45),fill=factor(label)),segment.color="black",size=3, nudge_x=3,direction=c("x") ) +
  theme(legend.position="none")


p2
p3<-plot_grid(p, p2, labels = c('A', 'B'))


ggsave(filename="TypeFunAB_nice.pdf", plot=p3, width=9, height=6, units="in")


