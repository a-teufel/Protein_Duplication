#script to calcuate the probability of an interface sub occuring during the first 250 subs and 
#graph it for each protein


rm(list=ls())
library(tidyverse)
library(stringr)
library(cowplot)
library(broom)
library(dplyr)
library(grid)
library(ggridges)
library(colorspace)


A_sum_WT<-0
B_sum_WT<-0
C_sum_WT<-0

A_sum_nobind<-0
B_sum_nobind<-0
C_sum_nobind<-0


#site numbers that are in the interface. these numbers come from Spider
A_sites<-c(20,21,23,25,27,29,30,33,67,68,69)
B_sites<-c(88,95,97,98,101,135,136,137,138,139)-77


#detects where two sequences differ
detect_subs<-function(prev_sequences, current_sequences){
  
  subs<-NULL
  
  for(k in 1:length(current_sequences)){
    split<-current_sequences[k]
    split_prev<-prev_sequences[k]
    
    temp<-data.frame(New=strsplit(split,"")[[1]],Old=strsplit(split_prev,"")[[1]],stringsAsFactors=F)
    
    subs<-c(subs,which(temp$New != temp$Old))
  }
  
  return(subs)
}




# Import the first wt sequences
name_WT<-paste("C:/Users/User/Documents/dup_alignments/M_SIM1_4gvb_csv_files_first_500/M_SIM1_4gvb_1.csv",sep="")  
wt_sequences = read_csv( name_WT, col_names = c("simulation", "proteinA", "proteinB1", "proteinB2", "proteinC1","proteinC2","empty"))
proteinB<-paste(wt_sequences$proteinB1,wt_sequences$proteinB2,sep="")
proteinC<-paste(wt_sequences$proteinC1,wt_sequences$proteinC2,sep="")
wt_sequences$proteinB=proteinB
wt_sequences$proteinC=proteinC

WT_prev<- wt_sequences


# Import the first no bind sequences
name_EXP<-paste("C:/Users/User/Documents/dup_alignments/M_SIM2_4gvb_csv_files_first_500/M_SIM2_4gvb_1.csv",sep="")  
nobind_sequences = read_csv( name_EXP, col_names = c("simulation", "proteinA", "proteinB1", "proteinB2", "proteinC1","proteinC2","empty"))
proteinB<-paste(nobind_sequences$proteinB1,nobind_sequences$proteinB2,sep="")
proteinC<-paste(nobind_sequences$proteinC1,nobind_sequences$proteinC2,sep="")
nobind_sequences$proteinB=proteinB
nobind_sequences$proteinC=proteinC


nobind_prev<-nobind_sequences

#set up tables to save the counts of where the differences are
A_table<-NULL
B_table<-NULL
C_table<-NULL

A_table_no_bind<-NULL
B_table_no_bind<-NULL
C_table_no_bind<-NULL

#go over the first 250 
for(i in seq(11,250,10)){
  # Import & tidy wildtype sequencesseq
  name_WT<-paste("C:/Users/User/Documents/dup_alignments/M_SIM1_4gvb_csv_files_first_500/M_SIM1_4gvb_",i,".csv",sep="")  
  
  #get next wt_sequences 
  wt_sequences = read_csv( name_WT, col_names = c("simulation", "proteinA", "proteinB1", "proteinB2", "proteinC1","proteinC2","empty"))
  proteinB<-paste(wt_sequences$proteinB1,wt_sequences$proteinB2,sep="")
  proteinC<-paste(wt_sequences$proteinC1,wt_sequences$proteinC2,sep="")
  wt_sequences$proteinB=proteinB
  wt_sequences$proteinC=proteinC
  
  #differences in A
  diffsitesA<-detect_subs(WT_prev$proteinA,wt_sequences$proteinA)
  
  #differences in B
  diffsitesB<-detect_subs(WT_prev$proteinB,wt_sequences$proteinB)
  
  #differences in C
  diffsitesC<-detect_subs(WT_prev$proteinC,wt_sequences$proteinC)
  
  # up date the prev sequences
  WT_prev<-wt_sequences
  
  
  # make counts of which positions differed
  a<-table( factor( diffsitesA, levels = 1:77)) 
  b<-table( factor( diffsitesB, levels = 1:74)) 
  c<-table( factor( diffsitesC, levels = 1:74)) 
  
  # subs at interface in A
  interfaceA<-sum(a[A_sites])
  #all subs in A
  all_A<-sum(a)
  #table of the ration of interface subs over all subs
  A_table<-c(A_table,interfaceA/all_A)
  
  #B
  interfaceB<-sum(b[B_sites])
  all_B<-sum(b)
  B_table<-c(B_table,interfaceB/all_B)
  
  #C
  interfaceC<-sum(c[B_sites])
  all_C<-sum(c)
  C_table<-c(C_table,interfaceC/all_C)
  
  #update all counts
  A_sum_WT<-A_sum_WT+all_A
  B_sum_WT<-B_sum_WT+all_B
  C_sum_WT<-C_sum_WT+all_C
}



for(i in seq(11,250,10)){       
  # Import & tidy nobind sequences
  name_EXP<-paste("C:/Users/User/Documents/dup_alignments/M_SIM2_4gvb_csv_files_first_500/M_SIM2_4gvb_",i,".csv",sep="")  
  nobind_sequences = read_csv( name_EXP, col_names = c("simulation", "proteinA", "proteinB1", "proteinB2", "proteinC1","proteinC2","empty"))
  proteinB<-paste(nobind_sequences$proteinB1,nobind_sequences$proteinB2,sep="")
  proteinC<-paste(nobind_sequences$proteinC1,nobind_sequences$proteinC2,sep="")
  nobind_sequences$proteinB=proteinB
  nobind_sequences$proteinC=proteinC
  
  diffsitesA<-detect_subs(nobind_prev$proteinA,nobind_sequences$proteinA)
  
  diffsitesB<-detect_subs(nobind_prev$proteinB,nobind_sequences$proteinB)
  
  diffsitesC<-detect_subs(nobind_prev$proteinC,nobind_sequences$proteinC)
  
  nobind_prev<-nobind_sequences
  
  a<-table( factor( diffsitesA, levels = 1:77)) 
  b<-table( factor( diffsitesB, levels = 1:74)) 
  c<-table( factor( diffsitesC, levels = 1:74)) 
  
  interfaceA<-sum(a[A_sites])
  all_A<-sum(a)
  A_table_no_bind<-c(A_table_no_bind,interfaceA/all_A)
  
  interfaceB<-sum(b[B_sites])
  all_B<-sum(b)
  B_table_no_bind<-c(B_table_no_bind,interfaceB/all_B)
  
  interfaceC<-sum(c[B_sites])
  all_C<-sum(c)
  C_table_no_bind<-c(C_table_no_bind,interfaceC/all_C)
  
  
  A_sum_nobind<-A_sum_nobind+all_A
  B_sum_nobind<-B_sum_nobind+all_B
  C_sum_nobind<-C_sum_nobind+all_C
  
  
}


cols<-c("#D55E00","#0072B2")

#set up graphing for A
ratiosA<-data.frame(Inter=A_table,Simulation=rep("bind both",length(A_table)))

ratiosA_no_bind<-data.frame(Inter=A_table_no_bind,Simulation=rep("bind B and not B'",length(A_table_no_bind)))

bothA<- rbind(ratiosA,ratiosA_no_bind)

t.test(ratiosA$Inter,ratiosA_no_bind$Inter)
#graph A
A<-ggplot(bothA, aes(Inter,color=Simulation,fill=Simulation))+ geom_density(alpha=.5)+scale_fill_manual(values=cols)+ scale_colour_manual(values=cols)+ylab("Density")+xlab("Probability of Interface Substitution")+ scale_y_continuous(expand = c(0.001, 0))+ scale_x_continuous(limit=c(0,.2))

#set up graphing for B
ratiosB<-data.frame(Inter=B_table,Simulation=rep("bind both",length(B_table)))

ratiosB_no_bind<-data.frame(Inter=B_table_no_bind,Simulation=rep("bind B and not B'",length(B_table_no_bind)))

bothB<- rbind(ratiosB,ratiosB_no_bind)
t.test(ratiosB$Inter,ratiosB_no_bind$Inter)
var.test(ratiosB$Inter,ratiosB_no_bind$Inter)
#graph B
B<-ggplot(bothB, aes(Inter,color=Simulation,fill=Simulation))+ geom_density(alpha=.5)+scale_fill_manual(values=cols)+ 
  scale_colour_manual(values=cols)+xlab("Probability of Interface Substitution")+ylab("Density")+scale_y_continuous(expand = c(0.001, 0),limit=c(0,15))+ scale_x_continuous(limit=c(0,.55))


#set up graphing for C
ratiosC<-data.frame(Inter=C_table,Simulation=rep("bind both",length(C_table)))


ratiosC_no_bind<-data.frame(Inter=C_table_no_bind,Simulation=rep("bind B and not B'",length(C_table_no_bind)))

bothC<- rbind(ratiosC,ratiosC_no_bind)
t.test(ratiosC$Inter,ratiosC_no_bind$Inter)
#graph C
C<-ggplot(bothC, aes(Inter,color=Simulation,fill=Simulation))+geom_density(alpha=.5)+scale_fill_manual(values=cols)+ 
  scale_colour_manual(values=cols)+xlab("Probability of Interface Substitution")+ylab("Density")+scale_y_continuous(expand = c(0.001, 0),limit=c(0,15))+scale_x_continuous(limit=c(0,.55))


legend<-get_legend(C)
pp<- plot_grid( A+ theme(legend.position= c(0.1, 0.8)),
                B+ theme(legend.position="none"),
                C+ theme(legend.position="none"),
                align = 'vh',
                labels = c("A", "B","C"),
                hjust = -1,
                ncol = 1)
p<-plot_grid(pp)
p


print(A_sum_WT/(A_sum_WT+B_sum_WT+C_sum_WT))
print(B_sum_WT/(A_sum_WT+B_sum_WT+C_sum_WT))
print(C_sum_WT/(A_sum_WT+B_sum_WT+C_sum_WT))

print(A_sum_nobind/(A_sum_nobind+B_sum_nobind+C_sum_nobind))
print(B_sum_nobind/(A_sum_nobind+B_sum_nobind+C_sum_nobind))
print(C_sum_nobind/(A_sum_nobind+B_sum_nobind+C_sum_nobind))

print(C_sum_nobind)

p

ggplot2::ggsave("C:/Users/User/Documents/dup_alignments/Interface_prob_expand_other_prot.pdf",plot=p, width = 11, height = 9,units = "in")


#>     print(A_sum_WT/(A_sum_WT+B_sum_WT+C_sum_WT))
#[1] 0.8791442
#>     print(B_sum_WT/(A_sum_WT+B_sum_WT+C_sum_WT))
#[1] 0.06101371
#>     print(C_sum_WT/(A_sum_WT+B_sum_WT+C_sum_WT))
#[1] 0.05984204
#>     
#  >     print(A_sum_nobind/(A_sum_nobind+B_sum_nobind+C_sum_nobind))
#[1] 0.8321816
#>     print(B_sum_nobind/(A_sum_nobind+B_sum_nobind+C_sum_nobind))
#[1] 0.07542606
#>     print(C_sum_nobind/(A_sum_nobind+B_sum_nobind+C_sum_nobind))
#[1] 0.09239239