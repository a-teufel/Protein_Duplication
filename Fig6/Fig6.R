require(survival)
require(survminer)
require(latticeExtra)
require(ggfortify)
require(ggplot2)
require(GGally)
require(scales)
require(cowplot)
require(ggrepel)
require(reshape2)
require(ggforce)


rm(list = ls())


#script to draw figure that shows anc. binding and divergence

#at .0001
survival.value=c(-5.34086620541)



#get data for A protein, to see if it can still bind B
get.data <- function(this.folder) {
  start <- this.folder
  dirs <- list.files(start)
  
  survival.divergence <- c()
  status.divergence <- c()
  
  for(i in dirs) {
    dat <- read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
    
    
    divergence <- 1 - dat$identityA
    survived <- dat$ancestral_interaction_AB <= survival.value
    
    if(is.na(dat$interface_identity_AC[1])){
      divergence <- 1 - dat$identityB
      survived <- dat$ancestral_interaction_AB <=  survival.value
    }
    
    if(length(survived)!=0){
      
      index <- min(which(survived == FALSE))
      
      if(!is.infinite(index)){
        survived[index:length(survived)]<-FALSE
      }
      
      
      cutoff.divergence<-divergence[index]
      max_div<-max(divergence)
      
      cutoff.divergence[is.infinite(cutoff.divergence) | is.na(cutoff.divergence)] <- max_div
      
      survival.divergence <- append(survival.divergence, cutoff.divergence)
      
      status.divergence <- append(status.divergence, as.numeric(!cutoff.divergence == max_div))
    }
  }
  
  tmp.survival.data <- data.frame(survival.divergence=survival.divergence,
                                  status.divergence=status.divergence)
  
  return(tmp.survival.data)
}


#get data from B protein to see if it can bind A
get.dataB <- function(this.folder) {
  start <- this.folder
  dirs <- list.files(start)
  
  survival.divergence <- c()
  status.divergence <- c()
  
  for(i in dirs) {
    dat <- read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
    
    divergence <- 1 - dat$identityB
    
    survived <- dat$ancestral_interaction_BA <= survival.value
    
    if(is.na(dat$interface_identity_AC[1])){
      divergence <- 1 - dat$identityB
      survived <- dat$ancestral_interaction_BA <=  survival.value
    }
    
    if(length(survived)!=0){
      
      index <- min(which(survived == FALSE))
      
      if(!is.infinite(index)){
        survived[index:length(survived)]<-FALSE
      }
      
      cutoff.divergence<-divergence[index]
      max_div<-max(divergence)
      
      
      cutoff.divergence[is.infinite(cutoff.divergence) | is.na(cutoff.divergence)] <- max_div
      
      survival.divergence <- append(survival.divergence, cutoff.divergence)
      
      status.divergence <- append(status.divergence, as.numeric(!cutoff.divergence == max_div))
    }
  }
  
  tmp.survival.data <- data.frame(survival.divergence=survival.divergence,
                                  status.divergence=status.divergence)
  
  return(tmp.survival.data)
}


#get data from C protien to see if it can bind A
get.dataC <- function(this.folder, which.chain) {
  start <- this.folder
  dirs <- list.files(start)
  
  survival.divergence <- c()
  status.divergence <- c()
  
  for(i in dirs) {
    dat <- read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
    
    
    divergence <- 1 - dat$indentityC
    
    survived <- dat$ancestral_interaction_CA <= survival.value
    if(is.na(dat$interface_identity_AC[1])){
      divergence <- 1 - dat$identityB
      survived <- dat$ancestral_interaction_BA <=  survival.value
    }
    
    if(length(survived)!=0){
      
      index <- min(which(survived == FALSE))
      
      if(!is.infinite(index)){
        survived[index:length(survived)]<-FALSE
      }
      cutoff.divergence<-divergence[index]
      max_div<-max(divergence)
      
      cutoff.divergence[is.infinite(cutoff.divergence) | is.na(cutoff.divergence)] <- max_div
      
      survival.divergence <- append(survival.divergence, cutoff.divergence)
      
      status.divergence <- append(status.divergence, as.numeric(!cutoff.divergence == max_div))
    }
  }
  
  tmp.survival.data <- data.frame(survival.divergence=survival.divergence,
                                  status.divergence=status.divergence)
  
  return(tmp.survival.data)
}


Compare_get_div<-function(B,C){
  
  
  notbind<-NULL
  bind<-NULL
  
  
 for(i in 1:length(B$status.divergence)){
     print(paste(B$status.divergence[i],C$status.divergence[i]))
  if(B$status.divergence[i]!=C$status.divergence[i]){
    
     if(B$status.divergence[i] == 1){
        bind<-c(bind, B$survival.divergence[i] ) 
     }
    if(C$status.divergence[i] == 1){
      bind<-c(bind, C$survival.divergence[i] ) 
    }
    
    if(B$status.divergence[i] == 0){
      notbind<-c(notbind, B$survival.divergence[i] ) 
    }
    if(C$status.divergence[i] == 0){
      notbind<-c(notbind, C$survival.divergence[i] ) 
    }
    
  }
   
 }
  results <- data.frame(notbind=notbind,bind=bind)
  return(results)
  
  
}

Compare_get_div_SD<-function(B,C){
  
  
  notbind<-NULL
  bind<-NULL
  

  for(i in 1:length(B$status.divergence)){

    if(B$status.divergence[i]!=C$status.divergence[i]){
      
      if(B$status.divergence[i] == 1){
        bind<-c(bind, B$survival.divergence[i] ) 
      }
      if(C$status.divergence[i] == 1){
        bind<-c(bind, C$survival.divergence[i] ) 
      }
      
      if(B$status.divergence[i] == 0){
        notbind<-c(notbind, B$survival.divergence[i] ) 
      }
      if(C$status.divergence[i] == 0){
        notbind<-c(notbind, C$survival.divergence[i] ) 
      }
      
    }
    
  }
  results <- data.frame(notbind=sd(notbind)/sqrt(length(notbind)),bind=sd(bind)/sqrt(length(notbind)))
  return(results)
  
  
}


Compare_get_div_t<-function(B,C){
  
  
  notbind<-NULL
  bind<-NULL
  
  for(i in 1:length(B$status.divergence)){
    if(B$status.divergence[i]!=C$status.divergence[i]){
      
      if(B$status.divergence[i] == 1){
        bind<-c(bind, B$survival.divergence[i] ) 
      }
      if(C$status.divergence[i] == 1){
        bind<-c(bind, C$survival.divergence[i] ) 
      }
      
      if(B$status.divergence[i] == 0){
        notbind<-c(notbind, B$survival.divergence[i] ) 
      }
      if(C$status.divergence[i] == 0){
        notbind<-c(notbind, C$survival.divergence[i] ) 
      }
      
    }
    
  }
  

  test<-t.test(notbind,bind,paired=TRUE)
  results <- data.frame(p=test$p.value)
  return(results)
  
  
}


survival.data.WTC <- get.dataC(paste('data/M_SIM_WT_anc/',sep=""))
survival.data.Un1C <- get.dataC(paste('data/M_SIM1_anc/',sep=""))
survival.data.Un2C <- get.dataC(paste('data/M_SIM2_anc/',sep=""))
survival.data.Un3C <- get.dataC(paste('data/M_SIM3_anc/',sep=""))
survival.data.Un4C <- get.dataC(paste('data/M_SIM4_anc/',sep=""))
survival.data.Un5C <- get.dataC(paste('data/M_SIM_WT_nobind_anc/',sep=""))
survival.data.Un6C <- get.dataC(paste('data/M_SIM5_anc/',sep=""))

survival.data.WTB <- get.dataB(paste('data/M_SIM_WT_anc/',sep=""))
survival.data.Un1B <- get.dataB(paste('data/M_SIM1_anc/',sep=""))
survival.data.Un2B <- get.dataB(paste('data/M_SIM2_anc/',sep=""))
survival.data.Un3B <- get.dataB(paste('data/M_SIM3_anc/',sep=""))
survival.data.Un4B <- get.dataB(paste('data/M_SIM4_anc/',sep=""))
survival.data.Un5B <- get.dataB(paste('data/M_SIM_WT_nobind_anc/',sep=""))
survival.data.Un6B <- get.dataB(paste('data/M_SIM5_anc/',sep=""))



survival.data.WT <- get.dataB(paste('data/M_SIM_WT_anc/',sep=""))
survival.data.Un1 <- get.dataB(paste('data/M_SIM1_anc/',sep=""))
survival.data.Un2 <- get.dataB(paste('data/M_SIM2_anc/',sep=""))
survival.data.Un3 <- get.dataB(paste('data/M_SIM3_anc/',sep=""))
survival.data.Un4 <- get.dataB(paste('data/M_SIM4_anc/',sep=""))
survival.data.Un5 <- get.dataB(paste('data/M_SIM_WT_nobind_anc/',sep=""))
survival.data.Un6 <- get.dataB(paste('data/M_SIM5_anc/',sep=""))


#put it in a dataframe
survival.data <- data.frame(time=c(survival.data.WT$survival.divergence, 
                                   survival.data.Un1$survival.divergence,
                                   survival.data.Un2$survival.divergence,
                                   survival.data.Un3$survival.divergence,
                                   survival.data.Un4$survival.divergence,
                                   survival.data.Un5$survival.divergence,
                                   survival.data.Un6$survival.divergence), 
                            status=c(survival.data.WT$status.divergence, 
                                     survival.data.Un1$status.divergence,
                                     survival.data.Un2$status.divergence,
                                     survival.data.Un3$status.divergence,
                                     survival.data.Un4$status.divergence,
                                     survival.data.Un5$status.divergence,
                                     survival.data.Un6$status.divergence), 
                            replicate=c(rep('WT', length(survival.data.WT$status.divergence)), 
                                        rep('bind both', length(survival.data.Un1$status.divergence)),
                                        rep('bind B and not B\'', length(survival.data.Un2$status.divergence)),
                                        rep('bind max', length(survival.data.Un3$status.divergence)),
                                        rep('no bind', length(survival.data.Un4$status.divergence)),
                                        rep('WT no bind', length(survival.data.Un5$status.divergence)),
                                        rep('bind B', length(survival.data.Un6$status.divergence)))
)

#pick out just the data I want to graph
survival.data$replicate <- factor(survival.data$replicate, levels = c('bind both', 'bind B and not B\'','bind max','no bind','bind B'))


#fit surv curves to data
fit = survfit(Surv(time,status)~replicate, data=survival.data)


#figure out how long things are to float labels
temp = survfit(Surv(time,status)~replicate, data=survival.data[survival.data$replicate=='bind both', ])
bb1<-length(summary(temp,censored=TRUE)$time)

temp = survfit(Surv(time,status)~replicate, data=survival.data[survival.data$replicate=='bind B and not B\'', ])
bb2<-length(summary(temp,censored=TRUE)$time)


temp = survfit(Surv(time,status)~replicate, data=survival.data[survival.data$replicate=='bind max', ])
bb3<-length(summary(temp,censored=TRUE)$time)


temp = survfit(Surv(time,status)~replicate, data=survival.data[survival.data$replicate=='no bind', ])
bb4<-length(summary(temp,censored=TRUE)$time)

temp = survfit(Surv(time,status)~replicate, data=survival.data[survival.data$replicate=='bind B', ])
bb5<-length(summary(temp,censored=TRUE)$time)


cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7")
values = c("bind both" = "#E69F00", "bind B and not B\'" = "#0072B2", 'bind max'="#009E73",'no bind'="#CC79A7", 'bind B'="#56B4E9")


#make plot
g2 <- ggsurv(fit, plot.cens = TRUE,cens.shape = 16)  + scale_colour_manual(values=values) + scale_fill_manual(values=values)

#float labels
g2<-g2+geom_text_repel( 
  label=c( rep(NA,bb1),'bind both', rep(NA,bb2), 'bind B and not B\'',rep(NA,bb3),'bind max',rep(NA,bb4-3),'no bind',rep(NA,bb5+3),'bind B'),
  inherit.as=FALSE, 
  segment.size = 0.5,  force = 3, nudge_x = 0.1,
  box.padding = unit(0.5, 'lines'),
  # Add extra padding around each data point.
  point.padding = unit(1, 'lines'))


g2<-g2+scale_y_continuous(labels = scales::percent,limits=c(0,1.02),expand = c(0, 0))

g2<-g2+scale_x_continuous(labels = scales::percent,limits=c(0,1.02),expand = c(0, 0))

g2<-g2+xlab("Divergence (100 - % Amino Acid Identity)")+ylab("B \n% Binding Ancestor")

g2<-g2 + theme(legend.position = "none",
               legend.title=element_blank(), 
               legend.key = element_blank())


g2<-g2


survival.data.WT <- get.dataC(paste('data/M_SIM_WT_anc/',sep=""))
survival.data.Un1 <- get.dataC(paste('data/M_SIM1_anc/',sep=""))
survival.data.Un2 <- get.dataC(paste('data/M_SIM2_anc/',sep=""))
survival.data.Un3 <- get.dataC(paste('data/M_SIM3_anc/',sep=""))
survival.data.Un4 <- get.dataC(paste('data/M_SIM4_anc/',sep=""))
survival.data.Un5 <- get.dataC(paste('data/M_SIM_WT_nobind_anc/',sep=""))
survival.data.Un6 <- get.dataC(paste('data/M_SIM5_anc/',sep=""))


survival.data <- data.frame(time=c(survival.data.WT$survival.divergence, 
                                   survival.data.Un1$survival.divergence,
                                   survival.data.Un2$survival.divergence,
                                   survival.data.Un3$survival.divergence,
                                   survival.data.Un4$survival.divergence,
                                   survival.data.Un5$survival.divergence,
                                   survival.data.Un6$survival.divergence), 
                            status=c(survival.data.WT$status.divergence, 
                                     survival.data.Un1$status.divergence,
                                     survival.data.Un2$status.divergence,
                                     survival.data.Un3$status.divergence,
                                     survival.data.Un4$status.divergence,
                                     survival.data.Un5$status.divergence,
                                     survival.data.Un6$status.divergence), 
                            replicate=c(rep('WT', length(survival.data.WT$status.divergence)), 
                                        rep('bind both', length(survival.data.Un1$status.divergence)),
                                        rep('bind B and not B\'', length(survival.data.Un2$status.divergence)),
                                        rep('bind max', length(survival.data.Un3$status.divergence)),
                                        rep('no bind', length(survival.data.Un4$status.divergence)),
                                        rep('WT no bind', length(survival.data.Un5$status.divergence)),
                                        rep('bind B', length(survival.data.Un6$status.divergence)))
)

#keep just the data I want they didnt care about the WT sims in eds paper
survival.data$replicate <- factor(survival.data$replicate, levels = c('bind both', 'bind B and not B\'','bind max','no bind','bind B'))

#fit surival curve
fit = survfit(Surv(time,status)~replicate, data=survival.data)

#figure out lengths to float labels nicely
temp = survfit(Surv(time,status)~replicate, data=survival.data[survival.data$replicate=='bind both', ])
bb1<-length(summary(temp,censored=TRUE)$time)

temp = survfit(Surv(time,status)~replicate, data=survival.data[survival.data$replicate=='bind B and not B\'', ])
bb2<-length(summary(temp,censored=TRUE)$time)


temp = survfit(Surv(time,status)~replicate, data=survival.data[survival.data$replicate=='bind max', ])
bb3<-length(summary(temp,censored=TRUE)$time)


temp = survfit(Surv(time,status)~replicate, data=survival.data[survival.data$replicate=='no bind', ])
bb4<-length(summary(temp,censored=TRUE)$time)

temp = survfit(Surv(time,status)~replicate, data=survival.data[survival.data$replicate=='bind B', ])
bb5<-length(summary(temp,censored=TRUE)$time)


#make C plot
g3 <- ggsurv(fit, plot.cens = TRUE,cens.shape = 16)  + scale_colour_manual(values=values) + scale_fill_manual(values=values)

g3<-g3+geom_text_repel( 
  label=c( rep(NA,bb1),'bind both', rep(NA,bb2), 'bind B and not B\'',rep(NA,bb3),'bind max',rep(NA,bb4),'no bind',rep(NA,bb5),'bind B'),
  inherit.as=FALSE,  
  segment.size = 0.5,  force = 0, nudge_x = 0.1,
  box.padding = unit(0.5, 'lines'),
  # Add extra padding around each data point.
  point.padding = unit(1, 'lines'))


g3<-g3+scale_y_continuous(labels = scales::percent,limits=c(0,1.02),expand = c(0, 0))

g3<-g3+scale_x_continuous(labels = scales::percent,limits=c(0,1.02),expand = c(0, 0))

g3<-g3+xlab("Divergence (100 - % Amino Acid Identity)")+ylab("B'\n % Binding Ancestor")

g3<-g3 + theme(legend.position = "none",
               legend.title=element_blank(), 
               legend.key = element_blank())

g3<-g3


#make violin plots

exp1<-Compare_get_div(survival.data.Un1B,survival.data.Un1C)

exp2<-Compare_get_div(survival.data.Un2B,survival.data.Un2C)

exp3<-Compare_get_div(survival.data.Un3B,survival.data.Un3C)

exp4<-Compare_get_div(survival.data.Un4B,survival.data.Un4C)

exp6<-Compare_get_div(survival.data.Un6B,survival.data.Un6C)

exp1sd<-Compare_get_div_SD(survival.data.Un1B,survival.data.Un1C)

exp2sd<-Compare_get_div_SD(survival.data.Un2B,survival.data.Un2C)

exp3sd<-Compare_get_div_SD(survival.data.Un3B,survival.data.Un3C)

exp4sd<-Compare_get_div_SD(survival.data.Un4B,survival.data.Un4C)

exp6sd<-Compare_get_div_SD(survival.data.Un6B,survival.data.Un6C)

exp1t<-Compare_get_div_t(survival.data.Un1B,survival.data.Un1C)

exp2t<-Compare_get_div_t(survival.data.Un2B,survival.data.Un2C)

exp3t<-Compare_get_div_t(survival.data.Un3B,survival.data.Un3C)

exp4t<-Compare_get_div_t(survival.data.Un4B,survival.data.Un4C)

exp6t<-Compare_get_div_t(survival.data.Un6B,survival.data.Un6C)

t<-cbind(exp1t,exp2t,exp3t,exp4t,exp6t)


sd<-cbind(exp1sd,exp2sd,exp3sd,exp4sd,exp6sd)

##Get all of the data for A
q<-data.frame(exp1,exp2,exp3,exp4,exp6)

e1<-data.frame(exp1)
e1$label<-rep('bind both',length(e1$notbind))


e2<-data.frame(exp2)
e2$label<-rep('bind B and not B\'',length(e2$notbind))

e3<-data.frame(exp3)
e3$label<-rep('bind max',length(e3$notbind))

e4<-data.frame(exp4)
e4$label<-rep('no bind',length(e4$notbind))

e6<-data.frame(exp6)
e6$label<-rep('bind B',length(e6$notbind))



t<-rbind(e1,e3,e6,e2,e4)

replicate=c(
  rep('bind both', 2),
  rep('bind B and not B\'', 2),
  rep('bind max', 2),
  rep('no bind', 2),
  rep('bind B', 2))

qq<-melt(q)
qq$replicate<-replicate
qq$sd<-unname(t(sd))


cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7")


cbPalette <- c("#E69F00", "#009E73","#56B4E9","#0072B2",   "#CC79A7")


dl_Palette<-NULL
for(i in 1:length(cbPalette)){
  dl_Palette<-c(dl_Palette,colorspace::darken(cbPalette[i],amount=0.2),colorspace::lighten(cbPalette[i],amount=0.2))
  
}  
values = c("bind both notbind" = dl_Palette[1], "bind both bind"=dl_Palette[2],'bind max notbind'=dl_Palette[3],'bind max bind'=dl_Palette[4],
           "bind B notbind" = dl_Palette[5], "bind B bind" = dl_Palette[6], 'bind B and not B\' notbind' = dl_Palette[7], 'bind B and not B\' bind' = dl_Palette[8],
           'no bind notbind' = dl_Palette[9], 'no bind bind' = dl_Palette[10])

t<-melt(t)
t$label=factor(t$label, levels=c("bind both", "bind max", "bind B","bind B and not B\'","no bind"))
t$color<-paste(t$label,t$variable)

level_order<-c("bind both", "bind max", "bind B","bind B and not B\'","no bind")


g_bar<-ggplot(data=t, aes(x=factor(label, levels=level_order), y=value,  fill=color )) +
  geom_violin(position=position_dodge(),alpha=0.5,linetype="blank") + scale_colour_manual(values=values) + scale_fill_manual(values=values)+
  geom_sina(aes(colour = factor(color)),size=1)+
  ylab("Divergence \n (100 - % Amino Acid Identity)") + xlab("Selection Scheme")

g_bar<-g_bar+scale_y_continuous(labels = scales::percent,limits=c(0,.5),expand = c(0, 0))


g_bar<-g_bar+scale_x_discrete(expand = c(0, .5))

g_bar<-g_bar+ theme(legend.position = "none",
                    legend.title=element_blank(), 
                    legend.key = element_blank())

pp<-plot_grid(g2,g3,g_bar,labels=c('A','B','C'),nrow=3,ncol=1,align='v')

save_plot("Fig6.pdf",pp, nrow=3,ncol=2,   base_aspect_ratio = 1.3)
