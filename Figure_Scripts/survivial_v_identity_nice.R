require(survival)
require(survminer)
require(latticeExtra)
require(ggfortify)
require(ggplot2)
require(GGally)
require(scales)
require(cowplot)
require(ggrepel)

rm(list = ls())


#script to figure that shows the probability that a protein can still interact with its ancestal partner. 
#The ability to interact is if that binding is over the threshold set by the survival.value (fig s7)


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




##Get all of the data for A
survival.data.WT <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM_WT_beta1/',sep=""))
survival.data.Un1 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM1_beta1/',sep=""))
survival.data.Un2 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM2_beta1/',sep=""))
survival.data.Un3 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM3_beta1/',sep=""))
survival.data.Un4 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM4_beta1/',sep=""))
survival.data.Un5 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM_WT_no_bind_beta1/',sep=""))
survival.data.Un6 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM5_beta1/',sep=""))

#make it into one big data frame with labels
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

#just pick out the ones you want to graph (it doesnt make sense to include the WT)
survival.data$replicate <- factor(survival.data$replicate, levels = c('bind both', 'bind B and not B\'','bind max','no bind','bind B'))


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


#make custom color palet
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7")
values = c("bind both" = "#E69F00", "bind B and not B\'" = "#0072B2", 'bind max'="#009E73",'no bind'="#CC79A7", 'bind B'="#56B4E9")

g <- ggsurv(fit, plot.cens = TRUE,cens.shape = 16)  + scale_colour_manual(values=values) + scale_fill_manual(values=values)

#repel labels at resonable intervals so things look nice
g<-g+geom_text_repel( 
  label=c( rep(NA,bb1),'bind both', rep(NA,bb2), 'bind B and not B\'',rep(NA,bb3),'bind max',rep(NA,bb4-1),'no bind',rep(NA,bb5+1),'bind B'),
  inherit.as=FALSE,#arrow = arrow(length = unit(0.01, 'npc')),  
  segment.size = 0.5,  force = 0, nudge_x = 0.1,
  box.padding = unit(0.5, 'lines'),
  # Add extra padding around each data point.
  point.padding = unit(1, 'lines'))


g<-g+scale_y_continuous(labels = scales::percent,limits=c(0,1.02),expand = c(0, 0))

g<-g+scale_x_continuous(labels = scales::percent,limits=c(0,1.02),expand = c(0, 0))

g<-g+xlab("Divergence (100 - % Amino Acid Identity)")+ylab("% Binding Ancestor")

g<-g + theme(legend.position = "none",
             legend.title=element_blank(), 
             legend.key = element_blank())
g1<-g


#get data from B
survival.data.WT <- get.dataB(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM_WT_beta1/',sep=""))
survival.data.Un1 <- get.dataB(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM1_beta1/',sep=""))
survival.data.Un2 <- get.dataB(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM2_beta1/',sep=""))
survival.data.Un3 <- get.dataB(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM3_beta1/',sep=""))
survival.data.Un4 <- get.dataB(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM4_beta1/',sep=""))
survival.data.Un5 <- get.dataB(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM_WT_no_bind_beta1/',sep=""))
survival.data.Un6 <- get.dataB(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM5_beta1/',sep=""))


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

g2<-g2+xlab("Divergence (100 - % Amino Acid Identity)")+ylab("% Binding Ancestor")

g2<-g2 + theme(legend.position = "none",
               legend.title=element_blank(), 
               legend.key = element_blank())


g2<-g2



#get data from C
survival.data.WT <- get.dataC(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM_WT_beta1/',sep=""))
survival.data.Un1 <- get.dataC(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM1_beta1/',sep=""))
survival.data.Un2 <- get.dataC(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM2_beta1/',sep=""))
survival.data.Un3 <- get.dataC(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM3_beta1/',sep=""))
survival.data.Un4 <- get.dataC(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM4_beta1/',sep=""))
survival.data.Un5 <- get.dataC(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM_WT_no_bind_beta1/',sep=""))
survival.data.Un6 <- get.dataC(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM5_beta1/',sep=""))

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

#keep just the data I want
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

g3<-g3+xlab("Divergence (100 - % Amino Acid Identity)")+ylab("% Binding Ancestor")

g3<-g3 + theme(legend.position = "none",
               legend.title=element_blank(), 
               legend.key = element_blank())

g3<-g3



pp<-plot_grid(g1,g2,g3,labels=c('A','B','C'),nrow=3,ncol=1,align='v')

save_plot("sur_all_filpB_nice.pdf",pp, nrow=3,ncol=2,   base_aspect_ratio = 1.3)
