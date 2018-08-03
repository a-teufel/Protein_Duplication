require(splines)
require(quantreg)
require(latticeExtra)
require(cowplot)
require(vioplot)
require(ggrepel)


#script to graph the stability of protein A-B during simulations (fig 2)


rm(list = ls())

##Function that grabs the data from the files in the specified folder
get.data <- function(this.folder, which.chain) {
  start <- this.folder
  dirs <- list.files(start)
  ab.count <- c()
  an.binding <- c()
  name <- c()
  count <- 0
  
  for(i in dirs) {
    dat <- read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
    
    count <- count + 1
    
    ab.count <- append(ab.count, dat$count)
    #get the stabs, function was made to handle WT data too
    if(!is.na(dat$evolved_interaction_AC)){
      an.binding<-append(an.binding,dat$evolved_interaction_AB);
    }
    if(is.na(dat$evolved_interaction_AC)){
      an.binding<-append(an.binding,dat$evolved_interaction_AB);
    }	
    name <- append(name, dat$name)
    
  }
  
  
  tmp.data <- data.frame(names=name, ab.count=ab.count, an.binding=an.binding)
  return(tmp.data)
  
}


#read in data
WT <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM_WT_beta1/',sep=""))
Un1 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM1_beta1/',sep=""))
Un2 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM2_beta1/',sep=""))
Un3 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM3_beta1/',sep=""))
Un4 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM4_beta1/',sep=""))
Un5 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM5_beta1/',sep=""))
Un6 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM_WT_no_bind_beta1/',sep=""))



plot.data <- data.frame(x=c(WT$ab.count, 
                            Un1$ab.count,
                            Un2$ab.count,
                            Un3$ab.count,
                            Un4$ab.count,
                            Un5$ab.count,
                            Un6$ab.count), 
                        y=c(WT$an.binding, 
                            Un1$an.binding,
                            Un2$an.binding,
                            Un3$an.binding,
                            Un4$an.binding,
                            Un5$an.binding,
                            Un6$an.binding), 
                        id=c(rep('WT', length(WT$an.binding)), 
                             rep('Dup bind both', length(Un1$an.binding)),
                             rep('Dup bind one not other', length(Un2$an.binding)),
                             rep('Dup bind max', length(Un3$an.binding)),
                             rep('Dup no bind', length(Un4$an.binding)),
                             rep('Dup bind one bind', length(Un5$an.binding)),
                             rep('WT no bind', length(Un6$an.binding)))
)



plot.data <- plot.data[order(plot.data$x), ]


degree.freedom<-100
alpha<-1

##Start of WT
plot.data.wt <- plot.data[plot.data$id == 'WT', ]


#draw lines through points

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.wt)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.wt)
  y.fit.wt <- X %*% fit$coef
  plot.data.wt <- cbind(plot.data.wt, y.fit.wt)
}



plot.data.un1 <- plot.data[plot.data$id == 'Dup bind both', ]

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.un1)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.un1)
  y.fit.un1 <- X %*% fit$coef
  plot.data.un1<- cbind(plot.data.un1, y.fit.un1)
}


plot.data.un2 <- plot.data[plot.data$id == 'Dup bind one not other', ]

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.un2)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.un2)
  y.fit.un2 <- X %*% fit$coef
  plot.data.un2<- cbind(plot.data.un2, y.fit.un2)
}

plot.data.un3 <- plot.data[plot.data$id == 'Dup bind max', ]

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.un3)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.un3)
  y.fit.un3 <- X %*% fit$coef
  plot.data.un3<- cbind(plot.data.un3, y.fit.un3)
}


plot.data.un4 <- plot.data[plot.data$id == 'Dup no bind', ]

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.un4)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.un4)
  y.fit.un4 <- X %*% fit$coef
  plot.data.un4<- cbind(plot.data.un4, y.fit.un4)
}

plot.data.un5 <- plot.data[plot.data$id == 'Dup bind one bind', ]

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.un5)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.un5)
  y.fit.un5 <- X %*% fit$coef
  plot.data.un5<- cbind(plot.data.un5, y.fit.un5)
}


plot.data.un6 <- plot.data[plot.data$id == 'WT no bind', ]

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.un6)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.un6)
  y.fit.un6 <- X %*% fit$coef
  plot.data.un6<- cbind(plot.data.un6, y.fit.un6)
}

##Assemble data together for plotting

plot.data <- data.frame(x=c(plot.data.wt$x,
                            plot.data.un1$x,
                            plot.data.un2$x,
                            plot.data.un3$x,
                            plot.data.un4$x,
                            plot.data.un5$x,
                            plot.data.un6$x),
                        lab_place=c(rep(1960,length(WT$an.binding)),
                                    rep(1990,length(Un1$an.binding)),
                                    rep(350,length(Un2$an.binding)),
                                    rep(1750,length(Un3$an.binding)),
                                    rep(1750,length(Un4$an.binding)),
                                    rep(1990,length(Un5$an.binding)),
                                    rep(1350,length(Un6$an.binding))),
                        y=c(plot.data.wt$y,
                            plot.data.un1$y,
                            plot.data.un2$y,
                            plot.data.un3$y,
                            plot.data.un4$y,
                            plot.data.un5$y,
                            plot.data.un6$y),
                        ysmooth=c(plot.data.wt[, 5], 
                                  plot.data.un1[, 5], 
                                  plot.data.un2[, 5],
                                  plot.data.un3[, 5],
                                  plot.data.un4[, 5],
                                  plot.data.un5[, 5],
                                  plot.data.un6[, 5]),
                        ymin=c(plot.data.wt[, 4], 
                               plot.data.un1[, 4], 
                               plot.data.un2[, 4],
                               plot.data.un3[, 4],
                               plot.data.un4[, 4],
                               plot.data.un5[, 4],
                               plot.data.un6[, 4]),
                        ymax=c(plot.data.wt[, 6], 
                               plot.data.un1[, 6], 
                               plot.data.un2[, 6],
                               plot.data.un3[, 6],
                               plot.data.un4[, 6],
                               plot.data.un5[, 6],
                               plot.data.un6[, 6]),
                        id=c(rep('ND', length(WT$an.binding)), 
                             rep('bind both', length(Un1$an.binding)),
                             rep('bind B and not B\'', length(Un2$an.binding)),
                             rep('bind max', length(Un3$an.binding)),
                             rep('no bind', length(Un4$an.binding)),
                             rep('bind B', length(Un5$an.binding)),
                             rep('ND no bind', length(Un6$an.binding)))
)



##Put in levels so the legend is in the right order
plot.data$id <- factor(plot.data$id, levels = c('bind both', 'bind B','bind max','bind B and not B\'','no bind','ND no bind',  'ND'  ))


#custom colors
mycols <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7","#D55E00","#000000")

##Plot function 
survival.lines <- function(df) {
  require(ggplot2)
  require(grid)
  
  
  g <- ggplot(df, aes(x=x, y=y, color=id, fill=id,label=id)) + 
    geom_line(aes(y=ysmooth), size=1) + 
    scale_colour_manual(breaks= c('ND' ,'bind both', 'bind B and not B\'','bind max','no bind','bind B','ND no bind'),values=mycols,labels= c('ND' ,'bind both', 'bind B and not B\'','bind max','no bind','bind B','ND no bind')) +
    scale_fill_manual(breaks= c('ND' ,'bind both', 'bind B and not B\'','bind max','no bind','bind B','ND no bind'),labels= c('ND' ,'bind both', 'bind B and not B\'','bind max','no bind','bind B','ND no bind'),values=mycols)
  
  g <- g + ylab(expression(paste(Delta,"G of A-B Interface")))
  g <- g + xlab('Time')
  g<-  g + coord_cartesian(ylim=c(-12, 2),xlim=c(0, 2200))
  
  g<-g+scale_y_continuous(expand = c(0, 0))
  
  g<-g+scale_x_continuous(expand = c(0, 0))
  
  
  g <- g + theme(legend.position = 'none',
                 legend.title=element_blank(), 
                 legend.key = element_blank())
  
  
  
  keep<-c("id","ysmooth","lab_place","x")
  data1=df[keep]
  data2=data1[!duplicated(data1), ]
  
  
  g<-g+geom_text_repel(
    data = subset(data2, lab_place==x),
    aes(x=lab_place,y=ysmooth),
    segment.color="black", 
    segment.size = 0.35,  force = 5,nudge_x=375,
    box.padding = unit(0.5, 'lines'),
    # Add extra padding around each data point.
    point.padding = unit(.75, 'lines'),min.segment.length = unit(0, "lines")
  )
  
  save_plot("STAB_AB_with_controls_nice.pdf",g)
  return(g)
}

survival.lines(plot.data)


