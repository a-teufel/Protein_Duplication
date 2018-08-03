require(splines)
require(quantreg)
require(latticeExtra)
require(cowplot)
require(vioplot)
require(ggrepel)


rm(list = ls())
#script to graph the stability of protein A-B during simulations using the 4gvb protien (S5 in manscript)


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
    #get stabilities, was initally set up to be able to use proteins that don't involve the duplicate as well
    ab.count <- append(ab.count, dat$count)
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



##Get all of the data
Un1 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM1_beta1_4gvb/',sep=""))
Un2 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM2_beta1_4gvb/',sep=""))
Un5 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM5_beta1_4gvb/',sep=""))



plot.data <- data.frame(x=c(
  Un1$ab.count,
  Un2$ab.count,
  Un5$ab.count), 
  y=c(
    Un1$an.binding,
    Un2$an.binding,
    Un5$an.binding), 
  id=c(
    rep('Dup bind both', length(Un1$an.binding)),
    rep('Dup bind one not other', length(Un2$an.binding)),
    rep('Dup bind one bind', length(Un5$an.binding)))
)



plot.data <- plot.data[order(plot.data$x), ]


#average points and draw lines
degree.freedom <-100
alpha<-1
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


plot.data.un5 <- plot.data[plot.data$id == 'Dup bind one bind', ]

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.un5)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.un5)
  y.fit.un5 <- X %*% fit$coef
  plot.data.un5<- cbind(plot.data.un5, y.fit.un5)
}


##Assemble data together for plotting
plot.data <- data.frame(x=c(
  plot.data.un1$x,
  plot.data.un2$x,
  plot.data.un5$x),
  lab_place=c(
    rep(1990,length(Un1$an.binding)),
    rep(1900,length(Un2$an.binding)),
    rep(200,length(Un5$an.binding))),
  y=c(
    plot.data.un1$y,
    plot.data.un2$y,
    plot.data.un5$y),
  ysmooth=c(
    plot.data.un1[, 5], 
    plot.data.un2[, 5],
    
    plot.data.un5[, 5]),
  ymin=c(
    plot.data.un1[, 4], 
    plot.data.un2[, 4],
    plot.data.un5[, 4]),
  ymax=c(
    plot.data.un1[, 6], 
    plot.data.un2[, 6],
    plot.data.un5[, 6]),
  id=c(
    rep('bind both', length(Un1$an.binding)),
    rep('bind B and not B\'', length(Un2$an.binding)),
    rep('bind B', length(Un5$an.binding)))
)




##Put in levels so the legend is in the right order
plot.data$id <- factor(plot.data$id, levels = c('bind both', 'bind B and not B\'','bind B'))

mycols <- c("#E69F00", "#0072B2", "#56B4E9")


##Plot function 
survival.lines <- function(df) {
  require(ggplot2)
  require(grid)
  
  
  g <- ggplot(df, aes(x=x, y=y, color=id, fill=id,label=id)) + 
    geom_line(aes(y=ysmooth), size=1) + 
    scale_colour_manual(breaks= c('bind both', 'bind B and not B\'','bind B'),values=mycols,labels= c('bind both', 'bind B and not B\'','bind B')) +
    scale_fill_manual(breaks= c('bind both', 'bind B and not B\'','bind B'),labels= c('bind both', 'bind B and not B\'','bind A'),values=mycols)
  
  g <- g + ylab(expression(paste(Delta,"G of A-B Interface ")))
  g <- g + xlab('Time')
  g<-  g + coord_cartesian(ylim=c(-9, -5),xlim=c(0, 2200))
  
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
    segment.size = 0.35,  force = 15,nudge_y=c(-.5,.5,.5),
    box.padding = unit(0.5, 'lines'),
    # Add extra padding around each data point.
    point.padding = unit(.75, 'lines'),min.segment.length = unit(0, "lines")
  )
  
  save_plot("STAB_interfaceAB_with_controls_other_prot_nice.pdf",g)
  return(g)
}

survival.lines(plot.data)




