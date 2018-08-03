
# script to draw ggridge plots from sites shown to be sig. different. 
# Note that the B' protein is refered to as protein C. 
# I only wrote one script to draw both the stickiness and the mass,  you just have to change which function is called later in the code if you want to do one or the other

rm(list=ls())
library(tidyverse)
library(stringr)
library(cowplot)
library(broom)
library(dplyr)
library(grid)
library(ggridges)
library(colorspace)



sticky_by_site <- function(site_df) {
  Aminos<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  sticky<-c(0.0062, 1.0372, -0.7485, -0.7893, 1.2727, -0.1771, 0.1204, 1.1109, -1.1806, 0.9138, 1.0124, -0.2693, -0.1799, -0.4114, -0.0876, 0.1376, 0.1031, 0.7599, 0.7925, 0.8806)
  
  sticky_sum<-NULL
  for(i in 1:length(site_df[[1]])){
    amino_to_look<-site_df[[1]][i]
    sticky_index<-which(Aminos==amino_to_look)
    sticky_sum<-c(sticky_sum,sticky[sticky_index])
  }
  
  return(sticky_sum)
}


mass_by_site <- function(site_df) {
  #naming here isn't right, but it makes graphing the masses easier. 
  Aminos<-c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  sticky<-c(71.08,156.2,114.1,115.1,103.1,129.1,128.1,57.05,137.1,113.2,113.2,128.2,131.2,147.2,97.12,87.08,101.1,186.2,163.2,99.13)
  
  sticky_sum<-NULL
  for(i in 1:length(site_df[[1]])){
    amino_to_look<-site_df[[1]][i]
    sticky_index<-which(Aminos==amino_to_look)
    
    sticky_sum<-c(sticky_sum,sticky[sticky_index])
  }
  
  return(sticky_sum)
}


#Makes a dataframe of everything that happened in the C protein
sites_and_time_data_C<-function(file,sumC,sites_timeC,BindBoth){
  
  cnt=1
  Sticky_C<-data.frame()
  
  for(i in seq(1,2000,50)){
    #read in data
    name_EXP<-paste(file,i,".csv",sep="")  
    sequences = read_csv(name_EXP, col_names = c("simulation", "proteinA", "proteinB", "proteinC", "empty"))
    
    seq_C<-sequences$proteinC
    
    split_C<-strsplit(seq_C, "")
    
    sum_sticky<-NULL
    which_time<-NULL
    
    for(j in 1:length(split_C)){
      #caluclate the stickiness or mass of each amino acid in the sequences
      sticky_vector<-mass_by_site(split_C[j]) 
      
      
      #multiple by the vector that has 1 marked for the sites we want to keep
      clean_sticky<-as.vector(sticky_vector)*as.vector(sumC)
      
      
      sum_sticky<-rbind(sum_sticky,clean_sticky)
      
      #keep track of time too, so we know when the site changes
      which_time<-rbind(which_time,as.vector(sites_timeC[cnt,]))
    }
    
    #dataframe is organized by sites on the columns.
    colnames(sum_sticky)<-paste(seq(1,73))
    colnames(which_time)<-paste(seq(1,73))
    temp_A_sticky<-NULL
    #make a temp data frame of the distribution from a single generation
    if(BindBoth){
      temp_C_sticky<-data.frame(Gen=rep(i-1,nrow(sum_sticky)),  Sticky=sum_sticky, Exp=rep("Bind Both",nrow(sum_sticky)),Sel=which_time)
    }
    else{
      temp_C_sticky<-data.frame(Gen=rep(i-1,nrow(sum_sticky)),  Sticky=sum_sticky, Exp=rep("Bind A and Not A'",nrow(sum_sticky)),Sel=which_time)  
    }
    
    #add it to the larger dataframe
    Sticky_C<-rbind.data.frame(Sticky_C,temp_C_sticky)
    cnt<-cnt+1
  }
  return(Sticky_C)
}


#Makes a dataframe of everything that happened in the C protein
sites_and_time_data_B<-function(file,sumB,sites_timeB,BindBoth){
  
  cnt=1
  Sticky_B<-data.frame()
  
  for(i in seq(1,2000,50)){
    #read in data
    name_EXP<-paste(file,i,".csv",sep="")  
    sequences = read_csv(name_EXP, col_names = c("simulation", "proteinA", "proteinB", "proteinC", "empty"))

    seq_B<-sequences$proteinB
    
    split_B<-strsplit(seq_B, "")
    
    sum_sticky<-NULL
    which_time<-NULL
    
    for(j in 1:length(split_B)){
      #caluclate the stickiness or mass of each amino acid in the sequences
      sticky_vector<-mass_by_site(split_B[j]) 
      
      
      #multiple by the vector that has 1 marked for the sites we want to keep
      clean_sticky<-as.vector(sticky_vector)*as.vector(sumB)
      
      
      sum_sticky<-rbind(sum_sticky,clean_sticky)
      
      #keep track of time too, so we know when the site changes
      which_time<-rbind(which_time,as.vector(sites_timeB[cnt,]))
    }
    
    #dataframe is organized by sites on the columns.
    colnames(sum_sticky)<-paste(seq(1,154))
    colnames(which_time)<-paste(seq(1,154))
    temp_B_sticky<-NULL
    #make a temp data frame of the distribution from a single generation
    if(BindBoth){
      temp_B_sticky<-data.frame(Gen=rep(i-1,nrow(sum_sticky)),  Sticky=sum_sticky, Exp=rep("Bind Both",nrow(sum_sticky)),Sel=which_time)
    }
    else{
      temp_B_sticky<-data.frame(Gen=rep(i-1,nrow(sum_sticky)),  Sticky=sum_sticky, Exp=rep("Bind A and Not A'",nrow(sum_sticky)),Sel=which_time)  
    }
    
    #add it to the larger dataframe
    Sticky_B<-rbind.data.frame(Sticky_B,temp_B_sticky)
    cnt<-cnt+1
  }
  return(Sticky_B)
}

#Makes a dataframe of everything that happened in the A protien
sites_and_time_data_A<-function(file,sumA,sites_timeA,BindBoth){
  
  cnt=1
  Sticky_A<-data.frame()
  
  for(i in seq(1,2000,50)){
    #read in data
    name_EXP<-paste(file,i,".csv",sep="")  
    sequences = read_csv(name_EXP, col_names = c("simulation", "proteinA", "proteinB", "proteinC", "empty"))
    
    seq_A<-sequences$proteinA
    
    split_A<-strsplit(seq_A, "")
    sum_sticky<-NULL
    which_time<-NULL
    
    for(j in 1:length(split_A)){
      #caluclate the mass or stickiness of each amino acid in the sequences
      sticky_vector<-mass_by_site(split_A[j]) 
      
      #multiple by the vector that has 1 marked for the sites we want to keep
      clean_sticky<-as.vector(sticky_vector)*as.vector(sumA)
      
      
      sum_sticky<-rbind(sum_sticky,clean_sticky)
      
      #keep track of time too, so we know when the site changes
      which_time<-rbind(which_time,as.vector(sites_timeA[cnt,]))
    }
    
    #dataframe is organized by sites on the columns.
    colnames(sum_sticky)<-paste(seq(1,154))
    colnames(which_time)<-paste(seq(1,154))
    
    
    temp_A_sticky<-NULL
    #make a temp data frame of the distribution from a single generation
    if(BindBoth){
      temp_A_sticky<-data.frame(Gen=rep(i-1,nrow(sum_sticky)),  Sticky=sum_sticky, Exp=rep("Bind Both",nrow(sum_sticky)),Sel=which_time)
    }
    else{
      temp_A_sticky<-data.frame(Gen=rep(i-1,nrow(sum_sticky)),  Sticky=sum_sticky, Exp=rep("Bind A and Not A'",nrow(sum_sticky)),Sel=which_time)
      
    }
    #add it to the larger dataframe
    Sticky_A<-rbind.data.frame(Sticky_A,temp_A_sticky)
    cnt<-cnt+1
  }
  return(Sticky_A)
}



#read in data of diffs, these files are made by the chi_square_test.R script
site_timeA_all_time<-t(read.csv(file="C:/Users/User/Documents/dup_alignments/diffs_A_dupA.csv"))
site_timeB_all_time<-t(read.csv(file="C:/Users/User/Documents/dup_alignments/diffs_B_dupA.csv"))
site_timeC_all_time<-t(read.csv(file="C:/Users/User/Documents/dup_alignments/diffs_C_dupA.csv"))

#figure out if a site ever differed over time and mark which sites those were
sumA<-colSums(site_timeA_all_time)
sumB<-colSums(site_timeB_all_time)
sumC<-colSums(site_timeC_all_time)

#sites that ever differed over time get marked with a 1
sumA[sumA>0]<-1
sumB[sumB>0]<-1
sumC[sumC>0]<-1


#make dataframes from the bind both data
C_wt<-sites_and_time_data_C("C:/Users/User/Documents/dup_alignments/M_SIM1_DupA_csv_files/M_SIM1_DupA_", sumC,site_timeC_all_time,TRUE)

B_wt<-sites_and_time_data_B("C:/Users/User/Documents/dup_alignments/M_SIM1_DupA_csv_files/M_SIM1_DupA_", sumB,site_timeB_all_time,TRUE)

#make dataframes from the bind B and not B' data
C_nobind<-sites_and_time_data_C("C:/Users/User/Documents/dup_alignments/M_SIM2_DupA_csv_files/M_SIM2_DupA_",sumC, site_timeC_all_time,FALSE)

B_nobind<-sites_and_time_data_B("C:/Users/User/Documents/dup_alignments/M_SIM2_DupA_csv_files/M_SIM2_DupA_", sumB,site_timeB_all_time,FALSE)

#combine both the wt and nobind data into one dataframe
Sticky_C<-rbind.data.frame(C_wt,C_nobind)  
Sticky_B<-rbind.data.frame(B_wt,B_nobind)  



#cut out everything then doesnt matter, ie. all the data where we never observed a difference
dfB<-Sticky_B[, colSums(Sticky_B != 0) > 0]

#just do it by 100s's rather then 50s, this can be changed to get higher or lower grandulatity out of the figure
dfB<-dfB[which(dfB$Gen %% 100 == 0) ,]

df<-Sticky_C[, colSums(Sticky_C != 0) > 0]

df<-df[which(df$Gen %% 100 == 0),]

#set up levels and and numbering for graphing
s<-seq(0,2000,100)

levels<-rev(as.character(s))

dfB$Gen <- factor(dfB$Gen, levels = levels)
df$Gen <- factor(df$Gen, levels = levels)

#set up colors and custom theme
cols<-c("#D55E00","#0072B2")

cols1 <- readhex(file = textConnection(paste(cols, collapse = "\n")),
                 class = "RGB")
#transform to hue/lightness/saturation colorspace
cols1 <- as(cols1, "HLS")
cols2 <- cols1
#additive decrease of lightness
cols1@coords[, "L"] <- pmax(0, cols1@coords[, "L"] + 0.5)

cols1 <- as(cols1, "RGB")
cols1 <- hex(cols1)

cols<-c(cols,cols1)


#custom theme 
theme_ridges_AT <- function(font_size = 14, font_family = "", line_size = .5, grid = TRUE, center_axis_labels = TRUE, ylab_on = TRUE) {
  half_line <- font_size / 2
  small_rel <- 0.857
  small_size <- small_rel * font_size
  color <- "grey90"
  
  if (grid) {
    panel.grid.major <- element_line(colour = color, size = line_size)
    axis.ticks       <- element_line(colour = color, size = line_size)
    axis.ticks.y     <- axis.ticks
  }
  else {
    panel.grid.major <- element_blank()
    axis.ticks       <- element_line(colour = "black", size = line_size)
    axis.ticks.y     <- element_blank()
  }
  
  if (center_axis_labels) {
    axis_just <- 0.5
  }
  else {
    axis_just <- 1.0
  }
  
  if(!ylab_on){
    print("off")
    axis.title.y <- element_blank()
  }
  
  else{
    
    axis.title.y<- element_text(
      angle = 90,
      margin = ggplot2::margin(r = small_size / 2, l = small_size / 4),
      hjust = axis_just
    )
  }
  
  theme_grey(base_size = font_size, base_family = font_family) %+replace%
    theme(
      rect              = element_rect(fill = "transparent", colour = NA, color = NA, size = 0, linetype = 0),
      text              = element_text(family = font_family, face = "plain", colour = "black",
                                       size = font_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = .9,
                                       margin = ggplot2::margin(), debug = FALSE),
      axis.text         = element_text(colour = "black", size = small_size),
      #axis.title        = element_text(face = "bold"),
      axis.text.x       = element_text(margin = ggplot2::margin(t = small_size / 4), vjust = 1),
      axis.text.y       =  element_text(margin = ggplot2::margin(r = small_size / 4), hjust = .5, vjust = 0, angle = 90),
      axis.title.x      = element_text(
        margin = ggplot2::margin(t = small_size / 2, b = small_size / 4),
        hjust = axis_just
      ),
      axis.title.y      = axis.title.y,
      axis.ticks        = axis.ticks,
      axis.ticks.y      = axis.ticks.y,
      axis.line         = element_blank(),
      legend.key        = element_blank(),
      legend.key.size   = grid::unit(1, "lines"),
      legend.text       = element_text(size = rel(small_rel)),
      legend.justification = c("left", "center"),
      panel.background  = element_blank(),
      panel.border      = element_blank(),
      # make grid lines
      panel.grid.major  = panel.grid.major,
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(size = rel(small_rel)),
      strip.background  = element_rect(fill = "grey80", colour = "grey50", size = 0),
      plot.background   = element_blank(),
      plot.title        = element_text(face = "bold",
                                       size = font_size,
                                       margin = ggplot2::margin(b = half_line), hjust = 0),
      plot.subtitle     = element_text(size = rel(small_rel),
                                       hjust = 0, vjust = 1,
                                       margin = margin(b = half_line * small_rel)),
      plot.caption      = element_text(size = rel(small_rel),
                                       hjust = 1, vjust = 1,
                                       margin = margin(t = half_line * small_rel)),
      plot.margin       = ggplot2::margin(half_line, font_size, half_line, half_line),
      
      complete = TRUE
    )
}

#swap out sticky and mass labels if you want to make the mass or sticky figures
Sticky<-c(71.08,156.2,114.1,115.1,103.1,129.1,128.1,57.05,137.1,113.2,113.2,128.2,131.2,147.2,97.12,87.08,101.1,186.2,163.2,99.13)
#Sticky<-c(0.0062, 1.0372, -0.7485, -0.7893, 1.2727, -0.1771, 0.1204, 1.1109, -1.1806, 0.9138, 1.0124, -0.2693, -0.1799, -0.4114, -0.0876, 0.1376, 0.1031, 0.7599, 0.7925, 0.8806)

lab="Mass"
#lab="Stickiness"
theme_tA<-theme_ridges_AT() 
theme_t<-theme_ridges_AT(ylab_on=FALSE)
#scale for mass
scaler_var<-1.75
#band for mass
band<-10

#scale for sticky
#scaler_var<-1.75
#band<-0.25

justleg<-ggplot(dfB, aes(x = dfB$Sticky.15, y = dfB$Gen, fill=paste(dfB$Exp))) + geom_density_ridges(rel_min_height = .001,color = "white",panel_scaling=FALSE, bandwidth=band, alpha=.5,scale=scaler_var,from=min(Sticky), to=max(Sticky))+ylab("Substitutions")+xlab(lab)

B15<-ggplot(dfB, aes(x = dfB$Sticky.15, y = dfB$Gen, fill=paste(dfB$Exp,dfB$Sel.15))) + geom_density_ridges(rel_min_height = .001,color = "white",panel_scaling=FALSE, bandwidth=band, alpha=.5, scale=scaler_var,from=min(Sticky), to=max(Sticky))+ylab("Substitutions")+xlab(lab)

theme_tA<-theme_ridges_AT() 
theme_t<-theme_ridges_AT(ylab_on=FALSE)

dis<-scale_y_discrete(breaks=c(seq(0,2000,200)), expand = c(0.01, 0))

dis2<-scale_x_continuous(breaks=c(round(min(Sticky)),round(mean(Sticky)),round(max(Sticky))),expand = c(0.01, 0))


breaks<-scale_fill_cyclical(breaks = c("Bind Both 0","Bind Both 1", "Bind A and Not A' 0","Bind A and Not A' 1"),
                            values = c(  "#B2D4FF","#0072B2" ,"#FFDDD5","#D55E00"))

breaks_leg<-scale_fill_cyclical(labels=c("bind both", "bind A and not A'"),breaks = c("Bind Both","Bind A and Not A'"),values = c("#0072B2","#D55E00"),name = "Simulation",guide = "legend")


B15<-B15+dis+breaks+dis2
B15<-B15+theme_tA


justleg<-justleg+breaks_leg+dis2 +guides(fill=guide_legend(nrow=1,byrow=TRUE))+ theme(legend.title.align=0.5) 



B16<-ggplot(dfB, aes(x = dfB$Sticky.16, y = dfB$Gen, fill=paste(dfB$Exp,dfB$Sel.16))) + geom_density_ridges(rel_min_height = .001,color = "white",panel_scaling=FALSE, bandwidth=band, alpha=.5, scale=scaler_var,from=min(Sticky), to=max(Sticky))
B16<-B16+theme_t+dis+breaks+xlab(lab)+dis2


B20<-ggplot(dfB, aes(x = dfB$Sticky.20, y = dfB$Gen, fill=paste(dfB$Exp,dfB$Sel.20))) + geom_density_ridges(rel_min_height = .001,color = "white",panel_scaling=FALSE, bandwidth=band, alpha=.5, scale=scaler_var,from=min(Sticky), to=max(Sticky))
B20<-B20+theme_t+dis+breaks+xlab(lab)+dis2


C61<-ggplot(df, aes(x = df$Sticky.61, y = df$Gen, fill=paste(df$Exp,df$Sel.61))) + geom_density_ridges(rel_min_height = .001,color = "white",panel_scaling=FALSE, bandwidth=band, alpha=.5,scale=scaler_var,from=min(Sticky), to=max(Sticky))

C61<-C61+theme_t+dis+breaks+xlab(lab)+dis2

C64<-ggplot(df, aes(x = df$Sticky.64, y = df$Gen, fill=paste(df$Exp,df$Sel.64))) + geom_density_ridges(rel_min_height = .001,color = "white",panel_scaling=FALSE, bandwidth=band, alpha=.5,scale=scaler_var,from=min(Sticky), to=max(Sticky))

C64<-C64+theme_t+dis+breaks+xlab(lab)+dis2

C70<-ggplot(df, aes(x = df$Sticky.70, y = df$Gen, fill=paste(df$Exp,df$Sel.70))) + geom_density_ridges(rel_min_height = .001,color = "white",panel_scaling=FALSE, bandwidth=band, alpha=.5,scale=scaler_var,from=min(Sticky), to=max(Sticky))
C70<-C70+theme_t+dis+breaks+xlab(lab)+dis2






legend<-get_legend(justleg)


pp<- plot_grid( B15+ theme(legend.position="none"),
                B16,
                B20,
                C61,
                C64,
                C70,
                align = 'vh',
                labels = c("A", "B","C","D","E","F"),
                hjust = -1,
                nrow = 1)

p<-plot_grid(pp,legend, rel_heights = c(3, .2),nrow=2)
print(p)

ggplot2::ggsave("C:/Users/User/Documents/dup_alignments/Mass_All_dup_A_grid_nice.pdf",plot=p, width = 11, height = 9,units = "in")