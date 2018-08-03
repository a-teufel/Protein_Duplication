# script to figure out which sites differ significantly based on a chi-square test
# WT refers to simulations that bind both, and no bind is the bind B and not B' simuations
# For the simulations of proteins of different lengths the internal values of this were changed to reflect those proteins lengths
# protein subunit A = 77 sites, 
# protein subunit B = 74 sites,
# protein subunit C (B') is the same as B
# Note, this script will print a bunch of troubleshooting stuff while running. 
# The inital version of this was written by M. Johnson and was later modified by A. Teufel

rm(list=ls())
library(tidyverse)
library(stringr)
library(cowplot)
library(broom)
library(dplyr)
library(grid)
library(ggridges)
library(colorspace)



makeA_sites<-function(site_test){
  site_vector<-NULL
  
  for(i in 1:77){
    if(i %in% site_test){
      site_vector<-c(site_vector,1)
    }
    else{
      site_vector<-c(site_vector,0)
    }
  }
  
  return(site_vector)
}

makeB_sites<-function(site_test){
  site_vector<-NULL
  for(i in 1:74){
    if(i %in% site_test){
      site_vector<-c(site_vector,1)
    }
    else{
      site_vector<-c(site_vector,0)
    }
  }
  
  return(site_vector)
}




# protein subunit A = 154 sites, 100 reps
# protein subunit B = 73 sites, 100 reps
sites_timeA<-data.frame()
sites_timeB<-data.frame()
sites_timeC<-data.frame()


timeA<-NULL
timeB<-NULL
timeC<-NULL


for(i in seq(1,2000,50)){
  # Import & tidy wildtype sequencesseq
  name_WT<-paste("C:/Users/User/Documents/dup_alignments/M_SIM1_4gvb_csv_files/M_SIM1_4gvb_",i,".csv",sep="")  
  
  # Import & tidy wildtype sequencesseq
  wt_sequences = read_csv( name_WT, col_names = c("simulation", "proteinA", "proteinB1", "proteinB2", "proteinC1","proteinC2","empty"))
  proteinB<-paste(wt_sequences$proteinB1,wt_sequences$proteinB2,sep="")
  proteinC<-paste(wt_sequences$proteinC1,wt_sequences$proteinC2,sep="")
  wt_sequences$proteinB=proteinB
  wt_sequences$proteinC=proteinC
  
  
  
  wt_sequences %>% select(-empty) -> wt_sequences
  
  wt_sequences %>% mutate(replicate = as.numeric(str_extract(wt_sequences$simulation, "[0-9]+"))) -> wt_sequences
  
  wt_sequences %>% 
    mutate(a = map(proteinA, 
                   ~ data.frame(aa = strsplit(., "")[[1]], 
                                site = 1:nchar(.), 
                                stringsAsFactors = FALSE)),
           b = map(proteinB, 
                   ~ data.frame(aa = strsplit(., "")[[1]], 
                                site = 1:nchar(.), 
                                stringsAsFactors = FALSE)),
           c = map(proteinC, 
                   ~ data.frame(aa = strsplit(., "")[[1]], 
                                site = 1:nchar(.), 
                                stringsAsFactors = FALSE)))%>%
    
    select(replicate, a, b, c) %>%
    gather(protein, aa_data, -replicate) %>%
    unnest() -> tidy_wt
  
  print(wt_sequences)
  # Import & tidy nobind sequences
  name_EXP<-paste("C:/Users/User/Documents/dup_alignments/M_SIM2_4gvb_csv_files/M_SIM2_4gvb_",i,".csv",sep="")  
  
  nobind_sequences = read_csv( name_EXP, col_names = c("simulation", "proteinA", "proteinB1", "proteinB2", "proteinC1","proteinC2","empty"))
  proteinB<-paste(nobind_sequences$proteinB1,nobind_sequences$proteinB2,sep="")
  proteinC<-paste(nobind_sequences$proteinC1,nobind_sequences$proteinC2,sep="")
  nobind_sequences$proteinB=proteinB
  nobind_sequences$proteinC=proteinC
  
  nobind_sequences %>% select(-empty) -> nobind_sequences
  
  nobind_sequences %>% mutate(replicate = as.numeric(str_extract(nobind_sequences$simulation, "[0-9]+"))) -> nobind_sequences
  
  nobind_sequences %>% 
    mutate(a = map(proteinA, 
                   ~ data.frame(aa = strsplit(., "")[[1]], 
                                site = 1:nchar(.), 
                                stringsAsFactors = FALSE)),
           b = map(proteinB, 
                   ~ data.frame(aa = strsplit(., "")[[1]], 
                                site = 1:nchar(.), 
                                stringsAsFactors = FALSE)),
           c = map(proteinC, 
                   ~ data.frame(aa = strsplit(., "")[[1]], 
                                site = 1:nchar(.), 
                                stringsAsFactors = FALSE)))%>%
    
    select(replicate, a, b,c) %>%
    gather(protein, aa_data, -replicate) %>%
    unnest() -> tidy_nobind
  print(nobind_sequences)
  
  # Perform a chi-squared test of homogeneity on wildtype and nobind proteins
  # 1. Function that can take a data frame with amino acids at one site, and return a data frame with counts, including zeros.
  empty_counts <- data_frame(aa = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T',
                                    'V', 'W', 'Y'), 
                             count = rep(0, 20))
  
  count_by_site <- function(site_df) {
    site_df %>%
      group_by(aa) %>%
      summarize(count = n()) -> counts
    suppressMessages(full_join(counts, anti_join(empty_counts, counts, by = 'aa')))
  }
  
  
  # 2. Call that function for all sites and combine the resulting data frames.
  # wildtype protein A
  tidy_wt %>% 
    filter(protein == "a") %>%
    nest(-site) %>%
    mutate(counts = map(data, count_by_site)) %>%
    select(-data) %>%
    unnest() -> dist_wt_a
  colnames(dist_wt_a) <- c("site", "aa", "countWT")
  print(dist_wt_a)
  
  
  # wildtype protein B
  tidy_wt %>% 
    filter(protein == "b") %>%
    nest(-site) %>%
    mutate(counts = map(data, count_by_site)) %>%
    select(-data) %>%
    unnest() -> dist_wt_b
  colnames(dist_wt_b) <- c("site", "aa", "countWT")
  
  # wildtype protein C
  tidy_wt %>% 
    filter(protein == "c") %>%
    nest(-site) %>%
    mutate(counts = map(data, count_by_site)) %>%
    select(-data) %>%
    unnest() -> dist_wt_c
  colnames(dist_wt_c) <- c("site", "aa", "countWT")
  #head(dist_wt_b)
  
  # nobind protein A
  tidy_nobind %>% 
    filter(protein == "a") %>%
    nest(-site) %>%
    mutate(counts = map(data, count_by_site)) %>%
    select(-data) %>%
    unnest() -> dist_nobind_a
  colnames(dist_nobind_a) <- c("site", "aa", "countNB")
  #head(dist_nobind_a)
  
  # nobind protein B
  tidy_nobind %>% 
    filter(protein == "b") %>%
    nest(-site) %>%
    mutate(counts = map(data, count_by_site)) %>%
    select(-data) %>%
    unnest() -> dist_nobind_b
  colnames(dist_nobind_b) <- c("site", "aa", "countNB")
  #head(dist_nobind_b)
  
  
  # nobind protein C
  tidy_nobind %>% 
    filter(protein == "c") %>%
    nest(-site) %>%
    mutate(counts = map(data, count_by_site)) %>%
    select(-data) %>%
    unnest() -> dist_nobind_c
  colnames(dist_nobind_c) <- c("site", "aa", "countNB")
  #head(dist_nobind_b)
  
  # 3. Perform a chi-squared test at each site - are amino acids distributed equally under wildtype and nobind conditions?
  # wildtype vs nobind protein A
  protein_a <- left_join(dist_wt_a, dist_nobind_a) # create one data frame
  
  # chi square test modifications:
  # instead of x and y, cbind(wt and nb) to make contingency table - cbind = vectors as columns
  # add one to all counts, bc chisq.test() ignores counts that are zero
  # adjust the p-values to account for false discovery
  # the function automatically does a Yates correction to account for small count values
  protein_a %>%
    group_by(site) %>%
    do(tidy(chisq.test(cbind(.$countWT + 1, .$countNB + 1)))) %>%
    ungroup() %>%
    mutate(p.adjusted = p.adjust(p.value,method= "fdr")) -> chisquare_a
  # signficant sites: 16 and 20
  print(chisquare_a)
  print(which(chisquare_a$p.adjusted < 0.05))
  if(length(which(chisquare_a$p.adjusted < 0.05)) == 0){
    #df<-append(df,as.data.frame(site_vector))
    #sites_timeA<-c(sites_timeA,0)
    site_vector<-c(rep(0,77))
    sites_timeA<-append(sites_timeA,as.data.frame(site_vector))
  }
  else{
    timeA<-c(timeA, i)
    site_vector<-makeA_sites(which(chisquare_a$p.adjusted < 0.05))
    sites_timeA<-append(sites_timeA,as.data.frame(site_vector))
  }
  
  # wildtype vs nobind protein B
  protein_b <- left_join(dist_wt_b, dist_nobind_b)  # create one data frame
  
  protein_b %>%
    group_by(site) %>%
    do(tidy(chisq.test(cbind(.$countWT + 1, .$countNB + 1)))) %>%
    ungroup() %>%
    mutate(p.adjusted = p.adjust(p.value,method= "fdr")) -> chisquare_b
  
  print(which(chisquare_b$p.adjusted < 0.05))
  
  if(length(which(chisquare_b$p.adjusted < 0.05)) == 0){
    site_vector<-c(rep(0,74))
    sites_timeB<-append(sites_timeB,as.data.frame(site_vector))
  }
  else{
    timeB<-c(timeB, i)
    site_vector<-makeB_sites(which(chisquare_b$p.adjusted < 0.05))
    sites_timeB<-append(sites_timeB,as.data.frame(site_vector))
  }
  protein_c <- left_join(dist_wt_c, dist_nobind_c)  # create one data frame
  
  protein_c %>%
    group_by(site) %>%
    do(tidy(chisq.test(cbind(.$countWT + 1, .$countNB + 1)))) %>%
    ungroup() %>%
    mutate(p.adjusted = p.adjust(p.value,method= "fdr")) -> chisquare_c
  
  print(which(chisquare_c$p.adjusted < 0.05))
  if(length(which(chisquare_c$p.adjusted < 0.05)) == 0){
    site_vector<-c(rep(0,74))
    sites_timeC<-append(sites_timeC,as.data.frame(site_vector))
  }
  else{
    timeC<-c(timeC, i)
    site_vector<-makeB_sites(which(chisquare_c$p.adjusted < 0.05))
    sites_timeC<-append(sites_timeC,as.data.frame(site_vector))
  }
  
  
}


print(sites_timeA)
print(sites_timeB)
print(sites_timeC)

# writes a file where the rows are the sites and the columns are the generations mesured.
write.csv(file="C:/Users/User/Documents/dup_alignments/diffs_A_other.csv", sites_timeA,row.names=FALSE)
write.csv(file="C:/Users/User/Documents/dup_alignments/diffs_B_other.csv", sites_timeB,row.names=FALSE)
write.csv(file="C:/Users/User/Documents/dup_alignments/diffs_C_other.csv", sites_timeC,row.names=FALSE)