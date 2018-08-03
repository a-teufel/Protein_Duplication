#figures out different functional types based on ability to bind (have to change the ancestral binding and survival cut off for the other protein and if you want to do B' or A) 
#this script is needed to make the mosaic plots
get.data <- function(this.folder) {
  start <- this.folder
  dirs <- list.files(start)
  
  ab.count <- c()
  ev.binding <- c()
  an.binding <- c()
  name <- c()
  replicate <- c()
  count <- 0
  survival=c(-5.34086620541)
  #survival=c(-3.9399600698172623) # if using other protien
  for(i in dirs) {
    dat <- read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
    
    count <- count + 1
    
    
    ab.count <- append(ab.count, dat$count)
    
    ev.binding <-append(ev.binding, dat$evolved_interaction_AB <= survival)
    
    an.binding<-append(an.binding, dat$ancestral_interaction_BA <= survival)
    
    name <- append(name, dat$name)
    
    replicate <-append(replicate, rep(count,length(dat$name)))
    
  }
  
  
  tmp.data <- data.frame(rep=replicate, names=name, ab.count=ab.count, ev.binding=ev.binding, an.binding=an.binding,dirs=length(dirs))
  return(tmp.data)
  
}

#read in data frile
Un2 <- get.data(paste('/home/ateufel/Marcott_Sim/beta1_data/M_SIM2_beta1/',sep=""))


Change_evo<-NULL
Change_anc<-NULL
Change_iso<-NULL


#figure out when the changes happen
for(i in 1:Un2$dirs){
  rep<-Un2[Un2$rep==i,]
  rep_lag<-rep[1, ]
  Change_evo<-c(Change_evo,rep$ev.binding-rep_lag$ev.binding)
  Change_anc<-c(Change_anc,rep$an.binding-rep_lag$an.binding)
  Change_iso<-c(Change_iso,rep$ev.binding-rep$an.binding)
}

#Add the this marks of change to the data frame
Un2$Evo_Change<-Change_evo
Un2$Anc_Change<-Change_anc
Un2$Iso_Change<-Change_iso


Defunc_cnt<-NULL
Refunc_cnt<-NULL
Isofunc_cnt<-NULL
Retfunc_cnt<-NULL

Iso_refunc_cnt<-NULL
Iso_some_refunc_cnt<-NULL

Functional_reaquire<-NULL


#loop to dected different conditions of functionalization
for(i in 1:Un2$dirs){
  
  rep<-Un2[Un2$rep==i,]
  
  
  Defunc_cnt<-c(Defunc_cnt,any(rep$ev.binding==FALSE))
  
  Isofunc_cnt<-c(Isofunc_cnt,any(rep$an.binding==FALSE) & all(rep$ev.binding==TRUE))
  
  Retfunc_cnt<-c(Retfunc_cnt,all(rep$ev.binding==TRUE) & all(rep$an.binding==TRUE))
  
  
  Lose_func<-which(rep$Evo_Change==-1)
  
  Gain_func<-which(rep$Evo_Change==0)
  
  
  if( (max(Gain_func)>max(Lose_func)) & length(Lose_func)!=0 ){
    
    Refunc_cnt<-c(Refunc_cnt,TRUE)
    
    print("Function Regained!")
    
    After_gain<-tail(rep, -(max(Lose_func)+1))
    Iso_refunc_cnt<-c(Iso_refunc_cnt,all(After_gain$an.binding==FALSE) & all(After_gain$ev.binding==TRUE))
    
    Iso_some_refunc_cnt<-c(Iso_some_refunc_cnt,any(After_gain$an.binding==FALSE) & all(After_gain$ev.binding==TRUE))
    
    Functional_reaquire<-c(Functional_reaquire,all(After_gain$an.binding==TRUE) & all(After_gain$ev.binding==TRUE))
  }
  
  else{
    Refunc_cnt<-c(Refunc_cnt,FALSE)
  }
}



print(sum(Defunc_cnt))
print(sum(Refunc_cnt))
print(sum(Isofunc_cnt))
print(sum(Retfunc_cnt))
print(sum(Iso_refunc_cnt))
print(sum(Iso_some_refunc_cnt))


print("just Defun")
true_def<-sum(Defunc_cnt)-sum(Refunc_cnt)
print(true_def)
print("true iso")
true_iso<-sum(Isofunc_cnt)
print(true_iso)
print("All Refuncs")
true_refun<-sum(Refunc_cnt)
print(true_refun)
print("Hard Iso refuncs")
true_refun_iso_hard<-sum(Iso_refunc_cnt)
print(true_refun_iso_hard)
print("Soft Iso refuncs")
true_refun_iso_soft<-sum(Iso_some_refunc_cnt)-sum(Iso_refunc_cnt)
print(true_refun_iso_soft)

print("Functional reaquizition")
print(sum(Functional_reaquire))

print("total iso regains")
print(true_refun_iso_hard+true_refun_iso_soft)



print((Defunc_cnt))
print((Refunc_cnt))
print((Isofunc_cnt))
print((Retfunc_cnt))
print((Iso_refunc_cnt))
