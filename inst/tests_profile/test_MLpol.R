rm(list=objects())
library(opera)
library(magrittr)


####simulation  une gaussienne avec des sauts dans la moyenne et  des prédicteurs constants
#mu, correspond à l'espérance de chacun des régimes
#n, nb de donnnées par période, la taille totale de la série est donc  de length(mu)*n
#M, nb d'experts
#sd, standard deviation  of the error
set.seed(3)

simu_data <- function(n, M, sd=0.5, mu=c(-0.5,0,0.5,-0.5,0,0.5), expert_range=c(-1,1))
{
  y <- lapply(mu, rnorm, n=n, sd=sd)%>%unlist()
  experts <- lapply(seq(expert_range[1], expert_range[2], length=M), 
                    rep, times=length(mu)*n)%>%unlist()%>%matrix(, byrow=F, ncol=M)
  l <- list(y=y, experts=experts)
}

computeloss <- function(n,M,mod,lt,lg)
{
  set.seed(3)
  nsample<-5
  elapsed<-NA
  for (i in 1:nsample) {
    sim  <- simu_data(n=n, M=M, sd=0.5, mu=c(-0.5,0,0.5,-0.5,0,0.5), expert_range=c(-1,1))
    start_time <- Sys.time()
    m <- mixture(Y=sim$y, experts=sim$experts, model=mod,loss.type=lt,loss.gradient = lg)
    end_time <- Sys.time()
    if (i==1){
      loss<-m$loss
    }
    current_elapsed=as.double(end_time-start_time)
    elapsed=min(current_elapsed,elapsed,na.rm=TRUE)
    if (current_elapsed>0.2 && i>1){
      break
    }
    elapsed=min(as.double(end_time-start_time),elapsed,na.rm=TRUE)
  }
  
  return(list(loss,elapsed))
}



ns<-rbind(50)
#ns<-rbind(50)
#Ms<-rbind(2,4,10,100,1000)
Ms<-rbind(10)
#Ms<-rbind(5)
#mods <- c("Ridge","EWA","MLpol","MLprod","BOA")
mods <- c("MLpol","MLprod","BOA","EWA","Ridge")
ls <- c("absolute", "square", "pinball", "percentage")
#ls <- c("square")
gs <- c(TRUE,FALSE)
#gs <- c(TRUE)



l_profis <- list()
nm<-length(ns)*length(Ms)*length(ls)*length(gs)*length(mods)
tm<-data.frame(n=rep(NA,nm),
               M=rep(NA,nm),
               Meth=rep("",nm),
               LossF=rep("",nm),
               Gr=rep(TRUE,nm),
               CLoss=rep(NA,nm),
               RLoss=rep(NA,nm),
               Rdiff=rep(NA,nm),
               Ctime=rep(NA,nm),
               Rtime=rep(NA,nm),
               SpUp=rep(NA,nm),
               stringsAsFactors = FALSE)

i<-1
for (n in ns) {
  for (M in Ms) {
    for (mod in mods){
      for (l in ls){
        for (g in gs){
          
          cat("computation ",i," n=",n,"M=",M,"mod=",mod,"l=",l,"gs=", g,"\n")
          
          if ((mod=="Ridge") && (l!="square")){#Ridge only allows for square loss
            break;
          }
          
          options("opera_use_cpp" = F)
          {
            result<-computeloss(n,M,mod,l,g)
            rloss<-result[[1]]
            relapsed<-result[[2]]
          }
          options("opera_use_cpp" = T)
          {
            result<-computeloss(n,M,mod,l,g)
            closs<-result[[1]]
            celapsed<-result[[2]]
          }
          
          lossdiff<-(closs-rloss)/max(closs,rloss)
          spup=relapsed/celapsed  
          tm[i,]<-list(n,M,mod,l,g,closs,rloss,lossdiff,celapsed,relapsed,spup)
          i<-i+1
        }
      }
    }
  }
}

tm<-tm[1:i-1,]

print(tm)







