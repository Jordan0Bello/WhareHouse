###############################################################
##  script fonctions a appeler par source()                ####
##                                                         ####
##  Bello Jordan                   le 12/04/2019           ####
##                                                         ####
##  maj le          --/--/----                             ####
##                                                         ####
##                                                         ####
###############################################################


##################################################################################################################################################
# ***********************************   Moyenne et erreur standard etc   *************************************************************************
##################################################################################################################################################

sem<-function(x,digits= 4,na.rm=FALSE){x<-x[!is.na(x)];return(round(sd(x)/sqrt(length(x)),digits))}
sumfun <- function(x, ...){c(mean=mean(x, ...), sem=sem(x, ...),sd=sd(x, ...),max=max(x,...),min=min(x,...))}
funMSE<-function(x, ...){c(mean=mean(x, ...), sem=sem(x, ...))}
funMSD<- function(x, ...){c(mean=mean(x, ...), sd=sd(x, ...))}

library("data.table")

MeanSD <-function(data, VbleRep=VbleRep, VbleEx=VbleEx, Function=Function) 
{
DT <- data.table(data)
DT[,.(mean(VbleRep,na.rm=T),sd(VbleRep,na.rm=T)), by=list(as.character(VbleEx))]}



##################################################################################################################################################
# **************************************************   Graphs   **********************************************************************************
##################################################################################################################################################

Distrib<-function(x,Vble=Vble){
  par(mfrow=c(1,2))
  hist(x[,Vble],col="#8bbaea",xlab="",ylab="effectifs",main=Vble,prob=TRUE)
  lines(density(x[,Vble],na.rm=TRUE),col="red",lwd=2)
  qqnorm(x[,Vble]);qqline(x[,Vble])   
}

##

PlotActEnzy<-function(x,Enzyme=Enzyme,na.rm=T){
  plot(x[,2],x[,Enzyme],
       main=Enzyme,ylab="Activite Enzymatique",
       ylim = c(0,250))}

##

PlotActEnzyOcSOL<-function(x,Abscisse=Abscisse,Enzyme=Enzyme,na.rm=T){
  plot(x[,Abscisse],x[,Enzyme],
       main=Enzyme,ylab="Activite Enzymatique")}

##

PlotOcSOL<-function(x,Abscisse=Abscisse,Enzyme=Enzyme,na.rm=T){
  plot(x[,Abscisse],x[,Enzyme],
       main=Enzyme,ylab="Activite Enzymatique")}

##



##################################################################################################################################################
# **************************************************   Outliers   ********************************************************************************
##################################################################################################################################################
# alpha<-0.01
# x<-tmp
# Variable<-"URE"

TauThompson<-function(x,Variable=Variable,alpha=alpha){
  
  
  xTMP<-x[!(is.na(x[,Variable])),]
  Outliers<-NULL
  # DataClean<-NULL
  
    if (length(xTMP[,Variable])<3){return(c(Variable,"IMPOSSIBLE_NB_ECH"))}
    
    if (alpha<0 || alpha>0.05){return("WARNING : aplha is outside the allowed range")}  
  
  
    if (length(xTMP[,Variable])>2){
  
      n<-as.numeric(length(xTMP[,Variable]))
      t<-as.numeric(unlist(t.test(xTMP[,Variable],conf.level = 1-alpha))["statistic.t"])
   
      Tau<- t*(n-1)/(sqrt(n)*sqrt(n-2+(t*t)))
    
      laMean<-mean(xTMP[,Variable])
      xTMP$sigma<-xTMP[,Variable]-laMean
      Lechantillon<-xTMP[abs(xTMP$sigma)==max(abs(xTMP$sigma)),]
  
          
        while (abs(Lechantillon$sigma)>Tau*sd(xTMP[,Variable]) & length(xTMP[,Variable])>2) {
        
          # return(c(Variable,"YES"))
        
          xTMP<-xTMP[!(abs(xTMP$sigma)==max(abs(xTMP$sigma))),]
          
          
          n<-as.numeric(length(xTMP[,Variable]))
          t<-as.numeric(unlist(t.test(xTMP[,Variable],conf.level = 1-alpha))["statistic.t"])
          
          
          Tau<- t*(n-1)/(sqrt(n)*sqrt(n-2+(t*t)))
          
          laMean<-mean(xTMP[,Variable])
          xTMP$sigma<-xTMP[,Variable]-laMean
          Lechantillon<-xTMP[abs(xTMP$sigma)==max(abs(xTMP$sigma)),]
          Outliers<-rbind(Outliers,Lechantillon)
          
          if (Lechantillon$sigma<Tau*sd(xTMP[,Variable])){
          DataClean<-rbind(DataClean,xTMP)
          
          }
        }
  return(tail(DataClean))  
  # return(data.frame(Outliers))
  # return(c(length(x),length(DataClean), length(Outliers)))
      }
  # return(as.data.frame(DataClean))
}




