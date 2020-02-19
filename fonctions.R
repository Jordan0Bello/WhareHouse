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



##################################################################################################################################################
# **************************************************   Graphs   **********************************************************************************
##################################################################################################################################################

Distrib<-function(x,Vble=Vble){
  hist(x[,Vble],col="#8bbaea",xlab="",ylab="effectifs",main=Vble,prob=TRUE)
  lines(density(x[,Vble],na.rm=TRUE),col="red",lwd=2)
}

##
DistribNorm<-function(x,Vble=Vble){
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
PlotHAP <-function(x,HAP=HAP,Enzyme=Enzyme,na.rm=T){
  plot(x[,HAP],x[,Enzyme],
       main="Activite Enzymatique",
       ylab=Enzyme, xlab=HAP)}




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


###############################################################################################
#
# Test Critere de Peirce    ###################################################################
#
###############################################################################################
findx <- function(N,k,m){ 
  # method by K. Thomsen (2008)
  # written by C. Dardis and S. Muller (2012)
  # Available online: https://r-forge.r-project.org/R/?group_id=1473
  #
  # Variable definitions:
  # N :: number of observations
  # k :: number of potential outliers to be removed
  # m :: number of unknown quantities
  #
  # Requires the complementary error function, erfc:
  erfc <- function(x) 2 * pnorm(x * sqrt(2), lower.tail = FALSE)
  #
  x <- 1
  if ((N - m - k) <= 0) {
    return(NaN)
    print(NaN)
  }  else {
    x    <- min(x, sqrt((N - m)/k) - 1e-10)
    #
    # Log of Gould's equation B:
    LnQN <- k * log(k) + (N - k) * log(N - k) - N * log(N)
    #
    # Gould's equation D:
    R1   <- exp((x^2 - 1)/2) * erfc(x/sqrt(2))
    #
    # Gould's equation A' solved for R w/ Lambda substitution:
    R2   <- exp( (LnQN - 0.5 * (N - k) * log((N-m-k*x^2)/(N-m-k)) )/k )
    #
    # Equate the two R equations:
    R1d  <- x * R1 - sqrt(2/pi/exp(1))
    R2d  <- x * (N - k)/(N - m - k * x^2) * R2
    #
    # Update x:
    oldx <- x
    x    <- oldx - (R1 - R2)/(R1d - R2d)
    #
    # Loop until convergence:
    while (abs(x - oldx) >= N * 2e-16){
      R1   <- exp((x^2 - 1)/2) * erfc(x/sqrt(2))
      R2   <- exp( (LnQN - 0.5 * (N - k) * log((N-m-k*x^2)/(N-m-k)) )/k )
      R1d  <- x * R1 - sqrt(2/pi/exp(1))
      R2d  <- x * (N - k)/(N - m - k * x^2) * R2
      oldx <- x
      x    <- oldx - (R1 - R2)/(R1d - R2d)
    }
  }
  return(x)
}


##
#
# 2e partie
#
##

peirce_dev <- function(N, n, m){
  # N :: total number of observations
  # n :: number of outliers to be removed
  # m :: number of model unknowns (e.g., regression parameters)
  #
  # Check number of observations:
  if (N > 1) {
    # Calculate Q (Nth root of Gould's equation B):
    Q = (n^(n/N) * (N-n)^((N-n)/N))/N
    #
    # Initialize R values:
    Rnew = 1.0
    Rold = 0.0  # <- Necessary to prompt while loop
    #
    while(abs(Rnew-Rold) > (N*2.0e-16)){
      # Calculate Lamda (1/(N-n)th root of Gould's equation A'):
      ldiv = Rnew^n
      if (ldiv == 0){
        ldiv = 1.0e-6
      }
      Lamda = ((Q^N)/(ldiv))^(1.0/(N-n))
      #
      # Calculate x-squared (Gould's equation C):
      x2 = 1.0 + (N-m-n)/n * (1.0-Lamda^2.0)
      #
      # If x2 goes negative, set equal to zero:
      if (x2 < 0){
        x2 = 0
        Rold = Rnew
      } else {
        #
        # Use x-squared to update R (Gould's equation D):
        # NOTE: error function (erfc) is replaced with pnorm (Rbasic):
        # source: 
        # http://stat.ethz.ch/R-manual/R-patched/library/stats/html/Normal.html
        Rold = Rnew
        Rnew = exp((x2-1)/2.0)*(2*pnorm(sqrt(x2)/sqrt(2)*sqrt(2), lower=FALSE))
      }
    }
  } else {
    x2 = 0
  }
  x2
}


###############################################################################################
#
# Diagramme gradient indicateur pH   ##########################################################
#
###############################################################################################

# x.shade<-seq(0,14,0.01)
# y.gradient<-rep(1,length(x.shade))
# Gradient<-as.data.frame(cbind(x.shade,y.gradient))


gradient<-function(data=data, value=value){
p<-ggplot(data) +
  geom_tile(aes(x = seq(0,14,0.01), y=1, fill = seq(0,14,0.01))) +
  scale_y_continuous(breaks=NULL,name = NULL)+
  scale_x_continuous(breaks=seq(0,14,1),name = "pH")+
  scale_fill_gradient2(low = 'red', mid = 'green', high = 'red', midpoint = 7) +
  geom_point(aes(x=value,y=1),shape=23, fill="blue",color="black",size=5)+
  labs(fill="Tolerance")+
  theme(panel.background = element_blank())
p
}




