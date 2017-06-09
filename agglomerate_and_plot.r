# Code from Clavel et al. (2014) Syst. Bio. - Supplementary Material S2

## Combination of the results

#Dependencies
library(shapes)
library(mice)
library(Amelia)
library(missMDA)
library(Hmisc)
library(norm)

# agglomerate.data (data,imp,Mimp,Method="mice")
#
# plot.MI(IM, symmetric=TRUE, DIM=c(1,2),web=FALSE,ellipses=TRUE)
#
# data= dataset with missing values
#
# imp=Imputed datasets
#
# Mimp=Number of MI
#
# Method= MI method used, must be "mice","amelia","hmisc","missmda" or “norm”
#
# Julien Clavel - 2013
#################################################
#################################################
## Agglomerate.data. Function to store the m imputed datasets and to compute the averaged dataset

agglomerate.data<-function(data,imp,Mimp,Method="mice"){
  
  Moy<-Mimp+1
  redata<-as.matrix(data)
  ximp<-array(redata,dim=c(nrow(redata),ncol(redata),Moy))
  #####################
  
  if(any(is.na(redata))==TRUE){
    if(Method=="mice" || Method=="amelia" || Method=="missmda" || Method=="hmisc" || Method=="norm"){
      #####################MICE
      if(Method=="mice"){
        
        for(i in 1:Mimp){
          ximp[,,i]<-as.matrix(complete(imp,i))
          
        }
        ##Averaged dataset
        ximp[,,Moy]<-apply(ximp[,,1:Mimp],c(1,2),mean)
      }
      #####################
      #####################Amelia
      if(Method=="amelia"){
        for(i in 1:Mimp){
          ximp[,,i]<-as.matrix(imp$imputations[[i]])
          
        }
        ##Averaged dataset
        ximp[,,Moy]<-apply(ximp[,,1:Mimp],c(1,2),mean)
      }
      #####################
      #####################
      #####################NORM
      if(Method=="norm"){
        for(i in 1:Mimp){
          ximp[,,i]<-as.matrix(imp[[i]])
          
        }
        ##Averaged dataset
        ximp[,,Moy]<-apply(ximp[,,1:Mimp],c(1,2),mean)
      }
      #####################
      #####################MDA
      if(Method=="missmda"){
        
        for(i in 1:Mimp){
          ximp[,,i]<-as.matrix(imp$res.MI[,,i])
          
        }
        ##Averaged dataset
        ximp[,,Moy]<-apply(ximp[,,1:Mimp],c(1,2),mean)
      }
      ####################
      ####################Hmisc
      if(Method=="hmisc"){
        ##Extract the m data imputed for each variables
        ximp<-array(redata, dim=c(nrow(redata),ncol(redata),Moy))
        col<-1:ncol(redata)
        for(j in 1:ncol(redata)){
          if(sum(is.na(redata[,j]))==0){
            col<-col[-which(col==j)]
            next
          }
        }
        for(m in 1:Mimp){
          for(g in col){
            ximp[,,m][!complete.cases(ximp[,,m][,g]),g]<-imp$imputed[[g]][,m]
          }
        }
        ##Averaged dataset
        ximp[,,Moy]<-apply(ximp[,,1:Mimp],c(1,2),mean)
        
      }
    }else{
      ####################
      ##Warning messages
      cat("Error! You must indicate if you are using Mice, Amelia, missMDA, NORM, or Hmisc package","\n")
    }}else{
      ## Warning messages
      cat("There is no missing value in your dataset","\n")
    }
  #return(ximp)
  tabM<-ximp[,,Moy]
  colnames(tabM)<-colnames(redata)
  list("ImpM"=tabM,"Mi"=ximp[,,1:Mimp],"nbMI"=Mimp, "missing"=as.data.frame(redata))
}#End
###################################
########################################
###############################################
######################################################
##Function to draw confidence ellipses
ELLI<-function(x,y,conf=0.95,np)
{centroid<-apply(cbind(x,y),2,mean)
ang <- seq(0,2*pi,length=np)
z<-cbind(cos(ang),sin(ang))
radiuscoef<-qnorm((1-conf)/2, lower.tail=F)
vcvxy<-var(cbind(x,y))
r<-cor(x,y)
M1<-matrix(c(1,1,-1,1),2,2)
M2<-matrix(c(var(x), var(y)),2,2)
M3<-matrix(c(1+r, 1-r),2,2, byrow=T)
ellpar<-M1*sqrt(M2*M3/2)
t(centroid + radiuscoef * ellpar %*% t(z))}
############################### Function to Plot MI confidence ellipses using procruste superimposition
##################################
####################################
######################################

plot.MI<-function(IM,symmetric=FALSE,DIM=c(1,2),scale=FALSE,web=FALSE,ellipses=TRUE,...){
  if(any(is.na(IM$ImpM)==TRUE))
  { cat("There is still missing values in the imputed dataset, please check your imputation")
    break
  }else{
    Mo<-IM$nbMI+1
    pcaM<-princomp(IM$ImpM)
    cpdimM<-as.matrix(pcaM$scores[,DIM])   
    opa<-array(cpdimM,dim=c(nrow(cpdimM),ncol(cpdimM),Mo)) 
    for(i in 1:IM$nbMI){
      pca<-princomp(IM$Mi[,,i])
      opa[,,i]<-as.matrix(pca$scores[,DIM])
    }
    if(symmetric==TRUE){
      for (i in 1:IM$nbMI+1){
        trace<-sum(opa[,,i]^2)
        opa[,,i]<-opa[,,i]/sqrt(trace) 
      }
    }
    ############################ Ordinary Procrustes Analysis (library(shapes))
    for(k in 1:IM$nbMI){
      analyse<-procOPA(opa[,,Mo],opa[,,k], reflect=TRUE)
      opa[,,k]<-analyse$Bhat
    }
    opa[,,Mo]<-analyse$Ahat
    ######################## Principal component explained variance
    pvar<-pcaM$sdev^2
    tot<-sum(pvar)
    valX<-pvar[DIM[1]]
    valY<-pvar[DIM[2]]
    valX<-round(valX*100/tot,digits=2)
    valY<-round(valY*100/tot, digits=2)
    ######################## Plot function
    op <- par(no.readonly=TRUE)
    if(scale==TRUE){
      plot(opa[,1,Mo],opa[,2,Mo], type="p", pch=3, col=c(as.factor(ifelse(complete.cases(IM$missing) ==T, 1, 5))),lwd=1,xlim=range(opa[,1,Mo]),ylim=range(opa[,1,Mo]),xlab=paste("DIM",DIM[1],valX,"%",sep=" "),ylab=paste("DIM",DIM[2],valY,"%",sep=" "))
    }
    if(scale==FALSE){
      plot(opa[,1,Mo],opa[,2,Mo], type="p", pch=3, col=c(as.factor(ifelse(complete.cases(IM$missing) ==T, 1, 5))),lwd=1,xlab=paste("DIM",DIM[1],valX,"%",sep=" "),ylab=paste("DIM",DIM[2],valY,"%",sep=" "))
    }
    title("MI effect on Multivariate Analysis", font.main=3, adj=1)
    ## Store row names
    NR<-IM$missing
    rownames(IM$missing)<-NULL
    ##
    if(ellipses==TRUE){                                      
      coul<-as.numeric(rownames(IM$missing[complete.cases(IM$missing),]))
      for (j in coul){
        lines(ELLI(opa[j,1,],opa[j,2,],np=Mo), col="black", lwd=1)}
      coul<-as.numeric(rownames(IM$missing[!complete.cases(IM$missing),]))
      for (j in coul){
        lines(ELLI(opa[j,1,],opa[j,2,],np=Mo), col="red", lwd=1)}
    }else{ points(opa[,1,],opa[,2,],cex=0.5) }
    if(web==TRUE){
      coul<-as.numeric(rownames(IM$missing[complete.cases(IM$missing),]))
      for (j in coul){ 
        for(f in 1:IM$nbMI){
          segments(opa[j,1,Mo],opa[j,2,Mo], opa[j,1,f],opa[j,2,f], col="black", lwd=1) }
      } 
      coul<-as.numeric(rownames(IM$missing[!complete.cases(IM$missing),]))
      for (j in coul){ 
        for(f in 1:IM$nbMI){
          segments(opa[j,1,Mo],opa[j,2,Mo], opa[j,1,f],opa[j,2,f], col="red", lwd=1)}
      } 
      points(opa[,1,],opa[,2,],cex=0.5) 
    } 
    nom<-rownames(NR)
    text(opa[,1,Mo],opa[,2,Mo],nom, pos=1)
    
    abline(h=0,v=0, lty=3)
    par(xpd=TRUE)  # Do not clip to the drawing area
    lambda <- .025
    legend(par("usr")[1], (1 + lambda) * par("usr")[4] - lambda * par("usr")[3],c("Complete", "Missing"), xjust = 0, yjust = 0,lwd=3, lty=1, col=c(par('fg'), 'red'))
    par(op)      
  }
}
