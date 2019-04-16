install.packages("plm")
library(plm)
KT<-read.csv("/Users/davidch91/Desktop/613/HW4/Koop-Tobias.csv")

#EXERCISE 1
KT_wages<-subset(KT,select=c(PERSONID,LOGWAGE))
set.seed(11)
x<-sample(unique(KT_wages$PERSONID),5)
KT_subset<-filter(KT_wages,KT_wages$PERSONID %in% x)

#Setting up time variable (t)
t<-matrix(1,nrow=17919,ncol=1)
as.vector(t)
KT<-cbind(KT,t)
for(i in 2:17919){
t[i]<-ifelse(KT$PERSONID[i]==KT$PERSONID[i-1],t[i-1]+1,1)
}

#EXERCISE 2 - RANDOM EFFECTS

#Telling R to work with panel data
KTpanel<-pdata.frame(KT,index=c("PERSONID","t"))
X <- cbind(KTpanel$EDUC, KTpanel$POTEXPER)

randomkt <- plm(LOGWAGE ~ X, data=KTpanel, model="random")
summary(randomkt)

#EXERCISE 3 - FIXED EFFECTS MODEL

#Between Estimator Model
fixedbekt <- plm(LOGWAGE ~ X, data=KTpanel, model="between")
summary(fixedbekt)

#Within Estimator Model
fixedwekt <- plm(LOGWAGE ~ X, data=KTpanel, model="within")
summary(fixedwekt)

#First Difference Model
fixedfdkt <- plm(LOGWAGE ~ X, data=KTpanel, model="fd")
summary(fixedfdkt)

#EXERCISE 4 - UNDERSTANDING FIXED EFFECTS

#Randomly selecting 100 individuals into our regression
set.seed(182)
s<-sample(unique(KT$PERSONID),100)
KT_su100<-filter(KT,KT$PERSONID %in% s)

IDsort<-matrix(1,nrow(882),1)
as.vector(IDsort)
KT_su100<-cbind(IDsort,KT_su100)

for(i in 2:882){
  KT_su100$IDsort[i]<-ifelse(KT_su100$PERSONID[i]==KT_su100$PERSONID[i-1],KT_su100$IDsort[i-1],KT_su100$IDsort[i-1]+1)
}

Xs<-cbind(KT_su100$EDUC, KT_su100$POTEXPER)

e<-matrix(0,nrow=100,ncol=100)
diag(e)<-1
set.seed(80)
indivfe<-matrix(runif(100,0.45,0.75),nrow=100,ncol=1)
as.vector(indivfe)
beta<-matrix(0.1,nrow=2,ncol=1)
as.vector(beta)
param<-rbind(indivfe,beta)

MLEpanel<-function(param){
  interi<-matrix(1,nrow=882,ncol=1)
  as.vector(interi)
  for(n in 1:882){
    interi[n]<-replace(interi[n],interi[n]>0,indivfe[(KT_su100$IDsort[n]),])  
  }
  yhat<-interi+Xs%*%param[101:102]
  error<-(KT_su100$LOGWAGE-yhat)/var(KT_su100$LOGWAGE)
  errorsq<-(error)^2
  LL<-(-882/2)*log(var(KT_su100$LOGWAGE))-(882/2*log(2*3.14))-0.5*sum(errorsq)
  return(-LL)
}
optparam<-optim(par=param, fn=MLEpanel)
answer<-optparam$par
est_indiv_fe<-answer[1:100]

estif<-matrix(1,nrow=882,ncol=1)
as.vector(estif)
for(n in 1:882){
  estif[n]<-replace(estif[n],estif[n]>0,answer[(KT_su100$IDsort[n]),])  
}
KT_su100<-cbind(KT_su100,estif)

X_ti <- cbind(KT_su100$ABILITY,KT_su100$MOTHERED,KT_su100$FATHERED,KT_su100$BRKNHOME,KT_su100$SIBLINGS)
IFE<-glm(KT_su100$estif ~ X_ti)
summary(IFE)