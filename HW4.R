install.packages("plm")
library(plm)
KT<-read.csv("/Users/davidch91/Desktop/613/HW4/Koop-Tobias.csv")

#EXERCISE 1

#Randomly selecting wages profile of 5 individuals
KT_wages<-subset(KT,select=c(PERSONID,LOGWAGE))
set.seed(11)
x<-sample(unique(KT_wages$PERSONID),5)
KT_subset <- subset(KT_wages, PERSONID %in% x)

#EXERCISE 2 - RANDOM EFFECTS

#Telling R to work with panel data
KTpanel<-pdata.frame(KT,index=c("PERSONID","TIMETRND"))
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
KT_su100 <- subset(KT, PERSONID %in% s)

#Creating ID_sorting variable
IDsort<-matrix(1,882,1)
as.vector(IDsort)
KT_su100<-cbind(IDsort,KT_su100)
for(i in 2:882){
  KT_su100$IDsort[i]<-ifelse(KT_su100$PERSONID[i]==KT_su100$PERSONID[i-1],KT_su100$IDsort[i-1],KT_su100$IDsort[i-1]+1)
}
Xs<-cbind(KT_su100$EDUC, KT_su100$POTEXPER)

#e<-matrix(0,nrow=100,ncol=100)
#diag(e)<-1

#We want to estimate 100 alphai's (individual fixed effects) for 100 random individuals.

#We start from generating 100 random initial values
set.seed(80)
indivfe<-matrix(runif(100,0.15,3.75),nrow=100,ncol=1)
as.vector(indivfe)

#We also create two betas for EDUC and POTEXPR and merge them with individual FE's
beta<-matrix(0.1,nrow=2,ncol=1)
as.vector(beta)
param<-rbind(indivfe,beta)

#Next, we try to optimize param that maximizes the following likelihood function

#Since there are 882 observations,
MLEpanel<-function(param){
  interi<-matrix(1,nrow=882,ncol=1)
  as.vector(interi)
  for(n in 1:882){
    interi[n]<-replace(interi[n],interi[n]>0,indivfe[(KT_su100$IDsort[n]),])  
  }
  yhat<-interi+Xs%*%param[101:102]
  error<-(KT_su100$LOGWAGE-yhat)/var(KT_su100$LOGWAGE)
  errorsq<-(error)^2
  #The following is the formula for likelihood function
  LL<-(-882/2)*log(var(KT_su100$LOGWAGE))-(882/2*log(2*3.14))-0.5*sum(errorsq)
  return(-LL)
}
optparam<-optim(par=param, fn=MLEpanel)
answerex4a<-optparam$par
est_indiv_fe<-answerex4a[1:100]
est_indiv_fe

#est_indiv_fe above is the estimated 100 individual FE's we are looking for in Exercise 4a.

#Next, we regress this 100 individual FE's on other time-invariant variables
estif<-matrix(1,882,1)
as.vector(estif)
for(n in 1:882){
  estif[n]<-replace(estif[n],estif[n]>0,answerex4a[(KT_su100$IDsort[n]),])  
}
KT_su100n<-cbind(KT_su100,estif)

#We combine all time-invariant variables into Xi
Xi <- cbind(KT_su100n$ABILITY,KT_su100n$MOTHERED,KT_su100n$FATHERED,KT_su100n$BRKNHOME,KT_su100n$SIBLINGS)

#and regress individual fixed effects on Xi
IFE<-glm(KT_su100n$estif ~ Xi)
summary(IFE)
#the following is the standard errors
ses<-summary(IFE)$coefficients[,2]

#Correcting Standard Errors using Bootstrap at Both Levels of Regression

#Using 49 replications bootstrap
#In each replication, we sample with replacement 100 random individuals from our initial sample KT_100su,
#then optimize for individual fixed effects under that particular bootstrap sample, then
#regress it on other time-invariant variables
#and then we store the estimates of Std Errors of this particular regression from this particular bootstrap sample.
#Lastly, we will take average of 49 standard errors estimates.

#Bootstrap calculation took about 2-3 minutes to complete in my computer

set.seed(1650)
nboot<-49
seboot<-matrix(NA,49,6)
for(i in 1:nboot){
  sboot<-sample(unique(KT$PERSONID),100,replace=TRUE)
  KT_su100b <- subset(KT, PERSONID %in% sboot)
  IDsortb<-matrix(1,nrow(KT_su100b),1)
  as.vector(IDsortb)
  KT_su100b<-cbind(IDsortb,KT_su100b)
  for(k in 2:nrow(KT_su100b)){
    KT_su100b$IDsortb[k]<-ifelse(KT_su100b$PERSONID[k]==KT_su100b$PERSONID[k-1],KT_su100b$IDsortb[k-1],KT_su100b$IDsortb[k-1]+1)
  }
  Xsb<-cbind(KT_su100b$EDUC, KT_su100b$POTEXPER)
  indivfeb<-matrix(runif(100,0.45,0.75),nrow=100,ncol=1)
  as.vector(indivfeb)
  betab<-matrix(0.1,nrow=2,ncol=1)
  as.vector(betab)
  paramb<-rbind(indivfeb,betab)
  MLEpanelb<-function(paramb){
    interib<-matrix(1,nrow(KT_su100b),1)
    as.vector(interib)
    for(n in 1:nrow(KT_su100b)){
      interib[n]<-replace(interib[n],interib[n]>0,indivfeb[(KT_su100b$IDsortb[n]),])  
    }
    yhatb<-interib+Xsb%*%paramb[101:102]
    errorb<-(KT_su100b$LOGWAGE-yhatb)/var(KT_su100b$LOGWAGE)
    errorsqb<-(errorb)^2
    LLb<-(-nrow(KT_su100b)/2)*log(var(KT_su100b$LOGWAGE))-(nrow(KT_su100b)/2*log(2*3.14))-0.5*sum(errorsqb)
    return(-LLb)
  }
  optparamb<-optim(par=paramb, fn=MLEpanelb)
  answerex4c<-optparamb$par
  est_indiv_fe_b<-answerex4c[1:100]
  estifb<-matrix(1,nrow(KT_su100b),1)
  as.vector(estifb)
  for(n in 1:nrow(KT_su100b)){
    estifb[n]<-replace(estifb[n],estifb[n]>0,answerex4c[(KT_su100b$IDsortb[n]),])  
  }
  KT_su100nb<-cbind(KT_su100b,estifb)
  Xib <- cbind(KT_su100nb$ABILITY,KT_su100nb$MOTHERED,KT_su100nb$FATHERED,KT_su100nb$BRKNHOME,KT_su100nb$SIBLINGS)
  IFEb<-glm(KT_su100nb$estif ~ Xib)
  stderror<-summary(IFEb)$coefficients[,2]
  seboot[i,] <- stderror
}
sebootstrap <- colMeans(seboot)
sebootstrap