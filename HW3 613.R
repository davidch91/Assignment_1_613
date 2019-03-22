install.packages("bayesm")
library(bayesm)
data(margarine)

demos <- read.csv(file="/Users/davidch91/Desktop/demos.csv", header=TRUE, sep=",")
products <- read.csv(file="/Users/davidch91/Desktop/product.csv", header=TRUE, sep=",")

#merging both data frames
total <- merge(demos,products,by="hhid")

#Report of the results and the interpretations of coefficients can be found in a separate pdf document uploaded also on Github

#--------------------------------------------------------------
#EXERCISE 1 - Data Description
#--------------------------------------------------------------

#average and dispersion of prod char (1a) - average and stdev of price of each ?
price <- total[,12:21]
average_price <- colMeans(price)
dispersion_price <- apply(price,2,sd)

#There are 7 brands and 2 product types
#market share of each brand (1b) - tabulate choice
R<-tabulate(products$choice) / 4470 * 100
Pk<-R[1]+R[8]
Fl<-R[3]+R[9]
Hse<-R[4]+R[10]

#market share by product characteristics (stick 89.69% vs tub 10.31%)
Stk<-sum(R[1:7])
Tub<-100-Stk

#crosstab observed attributes (of the households) vs choice (1c)

#Observed Attributes here are income,Fs3_4,Fs5.,Fam_Size,college,whtcollar,retired, 
#We want to cross tabulate them all with choice

#choice by Income:
CT1<-table(total$Income,total$choice)
#choice by Fs3_4:
CT2<-table(total$Fs3_4,total$choice)
#choice by Fs5.:
CT3<-table(total$Fs5.,total$choice)
#choice by Fam_Size:
CT4<-table(total$Fam_Size,total$choice)
#choice by college:
CT5<-table(total$college,total$choice)
#choice by whtcollar:
CT6<-table(total$whtcollar,total$choice)
#choice by retired:
CT7<-table(total$retired,total$choice)

#--------------------------------------------------------------
#EXERCISE 2 - First Model - Conditional Logit
#--------------------------------------------------------------

#We start from a random initial guess of beta's and alpha's values and combine them into theta
betaclogit <- c(0.5)
set.seed(11)
alpha<-c(rnorm(9, mean=0.2, sd=0.05))
theta<-c(betaclogit,alpha)
price <- total[,12:21]

#We set choice 1 as base and change first column into 0
baseprice1<-price[,1]
pricetilda<-price-baseprice1

#Setting up the function that returns us negative log likelihood clogit
CLOGIT_EX2 <- function(theta){
  vij <- pricetilda*theta[1]
  vij <- vij[,2:10]+theta[2:10]
  vij <- cbind(0,vij)
  denominator <- apply(exp(vij),1,sum)
  probij <- exp(vij)/denominator
  pij <- log(probij)
  #The following is equivalent to set dij=1 and multiply it to pij
  LLi <- ifelse(total$choice==1, pij[,1], ifelse(total$choice==2, pij[,2], ifelse(total$choice==3, pij[,3], ifelse(total$choice==4, pij[,4], ifelse(total$choice==5, pij[,5], ifelse(total$choice==6, pij[,6], ifelse(total$choice==7, pij[,7], ifelse(total$choice==8, pij[,8], ifelse(total$choice==9, pij[,9],ifelse(total$choice==10, pij[,10], NA))))))))))
  return(-sum(LLi))
}

#Optimization of theta
answer_clogit<-optim(par=theta,fn=CLOGIT_EX2)$par
answer_clogit

#--------------------------------------------------------------
#EXERCISE 3 - Second Model - Multinomial Logit
#--------------------------------------------------------------

#We want 2x9 parameters (alpha2-10 and beta2-10). We we use choice 1 as base, of which alpha1=beta1=0)

#Again we set a random initial guess of betas and alphas and combine them into thetam
set.seed(14)
betamlogit <- c(rnorm(9, mean=0.01, sd=0.0045))
set.seed(12)
alpham<-c(rnorm(9, mean=-2.4, sd=0.16))
thetam<-c(alpham,betamlogit)
incomei<-total[,3]

#Setting up the function that returns us the negative log likelihood of mlogit
MLOGIT_EX3 <- function(thetam){
  betaj<-thetam[10:18]
  dim(betaj)<-c(1,9)
  alphaj<-thetam[1:9]
  dim(alphaj)<-c(1,9)
  vij<-(incomei%*%betaj)
  vijmlogit <- data.frame(matrix(nrow=4470, ncol=9))
  for(i in 1:9){
    vijint<-vij[,i]+alphaj[,i]
    vijmlogit[,i]<-vijint
  }
  #First column covering base choice 1 will be = 0, since a1=b1=0. Hence, we add 0 to the left of the other 9 columns
  vijmlogit<-cbind(0,vijmlogit)
  denominator <- apply(exp(vijmlogit),1,sum)
  probij <- exp(vijmlogit)/denominator
  pij <- log(probij)
  LLi <- ifelse(total$choice==1, pij[,1], ifelse(total$choice==2, pij[,2], ifelse(total$choice==3, pij[,3], ifelse(total$choice==4, pij[,4], ifelse(total$choice==5, pij[,5], ifelse(total$choice==6, pij[,6], ifelse(total$choice==7, pij[,7], ifelse(total$choice==8, pij[,8], ifelse(total$choice==9, pij[,9],ifelse(total$choice==10, pij[,10], NA))))))))))
  return(-sum(LLi))
}

#Optimization
answermlogit<-optim(par=thetam,fn=MLOGIT_EX3)$par
answermlogit
                 
#--------------------------------------------------------------
#EXERCISE 4 - Marginal Effects
#--------------------------------------------------------------

#Marginal Effect for Conditional Logit

#We repeat our process earlier to get the pij from optimal parameter
betaclogit <- c(0.5)
set.seed(11)
alpha<-c(rnorm(9, mean=0.2, sd=0.05))
theta<-c(betaclogit,alpha)
price <- total[,12:21]
baseprice1<-price[,1]
pricetilda<-price-baseprice1

#Setting up a function that will return us the optimal pij
PIJ_CLOGIT <- function(theta){
  vij <- pricetilda*theta[1]
  vij <- vij[,2:10]+theta[2:10]
  vij <- cbind(0,vij)
  denominator <- apply(exp(vij),1,sum)
  probij <- exp(vij)/denominator
  return(probij)
}
probij<-PIJ_CLOGIT(answer_clogit)
as.matrix(probij)

#Setting up dummy variable
deltaijk<-matrix(0,nrow=4470,ncol=10)
for (i in 1:10){
  deltaijk[,i]<-ifelse(total$choice==i,1,0)
}

#Setting up the multiplier of pij
center<-deltaijk-probij

#Multiply it again by the optimal beta
ME_clogit<-(probij*center)*answer_clogit[1]

#Average ME of all households
avgME_clogit<-colMeans(ME_clogit)
avgME_clogit

#Marginal Effect for Multinomial Logit

#First, find what the betai_bar is from our optimal solution for beta_mlogit in Ex3

beta_optimal_mlogit <- answermlogit[10:18]

#optimal beta for choice1 will be 0 as choice 1 is designated as base. Hence, we add 0 to the left.
beta_optimal_mlogit <- c(0,beta_optimal_mlogit)
as.matrix(beta_optimal_mlogit)
beta_optimal_mlogit<-t(beta_optimal_mlogit)

#Next we want to find the pij of the optimal beta
set.seed(14)
betamlogit <- c(rnorm(9, mean=0.01, sd=0.0045))
set.seed(12)
alpham<-c(rnorm(9, mean=-2.4, sd=0.16))
thetam<-c(alpham,betamlogit)
incomei<-total[,3]

#Setting up a function that returns us optimal pij
PROBIJ_MLOGIT <- function(thetam){
  betaj<-thetam[10:18]
  dim(betaj)<-c(1,9)
  alphaj<-thetam[1:9]
  dim(alphaj)<-c(1,9)
  vij<-(incomei%*%betaj)
  vijmlogit <- data.frame(matrix(nrow=4470, ncol=9))
  for(i in 1:9){
    vijint<-vij[,i]+alphaj[,i]
    vijmlogit[,i]<-vijint
  }
  vijmlogit<-cbind(0,vijmlogit)
  denominator <- apply(exp(vijmlogit),1,sum)
  probij <- exp(vijmlogit)/denominator
  return(probij)
}
probij_mlogit<-PROBIJ_MLOGIT(answermlogit)
as.matrix(probij_mlogit)
betai_bar<-probij_mlogit * beta_optimal_mlogit

#summation by j, to get betai_bar of Nx1
betai_bar<-apply(betai_bar,1,sum)

#replicating beta_optimal_mlogit for 4470 times to get "betaj" in the formula
df<-data.frame(beta_optimal_mlogit)
as.matrix(df)
df <- as.data.frame(lapply(df, rep, 4470))
beta_optimal_mlogit<-df

#betaj minus betai_bar
parenthesis<-beta_optimal_mlogit-betai_bar
ME_mlogit<-probij_mlogit*parenthesis

#Averaging the marginal effects of income for all households
avgME_mlogit<-colMeans(ME_mlogit)

#If we want to average the household's ME of income across all choices
avgME_mlogit<-mean(avgME_mlogit)

#--------------------------------------------------------------
#EXERCISE 5 - IIA using Mixed Logit
#--------------------------------------------------------------

#First we set up the mixed logit optimization using full model:
betamxl <- c(-1)
set.seed(19)
gammaj<-c(rnorm(9, mean=0.01, sd=0.005))
set.seed(18)
alphaj<-c(rnorm(9, mean=-2, sd=0.15))
thetamx<-c(betamxl,gammaj,alphaj)

MXLOGIT_EX5a <- function(thetamx){
  #supplying conditional part of vij
  vij<-(pricetilda*thetamx[1])
  vij <- vij[,2:10]
  vij <- cbind(0,vij)
  
  #we will add the choice-specific intercepts later
  
  #supplying the multinomial part to the vij and the choice-specific intercept alphaj
  gamma<-thetamx[2:10]
  dim(gamma)<-c(1,9)
  alpha<-thetamx[11:19]
  dim(alpha)<-c(1,9)
  vijml<-incomei%*%gamma
  vijmlogit <- data.frame(matrix(nrow=4470, ncol=9))
  for(i in 1:9){
    vijint<-vijml[,i]+alpha[,i]
    vijmlogit[,i]<-vijint
  }
  vijmlogit<-cbind(0,vijmlogit)
  
  #adding conditional+multinomial part together to vij
  vijmixed<-vij+vijmlogit
  
  #Finding pij
  denominator <- apply(exp(vijmixed),1,sum)
  pij <- exp(vijmixed)/denominator
  pij <- log(pij)
  LLi <- ifelse(total$choice==1, pij[,1], ifelse(total$choice==2, pij[,2], ifelse(total$choice==3, pij[,3], ifelse(total$choice==4, pij[,4], ifelse(total$choice==5, pij[,5], ifelse(total$choice==6, pij[,6], ifelse(total$choice==7, pij[,7], ifelse(total$choice==8, pij[,8], ifelse(total$choice==9, pij[,9],ifelse(total$choice==10, pij[,10], NA))))))))))
  return(-sum(LLi))
}

answermixedfull<-optim(par=thetamx, fn=MXLOGIT_EX5a)$par
answermixedfull

#To test IIA, we run regression once more, removing choice==10 for example. Then we will compare both sets of parameters.

#Removing choice 10 in our price
pricesubset <- total[,12:20]
baseprice1sub<-pricesubset[,1]
pricetildasubset<-pricesubset-baseprice1sub

betamxls <- c(-1)
set.seed(19)
gammajs<-c(rnorm(8, mean=0.01, sd=0.005))
set.seed(18)
alphajs<-c(rnorm(8, mean=-2, sd=0.15))
thetamxs<-c(betamxls,gammajs,alphajs)

MXLOGIT_EX5b <- function(thetamxs){
  #supplying conditional part of vij
  vijs<-pricetildasubset*thetamxs[1]
  vijs<- vijs[,2:9]
  vijs<- cbind(0,vijs)
  
  #supplying the multinomial part to the vij and choice-specific intercepts
  gammas<-thetamxs[2:9]
  dim(gammas)<-c(1,8)
  alphas<-thetamxs[10:17]
  dim(alphas)<-c(1,8)
  vijml<-incomei%*%gammas
  vijmxsub <- data.frame(matrix(nrow=4470, ncol=8))
  for(i in 1:8){
    vijints<-vijml[,i]+alphas[,i]
    vijmxsub[,i]<-vijints
  }
  vijmxsub<-cbind(0,vijmxsub)
  
  #adding conditional+multinomial part together to vij
  vijmixed<-vijs+vijmxsub
  denominator <- apply(exp(vijmixed),1,sum)
  pij <- exp(vijmixed)/denominator
  pij <- log(pij)
  LLi <- ifelse(total$choice==1, pij[,1], ifelse(total$choice==2, pij[,2], ifelse(total$choice==3, pij[,3], ifelse(total$choice==4, pij[,4], ifelse(total$choice==5, pij[,5], ifelse(total$choice==6, pij[,6], ifelse(total$choice==7, pij[,7], ifelse(total$choice==8, pij[,8], ifelse(total$choice==9, pij[,9], 0)))))))))
  return(-sum(LLi))
}
answermixedsubset<-optim(par=thetamxs,fn=MXLOGIT_EX5b)$par
answermixedsubset

#Now to test for IIA, we compare the log likelihood from both sets of parameters we found

#LL of full model and restricted model respectively:
BF<--MXLOGIT_EX5a(answermixedfull)
BR<--MXLOGIT_EX5b(answermixedsubset)

#McFadden Train Tye Test
MTT<--2*(BF-BR)

#We get MTT value of 693.02. Next we compare this with critical value of chi distribution
#Degree of freedom is equal to the number of parameters in the restricted model, i.e. 17
#Suppose we use significance level 95%, for which the critical value is 27.5871
#Since our MTT is well past the critical value, it is firmly in the rejection area.
#Conclusion is that IIA assumption is violated in our model.
