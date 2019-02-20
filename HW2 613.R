setwd("/Users/davidch91/Desktop/613/HW2")

#Exercise 1 - Data Creation

#Setting the seed and draw
seed <- 15
draw <- 10000
set.seed(15)

#Setting up the vectors
x1 <- runif(10000,1,3)
x2 <- rgamma(10000,shape=3,scale=2)
x3 <- rbinom(10000,size=1,prob=0.3)
eps <- rnorm(10000,2,1)

#Setting up the variables
y <- 0.5+1.2*x1-0.9*x2+0.1*x3+eps
ydum <- replicate(n=10000,0)
ydum[y>mean(y)] <- 1

#Creating a vector of ones
constant <- replicate(n=10000,1)

#Setting up a dataframe to combine all vectors
OLS = data.frame(y,constant,x1,x2,x3)

#Exercise 2

#To compute correlation between Y and X1, we divide Cov(X1,Y) by (VarX1)^0.5 * (VarY)^0.5.
#If we do it without built-in command, Cov(X1,Y) can be obtained by summing (X1i-X1bar)(Yi-Ybar) and divide it by n-1
#Since correlation could only take up a value of 0-1, it is indeed quite different from 1.2.

corrx1y <- (sum((x1-mean(x1))*(y-mean(y)))/9999)/(sqrt(var(x1))*(sqrt(var(y))))
print(corrx1y)

#Here the result is 0.4147/1.9342, which roughly is equal to 0.2144
#To make sure, we get the same result with the built-in command in R using cor.test
cor.test(x1,y)

#To compute the coefficient of the regression, B=(B0,B1,B2,B3), the formula for OLS is B=(X'X)^-1X'Y
X <- cbind(constant,x1,x2,x3)
as.matrix(X)
as.matrix(y)

#Computing Beta OLS
betaols <- solve(t(X)%*%(X)) %*% t(X) %*% y
#We get the result that B0=2.462, B1=1.220, B2=-0.899, B3=0.074

#To calculate the Standard Error of OLS we obtain it from the square root of the sum of error square (ei^2) divided by n

#Creating sum of error square by substracting predicted y from y
error2 <- (y-(X %*% betaols))^2

#Divide sum of error square by n-k to get sigma square hat
sigmasqhat <- (sum(error2))/(nrow(X)-ncol(X))

#Obtain variance-covariance matrix & standard errors
vcvmatrix <- sigmasqhat[1] * solve(t(X) %*% X)
seols <- sqrt(diag(vcvmatrix))
print(seols)

#Using BOOTSTRAP to compute std errors
set.seed(15)
bs_sample <- X[sample(nrow(X),size=10000, replace=TRUE),]
nboot <- 49
se_ols_bs49 <- matrix(ncol=4, nrow=nboot)
for(i in 1:nboot){
  bs_sample <- X[sample(nrow(X),size=10000, replace=TRUE),]
  beta_bs <- solve(t(bs_sample) %*% bs_sample) %*% t(bs_sample) %*% y
  error_bs <- y-(X %*% beta_bs)
  sigmasqhat_bs <- (t(error_bs) %*% error_bs)/(nrow(bs_sample)-ncol(bs_sample))
  vcvmatrix_bs <- sigmasqhat_bs[1,] * solve(t(bs_sample) %*% bs_sample)
  seols_bs49 <- sqrt(diag(vcvmatrix_bs))
  se_ols_bs49[i,] <- seols_bs
}
avgse_ols_bs49 <- colMeans(se_ols_bs49)
print(avgse_ols_bs49)

set.seed(15)
bs_sample <- X[sample(nrow(X),size=10000, replace=TRUE),]
nboot <- 499
se_ols_bs499 <- matrix(ncol=4, nrow=nboot)
for(i in 1:nboot){
  bs_sample <- X[sample(nrow(X),size=10000, replace=TRUE),]
  beta_bs <- solve(t(bs_sample) %*% bs_sample) %*% t(bs_sample) %*% y
  error_bs <- y-(X %*% beta_bs)
  sigmasqhat_bs <- (t(error_bs) %*% error_bs)/(nrow(bs_sample)-ncol(bs_sample))
  vcvmatrix_bs <- sigmasqhat_bs[1,] * solve(t(bs_sample) %*% bs_sample)
  seols_bs499 <- sqrt(diag(vcvmatrix_bs))
  se_ols_bs499[i,] <- seols_bs499
}
avgse_ols_bs499 <- colMeans(se_ols_bs499)
print(avgse_ols_bs499)

#Exercise 3 - Numerical Optimization

#We want to find a set of betas that maximize the objective fn, which in this case is the log likelihood of probit. Let's start from an arbitrarily define beta
beta <- c(0.16,-0.0060,-0.0040,0.01)

#This is the function we want to maximize. This function returns total log likelihood given a set of betas
probit_likelihood <- function(beta) {
  yhat <- X %*% beta
  sllp<-sum((ydum*log(pnorm(yhat[,1]))+(1-ydum)*log(1-(pnorm(yhat[,1])))))
  return(sllp)
}

#In this exercise, we're using steepest ascent to find the optimum beta
#We start from an initial guess (ig) of a vector of betas (B0,B1,B2,B3)
betaig <- c(2.1,1.0,-0.7,0.01)

#We set a small number h as a finite difference to compute derivative
h = 0.00005

#Next, we perform steepest ascent procedures, the coding of which are
betatmp <- betaig
print(cat("a: ",betatmp))
llprev <- 0
lllatest <- probit_likelihood(betaig)
i = 0
while(lllatest - llprev > 0.0001 || lllatest - llprev < -1000) {
  #Computing gradient consisting of 4 partial derivatives. Could've used loops but couldn't figure out the code.
  betatmppos <- betaig
  betatmpneg <- betaig
  betatmppos[1] <- betaig[1]+h
  betatmpneg[1] <- betaig[1]-h
  gk0 <- (probit_likelihood(betatmppos)-probit_likelihood(betatmpneg))/(2*h)
  betatmppos <- betaig
  betatmpneg <- betaig
  betatmppos[2] <- betaig[2]+h
  betatmpneg[2] <- betaig[2]-h
  gk1 <- (probit_likelihood(betatmppos)-probit_likelihood(betatmpneg))/(2*h)
  betatmppos <- betaig
  betatmpneg <- betaig
  betatmppos[3] <- betaig[3]+h
  betatmpneg[3] <- betaig[3]-h
  gk2 <- (probit_likelihood(betatmppos)-probit_likelihood(betatmpneg))/(2*h)
  betatmppos <- betaig
  betatmpneg <- betaig
  betatmppos[4] <- betaig[4]+h
  betatmpneg[4] <- betaig[4]-h
  gk3 <- (probit_likelihood(betatmppos)-probit_likelihood(betatmpneg))/(2*h)
  g <- c(gk0,gk1,gk2,gk3)
  alpha <- 0.000005
  alphastar <- alpha
  
  #The following is a function of negative sum of LL of beta from next iteration, which we want to minimize with optim
  NLLnextit<-function(alphastar){
    betanext <- betaig+(alphastar*g)
    yhatnext <- X %*% beta
    ss <- sum((ydum*log(pnorm(yhatnext[,1]))+(1-ydum)*log(1-(pnorm(yhatnext[,1])))))
    return(-ss)
  }
  optimal<-optim(par=alphastar,fn=NLLnextit)
  optimalphastar <- optimal$par
  print(optimalphastar)
  
  #Updating position of the beta for next iteration:
  betaig<-betaig+(optimalphastar*g)
  llprev <- lllatest
  lllatest <- probit_likelihood(betaig)
  print(lllatest-llprev)
  print(betaig)
  i <- i + 1
}
#Showing the final beta after 2-3 minutes of iterating
print(betaig)

#The result we get from steepest ascent optimization is indeed quite close to the true parameters, except for the intercept, and slightly less so for x3.


#EXERCISE 4: DISCRETE CHOICE:

#PROBIT MODEL

#For probit, we essentially use the same function written in Exercise 3, i.e. probit_likelihood()
#But we change it to negative since command optim minimizes function by default.
#It will return us a positive number, i.e. minus of total log likelihood of probit

beta <- betaig
probit_neglikelihood <- function(beta) {
  yhat <- X %*% beta
  llprob <- sum((ydum*log(pnorm(yhat[,1]))+(1-ydum)*log(1-(pnorm(yhat[,1])))))
  nllprob <- -llprob
  return(nllprob)
}

#Command optim finds beta that minimizes this function
beta_probit_optim <- (optim(par=beta, fn=probit_neglikelihood)$par)
print(beta_probit_optim)

#Computing standard error of Betas Probit obtained from optim using Bootstrap 49 & Significance

seed<-15
set.seed(15)
nboot <- 49
beta_probit_bs <- matrix(ncol=4, nrow=nboot)
for(i in 1:nboot){
  bs_sample <- X[sample(nrow(X),size=10000, replace=TRUE),]
  beta_bs <- betaig
  probit_neglikelihood <- function(beta_bs) {
    yhat <- bs_sample %*% beta_bs
    llprob <- sum((ydum*log(pnorm(yhat[,1]))+(1-ydum)*log(1-(pnorm(yhat[,1])))))
    nllprob <- -llprob
    return(nllprob)
  }
  optim_beta_probit_bs <- optim(par=beta, fn=probit_neglikelihood)
  optimalbeta_bs <- optim_beta_probit_bs$par
  beta_probit_bs[i,] <- optimalbeta_bs
}
print(outputprobit)

#We compute standard error using two different methods

#This is simply taking the sd of each column of 49x4 matrix of betas obtained from Bootstrap
sd_beta_probit_bs <- apply(beta_probit_bs,2,sd)

#This is by obtaining variance-covariance matrix
vcm_beta_probit_bs <- cov(beta_probit_bs)
se_beta_probit_bs <- sqrt(diag(vcm_beta_probit_bs))

#The results from both are the same
comparepro <- rbind(sd_beta_probit_bs,se_beta_probit_bs)
print(comparepro)

#Computing t-statistic and significance of optimum beta in Probit
tstat_probit <- beta_probit_optim/se_beta_probit_bs
print(tstat_probit)
#All of X0, X1, X2 are significant statistically, while X3 isn't in this particular case.

#LOGIT MODEL
#First we need to write the log likelihood function to maximize, i.e. total log likelihood logit
#Form is similar with probit, but with different cdf used. Here we use logistic instead of normal dist.
#Again we create negative log likelihood of logit that will return us a positive number, because command optim minimizes by default

beta <- betaig
logit_neglikelihood <- function(beta){
  yhat <- X %*% beta
  ll <- sum(ydum*log(exp(yhat[,1])/(1+exp(yhat[,1])))+(1-ydum)*log(1-(exp(yhat[,1])/(1+exp(yhat[,1])))))
  return(-ll)
}

#Next find beta that minimizes total negative log likelihood logit:
beta_logit_optim <- (optim(par=beta, fn=logit_neglikelihood)$par)
print(beta_logit_optim)

#Computing standard error of Betas Logit obtained from optim using Bootstrap 49 & Significance (tstat)

seed<-15
set.seed(15)
nboot <- 49
beta_logit_bs <- matrix(ncol=4, nrow=nboot)
for(i in 1:nboot){
  bs_sample <- X[sample(nrow(X),size=10000, replace=TRUE),]
  beta_bs <- betaig
  logit_neglikelihood <- function(beta_bs) {
    yhat <- bs_sample %*% beta_bs
    lllogit <- sum(ydum*log(exp(yhat[,1])/(1+exp(yhat[,1])))+(1-ydum)*log(1-(exp(yhat[,1])/(1+exp(yhat[,1])))))
    return(-lllogit)
  }
  optim_beta_logit_bs <- optim(par=beta_bs, fn=logit_neglikelihood)
  optimalbeta_logit_bs <- optim_beta_logit_bs$par
  beta_logit_bs[i,] <- optimalbeta_logit_bs
}
#Again, we compute standard error using two different methods

#This is simply taking the sd of each column of 49x4 matrix of betas obtained from Bootstrap
sd_beta_logit_bs <- apply(beta_logit_bs,2,sd)

#This is by obtaining variance-covariance matrix
vcm_beta_logit_bs <- cov(beta_logit_bs)
se_beta_logit_bs <- sqrt(diag(vcm_beta_logit_bs))

#The results from both are the same
comparelog <- rbind(sd_beta_logit_bs,se_beta_logit_bs)
print(comparelog)

#Computing t-statistic and significance of optimum beta in Probit
tstat_logit <- beta_logit_optim/se_beta_logit_bs
print(tstat_logit)
#Again, just like probit, all explanatory vars and intercept are significant statistically, except for X3.

#LINEAR PROBABILITY MODEL

#Having set up LPM in max likelihood setting, the objective fn to maximize is total log likelihood lpm.
#As before, we change it into negative due to the way 'optim' command works.

#Setting up the objective function:
beta <- betaig
lpm_neglikelihood <- function(beta){
  yhat <- X %*% beta
  zz <- sum((ydum*yhat)+(1-ydum)*(1-yhat))
  return(-zz)
}

#Find beta that minimizes negative log likelihood lpm using optim
beta_lpm_optim <- (optim(par=beta, fn=lpm_neglikelihood)$par)
print(beta_lpm_optim)

#Computing standard error of Betas LPM obtained from optim using Bootstrap 49 & Significance (tstat)

seed<-15
set.seed(15)
nboot <- 49
beta_lpm_bs <- matrix(ncol=4, nrow=nboot)
for(i in 1:nboot){
  bs_sample <- X[sample(nrow(X),size=10000, replace=TRUE),]
  beta_bs <- betaig
  lpm_neglikelihood <- function(beta_bs) {
    yhat <- bs_sample %*% beta_bs
    lllpm<-sum((ydum*yhat)+(1-ydum)*(1-yhat))
    return(-lllpm)
  }
  optim_beta_lpm_bs <- optim(par=beta_bs, fn=lpm_neglikelihood)
  optimalbeta_lpm_bs <- optim_beta_lpm_bs$par
  beta_lpm_bs[i,] <- optimalbeta_lpm_bs
}

#Again, we compute standard error using two different methods

#This is simply taking the sd of each column of 49x4 matrix of betas obtained from Bootstrap
sd_beta_lpm_bs <- apply(beta_lpm_bs,2,sd)

#This is by obtaining variance-covariance matrix
vcm_beta_lpm_bs <- cov(beta_lpm_bs)
se_beta_lpm_bs <- sqrt(diag(vcm_beta_lpm_bs))

#The results from both are the same
comparelpm <- rbind(sd_beta_lpm_bs,se_beta_lpm_bs)
print(comparelpm)

#Computing t-statistic and significance of optimum beta in Probit
tstat_lpm <- beta_lpm_optim/se_beta_lpm_bs
print(tstat_lpm)
#If we use 1.96 critical values, only X2 is significant statistically under LPM


#Comparing estimated coefficients from Probit, Logit, LPM

trueparameters<-c(.5,1.2,-.9,.1)
compareall<-rbind(trueparameters,beta_probit_optim,beta_logit_optim,beta_lpm_optim)
colnames(compareall)<-c("Inter","X1","X2","X3")
print(compareall)

#Additionally we can also compare estimated coefficients with the ones obtained from built-in command like GLM:

logit <- glm(ydum ~ X, family=binomial (link="logit"))
summary(logit)
probit <- glm(ydum ~ X, family=binomial (link="probit"))
summary(probit)


#Exercise 5 - Marginal Effects

#PROBIT MODEL

#Here we use optimum betas we obtained from Ex4 to compute marginal effects
#For probit, optimum betas obtained previously are (2.91999344,1.21102642,-0.89133655,0.03401492)

#To compute Marginal Effects of Probit evaluated at each i
betaoptimum <- c(2.91999344,1.21102642,-0.89133655,0.03401492)
yhat <- X %*% betaoptimum
pdfprobiti <- dnorm(yhat)
meprobit <- pdfprobiti %*% betaoptimum

#Taking average of them all, we get average marginal effects of probit
ameprobit <- colMeans(meprobit)
print(ameprobit)

#Next, we want to compute standard deviations of marginal effects of probit using Delta Method
#Using finite differences to compute the Jacobian matrix of ME Probit:
#Jacobian matrix (4x4) contains first partial derivative of each ME with respect to each of the optimum betas

#Here we construct up a loop to set up the Jacobian matrix for probit

jacprobit <- data.frame(matrix(nrow=4, ncol=4))

for(t in 1:4){
  marginal_effect_probitt <- function(betaoptimum) {
  yhat <- X %*% betaoptimum
  pdfprobiti <- dnorm(yhat)
  meprobit <- pdfprobiti %*% betaoptimum
  avgmeprobit <- colMeans(meprobit)
  return(avgmeprobit[t])
  }
  print(marginal_effect_probitt(betaoptimum))
  for(k in 1:4){
    betaoptimumpos <- betaoptimum
    betaoptimumneg <- betaoptimum
    betaoptimumpos[k] <- betaoptimum[k]+h
    betaoptimumneg[k] <- betaoptimum[k]-h
    jacobiantk <- (marginal_effect_probitt(betaoptimumpos)-marginal_effect_probitt(betaoptimumneg))/(2*h)
    jacprobit[t,k] <- jacobiantk
    }
}

#Next we construct variance-covariance matrix of marginal effects of probit by Delta Method
vcm_me_probit <- t(jacprobit) %*% vcm_beta_probit_bs %*% as.matrix(jacprobit)
se_me_probit_delta <- sqrt(diag(vcm_me_probit))
print(se_me_probit_delta)

#Computing std deviations of Marginal Effects of Probit by Bootstrap 49:

seed<-15
set.seed(15)
nboot <- 49
me_probit_bs <- matrix(ncol=4, nrow=nboot)
for(i in 1:nboot){
  bs_sample <- X[sample(nrow(X),size=10000, replace=TRUE),]
  beta_bs <- betaig
  probit_neglikelihood <- function(beta_bs) {
    yhat <- bs_sample %*% beta_bs
    llprob <- sum((ydum*log(pnorm(yhat[,1]))+(1-ydum)*log(1-(pnorm(yhat[,1])))))
    nllprob <- -llprob
    return(nllprob)
  }
  optim_beta_probit_bs <- optim(par=beta, fn=probit_neglikelihood)
  optimalbeta_bs <- optim_beta_probit_bs$par
  outputprobit[i,] <- optimalbeta_bs
  yhati <- X %*% outputprobit[i,]
  pdfprobiti <- dnorm(yhat)
  meprobiti <- pdfprobiti %*% outputprobit[i,]
  ameprobiti <- colMeans(meprobiti)
  me_probit_bs[i,] <- ameprobiti
}
print(me_probit_bs)

#Next we compute the standard deviations of ME Probit using variance-covariance matrix
vcm_me_probit_bs <- cov(me_probit_bs)
se_me_probit_bs <- sqrt(diag(vcm_me_probit_bs))
print(se_me_probit_bs)

#LOGIT MODEL

#For logit, we plug in the optimum betas obtained previously in ex4

#To compute Marginal Effects of Probit evaluated at each i
betaoptimum <- c(5.1421674,2.1666104,-1.5818825,0.0572931)
yhat <- X %*% betaoptimum
pdflogiti <- ((betaoptimum*exp(yhat))/(exp(yhat)+1)^2)
melogit <- pdflogiti %*% betaoptimum

#Taking average of them all, we get average marginal effects of probit
amelogit <- colMeans(melogit)
print(amelogit)

#Computing std deviations of ME Logit by Delta Method

#Using finite differences to compute the Jacobian matrix of ME Logit:
#Again here we use a loop to set up the Jacobian matrix for logit

jaclogit <- data.frame(matrix(nrow=4, ncol=4))

for(t in 1:4){
  marginal_effect_logitt <- function(betaoptimum) {
    yhat <- X %*% betaoptimum
    pdflogit <- ((betaoptimum*exp(yhat))/(exp(yhat)+1)^2)
    melogit <- pdflogit %*% betaoptimum
    avgmelogit <- colMeans(melogit)
    return(avgmelogit[t])
  }
  print(marginal_effect_logitt(betaoptimum))
  for(k in 1:4){
    betaoptimumpos <- betaoptimum
    betaoptimumneg <- betaoptimum
    betaoptimumpos[k] <- betaoptimum[k]+h
    betaoptimumneg[k] <- betaoptimum[k]-h
    jacobiantk <- (marginal_effect_logitt(betaoptimumpos)-marginal_effect_logitt(betaoptimumneg))/(2*h)
    jaclogit[t,k] <- jacobiantk
  }
}

#Next we construct variance-covariance matrix of marginal effects of probit by Delta Method
vcm_me_logit <- t(jaclogit) %*% vcm_beta_logit_bs %*% as.matrix(jaclogit)
se_me_logit_delta <- sqrt(diag(vcm_me_logit))
print(se_me_logit_delta)

#Computing std deviations of Marginal Effects of Logit by Bootstrap 49

seed<-15
set.seed(15)
nboot <- 49
me_logit_bs <- matrix(ncol=4, nrow=nboot)
for(i in 1:nboot){
  bs_sample <- X[sample(nrow(X),size=10000, replace=TRUE),]
  beta_bs <- betaig
  logit_neglikelihood <- function(beta_bs) {
    yhat <- bs_sample %*% beta_bs
    lllog <- sum(ydum*log(exp(yhat[,1])/(1+exp(yhat[,1])))+(1-ydum)*log(1-(exp(yhat[,1])/(1+exp(yhat[,1])))))
    return(-lllog)
  }
  optim_beta_logit_bs <- optim(par=beta_bs, fn=logit_neglikelihood)
  optimalbeta_bs <- optim_beta_logit_bs$par
  optimum_beta_logit_bs[i,] <- optimalbeta_bs
  yhati <- X %*% optimum_beta_logit_bs[i,]
  pdflogiti <- ((optimum_beta_logit_bs[i,]*exp(yhati))/(exp(yhati)+1)^2)
  melogiti <- pdflogiti %*% optimum_beta_logit_bs[i,]
  amelogiti <- colMeans(melogiti)
  me_logit_bs[i,] <- amelogiti
}
print(me_logit_bs)

#Next we compute the standard deviations of ME Logit using variance-covariance matrix
vcm_me_logit_bs <- cov(me_logit_bs)
se_me_logit_bs <- sqrt(diag(vcm_me_logit_bs))
print(se_me_logit_bs)