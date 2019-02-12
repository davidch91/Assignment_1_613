setwd("/Users/davidch91/Desktop/613/HW2")

#Exercise 1

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
OLS <- data.frame(y,constant,x1,x2,x3,eps)

#Exercise 2

#To compute correlation between Y and X1, we divide Cov(X1,Y) by (VarX1)^0.5 * (VarY)^0.5.
#If we do it without built-in command, Cov(X1,Y) can be obtained by summing (X1i-X1bar)(Yi-Ybar) and divide it by n-1
#Since correlation could only take up a value of 0-1, it is indeed quite different from 1.2.

covarx1 <- {x1-mean(x1)}*{y-mean(y)}
sum(covarx1)/9999
{{var(x1)}^0.5}*{{var(y)}^0.5}

#Here the result is 0.4147/1.9342, which roughly is equal to 0.2144
#To make sure, we get the same result with the built-in command in R using cor.test
cor.test(x1,y)

#To compute the coefficient of the regression, B=(B0,B1,B2,B3), the formula for OLS is B=(X'X)^-1X'Y
X <- cbind(constant,x1,x2,x3)
as.matrix(X)
as.matrix(y)

#Matrix multiplication (X'X)
xprimex <- t(X) %*% X

#Creating inverse of (X'X)
xpxinv <- solve(xprimex)

#Multiply (X'X)^-1 with X'y to get B=(B0,B1,B2,B3)
betaols <- xpxinv %*% t(X) %*% y
#We get the result that B0=2.462, B1=1.220, B2=-0.899, B3=0.074

b0 <- betaols[1,]
b1 <- betaols[2,]
b2 <- betaols[3,]
b3 <- betaols[4,]

#To calculate the Standard Error of OLS we obtain it from the square root of the sum of error square (ei^2) divided by n

#Creating sum of error square by substracting predicted y from y
ei2 <- {y-{(b0*1)+(b1*x1)+(b2*x2)+(b3*x3)}}^2
se<- {sum(ei2)/10000}^0.5
#We get the result that SE from OLS is about 0.9928

#Next, creating a bootstrap with 49 replications:
#We have to loop this whole process 49 times
#This will result in SE49 vector consisting of 49 SE's from 49 different coefficients/samples)

nboot <- 49
bootSE49 <- rep(NA,nboot)
set.seed(15)
for(i in 1:nboot){
    x1 <- runif(10000,1,3)
    x2 <- rgamma(10000,shape=3,scale=2)
    x3 <- rbinom(10000,size=1,prob=0.3)
    eps <- rnorm(10000,2,1)
    y <- 0.5+1.2*x1-0.9*x2+0.1*x3+eps
    constant <- replicate(n=10000,1)
    OLS <- data.frame(y,constant,x1,x2,x3,eps)
    X <- cbind(constant,x1,x2,x3)
    as.matrix(X)
    as.matrix(y)
    xprimex <- t(X) %*% X
    xpxinv <- solve(xprimex)
    betaols <- xpxinv %*% t(X) %*% y
    b0 <- betaols[1,]
    b1 <- betaols[2,]
    b2 <- betaols[3,]
    b3 <- betaols[4,]
    ei2 <- {y-{(b0*1)+(b1*x1)+(b2*x2)+(b3*x3)}}^2
    bootSE49[i] <- {sum(ei2)/10000}^0.5
}

#Next, create another loop for 499 times (resulting in SE499 vector of 499 SE's from 499 different coefficients/samples)

nboot <- 499
bootSE499 <- rep(NA,nboot)
set.seed(15)
for(i in 1:nboot){
  x1 <- runif(10000,1,3)
  x2 <- rgamma(10000,shape=3,scale=2)
  x3 <- rbinom(10000,size=1,prob=0.3)
  eps <- rnorm(10000,2,1)
  y <- 0.5+1.2*x1-0.9*x2+0.1*x3+eps
  constant <- replicate(n=10000,1)
  OLS <- data.frame(y,constant,x1,x2,x3,eps)
  X <- cbind(constant,x1,x2,x3)
  as.matrix(X)
  as.matrix(y)
  xprimex <- t(X) %*% X
  xpxinv <- solve(xprimex)
  betaols <- xpxinv %*% t(X) %*% y
  b0 <- betaols[1,]
  b1 <- betaols[2,]
  b2 <- betaols[3,]
  b3 <- betaols[4,]
  ei2 <- {y-{(b0*1)+(b1*x1)+(b2*x2)+(b3*x3)}}^2
  bootSE499[i] <- {sum(ei2)/10000}^0.5
}

#Comparing standard errors from Bootstrap of 49, 499 and OLS
mean(bootSE49)
mean(bootSE499)

#From boot49 we get 0.999699, from boot499 we get 0.9995551, and from OLS we get 0.9928
#The more number of replication the bootstrap has, the closer their standard errors get to that of OLS in this particular case.

#Exercise 3 - Numerical Optimization

#We want to set up or write the probit function to estimate ydum given X
#The formula for the likelihood fn in probit involves the cdf of standard normal distribution, i.e. chi(X)
#F(XB)=Chi(XB) where B=(B0,B1,B2,B3) from OLS

#Create an empty 4 rows vector named betadiscrete, consisting of B0, B1, B2, B3 that we hope optimization can give us.
betadiscrete <- rep(NA,4)

#Set up a vector yhat, as equal to X multiplied with betadiscrete (XiB)
yhat <- X %*% betadiscrete

#The cdf of standard normal distribution, F(XiB)=F(yhat)=Chi(yhat), is equal to cdfprobit below:
int <- function(yhat) (exp(0.5*(yhat^2)))
cdfprobit <- function(yhat) (1/sqrt(2*pi))*integrate(int, lower=-Inf, upper=yhat)

#Next, we set up the log likelihood fn (from which to maximize). Obtain by summation of (F(XiB)^ydumi)(1-F(XiB)^1-ydumi)
#We call this function likeprobit

likeprobit <- function(cdfprobit) {
  a=sum((cdfprobit*ydum)+((1-cdfprobit)*(1-ydum))
  return(a)
}

#We now apply optimization to maximize likeprobit using steepest ascent method


#Exercise 4 - Discrete Choice (I'm sorry the codes are lousy, I got lost constructing the function, codes and their specific options)

#Creating the dataframe
Probit <- data.frame(ydum,yhat,y,x1,x2,x3)

#Writing the functions
#For Probit model is essentially the same with Ex 3, shown by likeprobit
#However, if we change it into negative log likelihood function, the function is as follows:

nllprobit <- function(cdfprobit) {
  ap <- sum(ydum*log(cdfprobit)+(1-ydum)*log(1-cdfprobit))
  return(-ap)
}

#For Logit model, we call it nlllogit, or negative log likelihood logit
#We get from first setting up logistic function F(XiB) first
cdflogit <- ((exp^yhat)/(1+exp^yhat))
nlllogit <- function(cdflogit) {
  al <- sum(ydum*log(cdflogit)+(1-ydum)*log(1-cdflogit))
  return(-al)
}

#We can use optim to maximize likeprobit and likelogit.
#Since optim command by design minimizes a function,
#to maximize llprobit and lllogit, we need to minimize -likeprobit & -likelogit

optim(par=betadiscrete, fn=nllprobit)
optim(par=betadiscrete, fn=nlllogit)

#Optimal value of betadiscrete that correspond with highest llprobit or lllogit are:

#Exercise 5 - Marginal Effects

#To compute marginal effect on probit/logit, we need to take derivative of conditional probability of Y equals to 1, with respect to X.
