clear all
set more off, perm
set scrollbufsize 2000000
cd "/Users/davidch91/Desktop/613/HW5"
pwd
log using "HW5_result.txt", replace

*HW2 OLS and Discrete Choice

*Ex1:Data Creation
set seed 143
set obs 10000
gen X1 = runiform(1,3)
gen X2 = rgamma(3,2)
gen X3 = rbinomial(10,0.3)
gen eps = rnormal(2,1)
gen Y = 0.5+(1.2*X1)-(0.9*X2)+(0.1*X3)+eps
gen ydum=0
replace ydum=1 if Y>-0.14295

*Ex2: OLS
corr Y X1
reg Y X1 X2 X3

*Ex4&5: Discrete Choice & Marginal Effects:
logit ydum X1 X2 X3
prchange
probit ydum X1 X2 X3
prchange


*HW3 Multinomial Choices





*HW4 Linear Panel Data



use "/Users/davidch91/Desktop/613/HW5/datstu.dta"
foreach x of varlist score agey schoolcode1-rankplace {
replace `x'="" if(`x'=="NA")
}
destring score agey schoolcode1-rankplace, replace



