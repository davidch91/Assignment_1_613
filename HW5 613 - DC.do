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
eststo OLS: reg Y X1 X2 X3

net install spost9_ado, from(http://www.indiana.edu/~jslsoc/stata) replace

*Ex4&5: Discrete Choice & Marginal Effects:
eststo LPM: reg ydum X1 X2 X3
eststo logit: logit ydum X1 X2 X3
eststo MElogit: prchange
eststo probit: probit ydum X1 X2 X3
eststo MEprobit: prchange
esttab OLS LPM logit probit, nodepvars nogaps noeqlines se(3) title(Comparison of Coefficients Across Models)
esttab MElogit MEprobit, se(3) nodepvars nogaps noeqlines title(Marginal Effects of Logit & Probit)

*HW3 Multinomial Choices

*Having exported the margarine data from R and merged them
*We first load data(margarine) in R, then use write.dta function in R

*Merging both choicePrice and demos in Stata
use "/Users/davidch91/Desktop/613/HW5/choicePrice.dta", clear
merge m:m hhid using "/Users/davidch91/Desktop/613/HW5/demos.dta"
save "/Users/davidch91/Desktop/613/HW5/margarine.dta", replace

generate chosen= "Pk_Stk"
replace chosen="BB_Stk" if choice==2
replace chosen="Fl_Stk" if choice==3
replace chosen="Hse_Stk" if choice==4
replace chosen="Gen_Stk" if choice==5
replace chosen="Imp_Stk" if choice==6
replace chosen="SS_Tub" if choice==7
replace chosen="Pk_Tub" if choice==8
replace chosen="Fl_Tub" if choice==9
replace chosen="Hse_Tub" if choice==10
encode(chosen), gen(chosenprod)

*Ex1-Data Descriptions

*Average and dispersion of each brand's price
summarize PPk_Stk-PHse_Tub

*Creating brand & type variable
generate brand="Pk"
replace brand="BB" if choice==2
replace brand="Fl" if choice==3 | choice==9
replace brand="Hse" if choice==4 | choice==10
replace brand="Gen" if choice==5
replace brand="Imp" if choice==6
replace brand="SS" if choice==7
generate type="Stick"
replace type="Tub" if choice>6

*Market share by brand and by product type
tab brand
tab type

*Mapping between observed attributes and choices
qui tab Income choice
qui tab Fs3_4 choice
qui tab Fs5_ choice
qui tab Fam_Size choice
qui tab college choice
qui tab whtcollar choice
qui tab retired choice

*Ex2-Conditional Logit & Marginal Effects

*Expanding our dataset to create 10 (because 10 options) binary choice outcome for each individuals
gen d=_n
expand 10
bysort d: gen group_index = _n
gen choice_decision=0
gen price=0
replace price= PPk_Stk if group_index==1
replace price= PBB_Stk if group_index==2
replace price= PFl_Stk if group_index==3
replace price= PHse_Stk if group_index==4
replace price= PGen_Stk if group_index==5
replace price= PImp_Stk if group_index==6
replace price= PSS_Tub if group_index==7
replace price= PPk_Tub if group_index==8
replace price= PFl_Tub if group_index==9
replace price= PHse_Tub if group_index==10
replace choice_decision=1 if group_index==choice

gen altname="Pk_Stk"
replace altname="BB_Stk" if group_index==2
replace altname="Fl_Stk" if group_index==3
replace altname="Hse_Stk" if group_index==4
replace altname="Gen_Stk" if group_index==5
replace altname="Imp_Stk" if group_index==6
replace altname="SS_Tub" if group_index==7
replace altname="Pk_Tub" if group_index==8
replace altname="Fl_Tub" if group_index==9
replace altname="Hse_Tub" if group_index==10

global y choice_decision
global z price
global id d
global alternative altname
global basealternative1 Pk_Stk
eststo Cond_Logit: asclogit $y $z, case($id) alternatives($alternative) basealternative($basealternative1)
esttab Cond_Logit, se(3) noeqlines nodepvars compress

*Marginal Effects of price for each choice
eststo MEclogit: qui estat mfx, varlist($z)

*Ex3-Multinomial Logit & Marginal Effects
eststo mlogit: mlogit chosenprod Income
esttab Cond_Logit mlogit, se(3) nogaps noeqlines nodepvars compress wide title(Effects of Price, Income on Demand: Conditional & Multinomial Logit Regression)

*Marginal Effects of Income for each choice
mfx, predict(pr outcome(1))
estimates store MEm1
mfx, predict(pr outcome(2))
estimates store MEm2
mfx, predict(pr outcome(3))
estimates store MEm3
mfx, predict(pr outcome(4))
estimates store MEm4
mfx, predict(pr outcome(5))
estimates store MEm5
mfx, predict(pr outcome(6))
estimates store MEm6
mfx, predict(pr outcome(7))
estimates store MEm7
mfx, predict(pr outcome(8))
estimates store MEm8
mfx, predict(pr outcome(9))
estimates store MEm9
mfx, predict(pr outcome(10))
estimates store MEm10
esttab MEm1 MEm2 MEm3 MEm4 MEm5 MEm6 MEm7 MEm8 MEm9 MEm10, se(3) nocons nogaps noeqlines not nopa title(Marginal Effects of Income at Each Choice)

*Ex5-Mixed Logit & IIA
asmixlogit $y, random($z) casevars($x) alternatives($alternative) case($id)
eststo mixlogitfull

*To test for IIA, we remove one choice from our data and later compare likelihood score of both regressions:
drop if group_index==10
asmixlogit $y, random($z) casevars($x) alternatives($alternative) case($id)
eststo mixlogitpartial
hausman mixlogitfull mixlogitpartial, alleqs constant

*The result of Hausman test is -4.68, meaning that H0 is rejected and IIA is violated

*HW4 Linear Panel Data
use "/Users/davidch91/Desktop/613/HW5/KoopTobias.dta", clear

*Creating time variable and establishing panel dimension
gen t=1
replace t=t[_n-1]+1 if personid[_n]==personid[_n-1]
xtset personid timetrnd

*Random Effects
eststo Random: xtreg logwage educ potexper, re
eststo Fixed_B: xtreg logwage educ potexper, be
eststo Fixed_W: xtreg logwage educ potexper, fe
eststo Pool_LS: reg
esttab Random Fixed_B Fixed_W Pool_LS, se(3) nocons nogaps noeqlines not nopa title(Estimated Coefficients of Different Panel Data Models)

log close


