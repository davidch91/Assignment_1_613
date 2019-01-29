setwd("/Users/davidch91/Desktop/613")
datstu <- read.csv("datstu.csv", header=T)
attach(datstu)
datjss <- read.csv("datjss.csv", header=T)
attach(datjss)
datsss <- read.csv("datsss.csv", header=T)
attach(datsss)
install.packages("dplyr")
library(dplyr)

#Exercise 1

#number of students
length(datstu$X)

#number of schools
all_sc <- c(as.matrix(datstu[,5:10]))
length(unique(all_sc))

#number of programs
all_pgm <- c(as.matrix(datstu[,11:16]))
length(unique(all_pgm))

#number of choices
all_scpgm <- cbind(all_sc, all_pgm)
length(unique(all_scpgm))

#missing test score
sum(is.na(datstu$score))

#apply to the same school
#figure out later

#apply to less than 6 schools
sum(is.na(datstu$schoolcode6))

#Exercise 2
names(datsss)[3] <- "all_sc"
datsss$X <- NULL
datsssuni <- unique(datsss)
datstusss <- merge(uniqscpgm,datsssuni, by="all_sc")

admitstu <- subset(datstu, rankplace<7)
admit1 <- subset(admitstu, rankplace==1, select=c(X,score,agey,male,schoolcode1,choicepgm1,jssdistrict,rankplace))
names(admit1)[5] <- "schoolcode"
names(admit1)[6] <- "choicepgm"
admit2 <- subset(admitstu, rankplace==2, select=c(X,score,agey,male,schoolcode2,choicepgm2,jssdistrict,rankplace))
names(admit2)[5] <- "schoolcode"
names(admit2)[6] <- "choicepgm"
admit3 <- subset(admitstu, rankplace==3, select=c(X,score,agey,male,schoolcode3,choicepgm3,jssdistrict,rankplace))
names(admit3)[5] <- "schoolcode"
names(admit3)[6] <- "choicepgm"
admit4 <- subset(admitstu, rankplace==4, select=c(X,score,agey,male,schoolcode4,choicepgm4,jssdistrict,rankplace))
names(admit4)[5] <- "schoolcode"
names(admit4)[6] <- "choicepgm"
admit5 <- subset(admitstu, rankplace==5, select=c(X,score,agey,male,schoolcode5,choicepgm5,jssdistrict,rankplace))
names(admit5)[5] <- "schoolcode"
names(admit5)[6] <- "choicepgm"
admit6 <- subset(admitstu, rankplace==6, select=c(X,score,agey,male,schoolcode6,choicepgm6,jssdistrict,rankplace))
names(admit6)[5] <- "schoolcode"
names(admit6)[6] <- "choicepgm"
admit_all <- rbind(admit1,admit2,admit3,admit4,admit5,admit6)
admit_all %>% group_by("schoolcode","choicepgm")

#Here, admit_all is database for all admitted students and the codes of their respective admitted school/program. This is key reference point for the next sections, because cutoff and quality only care about stat of admitted students.

#Cutoff, Average, Size
aggregate(admit_all[,2:2], list(admit_all$schoolcode, admit_all$choicepgm), min)
aggregate(admit_all[,2:2], list(admit_all$schoolcode, admit_all$choicepgm), mean)
aggregate(admit_all, schoolcode~choicepgm, function(x) length(unique)) #return later - size

#Exercise 3

#Since there are six senior school choices, each student has 6 different distances from junior HS to senior HS, unless for those students who are also applying for different programs in the same school)

#First, we merge datstu with datjss via jssdistrict
datstuj <- merge(datstu,datjss, by="jssdistrict")

#Then, try to consolidate datstuj with the corresponding coordinates from datsss based on each of the 6 schoolcodes applied by each student
datstuj$senlat1=latitude1
senlat1=latitude1
datstuj$senlat2=latitude2
senlat2=latitude2
datstuj$senlat3=latitude3
senlat3=latitude3
datstuj$senlat4=latitude4
senlat4=latitude4
datstuj$senlat5=latitude5
senlat5=latitude5
datstuj$senlat6=latitude6
senlat6=latitude6

datstuj$senlong1=longitude1
senlong1=longitude1
datstuj$senlong2=longitude2
senlong2=longitude2
datstuj$senlong3=longitude3
senlong3=longitude3
datstuj$senlong4=longitude4
senlong4=longitude4
datstuj$senlong5=longitude5
senlong5=longitude5
datstuj$senlong6=longitude6
senlong6=longitude6

#Exercise 4 Cutoff Quality Distance (Avg Sd)
aggregate(admit_all[,2:2], list(admit_all$rankplace), min)
aggregate(admit_all[,2:2], list(admit_all$rankplace), mean)
aggregate(admit_all[,2:2], list(admit_all$rankplace), sd)
byranked <- group_by(admit_all, rankplace)
summarize(byranked, cutoff_avg = min(score, na.rm=TRUE))
summarize(byranked, cutoff_sd = sd(score, na.rm=TRUE))
summarize(byranked, cutoff_avg = mean(score, na.rm=TRUE))
aggregate(admit_all[,2:2], list(admit_all$rankplace), mean)

#4b: by student test score quintiles
scorequintiles <- cut(admit_all$score,5)
admit_allquint <- cbind(admit_all, scorequintiles, na.rm=TRUE)
byquintiles <- group_by(admit_allquint, scorequintiles)
summarize(byquintiles, cutoff_avg = min(score, na.rm=TRUE))
summarize(byquintiles, quality_avg = mean(score, na.rm=TRUE))
summarize(byquintiles, quality_avg = sd(score, na.rm=TRUE))

#Exercise 5
#grouping schools by decile of selectivity (cutoffs)
by_school <- group_by(admit_all, schoolcode)
cutoffsc <- cut(by_school$score,10)
scselectivity <- cbind(by_school, cutoffsc)

#grouping by student score quintiles
byquintiles <- group_by(admit_allquint, scorequintiles)
#need to merge by_school with admit_allquint, e.g. into by_scquint, before doing the same thing.