###Replication code for the article "Development from representation?"

library(foreign)
library(catspec)
library(car)
library(xtable)
library(sandwich)
library(lmtest)
library(exactRankTests)
library(rms)

rm(list = ls())

# Dataset of development indicators for 3,134 of India’s state assembly
# constituencies (ACs) in 1971 and 2001

load("devDTA.Rdata")

# Variables are described in codebook: ReadMe.rtf
names(devDTA)

#################################################################
############ SUMMARY STATS ######################################

#Sample size including ST constituencies

# Summary of reservation status from 1974-2000
# SC: Scheduled Caste
# ST: Scheduled Tribe
# GEN: General Constituency

table(devDTA$AC_type_1976)

#Sample size Census PCA variables without ST constituencies

# Same as AC_type_1976 bt excluding Scheduled Tribes, as analysis is limited to
# comparing Scheduled Caste Constituencies and General Consituencies

table(devDTA$AC_type_noST)

#Sample size for Village Directory variables without ST constituencies

# Subsetting AC_type_noST for assembly constituencies with data for percentage
# of rural population living in a village with electricity

table(devDTA$AC_type_noST[complete.cases(devDTA$P_elecVD01)])

###Creating clustered SEs

# Creating a generalized function for clustering standard errors by state,
# however Jensenius reports that "coefficients are similarly insignificant if
# the standard errors are clustered at the district of state assembly
# constituent levels" (p.33, 2015)

clusterSE<-function(model, data, cluster){
require(sandwich, quietly = TRUE)
require(lmtest, quietly = TRUE)
cluster<-as.factor(as.character(data[as.numeric(row.names(model.matrix(model))),cluster]))
M <- length(unique(cluster))
N <- length(cluster)
K <- model$rank
dfc <- (M/(M-1))*((N-1)/(N-K))
u.clust<-apply(estfun(model),2, function(x) tapply(x, cluster, sum))
cl.vcov<- dfc*sandwich(model, meat=crossprod(u.clust)/N)
print(coeftest(model, cl.vcov))
return(cl.vcov)
}

# List of outcome variables of interest, where the number after the comma refers
# to the column number of the variable, as delineated in the codebook

##Outcome variables 

#Psc, 35

#Plit_7, 60
#P_W, 52
#P_al, 58

#P_elecVD01, 64
#P_educVD01, 68
#P_medicVD01, 72
#P_commVD01, 76

#Plit_7_gap, 63, 
#P_W_gap, 53
#P_al_gap, 59

#P_elecVD01_gap, 67
#P_educVD01_gap, 71
#P_medicVD01_gap, 75
#P_commVD01_gap, 79

# Generating a vector of the index values of outcome variables of interest to
# allow for iterative regression over the same explanatory variables below

outcomeindex<-c(35, 60, 52, 58, 64, 68, 72, 76, 63, 53, 59, 67, 71, 75, 79)

##Values for Table 1
##Remove # from lines with alternative SE specification to check robustness to other SEs.

# Creating a matrix to store regression results with one row for each outcome
# variables of interest

mymatrix<-matrix(nrow=length(outcomeindex), ncol=4)

# Regression for loop begins 

for(i in 1:length(outcomeindex)){
devDTAer<-devDTA[complete.cases(devDTA[,outcomeindex[i]], devDTA$AC_type_noST),]

mymatrix[i,c(1,2)]<-round(tapply(devDTAer[, outcomeindex[i]], devDTAer$AC_type_noST, mean, na.rm=T),1)

myOLS<-lm(devDTAer[, outcomeindex[i]]~ devDTAer$AC_type_noST)

mymatrix[i,c(3)]<-round(myOLS$coef[2],1)

#SEs clustered at state level 
mySE<-clusterSE(myOLS, data=devDTAer, cluster="State_no_2001_old")
mymatrix[i,c(4)]<-ifelse(coeftest(myOLS, mySE)[2,4]<0.01, "<0.01", 	round(coeftest(myOLS, mySE)[2,4],2))

#Bootstrapped SEs clustered at state level using rms
#cluster<-devDTAer$State_no_2001_old
#myRMS<-ols(devDTAer[, outcomeindex[i]]~ devDTAer$AC_type_noST, x=TRUE, y=TRUE)
#mySEboot<- bootcov(myRMS, cluster, B=1000)$var
#mymatrix[i,c(4)]<-ifelse(coeftest(myRMS, mySEboot)[2,4]<0.01, "<0.01", 	round(coeftest(myRMS, mySEboot)[2,4],2))

# SEs from permutation test	
#mymatrix[i,c(4)]<-ifelse(perm.test(devDTAer[,outcomeindex[i]]~ devDTAer$AC_type_noST)$p.value<0.01, "<0.01", 
#	round(perm.test(devDTAer[,outcomeindex[i]]~ devDTAer$AC_type_noST)$p.value,2))

# Regression for loop ends
}

row.names(mymatrix)<-names(devDTA[outcomeindex])
colnames(mymatrix)<-c("Mean general", "Mean reserved", "Difference", "P-value")

row.names(mymatrix)<-c("Percentage of SCs", "Literacy rate", " Employment Rate", "Agricultural laborers", "Electricity in village", "School in village ","Medical facility in village","Comm. channel in village",
"Literacy gap", " Employment gap", "Agricultural laborers gap",   "Electricity in village gap", "School in village gap","Medical facility in village gap","Comm. channel in village gap")

library(xtable)
xtable(mymatrix, align=c("l", "r", "r", "r", "r"), caption="Difference in general and SC-reserved constituencies in 2001")

######ILLUSTRATION OF EDUC CHANGE
###EDUCATION CHANGE
attach(devDTA)

myplot_gen<-rbind(mean(Plit71_nonSC[AC_type_noST=="GEN"], na.rm=T), mean(Plit_nonSC_7[AC_type_noST=="GEN"], na.rm=T))
myplot_sc<-rbind(mean(Plit71_nonSC[AC_type_noST=="SC"], na.rm=T), mean(Plit_nonSC_7[AC_type_noST=="SC"], na.rm=T))
myplot2_gen<-rbind(mean(Plit71_SC[AC_type_noST=="GEN"], na.rm=T), mean(Plit_SC_7[AC_type_noST=="GEN"], na.rm=T))
myplot2_sc<-rbind(mean(Plit71_SC[AC_type_noST=="SC"], na.rm=T), mean(Plit_SC_7[AC_type_noST=="SC"], na.rm=T))

#pdf(file="Figures/Fig_educ_change.pdf", height = 4, width=8)
par(mar=c(4,4,2,2))
par(mfrow=c(1,2))

plot(myplot_gen, type="l", lty=c(2), col="#00688B", ylab="Percentage", ylim=c(10, 70), xaxt="n", main="Non-SC population", xlab="Year", las=1)
lines(myplot_sc, col="#FF1493")
axis(1, at=c(1,2), labels=c("1971", "2001"))
legend("bottomright", c(paste("General (N=", summary(devDTA$AC_type_noST)[1], ")", sep=""), paste("Reserved (N=", summary(devDTA$AC_type_noST)[2], ")", sep="")), lty=c(2,1), col=c("#00688B", "#FF1493"), cex=.8)

#Adding line and text for general line
abline(h=myplot_gen[2], lty=3, col="#00688B")
arrows(x0=1.02, y0=myplot_gen[1]+2, x1=1.02, y1=myplot_gen[2]-2, code=3, length=.05)
text(x=1.02, y=myplot_gen[2]-9,  labels=paste("Change\ngeneral:\n", round(myplot_gen[2]-myplot_gen[1],2), sep=""), pos=4, cex=0.8, col = "black")

#Adding line and text for reserved line
abline(h=myplot_sc[1], lty=3, col="#FF1493")
arrows(x0=1.98, y0=myplot_sc[1]+2, x1=1.98, y1=myplot_sc[2]-2, code=3, length=.05)
text(x=1.98, y=myplot_sc[1]+5,  labels=paste("Change\nreserved:\n", round(myplot_sc[2]-myplot_sc[1],2), sep=""), pos=2, cex=0.8, col = "black")

plot(myplot2_gen, type="l", lty=c(2), col="#00688B", ylab="Percentage", ylim=c(10, 70), xaxt="n", main="SC population", xlab="Year", las=1)
lines(myplot2_sc, col="#FF1493")
axis(1, at=c(1,2), labels=c("1971", "2001"))

#Adding line and text for general line
abline(h=myplot2_gen[2], lty=3, col="#00688B")
arrows(x0=1.02, y0=myplot2_gen[1]+2, x1=1.02, y1=myplot2_gen[2]-2, code=3, length=.05)
text(x=1.02, y=myplot2_gen[2]-9,  labels=paste("Change\ngeneral:\n", round(myplot2_gen[2]-myplot2_gen[1],2), sep=""), pos=4, cex=0.8, col = "black")

#Adding line and text for reserved line
abline(h=myplot2_sc[1], lty=3, col="#FF1493")
arrows(x0=1.98, y0=myplot2_sc[1]+2, x1=1.98, y1=myplot2_sc[2]-2, code=3, length=.05)
text(x=1.98, y=myplot2_sc[1]+5,  labels=paste("Change\nreserved:\n", round(myplot2_sc[2]-myplot2_sc[1],2)), pos=2, cex=0.8, col = "black")
dev.off()

###############################################
##MATCHING MODELS
###############################################


for (i in 1:2) {

detach()

matchdta<-devDTA[complete.cases(devDTA$SC_percent71_true, devDTA$State_no_2001_old, devDTA$AC_type_noST, devDTA$Plit71_SC, devDTA$Plit_SC),] 

matchdta$SC_percent71_true<-as.numeric(as.character(matchdta$SC_percent71_true))

detach()
attach(matchdta)
dim(matchdta)
names(matchdta)
Tr<-ifelse(AC_type_noST=="SC", 1, 0)

X<-as.data.frame(cbind(as.numeric(State_no_2001_old), as.numeric(DELIM_district_no), as.numeric(PC_no_1976), SC_percent71_true))

library(Matching)

if (i==1) {
Matched_norep<-Match(Y=Plit, Tr=Tr, X=X, estimand="ATT", exact=c(TRUE, TRUE, TRUE, FALSE), replace=FALSE)
} else if (i==2) {
Matched_norep<-Match(Y=Plit, Tr=Tr, X=X, estimand="ATT", exact=c(TRUE, TRUE, TRUE, FALSE), replace=FALSE, caliper=c(0,0,0,.5))
}

summary(Matched_norep)

#####Reporting balance
if (i==1) {
bal_SC_norep1<-MatchBalance(Tr~SC_percent71_true, match.out=Matched_norep, nboots=1000, data=matchdta)
} else if (i==2) {
bal_SC_norep2<-MatchBalance(Tr~SC_percent71_true, match.out=Matched_norep, nboots=1000, data=matchdta)
}	

bal.out_norep<-MatchBalance(Tr~Pop_tot1971+ P_ST71 +Plit71_nonSC+Plit71_SC+ P_W71_nonSC + P_W71_SC +P_al71_nonSC+P_al71_SC, match.out=Matched_norep, nboots=1000, data=matchdta)

covariates<-as.data.frame(cbind(Pop_tot1971, P_ST71, Plit71_nonSC, Plit71_SC, P_W71_nonSC , P_W71_SC, P_al71_nonSC, P_al71_SC))

names(covariates)<-c("Population size", "Percentage of STs", "Literacy rate (non-SCs)", "Literacy rate (SCs)", "Employment (non-SCs)", "Employment (SCs)", "Agricultural laborers (non-SCs)", "Agricultural laborers (SCs)")

balanceTable <- function(covariates, bal.out){

  cat("\\begin{table}[ht] \n")
  cat("\\caption{Difference in means for treated and control and Balance output from matches} \n")
  cat(" \\begin{tabular}{lrrcrr} \\hline \\hline \n")
  cat("Covariate	&\\multicolumn{2}{c}{Before matching}",
    "&&\\multicolumn{2}{c}{After matching}", "\\", "\\", "\\cline{2-3} \\cline{5-6} \n",
    sep="")
  cat("& \\emph{t p}-value &KS \\emph{p}-value &&",
    "\\emph{t p}-value &KS \\emph{p}-value", "\\", "\\", "\n",
    sep="")
  z <- sapply(1:dim(covariates)[2], function(x){
    cat(names(covariates)[x], "&",
    round(bal.out$BeforeMatching[[x]]$tt$p.value,2), "&",
    ifelse(is.null(bal.out$BeforeMatching[[x]]$ks$ks.boot.pvalue) == 0,
      round(bal.out$BeforeMatching[[x]]$ks$ks.boot.pvalue,2), "---"), "&&",
    round(bal.out$AfterMatching[[x]]$tt$p.value,2), "&",
    ifelse(is.null(bal.out$AfterMatching[[x]]$ks$ks.boot.pvalue) == 0,
      round(bal.out$AfterMatching[[x]]$ks$ks.boot.pvalue,2), "---"), "\\", "\\", "\n",
      sep="")
  })
  cat("\\end{tabular} \\end{table} \n")
  }

balanceTable(covariates, bal.out_norep)

##output
treatedDTA<-matchdta[Matched_norep$index.treated,]
controlDTA<-matchdta[Matched_norep$index.control,]
treatedDTA$index.match<-c(1:dim(treatedDTA)[1])
controlDTA$index.match<-c(1:dim(controlDTA)[1])

if(i==1){
matched1<-rbind(treatedDTA, controlDTA)	
} else {
matched2<-rbind(treatedDTA, controlDTA)		
}
}

##Table 2 combines the balance output from the two matching models

###############################################
##BALANCE FIGURE
###############################################

detach(matchdta)

##Balance on SC percentage

t.test(devDTA$SC_percent71_true~devDTA$AC_type_noST)
t.test(matched1$SC_percent71_true~matched1$AC_type_noST)
t.test(matched2$SC_percent71_true~matched2$AC_type_noST)

###Figure showing balance on percentage SC
#pdf("Figures/Fig_SC_percent71_balance_AEJ.pdf", width=7, height=3)
par(mai = c(0.8, 0.3, 0.3, 0.1))
par(mfrow=c(1,3))
plot(density(devDTA$SC_percent71_true[devDTA$AC_type_noST=="SC" & complete.cases(devDTA$AC_type_noST)]), col="#FF1493", xlim=c(0,60), ylim=c(0, 0.1), main="Before matching", xlab="Percentage of SCs in constituency", las=1)
lines(density(devDTA$SC_percent71_true[devDTA$AC_type_noST=="GEN" & complete.cases(devDTA$AC_type_noST)]), col="#00688B", lty=2)
legend("topright", c(paste("General (N=", summary(devDTA$AC_type_noST)[1], ")", sep=""), paste("Reserved (N=", summary(devDTA$AC_type_noST)[2], ")", sep="")), lty=c(2,1), col=c("#00688B", "#FF1493"))

plot(density(matched1$SC_percent71_true[matched1$AC_type_noST=="SC"]), col="#FF1493", xlim=c(0,60), ylim=c(0, 0.1), main="After matching", xlab="Percentage of SCs in constituency" , las=1)
lines(density(matched1$SC_percent71_true[matched1$AC_type_noST=="GEN"]), col="#00688B", lty=2)
legend("topright", c(paste("General (N=", summary(matched1$AC_type_noST)[1], ")", sep=""), paste("Reserved (N=", summary(matched1$AC_type_noST)[2], ")", sep="")), lty=c(2,1), col=c("#00688B", "#FF1493"))

plot(density(matched2$SC_percent71_true[matched2$AC_type_noST=="SC"]), col="#FF1493", xlim=c(0,60), ylim=c(0, 0.1), main="After matching with caliper", xlab="Percentage of SCs in constituency", las=1)
lines(density(matched2$SC_percent71_true[matched2$AC_type_noST=="GEN"]), col="#00688B", lty=2)
legend("topright", c(paste("General (N=", summary(matched2$AC_type_noST)[1], ")", sep=""), paste("Reserved (N=", summary(matched2$AC_type_noST)[2], ")", sep="")), lty=c(2,1), col=c("#00688B", "#FF1493"))
dev.off()

###############################################
##MATCHING ESTIMATES
###############################################

##Calculating matching estimates for Table 3 and Figure 4
##Remove # from lines with alternative SE specification to check robustness to other SEs.

mymatrix<-matrix(nrow=length(outcomeindex), ncol=12)

for(i in 2:length(outcomeindex)){

matched1_smaller<-matched1[complete.cases(matched1[outcomeindex[i]], matched1$AC_type_noST),]
matched2_smaller<-matched2[complete.cases(matched2[outcomeindex[i]], matched2$AC_type_noST),]

mymodel<-lm(matched1_smaller[, outcomeindex[i]]~ matched1_smaller$AC_type_noST)

#SEs clustered at state level
mySE<-clusterSE(mymodel, data= matched1_smaller, cluster="State_no_2001_old")

mymatrix[i,1]<-round(mymodel$coef[2],2)
mymatrix[i,2]<-round(mymodel$coef[2]+qnorm(.975)*coeftest(mymodel, mySE)[2,2],2)
mymatrix[i,3]<-round(mymodel$coef[2]-qnorm(.975)*coeftest(mymodel, mySE)[2,2],2)
mymatrix[i,4]<-ifelse(coeftest(mymodel, mySE)[2,4]<0.01, "<0.01", 	round(coeftest(mymodel, mySE)[2,4],2))

# Naive SEs from OLS
#mymatrix[i,c(2:3)]<-round(confint(mymodel)[2,],2)
#mymatrix[i,4]<-ifelse(coef(summary.lm(mymodel))[2,4]<0.01, "<0.01", round(coef(summary.lm(mymodel))[2,4],2))

#And now with caliper
mymodel<-lm(matched2_smaller[, outcomeindex[i]]~ matched2_smaller$AC_type_noST)

#SEs clustered at state level
mySE<-clusterSE(mymodel, data= matched2_smaller, cluster="State_no_2001_old")

mymatrix[i,5]<-round(mymodel$coef[2],2)
mymatrix[i,6]<-round(mymodel$coef[2]+qnorm(.975)*coeftest(mymodel, mySE)[2,2],2)
mymatrix[i,7]<-round(mymodel$coef[2]-qnorm(.975)*coeftest(mymodel, mySE)[2,2],2)
mymatrix[i,8]<-ifelse(coeftest(mymodel, mySE)[2,4]<0.01, "<0.01", 	round(coeftest(mymodel, mySE)[2,4],2))

# Naive SEs from OLS
#mymatrix[i,c(6:7)]<-round(confint(mymodel)[2,],2)
#mymatrix[i,8]<-ifelse(coef(summary.lm(mymodel))[2,4]<0.01, "<0.01", round(coef(summary.lm(mymodel))[2,4],2))

#And now with bias adjust

SCpop<-matched2_smaller$SC_pop71_true
mymodel<-lm(matched2_smaller[, outcomeindex[i]]~ matched2_smaller$AC_type_noST+ SCpop)

#SEs clustered at state level
mySE<-clusterSE(mymodel, data= matched2_smaller, cluster="State_no_2001_old")

mymatrix[i,9]<-round(mymodel$coef[2],2)
mymatrix[i,10]<-round(mymodel$coef[2]+qnorm(.975)*coeftest(mymodel, mySE)[2,2],2)
mymatrix[i,11]<-round(mymodel$coef[2]-qnorm(.975)*coeftest(mymodel, mySE)[2,2],2)
mymatrix[i,12]<-ifelse(coeftest(mymodel, mySE)[2,4]<0.01, "<0.01", 	round(coeftest(mymodel, mySE)[2,4],2))

# Naive SEs from OLS
#mymatrix[i,c(10:11)]<-round(confint(mymodel)[2,],2)
#mymatrix[i,12]<-ifelse(coef(summary.lm(mymodel))[2,4]<0.01, "<0.01", round(coef(summary.lm(mymodel))[2,4],2))
}	
	
row.names(mymatrix)<-names(devDTA[outcomeindex])
colnames(mymatrix)<-c("Difference", "Conf.int min", "Conf.int max", "P-value", "Difference", "Conf.int min", "Conf.int max", "P-value", "Difference", "Conf.int min", "Conf.int max", "P-value")

row.names(mymatrix)<-c("Percentage SCs", "Literacy rate ", "Employment rate ", "Agricultural laborers", "Electricity in village", "School in village ","Medical facility in village","Comm. channel in village",
"Literacy gap", "Employment gap", "Agricultural laborers gap",   "Electricity in village gap", "School in village gap","Medical facility in village gap","Comm. channel in village gap")

figurematrix<-mymatrix[-1,-c(1:4)]
articlematrix<-mymatrix[-1,c(1,4,NA, 5,8,NA, 9,12)]

library(xtable)
xtable(articlematrix)

#pdf(file="Figures/Fig_matching_est.pdf", height=7, width =7)
par(mar = c(4, 12, 1, .5))

plot(x=NULL,axes=F, xlim=c(-10, 10), ylim=c(1,14),xlab="Difference in percentage points in 2001 (SC-GEN)", ylab="", cex.main=2)
# add the 0, vertical lines

abline(v=0, lty=3)
axis(side=1,tick=TRUE, las=1, cex.axis=1)
axis(side=2,at=c(14:1), labels=row.names(figurematrix), cex.axis=1, las=1, tick=F)
#Writing in variable names

for (i in 1:nrow(figurematrix)) {
arrows(x0=as.numeric(figurematrix[i,2]), x1= as.numeric(figurematrix[i,3]), y0=(15-i+.1), y1=(15-i+.1), angle=90, length=.025, code=3)

arrows(x0=as.numeric(figurematrix[i,6]), x1= as.numeric(figurematrix[i,7]), y0=(15-i-.1), y1=(15-i-.1), angle=90, length=.025, code=3)

points(x=as.numeric(figurematrix[i,1]), y=15-i+.1, pch=19, cex=.7)

points(x=as.numeric(figurematrix[i,5]), y=15-i-.1, pch=1, cex=.7 ,  col="gray40")
text(x=95, y=i, figurematrix[i,4], cex=.8)

}
legend(x=2.2, y=4, pch=c(19, 1, NA), lty=c(NA, NA, 1), legend=c("Matching est.", "Bias-adjusted est.", "95% conf. interval"), cex=.9, bg="white", merge=T)
dev.off()


#######################################################################################
##ROBUSTNESS CHECKS ON MATCHED DATA 2
#######################################################################################

#########Checking development variation by SC percent on
#DTA$Plit_71gap
#DTA$P_al_71gap
#DTA$P_W71_gap

PropSC<-matched2$SC_percent71_true 
educ_lag<-matched2$Plit71_SC
worker_lag<-matched2$P_W71_SC 
agr_lag<-matched2$P_al71_SC
stateFE<-as.factor(matched2$State_no_2001_old)

model1lm<-lm(matched2$Plit_SC_7 ~ matched2$AC_type_noST* PropSC)
model3lm<-lm(matched2$Plit_SC_7 ~ educ_lag +matched2$AC_type_noST* PropSC + stateFE)

model4lm<-lm(matched2$P_W_SC ~ matched2$AC_type_noST* PropSC)
model6lm<-lm(matched2$P_W_SC ~ worker_lag +matched2$AC_type_noST* PropSC + stateFE)

model7lm<-lm(matched2$P_al_SC ~ matched2$AC_type_noST* PropSC)
model9lm<-lm(matched2$P_al_SC ~ agr_lag +matched2$AC_type_noST* PropSC + stateFE)

summary(model1lm)
summary(model3lm)
summary(model4lm)
summary(model6lm)
summary(model7lm)
summary(model9lm)

model1lm$se<-clusterSE(model1lm, data=matched2, cluster="State_no_2001_old")
model3lm$se<-clusterSE(model3lm, data=matched2, cluster="State_no_2001_old")
model4lm$se<-clusterSE(model4lm, data=matched2, cluster="State_no_2001_old")
model6lm$se<-clusterSE(model6lm, data=matched2, cluster="State_no_2001_old")
model7lm$se<-clusterSE(model7lm, data=matched2, cluster="State_no_2001_old")
model9lm$se<-clusterSE(model9lm, data=matched2, cluster="State_no_2001_old")


library(apsrtable)
table_OLS<-apsrtable(model1lm, model3lm, model4lm, model6lm, model7lm, model9lm, se="robust", omitcoef=c(6:19),  coef.names=c("Intercept", "SC reserved", "Percentage SC", "SC reserved * Percentage SC", "Literacy SC in 1971", "Worker SC in 1971", "Agr. laborer SC in 1971"))

#Table 4
table_OLS


####Look at within constituency patterns

#load full village dataset
load("Vill_AC.RData")
dim(vill_con) 
names(vill_con)

vill_con$VD01_state_id<-as.numeric(as.character(vill_con$VD01_state_id))
summary(vill_con$VD01_state_id)
summary(vill_con$VD01_AC_id)
summary(matched2$State_number_2001)
summary(matched2$AC_no_2001)
##reducde to ACs that are in matched2
vill<-merge(vill_con, matched2[,c(1:2, 28)], by.x=c("VD01_state_id", "VD01_AC_id"), by.y=c("State_number_2001", "AC_no_2001"))

names(vill)
dim(vill)

PropSC_vill<-(as.numeric(as.character(vill$VD01_sc_p))/as.numeric(as.character(vill$VD01_t_p)))
states<-as.factor(vill$VD01_state_id)

model1glm<-glm(vill$VD01_power_supl ~ vill$AC_type_noST*PropSC_vill +as.factor(vill$VD01_state_id), family=binomial(link = "logit"))

model2glm<-glm(vill$VD01_educ ~ vill$AC_type_noST*PropSC_vill +as.factor(vill$VD01_state_id), family=binomial(link = "logit"))

model3glm<-glm(vill$VD01_medic ~ vill$AC_type_noST*PropSC_vill +as.factor(vill$VD01_state_id), family=binomial(link = "logit"))

model4glm<-glm(vill$VD01_comm ~ vill$AC_type_noST*PropSC_vill +as.factor(vill$VD01_state_id), family=binomial(link = "logit"))

##Try to cluster errors by AC, district, and state
#SEs clustered at state level. Remove # to check clustered at other levels
#cluster<-"VD01_uniqueAC"
#vill$VD01_district_unique<-paste(vill$VD01_state_id, "-", vill$VD01_district_id, sep="")
#cluster<-"VD01_district_unique"
cluster<-"VD01_state_id"

model1glm$se<-clusterSE(model1glm, data=vill, cluster=cluster)
model2glm$se<-clusterSE(model2glm, data=vill, cluster=cluster)
model3glm$se<-clusterSE(model3glm, data=vill, cluster=cluster)
model4glm$se<-clusterSE(model4glm, data=vill, cluster=cluster)

library(apsrtable)
table_logit<-apsrtable(model1glm, model2glm, model3glm , model4glm, se="both", stars=1, omitcoef=c(4:19), coef.names=c("Intercept", "SC reserved", "Proportion SC", "SC reserved * Proportion SC"))

##Table 5
table_logit
