


###############################################
###   Definition / inclusion of functions   ###
###############################################

library(corrgram)   # Correlogram
library(boot)       # Drawing bootstrap samples 
library(GAMBoost)   # For boosting 
library(rms)        # For 'fast backward selection' in investigation of subsamples and boosting

source("screening.step.glm.R")   # Variable selection based on BIFs 
source("fma.step.glm.R")   # Model averaging 
source("fma.pred.glm.R")   # Predictions for model averaging 



## Backward elimination using P-values to delete predictors one-at-a-time
## 1. START with fitting full model
## 2. IDENTIFY the predictor (if any) with the largest P-value
## 3. If p > alpha then ELIMINATE the corresponding predictor and refit model
## 4. Repeat Steps 2 and 3 if predictor was identified, or 
##    STOP stepwise regression if all remaining P-values are below alpha or the null model is reached
select.by.alpha <- function(response, candidates, alpha, dat, outp=FALSE){

  tmp.maxp <- 1

  # Start elimination of variables
  while (tmp.maxp>alpha & !is.null(candidates)) {
    tmp.formula <- as.formula(paste(response, " ~ ", paste(candidates,collapse=" + ")))
    tmp.model <- lm(tmp.formula, data=dat)
    tmp.p <- summary(tmp.model)$coefficients[,"Pr(>|t|)"]
    tmp.maxp <- max(tmp.p)
    if (tmp.maxp>alpha) {
      tmp.drop <- names(tmp.p)[round(tmp.p,7)==round(max(tmp.p),7)]
      if (outp==TRUE){
        print("Current p values:")
        print(round(sort(tmp.p, decreasing=TRUE), 4))
        print(paste("-> Drop ", tmp.drop, " (p=", round(tmp.maxp,4), ")", sep=""))
      }
      candidates <- candidates[-which(candidates==tmp.drop)]
    }else{     
      break
    }  
  }
  return(final=tmp.model)
}


##  K-fold cross validation
cross.split <- function(n, K){
  seed <- .Random.seed
  out <- NULL
  if((K > n) || (K <= 1))
    stop("K outside allowable range")
  K.o <- K
  K <- round(K)
  kvals <- unique(round(n/(1:floor(n/2))))
  temp <- abs(kvals - K)
  if(!any(temp == 0))
    K <- kvals[temp == min(temp)][1]
  if(K != K.o)
    warning(paste("K has been set to", K))
  f <- ceiling(n/K)
  s <- sample(rep(1:K, f), n)
  n.s <- table(s)
  return(s)
} # end of cross.split





########################################################################################





############################
###   General settings   ### 
############################


### Margin and axis settings

## Scatter plots
xlim.sc <- ylim.sc <- c(4, 48.3)

## Bland-Altman plots
mar.ba <- c(5, 4+1, 4-2, 2) + 0.1   # (bottom, left, top, right)
xlim.ba <- c(3.6, 47.6)
ylim.ba <- c(-5.6, 3.8)

## Boxplots of SE ratios
mar.ser <- c(5-3, 4+1, 4-2, 2) + 0.1   # (bottom, left, top, right)
ylim.ser <- c(0.9,4.9)





########################################################################################





#########################################
###   Creation of Analysis Data Set   ###
#########################################

BF <- read.csv(file="bodyfatfyp.csv")
dim(BF)

varnames <- c("Sex", "Age", "Weight", "Height", "Neck", "Chest",
              "Abdomen", "Hip", "Thigh", "Knee", "Ankle",
              "Biceps", "Forearm", "Wrist")
response <- "BodyFat"





########################################################################################





######################################################
###   Section 2 Example: Respiratory health data   ###
######################################################


(sp.cor <- cor(BF[,varnames], method="spearman"))
sp.cor[upper.tri(sp.cor, diag=TRUE)] <- NA   # delete 'duplicate entries' in upper triangle of correlation matrix

vlarge <- mlarge <- NULL
for (i in 1:nrow(sp.cor)){
  for (j in 1:ncol(sp.cor)){
    if (!is.na(sp.cor[i,j])){
      if (abs(sp.cor[i,j])>=0.5) vlarge <- rbind(vlarge, c(rownames(sp.cor)[i], colnames(sp.cor)[j], round(sp.cor[i,j], 2)))
      if (abs(sp.cor[i,j])>=0.4 & abs(sp.cor[i,j])<0.5) mlarge <- rbind(mlarge, c(rownames(sp.cor)[i], colnames(sp.cor)[j], round(sp.cor[i,j],2)))
    }
  }
}

vstr <- apply(vlarge, 1, function(x){paste(x[1], " and ", x[2], " (", x[3], ")", sep="")})
paste("Strong correlations (>=0.5) are observed between", paste(vstr,collapse=", "))

mstr <- apply(mlarge, 1, function(x){paste(x[1], " and ", x[2], sep="")})
paste("Furthermore, the pairs", paste(mstr,collapse=", "), "show medium strong correlations (between 0.4 and 0.5).")



##############
## WebFigure 4
postscript(file="../results/WebFigure_4.eps", width=10, height=10)
corrgram(cor(BF[,varnames], method="spearman"), type="cor", order=FALSE, lower.panel=panel.shade,
  upper.panel=panel.pie, text.panel=panel.txt,
  main="")
dev.off()





########################################################################################





##########################################################
###   Section 4.1 Variable selection and predictions   ###
##########################################################


## Prepare result matrix
table1 <- matrix(NA, nrow=18, ncol=8)
rownames(table1) <- c("(Intercept)", varnames, "R2", "adj.R2", "nvar")
colnames(table1) <- c("full.beta", "full.SE", "full.p", "BE.AIC", "BE.BIC", "BE.alpha157", "BE.alpha05", "BE.alpha013")


## 'Full' model
full.model <- as.formula(paste(response, " ~ ", paste(varnames,collapse=" + ")))
(full <- lm(full.model, data=BF))
lp.full <- predict(full, type="response", se.fit=T)
table1["nvar","full.beta"] <- length(full$coefficients)-1   # subtract 1 for the intercept
table1["R2","full.beta"] <- summary(full)$r.squared
table1["adj.R2","full.beta"] <- summary(full)$adj.r.squared
table1[rownames(summary(full)$coefficients),c("full.beta", "full.SE")] <- summary(full)$coefficients[,c(1,2)]
table1[rownames(summary(full)$coefficients),"full.p"] <- summary(full)$coefficients[,4]


## BE(AIC)
(BE.AIC <- step(full,direction="backward", data=BF, k=2, trace=0))
lp.BE.AIC <- predict(BE.AIC, type="response", se.fit=T)
table1["nvar","BE.AIC"] <- length(BE.AIC$coefficients)-1   # subtract 1 for the intercept
table1["R2","BE.AIC"] <- summary(BE.AIC)$r.squared
table1["adj.R2","BE.AIC"] <- summary(BE.AIC)$adj.r.squared
table1[names(full$coefficients),"BE.AIC"] <- as.integer(names(full$coefficients)%in%names(BE.AIC$coefficients))


## BE(BIC)
(BE.BIC <- step(full,direction="backward", data=BF, k=log(nrow(BF)), trace=0))
lp.BE.BIC <- predict(BE.BIC, type="response", se.fit=T)
table1["nvar","BE.BIC"] <- length(BE.BIC$coefficients)-1   # subtract 1 for the intercept
table1["R2","BE.BIC"] <- summary(BE.BIC)$r.squared
table1["adj.R2","BE.BIC"] <- summary(BE.BIC)$adj.r.squared
table1[names(full$coefficients),"BE.BIC"] <- as.integer(names(full$coefficients)%in%names(BE.BIC$coefficients))


## Backward elimination based on p value 15.7% (corresponding to AIC)
BE.alpha157 <- select.by.alpha(response, varnames, 0.157, BF)
lp.BE.alpha157 <- predict(BE.alpha157, type="response", se.fit=T)
table1["nvar","BE.alpha157"] <- length(BE.alpha157$coefficients)-1   # subtract 1 for the intercept
table1["R2","BE.alpha157"] <- summary(BE.alpha157)$r.squared
table1["adj.R2","BE.alpha157"] <- summary(BE.alpha157)$adj.r.squared
table1[names(full$coefficients),"BE.alpha157"] <- as.integer(names(full$coefficients)%in%names(BE.alpha157$coefficients))


## Backward elimination based on p value 5%
(BE.alpha05 <- select.by.alpha(response, varnames, 0.05, BF))
lp.BE.alpha05 <- predict(BE.alpha05, type="response", se.fit=T)
table1["nvar","BE.alpha05"] <- length(BE.alpha05$coefficients)-1   # subtract 1 for the intercept
table1["R2","BE.alpha05"] <- summary(BE.alpha05)$r.squared
table1["adj.R2","BE.alpha05"] <- summary(BE.alpha05)$adj.r.squared
table1[names(full$coefficients),"BE.alpha05"] <- as.integer(names(full$coefficients)%in%names(BE.alpha05$coefficients))


## Backward elimination based on p value 1.3% (corresponding to BIC)
(p <- 1-pchisq(q=log(nrow(BF)), df=1))
(BE.alpha013 <- select.by.alpha(response, varnames, 0.013, BF))
lp.BE.alpha013 <- predict(BE.alpha013, type="response", se.fit=T)
table1["nvar","BE.alpha013"] <- length(BE.alpha013$coefficients)-1   # subtract 1 for the intercept
table1["R2","BE.alpha013"] <- summary(BE.alpha013)$r.squared
table1["adj.R2","BE.alpha013"] <- summary(BE.alpha013)$adj.r.squared
table1[names(full$coefficients),"BE.alpha013"] <- as.integer(names(full$coefficients)%in%names(BE.alpha013$coefficients))



###########
## Table 1
table1
write.csv2(data.frame(table1), file="../results/Table1.csv", eol = "\r\n")


###########
## Figure 1 
postscript(file="../results/Figure_1.eps", width=10, height=5)
ref.model <- lp.full
comp.model <- lp.BE.BIC
par(mfrow=c(1,2))
plot(ref.model$fit, comp.model$fit,type = "points", col="grey60", main="", xlab="full", ylab="BE(BIC)", xlim=xlim.sc, ylim=ylim.sc)
abline(a=0, b=1)
plot((comp.model$fit+ref.model$fit)/2, comp.model$fit-ref.model$fit, col="grey60", main="", xlab="( BE(BIC) + full ) / 2", ylab="BE(BIC) - full", xlim=xlim.ba, ylim=ylim.ba)
abline(a=0, b=0, lty=2)
lines(lowess(x=(comp.model$fit+ref.model$fit)/2, y=comp.model$fit-ref.model$fit))
dev.off()


############
## Figure 2a
se.ratio <- lp.full$se.fit/lp.BE.BIC$se.fit
postscript(file="../results/Figure_2a.eps", width=5, height=5)
par(mfrow=c(1,2), mar=mar.ba)
boxplot(se.ratio, ylab="SE( full ) / SE( BE(BIC) )", ylim=ylim.ser)
points(1, mean(se.ratio), pch=4, lwd=2)
lines(c(0.95,1.05), quantile(se.ratio, probs=rep(0.05, 2)))
text(x=1.05, y=quantile(se.ratio,0.05), label="5% quantile", pos=4)
lines(c(0.95,1.05), quantile(se.ratio, probs=rep(0.95, 2)))
text(x=1.05, y=quantile(se.ratio,0.95), label="95% quantile", pos=4)
dev.off()





########################################################################################





####################################################
###   Section 4.2 Variable selection stability   ###
####################################################


## General settings for bootstrapMA
nboot <- 500   # how many boostrap replications
nxvar <- length(varnames)
alpha <- 0.05

full.f <- full.model
lower.f <- as.formula(paste(response, "~ 1"))

build <- BF[, c(response, varnames)]
test <- BF[, c(response, varnames)]


## Bootstrap indices
set.seed(538785)
build.boot <- boot(build, corr, R=nboot, stype="w")
boot.index <- boot.array(build.boot,indices=T)

set.seed(863786)
build.boot <- boot(build, corr, R=nboot, stype="w")
boot.index2 <- boot.array(build.boot,indices=T)


## Bootstrap variable selection AIC
print("Bootstrap variable selection AIC")
boot.AIC <- screening.step.glm(dset=build, newdata=NA, fam="gaussian",   
                               predtype="response", index=boot.index, varnames=varnames,
                               full.formula=full.f, lower.formula=lower.f,
                               crit.penal=2, direct="backward",
                               alpha.c=alpha, alpha.p=alpha, nval=1)
boot.AIC$screen


## Bootstrap variable selection BIC
print("Bootstrap variable selection BIC")
boot.BIC <- screening.step.glm(dset=build, newdata=NA, fam="gaussian",
                               predtype="response", index=boot.index, varnames=varnames,
                               full.formula=full.f, lower.formula=lower.f,
                               crit.penal=log(ncol(boot.index)), direct="backward",
                               alpha.c=alpha, alpha.p=alpha, nval=1)
boot.BIC$screen





########################################################################################





#########################################
###   Section 4.3 Model uncertainty   ###
#########################################


## Settings for bootstrapMA  
nboot <- 500   # number of bootstrap replications
nxvar <- length(varnames)
alpha <- 0.05

full.f <- full.model
lower.f <- as.formula(paste(response, "~ 1"))

build <- BF[, c(response, varnames)]
test <- BF[, c(response, varnames)]


## Bootstrap indices
set.seed(538785)
build.boot <- boot(build, corr, R=nboot, stype="w")
boot.index <- boot.array(build.boot,indices=T)

set.seed(863786)
build.boot <- boot(build, corr, R=nboot, stype="w")
boot.index2 <- boot.array(build.boot,indices=T)


## BootstrapMA with BE(AIC) and C=40%
print("Bootstrap model averaging AIC, 40% level")
models.BIC.40 <- fma.step.glm(dset=build, newdata=test, fam="gaussian", predtype="response",
                              index=boot.index2, varnames, full.f, lower.f, crit.penal=log(n),
                              direct="backward", screen=T, screen.lev=0.4, screen.res=boot.BIC,
                              alpha.c=alpha, alpha.p=alpha, nval=1)
preds.BIC.40 <- fma.pred.glm(dset=build, newdata=test, fam="gaussian", obj=models.BIC.40, best=F,
                             se.err=lp.full$residual.scale, alpha.c=alpha, alpha.p=alpha)



## BootstrapMA with BE(AIC) and C=30%
print("Bootstrap model averaging AIC, 30% level")
models.BIC.30 <- fma.step.glm(dset=build, newdata=test, fam="gaussian", predtype="response",
                              index=boot.index2, varnames, full.f, lower.f, crit.penal=log(n),
                              direct="backward", screen=T, screen.lev=0.3, screen.res=boot.BIC,
                              alpha.c=alpha, alpha.p=alpha, nval=1)
preds.BIC.30 <- fma.pred.glm(dset=build, newdata=test, fam="gaussian", obj=models.BIC.30, best=F,
                             se.err=lp.full$residual.scale, alpha.c=alpha, alpha.p=alpha)


## BootstrapMA with BE(BIC) and C=30%
print("Bootstrap model averaging BIC, 30% level")
models.BIC.30 <- fma.step.glm(dset=build, newdata=test, fam="gaussian", predtype="response",
                              index=boot.index2, varnames, full.f, lower.f, crit.penal=log(ncol(boot.index)),
                              direct="backward", screen=T, screen.lev=0.3, screen.res=boot.BIC,
                              alpha.c=alpha, alpha.p=alpha, nval=1)
preds.BIC.30 <- fma.pred.glm(dset=build, newdata=test, fam="gaussian", obj=models.BIC.30, best=F,
                             se.err=lp.full$residual.scale, alpha.c=alpha, alpha.p=alpha)



##########
## Table 2
nrow(models.BIC.40$m)  # no. of different models selected in at least on bootstrap sample
top10 <- models.BIC.40$m[1:10,-1]

table2 <- cbind(round(100*models.BIC.40$screen,1), t(top10))
colnames(table2) <- c("BIF_BIC", paste("M", 1:10, sep=""))
o <- rev(order(table2[,1]))
table2 <- table2[o,]

hMj.C40 <- round(100*models.BIC.40$m[1:10,1], 1)
table2 <- rbind(table2, c(NA, hMj.C40))
rownames(table2)[nrow(table2)] <- "h(Mj): BIC, C=40%"

hMj.C30 <- round(100*models.BIC.30$m[1:10,1], 1)
table2 <- rbind(table2, c(NA, hMj.C30))   # add bootstrap model selection frequencies (step 2) of bootstrapMA-AIC (30%)
## Note: The 10 top models are NOT identical to those selected with C=40%
rownames(table2)[nrow(table2)] <- "h(Mj): BIC, C=30%"

table2
write.csv2(data.frame(table2), file="../results/Table2.csv", eol = "\r\n")


##########
## Top model using BIC, C=30% as described in text
dim(models.BIC.40$m)
top10 <- models.BIC.40$m[1:10,-1]

extra.tab <- cbind(round(100*models.BIC.40$screen,1), t(top10))
colnames(extra.tab) <- c("BIF_BIC", paste("M", 1:10, sep=""))
o <- rev(order(extra.tab[,1]))
extra.tab <- extra.tab[o,]

hMj.BIC.C40 <- round(100*models.BIC.40$m[1:10,1], 1)
extra.tab <- rbind(extra.tab, c(NA, hMj.BIC.C40))
rownames(extra.tab)[nrow(extra.tab)] <- "h(Mj): BIC, C=40%"

extra.tab



###############################################



## BootstrapMA with 20-fold cross-validation


## Cross-validation specification
K <- 20   ## number of crossvalidations
set.seed(874568)
cv.index  <- cross.split(nrow(BF),K)
set.seed(54442)
cv.seedvec1 <- sample(1:1000000, K)
cv.seedvec2 <- sample(1:1000000, K)


## preparation of matrices 
dnames <- c("full", "BE.BIC", "bootMA.BIC.40")
cv.prediction <- matrix(NA,nrow=nrow(BF), ncol=length(dnames))
colnames(cv.prediction) <- dnames
cv.se <- matrix(NA,nrow=nrow(BF), ncol=length(dnames))
colnames(cv.se) <- dnames

for (j in 1:K){
  
  print(paste("j:",j))

  cv.test.ind <- which(cv.index==j)
  cv.build <- BF[-cv.test.ind, c(response, varnames)]
  cv.test <- BF[cv.test.ind, c(response, varnames)]

  
  ## full model
  print("full model")
  cv.full <- glm(full.f, family="gaussian", data=cv.build)
  cv.lp.full <- predict(cv.full, newdata=cv.test, type="response", se.fit=T)
  cv.prediction[cv.test.ind, "full"] <- cv.lp.full$fit
  cv.se[cv.test.ind, "full"] <- cv.lp.full$se.fit

    
  ## BE(BIC)
  print("BE(BIC)")
  cv.BE.BIC <- step(cv.full, scope = list(upper = full.f, lower = lower.f),
                     k=log(nrow(cv.build)), trace=0, direction="backward")
  cv.lp.BE.BIC <- predict(cv.BE.BIC, newdata=cv.test, type="response", se.fit=T)
  cv.prediction[cv.test.ind, "BE.BIC"] <- cv.lp.BE.BIC$fit
  cv.se[cv.test.ind, "BE.BIC"] <- cv.lp.BE.BIC$se.fit
  

  ## bootstrap matrices
  set.seed(cv.seedvec1[j])
  cv.build.boot <- boot(cv.build, corr, R=nboot, stype="w")
  cv.index1 <- boot.array(cv.build.boot,indices=T)
  set.seed(cv.seedvec2[j])
  cv.build.boot <- boot(cv.build, corr, R=nboot, stype="w")
  cv.index2 <- boot.array(cv.build.boot,indices=T)
  
  ## Bootstrap variable selection with BE(BIC)
  print("Bootstrap variable selection BIC")
  cv.v.BIC <- screening.step.glm(dset=cv.build, newdata=NA, fam="gaussian",
                                 predtype="response", index=cv.index1, varnames=varnames,
                                 full.formula=full.f, lower.formula=lower.f,
                                 crit.penal=log(nrow(cv.build)), direct="backward",
                                 alpha.c=alpha, alpha.p=alpha, nval=1)

  ## BootstrapMA BIC 40%
  print("Bootstrap model averaging BIC, 40% level")
  cv.models.bootMA.BIC.40 <- fma.step.glm(dset=cv.build, newdata=cv.test, fam="gaussian", predtype="response",
                                          index=cv.index2, varnames, full.f, lower.f, crit.penal=log(nrow(cv.build)),
                                          direct="backward", screen=T, screen.lev=0.4, screen.res=cv.v.BIC,
                                          alpha.c=alpha, alpha.p=alpha, nval=1)
  cv.lp.bootMA.BIC.40 <- fma.pred.glm(dset=cv.build, newdata=cv.test, fam="gaussian", obj=cv.models.bootMA.BIC.40, best=F,
                                      se.err=cv.lp.full$residual.scale, alpha.c=alpha, alpha.p=alpha)
  cv.prediction[cv.test.ind, "bootMA.BIC.40"] <- cv.lp.bootMA.BIC.40$pred
  cv.se[cv.test.ind, "bootMA.BIC.40"] <- cv.lp.bootMA.BIC.40$se.fit

}

## Figure 2a
se.ratio <- lp.full$se.fit/lp.BE.BIC$se.fit
postscript(file="../results/Figure_2a.eps", width=5, height=5)
par(mfrow=c(1,2), mar=mar.ba)
boxplot(se.ratio, ylab="SE( full ) / SE( BE(BIC) )", ylim=ylim.ser)
points(1, mean(se.ratio), pch=4, lwd=2)
lines(c(0.95,1.05), quantile(se.ratio, probs=rep(0.05, 2)))
text(x=1.05, y=quantile(se.ratio,0.05), label="5% quantile", pos=4)
lines(c(0.95,1.05), quantile(se.ratio, probs=rep(0.95, 2)))
text(x=1.05, y=quantile(se.ratio,0.95), label="95% quantile", pos=4)
dev.off()


############
## Figure 2b
se.ratio <- cv.se[,"full"]/cv.se[,"BE.BIC"]
summary(se.ratio)
postscript(file="../results/Figure_2b.eps", width=5, height=5)
par(mfrow=c(1,2), mar=mar.ba)
boxplot(se.ratio, ylab="SE( full ) / SE( BE(BIC) )", ylim=ylim.ser)
points(1, mean(se.ratio), pch=4, lwd=2)
lines(c(0.95,1.05), quantile(se.ratio, probs=rep(0.05, 2)))
text(x=1.05, y=quantile(se.ratio,0.05), label="5% quantile", pos=4)
lines(c(0.95,1.05), quantile(se.ratio, probs=rep(0.95, 2)))
text(x=1.05, y=quantile(se.ratio,0.95), label="95% quantile", pos=4)
dev.off()


###########
## Figure 3
postscript(file="../results/Figure_3a.eps", width=10, height=5)
ref.model <- list(fit=cv.prediction[,"BE.BIC"], se.fit=cv.se[,"BE.BIC"])
comp.model <- list(fit=cv.prediction[,"bootMA.BIC.40"], se.fit=cv.se[,"bootMA.BIC.40"])
par(mfrow=c(1,2), mar=mar.ba)
plot(ref.model$fit, comp.model$fit, col="grey60", main="", xlab="BE(BIC)", ylab="bootMA(BIC, C=40%)", xlim=xlim.sc, ylim=ylim.sc)
abline(a=0, b=1)
plot((comp.model$fit+ref.model$fit)/2, comp.model$fit-ref.model$fit, col="grey60", main="", xlab="( bootMA(BIC, C=40%) + BE(BIC) ) / 2", ylab="bootMA(BIC, C=40%) - BE(BIC)", xlim=xlim.ba, ylim=ylim.ba)
abline(a=0, b=0, lty=2)
lines(lowess(x=(comp.model$fit+ref.model$fit)/2, y=comp.model$fit-ref.model$fit))
dev.off()

se.ratio <- cv.se[,"bootMA.BIC.40"]/cv.se[,"BE.BIC"]
summary(se.ratio)
postscript(file="../results/Figure_3b_right.eps", width=5, height=5)
par(mfrow=c(1,2), mar=mar.ba)
boxplot(se.ratio, ylab="SE( bootMA(BIC, C=40%) ) / SE( BE(BIC) )", ylim=ylim.ser)
points(1, mean(se.ratio), pch=4, lwd=2)
lines(c(0.95,1.05), quantile(se.ratio, probs=rep(0.05, 2)))
text(x=1.05, y=quantile(se.ratio,0.05), label="5% quantile", pos=4)
lines(c(0.95,1.05), quantile(se.ratio, probs=rep(0.95, 2)))
text(x=1.05, y=quantile(se.ratio,0.95), label="95% quantile", pos=4)
dev.off()

postscript(file="../results/Figure_3b_left.eps", width=5, height=5)
comp.model <- list(fit=cv.prediction[,"BE.BIC"], se.fit=cv.se[,"BE.BIC"])
ref.model <- list(fit=cv.prediction[,"bootMA.BIC.40"], se.fit=cv.se[,"bootMA.BIC.40"])
par(mar=mar.ba)
plot(ref.model$se.fit, comp.model$se.fit, col="grey60", main="", xlab="bootMA(BIC, C=40%)", ylab="BE(BIC)",xlim = c(0,3),ylim = c(0,2))
abline(a=0, b=1)
text(1, 1, "angle bisector", srt=58, cex=1.1)
dev.off()


##############
## WebFigure 1
postscript(file="../results/WebFigure_1.eps", width=10, height=5)
ref.model <- list(fit=cv.prediction[,"full"], se.fit=cv.se[,"full"])
comp.model <- list(fit=cv.prediction[,"BE.BIC"], se.fit=cv.se[,"BE.BIC"])
par(mfrow=c(1,2), mar=mar.ba)
plot(ref.model$fit, comp.model$fit, col="grey60", main="", xlab="full", ylab="BE(BIC)", xlim=xlim.sc, ylim=ylim.sc)
abline(a=0, b=1)
plot((comp.model$fit+ref.model$fit)/2, comp.model$fit-ref.model$fit, col="grey60", main="", xlab="( BE(BIC) + full ) / 2", ylab="BE(BIC) - full", xlim=xlim.ba, ylim=ylim.ba)
abline(a=0, b=0, lty=2)
lines(lowess(x=(comp.model$fit+ref.model$fit)/2, y=comp.model$fit-ref.model$fit))
dev.off()





########################################################################################





########################################
###   Section 4.4 Stronger effects   ###
########################################


#############
## WebTable 3

dom.v <- c("Abdomen", "Sex")
change.v <- c("Hip", "Height", "Forearm", "Knee")
dom.cor <- matrix(NA, nrow=length(dom.v), ncol=length(change.v))
rownames(dom.cor) <- dom.v
colnames(dom.cor) <- change.v

for (i in 1:length(dom.v)){
  for (j in 1:length(change.v)){
    dom.cor[i, j] <- round(cor(BF[,dom.v[i]], BF[,change.v[j]], method="spearman"),2)
  }
}
dom.cor

write.csv2(data.frame(dom.cor), file="../results/WebTable_3.csv", eol = "\r\n")



#############################################################
## Variable selection


dom.resmat <- dom.full <- dom.lp.full <- dom.BE.AIC <- dom.lp.BE.AIC <- dom.BE.BIC <- dom.lp.BE.BIC <- list(NULL)


## Fit models without dominating variables (one at a time)

for (x in c("Abdomen", "Sex")){ 

  ## General settings
  tmp.varnames <- varnames[-which(varnames==x)]
  tmp.full.model <- as.formula(paste(response, " ~ ", paste(tmp.varnames,collapse=" + ")))

  ## Prepare result matrix
  dom.resmat[[x]] <- matrix(NA, nrow=17, ncol=5)
  rownames(dom.resmat[[x]]) <- c("(Intercept)", tmp.varnames, "R2", "adj.R2", "nvar")
  colnames(dom.resmat[[x]]) <- c("full.beta", "full.SE", "full.p", "BE.AIC", "BE.BIC")

  ## 'Full' model
  (dom.full[[x]] <- lm(tmp.full.model, data=BF))
  dom.lp.full[[x]] <- predict(dom.full[[x]], type="response", se.fit=T)
  dom.resmat[[x]]["nvar","full.beta"] <- length(dom.full[[x]]$coefficients)-1   # subtract 1 for the intercept
  dom.resmat[[x]]["R2","full.beta"] <- summary(dom.full[[x]])$r.squared
  dom.resmat[[x]]["adj.R2","full.beta"] <- summary(dom.full[[x]])$adj.r.squared
  tmp <- summary(dom.full[[x]])$coefficients
  dom.resmat[[x]][rownames(tmp),1:2] <- round(tmp[,c(1,2)],2)
  dom.resmat[[x]][rownames(tmp),3] <- round(tmp[,4],4)

  ## BE(AIC)
  (dom.BE.AIC[[x]] <- step(dom.full[[x]],direction="backward", data=BF, k=2, trace=0))
  dom.lp.BE.AIC[[x]] <- predict(dom.BE.AIC[[x]], type="response", se.fit=T)
  dom.resmat[[x]]["nvar","BE.AIC"] <- length(dom.BE.AIC[[x]]$coefficients)-1   # subtract 1 for the intercept
  dom.resmat[[x]]["R2","BE.AIC"] <- summary(dom.BE.AIC[[x]])$r.squared
  dom.resmat[[x]]["adj.R2","BE.AIC"] <- summary(dom.BE.AIC[[x]])$adj.r.squared
  dom.resmat[[x]][names(dom.full[[x]]$coefficients),"BE.AIC"] <- as.integer(names(dom.full[[x]]$coefficients)%in%names(dom.BE.AIC[[x]]$coefficients))

  ## BE(BIC)
  n <- nrow(BF)
  (dom.BE.BIC[[x]] <- step(dom.full[[x]],direction="backward", data=BF, k=log(n), trace=0))
  dom.lp.BE.BIC[[x]] <- predict(dom.BE.BIC[[x]], type="response", se.fit=T)
  dom.resmat[[x]]["nvar","BE.BIC"] <- length(dom.BE.BIC[[x]]$coefficients)-1   # subtract 1 for the intercept
  dom.resmat[[x]]["R2","BE.BIC"] <- summary(dom.BE.BIC[[x]])$r.squared
  dom.resmat[[x]]["adj.R2","BE.BIC"] <- summary(dom.BE.BIC[[x]])$adj.r.squared
  dom.resmat[[x]][names(dom.full[[x]]$coefficients),"BE.BIC"] <- as.integer(names(dom.full[[x]]$coefficients)%in%names(dom.BE.BIC[[x]]$coefficients))
 
  ## Complete result matrix
  dom.resmat[[x]][16:17,] <- round(dom.resmat[[x]][16:17,], 2)
  print(dom.resmat[[x]])
  
}


##################################################################
## Fit models without all two dominating variables (together) 

## General settings
x <- "wo2"
tmp.varnames <- varnames[-which(varnames%in%c("Abdomen", "Sex"))]
tmp.full.model <- as.formula(paste(response, " ~ ", paste(tmp.varnames,collapse=" + ")))

## Prepare result matrix
dom.resmat[[x]] <- matrix(NA, nrow=16, ncol=5)
rownames(dom.resmat[[x]]) <- c("(Intercept)", tmp.varnames, "R2", "adj.R2", "nvar")
colnames(dom.resmat[[x]]) <- c("full.beta", "full.SE", "full.p", "BE.AIC", "BE.BIC")

## 'Full' model
(dom.full[[x]] <- lm(tmp.full.model, data=BF))
dom.lp.full[[x]] <- predict(dom.full[[x]], type="response", se.fit=T)
dom.resmat[[x]]["nvar","full.beta"] <- length(dom.full[[x]]$coefficients)-1   # subtract 1 for the intercept
dom.resmat[[x]]["R2","full.beta"] <- summary(dom.full[[x]])$r.squared
dom.resmat[[x]]["adj.R2","full.beta"] <- summary(dom.full[[x]])$adj.r.squared
tmp <- summary(dom.full[[x]])$coefficients
dom.resmat[[x]][rownames(tmp),1:2] <- round(tmp[,c(1,2)],2)
dom.resmat[[x]][rownames(tmp),3] <- round(tmp[,4],4)

## BE(AIC)
(dom.BE.AIC[[x]] <- step(dom.full[[x]],direction="backward", data=BF, k=2, trace=0))
dom.lp.BE.AIC[[x]] <- predict(dom.BE.AIC[[x]], type="response", se.fit=T)
dom.resmat[[x]]["nvar","BE.AIC"] <- length(dom.BE.AIC[[x]]$coefficients)-1   # subtract 1 for the intercept
dom.resmat[[x]]["R2","BE.AIC"] <- summary(dom.BE.AIC[[x]])$r.squared
dom.resmat[[x]]["adj.R2","BE.AIC"] <- summary(dom.BE.AIC[[x]])$adj.r.squared
dom.resmat[[x]][names(dom.full[[x]]$coefficients),"BE.AIC"] <- as.integer(names(dom.full[[x]]$coefficients)%in%names(dom.BE.AIC[[x]]$coefficients))

## BE(BIC)
n <- nrow(BF)
(dom.BE.BIC[[x]] <- step(dom.full[[x]],direction="backward", data=BF, k=log(n), trace=0))
dom.lp.BE.BIC[[x]] <- predict(dom.BE.BIC[[x]], type="response", se.fit=T)
dom.resmat[[x]]["nvar","BE.BIC"] <- length(dom.BE.BIC[[x]]$coefficients)-1   # subtract 1 for the intercept
dom.resmat[[x]]["R2","BE.BIC"] <- summary(dom.BE.BIC[[x]])$r.squared
dom.resmat[[x]]["adj.R2","BE.BIC"] <- summary(dom.BE.BIC[[x]])$adj.r.squared
dom.resmat[[x]][names(dom.full[[x]]$coefficients),"BE.BIC"] <- as.integer(names(dom.full[[x]]$coefficients)%in%names(dom.BE.BIC[[x]]$coefficients))

## Complete result matrix
dom.resmat[[x]][14:16,] <- round(dom.resmat[[x]][14:16,], 3)
print(dom.resmat[[x]])



#############################
## Fit separate models by sex 

## General settings
tmp.data <- list("male"=BF[BF$Sex==0,], "female"=BF[BF$Sex==1,])
tmp.varnames <- varnames[-which(varnames=="Sex")]
tmp.full.model <- as.formula(paste(response, " ~ ", paste(tmp.varnames,collapse=" + ")))
x <- "by.Sex"
dom.resmat[[x]] <- dom.full[[x]] <- dom.lp.full[[x]] <- dom.BE.AIC[[x]] <- dom.lp.BE.AIC[[x]] <- dom.BE.BIC[[x]] <- dom.lp.BE.BIC[[x]] <- list(NULL)

## Prepare result matrix
dom.resmat[[x]] <- matrix(NA, nrow=18, ncol=10)
rownames(dom.resmat[[x]]) <- c("sex (0=male/1=female)", "(Intercept)", tmp.varnames, "R2", "adj.R2", "nvar")
colnames(dom.resmat[[x]]) <- c("full.beta0", "full.beta1", "full.SE0", "full.SE1", "full.p0", "full.p1", "BE.AIC0", "BE.AIC1", "BE.BIC0", "BE.BIC1")
dom.resmat[[x]][1,] <- rep(c(0,1), 5)

for (i in 0:1){   # loop over male/female
 
  ## 'Full' model
  (dom.full[[x]][[i+1]] <- lm(tmp.full.model, data=tmp.data[[i+1]]))
  dom.lp.full[[x]][[i+1]] <- predict(dom.full[[x]][[i+1]], type="response", se.fit=T)
  dom.resmat[[x]]["nvar",paste("full.beta", i, sep="")] <- length(dom.full[[x]][[i+1]]$coefficients)-1   # subtract 1 for the intercept
  dom.resmat[[x]]["R2",paste("full.beta", i, sep="")] <- summary(dom.full[[x]][[i+1]])$r.squared
  dom.resmat[[x]]["adj.R2",paste("full.beta", i, sep="")] <- summary(dom.full[[x]][[i+1]])$adj.r.squared
  tmp <- summary(dom.full[[x]][[i+1]])$coefficients
  dom.resmat[[x]][rownames(tmp),paste(c("full.beta", "full.SE"), i, sep="")] <- round(tmp[,c(1,2)],2)
  dom.resmat[[x]][rownames(tmp),paste("full.p", i, sep="")] <- round(tmp[,4],4)

  ## BE(AIC)
  (dom.BE.AIC[[x]][[i+1]] <- step(dom.full[[x]][[i+1]],direction="backward", data=tmp.data[[i+1]], k=2, trace=0))
  dom.lp.BE.AIC[[x]][[i+1]] <- predict(dom.BE.AIC[[x]][[i+1]], type="response", se.fit=T)
  dom.resmat[[x]]["nvar",paste("BE.AIC", i, sep="")] <- length(dom.BE.AIC[[x]][[i+1]]$coefficients)-1   # subtract 1 for the intercept
  dom.resmat[[x]]["R2",paste("BE.AIC", i, sep="")] <- summary(dom.BE.AIC[[x]][[i+1]])$r.squared
  dom.resmat[[x]]["adj.R2",paste("BE.AIC", i, sep="")] <- summary(dom.BE.AIC[[x]][[i+1]])$adj.r.squared
  dom.resmat[[x]][names(dom.full[[x]][[i+1]]$coefficients),paste("BE.AIC", i, sep="")] <- as.integer(names(dom.full[[x]][[i+1]]$coefficients)%in%names(dom.BE.AIC[[x]][[i+1]]$coefficients))

  ## BE(BIC)
  n <- nrow(tmp.data[[i+1]])
  (dom.BE.BIC[[x]][[i+1]] <- step(dom.full[[x]][[i+1]], direction="backward", data=tmp.data[[i+1]], k=log(n), trace=0))
  dom.lp.BE.BIC[[x]][[i+1]] <- predict(dom.BE.BIC[[x]][[i+1]], type="response", se.fit=T)
  dom.resmat[[x]]["nvar",paste("BE.BIC", i, sep="")] <- length(dom.BE.BIC[[x]][[i+1]]$coefficients)-1   # subtract 1 for the intercept
  dom.resmat[[x]]["R2",paste("BE.BIC", i, sep="")] <- summary(dom.BE.BIC[[x]][[i+1]])$r.squared
  dom.resmat[[x]]["adj.R2",paste("BE.BIC", i, sep="")] <- summary(dom.BE.BIC[[x]][[i+1]])$adj.r.squared
  dom.resmat[[x]][names(dom.full[[x]][[i+1]]$coefficients),paste("BE.BIC", i, sep="")] <- as.integer(names(dom.full[[x]][[i+1]]$coefficients)%in%names(dom.BE.BIC[[x]][[i+1]]$coefficients))

}  

## Complete result matrix
dom.resmat[[x]][17:18,] <- round(dom.resmat[[x]][17:18,], 3)
print(dom.resmat[[x]])



###########
## Figure 4
postscript(file="../results/Figure_4a.eps", width=10, height=5)
ref.model <- lp.full
comp.model <- dom.lp.full$Sex
par(mar=mar.ba, mfrow=c(1,2))
plot(ref.model$fit, comp.model$fit, col="grey60", main="", xlab="with Sex", ylab="w/o Sex", xlim=xlim.sc, ylim=ylim.sc)
abline(a=0, b=1)
plot((comp.model$fit+ref.model$fit)/2, comp.model$fit-ref.model$fit, col="grey60", main="", xlab="( with Sex + w/o Sex ) / 2", ylab="w/o Sex - with Sex", xlim=xlim.ba, ylim=ylim.ba)
abline(a=0, b=0, lty=2)
lines(lowess(x=(comp.model$fit+ref.model$fit)/2, y=comp.model$fit-ref.model$fit))
dev.off()

postscript(file="../results/Figure_4b.eps", width=10, height=5)
ref.model <- lp.BE.BIC
comp.model <- dom.lp.BE.BIC$Sex
par(mar=mar.ba, mfrow=c(1,2))
plot(ref.model$fit, comp.model$fit, col="grey60", main="", xlab="with Sex", ylab="w/o Sex", xlim=xlim.sc, ylim=ylim.sc)
abline(a=0, b=1)
plot((comp.model$fit+ref.model$fit)/2, comp.model$fit-ref.model$fit, col="grey60", main="", xlab="( with Sex + w/o Sex ) / 2", ylab="w/o Sex - with Sex", xlim=xlim.ba, ylim=ylim.ba)
abline(a=0, b=0, lty=2)
lines(lowess(x=(comp.model$fit+ref.model$fit)/2, y=comp.model$fit-ref.model$fit))
dev.off()



#################################################

##  Bootstrapping without Sex, Abdomen

## General settings
tmp.varnames <- varnames[-which(varnames%in%c("Abdomen", "Sex"))]
tmp.full.f <- as.formula(paste(response, " ~ ", paste(tmp.varnames,collapse=" + ")))
tmp.lower.f <- as.formula(paste(response, "~ 1"))
domboot.nxvar<- length(tmp.varnames)

## Generate bootstrap samples
set.seed(538785)
domboot.build.boot <- boot(build, corr, R=nboot, stype="w")
domboot.index1 <- boot.array(domboot.build.boot,indices=T)

## BIFs
print("Bootstrap variable selection AIC")
domboot.v.AIC <- screening.step.glm(dset=build, newdata=NA, fam="gaussian",
                                predtype="response", index=domboot.index1, varnames=tmp.varnames,
                                full.formula=tmp.full.f, lower.formula=tmp.lower.f,
                                crit.penal=2, direct="backward",
                                alpha.c=alpha, alpha.p=alpha, nval=1)
rev(sort(domboot.v.AIC$screen))

print("Bootstrap variable selection BIC")
domboot.v.BIC <- screening.step.glm(dset=build, newdata=NA, fam="gaussian",
                                predtype="response", index=domboot.index1, varnames=varnames,
                                full.formula=tmp.full.f, lower.formula=tmp.lower.f,
                                crit.penal=log(ncol(domboot.index1)), direct="backward",
                                alpha.c=alpha, alpha.p=alpha, nval=1)
rev(sort(domboot.v.BIC$screen))


##########
## Table 3
table3 <- matrix(NA, nrow=length(varnames), ncol=4)
rownames(table3) <- varnames
colnames(table3) <- c("AIC.14", "AIC.12", "BIC.14", "BIC.12")
table3[names(boot.AIC$screen),"AIC.14"] <- round(boot.AIC$screen*100, 1)   # BIFs from Section 4.2 
table3[names(domboot.v.AIC$screen),"AIC.12"] <- round(domboot.v.AIC$screen*100, 1)
table3[names(boot.BIC$screen),"BIC.14"] <- round(boot.BIC$screen*100, 1)   # BIFs from Section 4.2 
table3[names(domboot.v.BIC$screen),"BIC.12"] <- round(domboot.v.BIC$screen*100, 1)
o <- rev(order(table3[,1]))
(table3 <- table3[o,])
write.csv2(data.frame(table3), file="../results/Table3.csv", eol = "\r\n")


#############
## WebTable 2
rn <- c("all 14 variables", "w/o Abdomen", "w/o Sex", "w/o all 2")
cn <- c("full", "BE(AIC)", "BE(BIC)")

r.sq.all14 <- c(summary(full)$r.squared, summary(BE.AIC)$r.squared, summary(BE.BIC)$r.squared)
r.sq.woAbdomen <- c(summary(dom.full$Abdomen)$r.squared, summary(dom.BE.AIC$Abdomen)$r.squared, summary(dom.BE.BIC$Abdomen)$r.squared)
r.sq.woSex <- c(summary(dom.full$Sex)$r.squared, summary(dom.BE.AIC$Sex)$r.squared, summary(dom.BE.BIC$Sex)$r.squared)
r.sq.wo2 <- c(summary(dom.full$wo2)$r.squared, summary(dom.BE.AIC$wo2)$r.squared, summary(dom.BE.BIC$wo2)$r.squared)
r.sq.matrix <- rbind(r.sq.all14, r.sq.woAbdomen, r.sq.woSex, r.sq.wo2)
rownames(r.sq.matrix) <- rn
colnames(r.sq.matrix) <- cn

r.sq.change.raw <- t(apply(r.sq.matrix[-1,], 1, function(x){x - r.sq.matrix[1,]}))
r.sq.change.pc <- t(apply(r.sq.change.raw, 1, function(x){ (100/r.sq.matrix[1,]) * x}))
rownames(r.sq.change.raw) <- rownames(r.sq.change.pc) <- rn[-1]
colnames(r.sq.change.raw) <- colnames(r.sq.change.pc) <- cn

r.sq.matrix
r.sq.change.pc





########################################################################################


