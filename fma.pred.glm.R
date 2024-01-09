fma.pred.glm <- function(dset, newdata=NA, fam="binomial",
                         predtype="response", obj, best=F, se.err=1,
                         alpha.c=0.05, alpha.p=0.05) {

## ----------------------------------------------------------------------------
## Title: fma.pred.glm
## ----------------------------------------------------------------------------
## Authors: Willi Sauerbrei, Nicole Augustin, Norbert Hollaender, Anika Buchholz
## Center for Medical Biometry and Medical Informatics
## Stefan-Meier-Strasse 26, D-79104 Freiburg 
## http://www.imbi.uni-freiburg.de 
## ----------------------------------------------------------------------------
## Description: Predictions and confidence intervals for (generalised)
##              linear models as well as prediction intervals for
##              linear models. 
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: fma.pred.glm(dset, newdata=NA, fam="binomial",
##                     predtype="response", obj, best=F, se.err=1,
##                     alpha.c=0.05, alpha.p=0.05)
##
##        dset:      build data set
##        newdata:   new data to be predicted
##        fam:    character string specifying the regression
##        predtype:  the type of prediction required (see help to
##                   'predict.glm')
##        obj:       object derived by function 'fma.glm'
##        best:      logical value indicating whether the best
##                   model of the bootstrap procedure shall be
##                   used for prediction (best=T) or all models
##                   shall be used for model averaging (best=F,
##                   default)
##        se.err:    standard error used for calculation of prediction
##                   intervals (for linear models only)
##        alpha.c:   alpha level for confidence interval (to be
##                   used with generated, i.e. simulated, data),
##                   the default is alpha.c=0.05    
##        alpha.p:   alpha level for confidence interval (to be
##                   used with real data), the default is alpha.p=0.05    
## ----------------------------------------------------------------------------
## Value: pred:   prediction (of type 'predtype')
##        ci.vec: confidence interval for pred
##        pintv:  prediction interval for pred (linear regression only)
## ----------------------------------------------------------------------------


   if (best){
   formu <- obj$formu[1] 
   ### if best model of bootstrap procedure shall be used instead of
   ### model averaging
   form <- formu[1]
   fit <- glm(as.formula(form), data=dset, family=fam)
   if (length(fit$coef)==0) warning("Null model fitted!")
   if (!is.null(fit$coef)){
     if (!is.null(dim(newdata)))
       pred <- predict(fit,newdata=newdata,type=predtype, se.fit=T)
     else pred<-predict(fit,type=predtype,se.fit=T)
   }else{
     if (!is.null(dim(newdata)))
        pred <- rep(0, nrow(newdata))
      else pred <- rep(0, ncol(dset))
      warning("Null model fitted in best model. (prediction set to zero)")
   }
   
   ###
   ### obtain analytic confidence interval
   ### changes this to pred$fit -/+ 1.959964 * sqrt(1*sigma^2 +pred$se.fit^2) for
   ### prediction interval in the linear case and
   ### pred$fit -/+ 1.959964 * sqrt(1 +pred$se.fit^2) for logistic model (dispersion
   ### parameter is one).

   ## 1-alpha.c confidence interval
   ci.vec<-matrix(NA,nrow=length(pred$fit), ncol=2)
   ci.vec[,1] <- pred$fit - qnorm(1-(alpha.c/2)) * pred$se.fit
   ci.vec[,2] <- pred$fit + qnorm(1-(alpha.c/2)) * pred$se.fit

   ## calculation of prediction intervals for linear regression if wanted
   if (fam=="gaussian"){
     ## 1-alpha.p prediction interval
     pci.vec<-matrix(NA,nrow=length(pred$fit), ncol=2)
     pci.vec[,1] <- pred$fit - qnorm(1-(alpha.p/2)) * sqrt(pred$se.fit^2 + se.err^2)
     pci.vec[,2] <- pred$fit + qnorm(1-(alpha.p/2)) * sqrt(pred$se.fit^2 + se.err^2)
     return(list(pred=pred$fit, se.fit=pred$se.fit, ci=ci.vec, pintv=pci.vec))
   }
   else return(list(pred=pred$fit, se.fit=pred$se.fit, ci=ci.vec))
 } # end if(best)
 else{
   ### model averaging
   form <- obj$formu
   if (is.null(dim(newdata))) lpmat<-varmat<-matrix(0, nrow=nrow(dset), ncol=length(form))
   else lpmat<-varmat<-matrix(0, nrow=nrow(newdata), ncol=length(form))
   for (i in 1:length(form)){
     fitw <- glm(as.formula(form[i]), data=dset, family=fam)
     if (!is.null(fitw$coef)){ ## null model cannot be predicted
                               ## -> set prediction to zero
        if (!is.null(dim(newdata))) tmp <- predict(fitw, newdata=newdata, type=predtype, se.fit=T)
        else tmp <- predict(fitw, type=predtype, se.fit=T)
        lpmat[,i]<-tmp$fit
        varmat[,i]<-tmp$se.fit^2
      }else{
        lpmat[,i] <- 0
        varmat[,i] <- 0
        warning(paste("Null model fitted in formula ", i, ". (lp set to zero)", sep=""))
      }
     NULL
   }

   ## weighted average of predictions
   pred <- (lpmat %*% obj$m[,"weights"]) / sum(obj$m[,"weights"])
   ## calculate variance as equation (9) in Buckland et al 1997
   biasmat <- sweep(lpmat, 1, pred) # subtract predictor
   se.fit <- sqrt(varmat + biasmat^2) %*% obj$m[,"weights"]
   ci.vec <- matrix(NA,nrow=length(pred), ncol=2)
   ci.vec[,1] <- pred - qnorm(1-(alpha.c/2)) * se.fit
   ci.vec[,2] <- pred + qnorm(1-(alpha.c/2)) * se.fit

   ## calculation of prediction intervals for linear regression if wanted
   if (fam=="gaussian"){
     ## 1-alpha.p prediction interval
     pci.vec<-matrix(NA,nrow=length(pred), ncol=2)
     pci.vec[,1] <- pred - qnorm(1-(alpha.p/2)) * sqrt(se.fit^2 + se.err^2)
     pci.vec[,2] <- pred + qnorm(1-(alpha.p/2)) * sqrt(se.fit^2 + se.err^2)
     return(list(pred=pred, se.fit=se.fit, ci=ci.vec, pintv=pci.vec))
   } else return(list(pred=pred, se.fit=se.fit, ci=ci.vec))
 } # end else if(best)
 
} ## end of function fma.pred.glm

