screening.step.glm  <- function(dset, newdata=NA, fam="binomial",
                                predtype="response", index,
                                varnames, full.formula, lower.formula,
                                crit.penal=2, direct="backward",
                                alpha.c=0.05, alpha.p=0.05, nval=1) {
  ## -----------------------------------------------------------------------------
  ## Title: R-function screening.step.glm()
  ## -----------------------------------------------------------------------------
  ## Authors: Willi Sauerbrei, Nicole Augustin, Norbert Hollaender, Anika Buchholz
  ## Center for Medical Biometry and Medical Informatics
  ## Stefan-Meier-Strasse 26, D-79104 Freiburg,
  ## http://www.imbi.uni-freiburg.de
  ## -----------------------------------------------------------------------------
  ## Description: -This functions can be used as stand alone function in the screening step of
  ##               the two-step bootstrap model averaging approach (interesting output is SCREEN)
  ##               and it is also used in the model averaging step (step 2)
  ##              -performs variable selection in a GENERALIZED LINEAR MODEL (GLM)
  ##               within a set of bootstrap samples and counts bootstrap
  ##               selection frequencies for each covariate
  ##              -calculates a bagging (bootstrap averaging) estimator for 
  ##               the prediction of the response in a new data set   
  ##              -calculates confidence and prediction intervals 
  ## -----------------------------------------------------------------------------
  ## Required Packages: boot 
  ## -----------------------------------------------------------------------------
  ## Usage: screening.step.glm(dset, newdata=NA, fam="binomial",
  ##                           predtype="response", index,
  ##                           varnames, full.formula, lower.formula,
  ##                           crit.penal=2, direct="backward",
  ##                           alpha.c=0.05, alpha.p=0.05, nval=1)
  ##
  ##   dset               input data set
  ##   newdata            new data set for prediction
  ##   fam                family fam="gaussian" for linear model, fam="binomial" for logistic model
  ##   predtype           type of prediction (default is "response")
  ##   index              subsampling index array 
  ##   varnames           names of covariates (e.g. varnames=colnames(dset)[-1])
  ##   full.formula       formula of the full model
  ##   lower.formula      formula of smallest model (e.g. paste(colnames(dset)[1], "~1"))
  ##   crit.penal         penalty term (AIC:crit.penal=2, BIC:crit.penal=log(nrow(dset)))
  ##   direct             direction of variable selection procedure (e.g. direct="backward")
  ##   alpha.c            alpha level for confidence interval (to be used with generated (simulated) data) 
  ##   alpha.p            alpha level for prediction interval (to be used with real data) 
  ##   nval               number of new data sets (default is 1, for simulations nval may be larger than 1)
  ## -----------------------------------------------------------------------------
  ## Value: list of output values
  ##   m         selected models with index matrix of included covariates (STEP 2: model averaging step)    
  ##   screen    sampling selection frequencies of covariates (STEP 1: screening step)
  ##   sampw     sampling model weights (STEP 2)
  ##   formu     formulae of selected models (STEP 2)
  ##   selfreq   index matrix of selected covariates, used to calculate `screen' (STEP 1)
  ##   ci        (1-alpha.c) bootstrap/subsampling confidence interval (to be used with generated (simulated) data)
  ##   pi        (1-alpha.p) bootstrap/subsampling prediction interval (to be used with real data)  
  ##   bagg      bagging estimator for the prediction
  ##   se.bagg   estimated standard error of bagg
  ## -----------------------------------------------------------------------------

  
  if (nval==1 & !all(is.na(newdata))){
    dat <- newdata
    newdata <- list(NULL)
    newdata[[1]] <- dat
  }
  assign("newdata",newdata,1)
  assign("predtype",predtype,1)

  nboot <- nrow(index)
  sampsize <- ncol(index)
  
  sep <- "+"
  nvar <- length(varnames)
  n.col <- ifelse(attributes(terms(full.formula))$term.labels==".",
                  length(varnames),
                  length(attributes(terms(full.formula))$term.labels))
  if (attributes(terms(full.formula))$term.labels[1]=="."){
    col.nam <- varnames
  } else col.nam <- attributes(terms(full.formula))$term.labels
  
  selfreq <- matrix(0, nrow=nboot, 
  	            ##ncol=length(varnames), colnames=varnames)
		    ncol=n.col,
                    dimnames=list(NULL, col.nam))
  mods <- vector(mode="character", length=nboot)
  
  ## matrix to store fitted/predicted values for bootstrap confidence intervals
  if (is.null(dim(newdata[[1]]))){
    lpmat  <- matrix(nrow=sampsize, ncol=nboot)
    varmat <- matrix(nrow=sampsize, ncol=nboot) 
    ## bootstrap prediction interval  
    pimat <- matrix(nrow=sampsize, ncol=(2*nboot))
  } else {
    if (nval!=0){
      lpmat <- list(NULL)
      pimat <- list(NULL)
      varmat <- list(NULL)
      for (j in 1:nval){
        lpmat[[j]] <- matrix(0, nrow=nrow(newdata[[j]]), ncol= nboot)
        varmat[[j]] <- matrix(0, nrow=nrow(newdata[[j]]), ncol= nboot)        
        ## bootstrap prediction interval         
        pimat[[j]] <- matrix(0, nrow=nrow(newdata[[j]]), ncol= (2*nboot))
      }}
  }
  
  for (i in 1:nboot) {
    if (i%%50==0) print(paste(paste(i,"of"),nboot))
    ii <- index[i,]
    assign("dset.ii", dset[ii,], 1)
    if (direct=="backward") {
      b.full <- glm(full.formula, family=fam, data=dset.ii)
    } else b.full <- glm(lower.formula, family=fam, data=dset.ii)
    ## if full.formula is null model no need for model selection
    if (!is.null(dim(attr(terms(full.formula),"factors")))){
      b.step <- step(b.full, scope = list(upper = full.formula, lower = lower.formula),
                     k=crit.penal, trace=0, direction=direct)
    } else b.step <- b.full
    selfreq[i,attributes(b.step$terms)$term.labels] <- 1
    mods[i] <- paste(attributes(b.step$terms)$term.labels,collapse=sep)
    ## fitted values/predictions for bootstrap confidence intervals
    if (!is.null(dim(newdata[[1]])) & nval!=0){
      ## j is number of validation data sets
      for (j in 1:nval){
        pred <- predict(b.step, newdata=newdata[[j]],type=predtype,se.fit=T)
        lpmat[[j]][,i] <- pred$fit
        varmat[[j]][,i] <- pred$se.fit^2        
        pimat[[j]][,i]            <- pred$fit - qnorm(1-(alpha.p/2)) * sqrt(pred$se.fit^2 + pred$residual.scale^2)
        pimat[[j]][,(i + nboot)]  <- pred$fit + qnorm(1-(alpha.p/2)) * sqrt(pred$se.fit^2 + pred$residual.scale^2)     
      }
    } else if (nval!=0){
      pred <- predict(b.step,type=predtype, se.fit=T)
      lpmat[,i] <- pred$fit
      varmat[,i] <- pred$se.fit^2           
      pimat[,i]            <- pred$fit - qnorm(1-(alpha.p/2)) * sqrt(pred$se.fit^2 + pred$residual.scale^2)
      pimat[,(i + nboot)]  <- pred$fit + qnorm(1-(alpha.p/2)) * sqrt(pred$se.fit^2 + pred$residual.scale^2)   
    }
    NULL
  } # end for loop
  
  if (nval!=0){
    if (!is.null(dim(newdata[[1]]))){
      ci <- list(NULL)
      pci <- list(NULL)
      biasmat <- list(NULL)
      bagg <- list(NULL)
      varbagg <- list(NULL)
      se.bagg <- list(NULL)
      se.bagg.u <- list(NULL)
      ## bootstrap confidence intervals
      for (j in 1:nval){
        ci[[j]]  <- t(apply(lpmat[[j]],1,quantile,c(alpha.c/2,(1-alpha.c/2))))
        pci[[j]] <- t(apply(pimat[[j]],1,quantile,c(alpha.p/2,(1-alpha.p/2))))      
        ##  bagging estimate 
        ##  (this is the mean of the fitted values from all bootrap samples)
        bagg[[j]] <- apply(lpmat[[j]],1,mean)
        ## calculate variance estimation similar to that of Buckland et al. formula (9)
        biasmat[[j]] <- sweep(lpmat[[j]], 1, bagg[[j]]) # subtract predictor
        ## create vector with artifical weights 1/n.boot for each bootstrap sample
        weights <- rep((1/nboot), nboot) 
        se.bagg[[j]] <- sqrt(varmat[[j]] + biasmat[[j]]^2) %*% weights
        ## assuming independence between bootstrap results Buckland et al. formula (10)
        se.bagg.u[[j]] <- sqrt( (varmat[[j]] + biasmat[[j]]^2) %*% weights^2 )      
        ## leads to extrem underestimation of S.E.(bagg) 
      }
    } else{
      ci <- NA
      pci <- NA
      biasmat <- NA
      bagg <- NA
      varbagg <- NA
      se.bagg <- NA
      se.bagg.u <- NA
      ## bootstrap confidence intervals
      ci <- t(apply(lpmat,1,quantile,c(alpha.c/2,(1-alpha.c/2))))
      pci <- t(apply(pimat,1,quantile,c(alpha.p/2,(1-alpha.p/2))))      
      ##  bagging estimate 
      ##  (this is the mean of the fitted values from all bootrap samples)
      bagg <- apply(lpmat,1,mean)
      ##  calculate variance estimation similar to that of Buckland et al. formula (9)
      biasmat <- sweep(lpmat, 1, bagg) # subtract predictor
      ##  create vector with artifical weights 1/n.boot for each bootstrap sample
      weights <- rep((1/nboot), nboot) 
      se.bagg <- sqrt(varmat + biasmat^2) %*% weights
      ##  assuming independence between bootstrap results Buckland et al. formula (10)
      ##     se.bagg.u <- sqrt( (varmat + biasmat^2) %*% weights^2 )      
      ##  leads to extrem underestimation of S.E.(bagg) 
    } # end of else(!is.null(dim(...)))
  } # end of if(nval!=0)
  ## variabe inclusion frequencies
  screen <- apply(selfreq,2,sum)/nboot
  ## calculate bootstrap weights
  bootw <- table(mods)/sum(table(mods))
  formu <- names(bootw)
  formu[formu==""] <- 1
  y.var <- colnames(dset)[1]
  formu <- paste(paste(y.var,"~"),formu)
  ind <- rev(order(bootw))
  bootw <- bootw[ind]
  formu <- formu[ind]
  ## create matrix with result
  m <- matrix(0, nrow=length(formu),ncol=nvar)
  colnames(m) <- varnames
  for (i in 1:length(formu)) {
    whichvar <- all.vars(as.formula(formu[i]))[-1]
    m[i,whichvar] <- 1
    NULL
  }
  m <- cbind(bootw,m)
  colnames(m) <- c("weights", varnames)
  if (nval!=0){
    if (nval>1){
      if (fam=="gaussian")
        result <- list(m=m, screen=screen, weights=bootw, formu=formu, selfreq=selfreq, ci=ci, pi=pci, bagg=bagg,
                       se.bagg=se.bagg)
      else result <- list(m=m, screen=screen, weights=bootw, formu=formu, selfreq=selfreq, ci=ci, bagg=bagg,
                          se.bagg=se.bagg)
    }else{
      bci <- ci[[1]]
      bpi <- pci[[1]]
      
      if (fam=="gaussian")
        result <- list(m=m, screen=screen, weights=bootw, formu=formu,
                       selfreq=selfreq, ci=bci,
                       pi=bpi, bagg=unlist(bagg),
                       se.bagg=unlist(se.bagg))
       else result <- list(m=m, screen=screen, weights=bootw,
                           formu=formu, selfreq=selfreq,
                           ci=bci, bagg=unlist(bagg),
                           se.bagg=unlist(se.bagg))
    }
  }else  result <- list(m=m, screen=screen, weights=bootw, formu=formu, selfreq=selfreq)  
  class(result) <- c("bw")
  result
 
} # end of select.boot

