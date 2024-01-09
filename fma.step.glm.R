fma.step.glm <- function(dset, newdata=NA, fam="binomial",
                         predtype="response", index,
                         varnames, full.formula, lower.formula,
                         crit.penal=2, direct="backward",
                         screen=T, screen.lev=0.3, screen.res=NA,
                         alpha.c=0.05, alpha.p=0.05, nval=1) {
                  
## ----------------------------------------------------------------------------
## Title: fma.step.glm
## ----------------------------------------------------------------------------
## Authors: Willi Sauerbrei, Nicole Augustin, Norbert Hollaender, Anika Buchholz
## Center for Medical Biometry and Medical Informatics
## Stefan-Meier-Strasse 26, D-79104 Freiburg 
## http://www.imbi.uni-freiburg.de 
## ----------------------------------------------------------------------------
## Description: Step 2 of the model averaging process: model selection
##              for generalised linear models (e.g. linear and
##              logistic regression)
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: fma.step.glm(dset, newdata=NA, fam="binomial",
##                     predtype="response", index,
##                     varnames, full.formula, lower.formula,
##                     crit.penal=2, direct="backward",
##                     screen=T, screen.lev=0.3, screen.res=NA, nval=1)
##  
##        dset:          build data set
##        newdata:       new data to be predicted (default NA)
##        fam:           character string specifying the regression
##        predtype:      type of output desired, see 'predict.cph'
##                       (default NA)   
##        index:         index matrix for bootstrap samples / subsamples
##        varnames:      character vector containing the names of the
##                       covariates                   
##        full.formula:  formula specifying the full model
##        lower.formula: specifies lower formula for stepwise
##                       variable selection (default '"1"')
##                       model that shall be used, e.g. '"binomial"'
##                       for logistic regression (default) or
##                       '"gaussian"' for linear regression
##        crit.penal:    penalty criterion for stepwise variable
##                       selection, e.g. 'crit.penal = 2' for AIC
##                       (default) or 'crit.penal = log(n)' for BIC                   
##        direct:        direction of stepwise variable selection,
##                       being either '"backward"' (default),
##                       '"forward"' or '"both"'
##        screen:        'TRUE' if variable screening step has been
##                       performed beforehand (default), 'FALSE'
##                       otherwise
##        screen.lev:    selection level for variable screening
##                       (default 0.3)
##        screen.res:    object obtained from variable screening by
##                       function 'select.boot.cox', compulsory if
##                       screen=T
##        alpha.c:       alpha level for confidence interval (to be
##                       used with generated (simulated) data) 
##        alpha.p:       alpha level for prediction interval (to be
##                       used with real data)
##        nval:          number of new data sets (default is 1,
##                       for simulations nval may be larger than 1)
## ----------------------------------------------------------------------------
## Value: fma.res: list of results, as follows 
##         $m:         selected models with index matrix of included
##                     covariates     
##         $screen:    bootstrap / subsampling selection frequencies of
##                     covariates (obtained from STEP 1: screening step)
##         $weights:   bootstrap / subsampling model weights 
##         $formu:     formulae of selected models 
##         $selfreq:   index matrix of selected covariates
##         $ci:       (1-alpha.c) bootstrap / subsampling confidence interval
##                     (to be used with generated,  i.e. simulated, data)
##         $pi:       (1-alpha.p) bootstrap / subsampling prediction interval
##                     (to be used with real data)  
##         $bagg:      bagging estimator for the prediction
##         $se.bagg:   estimated standard error of bagg
##         $screen.lev:selection level used for variable screening
## ----------------------------------------------------------------------------

                                      
  sep <-"+"
  nvar <- length(varnames)
  assign("fam", fam, 1)

  nboot <- nrow(index)
  
  ### step 1: screening of variables
  if (screen) {
    s <- screen.res$screen
    # extract response from 'full.formula'                 
    response <- rownames(attr(terms(full.formula), "factors"))[1]
    ### maximum model after screening
    full.formula <- paste(names(s)[s>=screen.lev],collapse=sep)
    full.formula[full.formula==""] <- 1
    full.formula <- as.formula(paste(response,"~",full.formula))
  } # end if screen
  else print("NO SCREENING!!!")

  ### step 2: estimating model weights
  print("estimating model weigths")
  fma.res <- screening.step.glm(dset, newdata, fam, predtype,
                                index, varnames, full.formula,
                                lower.formula, crit.penal, direct,
                                alpha.c, alpha.p, nval)

  if (screen) {
    fma.res$screen.lev <- screen.lev
    fma.res$screen <- screen.res$screen    
  }
  else {
    fma.res$screen.lev <- 0
    fma.res$screen <- NA
  }

  fma.res

} ## end of function fma.step.glm

