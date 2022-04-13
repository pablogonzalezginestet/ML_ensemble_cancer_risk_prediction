library(eventglm) ## remotes::install_github("sachsmc/eventglm")
library(survival)
library(randomForestSRC)
library(CoxBoost)  ## remotes::install_github("binderh/CoxBoost")
library(riskRegression)
library(e1071)
library(splines)
library(glmnet)
library(class)

## models
bigform <- as.formula(paste("~", paste(colnames(XX), collapse = "+")))
## survival
############### Natively ############################
#' @description  Predictions based on those Machine Learning procedures in the library that allow for weights to be specified as an argument of the R function. No bagging occurs. This group of algorithms is denoted as  Native Weights
#' @param dset data set
#' @return  prediction
#' @rdname megalearner-internal
#'
#'
stepwise <- function(dset) {

  dset$id <- 1:nrow(dset)
  dset$yy <- Surv(dset$YY[, 1], ifelse(dset$YY[, 2] == 1, 1, 0))
  
  cfit <- coxph(yy ~ age_d1 + age_d2 + male, data = dset, id = id)
  cfin <- step(cfit, scope = bigform, trace = 0)
  function(valid) {
    predict(cfin, newdata = valid,
            type = "lp")
  }

}

coxglmnet <- function(dset) {
  
  yy <- Surv(dset$YY[, 1], ifelse(dset$YY[, 2] == 1, 1, 0))
  
  xx <- as.matrix(dset[, -1])
  cnfit <- glmnet(xx, yy, family = "cox",
                  alpha = .5, lambda = exp(-4.5))
  
  function(valid) {
    vx <- as.matrix(valid[, -1])
    predict(cnfit, newx = vx)[, 1]
    
  }
  
}

random.forest <- function(dset, time = 5) {

  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL

  cfit <- rfsrc(Surv(time, status) ~ ., data = dset)
  function(valid) {
    res <- predict(cfit, newdata = valid[, -1])
    res$cif[, max(which(res$time.interest <= time))[1], 2]
  }

}

coxboost <- function(dset) {

  mm <- model.matrix(YY ~ ., data = dset)[, -1]
  cfit <- CoxBoost(time = dset$YY[, "time"], status = dset$YY[, "status"],
           x = mm)
  function(valid) {
    mmv <- model.matrix(YY ~ ., data = valid)[, -1]
    res <- predict(cfit, newdata = mmv,
                   type = "lp")
    res[1, ]
  }

}

## binary

direct.binomial <- function(dset, time = 5) {
  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL
  dset$ybin <- 1.0 * (dset$time < time & dset$status == 1)
  dset$ybin[dset$time < time & dset$status == 0] <- NA

  sfit <- survfit(Surv(time, ifelse(status == 0, 1, 0)) ~ 1, data = dset)
  dset$weights <- 1 / summary(sfit, times = pmin(dset$time, time))$surv
  dset$time <- dset$status <- NULL

  mdex <- which(!is.na(dset$ybin))
  cfit <- glmnet(as.matrix(dset[mdex, -which(colnames(dset) %in% c("ybin", "weights"))]), dset$ybin[mdex], 
                    family = "binomial", weights = dset$weights[mdex], 
                  alpha = .25, lambda = exp(-6))
  
  function(valid) {
    vx <- as.matrix(valid[, -1])
    predict(cfit, newx = vx, type = "response")[, 1]
    
  }

}

## bagging ipcw

svm <- function(dset, time = 5){

  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL
  dset$ybin <- 1.0 * (dset$time < time & dset$status == 1)
  dset$ybin[dset$time < time & dset$status == 0] <- NA

  sfit <- survfit(Surv(time, ifelse(status == 0, 1, 0)) ~ 1, data = dset)

  dset <- dset[!is.na(dset$ybin),]
  wtmp <- 1 / summary(sfit, times = pmin(dset$time, time))$surv
  dset$samp.wts <- wtmp / sum(wtmp)
  dset$time <- dset$status <- NULL

  svmboot <- lapply(1:10, function(i) {

    dboot <- dset[sample(1:nrow(dset), nrow(dset),
                         replace = TRUE, prob = dset$samp.wts),]
    dboot$samp.wts <- NULL
    e1071::svm(update.formula(bigform, ybin ~ .), data = dboot)

  })


  function(valid) {

    vboots <- sapply(svmboot, function(cfit) {
      predict(cfit, newdata = valid)
    })
    rowMeans(vboots)

  }
}



knn <- function(dset, time = 5){

  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL
  dset$ybin <- 1.0 * (dset$time < time & dset$status == 1)
  dset$ybin[dset$time < time & dset$status == 0] <- NA

  sfit <- survfit(Surv(time, ifelse(status == 0, 1, 0)) ~ 1, data = dset)

  dset <- dset[!is.na(dset$ybin),]
  wtmp <- 1 / summary(sfit, times = pmin(dset$time, time))$surv
  dset$samp.wts <- wtmp / sum(wtmp)
  dset$time <- dset$status <- NULL


  function(valid) {

    vboots <- sapply(1:10, function(i) {

      dboot <- dset[sample(1:nrow(dset), nrow(dset),
                           replace = TRUE, prob = dset$samp.wts),]
      dboot$samp.wts <- NULL
      tclass <- as.factor(dboot$ybin)
      dboot$ybin <- NULL

      knboot <- class::knn(dboot, valid[, colnames(dboot)],
          cl = tclass, k = 5, prob = TRUE)

      ifelse(knboot == 0, 1 - attr(knboot, "prob"), attr(knboot, "prob"))

    })

    rowMeans(vboots)

  }
}


## pseudo obs

pseudo.glm <- function(dset, time = 5 ){

  dset$time <- dset$YY[, 1]
  dset$status <- factor(dset$YY[, 2])
  dset$YY <- NULL
  cfit <- cumincglm(update.formula(bigform, Surv(time, status) ~ .), time = time, data = dset, cause = "1")
  function(valid) {
    predict(cfit, newdata = valid, type = "response")
  }

}


pseudo.glm.age <- function(dset, time = 5 ){
  
  dset$time <- dset$YY[, 1]
  dset$status <- factor(dset$YY[, 2])
  dset$YY <- NULL
  cfit <- cumincglm(Surv(time, status) ~ bs(age_d1, knots = c(20, 40, 60, 80), Boundary.knots = c(-1, 100)) * male, 
                    time = time, data = dset, cause = "1")
  function(valid) {
    predict(cfit, newdata = valid, type = "response")
  }
  
}

pseudo.glmnet <- function(dset, time = 5 ) {

  dset$time <- dset$YY[, 1]
  dset$status <- factor(dset$YY[, 2])
  dset$YY <- NULL
  cfit <- cumincglm(update.formula(bigform, Surv(time, status) ~ .), time = time, data = dset, cause = "1")
  
  pfac <- rep(1, ncol(cfit$x)-1)
  pfac[c(1, 2, 3, 19, 20)] <- 0
  #cnfit <- cv.glmnet(cfit$x[, -1], cfit$y, penalty.factor = pfac, alpha = .25)
  cnfit <- glmnet(cfit$x[, -1], cfit$y, penalty.factor = pfac, alpha = .25, lambda = exp(-6))
  function(valid) {
    vx <- as.matrix(valid[, -1])
    predict(cnfit, newx = vx)[, 1]

  }

}


