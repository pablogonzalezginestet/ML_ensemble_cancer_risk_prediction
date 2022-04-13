library(pseudoloss) ## remotes::install_github("sachsmc/pseudoloss")
#source("00-model-specs.R")

Zout <- readRDS("Zout.rds")
fullfits <- readRDS("fullfits.rds")
folds <- readRDS("split.rds")
devel <- readRDS("devel.rds")

devel$time <- devel$YY[, 1]
devel$status <- factor(devel$YY[, 2])
devel$YY <- NULL
devel$PO <- cumincglm(Surv(time, status) ~ 1, data = devel, time = 5, cause = 1)$y
devel$PO10 <- cumincglm(Surv(time, status) ~ 1, data = devel, time = 10, cause = 1)$y

YYens <- devel$PO[unlist(folds)]
YYens10 <- devel$PO10[unlist(folds)]
Zmat <- do.call(rbind, Zout)


opt.fit2 <- glm(yy ~ . -1, data = data.frame(yy = YYens, Zmat), family = "gaussian")
opt.beta <- opt.fit2$coefficients

negloglik <- function(b) {
  b[-1] <- exp(b[-1])
  yhat <- cbind(1, Zmat) %*% b
  sum((YYens - yhat)^2) + 0 * abs(sum(b))

}

negloglik10 <- function(b) {
  b[-1] <- exp(b[-1])
  yhat <- cbind(1, Zmat) %*% b
  sum((YYens10 - yhat)^2) + 0 * abs(sum(b))

}

#nlm(negloglik, rep(0, 19), iterlim = 4000)
opt.beta <- exp(optim(rep(0, 20), fn = negloglik, method = "CG",
                      control = list(maxit = 500))$par)
opt.beta <- opt.beta[-1]

opt.beta10 <- exp(optim(rep(0, 20), fn = negloglik10, method = "CG",
                      control = list(maxit = 500))$par)

opt.beta10 <- opt.beta10[-1]
