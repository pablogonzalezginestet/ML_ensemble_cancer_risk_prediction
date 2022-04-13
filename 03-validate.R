library(pseudoloss)  ##
#source("00-model-specs.R")

library(xtable)
library(ggplot2)
library(patchwork)
library(broom)

fullfits <- readRDS("fullfits.rds")
cv.auc <- readRDS("cv-aucs-stack.rds")
opt.fit <- readRDS("opt-coeffs.rds")

valid <- validset
Zvalid <- lapply(fullfits, function(f) f(valid))

valid$time <- valid$YY[, 1]
valid$status <- factor(valid$YY[, 2])
valid$YY <- NULL
valid$PO <- cumincglm(Surv(time, status) ~ 1, data = valid, time = 5, cause = 1)$y
valid$PO10 <- cumincglm(Surv(time, status) ~ 1, data = valid, time = 10, cause = 1)$y

YYensvalid <- valid$PO
YYensvalid10 <- valid$PO10

OPO5 <- cbind(cumincglm(Surv(time, status) ~ 1, data = valid, time = 5, cause = 2)$y,
      cumincglm(Surv(time, status) ~ 1, data = valid, time = 5, cause = 3)$y)
OPO10 <- cbind(cumincglm(Surv(time, status) ~ 1, data = valid, time = 10, cause = 2)$y,
              cumincglm(Surv(time, status) ~ 1, data = valid, time = 10, cause = 3)$y)

aucalls <- lapply(Zvalid, function(z) {


  with(calc_roc(z, as.matrix(cbind(YYensvalid, OPO5))),
       pseudoloss::calc_auc(fpf, tpf))

})

aucalls10 <- lapply(Zvalid, function(z) {


  with(calc_roc(z, as.matrix(cbind(YYensvalid10, OPO10))),
       pseudoloss::calc_auc(fpf, tpf))

})

names.models <- c("Cox.stepwise", "Cox.lasso", "survival random forests time 5",
                 "survival random forests time 10",  "CoxBoost",
                  "IPCW logistic lasso time 5", "IPCW logistic lasso time 10",
                 "Bagged SVM time 5", "Bagged SVM time 10",
                 "Bagged KNN time 5", "Bagged KNN time 10",
                  "eventglm age, sex 4yr", "eventglm age, sex 5yr",
                 "eventglm age, sex 9yr", "eventglm age, sex, 10yr",
                  "LASSO eventglm 4yr", "LASSO eventglm 5yr",
                  "LASSO eventglm 9yr", "LASSO eventglm 10yr")


Zmat.valid <- do.call(cbind, Zvalid)

Zfin.opt <- (Zmat.valid %*% (opt.beta)/ sum(opt.beta))[, 1]
roc.opt <- calc_roc(Zfin.opt, as.matrix(cbind(YYensvalid, OPO5)))

Zfin.opt10 <- (Zmat.valid %*% (opt.beta10) / sum(opt.beta10))[, 1]
roc.opt10 <- calc_roc(Zfin.opt10, as.matrix(cbind(YYensvalid10, OPO10)))


with(roc.opt, pseudoloss::calc_auc(fpf,tpf))
with(roc.opt10, pseudoloss::calc_auc(fpf,tpf))

sumtable <- data.frame(R.package = c("survival", "survival",
                                     "randomForestSRC", "randomForestSRC",
                                     "CoxBoost", "glmnet", "glmnet",
                                     "e1071", "e1071", "class", "class", rep("eventglm", 4),
                                     rep("glmnet", 4)),
                       AUC = unlist(aucalls),
                       coefficient = (opt.beta) / sum(opt.beta),
                       AUC10 = unlist(aucalls10),
                       coefficient.10 = opt.beta10 / sum(opt.beta10))
rownames(sumtable) <- names.models


cat(print(xtable(sumtable, digits = 3)), file = "ctable.txt")




## roc and predictiveness curve

roc.opt$cut <- sort(unique(Zfin.opt))
roc.opt$cut <- NULL
roc.all <- rbind(cbind(time = 5, roc.opt), cbind(time = 10, roc.opt10))

p1<- ggplot(roc.all, aes(x = fpf, y = tpf, color = factor(time))) + geom_line() +
  theme_bw()+ plotROC::style_roc() + geom_abline(intercept = 0, slope = 1, color = "grey90") +
  theme(legend.position = "none")

Zopt.cut <- pmin(pmax(0, Zfin.opt), .99)


p2 <- ggplot(data.frame(pred.risk = c(Zopt.cut, Zopt.cut),
                        pos = c(YYensvalid, YYensvalid10),
                        time = rep(c(5, 10), each  = length(Zopt.cut))),
             aes(x = pred.risk, y = pos, color = factor(time))) +
  stat_smooth(method = "loess", span = 20, se = TRUE) +
  scale_x_continuous("Predicted risk") +
  scale_y_continuous("Estimated risk") +
  geom_rug(aes(x= pred.risk, y = NULL, color = NULL), sides = "b") +
  geom_hline(yintercept = mean(YYensvalid), linetype = 3, color = "salmon") +
  geom_hline(yintercept = mean(YYensvalid10), linetype = 3, color = "blue") +
  coord_cartesian(ylim = c(0, .35)) +
  theme_bw()
p1 + p2
ggsave("int-plot.pdf", width = 7.25, height = 3.75)

valid$Zopt.cut <- Zopt.cut
#concordance(Surv(time, status) ~ Zopt.cut, data = valid)


valid$group.Z <- cut(Zfin.opt,
                     breaks = c(-1, quantile(Zfin.opt, c(0.25, .5, .75)), 2),
                     labels = c("Q1", "Q2", "Q3", "Q4"))
valid$group.Z10 <- cut(Zfin.opt10,
                       breaks = c(-1, quantile(Zfin.opt10, c(0.25, .5, .75)), 2),
                       labels = c("Q1", "Q2", "Q3", "Q4"))

sfitg <- survfit(Surv(time, status) ~ group.Z, data = valid)
sfitg10 <- survfit(Surv(time, status) ~ group.Z10, data = valid)

meplot <- tidy(sfitg)
meplot2 <- tidy(sfitg10)
meplot2 <- subset(meplot2, state == "1")
meplot <- subset(meplot, state == "1")
meplot$strata <- gsub("group.Z=", "", meplot$strata, fixed = TRUE)
meplot2$strata <- gsub("group.Z10=", "", meplot2$strata, fixed = TRUE)

qplot <- rbind(cbind(meplot, stack.time = 5), cbind(meplot2, stack.time = 10))

ggplot(qplot, aes(x = time, y = estimate, color = strata)) +
  geom_step() +
  geom_step(aes(y = conf.low), linetype = 3) +
  geom_step(aes(y = conf.high), linetype = 3) +
  theme_classic() +
  scale_color_grey("Risk quartile", start = .7, end = .1) +
  xlab("Years since diagnosis") + ylab("Cumulative incidence of cancer") +
  facet_wrap(~ stack.time, labeller = "label_both") +
  theme(legend.position = "bottom")
ggsave("km-fig.pdf", width = 6.25, height = 4.75)

## distribution of top 5 features by risk quantile

fvars <- c("age_d2", 
           "age_d1", 
           "X.COUNT.drugs..",
           "male",
           "civil_statusOG")

library(table1)

x432 <- table1(  ~age_d2 + age_d1 + X.COUNT.drugs..+ factor(male) + factor(civil_statusOG) | quanti, 
         data = data.frame(quanti = valid$group.Z10, valid[, fvars]))


t1kable(x432, booktabs = FALSE, format = "latex")

## comparing bottom 10% to top 10%
hilo <- valid
hilo$group5 <- cut(Zfin.opt,
    breaks = c(-1, quantile(Zfin.opt, c(0.1, .5, .9)), 2),
    labels = c("bottom 10", "m1", "m2", "top 10"))
hilo$group10 <- cut(Zfin.opt10,
                   breaks = c(-1, quantile(Zfin.opt10, c(0.1, .5, .9)), 2),
                   labels = c("bottom 10", "m1", "m2", "top 10"))

shil5 <- summary(survfit(Surv(time, status) ~ group5,
        data = hilo),
        times = 5)

sprintf("%.3f (%.3f to %.3f)",
        shil5$pstate[c(1, 4), 2], shil5$lower[c(1,4), 2], shil5$upper[c(1,4), 2])

shil10 <- summary(survfit(Surv(time, status) ~ group10,
                         data = hilo),
                 times = 10)

sprintf("%.3f (%.3f to %.3f)",
        shil10$pstate[c(1, 4), 2], shil10$lower[c(1,4), 2], shil10$upper[c(1,4), 2])
