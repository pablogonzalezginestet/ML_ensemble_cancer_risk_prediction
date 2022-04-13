library(survival)

sarc <- as.data.frame(readRDS("training-file.rds"))

sarc$end_date[sarc$outcome == 0] <- as.Date("2019-12-31")
sarc <- subset(sarc, start_date >= as.Date("2006-01-01") & end_date <= as.Date("2019-12-31"))
sarc$time <- with(sarc, as.numeric(end_date - start_date)/(365.25))
sarc <- subset(sarc, time > 0 & age_d1 >= 18)

sfit <- survfit(Surv(time, outcome) ~ 1, data = sarc)
summary(sfit, times = c(5, 10))

pdf("cumincfig.pdf", width=5.5, height = 4.25)
plot(sfit, noplot = c("(s0)", "2", "3"), conf.int = TRUE, xlab = "Years since diagnosis", ylab = "Cumulative incidence of cancer")
dev.off()

YY = with(sarc, Surv(time, outcome))

sarc$nordic[sarc$nordic == 9] <- 0
sarc$`MONTH(end_date)` <- NULL
sarc$`TIME_SINCE(end_date)` <- NULL
sarc$`WEEKDAY(end_date)` <- NULL
sarc$inpt_visits <- NULL
sarc$outp_visits <- NULL
sarc$`TIME_SINCE(birth_d)` <- NULL

sarc[, grep("ATC(4|5)", colnames(sarc))] <- NULL

for(j in which(sapply(sarc, is.factor))[-4]) {

  t1 <- table(sarc[, j])
  dlevs <- which(t1 < 300)
  if(length(dlevs) > 0) {
    tmp <- as.character(sarc[[j]])

    for(d in names(dlevs)) {
      tmp[sarc[[j]] == d] <- "other"
    }
  sarc[[j]] <- factor(tmp, exclude = NULL)
  }
}

sarc$`FIRST(drugs.forpddd)`[is.na(sarc$`FIRST(drugs.forpddd)`)]<- 0


XX <- model.matrix(~ ., data = sarc[, -which(colnames(sarc) %in% c("lopnr", "cancer.type", "start_date", "end_date",
                                         "outcome", "time"))])
colnames(XX) <- make.names(colnames(XX))
XX <- XX[, -1]

alldat <- data.frame(YY, XX)


source("00-model-specs.R")

mymods <- list(stepwise, coxglmnet,
               function(xx) random.forest(xx, time = 5),
               function(xx) random.forest(xx, time = 10),
               coxboost,
               function(xx) direct.binomial(xx, time = 5),
               function(xx) direct.binomial(xx, time = 10),
               function(xx) svm(xx, time = 5),
               function(xx) svm(xx, time = 10),
               function(xx) knn(xx, time = 5),
               function(xx) knn(xx, time = 10),
               function(xx) pseudo.glm.age(xx, time = 4  ),
               function(xx) pseudo.glm.age(xx, time = 5  ),
               function(xx) pseudo.glm.age(xx, time = 9  ),
               function(xx) pseudo.glm.age(xx, time = 10  ),
               function(xx) pseudo.glmnet(xx, time = 4  ),
               function(xx) pseudo.glmnet(xx, time = 5  ),
               function(xx) pseudo.glmnet(xx, time = 9  ),
               function(xx) pseudo.glmnet(xx, time = 10  ))


set.seed(220317)
validdex <- sample(1:nrow(alldat), 3137)
validset <- alldat[validdex, ]

devel <- alldat[-validdex, ]

set.seed(420)
ndex <- 1:nrow(devel)
part <- as.factor(sample(1:10, length(ndex), replace = TRUE))
folds <- split(ndex, part)

Zout <- vector(mode = "list", length = 10)
for(j in 1:length(folds)) {

  training <- devel[unlist(folds[-j]), ]
  validation <- devel[folds[[j]], ]

  Zout[[j]] <- do.call(cbind, lapply(mymods, function(f){

    fhat <- f(training)
    fhat(validation)

  }))

  cat(j, "\n")

}


fullfits <- lapply(mymods, function(f) {

  f(devel)

})


saveRDS(fullfits, "fullfits.rds")
saveRDS(Zout, "Zout.rds")
saveRDS(folds, "split.rds")
saveRDS(validset, "holdout.rds")
