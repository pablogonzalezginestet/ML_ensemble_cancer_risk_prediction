## Table 1

sarc$educ_level <- addNA(sarc$educ_level)
x <- table1(~ age_d1 + age_d2 + factor(male) + cut(I(as.numeric(format(start_date, "%Y"))), c(2005, 2009, 2014, 2017, 2019), right = TRUE) +
         I(county == "01") + factor(nordic) + educ_level + civil_status + inpt_visits + outp_visits +
         `NUM_UNIQUE(drugs.ATC)`, data = sarc)

t1kable(x, booktabs = FALSE, format = "latex")

## list of all features

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

alldat <- data.frame(XX)

x2 <- table1( ~ ., data = alldat)
cat(t1kable(x2, booktabs = FALSE, format = "latex"), file = "longtable.txt")
