library(haven)
library(data.table)
library(reticulate)

cases <- subset(as.data.frame(read_sas("../inc_sarc_2visits_compar_0320.sas7bdat")),
                sarcoidosis == 1)
cases.lopnr <- cases$lopnr
outcome <- subset(as.data.frame(read_sas("../r_can__index11229_2020.sas7bdat")), LopNr %in% cases.lopnr)

outcome$cancer.date <- as.Date(outcome$DIADAT, format = "%Y%m%d")
outcome <- data.table(outcome, key = "LopNr")
outcome <- outcome[ben != "3"]
outcome[, first.diag := cancer.date == min(cancer.date), by = "LopNr"]
first.cancer <- unique(outcome[first.diag == TRUE, .(LopNr, cancer.date, ICD7)])
first.cancer[, cancer.type := paste(ICD7, collapse = ":"), by = "LopNr"]
first.cancer[, ICD7 := NULL]
first.cancer <- unique(first.cancer)
colnames(first.cancer)[1] <- "lopnr"

cases2 <- merge(cases, first.cancer, by = "lopnr", all.x = TRUE)

cases2$start_date <- cases$sarc_d1
cases2$end_date <- pmin(as.Date("2020-12-31"),
                        cases2$death_d, cases2$date_latest_emig,
                        cases2$cancer.date, na.rm = TRUE)

cases2$outcome <- factor(with(cases2, ifelse(!is.na(cancer.date) & end_date == cancer.date, 1,
                                             ifelse(!is.na(death_d) & end_date == death_d, 2,
                                                    ifelse(!is.na(date_latest_emig) & end_date == date_latest_emig, 3, 0)))),
                         levels = 0:3)

cases2 <- data.table(cases2, key = "lopnr")

### predictors

npr.out <- data.table(read_sas("../../../2021 Sarcoidosis Linkage/1. Raw data/Socialstyrelsen/National Patient Register/t_r_par_ov_index11229_2020.sas7bdat"))
npr.out <- npr.out[LopNr %in% cases.lopnr]
npr.out <- merge(npr.out, cases2, by.x = "LopNr", by.y = "lopnr")
npr.out[, npr.date := as.Date(INDATUMA, format = "%Y%m%d")]
npr.out <- npr.out[npr.date < start_date]
npr.out <- npr.out[HDIA != ""]
npr.out <- npr.out[, .(lopnr = LopNr, npr.date, HDIA, DIAGNOS, OP)]


npr.in <- data.table(read_sas("../../../2021 Sarcoidosis Linkage/1. Raw data/Socialstyrelsen/National Patient Register/t_r_par_sv_index11229_2020.sas7bdat"))
npr.in <- npr.in[LopNr %in% cases.lopnr]
npr.in <- merge(npr.in, cases2, by.x = "LopNr", by.y = "lopnr")
npr.in[, npr.date := as.Date(UTDATUMA, format = "%Y%m%d")]
npr.in <- npr.in[npr.date < start_date]
npr.in <- npr.in[HDIA != ""]
npr.in <- npr.in[, .(lopnr = LopNr, npr.date, HDIA, DIAGNOS, OP)]

drugs <- data.table(read_sas("../../../2021 Sarcoidosis Linkage/1. Raw data/Socialstyrelsen/Prescribed Drug Register/r_lmed__index11229_2020.sas7bdat"))
drugs <- drugs[LopNr %in% cases.lopnr]
drugs <- merge(drugs, cases2, by.x = "LopNr", by.y = "lopnr")
drugs <- drugs[EDATUM < start_date]
drugs <- drugs[, .(lopnr = LopNr, ATC, ATC5 = substr(ATC, 1, 5),
                   ATC4 = substr(ATC, 1, 4), ATC3 = substr(ATC, 1, 3), EDATUM, forpddd, lform)]
drugs[, drug.index := 1:.N]
npr.in[, npr.in.index := 1:.N]
npr.out[, npr.out.index := 1:.N]

cases3 <- data.table(cases2)
cases3 <- cases3[, .(lopnr, age_d1, age_d2, male, county, birth_d, nordic,
                     educ_level, civil_status, inpt_visits, outp_visits,
                     start_date, end_date, outcome)]
### featuretools

#virtualenv_create("feature", packages = "featuretools")
#virtualenv_install("feature", "featuretools")
use_virtualenv("feature")
ft <- import("featuretools")
pa <- import("pandas")
np <- import("numpy")

es <- ft$EntitySet()
es$add_dataframe(dataframe_name = "subjects", dataframe = cases3, index = "lopnr",
                 time_index = "start_date")
es$add_dataframe(dataframe_name = "npr_inpat", dataframe = npr.in, index = "npr.in.index",
                 time_index = "npr.date")
es$add_dataframe(dataframe_name = "npr_outpat", dataframe = npr.out, index = "npr.out.index",
                 time_index = "npr.date")
es$add_dataframe(dataframe_name = "drugs", dataframe = drugs, index = "drug.index",
                 time_index = "EDATUM")

es$add_relationship("subjects", "lopnr", "npr_inpat", "lopnr")
es$add_relationship("subjects", "lopnr", "npr_outpat", "lopnr")
es$add_relationship("subjects", "lopnr", "drugs", "lopnr")


deepfeatures <- ft$dfs(target_dataframe_name = "subjects",entityset = es,
                                        agg_primitives = list("mode", "trend", "percent_true", "mean",
                                                              "count", "num_true", "time_since_last",
                                                              "num_unique", "avg_time_between", "first", "last"),
                                        trans_primitives = list("time_since_previous", "time_since", "month", "weekday"),
                                        max_depth = 1)

featmat <- as.data.frame(deepfeatures[[1]])

conv <- which(sapply(featmat, function(x) class(x)[1]) == "pandas.core.arrays.integer.IntegerArray")
for(i in conv){

  featmat[[i]] <- unlist(np$asmatrix(featmat[[i]])[1,])

}

featmat$lopnr <- as.numeric(rownames(featmat))
outdat <- merge(cases2[, c("lopnr", "cancer.type", "start_date", "end_date")], featmat, by = "lopnr")

## fixing stuff

outdat2 <- subset(outdat, end_date > start_date)
outdat2$educ_level <- as.factor(outdat2$educ_level) |> addNA()
outdat2$`AVG_TIME_BETWEEN(npr_inpat.npr.date)`[is.na(outdat2$`AVG_TIME_BETWEEN(npr_inpat.npr.date)`)] <- max(outdat2$`AVG_TIME_BETWEEN(npr_inpat.npr.date)`, na.rm = TRUE)

for(j in grep("TIME_SINCE_LAST|AVG_TIME_BETWEEN", colnames(outdat2), value = TRUE)) {
  outdat2[[j]][is.na(outdat2[[j]])] <- max(outdat2[[j]], na.rm = TRUE)
}
for(j in grep("FIRST|LAST|MODE", colnames(outdat2), value = TRUE)){
  if(is.numeric(outdat2[[j]])) next
  outdat2[[j]] <- addNA(outdat2[[j]])
}
for(j in grep("NUM_UNIQUE|TREND|COUNT|MEAN", colnames(outdat2), value = TRUE)) {
  outdat2[[j]][is.na(outdat2[[j]])] <- 0
}

outdat2$`LAST(drugs.drug.index)` <- NULL
outdat2$`LAST(drugs.forpddd)`[is.na(outdat2$`LAST(drugs.forpddd)`)] <- 0

outdat2[, grep(".index)", colnames(outdat2), fixed = TRUE)] <- NULL

outdat2$`TIME_SINCE(start_date)` <- NULL
outdat2$`TIME_SINCE_PREVIOUS(start_date)` <- NULL

saveRDS(outdat2, "training-file.rds")
