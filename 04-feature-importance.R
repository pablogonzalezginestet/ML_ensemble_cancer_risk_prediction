library(iml)
library(randomForest)

pdat <- cbind(valid, pout = Zfin.opt)
pdat$PO10 <- pdat$PO <- pdat$Zopt.cut <-  NULL
pdat$group.Z <- pdat$group.Z10 <- pdat$time <- pdat$status <- NULL

surmod5 <- randomForest(pout ~ ., data = pdat, ntree = 50)

X <- pdat[, -which(names(pdat) == "pout")]
predictor <- Predictor$new(surmod5, data = X, y = pdat$pout)

imp <- FeatureImp$new(predictor, loss = "mae")
imp_results0 <- head(imp$results, 20)
imp_results <- data.frame(feature = imp_results0$feature,
                          `importance (95% CI)`=sprintf("%.3f (%.3f to %.3f)",
                                  imp_results0$importance,
                                  imp_results0$importance.05,
                                  imp_results0$importance.95),
                          error = imp_results0$permutation.error)

print(xtable::xtable(imp_results, digits = 2), include.rownames = FALSE)



set.seed(1288)
lime.explain <- LocalModel$new(predictor, 
                               x.interest = X[order(pdat$pout, decreasing = TRUE)[sample(1:40, 1)], ], 
                               k = 3)
lime.explain$results
p1 <- plot(lime.explain) + theme_bw()

lime.explain2 <- LocalModel$new(predictor, x.interest = X[order(pdat$pout, decreasing = FALSE)[sample(1:40, 1)], ], 
                                k = 3)
lime.explain2$results
p2 <- plot(lime.explain2) + theme_bw()

p1 + p2
pdf("hilo-risk-supp.pdf", width = 6.5, height = 5)
p1 + p2 + plot_layout(ncol = 1)
dev.off()


ale <- FeatureEffect$new(predictor, feature = "age_d1")
p1 <- ale$plot() + theme_bw()
ale$set.feature("age_d2")
p2 <- ale$plot() + theme_bw()
pdf("ages_effect-supp.pdf", width = 6.5, height = 4.5)
p1 + p2 + plot_layout(ncol = 2)
dev.off()


