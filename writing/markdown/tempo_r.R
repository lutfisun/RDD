######

library(learnr)
library(haven)
library(tidyverse)
library(stargazer)
library(estimatr)
library(rdd)
library(rdrobust)
library(rddensity)
library(knitr)
library(KernSmooth)
library(locfit)
library(gridExtra)

######

hansen <- read_dta("../../data/hansen_dwi.dta")
hansen$dui08 <- ifelse(hansen$bac1 >= 0.08, 1, 0)

######

p_mccrary = DCdensity(hansen$bac1, 
                      cutpoint = 0.08,
                      plot = FALSE)
density <- rddensity(hansen$bac1, c = 0.08, p=2, h=0.02)
rdplotdensity(density, 
              hansen$bac1,
              plotN = 10,
              histBreaks= seq(from = 0, to =0.3 , by = 0.0039),
              histFillShade = 0.8)

#####

density <- rddensity(hansen$bac1, c = 0.08, p=2, h=0.01)
rdplotdensity(density, 
              hansen$bac1,
              plotN = 10,
              histBreaks= seq(from = 0, to =0.3 , by = 0.0039),
              histFillShade = 0.8)

#####

hist(hansen$bac1,
     main="Histogram for BAC Levels Recorded in DUI Stops", 
     xlab="BAC",
     breaks = 40000)
abline(v = 0.08, col= 'red')

#####

lincb_male  <- RDestimate(male ~ bac1, data = hansen, 
                          cutpoint = 0.08, bw = 0.05, kernel = "rectangular")
lincb_white <- RDestimate(white ~ bac1, data = hansen, 
                          cutpoint = 0.08, bw = 0.05, kernel = "rectangular")
lincb_age   <- RDestimate(aged ~ bac1, data = hansen, 
                          cutpoint = 0.08, bw = 0.05, kernel = "rectangular")
lincb_acc   <- RDestimate(acc ~ bac1, data = hansen, 
                          cutpoint = 0.08, bw = 0.05, kernel = "rectangular")

late_male  = lincb_male$est[1]
late_white = lincb_white$est[1]
late_age   = lincb_age$est[1]
late_acc   = lincb_acc$est[1]

p_male  = lincb_acc$p[1]
p_white = lincb_white$p[1]
p_age   = lincb_age$p[1]
p_acc   = lincb_acc$p[1]

######

hansen$bac_centered <- hansen$bac1 - 0.08 
weights <- rdd::kernelwts(hansen$bac_centered, 
                          center = 0, 
                          bw = .05, 
                          kernel = "rectangular")

lmlin_male  <- lm_robust(male ~  bac_centered + dui08 + bac_centered*dui08, 
                         data = hansen, weights = weights)
lmlin_white <- lm_robust(white ~ bac_centered + dui08 + bac_centered*dui08, 
                         data = hansen, weights = weights)
lmlin_aged  <- lm_robust(aged ~ bac_centered + dui08 + bac_centered*dui08, 
                         data = hansen, weights = weights)
lmlin_acc   <- lm_robust(acc ~ bac_centered + dui08 + bac_centered*dui08, 
                         data = hansen, weights = weights)
texreg::screenreg(list(lmlin_male, lmlin_white, lmlin_aged, lmlin_acc),
                  custom.model.names = c("Male", "White", "Age", "Accident"))

######

rdplot(x = hansen$bac1, y = hansen$male, c = 0.08, p=1, kernel = "triangular", 
       title = "Male Cov Balance Linear",y.label = "Male", x.label = "BAC")
rdplot(x = hansen$bac1, y = hansen$male, c = 0.08, p=2, kernel = "triangular",
       title = "Male Cov Balance Quadratic",y.label = "Male", x.label = "BAC")

rdplot(x = hansen$bac1, y = hansen$white, c = 0.08, p=1, kernel = "triangular",
       title = "White Cov Balance Linear",y.label = "White", x.label = "BAC")
rdplot(x = hansen$bac1, y = hansen$white, c = 0.08, p=2, kernel = "triangular",
       title = "White Cov Balance Quadratic",y.label = "White", x.label = "BAC")


rdplot(x = hansen$bac1, y = hansen$aged, c = 0.08, p=1, kernel = "triangular",
       title = "Age Cov Balance Linear",y.label = "Age", x.label = "BAC")
rdplot(x = hansen$bac1, y = hansen$aged, c = 0.08, p=2, kernel = "triangular",
       title = "Age Cov Balance Quadratic",y.label = "Age", x.label = "BAC")

rdplot(x = hansen$bac1, y = hansen$acc, c = 0.08, p=1, kernel = "triangular",
       title = "Accident Cov Balance Linear",y.label = "Accident", x.label = "BAC")
rdplot(x = hansen$bac1, y = hansen$acc, c = 0.08, p=2, kernel = "triangular",
       title = "Accident Cov Balance Quadratic",y.label = "Accident", x.label = "BAC")

######

hansen$bac_centered2 <- hansen$bac_centered * hansen$bac_centered

sq_male  <- lm_robust(male ~  bac_centered + bac_centered2 + dui08 +
                         bac_centered*dui08 + bac_centered2*dui08, 
                       data = hansen, weights = weights)
sq_white <- lm_robust(white ~ bac_centered + bac_centered2 + dui08 +
                        bac_centered*dui08 + bac_centered2*dui08, 
                       data = hansen, weights = weights)
sq_age   <- lm_robust(aged ~ bac_centered + bac_centered2 + dui08 +
                        bac_centered*dui08 + bac_centered2*dui08, 
                       data = hansen, weights = weights)
sq_acc   <- lm_robust(acc ~ bac_centered + bac_centered2 + dui08 +
                        bac_centered*dui08 + bac_centered2*dui08, 
                       data = hansen, weights = weights)
texreg::screenreg(list(sq_male, sq_white, sq_age, sq_acc), 
                  custom.model.names = c("Male", "White", "Age", "Accident"))

#####

grid_bac <- hansen %>%
  mutate(bacbin = cut(bac1, breaks = 
                      seq(from = -0.01, to = 4.01, by = 0.01))) %>%
  group_by(bacbin) %>%
  summarize(bac1 = mean(bac1, na.rm = TRUE),
            mean_recid = mean(recidivism, na.rm = TRUE),
            mean_male = mean(male, na.rm = TRUE),
            mean_white = mean(white, na.rm = TRUE),
            mean_age = mean(aged, na.rm = TRUE),
            mean_acc = mean(acc, na.rm = TRUE))

grid_bac$dui08 <- ifelse(grid_bac$bac1 >= 0.08, 1, 0)
grid_bac$bac_centered  <- grid_bac$bac1 - 0.08
grid_bac$bac_centered2 <- grid_bac$bac_centered * grid_bac$bac_centered

pred_male <- predict(sq_male, grid_bac, interval="prediction")$fit
pred_white <- predict(sq_white, grid_bac, interval="prediction")$fit
pred_age <- predict(sq_age, grid_bac, interval="prediction")$fit
pred_acc <- predict(sq_acc, grid_bac, interval="prediction")$fit

#####
par(mfrow=c(2,2))

plot(grid_bac$bac1, grid_bac$mean_male,
     main = "Panel A. Male",
     xlab = "bac level", ylab = "percentage male")
polygon(c(rev(grid_bac$bac1), grid_bac$bac1), 
        c(rev(pred_male[,3]), pred_male[,2]), 
        col=adjustcolor("grey80",alpha.f=0.5), border = NA)
lines(grid_bac$bac1, pred_male[,1], col=4)
lines(grid_bac$bac1, pred_male[,2], lty = 'dashed', col=2)
lines(grid_bac$bac1, pred_male[,3], lty = 'dashed', col=2)
abline(v = 0.08, col= 1)

plot(grid_bac$bac1, grid_bac$mean_white,
     main = "Panel B. White",
     xlab = "bac level", ylab = "percentage white")
polygon(c(rev(grid_bac$bac1), grid_bac$bac1), 
        c(rev(pred_white[ ,3]), pred_white[ ,2]), 
        col=adjustcolor("grey80",alpha.f=0.5), border = NA)
lines(grid_bac$bac1, pred_white[,1], col=4)
lines(grid_bac$bac1, pred_white[,2], lty = 'dashed', col=2)
lines(grid_bac$bac1, pred_white[,3], lty = 'dashed', col=2)
abline(v = 0.08, col= 1)

plot(grid_bac$bac1, grid_bac$mean_age,
     main = "Panel C. Age",
     xlab = "bac level", ylab = "average age")
polygon(c(rev(grid_bac$bac1), grid_bac$bac1), 
        c(rev(pred_age[ ,3]), pred_age[ ,2]), 
        col=adjustcolor("grey80",alpha.f=0.5), border = NA)
lines(grid_bac$bac1, pred_age[,1], col=4)
lines(grid_bac$bac1, pred_age[,2], lty = 'dashed', col=2)
lines(grid_bac$bac1, pred_age[,3], lty = 'dashed', col=2)
abline(v = 0.08, col= 1)

plot(grid_bac$bac1, grid_bac$mean_acc,
     main = "Panel D. Accident",
     xlab = "bac level", ylab = "percentage male")
polygon(c(rev(grid_bac$bac1), grid_bac$bac1), 
        c(rev(pred_acc[ ,3]), pred_acc[ ,2]), 
        col=adjustcolor("grey80",alpha.f=0.5), border = NA)
lines(grid_bac$bac1, pred_acc[,1], col=4)
lines(grid_bac$bac1, pred_acc[,2], lty = 'dashed', col=2)
lines(grid_bac$bac1, pred_acc[,3], lty = 'dashed', col=2)
abline(v = 0.08, col= 1)

#######

plot(hansen2$bac1, hansen2$mean_acc)
lines(locpoly(hansen2$mean_acc, hansen2$bac1, bandwidth=0.05, degree=1), col=2)
lines(locpoly(hansen2$mean_acc, hansen2$bac1, bandwidth=0.05, degree=2), col=3)
lines(locfit(hansen$acc ~ lp(hansen$bac1, nn=0, h=0.05, deg=1)), col=5)

plot(hansen2$bac1, hansen2$mean_male)
lines(locpoly(hansen$male, hansen$bac1, bandwidth=0.05, degree=1), col=2)
lines(locpoly(hansen$male, hansen$bac1, bandwidth=0.05, degree=2), col=3)
lines(locfit(hansen$male ~ lp(hansen$bac1, nn=0, h=0.05, deg=1)), col=5)

#####

panelA_03 <- subset(hansen, bac1 >= 0.03  & bac1 <= 0.13)
panelB_05 <- subset(hansen, bac1 >= 0.055 & bac1 <= 0.105)

#Panel A

A_weights <- rdd::kernelwts(panelA_03$bac_centered, center = 0, 
                            bw = .05, kernel = "rectangular")

A_recid_lin0  <- lm_robust(recidivism ~  bac_centered, 
                           data = panelA_03, weights = A_weights)
A_recid_lin1 <- lm_robust(recidivism ~ bac_centered + dui08 + bac_centered*dui08,
                          data = panelA_03, weights = A_weights)
A_recid_sq   <- lm_robust(recidivism ~ bac_centered + bac_centered2 + dui08 +
                          bac_centered*dui08 + bac_centered2*dui08, 
                          data = panelA_03, weights = A_weights)

texreg::screenreg(list(A_recid_lin0, A_recid_lin1, A_recid_sq), 
                  custom.model.names = c("A. Linear", "A. Line + Cutoff", "A. Quad + Cutoff"))

# Panel B

B_weights <- rdd::kernelwts(panelB_05$bac_centered, center = 0, 
                            bw = .05, kernel = "rectangular")

B_recid_lin0  <- lm_robust(recidivism ~  bac_centered, 
                           data = panelB_05, weights = B_weights)
B_recid_lin1 <- lm_robust(recidivism ~ bac_centered + dui08 + bac_centered*dui08,
                          data = panelB_05, weights = B_weights)
B_recid_sq   <- lm_robust(recidivism ~ bac_centered + bac_centered2 + dui08 +
                            bac_centered*dui08 + bac_centered2*dui08, 
                          data = panelB_05, weights = B_weights)

texreg::screenreg(list(B_recid_lin0, B_recid_lin1, B_recid_sq), 
                  custom.model.names = c("B. Linear", "B. Line + Cutoff", "B. Quad + Cutoff"))

######

hansen_q8 <- subset(hansen, bac1 < 0.15)

rdplot(x = hansen_q8$bac1, y = hansen_q8$recidivism, c = 0.08, p=1, kernel = "triangular",
       title = "Impact of 0.08 BAC cutoff on Recidivism (Linear Trend)",y.label = "Recidivism", x.label = "BAC")

rdplot(x = hansen_q8$bac1, y = hansen_q8$recidivism, c = 0.08, p=2, kernel = "triangular",
       title = "Impact of 0.08 BAC cutoff on Recidivism (Quadratic Trend)",y.label = "Recidivism", x.label = "BAC")

######

grid_bac <- subset(grid_bac, bac1 < 0.15)

q8_weights <- rdd::kernelwts(hansen_q8$bac_centered, center = 0, 
                            bw = .05, kernel = "rectangular")

q8_recid_lin1 <- lm_robust(recidivism ~ bac_centered + dui08 + bac_centered*dui08,
                          data = hansen_q8, weights = q8_weights)

q8_recid_sq   <- lm_robust(recidivism ~ bac_centered + bac_centered2 + dui08 + 
                             bac_centered*dui08 + bac_centered2*dui08, 
                           data = hansen_q8, weights = q8_weights)

predlin_recid <- predict(q8_recid_lin1, grid_bac, interval="prediction")$fit
predsq_recid  <- predict(q8_recid_sq, grid_bac, interval="prediction")$fit


par(mfrow=c(2,1))

plot(grid_bac$bac1, grid_bac$mean_recid,
     main = "Panel A. Recidivism with Linear Trend",
     xlab = "bac level", ylab = "ratio of recidivism")
polygon(c(rev(grid_bac$bac1), grid_bac$bac1), 
        c(rev(predlin_recid[,3]), predlin_recid[,2]), 
        col=adjustcolor("grey80",alpha.f=0.5), border = NA)
lines(grid_bac$bac1, predlin_recid[,1], col=4)
lines(grid_bac$bac1, predlin_recid[,2], lty = 'dashed', col=2)
lines(grid_bac$bac1, predlin_recid[,3], lty = 'dashed', col=2)
abline(v = 0.08, col= 1)

plot(grid_bac$bac1, grid_bac$mean_recid,
     main = "Panel B. Recidivism with Quadratic Trend",
     xlab = "bac level", ylab = "ratio of recidivism")
polygon(c(rev(grid_bac$bac1), grid_bac$bac1), 
        c(rev(predsq_recid[ ,3]), predsq_recid[ ,2]), 
        col=adjustcolor("grey80",alpha.f=0.5), border = NA)
lines(grid_bac$bac1, predsq_recid[,1], col=4)
lines(grid_bac$bac1, predsq_recid[,2], lty = 'dashed', col=2)
lines(grid_bac$bac1, predsq_recid[,3], lty = 'dashed', col=2)
abline(v = 0.08, col= 1)

# I learned that not the slope b












  