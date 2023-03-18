### 5.3 random forest (but adds in parameter tuning for the RF and one using defaults)
### random split-sample and sample size
### sample sizes on performance variability
###

library(MASS)
library(doMC)
library(data.table)
library(glmnet)
library(ggh4x)
library(ranger)
library(tuneRanger)
library(rms)
library(ragg)
library(progress)
library(doParallel)
library(facetscales)
library(tidyverse)
library(dcurves)
library(apricom)
library(patchwork)
if(!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])


setwd('/Users/collinsg/Library/CloudStorage/OneDrive-Nexus365/CSM/research/Data/Gusto (Frank)/')
#setwd('/Users/gary/Library/CloudStorage/OneDrive-Nexus365/CSM/research/Data/Gusto (Frank)/')
DATA.ORIG <- read.csv('gusto.csv')
DATA.ORIG <- select(DATA.ORIG ,c("sex", "age", "hyp", "htn", "hrt", "pmi", "ste", "day30"))
names(DATA.ORIG)[8] <- 'y'

DATA.ORIG$sex2 <- rep(0, nrow(DATA.ORIG))
DATA.ORIG$sex2[DATA.ORIG$sex=='male'] <- 1

DATA.ORIG$pmi2 <- rep(0, nrow(DATA.ORIG))
DATA.ORIG$pmi2[DATA.ORIG$pmi == 'yes'] <- 1
DATA.ORIG <- select(DATA.ORIG ,c("sex2", "age", "hyp", "htn", "hrt", "pmi2", "ste", "y"))
set.seed(532523)
min.ss <- pmsampsize::pmsampsize(type = 'b', rsquared = 0.08, parameters=7, prevalence = 2851/40830)
DATA   <- DATA.ORIG %>% sample_n(min.ss$results_table[4,1])

NBOOT      <- 200
n.cal.grid <- 100

#################################
### Random Forest (DEFAULTS) ####
#################################
min.ss   <- pmsampsize::pmsampsize(type = 'b', rsquared = 0.08, parameters = 7, prevalence = 2851/40830)
DATA2    <- DATA
DATA2$y  <- factor(DATA$y, levels = c(0, 1))

fit.orig    <- ranger(y~., data = DATA2, probability = T, num.trees = 500)
pred.fit.rf <- predict(fit.orig, data = DATA2)$predictions[,2]
Sm.orig     <- loess(DATA$y~pred.fit.rf, degree = 2)
Sm.orig     <- data.frame(Sm.orig$x, Sm.orig$fitted)
Sm.orig     <- Sm.orig[order(Sm.orig$pred.fit.rf),]

########################################
### Random Forest (tree depth of 3) ####
########################################
fit.orig.d3    <- ranger(y~., data = DATA2, probability = T, num.trees = 500, max.depth = 3, min.node.size = 1)
pred.fit.rf.d3 <- predict(fit.orig.d3, data = DATA2)$predictions[,2]
Sm.orig.d3     <- loess(DATA$y~pred.fit.rf.d3, degree = 2)
Sm.orig.d3     <- data.frame(Sm.orig.d3$x, Sm.orig.d3$fitted)
Sm.orig.d3     <- Sm.orig.d3[order(Sm.orig.d3$pred.fit.rf),]

##############################
### Random Forest (TUNED) ####
##############################
DATA2.task   <- makeClassifTask(data = DATA2, target = "y")
DATA2.tuning <- tuneRanger(DATA2.task, 
                           measure         = list(auc), 
                           tune.parameters = c("mtry", "min.node.size"), 
                           show.info       = F, 
                           num.threads     = 16,
                           num.trees       = 500)

fit.orig.tune    <- DATA2.tuning$model
pred.fit.rf.tune <- predict(DATA2.tuning$model, newdata=DATA2)$data$prob.1
Sm.orig.tune     <- loess(DATA$y~pred.fit.rf.tune, degree = 2)
Sm.orig.tune     <- data.frame(Sm.orig.tune$x, Sm.orig.tune$fitted)
Sm.orig.tune     <- Sm.orig.tune[order(Sm.orig.tune$pred.fit.rf),]

### define some storage matrices
pred.rf            <- matrix(ncol = NBOOT, nrow = nrow(DATA))
pred.rf.d3         <- matrix(ncol = NBOOT, nrow = nrow(DATA))
pred.rf.tune       <- matrix(ncol = NBOOT, nrow = nrow(DATA))
cal.x.rf           <- matrix(ncol = NBOOT, nrow = 100)
cal.y.rf           <- matrix(ncol = NBOOT, nrow = 100)
cal.x.rf.d3        <- matrix(ncol = NBOOT, nrow = 100)
cal.y.rf.d3        <- matrix(ncol = NBOOT, nrow = 100)
cal.x.tune         <- matrix(ncol = NBOOT, nrow = 100)
cal.y.tune         <- matrix(ncol = NBOOT, nrow = 100)
cal.x.rf.tune      <- matrix(ncol = NBOOT, nrow = 100)
cal.y.rf.tune      <- matrix(ncol = NBOOT, nrow = 100)
cal.x.orig.rf      <- matrix(ncol = NBOOT, nrow = 100)
cal.y.orig.rf      <- matrix(ncol = NBOOT, nrow = 100)
cal.x.orig.rf.d3   <- matrix(ncol = NBOOT, nrow = 100)
cal.y.orig.rf.d3   <- matrix(ncol = NBOOT, nrow = 100)
cal.x.orig.rf.tune <- matrix(ncol = NBOOT, nrow = 100)
cal.y.orig.rf.tune <- matrix(ncol = NBOOT, nrow = 100)
cal.x.rf.orig      <- matrix(ncol = NBOOT, nrow = 100)
cal.y.rf.orig      <- matrix(ncol = NBOOT, nrow = 100)
cal.x.rf.orig.d3   <- matrix(ncol = NBOOT, nrow = 100)
cal.y.rf.orig.d3   <- matrix(ncol = NBOOT, nrow = 100)
cal.x.rf.tune.orig <- matrix(ncol = NBOOT, nrow = 100)
cal.y.rf.tune.orig <- matrix(ncol = NBOOT, nrow = 100)

pb <- progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = NBOOT, clear = FALSE, width = 60)

# random forest
for(i in 1:NBOOT){
  pb$tick()
  index <- sample(1:nrow(DATA2), nrow(DATA2), replace = T)
  DATA.boot <- DATA2[index,]
  
  if(sum(as.numeric(DATA.boot$y) - 1) == 0){ ### catch instances where bootstrapped data contains zero events
    i <- i - 1
  } else {
    ### parameter tuning in the bootstrap data
    DATA.task   <- makeClassifTask(data = DATA.boot, target = "y")
    DATA.tuning <- tuneRanger(DATA.task, 
                              measure         = list(auc), 
                              tune.parameters = c("mtry", "min.node.size"), 
                              show.info       = F, 
                              num.threads     = 16,
                              num.trees       = 500)
    
    ### fit bootstrap RF (defaults), bootstrap RF (max.depth=3), and tuned bootstrap RF
    fit      <- ranger(y~., data = DATA.boot, probability = T, num.trees = 500)
    fit.d3   <- ranger(y~., data = DATA.boot, probability = T, max.depth = 3, min.node.size = 1)
    fit.tune <- DATA.tuning$model
    
    ### get predictions from bootstrap RF (defaults), bootstrap RF (max.depth=3), and tuned bootstrap RF on original data
    pred.rf[,i]       <- predict(fit,    data = DATA2)$predictions[,2]
    pred.rf.d3[,i]    <- predict(fit.d3, data = DATA2)$predictions[,2]
    pred.rf.tune[,i]  <- predict(DATA.tuning$model, newdata=DATA2)$data$prob.1
    
    ### get calibration of the bootstrap RF (defaults) in original data
    Sm           <- lowess(pred.rf[,i], as.numeric(DATA2$y)-1, iter = 0)
    pp.rf        <- seq(min(pred.rf[,i]), max(pred.rf[,i]), length = n.cal.grid)
    Sm.rf        <- approx(Sm, xout = pp.rf, ties = mean)$y
    cal.x.rf.orig[,i] <- pp.rf      # calibration of bootstrap model on original data
    cal.y.rf.orig[,i] <- Sm.rf
    
    ### get calibration of the bootstrap RF (max.depth=3) in original data
    Sm.d3           <- lowess(pred.rf.d3[,i], as.numeric(DATA2$y)-1, iter = 0)
    pp.rf.d3        <- seq(min(pred.rf.d3[,i]), max(pred.rf.d3[,i]), length = n.cal.grid)
    Sm.rf.d3        <- approx(Sm.d3, xout = pp.rf.d3, ties = mean)$y
    cal.x.rf.orig.d3[,i] <- pp.rf.d3      # calibration of bootstrap model on original data
    cal.y.rf.orig.d3[,i] <- Sm.rf.d3
    
    ### get calibration of the tuned bootstrap RF in original data
    Sm           <- lowess(pred.rf.tune[,i], as.numeric(DATA2$y)-1, iter = 0)
    pp.rf        <- seq(min(pred.rf.tune[,i]), max(pred.rf.tune[,i]), length = n.cal.grid)
    Sm.rf        <- approx(Sm, xout = pp.rf, ties = mean)$y
    cal.x.rf.tune.orig[,i] <- pp.rf      # calibration of bootstrap model on original data
    cal.y.rf.tune.orig[,i] <- Sm.rf
    
    ### get predictions from bootstrap RF (defaults), bootstrap RF (max.depth=3), and tuned bootstrap RF on bootstrap data
    pred      <- predict(fit,    data = DATA.boot)$predictions[,2]
    pred.d3   <- predict(fit.d3, data = DATA.boot)$predictions[,2]
    pred.tune <- predict(DATA.tuning$model, newdata = DATA.boot)$data$prob.1
    
    ### get calibration of the bootstrap RF (defaults) in bootstrap data
    Sm           <- lowess(pred,      as.numeric(DATA.boot$y)-1, iter = 0)
    pp.rf        <- seq(min(pred), max(pred), length = n.cal.grid)
    Sm.rf        <- approx(Sm, xout = pp.rf, ties = mean)$y
    cal.x.rf[,i] <- pp.rf      # calibration of bootstrap model on bootstrap data
    cal.y.rf[,i] <- Sm.rf
    
    ### get calibration of the bootstrap RF (max.depth=3) in bootstrap data
    Sm.d3           <- lowess(pred.d3,   as.numeric(DATA.boot$y)-1, iter = 0)
    pp.rf.d3        <- seq(min(pred.d3), max(pred.d3), length = n.cal.grid)
    Sm.rf.d3        <- approx(Sm.d3, xout = pp.rf.d3, ties = mean)$y
    cal.x.rf.d3[,i] <- pp.rf.d3      # calibration of bootstrap model on bootstrap data
    cal.y.rf.d3[,i] <- Sm.rf.d3
    
    ### get calibration of the tuned bootstrap RF in bootstrap data
    Sm.tune           <- lowess(pred.tune, as.numeric(DATA.boot$y)-1, iter = 0)
    pp.rf.tune        <- seq(min(pred.tune), max(pred.tune), length = n.cal.grid)
    Sm.rf.tune        <- approx(Sm.tune, xout = pp.rf.tune, ties = mean)$y
    cal.x.rf.tune[,i] <- pp.rf.tune # calibration of bootstrap model on bootstrap data
    cal.y.rf.tune[,i] <- Sm.rf.tune
    
    # calibration of the original model (defaults) on the bootstrap data
    pred.val          <- predict(fit.orig, data = DATA.boot)$predictions[,2]
    Sm                <- lowess(pred.val, as.numeric(DATA.boot$y)-1, iter = 0)
    cal.x.orig.rf[,i] <- seq(min(pred.val), max(pred.val), length = 100)
    cal.y.orig.rf[,i] <- approx(Sm, xout = cal.x.orig.rf[,i], ties = mean)$y
    
    # calibration of the original model (max.depth=3) on the bootstrap data
    pred.val             <- predict(fit.orig.d3, data = DATA.boot)$predictions[,2]
    Sm                   <- lowess(pred.val, as.numeric(DATA.boot$y)-1, iter = 0)
    cal.x.orig.rf.d3[,i] <- seq(min(pred.val), max(pred.val), length = 100)
    cal.y.orig.rf.d3[,i] <- approx(Sm, xout = cal.x.orig.rf.d3[,i], ties = mean)$y
    
    # calibration of the original tuned RF on the bootstrap data
    #pred.val.tune          <- predict(fit.orig.tune, data = DATA.boot)$predictions[,2]
    pred.val.tune          <- predict(DATA2.tuning$model, newdata=DATA.boot)$data$prob.1
    Sm.tune                <- lowess(pred.val.tune, as.numeric(DATA.boot$y)-1, iter = 0)
    cal.x.orig.rf.tune[,i] <- seq(min(pred.val.tune), max(pred.val.tune), length = 100)
    cal.y.orig.rf.tune[,i] <- approx(Sm.tune, xout = cal.x.orig.rf.tune[,i], ties = mean)$y
  }
}

### Instability plot of individual predictions
x.rf      <- pred.fit.rf[order(pred.fit.rf)]
x.rf.d3   <- pred.fit.rf.d3[order(pred.fit.rf.d3)]
x.rf.tune <- pred.fit.rf.tune[order(pred.fit.rf.tune)]

y1.rf      <- apply(pred.rf,      1, function(x) quantile(x, probs = 0.025, na.rm = T))[order(pred.fit.rf)]
y2.rf      <- apply(pred.rf,      1, function(x) quantile(x, probs = 0.975, na.rm = T))[order(pred.fit.rf)]
y1.rf.d3   <- apply(pred.rf.d3,   1, function(x) quantile(x, probs = 0.025, na.rm = T))[order(pred.fit.rf.d3)]
y2.rf.d3   <- apply(pred.rf.d3,   1, function(x) quantile(x, probs = 0.975, na.rm = T))[order(pred.fit.rf.d3)]
y1.rf.tune <- apply(pred.rf.tune, 1, function(x) quantile(x, probs = 0.025, na.rm = T))[order(pred.fit.rf.tune)]
y2.rf.tune <- apply(pred.rf.tune, 1, function(x) quantile(x, probs = 0.975, na.rm = T))[order(pred.fit.rf.tune)]

OUT2 <- data.frame(x = c(x.rf,  x.rf,  x.rf.d3,  x.rf.d3,  x.rf.tune,  x.rf.tune), 
                   y = c(y1.rf, y2.rf, y1.rf.d3, y2.rf.d3, y1.rf.tune, y2.rf.tune))

xx1.rf      <- lowess(y1.rf~x.rf,           delta = 0.01, f = 0.1)
xx2.rf      <- lowess(y2.rf~x.rf,           delta = 0.01, f = 0.1)
xx1.rf.d3   <- lowess(y1.rf.d3~x.rf.d3,     delta = 0.01, f = 0.1)
xx2.rf.d3   <- lowess(y2.rf.d3~x.rf.d3,     delta = 0.01, f = 0.1)
xx1.rf.tune <- lowess(y1.rf.tune~x.rf.tune, delta = 0.01, f = 0.1)
xx2.rf.tune <- lowess(y2.rf.tune~x.rf.tune, delta = 0.01, f = 0.1)

OUT3 <- data.frame(x = c(xx1.rf$x, xx1.rf.d3$x, xx1.rf.tune$x, xx2.rf$x, xx2.rf.d3$x, xx2.rf.tune$x), 
                   y = c(xx1.rf$y, xx1.rf.d3$y, xx1.rf.tune$y, xx2.rf$y, xx2.rf.d3$y, xx2.rf.tune$y))

OUT3$N.DEV <- c(rep("rf",       length(xx1.rf$x)),
                rep("rf.d3",    length(xx1.rf.d3$x)),
                rep("tuned rf", length(xx1.rf.tune$x)),
                rep("rf",       length(xx2.rf$x)),
                rep("rf.d3",    length(xx2.rf.d3$x)),
                rep("tuned rf", length(xx2.rf.tune$x))) 

OUT3$N.DEV <- factor(OUT3$N.DEV, 
                     levels = c("rf", "rf.d3", "tuned rf"),
                     labels = c("Random Forest (defaults)", "Random Forest (max.depth=3)", "tuned Random Forest"))
OUT3$limit <- c(rep("lower", length(xx1.rf$x) + length(xx1.rf.d3$x) + length(xx1.rf.tune$x)), 
                rep("upper", length(xx2.rf$x) + length(xx2.rf.d3$x) + length(xx2.rf.tune$x)))
OUT3$limit <- factor(OUT3$limit, levels = c("lower", "upper"))

pred.rf.long      <- reshape2::melt(pred.rf)
pred.rf.d3.long   <- reshape2::melt(pred.rf.d3)
pred.rf.tune.long <- reshape2::melt(pred.rf.tune)
names(pred.rf.long)[2]      <- "B"
names(pred.rf.long)[3]      <- "boot.pred"
names(pred.rf.d3.long)[2]   <- "B"
names(pred.rf.d3.long)[3]   <- "boot.pred"
names(pred.rf.tune.long)[2] <- "B"
names(pred.rf.tune.long)[3] <- "boot.pred"

pred.rf.long      <- add_column(pred.rf.long,      mod.pred = rep(pred.fit.rf, NBOOT))
pred.rf.long      <- add_column(pred.rf.long,      N.DEV    = rep("rf", NBOOT * nrow(pred.rf)))
pred.rf.d3.long   <- add_column(pred.rf.d3.long,   mod.pred = rep(pred.fit.rf.d3, NBOOT))
pred.rf.d3.long   <- add_column(pred.rf.d3.long,   N.DEV    = rep("rf.d3", NBOOT * nrow(pred.rf.d3)))
pred.rf.tune.long <- add_column(pred.rf.tune.long, mod.pred = rep(pred.fit.rf.tune, NBOOT))
pred.rf.tune.long <- add_column(pred.rf.tune.long, N.DEV    = rep("tuned rf", NBOOT * nrow(pred.rf.tune)))

OUT       <- bind_rows(pred.rf.long, pred.rf.d3.long, pred.rf.tune.long)
OUT$N.DEV <- factor(OUT$N.DEV, 
                    levels = c("rf", "rf.d3", "tuned rf"),
                    labels = c("Random Forest (defaults)", "Random Forest (max.depth=3)", "tuned Random Forest"))

p1 <- OUT  %>% ggplot(aes(x = mod.pred, y = boot.pred)) +
  geom_point(size = 0.1, alpha = 0.5, colour='grey') +
  geom_line(data = OUT3, aes(x=x, y=y, group = limit), colour='black', linetype=2) + 
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(~N.DEV) + 
  xlim(0, 1) +
  ylim(0, 1) +
  xlab('Estimated risk from the developed model') +
  ylab('Estimated risk in the bootstrap samples') +
  theme_bw() +
  theme(axis.text = element_text(size = 6))


### Instability plot of MAPE values
par(mfrow = c(1, 2))
plot(pred.fit.rf,      apply(abs(pred.rf      - pred.fit.rf),      1, median, na.rm = T), pch = 20, xlab = 'Estimated risk from developed model', ylab = 'Instability index')
plot(pred.fit.rf.d3,   apply(abs(pred.rf.d3   - pred.fit.rf.d3),   1, median, na.rm = T), pch = 20, xlab = 'Estimated risk from developed model', ylab = 'Instability index')
plot(pred.fit.rf.tune, apply(abs(pred.rf.tune - pred.fit.rf.tune), 1, median, na.rm = T), pch = 20, xlab = 'Estimated risk from developed model', ylab = 'Instability index')

OUT.ii <- data.frame(pred  = c(pred.fit.rf, pred.fit.rf.d3, pred.fit.rf.tune),
                     index = c(apply(abs(pred.rf      - pred.fit.rf),      1, median, na.rm = T),
                               apply(abs(pred.rf.d3   - pred.fit.rf.d3),   1, median, na.rm = T),
                               apply(abs(pred.rf.tune - pred.fit.rf.tune), 1, median, na.rm = T)))

OUT.ii$label <- c(rep("rf",       length(pred.fit.rf)),
                  rep("rf.d3",    length(pred.fit.rf.d3)),
                  rep("tuned rf", length(pred.fit.rf.tune)))

OUT.ii$label <- factor(OUT.ii$label, 
                       levels = c("rf", "rf.d3", "tuned rf"),
                       labels = c("Random Forest (defauts)", "Random Forest (max.depth=3)", "tuned Random Forest"))

p3 <- OUT.ii  %>% ggplot(aes(x = pred, y = index)) +
  geom_point(size = 0.1, alpha = 0.5, colour = 'black') +
  facet_grid(~label) + 
  xlab('Estimated risk from the developed model') +
  ylab('Instability index') +
  theme_bw() +
  theme(axis.text = element_text(size = 6))

### MAPE_individual
MAPEi <- function(X, y){
  mean(abs(X-y))
}
mean(apply(pred.rf,      2, MAPEi, y = pred.fit.rf))
mean(apply(pred.rf.d3,   2, MAPEi, y = pred.fit.rf.d3))
mean(apply(pred.rf.tune, 2, MAPEi, y = pred.fit.rf.tune))

# plot calibration of the original model on the bootstrap data
cal.rf.orig.xx     <- reshape2::melt(cal.x.rf.orig)
cal.rf.orig.yy     <- reshape2::melt(cal.y.rf.orig)

cal.rf.d3.orig.xx   <- reshape2::melt(cal.x.rf.orig.d3)
cal.rf.d3.orig.yy   <- reshape2::melt(cal.y.rf.orig.d3)

cal.rf.tune.orig.xx <- reshape2::melt(cal.x.rf.tune.orig)
cal.rf.tune.orig.yy <- reshape2::melt(cal.y.rf.tune.orig)

cal.rf.orig.xx$N      <- rep("rf",       nrow(cal.rf.orig.xx))
cal.rf.d3.orig.xx$N   <- rep("rf.d3",    nrow(cal.rf.d3.orig.xx))
cal.rf.tune.orig.xx$N <- rep("tuned rf", nrow(cal.rf.tune.orig.xx))

cal <- rbind(cal.rf.orig.xx, cal.rf.d3.orig.xx, cal.rf.tune.orig.xx)

colnames(cal) <- c("Var1", "Sim", "value", "N")
cal$N         <- factor(cal$N, 
                        levels = c("rf", "rf.d3", "tuned rf"),
                        labels = c("Random Forest (defaults)", "Random Forest (max.depth=3)", "tuned Random Forest"))

OUT.cal <- data.table(cal)
OUT.cal[, y:=c(cal.rf.orig.yy$value,
               cal.rf.d3.orig.yy$value,
               cal.rf.tune.orig.yy$value)]

p2 <- OUT.cal %>% ggplot(aes(x = value, y = y, group = Sim)) +
  geom_line(alpha = 0.5, colour = '#999999') + 
  facet_grid(~N) + 
  xlim(0, 1) + 
  ylim(0, 1) +
  xlab('Estimated risk from the developed model') +
  ylab('Observed risk in the bootstrap samples') +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  theme(axis.text = element_text(size = 6))

p1/p2/p3

