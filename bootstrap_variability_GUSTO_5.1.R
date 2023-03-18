### 5.1 Unpenalised logistic regression forcing in 7 predictors
### random split-sample and sample size
### sample sizes on performance variability

rm(list = ls())
library(MASS)
library(doMC)
library(data.table)
library(glmnet)
library(ggh4x)
library(ranger)
library(rms)
library(ragg)
library(progress)
library(doParallel)
library(facetscales)
library(tidyverse)
library(dcurves)
library(apricom)
library(patchwork)
#library(doParallel)
#registerDoParallel(cl=4, cores=4)
#dev.off(dev.list()["RStudioGD"])

### define some functions
# evaluate model performance (c-statistic, R2, calibration slope and intercept)
val.perf     <- function(data = DATA.val, fit = fit.full){
  fit.val       <- lrm(data$y~predict(fit, data))
  cal.slope     <- as.numeric(coef(fit.val)[2])
  fit.cal       <- glm(data$y~offset(predict(fit, data)), family = 'binomial')
  cal.intercept <- as.numeric(coef(fit.cal)[1])
  c.val         <- as.numeric(fit.val$stats['C'])
  R2.val        <- as.numeric(fit.val$stats['R2'])
  return(list(c.val = c.val, 
              R2.val = R2.val, 
              cal.intercept = cal.intercept,
              cal.slope = cal.slope))
}
# apply the shrinkage factor to the model coefficients and re-estimate the intercept
recalibrate  <- function(shrinkage, fit, data = DATA){
  shrunk.coef   <- shrinkage * coef(fit)
  fit.formula   <- formula(fit)
  fit.shrunk    <- lrm(fit.formula, data = data, init = shrunk.coef, maxit = 1)
  lrm.offset    <- lrm.fit(y = data$y, offset = predict(fit.shrunk, data, type = 'lp'))
  new.intercept <- fit.shrunk$coefficients[1] + lrm.offset$coefficients[1]
  fit.shrunk    <- lrm(fit.formula, data = data, init = c(new.intercept, shrunk.coef[-1]), maxit = 1)
  return(fit.shrunk)
}
full_model   <- function(data.d = DATA, B = B, data.v = DATA.val, n.cal.grid = 100){
  ### fit maximum likelihood full model (no variable selection)
  fit.full <- lrm(y ~ ., data = data.d, x = T, y = T)
  
  if(!fit.full$fail){
    ### get apparent performance of full model (full.apparent)
    apparent <- as.numeric(val.perf(data = data.d, fit = fit.full))
    val      <- as.numeric(val.perf(data = data.v, fit = fit.full))
    
    # bootstrap internal validation of the full model to get the optimism-corrected bootstrap
    # calibration slope
    val.full  <- rms::validate(fit.full, bw = F, B = B, pr = F, estimates = F)
    xx        <- c(as.vector(val.full[c("Dxy", "R2", "Intercept","Slope"), 'index.corrected']))
    boot.opt  <- c((1 + xx[1])/2, xx[-1])
    
    ## apply the bootstrap corrected slope shrinkage factor (and re-estimate the intercept)
    fit.shrunk.full <- recalibrate(shrinkage = boot.opt[4], fit = fit.full, data = data.d)
    shrunk.apparent <- as.numeric(val.perf(data = data.d, fit = fit.shrunk.full))
    shrunk.val      <- as.numeric(val.perf(data = data.v, fit = fit.shrunk.full))
    
    pred.val <- predict(fit.full, newdata = data.v, type = 'fitted')
    Sm       <- lowess(pred.val, data.v$y, iter = 0)
    pp.full  <- seq(min(pred.val), max(pred.val), length = n.cal.grid)
    Sm.full  <- approx(Sm, xout = pp.full, ties = mean)$y
    
    pred.shrunk <- predict(fit.shrunk.full, newdata = data.v, type = 'fitted')
    Sm          <- lowess(pred.shrunk, data.v$y, iter = 0)
    pp.shrunk   <- seq(min(pred.shrunk), max(pred.shrunk), length = n.cal.grid)
    Sm.shrunk   <- approx(Sm, xout = pp.shrunk, ties = mean)$y 
  } else {
    apparent <- NA
    val      <- NA
    boot.opt  <- NA
    shrunk.apparent <- NA
    shrunk.val <- NA
    pp.full <- rep(NA, n.cal.grid)
    Sm.full <- rep(NA, n.cal.grid)
    pp.shrunk <- rep(NA, n.cal.grid)
    Sm.shrunk <- rep(NA, n.cal.grid)
    fit.full <- NA
    fit.shrunk.full <- NA
    boot.opt[4] <- NA
  }
  
  return(list(apparent        = apparent, 
              val             = val, 
              opt.cor         = boot.opt,
              shrunk.apparent = shrunk.apparent,
              shrunk.val      = shrunk.val,
              fit.cal.x       = pp.full,
              fit.cal.y       = Sm.full,
              shrunk.cal.x    = pp.shrunk,
              shrunk.cal.y    = Sm.shrunk,
              fit             = fit.full,
              fit.shrunk      = fit.shrunk.full, 
              shrinkage.f     = boot.opt[4]))
}

setwd('/Users/collinsg/Library/CloudStorage/OneDrive-Nexus365/CSM/research/Data/Gusto (Frank)/')
DATA.ORIG <- read.csv('gusto.csv')
DATA.ORIG <- DATA.ORIG[,c("sex", "age", "hyp", "htn", "hrt", "pmi", "ste", "day30")]
names(DATA.ORIG)[8] <- 'y'

NBOOT <- 200

### N=LARGE
DATA <- DATA.ORIG
fit.orig      <- full_model(data.d = DATA, B = 100, data.v = DATA)
pred.fit.large  <- plogis(predict(fit.orig$fit))

pred.large       <- matrix(ncol = NBOOT, nrow = nrow(DATA))
cal.x.large      <- matrix(ncol = NBOOT, nrow = 100)
cal.y.large      <- matrix(ncol = NBOOT, nrow = 100)
cal.x.orig.large <- matrix(ncol = NBOOT, nrow = 100)
cal.y.orig.large <- matrix(ncol = NBOOT, nrow = 100)

pb <- progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = NBOOT, clear = FALSE, width = 60)

for(i in 1:NBOOT){
  pb$tick()
  index <- sample(1:nrow(DATA), nrow(DATA), replace = T)
  DATA.boot <- DATA[index,]
  if(sum(DATA.boot$y) == 0){ ### catch instances where bootstrapped data contains zero events
    i <- i-1
  } else {
    fit <- full_model(data.d = DATA.boot, B = NBOOT, data.v = DATA)
    #if(!is.na(fit$fit)){
    if(any(sqrt(diag(vcov(fit$fit)))>30)) {
      i <- i - 1  
    } else {
      pred.large[,i]  <- plogis(predict(fit$fit, newdata = DATA))
      cal.x.large[,i] <- fit$fit.cal.x # calibraiton of bootstrap model on boostrap data
      cal.y.large[,i] <- fit$fit.cal.y
      
      # calibration of the original model on the bootstrap data
      pred.val <- predict(fit.orig$fit, newdata = DATA.boot, type = 'fitted')
      Sm       <- lowess(pred.val, DATA.boot$y, iter = 0)
      #Sm <- loess(DATA.boot$y~pred.val, degree = 2)
      #Sm <- data.frame(Sm$x, Sm$fitted)
      #Sm <- Sm[order(Sm$pred.val), ]
      cal.x.orig.large[,i] <- seq(min(pred.val), max(pred.val), length = 100)
      cal.y.orig.large[,i] <- approx(Sm, xout = cal.x.orig.large[,i], ties = mean)$y
    }
    #}
  }
}

### N=752 (MINIMUM)
set.seed(532523)
min.ss <- pmsampsize::pmsampsize(type = 'b', rsquared = 0.08, parameters=7, prevalence = 2851/40830)
DATA         <- DATA.ORIG %>% sample_n(min.ss$results_table[4,1])
fit.orig     <- full_model(data.d = DATA, B = 100, data.v = DATA)
pred.fit.min <- plogis(predict(fit.orig$fit))
Sm.orig      <- loess(DATA$y~pred.fit.min, degree=2)
Sm.orig      <- data.frame(Sm.orig$x, Sm.orig$fitted)
Sm.orig      <- Sm.orig[order(Sm.orig$pred.fit.min),]

pred.min       <- matrix(ncol = NBOOT, nrow = nrow(DATA))
cal.x.min      <- matrix(ncol = NBOOT, nrow = 100)
cal.y.min      <- matrix(ncol = NBOOT, nrow = 100)
cal.x.orig.min <- matrix(ncol = NBOOT, nrow = 100)
cal.y.orig.min <- matrix(ncol = NBOOT, nrow = 100)

pb <- progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = NBOOT, clear = FALSE, width = 60)

for(i in 1:NBOOT){
  pb$tick()
  index <- sample(1:nrow(DATA), nrow(DATA), replace = T)
  DATA.boot <- DATA[index,]
  if(sum(DATA.boot$y) == 0){ ### catch instances where bootstrapped data contains zero events
    i <- i-1
  } else {
    fit <- full_model(data.d = DATA.boot, B = NBOOT, data.v = DATA)
    if(any(sqrt(diag(vcov(fit$fit)))>30)) {
      i <- i - 1  
    } else {
      pred.min[,i]  <- plogis(predict(fit$fit, newdata = DATA))
      cal.x.min[,i] <- fit$fit.cal.x # calibration of bootstrap model on boostrap data
      cal.y.min[,i] <- fit$fit.cal.y
      
      # calibration of the original model on the bootstrap data
      pred.val <- predict(fit.orig$fit, newdata = DATA.boot, type = 'fitted')
      Sm       <- lowess(pred.val, DATA.boot$y, iter = 0)
      #Sm <- loess(DATA.boot$y~pred.val, degree=2)
      #Sm <- data.frame(Sm$x, Sm$fitted)
      #Sm <- Sm[order(Sm$pred.val),]
      cal.x.orig.min[,i] <- seq(min(pred.val), max(pred.val), length = 100)
      cal.y.orig.min[,i] <- approx(Sm, xout = cal.x.orig.min[,i], ties = mean)$y
    }
  }
}

### N=300 (SMALL)
set.seed(657222)
DATA           <- DATA.ORIG %>% sample_n(300)
fit.orig       <- full_model(data.d = DATA, B = 100, data.v = DATA)
pred.fit.small <- plogis(predict(fit.orig$fit))

pred.small       <- matrix(ncol = NBOOT, nrow = nrow(DATA))
cal.x.small      <- matrix(ncol = NBOOT, nrow = 100)
cal.y.small      <- matrix(ncol = NBOOT, nrow = 100)
cal.x.orig.small <- matrix(ncol = NBOOT, nrow = 100)
cal.y.orig.small <- matrix(ncol = NBOOT, nrow = 100)

pb <- progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = NBOOT, clear = FALSE, width = 60)

for(i in 1:NBOOT){
  pb$tick()
  index <- sample(1:nrow(DATA), nrow(DATA), replace = T)
  DATA.boot <- DATA[index,]
  if(sum(DATA.boot$y) == 0){ ### catch instances where bootstrapped data contains zero events
    i <- i-1
  } else {
    fit <- full_model(data.d = DATA.boot, B = NBOOT, data.v = DATA)
    #if(!is.na(fit$fit)){
    if(any(sqrt(diag(vcov(fit$fit)))>30)) {
      i <- i - 1  
    } else {
      pred.small[,i]  <- plogis(predict(fit$fit, newdata = DATA))
      pred.val <- predict(fit$fit, newdata = DATA.boot, type = 'fitted')
      Sm       <- lowess(pred.val, DATA.boot$y, iter = 0)
      pp.full  <- seq(min(pred.val), max(pred.val), length = 100)
      Sm.full  <- approx(Sm, xout = pp.full, ties = mean)$y
      cal.x.small[,i] <- fit$fit.cal.x # calibration of bootstrap model on boostrap data
      cal.y.small[,i] <- fit$fit.cal.y
      
      # calibration of the original model on the bootstrap data
      pred.val <- predict(fit.orig$fit, newdata = DATA.boot, type = 'fitted')
      Sm       <- lowess(pred.val, DATA.boot$y, iter = 0)
      #Sm <- loess(DATA.boot$y~pred.val, degree=2)
      #Sm <- data.frame(Sm$x, Sm$fitted)
      #Sm <- Sm[order(Sm$pred.val),]
      cal.x.orig.small[,i] <- seq(min(pred.val), max(pred.val), length = 100)
      cal.y.orig.small[,i] <- approx(Sm, xout = cal.x.orig.small[,i], ties = mean)$y
    }
    #}
  }
}


### Instability plot of individual predictions
x.large <- pred.fit.large[order(pred.fit.large)]
x.min   <- pred.fit.min[order(pred.fit.min)]
x.small <- pred.fit.small[order(pred.fit.small)]
y1.large <- apply(pred.large, 1, function(x) quantile(x, probs = 0.025, na.rm = T))[order(pred.fit.large)]
y1.min   <- apply(pred.min,   1, function(x) quantile(x, probs = 0.025, na.rm = T))[order(pred.fit.min)]
y1.small <- apply(pred.small, 1, function(x) quantile(x, probs = 0.025, na.rm = T))[order(pred.fit.small)]
y2.large <- apply(pred.large, 1, function(x) quantile(x, probs = 0.975, na.rm = T))[order(pred.fit.large)]
y2.min   <- apply(pred.min,   1, function(x) quantile(x, probs = 0.975, na.rm = T))[order(pred.fit.min)]
y2.small <- apply(pred.small, 1, function(x) quantile(x, probs = 0.975, na.rm = T))[order(pred.fit.small)]

OUT2 <- data.frame(x = c(x.large,  x.large,  x.min,  x.min,  x.small,  x.small),
                   y = c(y1.large, y2.large, y1.min, y2.min, y1.small, y2.small))

xx1.large <- lowess(y1.large~x.large, delta = 0.3)
xx1.min   <- lowess(y1.min~x.min,     delta = 0.3)
xx1.small <- lowess(y1.small~x.small, delta = 0.3)
xx2.large <- lowess(y2.large~x.large, delta = 0.3)
xx2.min   <- lowess(y2.min~x.min,     delta = 0.3)
xx2.small <- lowess(y2.small~x.small, delta = 0.3)

OUT3 <- data.frame(x = c(xx1.large$x, xx1.min$x, xx1.small$x,
                         xx2.large$x, xx2.min$x, xx2.small$x), 
                   y = c(xx1.large$y, xx1.min$y, xx1.small$y,
                         xx2.large$y, xx2.min$y, xx2.small$y))

OUT3$N.DEV <- c(rep("large",   length(xx1.large$x)), 
                rep("minimum", length(xx1.min$x)), 
                rep("small",   length(xx1.small$x)),
                rep("large",   length(xx2.large$x)),
                rep("minimum", length(xx2.min$x)), 
                rep("small",   length(xx2.small$x))) 

OUT3$N.DEV <- factor(OUT3$N.DEV, levels = c("large", "minimum", "small"))
OUT3$limit <- c(rep("lower", length(xx1.large$x) + length(xx1.min$x) + length(xx1.small$x)), 
                rep("upper", length(xx2.large$x) + length(xx2.min$x) + length(xx2.small$x)))
OUT3$limit <- factor(OUT3$limit, levels = c("lower", "upper"))

pred.large.long <- reshape2::melt(pred.large)
pred.min.long   <- reshape2::melt(pred.min)
pred.small.long <- reshape2::melt(pred.small)
names(pred.large.long)[2] <- "B"
names(pred.min.long)[2]   <- "B"
names(pred.small.long)[2] <- "B"
names(pred.large.long)[3] <- "boot.pred"
names(pred.min.long)[3]   <- "boot.pred"
names(pred.small.long)[3] <- "boot.pred"

pred.large.long <- add_column(pred.large.long, mod.pred = rep(pred.fit.large, NBOOT))
pred.min.long   <- add_column(pred.min.long,   mod.pred = rep(pred.fit.min,   NBOOT))
pred.small.long <- add_column(pred.small.long, mod.pred = rep(pred.fit.small, NBOOT))
pred.large.long <- add_column(pred.large.long, N.DEV = rep("large",   NBOOT * nrow(pred.large)))
pred.min.long   <- add_column(pred.min.long,   N.DEV = rep("minimum", NBOOT * nrow(pred.min)))
pred.small.long <- add_column(pred.small.long, N.DEV = rep("small",   NBOOT * nrow(pred.small)))

OUT <- bind_rows(pred.large.long, pred.min.long, pred.small.long)
OUT$N.DEV <- factor(OUT$N.DEV, levels = c("large", "minimum", "small"))
OUT  %>% ggplot(aes(x = mod.pred, y = boot.pred)) +
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
par(mfrow=c(1, 3))
plot(pred.fit.large, apply(abs(pred.large - pred.fit.large), 1, median, na.rm = T), pch = 20)
plot(pred.fit.min,   apply(abs(pred.min   - pred.fit.min),   1, median, na.rm = T), pch = 20)
plot(pred.fit.small, apply(abs(pred.small - pred.fit.small), 1, median, na.rm = T), pch = 20)

### MAPE_individual
MAPEi <- function(X, y){
  mean(abs(X-y))
}
mean(apply(pred.large, 2, MAPEi, y = pred.fit.large))
mean(apply(pred.min,   2, MAPEi, y = pred.fit.min))
mean(apply(pred.small, 2, MAPEi, y = pred.fit.small))

median(apply(pred.large, 2, MAPEi, y = pred.fit.large))
median(apply(pred.min,   2, MAPEi, y = pred.fit.min))
median(apply(pred.small, 2, MAPEi, y = pred.fit.small))


# plot calibration
cal.orig.large.xx <- reshape2::melt(cal.x.orig.large)
cal.orig.min.xx   <- reshape2::melt(cal.x.orig.min)
cal.orig.small.xx <- reshape2::melt(cal.x.orig.small)

cal.orig.large.yy <- reshape2::melt(cal.y.orig.large)
cal.orig.min.yy   <- reshape2::melt(cal.y.orig.min)
cal.orig.small.yy <- reshape2::melt(cal.y.orig.small)

cal.orig.large.xx$N <- rep("large",   nrow(cal.orig.large.xx))
cal.orig.min.xx$N   <- rep("minimum", nrow(cal.orig.min.xx))
cal.orig.small.xx$N <- rep("small",   nrow(cal.orig.small.xx))

cal <- rbind(cal.orig.large.xx,
             cal.orig.min.xx,
             cal.orig.small.xx)

colnames(cal) <- c("Var1", "Sim", "value", "N")
cal$N         <- factor(cal$N, labels = c("large", "minimum", "small"))

OUT.cal <- data.table(cal)
OUT.cal[, y:=c(cal.orig.large.yy$value,
               cal.orig.min.yy$value,
               cal.orig.small.yy$value)]


OUT.cal %>% ggplot(aes(x = value, y = y, group = Sim)) +
  geom_line(alpha = 0.5, colour = '#999999') + 
  facet_grid(~N) + 
  xlim(0, 1) + 
  ylim(0, 1) +
  xlab('predicted') +
  ylab('observed') +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  theme(axis.text = element_text(size = 6))

