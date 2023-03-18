### 5.2 Small sample size and mony noise variables
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
if(!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])


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
BE_model.p   <- function(fit.full = fit.full, data.d = data, B = B, data.v = DATA.val, n.cal.grid = 100){
  
  ### fit model using backwards elimination (then refit the model with only 
  ### those predictors retained)
  fit.fastbw <- fastbw(fit.full, rule = 'p', type = 'individual', sls = 0.05)
  
  # get the bootstrap-corrected slope for the BE model and apply and re-estimate the 
  # intercept, but need to check the model after backwards elimination is not an 
  # empty model
  if(is.null(fit.fastbw$factors.kept)) {
    apparent        <- rep(NA, 4)
    shrink.apparent <- rep(NA, 4)
    val             <- rep(NA, 4)
    shrink.val      <- rep(NA, 4)
    boot.opt        <- rep(NA, 4)
    return(list(apparent        = apparent, 
                val             = val, 
                opt.cor         = boot.opt,
                shrink.apparent = shrink.apparent,
                shrink.val      = shrink.val,
                fit.cal.x       = rep(NA, n.cal.grid),
                fit.cal.y       = rep(NA, n.cal.grid),
                shrink.cal.x    = rep(NA, n.cal.grid),
                shrink.cal.y    = rep(NA, n.cal.grid),
                fit             = NULL,
                fit.shrink      = NULL,
                shrinkage.f     = NULL))
  } else {
    BE.formula <- as.formula(paste("y~", paste(fit.fastbw$names.kept, collapse = "+"), sep = ''))
    fit.BE     <- rms::lrm(BE.formula, data = data.d, x = T, y = T) 
    val.BE     <- rms::validate(fit.full, bw = T, type = 'individual', rule = 'p', B = B, pr = F, estimates = F)
    xx         <- c(as.vector(val.BE[c("Dxy", "R2", "Intercept","Slope"), 'index.corrected']))
    boot.opt   <- c((1 + xx[1])/2, xx[-1])
    
    fit.shrink.BE <- recalibrate(shrinkage = boot.opt[4], fit = fit.BE, data = data.d)
    
    ## Assess (apparent and validation) performance of the BE model
    apparent        <- as.numeric(val.perf(data = data.d, fit = fit.BE))
    shrink.apparent <- as.numeric(val.perf(data = data.d, fit = fit.shrink.BE))
    val             <- as.numeric(val.perf(data = data.v, fit = fit.BE))
    shrink.val      <- as.numeric(val.perf(data = data.v, fit = fit.shrink.BE))
    
    pred.val <- predict(fit.BE, newdata = data.v, type = 'fitted')
    Sm       <- lowess(pred.val, data.v$y, iter = 0)
    pp.BE    <- seq(min(pred.val), max(pred.val), length = n.cal.grid)
    Sm.BE    <- approx(Sm, xout = pp.BE, ties = mean)$y
    
    pred.val  <- predict(fit.shrink.BE, newdata = data.v, type = 'fitted')
    Sm        <- lowess(pred.val, data.v$y, iter = 0)
    pp.shrink <- seq(min(pred.val), max(pred.val), length = n.cal.grid)
    Sm.shrink <- approx(Sm, xout = pp.shrink, ties = mean)$y
    
    return(list(apparent        = apparent, 
                val             = val, 
                opt.cor         = boot.opt,
                shrink.apparent = shrink.apparent,
                shrink.val      = shrink.val,
                fit.cal.x       = pp.BE,
                fit.cal.y       = Sm.BE,
                shrink.cal.x    = pp.shrink,
                shrink.cal.y    = Sm.shrink,
                fit             = fit.BE,
                fit.shrink      = fit.shrink.BE, 
                shrinkage.f     = boot.opt[4]))
  }
}
EN_model     <- function(data.d = data, data.v = data.v, mixing = 1, n.cal.grid = 100){
  P <- ncol(data.d)-1
  fit.glmnet.elasticnet <- cv.glmnet(x            = as.matrix(data.d[,1:P]), 
                                     y            = data.d$y, 
                                     family       = 'binomial', 
                                     alpha        = mixing, 
                                     nfolds       = 10,
                                     type.measure = 'deviance', 
                                     parallel     = F)
  
  fit.full <- lrm(y ~ ., data = data.d, x = T, y = T)
  ### Pull out lambda.min
  cv.elasticnet <- fit.glmnet.elasticnet$lambda.min
  
  ### Get betas for lambda.min
  optimal.beta.elasticnet <- as.numeric(predict(fit.glmnet.elasticnet, type = 'coefficients', s = "lambda.min"))
  
  ### create lrm object with coefficients from the glmnet
  fit.elasticnet <- lrm(y ~ ., data = data.d, init = optimal.beta.elasticnet, maxit = 1)
  
  if(all(coef(fit.elasticnet)[-1]==0)) { 
    apparent <- rep(NA, 4)
    val      <- rep(NA, 4)
    return(list(apparent        = apparent, 
                val             = val,
                EN.cal.x        = rep(NA, n.cal.grid),
                EN.cal.y        = rep(NA, n.cal.grid))) 
  } else {
    ### Get various performance measures (apparent)
    apparent <- as.numeric(val.perf(data = data.d, fit = fit.elasticnet))
    val      <- as.numeric(val.perf(data = data.v, fit = fit.elasticnet))
    
    pred  <- predict(fit.elasticnet, newdata = data.v, type = 'fitted')
    Sm    <- lowess(pred, data.v$y, iter = 0)
    pp.EN <- seq(min(pred), max(pred), length = 100)
    Sm.EN <- approx(Sm, xout = pp.EN, ties = mean)$y
    
    return(list(apparent        = apparent, 
                val             = val,
                EN.cal.x        = pp.EN,
                EN.cal.y        = Sm.EN,
                fit             = fit.full,
                fit.shrink      = fit.elasticnet)) 
  }
}

#setwd('/Users/collinsg/Library/CloudStorage/OneDrive-Nexus365/CSM/research/Data/Gusto (Frank)/')
setwd('/Users/gary/Library/CloudStorage/OneDrive-Nexus365/CSM/research/Data/Gusto (Frank)/')
DATA.ORIG <- read.csv('gusto.csv')
DATA.ORIG <- select(DATA.ORIG ,c("sex", "age", "hyp", "htn", "hrt", "pmi", "ste", "day30"))
names(DATA.ORIG)[8] <- 'y'

DATA.ORIG$sex2 <- rep(0, nrow(DATA.ORIG))
DATA.ORIG$sex2[DATA.ORIG$sex=='male'] <- 1

DATA.ORIG$pmi2 <- rep(0, nrow(DATA.ORIG))
DATA.ORIG$pmi2[DATA.ORIG$pmi == 'yes'] <- 1
DATA.ORIG <- select(DATA.ORIG, c("sex2", "age", "hyp", "htn", "hrt", "pmi2", "ste", "y"))

# simulate 20 noise predictors
for(i in 1:20){
  X <- rnorm(nrow(DATA.ORIG), 0, 1)
  DATA.ORIG[, ncol(DATA.ORIG) + 1] <- X
}
NBOOT <- 200

### get the minimum sample size and randomly sample from the GUSTO data
set.seed(723141)
min.ss <- pmsampsize::pmsampsize(type = 'b', rsquared = 0.08, parameters=7, prevalence = 2851/40830)
DATA   <- DATA.ORIG %>% sample_n(min.ss$results_table[4,1])

########################################
#### Logistic regression full model ####
########################################
fit.full       <- full_model(data.d = DATA, data.v = DATA, B = 100)
pred.fit.full  <- plogis(predict(fit.full$fit))
Sm.fit.full    <- loess(DATA$y~pred.fit.full, degree = 2)
Sm.fit.full    <- data.frame(Sm.fit.full$x, Sm.fit.full$fitted)
Sm.fit.full    <- Sm.fit.full[order(Sm.fit.full$pred.fit.full),]

pred.full       <- matrix(ncol = NBOOT, nrow = nrow(DATA))
cal.x.full      <- matrix(ncol = NBOOT, nrow = 100)
cal.y.full      <- matrix(ncol = NBOOT, nrow = 100)
cal.x.orig.full <- matrix(ncol = NBOOT, nrow = 100)
cal.y.orig.full <- matrix(ncol = NBOOT, nrow = 100)

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
    fit.boot <- full_model(data.d = DATA.boot, data.v = DATA, B = 100)
    if(any(sqrt(diag(vcov(fit.boot$fit)))>30)) {
      i <- i - 1  
    } else {
      pred.full[,i]  <- plogis(predict(fit.boot$fit, newdata = DATA))
      cal.x.full[,i] <- fit.boot$fit.cal.x # calibration of bootstrap model on boostrap data
      cal.y.full[,i] <- fit.boot$fit.cal.y
      
      # calibration of the original model on the bootstrap data
      pred.val <- predict(fit.full$fit, newdata = DATA.boot, type = 'fitted')
      Sm       <- lowess(pred.val, DATA.boot$y, iter = 0)
      cal.x.orig.full[,i] <- seq(min(pred.val), max(pred.val), length = 100)
      cal.y.orig.full[,i] <- approx(Sm, xout = cal.x.orig.full[,i], ties = mean)$y
    }
  }
}

###################################################
#### Logistic regression backwards elimination ####
###################################################
fit.full.BE <- lrm(y~., data = DATA, x = T, y = T)
fit.BE      <- BE_model.p(fit.full = fit.full.BE, data.d = DATA, data.v = DATA, B = 100)
pred.fit.BE <- plogis(predict(fit.BE$fit))
Sm.fit.BE   <- loess(DATA$y~pred.fit.BE, degree = 2)
Sm.fit.BE   <- data.frame(Sm.fit.BE$x, Sm.fit.BE$fitted)
Sm.fit.BE   <- Sm.fit.BE[order(Sm.fit.BE$pred.fit.BE), ]

pred.BE       <- matrix(ncol = NBOOT, nrow = nrow(DATA))
cal.x.BE      <- matrix(ncol = NBOOT, nrow = 100)
cal.y.BE      <- matrix(ncol = NBOOT, nrow = 100)
cal.x.orig.BE <- matrix(ncol = NBOOT, nrow = 100)
cal.y.orig.BE <- matrix(ncol = NBOOT, nrow = 100)

pb <- progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = NBOOT, clear = FALSE, width = 60)

for(i in 1:NBOOT){
  pb$tick()
  index <- sample(1:nrow(DATA), nrow(DATA), replace = T)
  DATA.boot <- DATA[index, ]
  if(sum(DATA.boot$y) == 0){ ### catch instances where bootstrapped data contains zero events
    i <- i-1
  } else {
    fit.full.boot <- lrm(y~., x = T, y = T, data = DATA.boot)
    fit.boot      <- BE_model.p(fit.full = fit.full.boot, data.d = DATA.boot, data.v = DATA, B = 100)
    
    if(is.null(fit.boot$fit)){
      i <- i - 1
      #if(any(sqrt(diag(vcov(fit.boot$fit))) > 30)) 
      # i <- i - 1  
    } else {
        pred.BE[,i]  <- plogis(predict(fit.boot$fit, newdata = DATA))
        cal.x.BE[,i] <- fit.boot$fit.cal.x # calibration of bootstrap model on bootstrap data
        cal.y.BE[,i] <- fit.boot$fit.cal.y
      
        # calibration of the original model on the bootstrap data
        pred.val <- predict(fit.BE$fit, newdata = DATA.boot, type = 'fitted')
        Sm       <- lowess(pred.val, DATA.boot$y, iter = 0)
        cal.x.orig.BE[,i] <- seq(min(pred.val), max(pred.val), length = 100)
        cal.y.orig.BE[,i] <- approx(Sm, xout = cal.x.orig.BE[,i], ties = mean)$y
    }
  }
}


################################################
#### Logistic regression with LASSO penalty ####
################################################
fit.lasso      <- EN_model(data.d = DATA, data.v = DATA)
pred.fit.lasso <- plogis(predict(fit.lasso$fit))
Sm.fit.lasso   <- loess(DATA$y~pred.fit.lasso, degree = 2)
Sm.fit.lasso   <- data.frame(Sm.fit.lasso$x, Sm.fit.lasso$fitted)
Sm.fit.lasso   <- Sm.fit.lasso[order(Sm.fit.lasso$pred.fit.lasso), ]

pred.lasso       <- matrix(ncol = NBOOT, nrow = nrow(DATA))
cal.x.lasso      <- matrix(ncol = NBOOT, nrow = 100)
cal.y.lasso      <- matrix(ncol = NBOOT, nrow = 100)
cal.x.orig.lasso <- matrix(ncol = NBOOT, nrow = 100)
cal.y.orig.lasso <- matrix(ncol = NBOOT, nrow = 100)

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
    fit.boot <- EN_model(data.d = DATA.boot, data.v = DATA)
    if(any(sqrt(diag(vcov(fit.boot$fit)))>30)) {
      i <- i - 1  
    } else {
      pred.lasso[,i]  <- plogis(predict(fit.boot$fit, newdata = DATA))
      cal.x.lasso[,i] <- fit.boot$EN.cal.x # calibration of bootstrap model on bootstrap data
      cal.y.lasso[,i] <- fit.boot$EN.cal.y
      
      # calibration of the original model on the bootstrap data
      pred.val <- predict(fit.lasso$fit, newdata = DATA.boot, type = 'fitted')
      Sm       <- lowess(pred.val, DATA.boot$y, iter = 0)
      cal.x.orig.lasso[,i] <- seq(min(pred.val), max(pred.val), length = 100)
      cal.y.orig.lasso[,i] <- approx(Sm, xout = cal.x.orig.lasso[,i], ties = mean)$y
    }
  }
}

### Instability plot of individual predictions
x.full  <- pred.fit.full[order(pred.fit.full)]
x.BE    <- pred.fit.BE[order(pred.fit.BE)]
x.lasso <- pred.fit.lasso[order(pred.fit.lasso)]

y1.full  <- apply(pred.full,  1, function(x) quantile(x, probs = 0.025, na.rm = T))[order(pred.fit.full)]
y2.full  <- apply(pred.full,  1, function(x) quantile(x, probs = 0.975, na.rm = T))[order(pred.fit.full)]
y1.BE    <- apply(pred.BE,    1, function(x) quantile(x, probs = 0.025, na.rm = T))[order(pred.fit.BE)]
y2.BE    <- apply(pred.BE,    1, function(x) quantile(x, probs = 0.975, na.rm = T))[order(pred.fit.BE)]
y1.lasso <- apply(pred.lasso, 1, function(x) quantile(x, probs = 0.025, na.rm = T))[order(pred.fit.lasso)]
y2.lasso <- apply(pred.lasso, 1, function(x) quantile(x, probs = 0.975, na.rm = T))[order(pred.fit.lasso)]

OUT2 <- data.frame(x = c(x.full,  x.full,  x.BE,  x.BE,  x.lasso,  x.lasso), 
                   y = c(y1.full, y2.full, y1.BE, y2.BE, y1.lasso, y2.lasso))

xx1.full  <- lowess(y1.full~x.full,   delta = 0.01, f = 0.1)
xx2.full  <- lowess(y2.full~x.full,   delta = 0.01, f = 0.1)
xx1.BE    <- lowess(y1.BE~x.BE,       delta = 0.01, f = 0.1)
xx2.BE    <- lowess(y2.BE~x.BE,       delta = 0.01, f = 0.1)
xx1.lasso <- lowess(y1.lasso~x.lasso, delta = 0.01, f = 0.1)
xx2.lasso <- lowess(y2.lasso~x.lasso, delta = 0.01, f = 0.1)

OUT3 <- data.frame(x = c(xx1.full$x, xx1.BE$x, xx1.lasso$x, xx2.full$x, xx2.BE$x, xx2.lasso$x), 
                   y = c(xx1.full$y, xx1.BE$y, xx1.lasso$y, xx2.full$y, xx2.BE$y, xx2.lasso$y))

OUT3$N.DEV <- c(rep("full",  length(xx1.full$x)), 
                rep("BE",    length(xx1.BE$x)),
                rep("lasso", length(xx1.lasso$x)),
                rep("full",  length(xx2.full$x)),
                rep("BE",    length(xx2.BE$x)),
                rep("lasso", length(xx2.lasso$x))) 

OUT3$N.DEV <- factor(OUT3$N.DEV, 
                     levels = c("full", "BE", "lasso"),
                     labels = c("Full model", "Backwards elimination", "LASSO"))
OUT3$limit <- c(rep("lower", length(xx1.full$x) + length(xx1.BE$x) + length(xx1.lasso$x)), 
                rep("upper", length(xx2.full$x) + length(xx2.BE$x) + length(xx2.lasso$x)))
OUT3$limit <- factor(OUT3$limit, levels = c("lower", "upper"))

pred.full.long  <- reshape2::melt(pred.full)
pred.BE.long    <- reshape2::melt(pred.BE)
pred.lasso.long <- reshape2::melt(pred.lasso)
names(pred.full.long)[2]  <- "B"
names(pred.BE.long)[2]    <- "B"
names(pred.lasso.long)[2] <- "B"
names(pred.full.long)[3]  <- "boot.pred"
names(pred.BE.long)[3]    <- "boot.pred"
names(pred.lasso.long)[3] <- "boot.pred"

pred.full.long  <- add_column(pred.full.long,  mod.pred = rep(pred.fit.full,  NBOOT))
pred.BE.long    <- add_column(pred.BE.long,    mod.pred = rep(pred.fit.BE,    NBOOT))
pred.lasso.long <- add_column(pred.lasso.long, mod.pred = rep(pred.fit.lasso, NBOOT))
pred.full.long  <- add_column(pred.full.long,  N.DEV = rep("full",  NBOOT * nrow(pred.full)))
pred.BE.long    <- add_column(pred.BE.long,    N.DEV = rep("BE",    NBOOT * nrow(pred.BE)))
pred.lasso.long <- add_column(pred.lasso.long, N.DEV = rep("lasso", NBOOT * nrow(pred.lasso)))

OUT       <- bind_rows(pred.full.long, pred.BE.long, pred.lasso.long)
OUT$N.DEV <- factor(OUT$N.DEV, 
                    levels = c("full", "BE", "lasso"),
                    labels = c("Full model", "Backwards elimination", "LASSO"))

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
par(mfrow = c(1, 3))
plot(pred.fit.full,  apply(abs(pred.full  - pred.fit.full),  1, median, na.rm = T), pch = 20, xlab='Estimated risk from developed model', ylab = 'Instability index')
plot(pred.fit.BE,    apply(abs(pred.BE    - pred.fit.BE),    1, median, na.rm = T), pch = 20, xlab='Estimated risk from developed model', ylab = 'Instability index')
plot(pred.fit.lasso, apply(abs(pred.lasso - pred.fit.lasso), 1, median, na.rm = T), pch = 20, xlab='Estimated risk from developed model', ylab = 'Instability index')

OUT.ii <- data.frame(pred  = c(pred.fit.full, pred.fit.BE, pred.fit.lasso),
                     index = c(apply(abs(pred.full  - pred.fit.full),  1, median, na.rm = T),
                               apply(abs(pred.BE    - pred.fit.BE),    1, median, na.rm = T),
                               apply(abs(pred.lasso - pred.fit.lasso), 1, median, na.rm = T)))

OUT.ii$label <- c(rep("full",  length(pred.fit.full)), 
                  rep("BE",    length(pred.fit.BE)), 
                  rep("lasso", length(pred.fit.lasso)))
OUT.ii$label <- factor(OUT.ii$label, 
                       levels = c("full", "BE", "lasso"),
                       labels = c("Full model", "Backwards elimination", "LASSO"))

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
mean(apply(pred.full,  2, MAPEi, y = pred.fit.full))
mean(apply(pred.BE,    2, MAPEi, y = pred.fit.BE))
mean(apply(pred.lasso, 2, MAPEi, y = pred.fit.lasso))

# plot calibration
cal.orig.full.xx  <- reshape2::melt(cal.x.orig.full)
cal.orig.BE.xx    <- reshape2::melt(cal.x.orig.BE)
cal.orig.lasso.xx <- reshape2::melt(cal.x.orig.lasso)

cal.orig.full.yy  <- reshape2::melt(cal.y.orig.full)
cal.orig.BE.yy    <- reshape2::melt(cal.y.orig.BE)
cal.orig.lasso.yy <- reshape2::melt(cal.y.orig.lasso)

cal.orig.full.xx$N  <- rep("full",  nrow(cal.orig.full.xx))
cal.orig.BE.xx$N    <- rep("BE",    nrow(cal.orig.BE.xx))
cal.orig.lasso.xx$N <- rep("lasso", nrow(cal.orig.lasso.xx))

cal <- rbind(cal.orig.full.xx,
             cal.orig.BE.xx,
             cal.orig.lasso.xx)

colnames(cal) <- c("Var1", "Sim", "value", "N")
cal$N         <- factor(cal$N, 
                        levels = c("full", "BE", "lasso"),
                        labels = c("Full model", "Backwards elimination", "LASSO"))

OUT.cal <- data.table(cal)
OUT.cal[, y:=c(cal.orig.full.yy$value,
               cal.orig.BE.yy$value,
               cal.orig.lasso.yy$value)]


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
