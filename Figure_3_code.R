### Figure 3 code (it's not pretty and it's not quick)
### looking at the stability of the stability index

library(MASS)
library(data.table)
library(glmnet)
library(rms)
library(ragg)
library(progress)
library(tidyverse)

if(!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])

### define some functions
generate_data <- function(NN, n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors, cor0 = 0, cor1 = 0, beta.0 = 0, var0 = 4){
  
  n.predictors <- n.true.predictors + n.noise.predictors
  mu0 <- rep(0, n.predictors)
  
  # Specify correlation matrix
  Sigma0 <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
  Sigma0[1:n.true.predictors, 1:n.true.predictors] <- cor0
  Sigma0[(n.true.predictors+1):n.predictors, (n.true.predictors+1):n.predictors] <- cor1
  diag(Sigma0) <- c(var0, rep(1.0, n.noise.predictors))
  
  x <- mvrnorm(NN, mu0, Sigma0)
  
  #beta <- c(0.5, 0.3, 0.3, 0.25, 0.25, rep(0, n.noise.predictors))
  beta <- c(1, rep(0, n.noise.predictors))
  
  y <- runif(NN) < 1 / (1 + exp(-beta.0 - x %*% beta))
  
  DATA   <- data.frame(x)
  DATA$y <- y * 1
  DATA
}
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
auroc <- function(p, y) {
  n1 <- sum(!y)
  n2 <- sum(y)
  U  <- sum(rank(p)[!y]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}
val.perf2 <- function(data = DATA.val, fit){
  pred          <- predict(fit, s = fit$lambda.min, newx=as.matrix(subset(data, select=-c(y))), type='response')
  fit.val       <- glm(data$y~binomial()$linkfun(pred), family = 'binomial')
  cal.slope     <- as.numeric(coef(fit.val)[2])
  
  fit.cal       <- glm(data$y~offset(binomial()$linkfun(pred)), family = 'binomial')
  cal.intercept <- as.numeric(coef(fit.cal)[1])
  auc <- auroc(p = pred, y = DATA$y)
  
  return(list(c.val         = auc,
              cal.intercept = cal.intercept,
              cal.slope     = cal.slope))
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
# full model fit
full_model   <- function(data.d = DATA, B = B, data.v = DATA.val, n.cal.grid = 100, perf = F){
  
  ### fit maximum likelihood full model (no variable selection)
  #fit.full <- lrm.fit(x = subset(data.d, select = -c(y)), y = data.d$y)
  fit.full <- tryCatch(glm(y~., data = data.d, family='binomial'), 
                       error   = function(e) T, 
                       warning = function(w) T)
  
  if(perf){
    if(!is.logical(fit.full)){
      ### get apparent performance of full model (full.apparent)
      apparent <- as.numeric(val.perf(data = data.d, fit = fit.full))
      val      <- as.numeric(val.perf(data = data.v, fit = fit.full))
      pred.val <- predict(fit.full, newdata = data.v, type = 'response')
      Sm       <- lowess(pred.val, data.v$y, iter = 0)
      pp.full  <- seq(min(pred.val), max(pred.val), length = n.cal.grid)
      Sm.full  <- approx(Sm, xout = pp.full, ties = mean)$y
      fail     <- FALSE
    } else {
      apparent <- NA
      val      <- NA
      pp.full  <- rep(NA, n.cal.grid)
      Sm.full  <- rep(NA, n.cal.grid)
      fit.full <- fit
      fail     <- TRUE
    }
    return(list(apparent        = apparent, 
                val             = val, 
                fit.cal.x       = pp.full,
                fit.cal.y       = Sm.full,
                fit             = fit.full,
                fail            = fail))
  } else {
    return(list(fit = fit.full, fail = FALSE))
  }
}
# ridge (mixing=0), elastic net (mixing=0.5), lasso (mixing=1) model fit
EN_model     <- function(data.d = data, data.v = data.v, mixing = 1, n.cal.grid = 100, perf = F){
  P <- ncol(data.d) - 1
  fit.elasticnet <- cv.glmnet(x                   = as.matrix(data.d[,1:P]), 
                                     y            = data.d$y, 
                                     family       = 'binomial', 
                                     alpha        = mixing, 
                                     nfolds       = 5,
                                     type.measure = 'deviance', 
                                     parallel     = F)

    ### Get various performance measures
    if(perf) {
      apparent <- as.numeric(val.perf2(data = data.d, fit = fit.elasticnet))
      val      <- as.numeric(val.perf2(data = data.v, fit = fit.elasticnet))
      fail     <- FALSE
      pred  <- predict(fit.elasticnet, s = fit.elasticnet$lambda.min, newx=as.matrix(subset(data.v, select = -c(y))), type='response')
      Sm    <- lowess(pred, data.v$y, iter = 0)
      pp.EN <- seq(min(pred), max(pred), length = 100)
      Sm.EN <- approx(Sm, xout = pp.EN, ties = mean)$y
  
    return(list(apparent        = apparent, 
                val             = val,
                EN.cal.x        = pp.EN,
                EN.cal.y        = Sm.EN,
                fail            = fail,
                fit.shrink      = fit.elasticnet)) 
    } else {
      return(list(fit.shrink = fit.elasticnet))
    }
  
}
### MAPE_individual
MAPEi <- function(X, y){
  mean(abs(X - y), na.rm = T)
}

### 'POPULATION' (i.e., large sample model)
n.true.predictors  <- 1
n.noise.predictors <- 10
N.POP <- 1000000
NSIM  <- 200
NBOOT <- 100

set.seed(723141)

DATA.POP <- generate_data(NN = N.POP, n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors, beta.0 = 0)
fit.POP  <- lrm(y ~ ., data = DATA.POP, x = T, y = T)

### get the minimum sample size and randomly sample from the GUSTO data
min.ss <- pmsampsize::pmsampsize(type = 'b', rsquared = as.numeric(fit.POP$stats['R2(1e+06)']), parameters = 11, prevalence = 0.5)
N.DEV <- sort(c(50, 100, 500, 1000, 5000, min.ss$results_table[4,1]))

########################################
#### Logistic regression full model ####
########################################
MAPE.full <- matrix(ncol = length(N.DEV), nrow = NSIM)

pb <- progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = length(N.DEV), clear = FALSE, width = 60)

for(k in 1:length(N.DEV)){
  pb$tick()
  for(j in 1:NSIM){
    DATA      <- generate_data(NN = N.DEV[k], n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors)
    
    fit.full  <- tryCatch(full_model(data.d = DATA, data.v = DATA, B = 100, perf = F), 
                          error   = function(e) T, 
                          warning = function(w) T)
    
    if(!is.logical(fit.full$fit)){
      pred.full     <- matrix(ncol = NBOOT, nrow = nrow(DATA))
      pred.fit.full <- predict(fit.full$fit, newdata = DATA, type = 'response')
    
      for(i in 1:NBOOT){
        index <- sample(1:nrow(DATA), nrow(DATA), replace = T)
        DATA.boot <- DATA[index,]
        
        fit.boot <- tryCatch(full_model(data.d = DATA.boot, data.v = DATA, B = 100, perf = F), 
                             error   = function(e) T, 
                             warning = function(w) T)
        if(!is.logical(fit.boot$fit)){
          pred.full[,i] <- predict(fit.boot$fit, newdata = DATA, type = 'response')
        } else {
          pred.full[,i] <- rep(NA, nrow(DATA))
        }
      }
    } else {
      pred.full <- matrix(NA, ncol = NBOOT, nrow = nrow(DATA))
      pred.fit.full <- rep(NA, length = nrow(DATA))
    }
    MAPE.full[j,k] <- mean(apply(pred.full, 2, MAPEi, y = pred.fit.full), na.rm = T)
  }
}

################################################
#### Logistic regression with LASSO penalty ####
################################################
MAPE.lasso <- matrix(ncol = length(N.DEV), nrow = NSIM)

pb <- progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = length(N.DEV), clear = FALSE, width = 60)

for(k in 1:length(N.DEV)){
  pb$tick()
  for(j in 1:NSIM){
    DATA   <- generate_data(NN = N.DEV[k], n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors)
    
    fit.lasso <- tryCatch(EN_model(data.d = DATA, data.v = DATA, perf = F), 
                         error   = function(e) T, 
                         warning = function(w) T)
    
    if(!is.logical(fit.lasso)){
      pred.fit.lasso <- predict(fit.lasso$fit.shrink, 
                                s = fit.lasso$fit.shrink$lambda.min, 
                                newx = as.matrix(subset(DATA, select = -c(y))),
                                type = 'response')
      pred.lasso     <- matrix(ncol = NBOOT, nrow = nrow(DATA))
      
      for(i in 1:NBOOT){
        index <- sample(1:nrow(DATA), nrow(DATA), replace = T)
        DATA.boot <- DATA[index,]
        
        fit.boot <- tryCatch(EN_model(data.d = DATA.boot, data.v = DATA, perf = F), 
                             error   = function(e) T, 
                             warning = function(w) T)
        
        if(!is.logical(fit.boot)){
            pred.lasso[,i] <- predict(fit.boot$fit.shrink, 
                                      s = fit.boot$fit.shrink$lambda.min, 
                                      newx = as.matrix(subset(DATA, select = -c(y))), 
                                      type = 'response')
        } else { 
          pred.lasso[,i] <- rep(NA, nrow(DATA))
        }
      }
    } else { 
      pred.lasso     <- matrix(NA, ncol = NBOOT, nrow = nrow(DATA))
      pred.fit.lasso <- rep(NA, nrow(DATA))
    }
    MAPE.lasso[j,k] <- mean(apply(pred.lasso, 2, MAPEi, y = pred.fit.lasso), na.rm=T)
  }
}

#### plotting stuff
apply(MAPE.full,   2, summary)
apply(MAPE.lasso,  2, summary)

MAPE.full.long   <- reshape2::melt(MAPE.full)
MAPE.lasso.long  <- reshape2::melt(MAPE.lasso)

MAPE.full.long <- add_column(MAPE.full.long, N = c(rep(50, NSIM), 
                                                   rep(100, NSIM), 
                                                   rep(min.ss$results_table[4, 1], NSIM), 
                                                   rep(500, NSIM), 
                                                   rep(1000, NSIM), 
                                                   rep(5000, NSIM)))

MAPE.lasso.long <- add_column(MAPE.lasso.long, N = c(rep(50, NSIM), 
                                                     rep(100, NSIM), 
                                                     rep(min.ss$results_table[4, 1], NSIM), 
                                                     rep(500, NSIM), 
                                                     rep(1000, NSIM), 
                                                     rep(5000, NSIM)))

MAPE.full.long$N   <- factor(MAPE.full.long$N,   levels = N.DEV)
MAPE.lasso.long$N  <- factor(MAPE.lasso.long$N,  levels = N.DEV)

MAPE.full.long   <- add_column(MAPE.full.long,   approach = rep("full model", nrow(MAPE.full.long)))
MAPE.lasso.long  <- add_column(MAPE.lasso.long,  approach = rep("lasso",      nrow(MAPE.lasso.long)))

OUT.MAPE          <- bind_rows(MAPE.full.long, MAPE.lasso.long)
OUT.MAPE$approach <- factor(OUT.MAPE$approach, levels = c("full model", "lasso"))

alpha_val <- 1
p <- ggplot(OUT.MAPE, aes(x = N, y = value, fill = approach, colour = approach)) + 
  geom_jitter(alpha = alpha_val, aes(color = approach, shape = approach), 
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8),
              size = 1) +
  stat_summary(
    fun      = median, 
    geom     = "errorbar",
    aes(ymax = ..y.., ymin = ..y..), 
    position = position_dodge(width = 0.8), 
    width    = 0.25,
    colour   = 'black') + 
  xlab("Sample size") + 
  ylab("MAPE") + 
  theme_bw() +
  scale_colour_grey(start = 0.4, end = 0.6) + 
  theme(legend.position = "bottom") + 
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text = element_text(size = 8)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
p  




