### Figure 2 code (it's not pretty and it's not quick)

### Looking at model stability in terms of variability in predictions 
### as a function of sample size.
### Backwards Elimination models with and without shrinkage
### compared against elastic net

rm(list = ls())
library(MASS)
#library(doMC)
library(data.table)
library(glmnet)
#library(ggh4x)
#library(ranger)
library(rms)
library(ragg)
library(progress)
#library(doParallel)
#library(facetscales)
library(tidyverse)
library(dcurves)
#library(apricom)
#library(patchwork)
if(!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])

### define some functions
# function to generate the data
generate_data_1_predictor_for_paper <- function(NN, n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors, cor0 = 0, cor1 = 0, beta.0 = 0, var0 = 4){
  
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
generate_data <- function(NN, n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors, cor0 = 0.1, cor1 = 0.05, beta.0 = 0){
  
  n.predictors <- n.true.predictors + n.noise.predictors
  mu0 <- rep(0, n.predictors)
  
  # Specify correlation matrix
  Sigma0 <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
  Sigma0[1:n.true.predictors, 1:n.true.predictors] <- cor0
  Sigma0[(n.true.predictors+1):n.predictors, (n.true.predictors+1):n.predictors] <- cor1
  diag(Sigma0) <- 1.0
  
  x <- mvrnorm(NN, mu0, Sigma0)
  
  beta <- c(0.5, 0.3, 0.3, 0.25, 0.25, rep(0, n.noise.predictors))
  
  y <- runif(NN) < 1 / (1 + exp(-beta.0 - x %*% beta))
  
  DATA   <- data.frame(x)
  DATA$y <- y * 1
  DATA
}
val.perf      <- function(data = DATA.val, fit = fit.full){
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
  
  ### get apparent performance of full model (full.apparent)
  apparent <- as.numeric(val.perf(data = data.d, fit = fit.full))
  val      <- as.numeric(val.perf(data = data.v, fit = fit.full))
  
  # bootstrap internal validation of the full model to get the optimism-corrected bootstrap
  # calibration slope
  val.full  <- rms::validate(fit.full, bw = F, B = B, pr = F, estimates = F)
  xx        <- c(as.vector(val.full[c("Dxy", "R2", "Intercept","Slope"), 'index.corrected']))
  boot.opt  <- c((1 +xx[1])/2, xx[-1])
  
  ## apply the bootstrap corrected slope shrinkage factor (and re-estimate the intercept)
  fit.shrunk <- recalibrate(shrinkage = boot.opt[4], fit = fit.full, data = data.d)
  shrunk.apparent <- as.numeric(val.perf(data = data.d, fit = fit.shrunk))
  shrunk.val      <- as.numeric(val.perf(data = data.v, fit = fit.shrunk))
  
  pred.val <- predict(fit.full, newdata = data.v, type = 'fitted')
  Sm       <- lowess(pred.val, data.v$y, iter = 0)
  pp.full  <- seq(min(pred.val), max(pred.val), length = n.cal.grid)
  Sm.full  <- approx(Sm, xout = pp.full, ties = mean)$y
  
  pred.shrunk <- predict(fit.shrunk, newdata = data.v, type = 'fitted')
  Sm          <- lowess(pred.shrunk, data.v$y, iter = 0)
  pp.shrunk   <- seq(min(pred.shrunk), max(pred.shrunk), length = n.cal.grid)
  Sm.shrunk   <- approx(Sm, xout = pp.shrunk, ties = mean)$y
  
  return(list(apparent        = apparent, 
              val             = val, 
              opt.cor         = boot.opt,
              shrunk.apparent = shrunk.apparent,
              shrunk.val      = shrunk.val,
              #fit.cal.x       = Sm.full$x,
              #fit.cal.y       = Sm.full$y,
              #shrunk.cal.x    = Sm.shrunk$x,
              #shrunk.cal.y    = Sm.shrunk$y,
              fit.cal.x       = pp.full,
              fit.cal.y       = Sm.full,
              shrunk.cal.x    = pp.shrunk,
              shrunk.cal.y    = Sm.shrunk,
              fit             = fit.full,
              fit.shrunk      = fit.shrunk, 
              shrinkage.f     = boot.opt[4]))
}
EN_model     <- function(data.d = data, data.v = data.v, mixing = 1, n.cal.grid = 100){
  P <- ncol(data.d)-1
  fit.glmnet.elasticnet <- cv.glmnet(x            = as.matrix(data.d[,1:P]), 
                                     y            = data.d$y, 
                                     family       = 'binomial', 
                                     alpha        = mixing, 
                                     nfolds       = 5,
                                     type.measure = 'deviance', 
                                     parallel     = F)
  
  fit.full <- lrm(y ~ ., data = data.d, x = T, y = T)
  ### Pull out lambda.min
  cv.elasticnet <- fit.glmnet.elasticnet$lambda.min
  
  ### Get betas for lambda.min
  optimal.fullta.elasticnet <- as.numeric(predict(fit.glmnet.elasticnet, type = 'coefficients', s = "lambda.min"))
  
  ### create lrm object with coefficients from the glmnet
  fit.elasticnet <- lrm(y ~ ., data = data.d, init = optimal.fullta.elasticnet, maxit = 1)
  
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
                fit.shrunk      = fit.elasticnet)) 
  }
}

NBOOT              <- 200
n.true.predictors  <- 1
n.noise.predictors <- 10
N.VAL              <- 5000
N.POP              <- 1000000

# fit model to a 'super' population (large sample model)
DATA.POP <- generate_data_1_predictor_for_paper(NN = N.POP, n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors)
fit.pop  <- lrm(y ~ ., data = DATA.POP, x = T, y = T)
phi      <- as.numeric(prop.table(table(DATA.POP$y))[2])
R2_CS    <- as.numeric(fit.pop$stats['R2']) * (1 - (phi^phi * (1 - phi)^(1 - phi))^2)

P <- ncol(DATA.POP) - 1
mixing <- 1

fit.glmnet.elasticnet <- cv.glmnet(x            = as.matrix(DATA.POP[,1:P]), 
                                   y            = DATA.POP$y, 
                                   family       = 'binomial', 
                                   alpha        = mixing, 
                                   nfolds       = 5,
                                   type.measure = 'deviance', 
                                   parallel     = F)
cv.elasticnet <- fit.glmnet.elasticnet$lambda.min
optimal.fullta.elasticnet <- as.numeric(predict(fit.glmnet.elasticnet, type = 'coefficients', s = "lambda.min"))
fit.EN.pop <- lrm(y ~ ., data = DATA.POP, init = optimal.fullta.elasticnet, maxit = 1)

# calculate required sample for model development
xx       <- pmsampsize::pmsampsize(type       = 'b', 
                                   parameters = n.true.predictors+n.noise.predictors,
                                   prevalence = as.numeric(prop.table(table(DATA.POP$y))[2]),
                                   rsquared   = R2_CS,
                                   shrinkage  = 0.9)

N.Riley <- xx$sample_size

DATA.val              <- generate_data_1_predictor_for_paper(NN = N.VAL, n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors)
N.SIM                 <- 500
n.cal.grid            <- 100
N                     <- sort(c(100, 200, 500, 1000, 5000, N.Riley))
p.full                <- array(0, dim = c(nrow(DATA.val), N.SIM, length(N)))
p.shrunk              <- array(0, dim = c(nrow(DATA.val), N.SIM, length(N)))
p.EN                  <- array(0, dim = c(nrow(DATA.val), N.SIM, length(N)))
cal.shrunk            <- matrix(ncol = length(N), nrow = N.SIM)
cal.EN                <- matrix(ncol = length(N), nrow = N.SIM)
full.val.cal.x        <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
full.val.cal.y        <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
full.val.shrunk.cal.x <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
full.val.shrunk.cal.y <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
EN.val.cal.x          <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
EN.val.cal.y          <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
dca.full              <- array(NA, dim = c(100, N.SIM, length(N)))
dca.shrunk            <- array(NA, dim = c(100, N.SIM, length(N)))
dca.EN                <- array(NA, dim = c(100, N.SIM, length(N)))

pb <- progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = length(N) * N.SIM, clear = FALSE, width = 60)

for(j in 1:length(N)){
  cat("\n loop = ", j)
  for(i in 1:N.SIM){
    pb$tick()
    DATA <- generate_data_1_predictor_for_paper(NN = N[j], n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors)
    
    # full model + shrinkage
    fit.full <- full_model(data.d = DATA, B = NBOOT, data.v = DATA.val)
    
    ### get the loess calibration on the validation data (full model)
    full.val.cal.x[,i,j] <- fit.full$fit.cal.x
    full.val.cal.y[,i,j] <- fit.full$fit.cal.y
    
    ### get the loess calibration on the validation data (shrunken model)
    full.val.shrunk.cal.x[,i,j] <- fit.full$shrunk.cal.x
    full.val.shrunk.cal.y[,i,j] <- fit.full$shrunk.cal.y
    
    if(length(class(fit.full$fit))>1){
      p.full[,i,j]    <- predict(fit.full$fit, newdata = DATA.val, type = 'fitted')
      p.shrunk[,i,j]  <- predict(fit.full$fit.shrunk, newdata = DATA.val, type = 'fitted')
      cal.shrunk[i,j] <- as.numeric(val.perf(data = DATA.val, fit = fit.full$fit.shrunk))[4]
      DAT.DCA         <- data.frame(pred.full   = predict(fit.full$fit,        newdata = DATA.val, type = 'fitted'),
                                    pred.shrunk = predict(fit.full$fit.shrunk, newdata = DATA.val, type = 'fitted'),
                                    y           = DATA.val$y)
      fit.dca.full     <- dca(y~pred.full, data = DAT.DCA)
      fit.dca.shrunk   <- dca(y~pred.shrunk, data = DAT.DCA)
      dca.full[,i,j]   <- fit.dca.full$dca$net_benefit[fit.dca.full$dca$label == 'pred.full']
      dca.shrunk[,i,j] <- fit.dca.shrunk$dca$net_benefit[fit.dca.shrunk$dca$label == 'pred.shrunk']
    } else {
      p.full[,i, j]        <- rep(NA, n.cal.grid)
      p.shrunk[,i, j] <- rep(NA, n.cal.grid)
      cal.shrunk[i, j]     <- NA
      dca.full[,i, j]      <- rep(NA, 100)
      dca.shrunk[,i, j]    <- rep(NA, 100)
    }
  
    # Elastic net
    fit.EN   <- EN_model(data.d = DATA, data.v = DATA.val)
    
    ### get the loess calibration on the validation data (elastic net model)
    EN.val.cal.x[,i,j] <- fit.EN$EN.cal.x
    EN.val.cal.y[,i,j] <- fit.EN$EN.cal.y
    
    if(length(class(fit.EN$fit))>1){
      p.EN[,i,j]   <- predict(fit.EN$fit, newdata = DATA.val, type = 'fitted')
      cal.EN[i,j]  <- as.numeric(val.perf(data = DATA.val, fit = fit.EN$fit.shrunk))[4]
      DAT.DCA      <- data.frame(pred.EN = predict(fit.EN$fit.shrunk, newdata = DATA.val, type = 'fitted'),y = DATA.val$y)
      fit.dca.EN   <- dca(y~pred.EN,     data = DAT.DCA)
      dca.EN[,i,j] <- fit.dca.EN$dca$net_benefit[fit.dca.EN$dca$label == 'pred.EN']
    } else {
      p.EN[,i, j]  <- rep(NA, n.cal.grid)
      cal.EN[i, j] <- NA
      dca.EN[,i,j] <- rep(NA, 100)
    }
  }
}

### get predictions of the 'population' model on the validation data
p.POP    <- predict(fit.pop, newdata = DATA.val, type = 'fitted')
p.POP.EN <- predict(fit.EN.pop, newdata = DATA.val, type = 'fitted')

# find the randomley chosen individual with large sample predictions of 0.1:0.9[0.1]
which.individual <- c(which(as.numeric(round(p.POP, 2)) == 0.1)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.1)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.2)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.2)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.3)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.3)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.4)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.4)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.5)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.5)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.6)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.6)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.7)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.7)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.8)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.8)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.9)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.9)))[1]])

which.individual.EN <- c(which(as.numeric(round(p.POP.EN, 2)) == 0.1)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.1)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.2)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.2)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.3)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.3)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.4)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.4)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.5)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.5)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.6)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.6)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.7)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.7)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.8)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.8)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.9)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.9)))[1]])

OUT.orig.1 <- gather(as_tibble(cbind(p.full[which.individual[1],,1], 
                                     p.full[which.individual[1],,2],
                                     p.full[which.individual[1],,3], 
                                     p.full[which.individual[1],,4], 
                                     p.full[which.individual[1],,5], 
                                     p.full[which.individual[1],,6])))

OUT.orig.2 <- gather(as_tibble(cbind(p.full[which.individual[2],,1],
                                     p.full[which.individual[2],,2],
                                     p.full[which.individual[2],,3], 
                                     p.full[which.individual[2],,4], 
                                     p.full[which.individual[2],,5], 
                                     p.full[which.individual[2],,6])))

OUT.orig.3 <- gather(as_tibble(cbind(p.full[which.individual[3],,1], 
                                     p.full[which.individual[3],,2],
                                     p.full[which.individual[3],,3], 
                                     p.full[which.individual[3],,4], 
                                     p.full[which.individual[3],,5], 
                                     p.full[which.individual[3],,6])))

OUT.orig.4 <- gather(as_tibble(cbind(p.full[which.individual[4],,1], 
                                     p.full[which.individual[4],,2],
                                     p.full[which.individual[4],,3], 
                                     p.full[which.individual[4],,4], 
                                     p.full[which.individual[4],,5], 
                                     p.full[which.individual[4],,6])))

OUT.orig.5 <- gather(as_tibble(cbind(p.full[which.individual[5],,1], 
                                     p.full[which.individual[5],,2],
                                     p.full[which.individual[5],,3], 
                                     p.full[which.individual[5],,4], 
                                     p.full[which.individual[5],,5], 
                                     p.full[which.individual[5],,6])))

OUT.orig.6 <- gather(as_tibble(cbind(p.full[which.individual[6],,1], 
                                     p.full[which.individual[6],,2],
                                     p.full[which.individual[6],,3], 
                                     p.full[which.individual[6],,4], 
                                     p.full[which.individual[6],,5], 
                                     p.full[which.individual[6],,6])))

OUT.orig.7 <- gather(as_tibble(cbind(p.full[which.individual[7],,1], 
                                     p.full[which.individual[7],,2],
                                     p.full[which.individual[7],,3], 
                                     p.full[which.individual[7],,4], 
                                     p.full[which.individual[7],,5], 
                                     p.full[which.individual[7],,6])))

OUT.orig.8 <- gather(as_tibble(cbind(p.full[which.individual[8],,1], 
                                     p.full[which.individual[8],,2],
                                     p.full[which.individual[8],,3], 
                                     p.full[which.individual[8],,4], 
                                     p.full[which.individual[8],,5], 
                                     p.full[which.individual[8],,6])))

OUT.orig.9 <- gather(as_tibble(cbind(p.full[which.individual[9],,1], 
                                     p.full[which.individual[9],,2],
                                     p.full[which.individual[9],,3], 
                                     p.full[which.individual[9],,4], 
                                     p.full[which.individual[9],,5], 
                                     p.full[which.individual[9],,6])))

OUT.shrunk.1 <- gather(as_tibble(cbind(p.shrunk[which.individual[1],,1], 
                                       p.shrunk[which.individual[1],,2],
                                       p.shrunk[which.individual[1],,3], 
                                       p.shrunk[which.individual[1],,4], 
                                       p.shrunk[which.individual[1],,5], 
                                       p.shrunk[which.individual[1],,6])))

OUT.shrunk.2 <- gather(as_tibble(cbind(p.shrunk[which.individual[2],,1],
                                       p.shrunk[which.individual[2],,2],
                                       p.shrunk[which.individual[2],,3], 
                                       p.shrunk[which.individual[2],,4], 
                                       p.shrunk[which.individual[2],,5], 
                                       p.shrunk[which.individual[2],,6])))

OUT.shrunk.3 <- gather(as_tibble(cbind(p.shrunk[which.individual[3],,1], 
                                       p.shrunk[which.individual[3],,2],
                                       p.shrunk[which.individual[3],,3], 
                                       p.shrunk[which.individual[3],,4], 
                                       p.shrunk[which.individual[3],,5], 
                                       p.shrunk[which.individual[3],,6])))

OUT.shrunk.4 <- gather(as_tibble(cbind(p.shrunk[which.individual[4],,1], 
                                       p.shrunk[which.individual[4],,2],
                                       p.shrunk[which.individual[4],,3], 
                                       p.shrunk[which.individual[4],,4], 
                                       p.shrunk[which.individual[4],,5], 
                                       p.shrunk[which.individual[4],,6])))

OUT.shrunk.5 <- gather(as_tibble(cbind(p.shrunk[which.individual[5],,1], 
                                       p.shrunk[which.individual[5],,2],
                                       p.shrunk[which.individual[5],,3], 
                                       p.shrunk[which.individual[5],,4], 
                                       p.shrunk[which.individual[5],,5], 
                                       p.shrunk[which.individual[5],,6])))

OUT.shrunk.6 <- gather(as_tibble(cbind(p.shrunk[which.individual[6],,1], 
                                       p.shrunk[which.individual[6],,2],
                                       p.shrunk[which.individual[6],,3], 
                                       p.shrunk[which.individual[6],,4], 
                                       p.shrunk[which.individual[6],,5], 
                                       p.shrunk[which.individual[6],,6])))

OUT.shrunk.7 <- gather(as_tibble(cbind(p.shrunk[which.individual[7],,1], 
                                       p.shrunk[which.individual[7],,2],
                                       p.shrunk[which.individual[7],,3], 
                                       p.shrunk[which.individual[7],,4], 
                                       p.shrunk[which.individual[7],,5], 
                                       p.shrunk[which.individual[7],,6])))

OUT.shrunk.8 <- gather(as_tibble(cbind(p.shrunk[which.individual[8],,1], 
                                       p.shrunk[which.individual[8],,2],
                                       p.shrunk[which.individual[8],,3], 
                                       p.shrunk[which.individual[8],,4], 
                                       p.shrunk[which.individual[8],,5], 
                                       p.shrunk[which.individual[8],,6])))

OUT.shrunk.9 <- gather(as_tibble(cbind(p.shrunk[which.individual[9],,1], 
                                       p.shrunk[which.individual[9],,2],
                                       p.shrunk[which.individual[9],,3], 
                                       p.shrunk[which.individual[9],,4], 
                                       p.shrunk[which.individual[9],,5], 
                                       p.shrunk[which.individual[9],,6])))

OUT.EN.1 <- gather(as_tibble(cbind(p.EN[which.individual.EN[1],,1], 
                                   p.EN[which.individual.EN[1],,2],
                                   p.EN[which.individual.EN[1],,3], 
                                   p.EN[which.individual.EN[1],,4], 
                                   p.EN[which.individual.EN[1],,5], 
                                   p.EN[which.individual.EN[1],,6])))

OUT.EN.2 <- gather(as_tibble(cbind(p.EN[which.individual.EN[2],,1],
                                   p.EN[which.individual.EN[2],,2],
                                   p.EN[which.individual.EN[2],,3], 
                                   p.EN[which.individual.EN[2],,4], 
                                   p.EN[which.individual.EN[2],,5], 
                                   p.EN[which.individual.EN[2],,6])))

OUT.EN.3 <- gather(as_tibble(cbind(p.EN[which.individual.EN[3],,1], 
                                   p.EN[which.individual.EN[3],,2],
                                   p.EN[which.individual.EN[3],,3], 
                                   p.EN[which.individual.EN[3],,4], 
                                   p.EN[which.individual.EN[3],,5], 
                                   p.EN[which.individual.EN[3],,6])))

OUT.EN.4 <- gather(as_tibble(cbind(p.EN[which.individual.EN[4],,1], 
                                   p.EN[which.individual.EN[4],,2],
                                   p.EN[which.individual.EN[4],,3], 
                                   p.EN[which.individual.EN[4],,4], 
                                   p.EN[which.individual.EN[4],,5], 
                                   p.EN[which.individual.EN[4],,6])))

OUT.EN.5 <- gather(as_tibble(cbind(p.EN[which.individual.EN[5],,1], 
                                   p.EN[which.individual.EN[5],,2],
                                   p.EN[which.individual.EN[5],,3], 
                                   p.EN[which.individual.EN[5],,4], 
                                   p.EN[which.individual.EN[5],,5], 
                                   p.EN[which.individual.EN[5],,6])))

OUT.EN.6 <- gather(as_tibble(cbind(p.EN[which.individual.EN[6],,1], 
                                   p.EN[which.individual.EN[6],,2],
                                   p.EN[which.individual.EN[6],,3], 
                                   p.EN[which.individual.EN[6],,4], 
                                   p.EN[which.individual.EN[6],,5], 
                                   p.EN[which.individual.EN[6],,6])))

OUT.EN.7 <- gather(as_tibble(cbind(p.EN[which.individual.EN[7],,1], 
                                   p.EN[which.individual.EN[7],,2],
                                   p.EN[which.individual.EN[7],,3], 
                                   p.EN[which.individual.EN[7],,4], 
                                   p.EN[which.individual.EN[7],,5], 
                                   p.EN[which.individual.EN[7],,6])))

OUT.EN.8 <- gather(as_tibble(cbind(p.EN[which.individual.EN[8],,1], 
                                   p.EN[which.individual.EN[8],,2],
                                   p.EN[which.individual.EN[8],,3], 
                                   p.EN[which.individual.EN[8],,4], 
                                   p.EN[which.individual.EN[8],,5], 
                                   p.EN[which.individual.EN[8],,6])))

OUT.EN.9 <- gather(as_tibble(cbind(p.EN[which.individual.EN[9],,1], 
                                   p.EN[which.individual.EN[9],,2],
                                   p.EN[which.individual.EN[9],,3], 
                                   p.EN[which.individual.EN[9],,4], 
                                   p.EN[which.individual.EN[9],,5], 
                                   p.EN[which.individual.EN[9],,6])))


OUT.orig <- bind_rows(OUT.orig.1, OUT.orig.2, OUT.orig.3,
                      OUT.orig.4, OUT.orig.5, OUT.orig.6,
                      OUT.orig.7, OUT.orig.8, OUT.orig.9)

OUT.shrunk <- bind_rows(OUT.shrunk.1, OUT.shrunk.2, OUT.shrunk.3,
                        OUT.shrunk.4, OUT.shrunk.5, OUT.shrunk.6,
                        OUT.shrunk.7, OUT.shrunk.8, OUT.shrunk.9)

OUT.EN <- bind_rows(OUT.EN.1, OUT.EN.2, OUT.EN.3,
                    OUT.EN.4, OUT.EN.5, OUT.EN.6,
                    OUT.EN.7, OUT.EN.8, OUT.EN.9)

OUT.orig   <- add_column(OUT.orig,   model = rep("full model", nrow(OUT.orig)))
OUT.shrunk <- add_column(OUT.shrunk, model = rep("shrunk LR",  nrow(OUT.shrunk)))
OUT.EN     <- add_column(OUT.EN,     model = rep("lasso",      nrow(OUT.EN)))

OUT     <- bind_rows(OUT.orig, OUT.shrunk, OUT.EN)
OUT$key <- factor(OUT$key, levels = paste("V", seq(1:length(N)), sep = ''), labels = as.character(N))
OUT     <- mutate(OUT, model = factor(model, levels = c("full model", "shrunk LR", "lasso")))
OUT     <- mutate(OUT, key = factor(key, levels = as.character(N)))
OUT     <- add_column(OUT, prob = c(rep("p=0.1", N.SIM * length(N)),
                                    rep("p=0.2", N.SIM * length(N)),
                                    rep("p=0.3", N.SIM * length(N)),
                                    rep("p=0.4", N.SIM * length(N)),
                                    rep("p=0.5", N.SIM * length(N)),
                                    rep("p=0.6", N.SIM * length(N)),
                                    rep("p=0.7", N.SIM * length(N)),
                                    rep("p=0.8", N.SIM * length(N)),
                                    rep("p=0.9", N.SIM * length(N)),
                                    rep("p=0.1", N.SIM * length(N)),
                                    rep("p=0.2", N.SIM * length(N)),
                                    rep("p=0.3", N.SIM * length(N)),
                                    rep("p=0.4", N.SIM * length(N)),
                                    rep("p=0.5", N.SIM * length(N)),
                                    rep("p=0.6", N.SIM * length(N)),
                                    rep("p=0.7", N.SIM * length(N)),
                                    rep("p=0.8", N.SIM * length(N)),
                                    rep("p=0.9", N.SIM * length(N)),
                                    rep("p=0.1", N.SIM * length(N)),
                                    rep("p=0.2", N.SIM * length(N)),
                                    rep("p=0.3", N.SIM * length(N)),
                                    rep("p=0.4", N.SIM * length(N)),
                                    rep("p=0.5", N.SIM * length(N)),
                                    rep("p=0.6", N.SIM * length(N)),
                                    rep("p=0.7", N.SIM * length(N)),
                                    rep("p=0.8", N.SIM * length(N)),
                                    rep("p=0.9", N.SIM * length(N))))
OUT <- mutate(OUT, prob = factor(prob))

OUT_mean <- OUT %>% filter(key == '5000') %>% group_by(prob) %>% summarize(mean_val = mean(value, na.rm = T))
OUT_mean$mean_val <- seq(0.1, 0.9, 0.1)

x1 <- OUT %>% filter(model == 'full model') %>% group_by(prob, key) %>% summarize(sd(value, na.rm = T))
x2 <- OUT %>% filter(model == 'shrunk LR') %>% group_by(prob, key) %>% summarize(sd(value, na.rm = T))
x3 <- OUT %>% filter(model == 'lasso') %>% group_by(prob, key) %>% summarize(sd(value, na.rm = T))
x1 <- as.data.frame(x1)
x2 <- as.data.frame(x2)
x3 <- as.data.frame(x3)
colnames(x1) <- c("prob", "key", "sd")
colnames(x2) <- c("prob", "key", "sd")
colnames(x3) <- c("prob", "key", "sd")
xx.sd <- cbind(x1, x2$sd, x3$sd)
colnames(xx.sd)<- c("prob", "N", "sd(full model)", "sd(shrunk)", "sd(lasso)")

### plot calibration
cal.full.xx <- reshape2::melt(full.val.cal.x)
cal.full.yy <- reshape2::melt(full.val.cal.y)

cal.shrunk.xx <- reshape2::melt(full.val.shrunk.cal.x)
cal.shrunk.yy <- reshape2::melt(full.val.shrunk.cal.y)

cal.EN.xx <- reshape2::melt(EN.val.cal.x)
cal.EN.yy <- reshape2::melt(EN.val.cal.y)

cal.full.xx$method   <- rep("full model", nrow(cal.full.xx))
cal.shrunk.xx$method <- rep("shrunk LR",   nrow(cal.shrunk.xx))
cal.EN.xx$method     <- rep("lasso", nrow(cal.EN.xx))

cal <- rbind(cal.full.xx,   
             cal.shrunk.xx,
             cal.EN.xx)

colnames(cal) <- c("Var1", "Sim", "N", "value", "method")
cal$method    <- factor(cal$method, levels = c('full model', 'shrunk LR', 'lasso'))
cal$N         <- factor(cal$N, levels = seq(1:length(N)), labels = paste("N=", N, sep=''))


OUT.cal <- data.table(cal)
OUT.cal[, y:=c(cal.full.yy$value,   
               cal.shrunk.yy$value,
               cal.EN.yy$value)]
  
p3 <- OUT.cal %>% dplyr::filter(method != 'shrunk LR') %>% ggplot(aes(x = value, y = y, group = Sim)) +
  geom_line(alpha = 0.05, colour = 'darkgrey') + 
  facet_grid(method~N) + 
  xlim(0, 1) + 
  ylim(0, 1) +
  xlab('estimated risk') +
  ylab('observed') +
  geom_abline(slope = 1, intercept = 0, colour = 'black') +
  theme_bw() +
  theme(axis.text = element_text(size = 6))
p3

