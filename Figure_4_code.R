### Figure 4 code (it's not pretty and it's not quick)
library(MASS)
library(rms)
library(ragg)
library(glmnet)
library(ggh4x)
library(tidyverse)
library(progress)
library(bde)
if(!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])

### define some functions
# function to generate the data
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
EN_model      <- function(data.d = data, data.v = data.v, mixing = 1, n.cal.grid = 100){
  # alpha = 1 -> lasso
  P <- ncol(data.d) - 1
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
                fit.shrunk      = fit.elasticnet)) 
  }
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

set.seed(723141)
min.ss <- pmsampsize::pmsampsize(type = 'b', cstatistic = 0.861, parameters = 11, prevalence = 0.5)

set.seed(24576)
n.true.predictors  <- 1
n.noise.predictors <- 10
N.SIM <- 1000
N.DEV <- sort(c(50, 100, 500, 1000, 5000, min.ss$results_table[4,1]))
N.VAL <- 100000
N.POP <- 1000000

### 'POPULATION' (i.e., large sample model)
DATA.POP <- generate_data(NN = N.POP, n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors, beta.0 = 0)
DATA.val <- generate_data(NN = N.VAL, n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors, beta.0 = 0)
fit.POP  <- lrm(y ~ ., data = DATA.POP, x = T, y = T)
p.POP    <- predict(fit.POP, newdata = DATA.val, type = 'fitted')
LP.POP   <- predict(fit.POP, newdata = DATA.val, type = 'lp')

# find a randomly chosen individual with large sample predictions of 0.1:0.9[0.1]
which.individual <- c(which(as.numeric(round(p.POP, 2)) == 0.1)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.1)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.2)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.2)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.3)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.3)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.4)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.4)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.5)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.5)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.6)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.6)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.7)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.7)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.8)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.8)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.9)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.9)))[1]])

### looking at mean predicted risk
p.EN     <- array(dim = c(nrow(DATA.val), N.SIM, length(N.DEV)))
X.coef   <- array(dim = c(N.SIM, n.true.predictors + n.noise.predictors + 1, length(N.DEV)))

pb <- progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = length(N.DEV) * N.SIM, clear = FALSE, width = 60)

for(j in 1:length(N.DEV)){
  for(i in 1:N.SIM){
    pb$tick()
    DATA   <- generate_data(NN = N.DEV[j], n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors)
    fit.EN <- EN_model(data.d = DATA, data.v = DATA.val)
    if(length(class(fit.EN$fit))>1){
      p.EN[,i,j]   <- predict(fit.EN$fit, newdata = DATA.val, type = 'fitted')
      X.coef[i,,j] <- as.numeric(coef(fit.EN$fit.shrunk))
    } else { # models that don't converge
      p.EN[,i,j]   <- rep(NA, nrow(DATA.val))
      X.coef[i,,j] <- rep(NA, (1+n.true.predictors+n.noise.predictors))
    }
  }
}

# Look at predictions for specific individuals with probabilities at 0.1:0.9[0.1]
OUT.1 <- gather(as_tibble(cbind(p.EN[which.individual[1],,1], 
                                p.EN[which.individual[2],,1],
                                p.EN[which.individual[3],,1], 
                                p.EN[which.individual[4],,1], 
                                p.EN[which.individual[5],,1], 
                                p.EN[which.individual[6],,1],   
                                p.EN[which.individual[7],,1],
                                p.EN[which.individual[8],,1],
                                p.EN[which.individual[9],,1])))

OUT.2 <- gather(as_tibble(cbind(p.EN[which.individual[1],,2], 
                                p.EN[which.individual[2],,2],
                                p.EN[which.individual[3],,2], 
                                p.EN[which.individual[4],,2], 
                                p.EN[which.individual[5],,2], 
                                p.EN[which.individual[6],,2],   
                                p.EN[which.individual[7],,2],
                                p.EN[which.individual[8],,2],
                                p.EN[which.individual[9],,2])))

OUT.3 <- gather(as_tibble(cbind(p.EN[which.individual[1],,3], 
                                p.EN[which.individual[2],,3],
                                p.EN[which.individual[3],,3], 
                                p.EN[which.individual[4],,3], 
                                p.EN[which.individual[5],,3], 
                                p.EN[which.individual[6],,3],   
                                p.EN[which.individual[7],,3],
                                p.EN[which.individual[8],,3],
                                p.EN[which.individual[9],,3])))

OUT.4 <- gather(as_tibble(cbind(p.EN[which.individual[1],,4], 
                                p.EN[which.individual[2],,4],
                                p.EN[which.individual[3],,4], 
                                p.EN[which.individual[4],,4], 
                                p.EN[which.individual[5],,4], 
                                p.EN[which.individual[6],,4],   
                                p.EN[which.individual[7],,4],
                                p.EN[which.individual[8],,4],
                                p.EN[which.individual[9],,4])))

OUT.5 <- gather(as_tibble(cbind(p.EN[which.individual[1],,5], 
                                p.EN[which.individual[2],,5],
                                p.EN[which.individual[3],,5], 
                                p.EN[which.individual[4],,5], 
                                p.EN[which.individual[5],,5], 
                                p.EN[which.individual[6],,5],   
                                p.EN[which.individual[7],,5],
                                p.EN[which.individual[8],,5],
                                p.EN[which.individual[9],,5])))

OUT.6 <- gather(as_tibble(cbind(p.EN[which.individual[1],,6], 
                                p.EN[which.individual[2],,6],
                                p.EN[which.individual[3],,6], 
                                p.EN[which.individual[4],,6], 
                                p.EN[which.individual[5],,6], 
                                p.EN[which.individual[6],,6],   
                                p.EN[which.individual[7],,6],
                                p.EN[which.individual[8],,6],
                                p.EN[which.individual[9],,6])))



OUT.1 <- add_column(OUT.1, N.DEV = rep(N.DEV[1], nrow(OUT.1)))
OUT.2 <- add_column(OUT.2, N.DEV = rep(N.DEV[2], nrow(OUT.2)))
OUT.3 <- add_column(OUT.3, N.DEV = rep(N.DEV[3], nrow(OUT.3)))
OUT.4 <- add_column(OUT.4, N.DEV = rep(N.DEV[4], nrow(OUT.4)))
OUT.5 <- add_column(OUT.5, N.DEV = rep(N.DEV[5], nrow(OUT.5)))
OUT.6 <- add_column(OUT.6, N.DEV = rep(N.DEV[6], nrow(OUT.6)))
OUT   <- bind_rows(OUT.1, OUT.2, OUT.3, OUT.4, OUT.5, OUT.6)

OUT <- add_column(OUT, prob = rep(c(rep("p=0.1", N.SIM),
                                    rep("p=0.2", N.SIM),
                                    rep("p=0.3", N.SIM),
                                    rep("p=0.4", N.SIM),
                                    rep("p=0.5", N.SIM),
                                    rep("p=0.6", N.SIM),
                                    rep("p=0.7", N.SIM),
                                    rep("p=0.8", N.SIM),
                                    rep("p=0.9", N.SIM)), length(N.DEV)))

OUT <- mutate(OUT, prob = factor(prob))
#OUT <- mutate(OUT, N.DEV = factor(N.DEV, levels = N.DEV))
OUT <- add_column(OUT, title = rep("Development data sample size", nrow(OUT)))

alpha_val <- 0.3
n_label <- as_labeller(c("50"   ="n[D]==50",
                         "100"  ="n[D]==100",
                         "385"  ="n[D]==385",
                         "500"  ="n[D]==500",
                         "1000" ="n[D]==1000",
                         "5000" ="n[D]==5000"),
                       default = label_parsed)

OUT %>% ggplot(aes(y = value, x = prob)) + 
  geom_jitter(alpha = alpha_val, 
              size = 1, colour = 'darkgrey') +
  #facet_nested(~title + N.DEV) + 
  facet_wrap(~N.DEV, nrow = 2, labeller = n_label) + 
  xlab("True risk for nine individuals’ ") +
  ylab("Predictions for the same individual") +
  stat_summary(
    fun      = median, 
    geom     = "errorbar",
    aes(ymax = ..y.., ymin = ..y..), 
    position = position_dodge(width = 0.8), 
    width    = 0.5,
    colour   = 'black',
    size     = 1.5) + 
  theme_bw() +
  theme(legend.position = "bottom") + 
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text = element_text(size = 7)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

