# rejection likelihood MCMC ASR.R

rm(list = ls()) # Clear environment
# ------------------------------------------------------------- Set up

set.seed(43210) # set seed for reproducibility
library("tidyverse") # import tidy function
library("MCMCpack") # MH sampler algorithm
library("brms") # Probability density functions
# Transform parameters to/from normal or ASR param space
source("transformations_for_asr.R")
source("priors_asr.R") # priors for ASR model
source("limits.r") # for plotting

# ------------------------------------------------------------- Import Data
Data <- read_table2("Conflict_Data_2.0 (1).txt",
  col_types = cols(
    Inc = col_integer(),
    SOA = col_integer(), Sub = col_integer()
  )
)
# Reshape Data
data_list <- list()
data_list[[1]] <- Data %>%
  dplyr::select(Inc, RT) %>%
  as.list()
Data <- data_list

# Density function for ASR
dasr <- function(t, theta, inc = 0, soa = 0) {
  alpha <- 1 / theta["alpha"] # Exponential scale for A
  beta <- 1 / theta["beta"] # Exponential scale for B
  mu <- theta["mu"] # Mean for C
  sigma <- theta["sigma"] # SD for C
  lambda <- theta["lambda"] # Incongruent delay
  p <- beta / (alpha + beta) # P(B<A)
  
  # Dexgaussian uses the scale, not the rate
  term_1 <- (1 - inc) * dexgaussian(t, mu + 1 / beta, sigma, 1 / beta)
  part_1 <- (p * dexgaussian(t, 
                             mu + 1 / (alpha + beta) + lambda, 
                             sigma, 
                             1 / (alpha + beta)))
  part_2 <- (1 - p) * ((1 + beta / alpha) * dexgaussian(t, 
                                                        mu + 1 / beta, 
                                                        sigma, 1 / beta) -
               (beta / alpha) * dexgaussian(t, 
                                            mu + 1 / (alpha + beta), 
                                            sigma, 
                                            1 / (alpha + beta)))
  
  
  term_2 <- inc * (part_1 + part_2)
  pdf <- term_1 + term_2
  return(pdf) # return density
}

# log posterior function
log_post <- function(theta, data) {
  names(theta) <- param_names # get param names

  ldp <- 0 # add log prior function here

  lds <- 0 # total log densities
  # transform theta from normal space to ASR param space
  theta_ <- transform_to_asr(theta)
  # check that prior has density
  if (!is.finite(ldp)) {
    # add lowest amount of log density in R for each observation
    lds <- lds + rep(-740, length(Data[[1]]$RT))
    lds <- lds + -740
  } else {
    # Evaluate likelihood for congruent and incongruent
    ld_con <- dasr(t = Data[[1]]$RT[Data[[1]]$Inc == 0], 
                   theta_, 
                   inc = 0, 
                   soa = 0)
    ld_incon <- dasr(t = Data[[1]]$RT[Data[[1]]$Inc == 1], 
                     theta_, 
                     inc = 1, 
                     soa = 0)
    
    # Be sure density is >= 0
    ld_incon <- pmax(ld_incon, 0)
    ld_con <- pmax(ld_con, 0)
    
    # Add log densities
    lds <- lds + log(ld_con) + log(ld_incon)
  }
  
  # Add all log densities and return
  out <- sum(lds + ldp)
  return(out)
}

# ------------------------------------------------------------- Initialize
theta_init <- mu_mean_vec # Initial starting point (at mean of prior)

# ------------------------------------------------------------- Sample
mcmc <- 100000 # number of mcmc draws
# Metropolis Hastings sampler
out <- MCMCmetrop1R(
  fun = log_post, # posterior function to sample
  theta.init = theta_init, # Starting point
  burnin = mcmc * .1, # number of draws to discard
  mcmc = mcmc, # number of mcmc draws
  thin = 1, # thinning interval used in the simulation
  tune = 1.3, # sd of proposal distribution
  seed = 1,  # seed
  logfun = T, 
  optim.method = "Nelder-Mead",
  force.samp = T,
  verbose = T
)

# --------------------------------------------- prior/posterior density plots


# Extract posteriors for each param
alpha <- out[,1]
beta <- out[,2]
mu <- out[,3]
sigma <- out[,4]
lambda <- out[,5]
# ground truth
theta <- c(.0075,.01,350,50,100)

# make plots
png("../MCMC_alpha.png",3300,1100)
par(mfrow=c(1,3),cex=2.5)
cs.alpha <- quantile(alpha,c(0,.025,.975,1))
range.alpha <- alpha.y[2]-alpha.y[1]
plot(density(alpha),xlim=alpha.x,ylim=alpha.y,
     xlab="Log Alpha",lwd=3,cex=1.5,main='',ylab="")
lines(rep(log(1/theta[1]),2),alpha.y,lty=3)
lines(cs.alpha[2:3],rep(-.005*range.alpha,2),lwd=5)
abline(h=1,lwd=1)

cs.beta <- quantile(beta,c(0,.025,.975,1))
range.beta <- beta.y[2]-beta.y[1]
plot(density(beta),xlim=beta.x,ylim=beta.y,
     xlab="Log Beta",lwd=3,cex=1.5,main='',ylab="")
lines(rep(log(1/theta[2]),2),beta.y,lty=3)
lines(cs.beta[2:3],rep(-.005*range.beta,2),lwd=5)
abline(h=1,lwd=1)

cs.lambda <- quantile(lambda,c(0,.025,.975,1))
range.lambda <- lambda.y[2]-lambda.y[1]
plot(density(lambda),xlim=lambda.x,ylim=lambda.y,
     xlab="Log Lambda",lwd=3,cex=1.5,main='',ylab="")
lines(rep(log(theta[5]),2),lambda.y,lty=3)
lines(cs.lambda[2:3],rep(-.005*range.lambda,2),lwd=5)
abline(h=1,lwd=1)
dev.off()

png("../MCMC_mu.png",2200,1100)
par(mfrow=c(1,2),cex=2.5)

cs.mu <- quantile(mu,c(0,.025,.975,1))
range.mu <- mu.y[2]-mu.y[1]
plot(density(mu),xlim=mu.x,ylim=mu.y,
     xlab="Log Mu",lwd=3,cex=1.5,main='',ylab="")
lines(rep(log(theta[3]),2),mu.y,lty=3)
lines(cs.mu[2:3],rep(-.005*range.mu,2),lwd=5)
abline(h=1,lwd=1)

cs.sigma <- quantile(sigma,c(0,.025,.975,1))
range.sigma <- sigma.y[2]-sigma.y[1]
plot(density(sigma),xlim=sigma.x,ylim=sigma.y,
     xlab="Log Sigma",lwd=3,cex=1.5,main='',ylab="")
lines(rep(log(theta[4]),2),sigma.y,lty=3)
lines(cs.sigma[2:3],rep(-.005*range.sigma,2),lwd=5)
abline(h=1,lwd=1)
dev.off()

