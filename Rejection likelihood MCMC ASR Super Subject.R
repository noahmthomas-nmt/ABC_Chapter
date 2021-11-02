rm(list = ls()) # Clear environment
# ------------------------------------------------------------- Set up

# Rejection based MCMC for ASR
set.seed(43210) # set seed for reproducibility
library("tidyverse") # import tidy function
library("compiler") # import R compiler
library("MCMCpack") # MH sampler algorithm
library("msm") # Probability density functions
library("brms")  # Probability density functions
source("asr_distribution.R") # ASR probability density functions
# functions to transform parameters to normal or ASR param space
source("transformations_for_asr.R")
source("priors_asr.R") # priors for ASR model
source("ABC_functions.R") # functions to run ABC
source("log_densities_asr.R") # likehood and prior functions
source("limits.r") # for plotting

# ------------------------------------------------------------- ASR Experiment Parameters

conds <- c(0, 1) # names of conditions
condition_table <- c(n, n) # make a table of the condition/stimuli
names(condition_table) <- conds

# ------------------------------------------------------------- Initialization
Data <- read_table2("Conflict_Data_2.0 (1).txt", 
                    col_types = cols(Inc = col_integer(), 
                                     SOA = col_integer(), Sub = col_integer())) # import data
data_list <- list()
data_list[[1]] <- Data %>% 
  dplyr::select(Inc, RT) %>% 
  as.list() 
Data <- data_list


dasr <- function(t, theta, inc=0, soa=0) {
  alpha <- 1/theta["alpha"]    # Exponential scale for A
  beta <- 1/theta["beta"]     # Exponential scale for B
  mu <- theta["mu"]       # Mean for C
  sigma <- theta["sigma"]    # SD for C
  lambda <- theta["lambda"]   # Incongruent delay
  p <- beta/(alpha + beta) # P(B<A)
  # Dexgaussian uses the scale, not the rate
  pdf <- (1-inc)*dexgaussian(t,mu+1/beta,sigma,1/beta) + inc*
    (   p*
          dexgaussian(t,mu+1/(alpha+beta)+lambda,sigma,1/(alpha+beta)) +
          (1-p)*
          ((1 + beta/alpha)*
             dexgaussian(t,mu+1/beta,sigma,1/beta) -
             (beta/alpha)*
             dexgaussian(t,mu+1/(alpha+beta),sigma,1/(alpha+beta))))
  return(pdf)
}

log_post <- function(theta, data){
  names(theta) <- param_names
  
  ldp <- 0 #dens_prior(theta, LOG = T)
  
  lds <- 0
  theta_ <- transform_to_asr(theta)
  if(!is.finite(ldp)){
    lds <-  lds + rep(-740, length(Data[[1]]$RT))
    lds <- lds + -740
  }else{
    # log likelihood
    ld_con <- dasr(t = Data[[1]]$RT[Data[[1]]$Inc == 0], theta_, inc=0, soa=0)
    ld_incon <- dasr(t = Data[[1]]$RT[Data[[1]]$Inc == 1], theta_, inc=1, soa=0)
    
    ld_incon <- pmax(ld_incon, 0)
    ld_con <- pmax(ld_con, 0)
    
    lds <- lds + log(ld_con) + log(ld_incon)  
  }
  
  out <- sum(lds + ldp)
  # print(theta_)
  # print(out)
  out
}

mcmc <- 100000

theta_init <- mu_mean_vec
out <- MCMCmetrop1R(fun = log_post, 
                    theta.init = theta_init, 
                    burnin = mcmc*.1, 
                    mcmc = mcmc, 
                    thin = 1, 
                    tune = 1.3, 
                    seed = 1, 
                    logfun = T, 
                    optim.method = "Nelder-Mead", 
                    force.samp = T, 
                    verbose = T)

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


