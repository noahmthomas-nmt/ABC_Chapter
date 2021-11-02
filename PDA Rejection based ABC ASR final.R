# PDA Rejection based ABC ASR.R

rm(list = ls()) # Clear environment

# ------------------------------------------------------------- Set up

set.seed(43210) # set seed for reproducibility
library("tidyverse") # import tidy function
library("MCMCpack") # MH sampler algorithm
library("compiler") # compiles R functions for speed
source("rasr.R") # simulate ASR
source("rconflict_asr.R") # Runs conflict model with experimental design
# Transform parameters to/from normal or ASR param space
source("transformations_for_asr.R")
source("priors_asr.R") # priors for ASR model
source("ABC_functions.R") # functions to run ABC
source(("PDA final.R")) # Calculates approximate likelihood
source("limits.r") # for plotting

# ---------------------------------------------------- Experiment Parameters

conds <- c(0, 1) # names of conditions

n <- 2500 # number of trials for cond (comp/incomp) and stim (left/right)
condition_table <- c(n, n) # make a table of the condition/stimuli
names(condition_table) <- conds

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

# log posterior function
log_post <- function(theta, data) {
  names(theta) <- param_names # get param names
  
  ldp <- 0 # add log prior function here
  
  lds <- 0 # total log densities
  # transform theta from normal space to ASR param space
  theta_ <- transform_to_asr(theta)
  # check that prior has density
  if(!is.finite(ldp)){
    # add lowest amount of log density in R for each observation
    lds <-  lds + rep(-740, length(Data[[1]]$RT))
    lds <- lds + -740
  }else{
    # Evaluate approximate likelihood for congruent
    tmp1 <- pda(condition_RT_data = Data[[1]]$RT[Data[[1]]$Inc == 0], 
        params = theta_, 
        num_samples = num_samples, 
        simulate_data_fn = rasr,  
        inc=0)
    
    # If pda returns only densities
    if (length(tmp1[[1]]) == 1) {
      tmp2 <- tmp1[[1]][[1]] # store densities
    } else {
      # store densities and RTs that have no postive density (coded -1)
      tmp2 <- c(tmp1[[1]][[1]], tmp1[[1]][[2]])
    }
    ld_con <- tmp2  # store densities for congruent condition
    
    # Evaluate approximate likelihood for congruent
    tmp1 <- pda(condition_RT_data = Data[[1]]$RT[Data[[1]]$Inc == 1], 
                params = theta_, 
                num_samples = num_samples, 
                simulate_data_fn = rasr,  
                inc=1)
    
    # If pda returns only densities
    if (length(tmp1[[1]]) == 1) {
      tmp2 <- tmp1[[1]][[1]] # store densities
    } else {
      # store densities and RTs that have no postive density (coded -1)
      tmp2 <- c(tmp1[[1]][[1]], tmp1[[1]][[2]])
    }
    ld_incon <- tmp2 # store densities for incongruent condition
    
    # Be sure density is >= 0
    ld_incon <- pmax(ld_incon, 0)
    ld_con <- pmax(ld_con, 0)
    
    # Add log densities
    lds <- lds + log(ld_con) + log(ld_incon)  
  }
  
  # Add all log densities and return
  out <- sum(lds + ldp)
  out
}

# ------------------------------------------------------------- Initialize
theta_init <- mu_mean_vec # Initial starting point (at mean of prior)

# ------------------------------------------------------------- Sample
# Covariance matrix for proposal distribution
Cov_mat <- matrix(c(5.2e-03,  3e-04, -2e-04, -4e-04, -6e-05,
           3e-04,  1.1e-03, -2e-04, -7e-04,  2.8e-04,
           -2e-04, -2e-04,  8.8e-05,  2.5e-04, -2.5e-04,
           -4.2e-04, -7.3e-04,  2.5e-04,  1.7e-03, -8e-04,
           -6e-05,  2.8e-04, -2.5e-04, -8e-04,  2.1e-03), 
       nrow = 5, 
       ncol = 5)


mcmc <- 1e6  # number of mcmc draws
num_samples <- 50000  # number of simulations for pda

debugonce(log_post)
out <- MCMCmetrop1R(fun = log_post, 
                    theta.init = theta_init, 
                    burnin = mcmc*1.1, 
                    mcmc = mcmc, 
                    thin = 1, 
                    V = Cov_mat,
                    seed = 1, 
                    logfun = T, 
                    optim.method = "Nelder-Mead", 
                    force.samp = T, 
                    verbose = T)


load("PDA_chains.RData") # save work

# --------------------------------------------- prior/posterior density plots

# Extract posteriors for each param
alpha <- out2[,1]
beta <- out2[,2]
mu <- out2[,3]
sigma <- out2[,4]
lambda <- out2[,5]
# ground truth
theta <- c(.0075,.01,350,50,100)

# make plots
png("../PDA_alpha.png",3300,1100)
par(mfrow=c(1,3),cex=2.5)
cs.alpha <- quantile(alpha,c(0,.025,.975,1))
range.alpha <- alpha.y[2]-alpha.y[1]
plot(density(alpha),xlim=alpha.x,ylim=alpha.y,
     xlab="Log Alpha",lwd=3,cex=1.5,main='',ylab="")
lines(rep(log(1/theta[1]),2),alpha.y,lty=3)
lines(cs.alpha[2:3],rep(-.005*range.alpha,2),lwd=5)
x <- seq(alpha.x[1], alpha.x[2], by = .001)
lines(x, dnorm(x, prior_list$alpha$mu[1], prior_list$alpha$mu[2]))

cs.beta <- quantile(beta,c(0,.025,.975,1))
range.beta <- beta.y[2]-beta.y[1]
plot(density(beta),xlim=beta.x,ylim=beta.y,
     xlab="Log Beta",lwd=3,cex=1.5,main='',ylab="")
lines(rep(log(1/theta[2]),2),beta.y,lty=3)
lines(cs.beta[2:3],rep(-.005*range.beta,2),lwd=5)
x <- seq(beta.x[1], beta.x[2], by = .001)
lines(x, dnorm(x, prior_list$beta$mu[1], prior_list$beta$mu[2]))

cs.lambda <- quantile(lambda,c(0,.025,.975,1))
range.lambda <- lambda.y[2]-lambda.y[1]
plot(density(lambda),xlim=lambda.x,ylim=lambda.y,
     xlab="Log Lambda",lwd=3,cex=1.5,main='',ylab="")
lines(rep(log(theta[5]),2),lambda.y,lty=3)
lines(cs.lambda[2:3],rep(-.005*range.lambda,2),lwd=5)
x <- seq(lambda.x[1], lambda.x[2], by = .001)
lines(x, dnorm(x, prior_list$lambda$mu[1], prior_list$lambda$mu[2]))
dev.off()

png("../PDA_mu.png",2200,1100)
par(mfrow=c(1,2),cex=2.5)

cs.mu <- quantile(mu,c(0,.025,.975,1))
range.mu <- mu.y[2]-mu.y[1]
plot(density(mu),xlim=mu.x,ylim=mu.y,
     xlab="Log Mu",lwd=3,cex=1.5,main='',ylab="")
lines(rep(log(theta[3]),2),mu.y,lty=3)
lines(cs.mu[2:3],rep(-.005*range.mu,2),lwd=5)
x <- seq(mu.x[1], mu.x[2], by = .001)
lines(x, dnorm(x, prior_list$mu$mu[1], prior_list$mu$mu[2]))

cs.sigma <- quantile(sigma,c(0,.025,.975,1))
range.sigma <- sigma.y[2]-sigma.y[1]
plot(density(sigma),xlim=sigma.x,ylim=sigma.y,
     xlab="Log Sigma",lwd=3,cex=1.5,main='',ylab="")
lines(rep(log(theta[4]),2),sigma.y,lty=3)
lines(cs.sigma[2:3],rep(-.005*range.sigma,2),lwd=5)
x <- seq(sigma.x[1], sigma.x[2], by = .001)
lines(x, dnorm(x, prior_list$sigma$mu[1], prior_list$sigma$mu[2]))
dev.off()
