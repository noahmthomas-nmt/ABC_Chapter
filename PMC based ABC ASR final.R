# PMC based ABC ASR.R

rm(list = ls()) # Clear environment

# ------------------------------------------------------------- Set up

set.seed(43210) # set seed for reproducibility
library("tidyverse") # import tidy function
source("rasr.R") # simulate ASR
source("rconflict_asr.R") # Runs conflict model with experimental design
# Transform parameters to/from normal or ASR param space
source("transformations_for_asr.R")
source("priors_asr.R") # priors for ASR model
source("ABC_functions.R") # functions to run ABC
source("limits.r") # for plotting

# --------------------------------------------------------- Algorithm Parameters
num_posts <- 1000 # number of posterior draws
num_particles <- 100 # number of particles in population
# array for posterior draws
theta <- array(NA, 
               dim = c(length(param_names), 
                       num_particles, 
                       num_posts)) 
rownames(theta) <- param_names 
log_weight <- matrix(NA, num_particles, num_posts) # array for log weights
eps <- seq(100, 25, length.out = num_posts) # vector of acceptable differences
# vector for variation of the population
sigma2_vec <- matrix(NA, length(param_names), num_posts) 

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


# ------------------------------------------------------------- Initialization
for(part_i in 1:num_particles){
  d <- Inf
  while(d > eps[1]){
    print(d)
    theta_star <- sample_prior() # sample from prior
    # check that the parameters have density under the prior
    if(is.finite(dens_prior(theta_star, LOG = T)) & 
       dens_prior(theta_star, LOG = F) > 0){
      # get the difference between simulated and observed data
      d <- get_d(params = theta_star, data = Data) 
    }
  }
  print(part_i)
  theta[, part_i, 1] <- theta_star # store result
}

log_weight[, 1] <- log(1/num_particles) # initialize log weights
# calulate variation of the population
sigma2_vec[, 1] <- 2 * apply(theta[,, 1], 1, var) 

# ------------------------------------------------------------- Sample
I <- 3 # initialize I
for(i in (I - 1):num_posts){
  print(i)
  for(part_i in 1:num_particles){
    d <- Inf
    while(d > eps[i]){ # while d is greater that acceptable difference
      # sample particles based on fitness
      theta_star_ind <-  sample(1:num_particles, 
                                1, 
                                prob = exp(log_weight[, i-1]) + exp(-740)) 
      theta_star <- theta[, theta_star_ind, i-1] # get particle that was sampled
      # generate proposal around that particle
      theta_ss <- rproposal(theta_star, sigma2_vec[, i-1]) 
      # check that the parameters have density under the prior
      if(is.finite(dens_prior(theta_ss, LOG = T)) & 
         dens_prior(theta_ss, LOG = F) > 0 ){
        # get the difference between simulated and observed data
        d <- get_d(params = theta_ss, data = Data) 
      }
    }
    theta[, part_i, i] <- theta_ss # store result
    
    # log(pi/sum(w*dens))
    log_weight[part_i, i] <- 
      log(
        dens_prior(theta[, part_i, i], LOG = F) / 
          sum(dnorm(theta[,, i-1], 
                    theta[, part_i, i], 
                    sqrt(sigma2_vec[, i-1]), 
                    log = F) %*% 
                exp(matrix(log_weight[, i-1], ncol = 1)), na.rm = T)
      )
  }
  # calculate variation in population
  sigma2_vec[, i] <- 2 * apply(theta[,, i], 1, var) 
}
I <- i-1 # store i

load("PMC_chains.Rdata") # save work

# --------------------------------------------- prior/posterior density plots

# Extract posteriors for each param
alpha <- theta[1,,]
beta <- theta[2,,]
mu <- theta[3,,]
sigma <- theta[4,,]
lambda <- theta[5,,]
# ground truth
theta <- c(.0075,.01,350,50,100)

# make plots
png("../PMC_alpha.png",3300,1100)
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

png("../PMC_mu.png",2200,1100)
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
