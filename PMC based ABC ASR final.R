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

# inspect posterior means and variances
apply(theta[, , I], MARGIN = 1, mean)
apply(theta[, , I], MARGIN = 1, var)

# inspect chains
ts.plot(t(exp(theta["alpha", , 1:I])), ylim = c(0, 200))
plot(density(exp(theta["alpha", , I])))

ts.plot(t(exp(theta["beta", ,1:I])))
plot(density(exp(theta["beta", , I])))

ts.plot(t(exp(theta["mu", ,1:I])))
plot(density(exp(theta["mu", , I])))

ts.plot(t(exp(theta["sigma", ,1:I])))
plot(density(exp(theta["sigma", , I])))

ts.plot(t(exp(theta["lambda", ,1:I])))
plot(density(exp(theta["lambda", , I])))

# load("ABC PMC ASR Super Subject 1 with Priors More Samples to eta 7 current.Rdata") # save work
