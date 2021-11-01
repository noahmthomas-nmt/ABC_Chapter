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

# Sample plots
plot(out)

# Chains
ts.plot(exp(out[1:mcmc, 1]))
ts.plot(exp(out[1:mcmc, 2]))
ts.plot(exp(out[1:mcmc, 3]))
ts.plot(exp(out[1:mcmc, 4]))
ts.plot(exp(out[1:mcmc, 5]))

# means
post_mean <- apply(out, 2, function(x) mean(exp(x)))
names(post_mean) <- param_names
post_mean

# joint distributions
psych::pairs.panels(out %>% as.matrix())
