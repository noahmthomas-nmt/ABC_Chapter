# Saved as PDA.R

# Creates an approximate log likelihood function

# Probability Density Approximation
#
#  Approximate a density function with a r sampling function

# data - observed data
# bandwidth - kernel bandwidth
# num_samples - number of samples for approximation
# rsampler - (psuedo)random number generator
approximate_dens <- function(data, num_samples, sample_data) {
  # Kernel density approximation
  dens <- try(density(sample_data, 
                      n = 1024, 
                      na.rm = T))
  # make sure densities are >= 0
  dens$y[dens$y < 0] <- 0
  # re-scale approximate pdf relative to probability of choice (one / one here)
  dens$y <- dens$y * length(sample_data) / num_samples
  
  out <- numeric(length(data))
  # check that observed data is within the range of the simulated data
  data_in_sample <- (data > dens$x[1]) & (data < dens$x[length(dens$x)])
  # approximate density for data within that range
  out[data_in_sample] <- approx(dens$x, dens$y, data[data_in_sample])$y
  # store 0 density for all other data
  out[is.na(out) | !is.finite(out)] <- 0
  ifelse(is.finite(out), out, -Inf) # handle strange value
}
approximate_dens <- cmpfun(approximate_dens) # compile for speed


# Creates an approximate log likelihood function
 
# Uses a kernel on simulated to approximate a likelihood
 
# condition_stim_data - list of rt and resp experimental data
# params - named vector of parameters to simulate simulate_data
# num_samples - number of samples to simulate for a kernel
# max_counter - maximum amount of steps to simulate simulate_data
# dt - dt of simulated process
# congruency - congruency of trials
# bandwidth - bandwith for kernel estimate
# simulate_data - simulation function of choice response time model
pda <- function(condition_RT_data,
                params,
                inc,
                num_samples,
                simulate_data_fn,
                transform_fn = function(x) {
                  (x)
                }) {
  
  PDF <- list() # initial approximate pdf list
  # simulate from model under params
  rt_samp_vec <- simulate_data_fn(n.obs = num_samples/2, theta = params, inc = inc)
  rt_obs_vec <- condition_RT_data # get observed data
  
  
  if (length(rt_obs_vec) == 0) {
    # if there are no observation, store no density
    PDF[[1]] <- NULL
  } else if (length(rt_samp_vec) < 2) {
    # if simulation does not produce enough samples, store 0 density
    PDF[[use_resp]] <- rep(0, length(rt_obs_vec))
  } else {
    # get range of observed data
    m <- min(rt_obs_vec)
    M <- max(rt_obs_vec)
    
    if ((min(rt_samp_vec) > M) | (max(rt_samp_vec) < m)) {
      # if the simulated data does not overlap with the observed data,
      # store 0 density
      PDF[[1]] <- rep(0, length(rt_obs_vec))
    } else {
      # approximate PDF
      PDF[[1]] <- approximate_dens(
        data = transform_fn(rt_obs_vec),
        num_samples = num_samples,
        sample_data = transform_fn(rt_samp_vec)
      )
    }
  }
  # return PDF for observations and 
  # number of simulations that failed to respond (none here)
  list(PDF, sum(sim_data$resp == -1))
}
pda <- cmpfun(pda) # compile for speed
