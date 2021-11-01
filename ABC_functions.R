# ABD_functions.R

# sample from prior
# HIER - logical, whether to sample hierarchically or not
sample_prior <- function(HIER = F) {
  # intialize prior samples
  prior_vec <- numeric(length(param_names))
  names(prior_vec) <- param_names
  
  if(HIER){
    # get sample from hierarchical structure
    for(p_name in param_names){
      prior <- prior_list[[p_name]] # get parameter of interest
      # sample from mean distribution of parameter
      mu_i <- rnorm(1, prior$mu[1], prior$mu[2])  
      # sample from standard deviation distribution of parameter
      sigma_i <- rgamma(1, prior$sigma[1], prior$sigma[2])  
      p_i <- rnorm(1, mu_i, sigma_i) # sample parameter
      prior_vec[p_name] <- p_i # store parameter
    }
  }else{
    # get sample from simple prior structure
    for(p_name in param_names){
      prior <- prior_list[[p_name]] # get parameter of interest
      # sample from mean distribution of parameter
      p_i <- rnorm(1, prior$mu[1], prior$mu[2])  
      prior_vec[p_name] <- p_i # store parameter
    }
  }
  prior_vec # return vector of parameters
}

# rho difference function (# here chi-squared similarity)
# x     - simulated data
# y     - observed data
# scale - proportion of simulated data with this response
rho <- function(x, y, scale) { 
  # determine deciles of simulated data
  deciles <- seq(.05, .95, by = .1) 
  decile_breaks <- sort(x)[length(x)*deciles]
  
  # no deciles
  if (length(x) < 1){
    x2_star <- NA
  }else if(length(x) > length(deciles)+1){ 
    # If all deciles are unique compute similarity for each bin
    # .0 to .05 quantile
    x2_star <- (sum(y < decile_breaks[1]) - 
                  .05 * length(x) * scale)^2 / (.05 * length(x) * scale)
    # .05 to .95 deciles
    for(d in seq_along(decile_breaks)[-1]){
      x2_star <- x2_star + 
        (sum(y < decile_breaks[d] & y > decile_breaks[d-1]) - 
           .1 * length(x) * scale)^2 / (.1 * length(x) * scale)
    }
    # .95 to 1.0 quantile
    x2_star <- x2_star + 
      (sum(y > decile_breaks[length(decile_breaks)]) - 
         .05 * length(x) * scale)^2 / (.05 * length(x) * scale)
  }else{
    # If all deciles are not unique -> bins may shift
    x2_star <- (sum(y < decile_breaks[1]) - 1)^2 / scale
    
    for(d in seq_along(decile_breaks)[-1]){
      x2_star <- x2_star + 
        (sum(y < decile_breaks[d] & y > decile_breaks[d-1]) - 1 * scale)^2 / 
        scale 
    }
    
    x2_star <- x2_star + 
      (sum(y > decile_breaks[length(decile_breaks)]) - 1 * scale)^2 / scale 
  }
  x2_star # return chi-squared similarity
}

# Get difference between parameter's simulated data and observed data
get_d <- function(params, data) {
  # simulate data 
  sim_data <<-  try(rconflict(condition_table = condition_table, 
                              params = params, 
                              param_names = param_names, 
                              params_fixed = params_fixed, 
                              param_names_fixed = param_names_fixed, 
                              simulate_fn = rasr, 
                              transform_fn = transform_to_asr
  )
  )
  
  # return Inf for computational issues: too large/small parameters, NAs, ...
  if(is.character(sim_data)){return(Inf)}
  d <- 0 # initiate difference
  for(cond in conds){
    # get simulated data that matches this condition and stimuli
    filt_sim_data <- sim_data %>% 
      as_data_frame() %>% 
      filter(Inc == cond)
    
    # get observed data that matches this condition and stimuli
    filt_data <- data %>% 
      as_data_frame() %>% 
      filter(Inc == cond)
    
    # ---------------------------------- Get vector of times for right responses
    filt_data <- filt_data %>% 
      pull(RT)
    
    filt_sim_data <- filt_sim_data %>% 
      pull(RT)
    
    #  -------------------- Compare simulated v. observed data for each response
    rho <- rho(filt_sim_data, filt_data, length(filt_sim_data)/
                 (length(filt_sim_data) + length(filt_sim_data)))
    
    d <- d + rho # compute the the joint differences
    # here we add the rescaled chi-squared similarities
  }
  if(!is.finite(d)){d <-  Inf}
  d
}

# Creates a new proposed theta vector based on theta_star
rproposal <- function(theta_star, sigma_vec) {
  sample_vec <- numeric(length(theta_star)) # initiate vector
  for(p in seq_along(theta_star)){
    # sample from distribution centered at theta and 
    # with standard deviation sigma
    # note: you may have different sigma for each parameter
    sample_vec[p] <- rnorm(1, theta_star[p], sigma_vec[p]) 
  }
  sample_vec # return new theta
}

# Calculate density of paramaters under prior
dens_prior <- function(params, LOG = T) {
  names(params) <- param_names
  if(LOG){
    out <- 0 # initiate log density 
    for (p in param_names) {
      out <- out + dnorm(params[p], 
                         prior_list[[p]]$mu[1], # mean
                         prior_list[[p]]$mu[2], # sd
                         log = LOG
      )
    }
  }else{
    out <- 1 # initiate density 
    for (p in param_names) {
      out <- out * dnorm(params[p], 
                         prior_list[[p]]$mu[1], 
                         prior_list[[p]]$mu[2],
                         log = LOG
      )
    }
  }
  out # return (log) density
}

# function to compute the metropolis hastings step
mh_step <- function(new_weight, old_weight) {
  out <- "reject"
  
  # check for bad weights
  if (is.null(new_weight)) new_weight <- -Inf
  if (!is.na(new_weight) & !is.finite(new_weight) & new_weight < 0){
    new_weight <- -Inf
  } 
  if (is.na(new_weight)) new_weight <- -Inf
  
  if (is.null(old_weight)) old_weight <- -Inf
  if (!is.na(old_weight) & !is.finite(old_weight) & old_weight < 0){
    old_weight <- -Inf
  } 
  if (is.na(old_weight)) old_weight <- -Inf
  
  # perform MH step comparing the log weights
  mh <- exp(new_weight - old_weight)
  if (is.finite(mh)) {
    if (runif(1) < mh) {
      out <- "accept"
    }
  }else if (is.infinite(mh) & mh > 1){
    out <- "accept"
  }
  out # return decision
}

# compute the log density of proposal based on the former parameters
log_dens_proposal <- function(theta_0, theta_star){
  for(p in 1:length(theta_0)){
    out <- 0
    out <- out + dnorm(theta_0[p], theta_star[p], sigma_proposal[p], log = T)
  }
  out
}
