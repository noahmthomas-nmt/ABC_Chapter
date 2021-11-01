# Saved as Priors.R

# helper function to make a list for priors
update_prior_bounds_lists <- function(param_name,
                                      mu_mean,
                                      mu_sd,
                                      sigma_shape,
                                      sigma_rate,
                                      lb,
                                      ub,
                                      param_names,
                                      prior_list,
                                      lower_bound_list,
                                      upper_bound_list) {
  
  # get parameter
  tmp <- grep(param_name, param_names, value = TRUE)
  
  # fill in parameter to prior list
  for (n in 1:length(tmp)) {
    tmp2 <- tmp[n]
    prior_list[[tmp2]] <- list(mu = c(mu_mean, mu_sd), 
                               sigma = c(sigma_shape, sigma_rate))
    lower_bound_list[[tmp2]] <- lb
    upper_bound_list[[tmp2]] <- ub
  }
  
  # return list of priors, lower and upper bounds
  list(prior_list = prior_list, 
       lower_bound_list = lower_bound_list, 
       upper_bound_list = upper_bound_list)
}


# -------------------------------------------------------- Parameters
# means
alpha <- 100 
beta <- 100
mu <- 300
sigma <- 60
lambda <- 75

# names of parameters
param_names <- c("alpha", "beta", "mu", "sigma", "lambda")

# parameters to hold fixed (none here)
param_names_fixed <- NULL
params_fixed <- NULL

# -------------------------------------------------------- Priors
# initiate list
prior_list <- list()
upper_bound_list <- list()
lower_bound_list <- list()

# means (in normal space)
start_points <- mu_mean_vec <- log(c(alpha, beta, mu, sigma, lambda))  
# standard deviation of prior
start_sds <- mu_sd_vec <-             c(.75, .75, .3,    .3,     1)

# gamma parameterization for hierarchical models (not used here)
sigma_shape_vec  <- rep(1, length(param_names))
sigma_rate_vec   <- rep(1, length(param_names))

# Bounds on normal parameters (unbounded)
lower_bounds_vec <- rep(-Inf, length(param_names))
upper_bounds_vec <- rep(Inf, length(param_names))

# build prior list for each parameter
for (i in seq_along(param_names)) {
  out <- update_prior_bounds_lists(
    param_name = param_names[i],
    mu_mean = mu_mean_vec[i],
    mu_sd = mu_sd_vec[i],
    sigma_shape = sigma_shape_vec[i],
    sigma_rate = sigma_rate_vec[i],
    lb = lower_bounds_vec[i],
    ub = upper_bounds_vec[i],
    param_names = param_names,
    prior_list = prior_list,
    lower_bound_list = lower_bound_list,
    upper_bound_list = upper_bound_list
  )
  prior_list <- out$prior_list
  upper_bound_list <- out$upper_bound_list
  lower_bound_list <- out$lower_bound_list
}
