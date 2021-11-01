# Saved as rconflict.R

# Simulates a conflict model

# Generates a synthetic experiment from a table of numbers for conditions

# condition_table - table containing the numbers of trials for each condition 
#   and stimulus direction
# params - vector of parameters to be freely estimated
# param_names - vector of names of parameters to be freely estimated
# params_fixed - vector of parameters to be held fixed
# param_names_fixed - vector of names of parameters to be held fixed
# simulate_fn  - function to simulate conflict model

rconflict <- function(condition_table, 
                      params, 
                      param_names, 
                      params_fixed, 
                      param_names_fixed, 
                      simulate_fn = simulate_DMC, 
                      transform_fn = transform_to_asr) {
  names(params) <- param_names
  out <- list(RT = NULL, Inc = NULL)
  
  for (cond in conds) {
      
      # get number of samples for this condition
      N <- condition_table[cond+1]
      
      # Transform parameters to ASR param space
      params_tmp <- transform_fn(params)
      # Add fixed parameters (None here)
      params_tmp[param_names_fixed] <- params_fixed
      
      # simulate ASR
      tmp <- simulate_fn(n.obs = N/2, theta = params_tmp, inc = cond)
      
      # get RT and condtion
      out$RT <- c(out$RT, tmp)
      out$Inc <- c(out$Inc, rep(cond, condition_table[cond+1]))
    }
  out # return list of RTs and conditions
}