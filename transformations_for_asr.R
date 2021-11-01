# transformations_for_asr.R

# transform parameters from normal space to ASR param space
transform_to_asr <- function(params){
  params_new <- exp(params) # exponentiate all parameters
  names(params_new) <- param_names
  params_new
}

# transform parameters from ASR param space to normal space
transform_to_normal <- function(params){
  params_new <- log(params)  # log all parameters
  names(params_new) <- param_names
  params_new
}
