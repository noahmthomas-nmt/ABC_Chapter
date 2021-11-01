# generate RTs for one subject with parameter vector theta
# n.obs is the number of trials in each condition (congruent or incongruent)
# inc is a vector of length 1 (all observations are from one condition) or 
#   length 2*n.obs (n.obs in each condition)
# theta is the model parameter vector as unpacked in the function
# soa is included for completeness but we won’t use it
rasr <- function(n.obs,theta,inc=0,soa=0) {
  if (length(inc)>1 & length(inc)!=2*n.obs){ 
    return("There’s a problem with the number of trials 
            or the length of the incongruent variable.")
  }
  
  # Pull out the different variables for readability
  alpha <- theta["alpha"]  # Exponential rate for A
  beta <- theta["beta"]   # Exponential rate for B
  mu <- theta["mu"]    # Mean for C
  sigma <- theta["sigma"]  # SD for C
  lambda <- theta["lambda"]  # Incongruent delay
  
  # Generate A, taking into account any SOA
  A <- rexp(2*n.obs)*alpha - soa
  
  # Generate B
  B <- rexp(2*n.obs)*beta
  
  # Generate C, with delay of lambda if A exceeds B for incongruent trials
  C <- rnorm(2*n.obs,mu + (B<A)*inc*lambda,sigma)
  
  # RT is equal to B + C, round to nearest whole number
  T <- round(B + C,0)
  return(T)
}
