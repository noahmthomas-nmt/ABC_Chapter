# Rejection base ABC MCMC ASR.R

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
num_posts <- 1000000 # number posterior draws
# create array for theta
theta <- array(NA, dim = c(length(param_names), num_posts)) 
rownames(theta) <- param_names

log_weight <- matrix(NA, num_posts) # create matrix for log_weights
eps <- 65 # determine maximum difference between simulated and observed data
sigma_proposal <- start_sds / 5 # determine sd or proposal distribution

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
d <- Inf
while(d >= eps){
  theta_star <- sample_prior() # sample from prior
  print(d)
  if(is.finite(dens_prior(theta_star, LOG = T))){
    # get difference between observed and similar data
    d <- get_d(params = theta_star, data = Data[[1]]) 
  }
}
# store the first accepted proposal
theta[, 1] <- theta_star

# ------------------------------------------------------------- Sample
I <- 3 # intialize iter
# sample posterior
for (i in (I-1):num_posts) { 
  theta_1 <- rproposal(theta[, i-1], sigma_proposal) # sample from proposal
  
  # get difference between simulated data and 
  # observed for this set of parameters
  d_new <- get_d(theta_1, Data) 
  if(d_new <= eps){
    
    # determine weight based on prior and other candidate
    new_weight <- dens_prior(theta_1, LOG = T) + 
      log_dens_proposal(theta[, i-1], theta_1)
    old_weight <- dens_prior(theta[, i-1], LOG = T) + 
      log_dens_proposal(theta_1, theta[, i-1])
    
    # preform MH step
    MH <- mh_step(new_weight = new_weight, old_weight = old_weight)
    if(MH == "accept"){
      theta[, i] <- theta_1 # store the accepted value  
    }else{
      theta[, i] <- theta[, i-1] # store the previous value
    }
  }else{
    theta[, i] <- theta[, i-1] # store the previous value
  }
  print(i)
}
I <- i # record where you leave off if interrupting

# load("chains.RData") # save work

# --------------------------------------------- prior/posterior density plots

# Extract posteriors for each param
alpha <- theta[1,]
beta <- theta[2,]
mu <- theta[3,]
sigma <- theta[4,]
lambda <- theta[5,]
# ground truth
theta <- c(.0075,.01,350,50,100)


# make plots
png("../reject_alpha.png",3300,1100)
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

png("../reject_mu.png",2200,1100)
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
