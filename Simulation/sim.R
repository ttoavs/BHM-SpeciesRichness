library(tidyverse)
library(jagsUI)
library(progress)



######################
##### Functions  #####
######################

# Chao function
Chao <- function(matrix) {
  data <- apply(matrix, 1, sum)
  t <- ncol(matrix)
  f1 <- sum(data==1)
  f2 <- sum(data==2)
  K <- (t-1)/t
  Sobs <- sum(data>0)
  f0 <- K*ifelse(f2==0,f1*(f1-1)/2,f1^2/(2*f2))
  S <- Sobs+f0
  f2 <- max(f2,1)
  
  func1 <- function(obs,f1,f2,t){obs+(t-1)/t*f1^2/(2*f2)}
  M1 <- -c(f0,f1,f2)%*%t(c(f0,f1,f2))/S+diag(c(f0,f1,f2)) 
  c1 <- -1
  c2 <- func1(Sobs,f1+1,f2,t)-func1(Sobs,f1,f2,t)
  c3 <- func1(Sobs,f1,f2+1,t)-func1(Sobs,f1,f2,t)
  V <- c(c1,c2,c3)
  sd <- (t(V)%*%M1%*%V)^(1/2)
  R <- exp(1.96*( log(1+sd^2/max(0.1^10,(S-Sobs)^2)))^(1/2))
  U <- Sobs+(S-Sobs)*R
  L <- Sobs+(S-Sobs)/R
  
  return(c(est=S,sd=sd,upper95=U,lower95=L))
}

# Jackknife function
Jack1=function(matrix){
  data <- apply(matrix, 1, sum)
  t <- ncol(matrix)
  f1=sum(data==1);Sobs=sum(data>0);
  
  func2=function(obs,f1,t){obs+(t-1)/t*f1}
  S=f1+Sobs;f0=S-Sobs;
  M2=-c(f0,f1)%*%t(c(f0,f1))/S+diag(c(f0,f1)) 
  c1=-1;
  c2=func2(Sobs,f1+1,t)-func2(Sobs,f1,t)
  V=c(c1,c2)
  sd=(t(V)%*%M2%*%V)^(1/2)
  R <- exp(1.96*( log(1+sd^2/max(0.1^10,(S-Sobs)^2)))^(1/2))
  U <- Sobs+(S-Sobs)*R
  L <- Sobs+(S-Sobs)/R
  
  return(c(est=S,sd=sd,upper95=U,lower95=L));
}

# BMM
sink("New Sim/BMM.txt")
cat("
model {

# Likelihood

for (i in 1:n) {
    y[i] ~ dnorm(y.hat[i], tau.y) 
    y.hat[i] <- a[site[i]] * x[i]/(b[site[i]] + x[i])
 }

# Priors

for (j in 1:J) {
    a[j] ~ dnorm(mu.a, tau.a)
    b[j] ~ dnorm(mu.b, tau.b)
 }
 
 mu.a ~ dnorm(0, 0.0001)
 tau.a <- 1 / (sigma.a * sigma.a)
 sigma.a ~ dunif(0, 1000)		

 mu.b ~ dnorm(0, 0.0001)
 tau.b <- 1 / (sigma.b * sigma.b)
 sigma.b ~ dunif(0, 1000)
 
 tau.y <- 1 / (sigma.y * sigma.y)
 sigma.y ~ dunif(0,1000)
}
",fill=TRUE)
sink()

# BNE
sink("New Sim/BNE.txt")
cat("
model {

# Likelihood
for (i in 1:n) {
    y[i] ~ dnorm(y.hat[i], tau.y) 
    y.hat[i] <- a[site[i]]*(1-exp(-b[site[i]]*x[i]))
}

# Priors
for (j in 1:J) {
    a[j] ~ dnorm(mu.a, tau.a)
    b[j] ~ dnorm(mu.b, tau.b)
}
 
 mu.a ~ dnorm(0, 0.0001)
 tau.a <- 1 / (sigma.a * sigma.a)
 sigma.a ~ dunif(0, 1000)		

 mu.b ~ dnorm(0, 0.0001)
 tau.b <- 1 / (sigma.b * sigma.b)
 sigma.b ~ dunif(0, 1000)
 
 tau.y <- 1 / (sigma.y * sigma.y)
 sigma.y ~ dunif(0,1000)
}
",fill=TRUE)
sink()

# create SAD and calculate species detection probabilities
get_probs <- function(x) {
  vec <- rnbinom(x, true_df$nbinom_size, true_df$nbinom_prob)
  vec <- ifelse(vec==0,1,vec)
  vec <- vec/sum(vec)
  return(vec)
}

# Sampling simulation function
sampling_func <- function(vec) {
  temp <- rmultinom(vec[1], size = vec[2], prob = vec[-(1:2)])
  temp <- ifelse(temp>0,1,0)
  return(temp)
}

# insert sampling function to repeat the desired number of times
replicate_sampling <- function(vec) {
  sample <- replicate(replicates, { sampling_func(vec) }, simplify = FALSE)
  return(sample)
}

# Convert formatting for BHMS
replicate_format <- function(vec) {
  check <- apply(vec, 2, function(vec) {
    return(which(vec > 0))
  }, simplify = FALSE)
  seq <- c(); agg <- c()
  for(i in 1:length(check)) {
    seq[i] = i
    if(i == 1) {
      agg[i] <- length(check[[i]])
      tot <- check[[i]]
    } else {
      tot <- union(tot, check[[i]])
      agg[i]<- length(tot)
    }
  }
  return(data.frame(cbind(seq=seq, agg=agg)))
}

###############################
##### Simulation Settings #####
###############################

site_num <- 30 # Number of sites in nested domain
richness_vec <- c(50,100,150)
#richness_vec <- c(100)
#samples_vec <- c(50,75,100,125,150)
samples_vec <- c(100)
#sampling_event_vec <- c(5,10,15)
sampling_event_vec <- c(5,10,15)
#nbinom_size = c(.5,1,1.5,2,2.5)
#nbinom_size = c(.5,1,1.5,2,2.5,3,3.5,4)
nbinom_size = c(1)
nbinom_prob = .1
iter_num <- length(richness_vec)*length(samples_vec)*+length(sampling_event_vec)*length(nbinom_size)
replicates <- 100

# set initials
inits <- function(){list(mu.a = rnorm(1, 50, 0),
                         mu.b = rlnorm(1),
                         sigma.y = runif(1,1,10),
                         sigma.a = runif(1,1,10),
                         sigma.b = runif(1,1,10)) }
# params to save
params1 <- c("mu.a", "mu.b","sigma.y","sigma.b","sigma.a","a","b")
nt <- 5 # thining
nb <- 1000 # burnin
nc <- 3 # chains

# Set progress bar
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = iter_num,
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE,    # If TRUE, clears the bar when finish
                       width = 100)      # Width of the progress bar



dat_final <- list()
iter <- 0

for(i in richness_vec) {
  for(j in samples_vec) {
    for(z in sampling_event_vec) {
      for(t in nbinom_size) {
        true_df <- data.frame(cbind(
          true=rpois(site_num, i), 
          site=seq(1,site_num,1),
          sampling_events = rep(z, site_num),
          samples_per_event = rep(j, site_num),
          nbinom_size = t,
          nbinom_prob = .1
        ))
       
        probs <- lapply(true_df$true, get_probs)
        for(f in 1:site_num){
          probs[[f]] <- append(true_df$samples_per_event[f], probs[[f]])
          probs[[f]] <- append(true_df$sampling_events[f], probs[[f]])
        }
        
        replicated_samples <- lapply(probs, replicate_sampling)
        
        reformated <- lapply(replicated_samples, function(x) lapply(x, replicate_format))
        reformated2 <- list()
        for(f in 1:site_num){
          temp <- reformated[[f]]
          reformated2[[f]] <- bind_rows(temp, .id = "replicate")
        }
        dat <- bind_rows(reformated2, .id = "site")
        
        chao_est <- lapply(replicated_samples, function(x) lapply(x, Chao))
        reformated2 <- list()
        for(f in 1:site_num){
          temp <- chao_est[[f]]
          reformated2[[f]] <- bind_rows(temp, .id = "replicate")
        }
        chao_est_df <- bind_rows(reformated2, .id = "site")
        chao_est_df$method <- "Chao"
        
        jack_est <- lapply(replicated_samples, function(x) lapply(x, Jack1))
        reformated2 <- list()
        for(f in 1:site_num){
          temp <- jack_est[[f]]
          reformated2[[f]] <- bind_rows(temp, .id = "replicate")
        }
        jack_est_df <- bind_rows(reformated2, .id = "site")
        jack_est_df$method <- "Jack"
        
        Naive_est_df <- data.frame(cbind(
          site <- seq()
        ))
        
        Naive_est_df <- dat %>%
          group_by(site, replicate) %>%
          summarise(est = max(agg)) %>%
          mutate(sd = NA, upper95 = NA, lower95 = NA, method = "Naive")
        
        
        
        MM_est_df <- list()
        NE_est_df <- list()
        
        
        for(f in 1:replicates) {
          mod.df <- dat %>%
            filter(replicate == f)
          data <- list(y = as.numeric(mod.df$agg),
                       x = as.numeric(mod.df$seq),
                       site = as.numeric(mod.df$site),
                       n = nrow(mod.df),
                       J = length(unique(mod.df$site)))
          out_MM <- autojags(data = data,
                             inits = inits,
                             parameters.to.save = params1,
                             model.file = "New Sim/BMM.txt",
                             n.chains = nc,
                             n.thin = nt,
                             n.burnin = nb,
                             verbose = FALSE)
          out_NE <- autojags(data = data,
                             inits = inits,
                             parameters.to.save = params1,
                             model.file = "New Sim/BNE.txt",
                             n.chains = nc,
                             n.thin = nt,
                             n.burnin = nb,
                             verbose = FALSE)
          MM_est_df[[f]] <- data.frame(cbind(
            est = out_MM$mean$a,
            sd = out_MM$sd$a,
            site = seq(1,site_num,1),
            replicate = f,
            upper95 = out_MM$q97.5$a,
            lower95 = out_MM$q2.5$a,
            method = "BMM"
          ))
          NE_est_df[[f]] <- data.frame(cbind(
            est = out_NE$mean$a,
            sd = out_NE$sd$a,
            site = seq(1,site_num,1),
            replicate = f,
            upper95 = out_NE$q97.5$a,
            lower95 = out_NE$q2.5$a,
            method = "BNE"
          ))
        }
        MM_est_df <- bind_rows(MM_est_df)
        NE_est_df <- bind_rows(NE_est_df)
        
        dat_all_df <- rbind(MM_est_df, NE_est_df, chao_est_df, jack_est_df, Naive_est_df)
        dat_all_df[,1:6] <- apply(dat_all_df[,1:6], 2, as.numeric)
        dat_all_df <- merge(dat_all_df, true_df[,1:2], by = "site")
        dat_all_df$richness_sim <- i
        dat_all_df$samples_sim <- j
        dat_all_df$events_sim <- z
        dat_all_df$nbinom_prob_sim <- t
        
        iter <- iter + 1
        dat_final[[iter]] <- dat_all_df
        pb$tick()
      }
    }
  }
}

test <- bind_rows(dat_final)


saveRDS(dat_final, "data/Created/sim_data.RDS")

