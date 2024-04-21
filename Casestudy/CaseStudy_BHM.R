# Load libraries
library(tidyverse)
library(progress)
library(sf)
library(jagsUI)

# Read in prepped data
BHM_data <- readRDS("Data/Created/CS_dat_BHM.RDS")

# Unique huc2s for main loop
HUC2s <- unique(BHM_data$HUC_2)

# list to store data
all_est <- list()

# Create progress bar object
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = length(HUC2s),
                       complete = "=",   
                       incomplete = "-", 
                       current = ">",    
                       clear = FALSE,   
                       width = 100)

# Main loop
for(i in 1:length(HUC2s)){
  # filter by huc2
  mod_df <- BHM_data %>%
    filter(HUC_2 == HUC2s[i])
  # create factor levels for huc8s
  mod_df$HUC_8_factor <- as.numeric(as.factor(mod_df$HUC_8))
  
  # create factor lookup
  HUC_check <- mod_df %>%
    dplyr::select(HUC_8, HUC_8_factor)
  
  # drop duplicates and arrange
  HUC_check <- HUC_check[!duplicated(HUC_check$HUC_8),] %>%
    arrange(HUC_8_factor)
  
  # number of huc8s
  J <- length(unique(mod_df$HUC_8_factor))
  
  # package data for JAGS
  data <- list(y = as.numeric(mod_df$agg), 
               x = as.numeric(mod_df$seq), 
               site = as.numeric(mod_df$HUC_8_factor),
               n = nrow(mod_df),
               J = J)
  
  # Set initals
  inits <- function(){list(mu.a = rnorm(1, 30, 0), 
                           mu.b = rlnorm(1),
                           sigma.y = runif(1,1,10),
                           sigma.a = runif(1,1,10),
                           sigma.b = runif(1,1,10)) }
  
  # MCMC settings
  params1 <- c("mu.a", "mu.b","sigma.y","sigma.b","sigma.a","a","b")
  nt <- 5
  nb <- 2000
  nc <- 3

  # run jags
  out_MM <- autojags(data = data, 
                     inits = inits, 
                     parameters.to.save = params1, 
                     model.file = "Data/Created/BMM.txt",
                     n.chains = nc, 
                     n.thin = nt, 
                     n.burnin = nb)
  
  # Grab value to store
  maxes <- mod_df %>%
    group_by(HUC_8) %>%
    summarise(naive = max(agg),
              events = max(seq))
  
  # create row of estimates
  est <- data.frame(cbind(
    est = out_MM$mean$a,
    sd = out_MM$sd$a,
    upper95 = out_MM$q97.5$a,
    lower95 = out_MM$q2.5$a,
    HUC_8 = HUC_check$HUC_8,
    method = "BMM"
  ))
  
  # merge estimates and stored values
  est <- merge(est, maxes, by = "HUC_8")
  
  # Store in list
  all_est[[i]] <- est
  
  # tick progress bar
  pb$tick()
}

# Set huc names
names(all_est) <- HUC2s

# bind into dataframe
final <- bind_rows(all_est, .id = "HUC_2")

# Save output
saveRDS(final, "Data/Created/out_JAGS.RDS")