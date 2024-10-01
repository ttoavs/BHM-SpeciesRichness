# Load libraries
library(tidyverse)
library(progress)
library(sf)
library(iNEXT)

# Read in prepped data
df <- readRDS("Data/Created/CS_dat_con.RDS")


# Estimator functions ----
# Chao
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
# Jackknife
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


# Set up

# DF to store the estimates
estimates.df <- list() 
# Vector of HUC8s
huc_vec <- unique(df$HUC_8)
# Progress bar for main function
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = length(huc_vec),
                       complete = "=",   
                       incomplete = "-", 
                       current = ">",    
                       clear = FALSE,   
                       width = 100)

# Main loop
for(i in 1:length(huc_vec)){ # Loop with each iteration corresponding to each HUC8
  dat <- df %>%
    filter(HUC_8 %in% huc_vec[i])
  # Filter data set by HUC8
  
  dat$species <- factor(dat$Scientific_Name)
  # Create factorized column for species
  
  dat <- dat %>% 
    group_by(species, eventNum, .drop = FALSE) %>% # Group data
    tally() %>% # Tally up counts
    pivot_wider(names_from = species, values_from = n) %>% # Pivot so there is a column for each species
    arrange(eventNum) # arrange by event sequence
  
  dat <- dat[,-1] # drop column 1
  dat <- data.frame(ifelse(is.na(dat), 0, 1)) # change abundance counts to incidence value
  
  chao_est <- data.frame(t(Chao(t(dat))))
  chao_est$method <- "Chao"
  chao_est$HUC_8 <- huc_vec[i]
  chao_est$naive <- ncol(dat)
  chao_est$events <- nrow(dat)
  
  jack_est <- data.frame(t(Jack1(t(dat))))
  jack_est$method <- "Jack"
  jack_est$HUC_8 <- huc_vec[i]
  jack_est$naive <- ncol(dat)
  jack_est$events <- nrow(dat)

  
  CMAX_est <- estimateD(t(dat), datatype = "incidence_raw", q = 0) %>%
    rename(est = qD, upper95 = qD.UCL, lower95 = qD.LCL) %>%
    mutate(method = "CMAX", HUC_8 = huc_vec[i], naive = ncol(dat), events = nrow(dat), sd = (upper95-lower95)/3.92) %>%
    dplyr::select(est, sd, upper95, lower95, method, HUC_8, naive, events)

  
  jack_est <- data.frame(t(Jack1(t(dat))))
  jack_est$method <- "Jack"
  jack_est$HUC_8 <- huc_vec[i]
  jack_est$naive <- ncol(dat)
  jack_est$events <- nrow(dat)
  
  est <- rbind(chao_est, jack_est, CMAX_est)
  
  estimates.df[[i]] <- est # Store estimates
  
  pb$tick() # tick over progress bar
}

final <- bind_rows(estimates.df)

#saveRDS(final, "Data/Created/out_conventional_revision.RDS")
