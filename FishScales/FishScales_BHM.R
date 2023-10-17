# Load libraries ----
library(tidyverse)
library(progress)
library(sf)
library(jagsUI)

# Read in the data ----
fish <- readRDS("Data/Downloaded/dataset_complete_2023-08-16.rds")


# Data manipulation ----
df <- fish %>%
  drop_na(Date, X, Y, HUC_8) %>% # filter NA's
  filter(Year >= 2010) %>% # filter date 
  mutate(sampleId = paste(COMID, Date, Prim_Gear, sep = "_")) %>% # create sampling ID
  group_by(HUC_8) %>% # Group data by HUC8 
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  arrange(Date) %>% # arranging in order of Data
  mutate(eventNum = match(sampleId, unique(sampleId))) # Count by unique sampleID in each HU8

huc_vec <- unique(df$HUC_8)


filter_func <- function(sample){
  dat <- df %>%
    filter(HUC_8 %in% sample)
  dat$species <- factor(dat$Scientific_Name) 
  dat <- dat %>% 
    group_by(species, eventNum, .drop = FALSE) %>% # Group data
    tally() %>% # Tally up counts
    pivot_wider(names_from = species, values_from = n) %>% # Pivot so there is a column for each species
    arrange(eventNum) # arrange by event sequence
  
  dat <- dat[,-1] # drop column 1
  dat <- data.frame(ifelse(is.na(dat), 0, 1)) # change abundance counts to incidence value
  return(dat)
}

manipulate_accumulation <- function(matrix) {
  agg <- c()
  if(ncol(matrix) == 1){
    return(data.frame(cbind(seq=seq(1,nrow(matrix),1), agg=1)))
  }
  else{
    for(i in 1:nrow(matrix)){
      agg[i] <- length(which(apply(matrix[1:i,], 2, sum)>0))
    }
    return(data.frame(cbind(seq=seq(1,nrow(matrix),1), agg=agg)))
  }
}

# test <- filter_func("10030105")
# test2 <- manipulate_accumulation(test)

sampling_list <- lapply(huc_vec, filter_func) # Create sampling matrices
BHM_data <- lapply(sampling_list, manipulate_accumulation) # Create seq data
names(BHM_data) <- huc_vec # Add huc8 data
BHM_data <- bind_rows(BHM_data, .id = "HUC_8") # Create df


filter_df <- BHM_data %>%
  group_by(HUC_8) %>%
  summarise(max = max(seq)) %>%
  filter(max >= 5)


BHM_data <- BHM_data %>%
  filter(HUC_8 %in% unique(filter_df$HUC_8)) %>%
  mutate(HUC_2 = substr(HUC_8, 1, 2))

HUC2s <- unique(BHM_data$HUC_2)

all_est <- list()

pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = length(HUC2s),
                       complete = "=",   
                       incomplete = "-", 
                       current = ">",    
                       clear = FALSE,   
                       width = 100)

for(i in 1:length(HUC2s)){
  mod_df <- BHM_data %>%
    filter(HUC_2 == HUC2s[i])
  mod_df$HUC_8_factor <- as.numeric(as.factor(mod_df$HUC_8))
  
  
  HUC_check <- mod_df %>%
    dplyr::select(HUC_8, HUC_8_factor)
  
  HUC_check <- HUC_check[!duplicated(HUC_check$HUC_8),] %>%
    arrange(HUC_8_factor)
  
  J <- length(unique(mod_df$HUC_8_factor))
  
  data <- list(y = as.numeric(mod_df$agg), 
               x = as.numeric(mod_df$seq), 
               site = as.numeric(mod_df$HUC_8_factor),
               n = nrow(mod_df),
               J = J)
  
  inits <- function(){list(mu.a = rnorm(1, 30, 0), 
                           mu.b = rlnorm(1),
                           sigma.y = runif(1,1,10),
                           sigma.a = runif(1,1,10),
                           sigma.b = runif(1,1,10)) }
  
  params1 <- c("mu.a", "mu.b","sigma.y","sigma.b","sigma.a","a","b")
  nt <- 5
  nb <- 2000
  nc <- 3

  out_MM <- autojags(data = data, 
                     inits = inits, 
                     parameters.to.save = params1, 
                     model.file = "New Sim/BMM.txt",
                     n.chains = nc, 
                     n.thin = nt, 
                     n.burnin = nb)
  
  maxes <- mod_df %>%
    group_by(HUC_8) %>%
    summarise(naive = max(agg),
              events = max(seq))
  
  
  est <- data.frame(cbind(
    est = out_MM$mean$a,
    sd = out_MM$sd$a,
    upper95 = out_MM$q97.5$a,
    lower95 = out_MM$q2.5$a,
    HUC_8 = HUC_check$HUC_8,
    method = "BMM"
  ))
  
  est <- merge(est, maxes, by = "HUC_8")
  
  
  all_est[[i]] <- est
  
  pb$tick()
}


names(all_est) <- HUC2s

final <- bind_rows(all_est, .id = "HUC_2")


saveRDS(final, "Data/Created/out_JAGS.RDS")