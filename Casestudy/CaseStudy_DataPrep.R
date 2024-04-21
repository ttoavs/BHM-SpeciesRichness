# Load libraries
library(tidyverse)
library(progress)
library(sf)

# Read in the downloaded data
fish <- readRDS("Data/Downloaded/dataset_complete_2023-08-16.rds")


# General data prep
df <- fish %>%
  drop_na(Date, X, Y, HUC_8) %>% # filter NA's
  filter(Year >= 2010) %>% # filter date 
  mutate(sampleId = paste(COMID, Date, Prim_Gear, sep = "_")) %>% # create sampling ID
  group_by(HUC_8) %>% # Group data by HUC8 
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  arrange(Date) %>% # arranging in order of Data
  mutate(eventNum = match(sampleId, unique(sampleId))) # Count by unique sampleID in each HU8


# Data prep for BHMs
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


saveRDS(BHM_data, "Data/Created/CS_dat_BHM.RDS")


# Data prep for conventional estimators
dat <- df %>%
  select(HUC_8, Scientific_Name, eventNum)

  
saveRDS(dat, "Data/Created/CS_dat_con.RDS")