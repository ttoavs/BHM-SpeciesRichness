library(tidyverse)
library(gt)

# Code to prepare table for publication


all_sim <- readRDS("data/Created/sim_revision_data.RDS")

dat <- bind_rows(all_sim)



dat_means <- dat %>%
  group_by(method, true, site, richness_sim, samples_sim, events_sim, nbinom_prob_sim) %>%
  mutate(recover = ifelse(true < lower95 | true > upper95, 0, 1)) %>%
  summarise(mean_est = mean(est), 
            mean_sd = mean(sd),
            mean_upper = mean(upper95),
            mean_lower = mean(lower95),
            rmse = mean((est-true)^2)^(1/2),
            percent_recover = sum(recover)/100)

dat_means <- dat_means %>%
  mutate(diff = mean_est - true,
         precision = mean_upper-mean_lower,
         w.diff = diff * (1/mean_sd))

dat_means <- dat_means %>%
  filter(events_sim != 20)


dat_means$richness_sim[dat_means$richness_sim == 150] <- "High"
dat_means$richness_sim[dat_means$richness_sim == 100] <- "Medium"
dat_means$richness_sim[dat_means$richness_sim == 50] <- "Low"

dat_means$events_sim[dat_means$events_sim == 15] <- "High"
dat_means$events_sim[dat_means$events_sim == 10] <- "Medium"
dat_means$events_sim[dat_means$events_sim == 5] <- "Low"


dat_means$method <- factor(dat_means$method, levels = c("BMM", "BNE", "Chao", "Jack","CMAX","Naive"))
dat_means$richness_sim <- factor(dat_means$richness_sim, levels = c('Low','Medium','High'))
dat_means$events_sim <- factor(dat_means$events_sim, levels = c('Low','Medium','High'))

table <- dat_means %>%
  group_by(richness_sim, events_sim, method) %>%
  summarise(Diffence = mean(diff),
            Precision = mean(precision),
            `Weighted Difference` = mean(w.diff),
            `Mean Percent Recovered` = mean(percent_recover),
            RMSE = mean(rmse))

ROUND <- function(vec) {
  round(vec, 2)
}

table[,4:8] <- apply(table[,4:8], 2, ROUND)

colnames(table) <- c("Richness Scenario","Events Scenario","Estimator","Difference from truth","Error range","Weighted Difference","Percent recovered","RMSE")

#write.csv(table, "Data/Created/sim_table_revision.csv", row.names = FALSE)



# Metric evaluation

# Accuracy
accuracy_by_sim <- table %>%
  group_by(`Richness Scenario`, `Events Scenario`) %>%
  slice(which.min(abs(`Difference from truth`))) %>%
  ungroup()

# table with the count of which estimators were closest to the truth 
table(accuracy_by_sim$Estimator, accuracy_by_sim$`Richness Scenario`, accuracy_by_sim$`Events Scenario`)


weighted_diff_by_sim <- table %>%
  group_by(`Richness Scenario`, `Events Scenario`) %>%
  slice(which.min(abs(`Weighted Difference`))) %>%
  ungroup()


rmse_by_sim <- table %>%
  group_by(`Richness Scenario`, `Events Scenario`) %>%
  slice(which.min(abs(RMSE))) %>%
  ungroup()







