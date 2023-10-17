library(tidyverse)
library(gt)

all_sim <- readRDS("data/Created/sim_data.RDS")

dat <- bind_rows(all_sim)

dat_means <- dat %>%
  group_by(method, true, site, richness_sim, samples_sim, events_sim, nbinom_prob_sim) %>%
  mutate(recover = ifelse(true < lower95 | true > upper95, 0, 1)) %>%
  summarise(mean_est = mean(est), 
            mean_sd = mean(sd),
            mean_upper = mean(upper95),
            mean_lower = mean(lower95),
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


dat_means$method <- factor(dat_means$method, levels = c("BMM", "BNE", "Chao", "Jack","Naive"))
dat_means$richness_sim <- factor(dat_means$richness_sim, levels = c('Low','Medium','High'))
dat_means$events_sim <- factor(dat_means$events_sim, levels = c('Low','Medium','High'))

table <- dat_means %>%
  group_by(richness_sim, events_sim, method) %>%
  summarise(Diffence = mean(diff),
            Precision = mean(precision),
            `Weighted Difference` = mean(w.diff),
            `Mean Percent Recovered` = mean(percent_recover))

ROUND <- function(vec) {
  round(vec, 2)
}

table[,4:7] <- apply(table[,4:7], 2, ROUND)

colnames(table) <- c("Richness Scenario","Events Scenario","Estimator","Difference from truth","Error range","Weighted Difference","Percent recovered")

write.csv(table, "Data/Created/sim_table.csv", row.names = FALSE)
