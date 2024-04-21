# Included is the code to analyze the simulation results.
# Also included is code to prepare figures for the manuscript.

library(tidyverse)

# read in results
all_sim <- readRDS("data/Created/sim_data.RDS")

dat <- bind_rows(all_sim)

naive <- dat %>%
  ungroup() %>%
  select(true, site, replicate, richness_sim, samples_sim, events_sim, nbinom_prob_sim) %>%
  mutate(est = true, method = "Naive")

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

dat_means$method <- factor(dat_means$method, levels = c("BMM", "BNE", "Chao", "Jack","Naive"))
dat_means$site <- factor(dat_means$site)
dat_means$richness_sim <- factor(dat_means$richness_sim, labels = c("Low-Richness","Medium-Richness","High-Richness"))
dat_means$samples_sim <- factor(dat_means$samples_sim)
dat_means$events_sim <- factor(dat_means$events_sim, labels = c("Low-Sampling","Medium-Sampling","High-Sampling"))
dat_means$nbinom_prob_sim <- factor(dat_means$nbinom_prob_sim)


plot.df <- dat_means

 mean_diff_samplesize <- plot.df %>%
  filter(method != "Naive") %>%
  group_by(events_sim, richness_sim) %>%
  summarise(mean = mean(diff))

 ggplot(data = plot.df, aes(x=events_sim, y=diff)) +
   geom_boxplot() 

 mean_prec_samplesize <- dat_means %>%
   filter(method != "Naive") %>%
   group_by(events_sim, richness_sim) %>%
   summarise(mean = mean(as.numeric(precision)))
 
 ggplot(data = plot.df, aes(x=events_sim, y=precision)) +
   geom_boxplot() 


###############
## Bar plots ##
###############

library(ggpubr)

dat_recover <- dat_means %>%
  group_by(method, richness_sim, samples_sim, events_sim, nbinom_prob_sim) %>%
  summarise(mean_recover = mean(percent_recover))


ggplot(data = dat_recover) +
  geom_bar(aes(y = mean_recover, x = method), stat = "identity") +
  facet_wrap(~richness_sim+events_sim) +
  theme_classic()

###############
###############
###############

#####################
## Panel plot code ##
#####################







plot.df <- dat_means

levels(plot.df$richness_sim) <- c("Low","Medium","High")

plot1 <- ggplot(data = plot.df, aes(x=richness_sim, y=diff, fill = method, color = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  facet_wrap(~events_sim) +
  scale_color_manual(values=c("steelblue4", "steelblue", "red4", "red3", "grey50")) +
  scale_fill_manual(values=c("steelblue2", "steelblue1", "red2", "red1", "grey")) +
  theme_classic(base_size = 16) +
  labs(x="", y = "Difference from truth", fill = "Estimator") +
  guides(color="none") +
  theme(axis.text.x = element_blank())

plot1b <- ggplot(data = plot.df, aes(x=richness_sim, y=diff, fill = method, color = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  facet_wrap(~events_sim) +
  scale_color_manual(values=c("#4B4B61", "#49585F", "#564444", "#8C594A", "grey50")) +
  scale_fill_manual(values=c("#0000FF", "#00B7FF", "#A50000", "#FF9D81", "grey")) +
  theme_classic(base_size = 16) +
  labs(x="", y = "Difference from truth", fill = "Estimator") +
  guides(color="none") +
  theme(axis.text.x = element_blank())

#9AEAFF

plot2 <- ggplot(data = plot.df, aes(x=richness_sim, y=precision, fill = method, color = method)) +
  geom_boxplot() +
  facet_wrap(~events_sim) +
  scale_color_manual(values=c("steelblue4", "steelblue", "red4", "red3"), ) +
  scale_fill_manual(values=c("steelblue2", "steelblue1", "red2", "red1")) +
  theme_classic(base_size = 16) +
  labs(x="", y="Error Range", fill = "Estimator") +
  theme(axis.text.x = element_blank(), legend.position="none", strip.text.x = element_blank())

plot2b <- ggplot(data = plot.df, aes(x=richness_sim, y=precision, fill = method, color = method)) +
  geom_boxplot() +
  facet_wrap(~events_sim) +
  scale_color_manual(values=c("#4B4B61", "#49585F", "#564444", "#8C594A")) +
  scale_fill_manual(values=c("#0000FF", "#00B7FF", "#A50000", "#FF9D81")) +
  theme_classic(base_size = 16) +
  labs(x="", y="Error Range", fill = "Estimator") +
  theme(axis.text.x = element_blank(), legend.position="none", strip.text.x = element_blank())

plot3 <- ggplot(data = plot.df, aes(x=richness_sim, y=w.diff, fill = method, color = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  facet_wrap(~events_sim) +
  scale_color_manual(values=c("steelblue4", "steelblue", "red4", "red3")) +
  scale_fill_manual(values=c("steelblue2", "steelblue1", "red2", "red1")) +
  theme_classic(base_size = 16) +
  labs(x="Richness", y = "Weighted Difference", fill = "Estimator") +
  theme(legend.position="none", strip.text.x = element_blank())

plot3b <- ggplot(data = plot.df, aes(x=richness_sim, y=w.diff, fill = method, color = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  facet_wrap(~events_sim) +
  scale_color_manual(values=c("#4B4B61", "#49585F", "#564444", "#8C594A")) +
  scale_fill_manual(values=c("#0000FF", "#00B7FF", "#A50000", "#FF9D81")) +
  theme_classic(base_size = 16) +
  labs(x="Richness", y = "Weighted Difference", fill = "Estimator") +
  theme(legend.position="none", strip.text.x = element_blank())

panel_plot1 <- ggarrange(plot1, plot2, plot3, ncol = 1, common.legend = TRUE, legend="right", labels = c("A","B","C"), hjust = -2.25, font.label = list(face = "bold", color = "black"), align = "v"); panel_plot1

panel_plot1b <- ggarrange(plot1b, plot2b, plot3b, ncol = 1, common.legend = TRUE, legend="right", labels = c("A","B","C"), hjust = -2.25, font.label = list(face = "bold", color = "black"), align = "v"); panel_plot1

#ggsave(panel_plot1, device="tiff", dpi=300, filename = "Figures/sim_panel_plot1_V2.tiff", height = 8.25, width = 7.25)

#ggsave(panel_plot1b, device="tiff", dpi=300, filename = "Figures/sim_panel_plot1_V3.tiff", height = 8.25, width = 7.25)

#####################
#####################
#####################





