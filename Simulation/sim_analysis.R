library(tidyverse)

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


# dat_means <- dat_means %>%
#   filter(events_sim != 20)






plot.df <- dat_means

# plot.df <- plot.df %>%
#   filter(richness_sim == "100")

# plot.df <- plot.df %>%
#   filter(events_sim == "10")

# plot.df <- plot.df %>%
#  filter(samples_sim == "100")

# plot.df <- plot.df %>%
#   filter(nbinom_prob_sim == "1")



library(ggpubr)

dat_recover <- dat_means %>%
  group_by(method, richness_sim, samples_sim, events_sim, nbinom_prob_sim) %>%
  summarise(mean_recover = mean(percent_recover))


ggplot(data = dat_recover) +
  geom_bar(aes(y = mean_recover, x = method), stat = "identity") +
  facet_wrap(~richness_sim+events_sim) +
  theme_classic()


plot1 <- ggplot(data = plot.df, aes(x=richness_sim, y=diff, fill = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  facet_wrap(~events_sim) +
  scale_fill_manual(values=c("steelblue1", "steelblue", "red1", "red3","grey")) +
  theme_classic(base_size = 16) +
  labs(x="", y = "Difference from truth", fill = "Estimator") +
  theme(axis.text.x = element_blank())



plot2 <- ggplot(data = plot.df, aes(x=richness_sim, y=precision, fill = method)) +
  geom_boxplot() +
  facet_wrap(~events_sim) +
  scale_fill_manual(values=c("steelblue1", "steelblue", "red1", "red3")) +
  theme_classic(base_size = 16) +
  labs(x="", y="Error Range", fill = "Estimator") +
  theme(axis.text.x = element_blank(), legend.position="none", strip.text.x = element_blank())

plot3 <- ggplot(data = plot.df, aes(x=richness_sim, y=w.diff, fill = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  facet_wrap(~events_sim) +
  scale_fill_manual(values=c("steelblue1", "steelblue", "red1", "red3")) +
  theme_classic(base_size = 16) +
  labs(x="", y = "Weighted Difference", fill = "Estimator") +
  theme(legend.position="none", strip.text.x = element_blank(), axis.text.x = element_text(angle = 25, vjust = .6, hjust=.50)); plot3

panel_plot1 <- ggarrange(plot1, plot2, plot3, ncol = 1, common.legend = TRUE, legend="right")

ggsave(panel_plot1, device="tiff", dpi=300, filename = "Figures/sim_panel_plot1.tiff", height = 8, width = 7)







