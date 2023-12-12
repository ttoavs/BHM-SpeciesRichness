library(tidyverse)

out_con <- readRDS("Data/Created/out_conventional.RDS")
out_jags <-readRDS("Data/Created/out_JAGS.RDS")

out_jags <- out_jags %>%
  select(-HUC_2)

out_con <- out_con %>%
  filter(HUC_8 %in% unique(out_jags$HUC_8))


dat <- rbind(out_con, out_jags)

dat <- dat %>%
  filter(method %in% c("BMM","Chao"))

dat <- dat %>%
  select(method, HUC_8, everything())

dat_wide <- dat %>%
  pivot_wider(names_from = method, values_from = c(3:6))

dat[,3:8] <- apply(dat[,3:8], 2, as.numeric)

dat_wide[,2:11] <- apply(dat_wide[,2:11], 2, as.numeric)

dat <- dat %>% 
  mutate(error_range = upper95-lower95)

dat_wide <- dat_wide %>%
  mutate(diff = est_BMM - est_Chao,
         direction = ifelse(diff>1,1,0))

mean(dat_wide$est_Chao)
mean(dat_wide$sd_Chao)
mean(dat_wide$est_BMM)
mean(dat_wide$sd_BMM)





plot1 <- ggplot(data = dat) +
  geom_point(aes(x=events, y=est, color = method, alpha = .15)) +
  facet_wrap(~method) +
  theme_classic(base_size = 16) +
  lims(x=c(0,500)) +
  scale_color_manual(values=c("steelblue2","red3")) +
  guides(alpha = "none") +
  labs(color = "Estimator", x= "", y = "Richness Estimate")
  

plot2 <- ggplot(data = dat) +
  geom_point(aes(x=as.numeric(naive), y=as.numeric(est), color = method, alpha = .15)) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  facet_wrap(~method) +
  theme_classic(base_size = 16) +
  scale_color_manual(values=c("steelblue1","red3")) +
  guides(alpha = "none") +
  labs(x="", y="Richness Estimate", color = "Estimator") +
  theme()

plot3 <- ggplot(data = dat) +
  geom_point(aes(x=as.numeric(events), y=as.numeric(error_range), color = method, alpha = .15)) +
  facet_wrap(~method) +
  theme_classic(base_size = 16) +
  lims(x=c(0,500)) +
  scale_color_manual(values=c("steelblue2","red3")) +
  labs(x="", y="Error Range") +
  theme(legend.position="none")


plot4 <- ggplot(data = dat) +
  geom_point(aes(x=as.numeric(naive), y=as.numeric(error_range), color = method, alpha = .15)) +
  facet_wrap(~method) +
  theme_classic(base_size = 16) +
  scale_color_manual(values=c("steelblue2","red3")) +
  labs(x="", y="Error Range") +
  theme(legend.position="none")

plot5 <- ggplot(data = dat_wide, aes(x=events, y=diff, color = factor(direction), alpha = .15)) +
  geom_point() +
  theme_classic(base_size = 16) +
  lims(x=c(0,500)) +
  scale_color_manual(values=c("red3","steelblue2")) +
  labs(x="Sample Size", y="Estimator Difference") +
  theme(legend.position="none")

plot6 <- ggplot(data = dat_wide, aes(x=naive, y=diff, color = factor(direction), alpha = .15)) +
  geom_point() +
  theme_classic(base_size = 16) +
  scale_color_manual(values=c("red3","steelblue2")) +
  labs(x="Naive Estimate", y="Estimator Difference") +
  theme(legend.position="none")

library(ggpubr)

ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 2, nrow = 3, common.legend = TRUE, legend="bottom")

panel_plot1 <- ggarrange(plot1, plot3, plot5, ncol = 1, nrow = 3, common.legend = TRUE, legend="right", labels = c("A","B","C"), hjust = -2.25)


panel_plot2 <- ggarrange(plot2, plot4, plot6, ncol = 1, nrow = 3, common.legend = TRUE, legend="right", labels = c("A","B","C"), hjust = -2.25)

ggsave(panel_plot1, device="tiff", dpi=300, filename = "Figures/fishscales_panel_plot_v2.tiff", height = 8.25, width = 7)

ggsave(panel_plot2, device="tiff", dpi=300, filename = "Figures/fishscales_panel_plot2_v2.tiff", height = 8.25, width = 7)




