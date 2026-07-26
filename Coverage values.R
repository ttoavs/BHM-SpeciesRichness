################################################################################
# Recover realized sample coverage (Cmax) for the coverage-based richness
# estimator used in the simulation study.
#
# This reproduces ONLY the data-generating and coverage-estimation portions of
# the original simulation script. The JAGS models, Chao2, Jackknife, and naive
# estimators are all omitted, since none of them are needed to recover the
# coverage values. Runtime should be minutes rather than hours.
#
# Output: a tidy data frame of the coverage level at which richness was
# estimated, for every site in every simulation scenario, plus summary tables
# suitable for reporting in the Methods or Results.
################################################################################

library(tidyverse)
library(iNEXT)

# The original simulation script called estimateD() but never loaded iNEXT,
# so it must have been attached earlier in the session. Loading it explicitly
# here.

set.seed(20260726)  # see note in the accompanying discussion about seeds


######################
##### Functions  #####
######################

# Unchanged from the original simulation script.

# Create SAD and calculate species detection probabilities
get_probs <- function(x) {
  vec <- rnbinom(x, true_df$nbinom_size, true_df$nbinom_prob)
  vec <- ifelse(vec == 0, 1, vec)
  vec <- vec / sum(vec)
  return(vec)
}

# Sampling simulation function
sampling_func <- function(vec) {
  temp <- rmultinom(vec[1], size = vec[2], prob = vec[-(1:2)])
  temp <- ifelse(temp > 0, 1, 0)
  return(temp)
}

# Repeat the sampling function the desired number of times
replicate_sampling <- function(vec) {
  sample <- replicate(replicates, { sampling_func(vec) }, simplify = FALSE)
  return(sample)
}


###############################
##### Simulation Settings #####
###############################

# Identical to the original script.

site_num <- 30
richness_vec <- c(50, 100, 150)
samples_vec <- c(100)
sampling_event_vec <- c(5, 10, 15)
nbinom_size <- c(1)
nbinom_prob <- 0.1
replicates <- 100


###############################
##### Coverage Extraction #####
###############################

coverage_out <- list()
iter <- 0

start <- Sys.time()

for (i in richness_vec) {
  for (j in samples_vec) {
    for (z in sampling_event_vec) {
      for (t in nbinom_size) {
        
        # Create df of true values
        true_df <- data.frame(cbind(
          true = rpois(site_num, i),
          site = seq(1, site_num, 1),
          sampling_events = rep(z, site_num),
          samples_per_event = rep(j, site_num),
          nbinom_size = t,
          nbinom_prob = 0.1
        ))
        
        # Calc probabilities
        probs <- lapply(true_df$true, get_probs)
        for (f in 1:site_num) {
          probs[[f]] <- append(true_df$samples_per_event[f], probs[[f]])
          probs[[f]] <- append(true_df$sampling_events[f], probs[[f]])
        }
        
        # Simulate replicated sampling
        replicated_samples <- lapply(probs, replicate_sampling)
        
        # Estimate coverage-based richness, one call per site.
        # nboot = 0 skips the bootstrap confidence intervals, which are not
        # needed here and are by far the slowest part of estimateD().
        for (f in 1:site_num) {
          
          est <- estimateD(replicated_samples[[f]],
                           datatype = "incidence_raw",
                           q = 0,
                           nboot = 0)
          
          # estimateD() standardizes every assemblage in a single call to a
          # common coverage level, so all replicates for a given site share
          # one value of Cmax. Taking the unique value guards against that
          # assumption silently breaking in a future iNEXT version.
          sc_vals <- unique(est$SC)
          
          coverage_out[[length(coverage_out) + 1]] <- data.frame(
            richness_sim   = i,
            samples_sim    = j,
            events_sim     = z,
            nbinom_size    = t,
            site           = f,
            true_richness  = true_df$true[f],
            coverage       = sc_vals[1],
            n_unique_sc    = length(sc_vals)
          )
        }
        
        iter <- iter + 1
        cat("Finished scenario", iter, "of",
            length(richness_vec) * length(samples_vec) *
              length(sampling_event_vec) * length(nbinom_size), "\n")
      }
    }
  }
}

end <- Sys.time()
runtime <- end - start
print(runtime)

coverage_df <- bind_rows(coverage_out)

hist(coverage_df$coverage)

# Sanity check: this should be all TRUE. If any row reports more than one
# unique coverage value, the assumption that estimateD standardizes to a
# single level per call no longer holds and the summaries below need revising.
stopifnot(all(coverage_df$n_unique_sc == 1))


###############################
##### Reportable Summaries ####
###############################

# Overall: the single sentence for the Methods section
overall_summary <- coverage_df %>%
  summarise(
    n_sites   = n(),
    mean_cov  = mean(coverage),
    sd_cov    = sd(coverage),
    min_cov   = min(coverage),
    max_cov   = max(coverage),
    median_cov = median(coverage)
  )

# By scenario: use this if coverage varies enough to warrant a Results table
scenario_summary <- coverage_df %>%
  group_by(richness_sim, events_sim) %>%
  summarise(
    n_sites  = n(),
    mean_cov = mean(coverage),
    min_cov  = min(coverage),
    max_cov  = max(coverage),
    .groups  = "drop"
  ) %>%
  arrange(richness_sim, events_sim)

print(overall_summary)
print(scenario_summary, n = Inf)

saveRDS(coverage_df, "Data/Created/coverage_values.RDS")
write.csv(scenario_summary, "Data/Created/coverage_by_scenario.csv",
          row.names = FALSE)