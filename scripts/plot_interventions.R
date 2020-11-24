########################################
#
# PLOT: Interventions for decreasing fake news spread
#
########################################

####################
# Load packages
####################
library(ggplot2)
library(dplyr)
library(Bolstad)
source("scripts/_plot_themes/theme_ctokita.R")


####################
# Paramters for analysis: gpaths to data, paths for output, and filename
####################
# Paths to files/directories
intervention_dirs <- c('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/interventions/reduceviz0.5_t6/',
                       '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/interventions/reduceviz0.5_t5/',
                       '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/interventions/reduceviz0.5_t4/',
                       '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/interventions/reduceviz0.5_t3/',
                       '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/interventions/reduceviz0.5_t2/')
outpath <- 'output/interventions/'


####################
# Load data
####################

#################### Exposure data ####################

# Function to pre-first share dummy rows of pre-first share for plotting purposes
add_dummy_time_points <- function(exposure_data, intervention_time, intervention_amount) {
  # (1) Create dummy rows that make two "non tweets" in at two time points leading up to the first real tweet
  n_articles <- length(unique(exposure_data$total_article_number)) #number of unique articles
  dummy_rows <- exposure_data %>% 
    distinct(total_article_number, replicate, simulation_type)
  dummy_row_time <- data.frame(total_article_number = rep(unique(exposure_data$total_article_number), each = 2),
                               time = rep(c(-2, 0), n_articles), 
                               new_exposed_users = rep(c(0, 0), n_articles),
                               cumulative_exposed = rep(c(0, 0), n_articles))
  dummy_rows <- merge(dummy_rows, dummy_row_time, by = "total_article_number")
  # (2) Join together
  exposure_data <- exposure_data %>% 
    rbind(dummy_rows) %>% 
    arrange(total_article_number, replicate, time) %>% 
    mutate(intervention_time = intervention_time, 
           intervention_amount = intervention_amount)
  rm(dummy_rows, dummy_row_time)
  return(exposure_data)
}

# Load data
intervention_exposure <- data.frame()
for (dir in intervention_dirs) {
  
  intervention_t <- as.numeric( gsub(".*_t([0-9]+)/", "\\1", dir, perl = TRUE) )
  intervention_am <- as.numeric( gsub(".*reduceviz([.0-9]+)_.*", "\\1", dir, perl = TRUE) )
  intervention_files <- list.files(dir, full.names = TRUE)
  exposure_file <- intervention_files[grepl("_exposetime.csv", intervention_files)]
  for (file in exposure_file) {
    exposure <- read.csv(file)
    exposure <- add_dummy_time_points(exposure_data = exposure, intervention_time = intervention_t, intervention_amount = intervention_am)
    intervention_exposure <- rbind(intervention_exposure, exposure)
    rm(exposure)
  }
  
}

# Create measure of relative exposure (relative to actual tweet data)
intervention_exposure <- intervention_exposure %>% 
  group_by(total_article_number) %>% 
  mutate(relative_cumulative_exposed = cumulative_exposed / max(cumulative_exposed), 
         simulation_number = paste0(total_article_number, "-", replicate, "-", intervention_time, "-", intervention_amount),
         simulation_type = factor(simulation_type, levels = c("no intervention", "intervention")),
         intervention = paste0(intervention_time, "-", intervention_amount)) %>% 
  filter(time <= 55)

  
####################
# Plot raw counts of cumulative exposure
####################
gg_exposed_raw <- ggplot(intervention_exposure) +
  geom_rect(xmin = 6, xmax = 55, ymin = 0, ymax = 1, color = NA, fill = "#f1eef6", alpha = 0.3) +
  geom_line(aes(x = time, y = cumulative_exposed, 
                group = simulation_number,
                color = simulation_type, 
                alpha = simulation_type)) +
  scale_color_manual(values = c("black", "#f03b20")) +
  scale_alpha_manual(values = c(1, 0.2)) +
  # scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  ylab("Relative number of users exposed") +
  xlab("Time since first article share (hrs)") +
  theme_ctokita() +
  facet_wrap(~total_article_number, scales = "free_y")
gg_exposed_raw


####################
# Plot relative exposure by intervention
####################
# Filter to data of interest
exposure_decrease <-  intervention_exposure %>% 
  filter(time == 55,
         replicate != -1) %>% 
  mutate(decrease_in_exposure = 1-relative_cumulative_exposed)

# Estimate mean reduction in exposure by intervention
treatments <- unique(intervention_exposure$intervention)
exposure_reduction_estimate <- data.frame()
for (treatment in treatments) {
  data <- exposure_decrease %>% 
    filter(intervention == treatment)
  est <- Bolstad::normnp(x = data$relative_cumulative_exposed, m.x = 0.5, s.x = 5, quiet = TRUE, plot = FALSE)
  df_est <- data.frame(intervention_time = unique(data$intervention_time),
                       intervention_amount = unique(data$intervention_amount),
                       est_mean = est$mean,
                       ci_95_low = est$quantileFun(0.025),
                       ci_95_high = est$quantileFun(0.975))
  exposure_reduction_estimate <- rbind(exposure_reduction_estimate, df_est)
  rm(df_est, est)
}

# Plot
gg_relative_exposure <- ggplot(exposure_reduction_estimate, aes(x = intervention_time, color = intervention_time)) +
  # geom_point(data = exposure_decrease,
  #            aes(x = intervention_time, y = decrease_in_exposure),
  #            size = 0.8, stroke = 0, alpha = 0.3, position = position_jitter(width = 0.1)) +
  geom_errorbar(aes(ymin = ci_95_low, ymax = ci_95_high),
                width = 0, size = 0.5) +
  geom_point(aes(y = est_mean),
             size = 2) +
  scale_x_continuous(breaks = seq(2, 6, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.05), limits = c(0.6, 0.8)) +
  scale_colour_viridis_c(option = "plasma", direction = -1, end = 0.9, guide = NULL) +
  xlab("Intervention time (hrs)") + 
  ylab("Relative user exposure") +
  theme_ctokita()
gg_relative_exposure
ggsave(gg_relative_exposure, filename = paste0(outpath, "relexposure_0.5reduction_by_interventiontime.png"), width = 45, height = 45, units = "mm", dpi = 400)


####################
# Plot example intervention time series
####################
intervention_pal <- scales::viridis_pal(begin = 0, end = 0.9, direction = 1, option = "plasma")
intervention_pal <- intervention_pal(5)

gg_example_timeseries <- intervention_exposure %>% 
  filter(total_article_number == 28) %>% 
  mutate(intervention_time = ifelse(simulation_type == "no intervention", "No intervention", paste(intervention_time, "hr."))) %>% 
  mutate(intervention_time = factor(intervention_time, levels = c("No intervention", "6 hr.", "5 hr.", "4 hr.", "3 hr.", "2 hr."))) %>% 
  ggplot(., aes(x = time, y = cumulative_exposed, color = intervention_time, group = simulation_number, alpha = simulation_type)) +
  geom_line(size = 0.5) +
  scale_y_continuous(breaks = seq(0, 10000000, 2000000), 
                     limits = c(0, 10000000),
                     labels = scales::comma) +
  scale_x_continuous(breaks = seq(0, 48, 12)) +
  scale_color_manual(values = c("black", intervention_pal),
                     name = "Intervention time,\n50% vis. reduction") +
  scale_alpha_manual(values = c(1, 0.3), 
                    guide = NULL) +
  xlab("Time since first article share (hrs)") +
  ylab("Total users exposed") +
  theme_ctokita()
gg_example_timeseries
ggsave(gg_example_timeseries, filename = paste0(outpath, "exampleintervention_article28.png"), width = 110, height = 90, units = "mm", dpi = 400)

