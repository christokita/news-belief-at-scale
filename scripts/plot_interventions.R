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
path_to_interventions <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/interventions/'
intervention_dirs <- list.dirs(path_to_interventions)
intervention_dirs <- intervention_dirs[grepl('.*/reduce_', intervention_dirs)]
outpath <- 'output/interventions/'


####################
# Load data
####################

#################### Exposure data ####################

# Function to pre-first share dummy rows of pre-first share for plotting purposes
add_dummy_time_points <- function(exposure_data) {
  # (1) Create dummy rows that make two "non tweets" in at two time points leading up to the first real tweet
  n_articles <- length(unique(exposure_data$total_article_number)) #number of unique articles
  dummy_rows <- exposure_data %>% 
    distinct(total_article_number, intervention_time, visibility_reduction, sharing_reduction, replicate, simulation_type)
  dummy_row_time <- data.frame(total_article_number = rep(unique(exposure_data$total_article_number), each = 2),
                               time = rep(c(-2, 0), n_articles), 
                               new_exposed_users = rep(c(0, 0), n_articles),
                               cumulative_exposed = rep(c(0, 0), n_articles))
  dummy_rows <- merge(dummy_rows, dummy_row_time, by = "total_article_number")
  # (2) Join together
  exposure_data <- exposure_data %>% 
    rbind(dummy_rows) %>% 
    arrange(total_article_number, replicate, time)
  rm(dummy_rows, dummy_row_time)
  return(exposure_data)
}

# Load data
intervention_exposure <- data.frame()
for (dir in intervention_dirs) {
  intervention_files <- list.files(dir, full.names = TRUE)
  exposure_file <- intervention_files[grepl("_exposetime.csv", intervention_files)]
  for (file in exposure_file) {
    exposure <- read.csv(file)
    exposure <- add_dummy_time_points(exposure_data = exposure)
    intervention_exposure <- rbind(intervention_exposure, exposure)
    rm(exposure)
  }
  
}

# Create measure of relative exposure (relative to actual tweet data)
max_time_of_expsoure <-  max(intervention_exposure$time[intervention_exposure$new_exposed_users > 0]) # find where new users are no longer being exposed. 55hrs
intervention_exposure <- intervention_exposure %>% 
  group_by(total_article_number) %>% 
  mutate(relative_cumulative_exposed = cumulative_exposed / max(cumulative_exposed), 
         simulation_number = paste0(total_article_number, "-", replicate, "-", intervention_time, "-visibility", visibility_reduction, "-sharing", sharing_reduction),
         simulation_type = factor(simulation_type, levels = c("no intervention", "intervention")),
         intervention = paste0("t", intervention_time, "-visibility", visibility_reduction, "-sharing", sharing_reduction)) %>% 
  filter(time <= max_time_of_expsoure) 

  
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
# Plot exposure decrease by intervention
####################
# Filter to final exposure values (55 hours out from first share is enough)
exposure_decrease <- intervention_exposure %>% 
  filter(time == 55,
         replicate != -1) %>% 
  mutate(decrease_in_exposure = relative_cumulative_exposed - 1,
         intervention_amount = paste0("visibility", visibility_reduction, "-", "sharing", sharing_reduction)) %>% 
  group_by(intervention_amount, intervention_time) %>% 
  summarise(mean_exposure_decrease = mean(decrease_in_exposure))

plot_labels <- data.frame(intervention_amount = unique(exposure_decrease$intervention_amount),
                          intervention_label = c("Sharing friction (light)", "Sharing friction (heavy)", "Visibility reduction (light)", "Visibility reduction (heavy)")) %>% 
  mutate(intervention_label = gsub(" \\(", "\n\\(", intervention_label)) %>% 
  merge(exposure_decrease %>%  
          filter(intervention_time == 8) %>% 
          group_by(intervention_amount) %>% 
          select(intervention_amount, mean_exposure_decrease) %>% 
          rename(y_pos = mean_exposure_decrease), by = "intervention_amount")

# Plot
gg_exposure_decrease <- ggplot(exposure_decrease, aes(x = intervention_time, y = mean_exposure_decrease, color = intervention_amount)) +
  geom_segment(aes(xend = intervention_time, yend = 0), size = 0.6) +
  # geom_hline(aes(yintercept = 0), size = 0.5) +
  geom_point(size = 2.5,
             stroke = 0) +
  geom_text(data = plot_labels, aes(x = 12.5, y = y_pos, label = intervention_label), 
            fontface = "bold", 
            size = 2.5, 
            lineheight = 0.8,
            color = "black", 
            hjust = 1, 
            nudge_y = -0.2) +
  scale_x_continuous(breaks = seq(1, 12, 1),
                     limits = c(0, 13),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-1, 0, 0.2), 
                     limits = c(-0.8, 0), 
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#74c476", "#238b45", "#6baed6", "#2171b5"),
                      name = "Intervention",
                      labels = c("Sharing friction (light)",
                                 "Sharing friction (heavy)",
                                 "Visibility reduction (light)",
                                 "Visibility reduction (heavy)")) +
  xlab("Intervention time (hrs)") + 
  ylab("Mean reduction in user exposure to misinformation") +
  theme_ctokita() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA),
        legend.position = "none",
        strip.text = element_blank()) +
  facet_wrap(~intervention_amount,
             ncol = 2)
gg_exposure_decrease

ggsave(gg_exposure_decrease, filename = paste0(outpath, "interventions_exposurereduction.pdf"), width = 90, height = 90, units = "mm", dpi = 400)

####################
# Plot relative exposure by intervention
####################
# Filter to final exposure values (55 hours out from first share is enough)
exposure_decrease <- intervention_exposure %>% 
  filter(time == 55,
         replicate != -1) %>% 
  mutate(decrease_in_exposure = 1-relative_cumulative_exposed)

# Estimate mean reduction in exposure by intervention
treatments <- unique(exposure_decrease$intervention)
exposure_reduction_estimate <- data.frame()
for (treatment in treatments) {
  data <- exposure_decrease %>% 
    filter(intervention == treatment)
  est <- Bolstad::normnp(x = data$relative_cumulative_exposed, m.x = 0.5, s.x = 5, quiet = TRUE, plot = FALSE)
  df_est <- data.frame(intervention = unique(data$intervention),
                       intervention_time = unique(data$intervention_time),
                       visibility_reduction = unique(data$visibility_reduction),
                       sharing_reduction = unique(data$sharing_reduction),
                       est_mean = est$mean,
                       ci_95_low = est$quantileFun(0.005),
                       ci_95_high = est$quantileFun(0.995)) %>% 
    mutate(intervention_amount = paste0("visibility", visibility_reduction, "-", "sharing", sharing_reduction))
  exposure_reduction_estimate <- rbind(exposure_reduction_estimate, df_est)
  rm(df_est, est)
}

# Plot
gg_relative_exposure <- ggplot(exposure_reduction_estimate, aes(x = intervention_time, color = intervention_amount)) +
  # geom_point(data = exposure_decrease,
  #            aes(x = intervention_time, y = decrease_in_exposure),
  #            size = 0.8, stroke = 0, alpha = 0.3, position = position_jitter(width = 0.1)) +
  geom_errorbar(aes(ymin = ci_95_low, ymax = ci_95_high),
                width = 0, size = 0.5) +
  geom_point(aes(y = est_mean),
             size = 1) +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0.3, 1.0), expand = c(0,0)) +
  scale_colour_manual(values = c("#74c476", "#238b45", "#6baed6", "#2171b5"),
                      name = "Intervention",
                      labels = c("Sharing friction (light)",
                                 "Sharing friction (heavy)",
                                 "Visibility reduction (light)",
                                 "Visibility reduction (heavy)")) +
  xlab("Intervention time (hrs)") + 
  ylab("Relative users exposed") +
  theme_ctokita()
gg_relative_exposure
ggsave(gg_relative_exposure, filename = paste0(outpath, "relexposure_by_intervention.pdf"), width = 90, height = 45, units = "mm", dpi = 400)


####################
# Plot example intervention time series
####################
intervention_pal <- scales::viridis_pal(begin = 0, end = 0.9, direction = -1, option = "plasma")
intervention_pal <- intervention_pal(6)

gg_example_timeseries <- intervention_exposure %>% 
  filter(total_article_number == 28, 
         sharing_reduction == 0.75,
         visibility_reduction == 0,
         intervention_time %in% c(seq(0, 12, 2))) %>% 
  mutate(intervention_time = ifelse(simulation_type == "no intervention", "No intervention", paste(intervention_time, "hr."))) %>%
  mutate(intervention_time = factor(intervention_time, levels = c("No intervention", paste(seq(2, 12, 2) , "hr.")) )) %>%
  ggplot(., aes(x = time, y = cumulative_exposed, color = intervention_time, group = simulation_number, alpha = simulation_type)) +
  geom_line(size = 0.3) +
  scale_y_continuous(breaks = seq(0, 10000000, 2000000), 
                     limits = c(0, 10000000),
                     # expand = c(0, 0),
                     labels = scales::comma) +
  scale_x_continuous(breaks = seq(0, 48, 12)) +
  scale_color_manual(values = c("black", intervention_pal),
                     name = "Intervention time") +
  scale_alpha_manual(values = c(1, 0.2), 
                    guide = NULL) +
  xlab("Time since first article share (hrs)") +
  ylab("Total users exposed") +
  theme_ctokita() +
  theme(legend.position = "right",
        aspect.ratio = 0.5,
        legend.key.height = unit(0.5, 'mm'),
        legend.spacing = unit(0, 'mm'),
        legend.title = element_text(vjust = -1))
gg_example_timeseries
ggsave(gg_example_timeseries, filename = paste0(outpath, "exampleintervention_article28.pdf"), width = 120, height = 45, units = "mm", dpi = 400)

