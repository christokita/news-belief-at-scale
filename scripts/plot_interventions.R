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
    distinct(total_article_number, intervention_time, visibility_reduction, belief_reduction, sharing_reduction, replicate, simulation_type)
  dummy_row_time <- data.frame(total_article_number = rep(unique(exposure_data$total_article_number), each = 2),
                               time = rep(c(-2, 0), n_articles), 
                               new_exposed_users = rep(c(0, 0), n_articles),
                               new_believing_users = rep(c(0, 0), n_articles),
                               cumulative_exposed = rep(c(0, 0), n_articles),
                               cumulative_believing = rep(c(0, 0), n_articles))
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

# Add belief reduction column akin to scripts (USE THIS ONCE THEN DELETE this after all data is updated)
# intervention_exposure <- data.frame()
# for (dir in intervention_dirs) {
#   intervention_files <- list.files(dir, full.names = TRUE)
#   exposure_file <- intervention_files[grepl("_exposetime.csv", intervention_files)]
#   belief_impact <- as.numeric (gsub(".*belief([.0-9]+).*", "\\1", dir, perl = T))
#   for (file in exposure_file) {
#     exposure <- read.csv(file)
#     exposure <- exposure %>% 
#       mutate(belief_reduction = belief_impact) %>% 
#       relocate(belief_reduction, .after = sharing_reduction)
#     write.csv(exposure, file = file, row.names = FALSE)
#   }
#   
# }

# Create measure of relative exposure (relative to actual tweet data)
max_time_of_expsoure <-  max(intervention_exposure$time[intervention_exposure$new_exposed_users > 0]) # find where new users are no longer being exposed. 55hrs
intervention_exposure <- intervention_exposure %>% 
  group_by(total_article_number) %>% 
  mutate(max_article_exposure = max(cumulative_exposed),
         max_article_believing = max(cumulative_believing),
         relative_cumulative_exposed = cumulative_exposed / max_article_exposure, 
         relative_cumulative_believing = cumulative_believing / max_article_believing, 
         change_in_cumlative_exposed = (cumulative_exposed - max_article_exposure) / max_article_exposure,
         change_in_cumlative_believing = (cumulative_believing - max_article_believing) / max_article_believing,
         simulation_number = paste0(total_article_number, "-", replicate, "-", intervention_time, "-visibility", visibility_reduction, "-sharing", sharing_reduction),
         simulation_type = factor(simulation_type, levels = c("no intervention", "intervention")),
         intervention = paste0("t", intervention_time, "-visibility", visibility_reduction, "-sharing", sharing_reduction)) %>% 
  filter(time <= max_time_of_expsoure) %>% 
  filter(total_article_number > 10) #we don't use first 10 articles in analysis


####################
# Prep data on relative effect of interventions
####################
# Filter to final exposure/belief values for total impact
relative_effect <- intervention_exposure %>% 
  filter(time == max_time_of_expsoure,
         replicate != -1) %>% 
  mutate(intervention_amount = paste0("visibility", visibility_reduction, "-", "sharing", sharing_reduction)) %>% 
  group_by(intervention_amount, intervention_time) %>% 
  summarise(mean_exposure = mean(relative_cumulative_exposed),
            median_exosure = median(relative_cumulative_exposed),
            exposure_5 = quantile(relative_cumulative_exposed, c(0.05)),
            exposure_10 = quantile(relative_cumulative_exposed, c(0.1)),
            exposure_25 = quantile(relative_cumulative_exposed, c(0.25)),
            exposure_75 = quantile(relative_cumulative_exposed, c(0.75)),
            exposure_90 = quantile(relative_cumulative_exposed, c(0.90)),
            exposure_95 = quantile(relative_cumulative_exposed, c(0.95)),
            mean_belief = mean(relative_cumulative_believing),
            median_belief = median(relative_cumulative_believing),
            belief_5 = quantile(relative_cumulative_believing, c(0.05)),
            belief_10 = quantile(relative_cumulative_believing, c(0.1)),
            belief_25 = quantile(relative_cumulative_believing, c(0.25)),
            belief_75 = quantile(relative_cumulative_believing, c(0.75)),
            belief_90 = quantile(relative_cumulative_believing, c(0.90)),
            belief_95 = quantile(relative_cumulative_believing, c(0.95)))

write.csv(relative_effect, file = paste0(path_to_interventions, "interventions_relative_effect.csv"), row.names = FALSE) #write to file for use in crowd-sourced intervention analysis


# Prep plot labeling
relative_effect <- read.csv(file = paste0(path_to_interventions, "interventions_relative_effect.csv"))
relative_effect <- relative_effect %>% 
  mutate(intervention_type = ifelse(intervention_amount == "visibility0-sharing0.25", "Fact-check labeling",
                                    ifelse(intervention_amount == "visibility0-sharing0.75", "Sharing friction",
                                           ifelse(intervention_amount == "visibility0.25-sharing0", "Visibility reduction (light)",
                                                  ifelse(intervention_amount == "visibility0.75-sharing0", "Visibility reduction (heavy)", ""))))) %>% 
  mutate(intervention_type = factor(intervention_type, levels = c("Fact-check labeling",
                                                                  "Sharing friction",
                                                                  "Visibility reduction (light)",
                                                                  "Visibility reduction (heavy)")))

####################
# Plot relative exposure by intervention
####################
# Bar plot of mean
gg_exposure_decrease <- ggplot(relative_effect, aes(x = intervention_time, y = mean_exposure, fill = intervention_type)) +
  geom_bar(stat = "identity", 
           color = NA,
           width = 0.8) +
  scale_x_continuous(breaks = seq(1, 12, 1),
                     limits = c(0, 13),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1), 
                     expand = c(0,0)) +
  scale_fill_manual(values = c("#74c476", "#238b45", "#6baed6", "#2171b5"),
                      name = "Intervention",
                      labels = c("Fact-check labeling",
                                 "Sharing friction",
                                 "Visibility reduction (light)",
                                 "Visibility reduction (heavy)")) +
  xlab(expression( paste("Intervention delay ", italic(t[int]), " (hr)") )) + 
  ylab("Relative user exposure to misinformation") +
  theme_ctokita() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA),
        legend.position = "none",
        strip.text = element_text(vjust = -0.5)) +
  facet_wrap(~intervention_type,
             ncol = 2)
gg_exposure_decrease

ggsave(gg_exposure_decrease, filename = paste0(outpath, "interventions_relativeexposure.pdf"), width = 85, height = 100, units = "mm", dpi = 400)


# Point plot of mean and quantiles
gg_exposure_decrease_point <- ggplot(relative_effect, aes(x = intervention_time, y = mean_exposure, color = intervention_type, fill = intervention_type)) +
  geom_ribbon(aes(ymin = exposure_5,
                  ymax = exposure_95),
              color = NA,
              alpha = 0.25) +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 12, 2),
                     limits = c(0, 12),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1), 
                     expand = c(0,0)) +
  scale_color_manual(values = c("#74c476", "#238b45", "#6baed6", "#2171b5"),
                    name = "Intervention",
                    labels = c("Fact-check labeling",
                               "Sharing friction",
                               "Visibility reduction (light)",
                               "Visibility reduction (heavy)")) +
  scale_fill_manual(values = c("#74c476", "#238b45", "#6baed6", "#2171b5"),
                    name = "Intervention",
                    labels = c("Fact-check labeling",
                               "Sharing friction",
                               "Visibility reduction (light)",
                               "Visibility reduction (heavy)")) +
  xlab(expression( paste("Intervention delay ", italic(t[int]), " (hr)") )) + 
  ylab("Relative user exposure to misinformation") +
  theme_ctokita() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA),
        panel.spacing = unit(0.25, "cm"),
        legend.position = "none",
        strip.text = element_text(vjust = -0.5)) +
  facet_wrap(~intervention_type,
             ncol = 2)
gg_exposure_decrease_point

ggsave(gg_exposure_decrease_point, filename = paste0(outpath, "interventions_relativeexposure_detailed.pdf"), width = 85, height = 100, units = "mm", dpi = 400)


####################
# Plot relative belief by intervention
####################
# Bar plot of mean
gg_belief_decrease <- ggplot(relative_effect, aes(x = intervention_time, y = mean_belief, fill = intervention_type)) +
  geom_bar(stat = "identity", 
           color = NA,
           width = 0.8) +
  scale_x_continuous(breaks = seq(1, 12, 1),
                     limits = c(0, 13),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1), 
                     expand = c(0,0)) +
  scale_fill_manual(values = c("#74c476", "#238b45", "#6baed6", "#2171b5"),
                    name = "Intervention",
                    labels = c("Fact-check labeling",
                               "Sharing friction",
                               "Visibility reduction (light)",
                               "Visibility reduction (heavy)")) +
  xlab(expression( paste("Intervention delay ", italic(t[int]), " (hr)") )) + 
  ylab("Relative user belief of misinformation") +
  theme_ctokita() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA),
        legend.position = "none",
        strip.text = element_text(vjust = -0.5)) +
  facet_wrap(~intervention_type,
             ncol = 2)
gg_belief_decrease

ggsave(gg_belief_decrease, filename = paste0(outpath, "interventions_relativebelief.pdf"), width = 85, height = 100, units = "mm", dpi = 400)

# Point plot of mean and quantiles
gg_belief_decrease_point <- ggplot(relative_effect, aes(x = intervention_time, y = mean_belief, color = intervention_type, fill = intervention_type)) +
  geom_ribbon(aes(ymin = belief_5,
                  ymax = belief_95),
              color = NA,
              alpha = 0.25) +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 12, 2),
                     limits = c(0, 12),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1), 
                     expand = c(0,0)) +
  scale_color_manual(values = c("#74c476", "#238b45", "#6baed6", "#2171b5"),
                     name = "Intervention",
                     labels = c("Fact-check labeling",
                                "Sharing friction",
                                "Visibility reduction (light)",
                                "Visibility reduction (heavy)")) +
  scale_fill_manual(values = c("#74c476", "#238b45", "#6baed6", "#2171b5"),
                    name = "Intervention",
                    labels = c("Fact-check labeling",
                               "Sharing friction",
                               "Visibility reduction (light)",
                               "Visibility reduction (heavy)")) +
  xlab(expression( paste("Intervention delay ", italic(t[int]), " (hr)") )) + 
  ylab("Relative user belief of misinformation") +
  theme_ctokita() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA),
        legend.position = "none",
        strip.text = element_text(vjust = -0.5)) +
  facet_wrap(~intervention_type,
             ncol = 2)
gg_belief_decrease_point

ggsave(gg_exposure_decrease_point, filename = paste0(outpath, "interventions_relativeexposure_detailed.pdf"), width = 85, height = 100, units = "mm", dpi = 400)


####################
# Plot example exposure intervention time series
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
  scale_y_continuous(breaks = seq(0, 20000000, 2000000), 
                     limits = c(0, 12000000),
                     # expand = c(0, 0),
                     labels = scales::comma) +
  scale_x_continuous(breaks = seq(0, 48, 6),
                     limits = c(-0.1, 48)) +
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

