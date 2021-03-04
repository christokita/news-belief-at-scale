########################################
#
# PLOT: Crowdsouced interventions for decreasing fake news spread
#
########################################

####################
# Load packages
####################
library(ggplot2)
library(dplyr)
library(brms)
library(RColorBrewer)
library(viridis)
source("scripts/_plot_themes/theme_ctokita.R")


####################
# Paramters for analysis: gpaths to data, paths for output, and filename
####################
# Paths to files/directories
path_to_interventions <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/interventions/'
crowd_intervention_dirs <- list.dirs(path_to_interventions)
crowd_intervention_dirs <- crowd_intervention_dirs[grepl('crowdsourced_reduce', crowd_intervention_dirs)]
outpath <- 'output/interventions/'

# Palette
veracity_pal <- c("#F1A208",  "#858585", "#000004")

####################
# Load data
####################

#################### Exposure data ####################

# Function to pre-first share dummy rows of pre-first share for plotting purposes
add_dummy_time_points <- function(exposure_data) {
  # (1) Create dummy rows that make two "non tweets" in at two time points leading up to the first real tweet
  n_articles <- length(unique(exposure_data$total_article_number)) #number of unique articles
  dummy_rows <- exposure_data %>% 
    distinct(total_article_number, intervention_time, visibility_reduction, sharing_reduction, replicate, simulation_type, fc_rule)
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

# Load intervention data
intervention_exposure <- data.frame()
for (dir in crowd_intervention_dirs) {
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
         simulation_type = factor(simulation_type, levels = c("baseline", "no intervention", "intervention"))) %>% 
  filter(time <= max_time_of_expsoure) 

# Add article fc rating (professional fact checkers)
article_info <- read.csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/tweets_labeled.csv') %>% 
  select(total_article_number, article_fc_rating, source_type) %>% 
  distinct()

intervention_exposure <- merge(intervention_exposure, article_info, by = "total_article_number") %>% 
  mutate(article_fc_rating = gsub("No Mode!|CND", "Borderline", article_fc_rating), #group the non-T/F ratings
         article_fc_rating = factor(article_fc_rating, levels = c("FM", "Borderline", "T")))

####################
# Plot rate of being labeled false by crowd, broken out by crowd rule
####################
gg_labelingrate <- intervention_exposure %>% 
  # Prep data
  filter(replicate != -1) %>% 
  distinct(total_article_number, simulation_type, replicate, fc_rule, article_fc_rating) %>% 
  mutate(fc_rule = factor(fc_rule, levels = c("mean", "mode", "median", "majority", "unanimity"))) %>% 
  group_by(article_fc_rating, fc_rule) %>% 
  summarise(labeled_false_rate = sum(simulation_type == "intervention") / length(simulation_type)) %>% 
  # Plot
  ggplot(., aes(x = article_fc_rating, y = labeled_false_rate, color = article_fc_rating)) +
  geom_segment(aes(xend = article_fc_rating, yend = 0), size = 0.6) +
  geom_point(size = 2.5, stroke = 0) +
  scale_x_discrete(labels = c("F", "Unc", "T")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), 
                     limits = c(0, 0.701),
                     expand = c(0, 0)) +
  scale_color_manual(values = veracity_pal) +
  xlab("Article rating by fact checkers") +
  ylab("Rate flagged \"false\" by crowd") +
  theme_ctokita() +
  theme(aspect.ratio = NULL,
        legend.position = "none") +
  facet_grid(~fc_rule)
  
gg_labelingrate
ggsave(gg_labelingrate, filename = paste0(outpath, "crowd_false_rating_rate.pdf"), width = 90, height = 45, units = "mm")

####################
# Plot relative exposure for true and fake news
####################
# Filter to final exposure values (55 hours out from first share is enough)
exposure_decrease <- intervention_exposure %>% 
  filter(replicate != -1,
         time == max_time_of_expsoure,
         article_fc_rating %in% c("FM", "T")) %>% 
  mutate(decrease_in_exposure = 1-relative_cumulative_exposed) 

# exposure_decrease$relative_cumulative_exposed[exposure_decrease$relative_cumulative_exposed == 1] <- 0.5

# Plot raw
ggplot(exposure_decrease, aes(x = article_fc_rating, y = relative_cumulative_exposed)) +
  geom_point(stroke = 0, alpha = 0.4, position = position_jitter(width = 0.1, height = 0.01)) +
  theme_ctokita() +
  facet_grid(~fc_rule)

# Estimate mean reduction in exposure by intervention
# prior <- set_prior("Beta(1,1)", class = "b")
# prior <- get_prior(relative_cumulative_exposed ~ 0 + article_fc_rating + fc_rule, data = exposure_decrease, family = gaussian())
# blm_exposure <- brm(relative_cumulative_exposed ~ 0 + article_fc_rating + fc_rule, 
#                     data = exposure_decrease, 
#                     prior = prior,
#                     family = gaussian(), 
#                     warmup = 500, 
#                     iter = 1500, 
#                     chains = 1)
# 
# exposure_posterior <- posterior_samples(blm_exposure)
# posterior_summary(blm_exposure)

# Plot mean reduction
gg_relative_exposure <- exposure_decrease %>% 
  group_by(fc_rule, article_fc_rating) %>% 
  mutate(fc_rule = factor(fc_rule, levels = c("mean", "mode", "median", "majority", "unanimity"))) %>% 
  summarise(decrease_in_exposure = mean(decrease_in_exposure)) %>% 
  # Plot
  ggplot(., aes(x = article_fc_rating, y = decrease_in_exposure, color = article_fc_rating)) +
  geom_segment(aes(xend = article_fc_rating, yend = 0), size = 0.6) +
  geom_point(size = 2.5,
             stroke = 0) +
  scale_x_discrete(labels = c("F", "T")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), 
                     limits = c(0, 0.4001),
                     expand = c(0, 0)) +
  scale_color_manual(values = veracity_pal[c(1, 3)]) +
  xlab("Article rating by fact checkers") +
  ylab("Mean reduction in user expsoure") +
  theme_ctokita() +
  theme(aspect.ratio = NULL,
        legend.position = "none") +
  facet_grid(~fc_rule)
gg_relative_exposure

ggsave(gg_relative_exposure, filename = paste0(outpath, "crowd_exposure_reduction_byrule.pdf"), width = 90, height = 45, units = "mm")

