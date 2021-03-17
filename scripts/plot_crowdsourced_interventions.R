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
veracity_pal <- c("#F18805",  "#8BAAAD", "#495867")

####################
# Load data
####################
# Load intervention data
intervention_exposure <- data.frame()
for (dir in crowd_intervention_dirs) {
  intervention_files <- list.files(dir, full.names = TRUE)
  exposure_file <- intervention_files[grepl("_exposetime.csv", intervention_files)]
  for (file in exposure_file) {
    exposure <- read.csv(file)
    intervention_exposure <- rbind(intervention_exposure, exposure)
    rm(exposure)
  }
  
}

# Create measure of relative exposure (relative to actual tweet data)
max_time_of_expsoure <-  max(intervention_exposure$time[intervention_exposure$new_exposed_users > 0]) # find where new users are no longer being exposed. 55hrs
intervention_exposure <- intervention_exposure %>% 
  group_by(total_article_number) %>% 
  mutate(max_article_exposure = max(cumulative_exposed),
         relative_cumulative_exposed = cumulative_exposed / max_article_exposure, 
         change_in_cumlative_exposed = (cumulative_exposed - max_article_exposure) / max_article_exposure,
         simulation_number = paste0(total_article_number, "-", replicate, "-", intervention_time, "-visibility", visibility_reduction, "-sharing", sharing_reduction),
         simulation_type = factor(simulation_type, levels = c("baseline", "no intervention", "intervention"))) %>% 
  filter(time == max_time_of_expsoure) 

# Add article fc rating (professional fact checkers)
article_info <- read.csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/tweets_labeled.csv') %>% 
  select(total_article_number, article_fc_rating, source_type) %>% 
  distinct()

intervention_exposure <- merge(intervention_exposure, article_info, by = "total_article_number") %>% 
  mutate(article_fc_rating = gsub("No Mode!|CND", "Unclear", article_fc_rating), #group the non-T/F ratings
         article_fc_rating = factor(article_fc_rating, levels = c("FM", "Unclear", "T")))

####################
# Plot rate of being labeled false by crowd, broken out by crowd rule
####################
# FM and T articles only
gg_labelingrate <- intervention_exposure %>% 
  # Prep data
  filter(replicate != -1,
         article_fc_rating %in% c("FM", "T")) %>% 
  distinct(total_article_number, simulation_type, replicate, fc_rule, article_fc_rating) %>% 
  mutate(fc_rule = factor(fc_rule, levels = c("mean", "mode", "median", "majority", "unanimity"))) %>% 
  group_by(article_fc_rating, fc_rule) %>% 
  summarise(labeled_false_rate = sum(simulation_type == "intervention") / length(simulation_type)) %>% 
  # Plot
  ggplot(., aes(x = article_fc_rating, y = labeled_false_rate, color = article_fc_rating)) +
  geom_segment(aes(xend = article_fc_rating, yend = 0), size = 0.6) +
  geom_point(size = 2.5, stroke = 0) +
  scale_x_discrete(labels = c("F", "T")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), 
                     limits = c(0, 0.701),
                     expand = c(0, 0)) +
  scale_color_manual(values = veracity_pal[c(1, 3)]) +
  xlab("Article rating by professional fact-checkers") +
  ylab("Rate labeled \"false\" by crowd") +
  theme_ctokita() +
  theme(aspect.ratio = NULL,
        legend.position = "none") +
  facet_grid(~fc_rule)
  
gg_labelingrate
ggsave(gg_labelingrate, filename = paste0(outpath, "crowd_false_rating_rate.pdf"), width = 60, height = 45, units = "mm")

# All articles
gg_labelingrate_all <- intervention_exposure %>% 
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
  xlab("Article rating by professional fact-checkers") +
  ylab("Rate labeled \"false\" by crowd") +
  theme_ctokita() +
  theme(aspect.ratio = NULL,
        legend.position = "none") +
  facet_grid(~fc_rule)

gg_labelingrate_all
ggsave(gg_labelingrate_all, filename = paste0(outpath, "crowd_false_rating_rate_alltypes.pdf"), width = 90, height = 45, units = "mm")

# Broken out by source type
gg_labelingrate_sourcetype <- intervention_exposure %>% 
  # Prep data
  filter(replicate != -1) %>% 
  distinct(total_article_number, simulation_type, replicate, fc_rule, article_fc_rating, source_type) %>% 
  mutate(fc_rule = factor(fc_rule, levels = c("mean", "mode", "median", "majority", "unanimity")),
         source_type = gsub("$", " source", source_type)) %>% 
  group_by(article_fc_rating, fc_rule, source_type) %>% 
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
  xlab("Article rating by professional fact-checkers") +
  ylab("Rate flagged \"false\" by crowd") +
  theme_ctokita() +
  theme(aspect.ratio = NULL,
        legend.position = "none") +
  facet_grid(source_type~fc_rule, scales = "free_x")

gg_labelingrate_sourcetype
ggsave(gg_labelingrate_sourcetype, filename = paste0(outpath, "crowd_false_rating_rate_bysourcetype.pdf"), width = 90, height = 90, units = "mm")


####################
# Plot relative exposure for true and fake news
####################
# Filter to final exposure values (55 hours out from first share is enough)
exposure_decrease <- intervention_exposure %>% 
  filter(replicate != -1,
         time == max_time_of_expsoure,
         article_fc_rating %in% c("FM", "T"))

# Plot raw
ggplot(exposure_decrease, aes(x = article_fc_rating, y = relative_cumulative_exposed)) +
  geom_point(stroke = 0, alpha = 0.4, position = position_jitter(width = 0.1, height = 0.01)) +
  theme_ctokita() +
  facet_grid(~fc_rule)

# Plot mean reduction
gg_relative_exposure <- exposure_decrease %>% 
  group_by(fc_rule, article_fc_rating) %>% 
  mutate(fc_rule = factor(fc_rule, levels = c("mean", "mode", "median", "majority", "unanimity"))) %>% 
  summarise(decrease_in_exposure = mean(change_in_cumlative_exposed)) %>% 
  # Plot
  ggplot(., aes(x = article_fc_rating, y = decrease_in_exposure, color = article_fc_rating)) +
  geom_segment(aes(xend = article_fc_rating, yend = 0), size = 0.6) +
  geom_point(size = 2.5,
             stroke = 0) +
  scale_x_discrete(labels = c("F", "T")) +
  scale_y_continuous(breaks = seq(-1, 1, 0.1), 
                     limits = c(-0.4001, 0),
                     expand = c(0, 0)) +
  scale_color_manual(values = veracity_pal[c(1, 3)]) +
  xlab("Article rating by professional fact-checkers") +
  ylab("Mean reduction in user exposure") +
  theme_ctokita() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA),
        legend.position = "none",
        # strip.text = element_blank(),
        aspect.ratio = NULL) +
  facet_grid(~fc_rule)
gg_relative_exposure

ggsave(gg_relative_exposure, filename = paste0(outpath, "crowd_exposure_reduction_byrule.pdf"), width = 90, height = 45, units = "mm")



#################### Compare interventions based on crowd-sourced and professional fact-checking ####################

####################
# Load data
####################
professional_intervention_exposure <- read.csv(paste0(path_to_interventions, "interventions_mean_exposure_reduction.csv")) %>% 
  filter(intervention_amount == "visibility0.75-sharing0") #only focus on same intervention type


####################
# Prep crowd sourced data
####################
crowd_intervention_exposure <- exposure_decrease %>% 
  group_by(fc_rule, article_fc_rating) %>% 
  summarise(decrease_in_exposure = mean(change_in_cumlative_exposed)) %>% 
  filter(fc_rule == "mean",
         article_fc_rating == "FM")


####################
# Plot comparison
####################
# Calcualte relative performance
comparison_exposure <- professional_intervention_exposure %>% 
  mutate(relative_speed = intervention_time - 1, #crowd_sourced fact-checking is implemented at t_int = 1
         crowd_sourced_exposure = crowd_intervention_exposure$decrease_in_exposure,
         relative_performance_professional = (crowd_sourced_exposure - mean_exposure_decrease) / mean_exposure_decrease)

# Plot
gg_comparison <- ggplot(comparison_exposure, aes(x = intervention_time, y = relative_performance_professional)) +
  geom_hline(aes(yintercept = 0), size = 0.3, linetype = "dotted") +
  geom_point(size = 2, stroke = 0) +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(breaks = seq(0, 12, 4),
                     limits = c(0, 12),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-0.6, 1.2, 0.3),
                     limits = c(-0.6, 1.2),
                     expand = c(0, 0)) +
  xlab("Prof. fact-checker\nturnaround time (hr)") +
  ylab("Relative performance of\ncrowd-sourced interventions") +
  theme_ctokita() +
  theme(aspect.ratio = NULL)
gg_comparison

ggsave(gg_comparison, filename = paste0(outpath, "comparison_crowd_professional_interventions.pdf"), width = 30, height = 45, units = "mm")
