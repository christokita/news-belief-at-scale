########################################
#
# PLOT: Ideology check - self-reported ideology (categorical) vs estimated Twitter ideology (continuous)
#
########################################

####################
# Load packages
####################
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
source("scripts/_plot_themes/theme_ctokita.R")


####################
# Paramters for analysis: paths to data, paths for output, and filename
####################
ideology_scores_path <- '/Volumes/CKT-DATA/news-belief-at-scale/data/ideology_check/selfreported_vs_pabloscore.csv' 
outpath <- 'output/ideology_basis/'

# Color palette
plot_color <- "#495867"
ideol_pal_small <- rev(brewer.pal(5, "RdBu"))
ideol_pal_small[3] <- "#e0e0e0"

ideol_pal_large <- rev(brewer.pal(7, "RdBu"))
ideol_pal_large[4] <- "#e0e0e0"



####################
# Load data 
####################
# Read in self-reported vs. pablo score data
ideologies <- read.csv(ideology_scores_path, header = TRUE) %>% 
  mutate(ideology_category = ideo5)


####################
# PLOT: scatter of two ideology scores for each user
####################
# Calculate mean by ideology category
ideology_means <- ideologies %>% 
  group_by(ideology_category) %>% 
  summarise(mean_pablo_score = mean(pablo_score),
            sd_pablo_score = sd(pablo_score),
            se_pablo_score = sd(pablo_score) / sqrt(length(pablo_score)))

# Create composite mean for "somewhat X" categories
somewhat_liberal <- ideology_means %>% 
  filter(ideology_category %in% c("Liberal", "Moderate")) %>% 
  summarise(mean_pablo_score = mean(mean_pablo_score)) %>% 
  mutate(ideology_category = "Somewhat liberal",
         sd_pablo_score = NA,
         se_pablo_score = NA)

somewhat_conservative <- ideology_means %>% 
  filter(ideology_category %in% c("Conservative", "Moderate")) %>% 
  summarise(mean_pablo_score = mean(mean_pablo_score)) %>% 
  mutate(ideology_category = "Somewhat conservative",
         sd_pablo_score = NA,
         se_pablo_score = NA)

ideology_means <- rbind(ideology_means, somewhat_liberal)
ideology_means <- rbind(ideology_means, somewhat_conservative) %>% 
  mutate(ideology_category = factor(ideology_category, levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative", "Not sure"))) %>% 
  arrange(ideology_category)

ideologies <- ideologies %>% 
  mutate(ideology_category = factor(ideology_category, levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative", "Not sure")))

# Plot
gg_ideology_comparison <- ggplot() +
  geom_hline(yintercept = 0, linetype = 'dotted', size = 0.3) +
  # Have to double plot group means to force factors to plot in order
  geom_point(data = ideology_means, aes(x = ideology_category, y = mean_pablo_score, fill = ideology_category),
             size = 3, 
             shape = 21,
             color = 'white') +  
  geom_point(data = ideologies, aes(x = ideology_category, y = pablo_score, color = ideology_category),
             size = 1.5, 
             stroke = 0, 
             alpha = 0.15,
             position = position_jitter(width = 0.1, height = 0)) +
  geom_point(data = ideology_means, aes(x = ideology_category, y = mean_pablo_score, fill = ideology_category),
             size = 3,
             shape = 21,
             color = 'white') +
  xlab('Self-reported ideology') +
  ylab('Inferred ideology from Twitter') +
  scale_x_discrete(labels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative", "Not sure")) +
  scale_y_continuous(breaks = seq(-4, 4, 1), limits = c(-4, 4), expand = c(0, 0)) +
  scale_color_manual(values = c(ideol_pal_small, 'black')) +
  scale_fill_manual(values = c(ideol_pal_large, 'black')) +
  theme_ctokita() +
  theme(legend.position = 'none',
        aspect.ratio = NULL,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.03))

gg_ideology_comparison
ggsave(gg_ideology_comparison, filename = paste0(outpath, 'empirical_ideology_means.pdf'), width = 90, height = 70, units = 'mm', dpi = 400)

####################
# Calculate bins using these empirical categories
####################
ideology_means %>% 
  filter(ideology_category != "Not sure") %>% 
  mutate(bin_edge_upper = (mean_pablo_score + lead(mean_pablo_score)) / 2,
         bin_edge_lower = (mean_pablo_score + lag(mean_pablo_score)) / 2,
         bin_size = bin_edge_upper - bin_edge_lower)
