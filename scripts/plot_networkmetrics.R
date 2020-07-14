########################################
#
# PLOT: Network metrics of articles
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
network_metric_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/networks/article_network_metrics.csv'
outpath <- 'output/networks/'

# Color palette
plot_color <- "#495867"


####################
# Load data 
####################
# Read in data
network_metrics <- read.csv(network_metric_path, header = TRUE)



############################## Ideology in networks ##############################

####################
# Plot: Ideological diversity by article veracity
####################
gg_ideodiversity <- network_metrics %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
           ggplot(., aes(x = article_fc_rating, y = ideology_sd)) +
  geom_point(size = 1, stroke = 0, alpha = 0.5,
             position = position_jitter(width = 0.03)) +
  scale_x_discrete(labels = c("Fake news", "Real news")) +
  ylab("Ideol. diversity of article tweeters") +
  xlab("") +
  theme_ctokita()
gg_ideodiversity
ggsave(gg_ideodiversity, filename = paste0(outpath, "ideodiversity_byveracity.png"), width = 45, height = 45, units = "mm")



############################## Network structure ##############################

####################
# Plot: Network density by article veracity
####################
gg_veracitydensity <- network_metrics %>% 
  filter(article_fc_rating %in% c("FM", "T")) %>% 
  ggplot(., aes(x = article_fc_rating, y = network_density)) +
  geom_point(size = 1, stroke = 0, alpha = 0.5,
             position = position_jitter(width = 0.03)) +
  scale_x_discrete(labels = c("Fake news", "Real news")) +
  ylab("Ideol. diversity of article tweeters") +
  xlab("") +
  theme_ctokita()
gg_veracitydensity
ggsave(gg_veracitydensity, filename = paste0(outpath, "networkdensity_byveracity.png"), width = 45, height = 45, units = "mm", dpi = 400)
