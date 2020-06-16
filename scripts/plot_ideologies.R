########################################
#
# PLOT: Distribution of ideologies of tweeters
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
tweeter_score_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/labeled_tweets.csv' #path to fitness cascade data
outpath <- 'output/ideology/'


############################## Ideological distribution of tweeters ##############################

####################
# Load data 
####################
# Read in data
tweeter_scores <- read.csv(tweeter_score_path, header = TRUE) %>% 
  mutate(article_ideology = article_con_feel - article_lib_feel)


####################
# Plot: Tweeter ideology distribution by article source lean
####################
# Plot colors
pal <- c("#d54c54", "#006195", "#C5CBD3")

# All source types
gg_fmtweeters <- ggplot(data = tweeter_scores, aes(x = ideology_score, group = article_type, fill = article_type)) +
  geom_histogram(alpha = 0.6, position = 'identity') +
  xlab("Tweeter ideology") +
  ylab("Log count") +
  scale_fill_manual(values = pal, name = "Article type", labels = c("False, Conservative source", "False, Liberal source", "False, Unclear source")) +
  scale_x_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  scale_y_continuous(trans = "log10", limits = c(1, 12000), expand = c(0,0)) +
  theme_ctokita() +
  theme(aspect.ratio = 0.5)
gg_fmtweeters
ggsave(gg_fmtweeters, file = paste0(outpath, "FMtweeters_ideologies_bysource.png"), width = 120, height = 45, units = "mm", dpi = 400)

# Just liberal and conservative sources
gg_fmtweeters_libcon <- tweeter_scores %>% 
  filter(article_type %in% c("false-conservative", "false-liberal")) %>% 
  ggplot(., aes(x = ideology_score, group = article_type, fill = article_type)) +
  geom_histogram(alpha = 0.6, position = 'identity') +
  xlab("Tweeter ideology") +
  ylab("Count") +
  scale_fill_manual(values = pal, name = "Article type", labels = c("False, Conservative source", "False, Liberal source")) +
  scale_x_continuous(limits = c(-4, 4)) +
  theme_ctokita() +
  theme(aspect.ratio = 0.5)
gg_fmtweeters_libcon
ggsave(gg_fmtweeters_libcon, file = paste0(outpath, "FMtweeters_ideologies_lib-con-sources.png"), width = 120, height = 45, units = "mm", dpi = 400)


####################
# Plot: Tweeter ideology by article lean
####################
# Plot colors
pal <- c("#d54c54", "#006195", "#9D69A3", "#C5CBD3")

gg_articlelean <- ggplot(data = tweeter_scores, aes(x = ideology_score, group = article_lean, fill = article_lean)) +
  geom_histogram(position = 'identity') +
  xlab("Tweeter ideology") +
  ylab("Count") +
  scale_x_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  scale_fill_manual(values = pal, name = "Article lean", labels = c("Conservative", "Liberal", "Neutral", "Unclear")) +
  facet_wrap(~article_lean, scales = 'free', ncol = 1) +
  theme_ctokita() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.2)
gg_articlelean
ggsave(gg_articlelean, file = paste0(outpath, "FMtweeters_ideologies_byarticlelean.png"), width = 90, height = 70, units = "mm", dpi = 400)


gg_articleideol <- ggplot(data = tweeter_scores, aes(x = article_ideology, y = ideology_score, color = article_lean)) +
  geom_hline(aes(yintercept = 0), size = 0.3, linetype = "dotted") +
  geom_vline(aes(xintercept = 0), size = 0.3, linetype = "dotted") +
  geom_point(alpha = 0.6, position = position_jitter(0.1)) +
  xlab("Article ideology") +
  ylab("Tweeter ideology") +
  scale_x_continuous(limits = c(-4, 4)) +
  scale_y_continuous(limits = c(-4, 4)) +
  theme_ctokita() 
gg_articleideol
