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
tweeter_score_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/tweets_labeled.csv' 
outpath <- 'output/ideology/'

# For labeling facet plots
label_veracity <- c("T" = "True news", 
                    "FM" = "Fake news",
                    "CND" = "CND",
                    "No mode!" = "No mode!")

label_lean <- c("C" = "Conservative",
                "L" = "Liberal",
                "N" = "Neutral",
                "U" = "Unclear")

# Color palette
plot_color <- "#495867"
ideol_pal <- rev(brewer.pal(5, "RdBu"))
ideol_pal[3] <- "#e0e0e0"
source_pal <- c(ideol_pal[c(5,1)], "#9D69A3", "#C5CBD3")


####################
# Load data 
####################
# Read in data
tweeter_scores <- read.csv(tweeter_score_path, header = TRUE) %>% 
  mutate(article_ideology = article_con_feel - article_lib_feel) %>% 
  filter(!is.na(total_article_number))



############################## Ideological distribution of tweeters ##############################

####################
# Plot: Tweeter ideology distribution by article source lean
####################
# All source types
gg_fmtweeters <- tweeter_scores %>% 
  filter(article_fc_rating == "FM",
         !is.na(user_ideology)) %>% 
  ggplot(., aes(x = user_ideology, group = source_lean, fill = source_lean)) +
  geom_histogram(alpha = 0.6, position = 'identity', binwidth = 0.25) +
  xlab("Tweeter ideology") +
  ylab("Log number of tweets") +
  scale_fill_manual(values = source_pal[c(1,2,4)], name = "Article type", labels = c("False, Conservative source", "False, Liberal source", "False, Unclear source")) +
  scale_x_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  scale_y_continuous(trans = "log10", limits = c(1, 12000), expand = c(0,0)) +
  theme_ctokita() +
  theme(aspect.ratio = 0.5)
gg_fmtweeters
ggsave(gg_fmtweeters, file = paste0(outpath, "FMtweeters_ideologies_bysource.png"), width = 120, height = 45, units = "mm", dpi = 400)

# Just liberal and conservative sources
gg_fmtweeters_libcon <- tweeter_scores %>% 
  filter(source_lean %in% c("C", "L"),
         article_fc_rating == "FM") %>% 
  ggplot(., aes(x = user_ideology, group = source_lean, fill = source_lean)) +
  geom_histogram(alpha = 0.6, position = 'identity', binwidth = 0.25) +
  xlab("Tweeter ideology") +
  ylab("Number of tweets") +
  scale_fill_manual(values = source_pal[c(1,2,4)], name = "Article type", labels = c("False, Conservative source", "False, Liberal source")) +
  scale_x_continuous(limits = c(-4, 4)) +
  theme_ctokita() +
  theme(aspect.ratio = 0.5)
gg_fmtweeters_libcon
ggsave(gg_fmtweeters_libcon, file = paste0(outpath, "FMtweeters_ideologies_lib-con-sources.png"), width = 120, height = 45, units = "mm", dpi = 400)


####################
# Plot: Tweeter ideology by article lean
####################

###### All articles #####
# Distribution of user ideology by article lean
gg_articlelean <- ggplot(data = tweeter_scores, aes(x = user_ideology, group = article_lean, fill = article_lean)) +
  geom_histogram(position = 'identity', binwidth = 0.25) +
  xlab("Tweeter ideology") +
  ylab("Number of tweets") +
  scale_x_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  scale_fill_manual(values = source_pal, name = "Article lean", labels = c("Conservative", "Liberal", "Neutral", "Unclear")) +
  facet_wrap(~article_lean, scales = 'free', ncol = 1) +
  theme_ctokita() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.2)
gg_articlelean
ggsave(gg_articlelean, file = paste0(outpath, "ideologies_byarticlelean.png"), width = 90, height = 70, units = "mm", dpi = 400)

# Tweeter ideology vs article ideology
gg_articleideol <- tweeter_scores %>% 
  filter(article_fc_rating %in% c("FM", "T"),
         !is.na(user_ideology)) %>% 
  ggplot(., aes(x = user_ideology, y = article_ideology, color = article_lean)) +
  geom_hline(aes(yintercept = 0), size = 0.3, linetype = "dotted") +
  geom_vline(aes(xintercept = 0), size = 0.3, linetype = "dotted") +
  geom_point(size = 0.5,
             alpha = 0.4, 
             stroke = 0, 
             position = position_jitter(0.1)) +
  xlab("Tweeter ideology") +
  ylab("Article ideology") +
  scale_color_manual(values = source_pal, name = "Article lean", labels = c("Conservative", "Liberal", "Neutral", "Unclear")) +
  scale_x_continuous(breaks = seq(-6, 6, 2), limits = c(-4.7, 4.7)) +
  scale_y_continuous(breaks = seq(-6, 6, 2), limits = c(-4.6, 4.6)) +
  theme_ctokita() +
  facet_wrap(~article_fc_rating, scales = "free")
gg_articleideol
ggsave(gg_articleideol, file = paste0(outpath, "ideology_vs_articleideology.png"), width = 100, height = 45, units = "mm", dpi = 400)


###### Just FM articles #####
gg_articlelean <- tweeter_scores %>% 
  filter(article_fc_rating == "FM") %>% 
  ggplot(., aes(x = user_ideology, group = article_lean, fill = article_lean)) +
  geom_histogram(position = 'identity', binwidth = 0.25) +
  xlab("Tweeter ideology") +
  ylab("Number of tweets") +
  scale_x_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  scale_fill_manual(values = source_pal, name = "Article lean", labels = c("Conservative", "Liberal", "Neutral", "Unclear")) +
  facet_wrap(~article_lean, scales = 'free', ncol = 1) +
  theme_ctokita() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.2)
gg_articlelean
ggsave(gg_articlelean, file = paste0(outpath, "FMtweeters_ideologies_byarticlelean.png"), width = 90, height = 70, units = "mm", dpi = 400)



############################## Summary of ideology by article veracity ##############################

####################
# Plot: Ideology by article veracity
####################
# Summarise data
tweeter_veracity <- tweeter_scores %>% 
  group_by(article_fc_rating) %>% 
  summarise(ideological_extremity = mean(user_ideol_extremity, na.rm = TRUE)) %>% 
  filter(article_fc_rating %in% c("FM", "T"))

# Plot ideological distribution
gg_veracityideo <- tweeter_scores %>% 
  filter(!is.na(user_ideology),
         article_fc_rating %in% c("FM", "T")) %>% 
  ggplot(., aes(x = user_ideology, fill = ..x..)) +
  geom_histogram(position = 'identity', binwidth = 0.25) +
  xlab("Tweeter ideology") +
  ylab("Number of tweets") +
  scale_x_continuous(limits = c(-4.5, 4.5), expand = c(0, 0), breaks = seq(-5, 5, 1)) +
  scale_fill_gradientn(colors = ideol_pal, limits = c(-2, 2), oob = scales::squish) +
  theme_ctokita() +
  theme(aspect.ratio = 0.2, 
        legend.position = "none") +
  facet_wrap(~article_fc_rating, 
             scales = 'free', 
             ncol = 1, 
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity))
gg_veracityideo
ggsave(gg_veracityideo, file = paste0(outpath, "ideologies_byveracity.png"), width = 90, height = 70, units = "mm", dpi = 400)


# Plot ideological extremity
gg_veracityextr <- tweeter_scores %>% 
  filter(!is.na(user_ideology)) %>% 
  ggplot(., aes(x = user_ideol_extremity, fill = ..x..)) +
  geom_histogram(position = 'identity', binwidth = 0.25) +
  theme_ctokita() +
  facet_wrap(~article_fc_rating, 
             scales = 'free', 
             ncol = 1, 
             strip.position = "right")
gg_veracityextr
