########################################
#
# PLOT: Tweets over time
#
########################################

####################
# Load packages
####################
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
source("scripts/_plot_themes/theme_ctokita.R")

####################
# Paramters for analysis: paths to data, paths for output, and filename
####################
tweeter_score_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/labeled_tweets.csv' #path to fitness cascade data
outpath <- 'output/ideology/'

# For labeling facet plots
label_veracity <- c("T" = "Real", 
                    "FM" = "Fake")


############################## Plot time series of article shares ##############################

####################
# Load data 
####################
# Read in data
tweeter_scores <- read.csv(tweeter_score_path, header = TRUE) %>% 
  mutate(article_ideology = article_con_feel - article_lib_feel,
         tweet_time_text = tweet_time,
         tweet_time = as.POSIXct(tweet_time, format = "%a %b %d %H:%M:%S %z %Y")) %>% 
  filter(!is.na(total_article_number)) %>% 
  arrange(total_article_number, tweet_time)

# Calculate time since first sharing of the story
tweeter_scores <- tweeter_scores %>% 
  group_by(total_article_number) %>% 
  mutate(article_first_time = min(tweet_time)) %>% 
  mutate(tweet_number = 1:length(tweet_time), #order tweets for plotting purposes
         relative_tweet_time = as.numeric( (tweet_time - article_first_time) / (60*60) ) ) %>%  #time diff is in seconds, so convert to hours
  mutate(tweet_count_percentile = tweet_number / max(tweet_number),
         hour_bin = cut(relative_tweet_time, breaks = seq(0, 50, 1), include.lowest = TRUE, labels = seq(0, 49)))

tweet_perhour <- tweeter_scores %>% 
  group_by(article_fc_rating, total_article_number, hour_bin) %>% 
  count(.)


####################
# Plot tweeting of stories over time
####################
pal <- rev(brewer.pal(5, "RdYlBu"))

# New tweets over time
gg_tweettime <- tweet_perhour %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
  ggplot(., aes(x = hour_bin, y = n, group = total_article_number)) +
  geom_line(size = 0.2, alpha = 0.5, color = "#495867") +
  xlab("Time since first article share (hrs)") +
  ylab("Tweets") +
  scale_y_log10() +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity))
gg_tweettime

gg_tweettime <- tweet_perhour %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
  mutate(article = as.factor(as.character(total_article_number))) %>% 
  ggplot(., aes(x = hour_bin, y = article, fill = n)) +
  # geom_line(size = 0.2, alpha = 0.5, color = "#495867") +
  geom_tile() +
  xlab("Time since first article share (hrs)") +
  ylab("Tweets") +
  # scale_y_log10() +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity))
gg_tweettime

# cumulative tweets
gg_totaltweets <- tweeter_scores %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
  ggplot(., aes(x = relative_tweet_time, y = tweet_number, group = total_article_number)) +
  geom_vline(aes(xintercept = 24), 
             linetype = "dotted",
             size = 0.3,
             color = "grey60") +
  # geom_text(aes(x = 24, y = 100, label = "24 hr intervention\n"), angle = 270, color = "grey80") +
  geom_line(size = 0.2, alpha = 0.5, color = "#495867") +
  scale_y_log10() +
  # scale_x_log10() +
  xlab("Time since first article share (hrs)") +
  ylab("Log total tweets") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity))
gg_totaltweets
ggsave(gg_totaltweets, filename = "output/timeseries/timeseries_totaltweetcount.png", width = 90, height = 45, units = "mm", dpi = 400)

# Percentiage of tweets per story
gg_perctweets <- tweeter_scores %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
  ggplot(., aes(x = relative_tweet_time, y = tweet_count_percentile, group = total_article_number)) +
  geom_vline(aes(xintercept = 24), 
             linetype = "dotted",
             size = 0.3,
             color = "grey60") +
  # geom_text(aes(x = 24, y = 100, label = "24 hr intervention\n"), angle = 270, color = "grey80") +
  geom_line(size = 0.2, alpha = 0.5, color = "#495867") +
  xlab("Time since first article share (hrs)") +
  ylab("Proportion of story stweets") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity))
gg_perctweets
ggsave(gg_perctweets, filename = "output/timeseries/timeseries_percentagestorytweets.png", width = 90, height = 45, units = "mm", dpi = 400)

####################
# Plot saturation time of stories
####################
percentile <- 0.5
  
gg_saturationcount <- tweeter_scores %>% 
  filter(tweet_count_percentile >= percentile,
         article_fc_rating %in% c("T", "FM")) %>% 
  group_by(total_article_number) %>% 
  filter(tweet_count_percentile == min(tweet_count_percentile)) %>% 
  ggplot(., aes(x = relative_tweet_time)) +
  geom_histogram(aes(y = stat(density))) +
  xlab(paste0("Time to ", percentile*100, "% sharing saturation")) +
  theme_ctokita() +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity))
gg_saturationcount


############################## Plot time series of article exposure ##############################

####################
# Plot user exposed to stories over time
####################
# TEMP: Load data and bind
files <- list.files('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/timeseries/individual_articles/', full.names = TRUE)
for (file in files) {
  story_data <- read.csv(file) 
  story_data <- story_data %>%
    # add dummy rows of pre-tweet data for plotting purposes
    add_row(X = 0, time = c(-1, -0.01), tweet_number = c(-2, -1), user_id = min(story_data$user_id), new_exposed_users = 0, cumulative_exposed = 0, total_article_number = unique(story_data$total_article_number))
  if (!exists("exposure_data")) {
    exposure_data <- story_data
  } else {
    exposure_data <- rbind(exposure_data, story_data)
  }
}

# Merge in relevant article level data
article_data <- tweeter_scores %>% 
  select(user_id, total_article_number, source_lean, article_fc_rating, article_lean, user_ideology, article_ideology)
exposure_timeseries <- merge(exposure_data, article_data, by = c("user_id", "total_article_number"), all.x = TRUE) %>% 
  group_by(total_article_number) %>% 
  mutate(relative_cumulative_exposed = cumulative_exposed / max(cumulative_exposed))

# Total cumulative exposed
gg_exposuretime <- ggplot(exposure_timeseries, aes(x = time, y = cumulative_exposed, group = total_article_number)) +
  geom_line(size = 0.3, alpha = 0.5) +
  # scale_y_log10() +
  scale_y_continuous(labels = scales::comma) +
  xlab("Time since first article share (hrs)") +
  ylab("Total users exposed") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity))
gg_exposuretime

# Percentiage of tweets per story
gg_relexpostime <- ggplot(exposure_timeseries, aes(x = time, y = relative_cumulative_exposed, group = total_article_number)) +
  geom_line(size = 0.3, alpha = 0.5) +
  # scale_y_log10() +
  scale_y_continuous(labels = scales::comma) +
  xlab("Time since first article share (hrs)") +
  ylab("Proportion of total users exposed") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity))
gg_relexpostime



