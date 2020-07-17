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
library(scales)
source("scripts/_plot_themes/theme_ctokita.R")

####################
# Paramters for analysis: paths to data, paths for output, and filename
####################
tweeter_score_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/all_tweets_labeled.csv' #path to fitness cascade data
outpath <- 'output/ideology/'

# For labeling facet plots
label_veracity <- c("T" = "True news", 
                    "FM" = "Fake news",
                    "C" = "Conservative",
                    "L" = "Liberal",
                    "N" = "Neutral",
                    "U" = "Unclear")

# Color palette
line_color <- "#495867"
ideol_pal <- rev(brewer.pal(5, "RdBu"))
ideol_pal[3] <- "#e0e0e0"
ideol_dist_pal <- rev(brewer.pal(5, "PuOr"))
ideol_dist_pal[3] <- "#e0e0e0"



############################## Plot time series of article shares ##############################

####################
# Load data 
####################
# Read in data, Calculate time since first sharing of the story
tweeter_scores <- read.csv(tweeter_score_path, header = TRUE) %>% 
  mutate(article_ideology = article_con_feel - article_lib_feel,
         tweet_time_text = tweet_time,
         tweet_time = as.POSIXct(tweet_time, format = "%a %b %d %H:%M:%S %z %Y")) %>% 
  filter(!is.na(total_article_number)) %>% 
  arrange(total_article_number, tweet_time) %>% 
  group_by(total_article_number) %>% 
  mutate(article_first_time = min(tweet_time)) %>% 
  mutate(tweet_number = 1:length(tweet_time), #order tweets for plotting purposes
         relative_tweet_time = as.numeric( (tweet_time - article_first_time) / (60*60) ) ) %>%  #time diff is in seconds, so convert to hours
  mutate(relative_tweet_count = tweet_number / max(tweet_number))

# Add dummy rows of pre-first share for plotting purposes
# (1) Create empty dataframe for dummy rows
n_articles <- length(unique(tweeter_scores$total_article_number)) #number of unique articles
dummy_rows <- data.frame(matrix(NA, ncol = ncol(tweeter_scores), nrow = 2*n_articles))  #create empty dataframe
names(dummy_rows) <- names(tweeter_scores) #give same column names
# (2) Create unique set of article IDs and fact-check rating to add to our dummy rows
unique_article_ratings <- tweeter_scores %>% 
  select(total_article_number, article_fc_rating) %>% 
  unique()
# (3) Join together
tweeter_scores <- dummy_rows %>% 
  select(-article_fc_rating) %>% 
  mutate(relative_tweet_time = rep(c(-2, -0.01), n_articles),
         tweet_number = 0,
         relative_tweet_count = 0,
         total_article_number = rep(unique(tweeter_scores$total_article_number), each = 2)) %>% 
  merge(., unique_article_ratings) %>% 
  rbind(tweeter_scores, .) %>% 
  mutate(hour_bin = cut(relative_tweet_time, breaks = seq(-2, 50, 1), include.lowest = TRUE, right = FALSE, labels = seq(-2, 49))) %>%  #bin by hour tweet appeared
  mutate(hour_bin = as.numeric(as.character(hour_bin))) #convert from factor to plain number
  
# Calculate new tweets per hour
tweet_perhour <- tweeter_scores %>% 
  group_by(article_fc_rating, total_article_number, hour_bin) %>% 
  count(.)
tweet_perhour$n[tweet_perhour$hour_bin %in% c(-2, -1)] <- 0 #zero out the count of dummy rows

####################
# Plot tweeting of stories over time
####################
# New tweets over time
gg_tweettime <- tweet_perhour %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
  ggplot(., aes(x = hour_bin, y = n, group = total_article_number)) +
  geom_line(size = 0.2, alpha = 0.5, color = line_color) +
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

# cumulative tweets
gg_totaltweets <- tweeter_scores %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
  ggplot(., aes(x = relative_tweet_time, y = tweet_number, group = total_article_number)) +
  geom_vline(aes(xintercept = 24), 
             linetype = "dotted",
             size = 0.3,
             color = "grey60") +
  # geom_text(aes(x = 24, y = 100, label = "24 hr intervention\n"), angle = 270, color = "grey80") +
  geom_step(size = 0.2, alpha = 0.5, color = line_color) +
  scale_y_continuous(breaks = c(10^seq(1, 5)),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     trans = scales::pseudo_log_trans(base = 10)) +
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
ggsave(gg_totaltweets, filename = "output/timeseries/total_tweet_count.png", width = 90, height = 45, units = "mm", dpi = 400)

# Percentiage of tweets per story
gg_perctweets <- tweeter_scores %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
  ggplot(., aes(x = relative_tweet_time, y = relative_tweet_count, group = total_article_number)) +
  geom_vline(aes(xintercept = 24), 
             linetype = "dotted",
             size = 0.3,
             color = "grey60") +
  # geom_text(aes(x = 24, y = 100, label = "24 hr intervention\n"), angle = 270, color = "grey80") +
  geom_step(size = 0.2, alpha = 0.5, color = line_color) +
  xlab("Time since first article share (hrs)") +
  ylab("Proportion of story stweets") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity))
gg_perctweets
ggsave(gg_perctweets, filename = "output/timeseries/percentage_story_tweets.png", width = 90, height = 45, units = "mm", dpi = 400)

####################
# Plot saturation time of stories
####################
percentile <- 0.5
  
gg_saturationcount <- tweeter_scores %>% 
  filter(relative_tweet_count >= percentile,
         article_fc_rating %in% c("T", "FM")) %>% 
  group_by(total_article_number) %>% 
  filter(relative_tweet_count == min(relative_tweet_count)) %>% 
  ggplot(., aes(x = relative_tweet_time)) +
  geom_histogram(aes(y = stat(density)), binwidth = 2, color = 'white', fill = line_color) +
  xlab(paste0("Time to ", percentile*100, "% sharing saturation (hrs.)")) +
  theme_ctokita() +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity))
gg_saturationcount
ggsave(gg_saturationcount, filename = paste0("output/timeseries/story_saturation", percentile*100, ".png"), width = 55, height = 90, units = "mm", dpi = 400)


####################
# Plot ideology distance of tweeters relative to article content
####################
# Calculate distance between ideological category of user and article
# Filter out users without ideological scores, bin by hour
gg_ideoldisttime <- tweeter_scores %>% 
  #data processing
  filter(!is.na(user_ideology)) %>% 
  mutate(article_ideol_category = ifelse(article_lean == "L", -1, ifelse(article_lean == "C", 1, 0))) %>% 
  mutate(ideol_diff = user_ideol_category - article_ideol_category) %>% 
  mutate(ideol_distance = sign(user_ideol_category - article_ideol_category)) %>%
  group_by(article_fc_rating, total_article_number, hour_bin) %>% 
  count(ideol_distance) %>% 
  mutate(freq_ideol_distance = n / sum(n)) %>% 
  group_by(article_fc_rating, hour_bin, ideol_distance) %>% 
  summarise(freq_ideol_distance = mean(freq_ideol_distance)) %>% 
  filter(article_fc_rating %in% c("FM", "T"),
         hour_bin >= 0) %>%
  #graph
  ggplot(., aes(x = hour_bin, y = freq_ideol_distance, fill = factor(ideol_distance, levels = c(1, 0, -1)))) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_fill_manual(values = rev(ideol_dist_pal[c(1,3,5)]),
                    name = NULL,
                    labels = c("User more conservative than article", 
                               "Same ideology", 
                               "User more liberal than article")) +
  xlab("Time since first article share (hrs)") +
  ylab("Prop. of tweeters") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm")) +
  facet_wrap(~article_fc_rating, 
             strip.position = "right",
             ncol = 1,
             labeller = labeller(article_fc_rating = label_veracity))
gg_ideoldisttime
ggsave(gg_ideoldisttime, filename = "output/timeseries/ideology_relative_tweeters.png", width = 120, height = 45, units = "mm", dpi = 400)


####################
# Plot ideology tweeters by article time
####################
# Calculate average ideological distribution of tweeters over time
# Bin ideologies, filter out users without ideological scores, bin by hour
gg_ideoltime <- tweeter_scores %>% 
  #data processing
  filter(!is.na(user_ideology)) %>% 
  mutate(ideol_bin = cut(user_ideology, breaks = c(-6, -1, 1, 6), labels = c(-1, 0, 1), right = FALSE, include.lowest = TRUE)) %>%
  group_by(article_fc_rating, article_lean, total_article_number, hour_bin) %>% 
  count(ideol_bin) %>% 
  mutate(freq_ideol_bin = n / sum(n)) %>% 
  group_by(article_fc_rating, article_lean, hour_bin, ideol_bin) %>% 
  summarise(freq_ideol_bin = mean(freq_ideol_bin)) %>% 
  filter(article_fc_rating %in% c("FM", "T")) %>%
  # graph
  ggplot(., aes(x = hour_bin, y = freq_ideol_bin, fill = factor(ideol_bin, levels = c(1, 0, -1)))) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  scale_x_continuous(breaks = seq(0, 48, 6), limits = c(-1, 32)) +
  scale_fill_manual(values = rev(ideol_pal[c(1,3,5)]),
                    name = "User ideology",
                    labels = c("Conservative", "Moderate", "Liberal")) +
  xlab("Time since first article share (hrs)") +
  ylab("Prop. of tweeters") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm")) +
  facet_grid(article_fc_rating~article_lean,
             labeller = labeller(article_fc_rating = label_veracity))
gg_ideoltime
ggsave(gg_ideoltime, filename = "output/timeseries/ideology_tweeters.png", width = 120, height = 45, units = "mm", dpi = 400)

gg_ideoltime_raw <- tweeter_scores %>% 
  filter(article_fc_rating %in% c("FM", "T"),
         !is.na(user_ideology)) %>%
  ggplot(., aes(x = relative_tweet_time, y = user_ideology, color = user_ideology)) +
  geom_point(size = 0.3, stroke = 0, position = position_jitter(width = 0.1, height = 0.1)) +
  scale_color_gradientn(colors = ideol_pal, limits = c(-2, 2), oob = scales::squish) +
  scale_x_continuous(breaks = seq(0, 48, 12)) +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  xlab("Time since first article share (hrs)") +
  ylab("User ideology") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.position = "none") +
  facet_grid(article_fc_rating~article_lean,
             labeller = labeller(article_fc_rating = label_veracity))
gg_ideoltime_raw
ggsave(gg_ideoltime_raw, filename = "output/timeseries/ideology_raw_tweeters.png", width = 90, height = 45, units = "mm", dpi = 400)

# Calculate average ideological distribution of tweeters over time, broken out by source and article veracity
# Bin ideologies, filter out users without ideological scores, bin by hour
gg_ideoltimesource <- tweeter_scores %>% 
  filter(!is.na(user_ideology),
         article_fc_rating %in% c("FM", "T")) %>% 
  mutate(ideol_bin = cut(user_ideology, breaks = c(-6, -1, 1, 6), labels = c(-1, 0, 1), right = FALSE, include.lowest = TRUE)) %>%
  group_by(article_fc_rating, source_type, article_fc_rating, total_article_number, hour_bin) %>% 
  count(ideol_bin) %>% 
  mutate(freq_ideol_bin = n / sum(n)) %>% 
  group_by(article_fc_rating, source_type, article_fc_rating, hour_bin, ideol_bin) %>% 
  summarise(freq_ideol_bin = mean(freq_ideol_bin)) %>% 
  ggplot(., aes(x = hour_bin, y = freq_ideol_bin, fill = factor(ideol_bin, levels = c(1, 0, -1)))) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  scale_x_continuous(breaks = seq(0, 48, 6), limits = c(-1, 32)) +
  scale_fill_manual(values = rev(ideol_pal[c(1,3,5)]),
                    name = "User ideology",
                    labels = c("Conservative", "Moderate", "Liberal")) +
  xlab("Time since first article share (hrs)") +
  ylab("Prop. of tweeters") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm")) +
  facet_grid(article_fc_rating~source_type,
             labeller = labeller(article_fc_rating = label_veracity))
gg_ideoltimesource
ggsave(gg_ideoltimesource, filename = "output/timeseries/ideology_tweeters_bysourceandveracity.png", width = 120, height = 45, units = "mm", dpi = 400)


############################## Plot time series of article exposure ##############################

####################
# Load data
####################
# TEMP: Load data and bind
files <- list.files('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/timeseries/individual_articles/', full.names = TRUE)
for (file in files) {
  story_data <- read.csv(file) 
  story_data <- story_data %>%
    mutate(tweet_number = tweet_number+1) %>% 
    # add dummy rows of pre-tweet data for plotting purposes
    add_row(time = c(-2, -0.01), tweet_number = c(-1, 0), user_id = min(story_data$user_id), new_exposed_users = 0, cumulative_exposed = 0, total_article_number = unique(story_data$total_article_number))
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
rm(exposure_data, story_data)

####################
# Plot user exposed to stories over time
####################
# Total cumulative exposed
gg_exposuretime <- exposure_timeseries %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
  ggplot(., aes(x = time, y = cumulative_exposed, group = total_article_number)) +
  geom_step(size = 0.3, alpha = 0.5, color = line_color) +
  # scale_y_log10() +
  scale_y_continuous(breaks = c(10^seq(1, 7, 2)),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     trans = scales::pseudo_log_trans(base = 10)) +
  xlab("Time since first article share (hrs)") +
  ylab("Log users exposed") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity))
gg_exposuretime
ggsave(gg_exposuretime, filename = "output/timeseries/total_exposure.png", width = 90, height = 45, units = "mm", dpi = 400)

# Percentiage of tweets per story
gg_relexpostime <- exposure_timeseries %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
  ggplot(., aes(x = time, y = relative_cumulative_exposed, group = total_article_number)) +
  geom_step(size = 0.3, alpha = 0.5, color = line_color) +
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
ggsave(gg_relexpostime, filename = "output/timeseries/timeseries_percentageexposure.png", width = 90, height = 45, units = "mm", dpi = 400)

####################
# Plot time tweet number vs exposure
####################
# Merge data to create relevant dataset
exposure_vs_tweet <- tweeter_scores %>% 
  select(total_article_number, tweet_number, relative_tweet_count) %>% 
  merge(exposure_timeseries, ., all.x = TRUE)

# Plot
gg_expVnum <- exposure_vs_tweet %>% 
  filter(article_fc_rating %in% c("T", "FM"),
         tweet_number >= 0) %>% 
  ggplot(., aes(x = relative_tweet_count, y = relative_cumulative_exposed, group = total_article_number)) +
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linetype = "dashed", size = 0.3)+
  geom_line(size = 0.3, alpha = 0.5, color = line_color) +
  xlab("Proportion of total tweets") +
  ylab("Proportion of total users exposed") +
  theme_ctokita() +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "right",
             labeller = labeller(article_fc_rating = label_veracity),
             scales = "free_x")
gg_expVnum
ggsave(gg_expVnum, filename = "output/timeseries/tweetcount_vs_exposure.png", width = 50, height = 90, units = "mm", dpi = 400)
