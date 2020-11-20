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
source("scripts/_plot_themes/theme_ctokita.R")


####################
# Paramters for analysis: gpaths to data, paths for output, and filename
####################
# Paths to files/directories
intervention_dir <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/interventions/reduceviz0.5_t6/'
tweet_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/tweets_labeled.csv'


####################
# Load data
####################
# Intervention tweets
intervention_tweets <- data.frame()
for (file in list.files(intervention_dir, full.names = TRUE)) {
  tweets <- read.csv(file, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
    rename(time = relative_tweet_time) %>% 
    mutate(tweet_number = tweet_number+1)
  if ("X" %in% names(tweets)) {
    tweets <- tweets %>% 
      select(-X)
  }
  intervention_tweets <- rbind(intervention_tweets, tweets)
  rm(tweets)
}

# Original tweets and exposure data
article_data <- read.csv(tweet_path, header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
  mutate(article_ideology = article_con_feel - article_lib_feel) %>% 
  select(tweet_id, total_article_number, source_type, source_lean, article_fc_rating, article_lean, user_ideology) 

exposure_data <- read.csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/exposure/estimated_users_exposed_over_time.csv', 
                          header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
  mutate(tweet_number = tweet_number+1) %>%  #python zero index
  rename(time = relative_time) %>% 
  arrange(total_article_number, tweet_number)

exposure_timeseries <- merge(exposure_data, article_data, by = c("tweet_id", "total_article_number"), all = TRUE) %>% 
  mutate(replicate = -1)
rm(article_data, exposure_data)

# Create full dataset
exposure_timeseries <- exposure_timeseries[ , names(exposure_timeseries) %in% names(intervention_tweets)] %>% 
  mutate(tweet_type = "No intervention") %>% 
  filter(total_article_number %in% unique(intervention_tweets$total_article_number))
intervention_tweets <- intervention_tweets[ , names(intervention_tweets) %in% names(exposure_timeseries)] %>% 
  mutate(tweet_type = "Intervention") %>% 
  rbind(exposure_timeseries) %>% 
  mutate(tweet_type = factor(tweet_type, levels = c("No intervention", "Intervention")))
rm(exposure_timeseries)

# Add dummy rows of pre-first share for plotting purposes
# (1) Create dummy rows that make two "non tweets" in at two time points leading up to the first real tweet
n_articles <- length(unique(intervention_tweets$total_article_number)) #number of unique articles
dummy_rows <- intervention_tweets %>% 
  select(-tweet_number, -time) %>% 
  distinct(total_article_number, tweet_type, replicate, .keep_all = TRUE) %>% 
  mutate(user_id = NA,
         tweet_id = NA,
         new_exposed_users = 0,
         cumulative_exposed = 0)
dummy_row_time <- data.frame(time = rep(c(-2, -0.01), n_articles), 
                             tweet_number = rep(c(-1, 0), n_articles),
                             total_article_number = rep(unique(intervention_tweets$total_article_number), each = 2))
dummy_rows <- merge(dummy_rows, dummy_row_time, by = "total_article_number")
# (2) Join together
intervention_tweets <- intervention_tweets %>% 
  rbind(dummy_rows) %>% 
  arrange(total_article_number, replicate, tweet_number) %>% 
  mutate(hour_bin = cut(time, breaks = seq(-2, 50, 1), include.lowest = TRUE, right = FALSE, labels = seq(-2, 49))) %>%  #bin by hour tweet appeared
  mutate(hour_bin = as.numeric(as.character(hour_bin)))  #convert from factor to plain number

# Create measure of relative exposure (relative to actual tweet data)
intervention_tweets <- intervention_tweets %>% 
  group_by(total_article_number) %>% 
  mutate(relative_cumulative_exposed = cumulative_exposed / max(cumulative_exposed),
         relative_tweet_count = tweet_number / max(tweet_number))


  
  
####################
# Plot raw counts of cumulative exposure
####################
gg_exposed_raw <- ggplot(intervention_tweets, aes(x = time, y = relative_cumulative_exposed, color = tweet_type, group = replicate, alpha = tweet_type)) +
  geom_line(size = 0.3) +
  ylab("Relative number of users exposed") +
  xlab("Time since first article share (hrs)") +
  scale_color_manual(values = c("black", "red")) +
  scale_alpha_manual(values = c(1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  theme_ctokita() +
  facet_wrap(~total_article_number)
gg_exposed_raw
