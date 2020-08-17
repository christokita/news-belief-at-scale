########################################
#
# PLOT: Exposure over time
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
library(ggridges)
source("scripts/_plot_themes/theme_ctokita.R")

####################
# Paramters for analysis: paths to data, paths for output, and filename
####################
tweet_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/tweets_labeled.csv' #path to fitness cascade data
exposure_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/timeseries/individual_articles/' #temporary path to individually processed articles
outpath <- 'output/exposure/'

# For labeling facet plots
label_veracity <- c("T" = "True news", 
                    "FM" = "Fake news")

# Color palette
density_pal <- c()
ideol_pal <- rev(brewer.pal(5, "RdBu"))


####################
# Load data
####################
# Read in tweet data, Calculate time since first sharing of the story
tweets <- read.csv(tweet_path, header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
  mutate(article_ideology = article_con_feel - article_lib_feel,
         tweet_time_text = tweet_time,
         tweet_time = as.POSIXct(tweet_time, format = "%a %b %d %H:%M:%S %z %Y")) %>% 
  arrange(total_article_number, tweet_time) %>% 
  group_by(total_article_number) %>% 
  mutate(article_first_time = min(tweet_time)) %>% 
  mutate(tweet_number = 1:length(tweet_time), #order tweets for plotting purposes
         relative_tweet_time = as.numeric( (tweet_time - article_first_time) / (60*60) ) ) %>%  #time diff is in seconds, so convert to hours
  mutate(relative_tweet_count = tweet_number / max(tweet_number))

# Load exposure data (TEMPORARY FORMAT)
exposure_files <- list.files(exposure_path, full.names = TRUE)
for (file in exposure_files) {
  data <- read.csv(file, colClasses = c("user_id"="character"))
  if (!exists("exposure_data")) {
    exposure_data <- data
  } else {
    exposure_data <- rbind(exposure_data, data)
  }
  rm(data)
}

# Shape exposure data
exposure_data <- exposure_data %>% 
  gather(key = "ideology_bin", value = "count", -time, -tweet_number, -tweet_id, -user_id, -new_exposed_users, -cumulative_exposed, -total_article_number) %>% 
  mutate(ideology_bin = gsub("ideol_", "", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("^\\.", "-", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("_\\.", "_-", ideology_bin)) %>% 
  separate(ideology_bin, c("lower", "upper"), sep = "_", convert = TRUE) %>% 
  mutate(ideology_bin = (lower + upper) / 2) %>% 
  select(-lower, -upper)

####################
# Plot exposure ideology distributions
####################
# Grab example story for now
story <- 70
example_story <- exposure_data %>% 
  filter(total_article_number == story)

# 
ggplot(data = example_story, aes(x = ideology_bin, y = tweet_number, fill = count)) +
  geom_tile() +
  scale_fill_gradientn(name = "count", trans = "log", colors = c("blue", "red"), na.value = "blue") +
  theme_ctokita()

# sample tweet
axis_labels <- as.character(unique(exposure_data$ideology_bin))
axis_labels[seq(2, length(axis_labels)+1, 2)] <- ""
gg_ideol_tweet <- exposure_data %>% 
  filter(total_article_number == 109,
         tweet_number == 4) %>% 
  # mutate(ideology_bin = as.factor(ideology_bin)) %>%
  ggplot(., aes(x = as.factor(ideology_bin), y = count, fill = ideology_bin)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = axis_labels) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-2, 2), oob = scales::squish) +
  xlab("Follower ideology") +
  ylab("Count") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL)
gg_ideol_tweet
ggsave(gg_ideol_tweet, filename = paste0(outpath, "example_ideolexposed_", story, ".png"), dpi = 400, width = 90, heigh = 45, units = "mm")

ggplot(test, aes(x = ideology_bin, y = count)) +
  geom_line() +
  geom_point()
