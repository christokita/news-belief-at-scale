#########################################
# Name: `plot_exposure.R`
# Author: Chris Tokita
# Purpose: Plot estimated exposure to news articles over time, comparing by article veracity or by news source type
# Details:
#   (These R scripts assume the use of the `.Rproj` at top of the news-belief-at-scale/ repo. Otherwise, set the working directory to one level above this script.)
#
#   The Variables at the beginning of the script that are in all caps need to be set by the user:
#     `DATA_DIRECTORY`: path to the data directory. (Copies of data are currently stored on external hard drive and high-performance cluster.)
#     `GROUPING`:       determines whether the plots will break out tweet belief according to article veracity ("article_fc_rating") or the source of the article ("source_type").
# 
# Data In:
# `<data storage location>/data_derived/tweets/tweets_labeled.csv`: article tweets with article and tweeter metadata.
# `<data storage location>/data_derived/exposure/estimated_users_exposed_over_time.csv: estimated exposure to each article tweet.
# `<data storage location>/data/articles/evaluations.csv`: article fact-check ratings for each unique news article.
# `<data storage location>/data/articles/daily_articles.csv`: article source information for each unique news article.
#
# Data Out: Plots written to output sub-folder depending on if we are comparing article veracity or news source type. 
# `output/exposure/veracity/`
# `output/exposure/source_type/`
# 
# Machine: Chris' laptop
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
# Set parameters for analysis
####################
# Choose location of data
DATA_DIRECTORY <- "/Volumes/CKT-DATA/news-belief-at-scale/"

# Choose GROUPING of interest. Options: 
#     (1) article veracity: "article_fc_rating"
#     (2) source: "source_type"
GROUPING <- "article_fc_rating"


####################
# Prepare for analysis: set paths to data, paths for output, and color palettes for plotting
####################
# Set paths for data
tweet_path <- paste0(DATA_DIRECTORY, "data_derived/tweets/tweets_labeled.csv") #tweets
exposure_path <- paste0(DATA_DIRECTORY, "data_derived/exposure/estimated_users_exposed_over_time.csv") #estimated exposure per tweet
article_evaluation_path <- paste0(DATA_DIRECTORY, "data/articles/evaluations.csv") #article fact-check rating
article_source_path <- paste0(DATA_DIRECTORY, "data/articles/daily_articles.csv") #article source type

# Set path for plots
outpath <- 'output/exposure/'
if (GROUPING == "article_fc_rating") {
  subdir_out <- 'veracity/'
} else if(GROUPING == "source_type") {
  subdir_out <- 'source_type/'
}

# Set color palette
line_color <- "#495867"
ideol_pal <- rev(brewer.pal(5, "RdBu"))
ideol_pal[3] <- "#e0e0e0"
ideol_dist_pal <- rev(brewer.pal(5, "PuOr"))
ideol_dist_pal[3] <- "#e0e0e0"
ideol_limit <- 3 #limit beyond which we squish the color palette


####################
# Load and prepare data 
####################
# Read in tweet data for article info
article_data <- read.csv(tweet_path, header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
  mutate(article_ideology = article_con_feel - article_lib_feel) %>% 
  select(tweet_id, total_article_number, source_type, source_lean, article_fc_rating, article_lean, user_ideology) 

# Count number of articles per GROUPING (useful for average distributions) and add to article info
 if (GROUPING == "article_fc_rating") {
  article_group_counts <- read.csv(article_evaluation_path, header = TRUE) %>% 
    select(article_num, mode.of.FC) %>% 
    filter(article_num > 10) %>% 
    group_by(mode.of.FC) %>% 
    count() %>% 
    as.data.frame() %>% 
    rename(article_fc_rating = mode.of.FC,
           n_articles_in_grouping = n)
} else if(GROUPING == "source_type") {
  article_group_counts <-  read.csv(article_source_path, header = TRUE) %>% 
    select(total.article.number, source) %>% 
    mutate(source = gsub("_con|_lib|_unclear", "", source),
           source = gsub("ct", "mainstream", source),
           source = gsub("rss", "fringe", source)) %>% 
    filter(total.article.number > 10) %>% 
    group_by(source) %>% 
    count() %>% 
    as.data.frame() %>% 
    rename(source_type = source,
           n_articles_in_grouping = n)
}

article_data <- merge(article_data, article_group_counts, by = GROUPING)

# Load exposure data 
#
# NOTE:
# - for the raw count of exposure of followers we have ideology scores for user: users_exposed_over_time.csv
# - for the estimated ideology of all exposed followers: estimated_users_exposed_over_time.csv
exposure_data <- read.csv(exposure_path, header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
  mutate(tweet_number = tweet_number+1) %>%  #python zero index
  rename(time = relative_time) %>% 
  arrange(total_article_number, tweet_number)

# Merge in relevant article level data
# NOTE: adds one extra row, check after double checking with new data
exposure_timeseries <- merge(exposure_data, article_data, by = c("tweet_id", "total_article_number"), all = TRUE) 

# Add dummy rows of pre-first share for plotting purposes
# NOTE: This will add rows to the dataframe---two extra rows per article---but it will not effect user counts
# (1) Create empty dataframe for dummy rows
n_articles <- length(unique(exposure_timeseries$total_article_number)) #number of unique articles
dummy_rows <- data.frame(matrix(NA, ncol = ncol(exposure_timeseries), nrow = 2*n_articles))  #create empty dataframe
names(dummy_rows) <- names(exposure_timeseries) #give same column names
# (2) Create unique set of article IDs and fact-check rating to add to our dummy rows
unique_article_ratings <- article_data %>% 
  select(source_type, source_lean, total_article_number, article_fc_rating) %>% 
  unique()
# (3) Join together
exposure_timeseries <- dummy_rows %>% 
  select(-article_fc_rating, -source_type, -source_lean) %>% 
  mutate(time = rep(c(-2, -0.01), n_articles),
         tweet_number = rep(c(-1, 0), n_articles),
         new_exposed_users = 0, 
         cumulative_exposed = 0, 
         total_article_number = rep(unique(article_data$total_article_number), each = 2)) %>% 
  merge(unique_article_ratings, by = "total_article_number") %>% 
  rbind(exposure_timeseries, .) %>% 
  mutate(hour_bin = cut(time, breaks = seq(-2, 24*14, 1), include.lowest = TRUE, right = FALSE, labels = seq(-2, 24*14-1))) %>%  #bin by hour tweet appeared
  mutate(hour_bin = as.numeric(as.character(hour_bin))) %>%  #convert from factor to plain number
  group_by(total_article_number) %>% 
  mutate(relative_cumulative_exposed = cumulative_exposed / max(cumulative_exposed),
         relative_tweet_count = tweet_number / max(tweet_number)) %>% 
  arrange(total_article_number, tweet_number)

rm(dummy_rows, article_data, exposure_data)

# If analyzing by veracity, drop out non-True/False articles
if (GROUPING == "article_fc_rating") {
  exposure_timeseries <- exposure_timeseries %>% 
    filter(article_fc_rating %in% c("T", "FM"))
}

# Clean up some labels
exposure_timeseries <- exposure_timeseries %>% 
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "False/Misleading news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream outlet", ifelse(source_type == "fringe", "Fringe outlet", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )



############################## Plot time series of article exposure ##############################

####################
# PLOT: Total cumulative exposed per article (each line is an article)
####################
gg_exposuretime <- exposure_timeseries %>% 
  ggplot(., aes(x = time/24, y = cumulative_exposed, group = total_article_number)) +
  geom_step(size = 0.3, alpha = 0.5, color = line_color) +
  scale_x_continuous(breaks = seq(0, 14, 1),
                     limits = c(-0.1, 7)) +
  scale_y_continuous(breaks = c(10^seq(1, 8, 2)),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     trans = scales::pseudo_log_trans(base = 10)) +
  xlab("Time since first article share (days)") +
  ylab("Total users exposed") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", GROUPING)), 
             ncol = 1,
             strip.position = "top")

gg_exposuretime
ggsave(gg_exposuretime, filename = paste0(outpath, subdir_out, "total_exposed_time.pdf"), width = 90, height = 45, units = "mm", dpi = 400)


####################
# PLOT: Cumulative relative exposure per article (i.e., % of users eventually exposed to the article)
####################
gg_relexpostime <- exposure_timeseries %>% 
  ggplot(., aes(x = time/24, y = relative_cumulative_exposed, group = total_article_number)) +
  geom_step(size = 0.3, alpha = 0.5, color = line_color) +
  # scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 14, 1),
                       limits = c(-0.1, 7)) +
  xlab("Time since first article share (days)") +
  ylab("Proportion of article exposures") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", GROUPING)), 
             ncol = 1,
             strip.position = "top")

gg_relexpostime
ggsave(gg_relexpostime, filename = paste0(outpath, subdir_out, "relative_exposure_time.pdf"), width = 90, height = 45, units = "mm", dpi = 400)


####################
# PLOT: Cumulative users exposed by tweet number
####################
gg_exposuretweet <- exposure_timeseries %>% 
  ggplot(., aes(x = tweet_number, y = cumulative_exposed, group = total_article_number)) +
  geom_step(size = 0.3, alpha = 0.5, color = line_color) +
  scale_x_continuous(breaks = c(0, 1, 10, 100, 1000, 10000, 100000),
                     labels = scales::comma_format(accuracy = 1),
                     trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(breaks = c(10^seq(1, 7, 2)),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     trans = scales::pseudo_log_trans(base = 10)) +
  xlab("Tweet number") +
  ylab("Total users exposed") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(as.formula(paste(GROUPING, "~.")), 
             ncol = 2)
gg_exposuretweet
ggsave(gg_exposuretweet, filename = paste0(outpath, subdir_out, "total_exposed_tweetnumber.pdf"), width = 90, height = 45, units = "mm", dpi = 400)


####################
# PLOT: Proportion of article tweets X Proportion of total article exposure
####################
gg_expVnum <- exposure_timeseries %>% 
  filter(tweet_number > -1) %>% 
  ggplot(., aes(x = relative_tweet_count, y = relative_cumulative_exposed, group = total_article_number, color = total_article_number)) +
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linetype = "dashed", size = 0.3)+
  geom_line(size = 0.3, alpha = 0.5, color = line_color) +
  xlab("Proportion of article tweets") +
  ylab("Proportion of article exposures") +
  theme_ctokita() +
  facet_wrap(as.formula(paste("~", GROUPING)), 
             ncol = 1,
             strip.position = "top",
             scales = "free_x")

gg_expVnum
ggsave(gg_expVnum, filename = paste0(outpath, subdir_out, "relative_tweet_vs_exposure.pdf"), width = 50, height = 90, units = "mm", dpi = 400)



############################## Plot article exposure by ideology ##############################

# Prep data (data is melted to make one ideological bin per tweet per row)
exposure_ideol <- exposure_timeseries %>% 
  filter(hour_bin >= 0) %>% 
  select(-source_lean, -relative_cumulative_exposed, -relative_tweet_count) %>% 
  gather(key = "ideology_bin", value = "count", 
         -time, -tweet_number, -tweet_id, -user_id, -user_ideology, -follower_count, -new_exposed_users, -cumulative_exposed, -total_article_number, -hour_bin, -source_type, -article_fc_rating, -article_lean, -n_articles_in_grouping) %>% 
  mutate(ideology_bin = gsub("ideol_", "", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("^\\.", "-", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("_\\.", "_-", ideology_bin)) %>% 
  separate(ideology_bin, c("lower", "upper"), sep = "_", convert = TRUE) %>% 
  mutate(ideology_bin = (lower + upper) / 2) %>% 
  select(-lower, -upper)


####################
# PLOT: Ideology of exposed users
####################
gg_ideol_total <- exposure_ideol %>% 
  group_by(!!sym(GROUPING), ideology_bin) %>% 
  summarise(count = sum(count, na.rm = TRUE),
            avg_count = sum(count, na.rm = TRUE) / length(unique(total_article_number))) %>% 
  ggplot(., aes(x = ideology_bin, y = count, fill = ideology_bin)) +
  geom_bar(stat = "identity", width = 0.5, color = NA) +
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 2)) +
  scale_y_continuous(expand = c(0, 0), 
                     labels = comma) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("User ideology") +
  ylab("Total number of exposed users") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", GROUPING)), 
             ncol = 1,
             strip.position = "top",
             scales = "free")

gg_ideol_total
ggsave(gg_ideol_total, filename = paste0(outpath, subdir_out, "ideol_total_exposed.pdf"), width = 45, height = 90, units = "mm", dpi = 400)

####################
# PLOT: Average number of users exposed per article, broken out by user ideology
####################
gg_ideol_avg <- exposure_ideol %>% 
  group_by(!!sym(GROUPING), ideology_bin) %>% 
  summarise(avg_count = sum(count, na.rm = TRUE) / unique(n_articles_in_grouping)) %>% 
  ggplot(., aes(x = ideology_bin, y = avg_count, fill = ideology_bin)) +
  geom_bar(stat = "identity", width = 0.5, color = NA) +
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 2)) +
  scale_y_continuous(expand = c(0, 0), 
                     labels = comma) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("User ideology") +
  ylab("Avg. number of exposed users per article") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", GROUPING)), 
             ncol = 1,
             strip.position = "top",
             scales = "free")

gg_ideol_avg
ggsave(gg_ideol_avg, filename = paste0(outpath, subdir_out, "ideol_avg_exposed.pdf"), width = 45, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Average ideological distribution of users exposed to an article
####################
gg_ideol_dist <- exposure_ideol %>% 
  # filter(total_article_number == 28) %>%
  # For each article, determine proportion exposed by ideology bin
  group_by(!!sym(GROUPING), ideology_bin, total_article_number) %>% 
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  group_by(total_article_number, n_articles_in_grouping) %>% 
  mutate(exposed_prop = count / sum(count),
         exposed_prop = ifelse( is.na(exposed_prop), 0, exposed_prop)) %>% 
  # Now determine average distribution shape by article GROUPING
  group_by(!!sym(GROUPING), ideology_bin, n_articles_in_grouping) %>% 
  summarise(avg_exposed_prop = sum(exposed_prop) / length(unique(total_article_number))) %>% 
  # Plot
  ggplot(., aes(x = ideology_bin, y = avg_exposed_prop, fill = ideology_bin)) +
  geom_bar(stat = "identity", width = 0.5, color = NA) +
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.05),
                     limits = c(0, 0.2),
                     expand = c(0, 0)) +
  scale_color_gradientn(colours = ideol_pal, 
                        limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_fill_gradientn(colours = ideol_pal, 
                       limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("Exposed user ideology") +
  ylab("Avg. proportion of article exposure per article") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", GROUPING)), 
             ncol = 1,
             strip.position = "top",
             scales = "free_x")

gg_ideol_dist  
ggsave(gg_ideol_dist, filename = paste0(outpath, subdir_out, "ideol_avg_exposure_distribution.pdf"), width = 45, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Proportion of newly exposed users (time series)
####################
gg_ideoltime <- exposure_ideol %>% 
  filter(hour_bin >= 0) %>% 
  group_by(!!sym(GROUPING), hour_bin, ideology_bin) %>% 
  summarise(count = sum(count)) %>%
  ggplot(., aes(x = hour_bin, y = count, fill = ideology_bin, color = ideology_bin)) +
  geom_bar(position = "fill", stat = "identity", width = 1, size = 0.01) +
  scale_fill_gradientn(colours = ideol_pal, 
                       name = "User\nideology",
                       limits = c(-2, 2), 
                       oob = squish) +
  scale_color_gradientn(colours = ideol_pal, 
                       name = "User\nideology",
                       limits = c(-2, 2), 
                       oob = squish) +  
  scale_x_continuous(breaks = seq(0, 72, 12),
                     expand = c(0, 0),
                     limits = c(-0.5, 72.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Time since first article share (hrs)") +
  ylab("Proportion of newly exposed users") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_rect(size = 0.6, fill = NA),
        axis.line = element_blank()) +
  facet_wrap(as.formula(paste("~", GROUPING)), 
             ncol = 1,
             strip.position = "top",
             scales = "free")

gg_ideoltime
ggsave(gg_ideoltime, filename = paste0(outpath, subdir_out, "ideol_exposed_hourbin.pdf"), width = 90, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Ideologies of users exposed to each of the false/misleading news articles
####################
gg_fake_exposure_article <- exposure_ideol %>% 
  filter(article_fc_rating == "False/Misleading news") %>% 
  group_by(total_article_number, ideology_bin) %>% 
  summarise(count = sum(count, na.rm = TRUE)) %>% 
  ggplot(., aes(x = ideology_bin, y = count, fill = ideology_bin)) +
  geom_bar(stat = "identity", width = 0.5, color = NA) +
  scale_x_continuous(breaks = seq(-3, 6, 3),
                     limits = c(-3, 6)) +
  scale_y_continuous(labels = comma) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("User ideology") +
  ylab("Total users exposed to article") +
  theme_ctokita() +
  theme(legend.position = "none") +
  facet_wrap(~total_article_number,
             scales = "free")

gg_fake_exposure_article
ggsave(gg_fake_exposure_article, filename = paste0(outpath, "veracity/fake_news_exposure_by_article.pdf"), width = 180, height = 180, units = "mm", dpi = 400)



############################## Plot cross-ideology exposure ##############################
# Here, we are interested in assessing when conservatives tweeted/retweeted an article and exposed liberals, and vice versa.
# Cross ideology exposure is an instance in which a user exposed another user with the opposite ideological sign (+/-).
# In this analysis, liberals are users with ideology <0, while conservatives are users with ideology >0.

####################
# PLOT: Average cross-ideology exposure per tweet, broken out by tweeter ideology
####################
# Re-load exposure data 
exposure_data <- read.csv('/Volumes/CKT-DATA/news-belief-at-scale/data_derived/exposure/estimated_users_exposed_over_time.csv', 
                          header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
  mutate(tweet_number = tweet_number+1) %>%  #python zero index
  rename(time = relative_time) %>% 
  arrange(total_article_number, tweet_number)

# Add in user_ideology
user_ideologies <- read.csv('/Volumes/CKT-DATA/news-belief-at-scale/data_derived/tweets/tweets_labeled.csv') %>% 
  select(user_id, user_ideology) %>% 
  distinct()

exposure_by_ideology <- user_ideologies %>% 
  merge(exposure_data, by = "user_id", all.y = TRUE) %>% 
  mutate(bin = cut(user_ideology, breaks = seq(-5.5, 5.5, 0.5)))%>% 
  mutate(lower_edge = as.numeric( gsub("^[^0-9]([-\\.0-9]+),.*", "\\1", bin, perl = TRUE) ),
         upper_edge = as.numeric( gsub(".*,([-\\.0-9]+)[^0-9]$", "\\1", bin, perl = TRUE) )) %>% 
  mutate(user_ideology_bin = (lower_edge + upper_edge) / 2) %>% 
  select(-bin, -lower_edge, -upper_edge)

# Calculate left-vs-right exposure by tweet
exposure_by_ideology <- exposure_by_ideology %>% 
  mutate(exposed_left = ideol_.3.0_.2.5 + ideol_.2.5_.2.0 + ideol_.2.0_.1.5 + ideol_.1.5_.1.0 + ideol_.1.0_.0.5 + ideol_.0.5_0.0,
         exposed_right = ideol_0.0_0.5 + ideol_0.5_1.0 + ideol_1.0_1.5 + ideol_1.5_2.0 + ideol_2.0_2.5 + ideol_2.5_3.0 + ideol_3.0_3.5 + ideol_3.5_4.0 + ideol_4.0_4.5 + ideol_4.5_5.0 + ideol_5.0_5.5) %>%
  mutate(cross_exposure = ifelse(user_ideology < 0, exposed_right, exposed_left),
         cross_exposure_percent = cross_exposure / (exposed_left + exposed_right)) 
  # select(time, tweet_number, tweet_id, user_id, user_ideology, user_ideology_bin, follower_count, new_exposed_users, cumulative_exposed, total_article_number, exposed_left, exposed_right, cross_exposure, cross_exposure_percent)

cross_exposure_data <- exposure_by_ideology %>% 
  group_by(user_ideology_bin) %>% 
  summarise(total_tweets = length(cross_exposure),
            cross_exposure_sum = sum(cross_exposure, na.rm = TRUE),
            cross_exposure_per_tweet = cross_exposure_sum / total_tweets,
            cross_exposure_percent_avg = mean(cross_exposure_percent, na.rm = TRUE))

# Plot
gg_crossexposure <- ggplot(cross_exposure_data %>% filter(!is.na(user_ideology_bin)), aes(x = user_ideology_bin, y = cross_exposure_per_tweet, fill = user_ideology_bin)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.3) +
  scale_x_continuous(limits = c(-3, 5), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 1)) +
  scale_y_continuous(limits = c(0, 12000),
                     expand = c(0, 0),
                     labels = scales::comma) +
  scale_fill_gradientn(colors = ideol_pal, limits = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("Tweeter ideology") +
  ylab("Avg. cross-ideology\nexposure per tweet") +
  theme_ctokita() +
  theme(legend.position = "none")

gg_crossexposure
ggsave(gg_crossexposure, filename = paste0(outpath, "avg_crossideology_exposure.pdf"), width = 55, height = 45, units = "mm", dpi = 400)


####################
# PLOT: Exposure diversity per tweet, broken out by tweeter ideology
####################
# 
# The idea here is to see what types of users (in terms of ideology) had the most diverse audience.
#

# Calculate standard deviation of user ideologies exposed to each tweet
exposure_diversity <- exposure_by_ideology %>% 
  select(user_ideology_bin, user_id, tweet_id, ideol_.3.0_.2.5:ideol_5.0_5.5) %>% 
  gather("bin", "count", -user_ideology_bin, -user_id, -tweet_id) %>% 
  mutate(bin = gsub("_\\.", "_-", bin)) %>% 
  mutate(lower_edge = as.numeric( gsub("^ideol_([-\\.0-9]+)_.*$", "\\1", bin, perl = TRUE) ),
         upper_edge = as.numeric( gsub("^ideol_[-\\.0-9]+_([-\\.0-9]+)$", "\\1", bin, perl = TRUE) )) %>% 
  mutate(exposed_bin = (lower_edge + upper_edge) / 2) %>% 
  select(-bin, -lower_edge, -upper_edge) %>% 
  mutate(weighted_ideol = exposed_bin * count) %>% 
  group_by(tweet_id, user_ideology_bin) %>% 
  summarise(exposure_mean = sum(weighted_ideol) / sum(count),
            exposure_sd = sqrt( sum((exposed_bin - exposure_mean)^2 * count) / sum(count) ) )

# Plot
gg_exposure_sd <- exposure_diversity %>% 
  filter(!is.na(user_ideology_bin)) %>% 
  ggplot(., aes(x = user_ideology_bin, y = exposure_sd, fill = user_ideology_bin, group = user_ideology_bin)) +
  geom_violin(size = 0, color = NA) +
  stat_summary(fun.y = mean, geom = "point", size = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.3) +
  scale_x_continuous(limits = c(-3, 5), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 1)) +
  scale_y_continuous(limits = c(0, 2),
                     breaks = seq(0, 3, 0.5),
                     expand = c(0, 0),
                     labels = scales::comma) +
  scale_fill_gradientn(colors = ideol_pal, limits = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("Tweeter ideology") +
  ylab("Standard deviation of\nideologies exposed per tweet") +
  theme_ctokita() +
  theme(legend.position = "none")

gg_exposure_sd
ggsave(gg_exposure_sd, filename = paste0(outpath, "exposure_ideology_diversity.pdf"), width = 55, height = 45, units = "mm", dpi = 400)
