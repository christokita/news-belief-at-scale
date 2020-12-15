########################################
#
# PLOT: Exposure to news articles over time, comparing by article veracity or by news source type
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
# Parameters for analysis: grouping of interest, paths to data, paths for output, and filename
####################
# Choose grouping of interest. Options: 
#     (1) article veracity: "article_fc_rating"
#     (2) source: "source_type"
grouping <- "source_type"

# Paths to files/directories
tweet_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/tweets_labeled.csv'
if (grouping == "article_fc_rating") {
  outpath <- 'output/exposure/veracity/'
} else if(grouping == "source_type") {
  outpath <- 'output/exposure/source_type/'
}

# Color palette
line_color <- "#495867"
ideol_pal <- rev(brewer.pal(5, "RdBu"))
ideol_pal[3] <- "#e0e0e0"
ideol_dist_pal <- rev(brewer.pal(5, "PuOr"))
ideol_dist_pal[3] <- "#e0e0e0"


####################
# Load and prep data 
####################
# Read in tweet data for article info
article_data <- read.csv(tweet_path, header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
  mutate(article_ideology = article_con_feel - article_lib_feel) %>% 
  select(tweet_id, total_article_number, source_type, source_lean, article_fc_rating, article_lean, user_ideology) 

# Load exposure data 
#
# NOTE:
# - for the raw count of exposure of followers we have ideology scores for user: users_exposed_over_time.csv
# - for the estimated ideology of all exposed followers: estimated_users_exposed_over_time.csv
exposure_data <- read.csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/exposure/estimated_users_exposed_over_time.csv', 
                          header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
  mutate(tweet_number = tweet_number+1) %>%  #python zero index
  rename(time = relative_time) %>% 
  arrange(total_article_number, tweet_number)

# Merge in relevant article level data
# NOTE: adds one extra row, check after double checking with new data
exposure_timeseries <- merge(exposure_data, article_data, by = c("tweet_id", "total_article_number"), all = TRUE) 

# Add dummy rows of pre-first share for plotting purposes
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
  mutate(hour_bin = cut(time, breaks = seq(-2, 50, 1), include.lowest = TRUE, right = FALSE, labels = seq(-2, 49))) %>%  #bin by hour tweet appeared
  mutate(hour_bin = as.numeric(as.character(hour_bin))) %>%  #convert from factor to plain number
  group_by(total_article_number) %>% 
  mutate(relative_cumulative_exposed = cumulative_exposed / max(cumulative_exposed),
         relative_tweet_count = tweet_number / max(tweet_number)) %>% 
  arrange(total_article_number, tweet_number)

rm(dummy_rows, article_data, exposure_data)

# If analyzing by veracity, drop out non-True/False articles
if (grouping == "article_fc_rating") {
  exposure_timeseries <- exposure_timeseries %>% 
    filter(article_fc_rating %in% c("T", "FM"))
}

# Clean up some labels
exposure_timeseries <- exposure_timeseries %>% 
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "Fake news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream", ifelse(source_type == "fringe", "Fringe", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )



############################## Plot time series of article exposure ##############################

####################
# Total cumulative exposed
####################
gg_exposuretime <- exposure_timeseries %>% 
  ggplot(., aes(x = time, y = cumulative_exposed, group = total_article_number)) +
  geom_step(size = 0.3, alpha = 0.5, color = line_color) +
  # scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_y_continuous(breaks = c(10^seq(1, 7, 2)),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     trans = scales::pseudo_log_trans(base = 10)) +
  xlab("Time since first article share (hrs)") +
  ylab("Total users exposed") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right")
gg_exposuretime
ggsave(gg_exposuretime, filename = paste0(outpath, "total_exposed_time.pdf"), width = 90, height = 45, units = "mm", dpi = 400)


####################
# Percentiage of tweets per story
####################
gg_relexpostime <- exposure_timeseries %>% 
  ggplot(., aes(x = time, y = relative_cumulative_exposed, group = total_article_number)) +
  geom_step(size = 0.3, alpha = 0.5, color = line_color) +
  # scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_y_continuous(labels = scales::comma) +
  xlab("Time since first article share (hrs)") +
  ylab("Proportion of total users exposed") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right")
gg_relexpostime
ggsave(gg_relexpostime, filename = paste0(outpath, "percentage_exposed_time.pdf"), width = 90, height = 45, units = "mm", dpi = 400)


####################
# Cumulative exposed by tweet number
####################
gg_exposuretweet <- exposure_timeseries %>% 
  ggplot(., aes(x = tweet_number, y = cumulative_exposed, group = total_article_number)) +
  geom_step(size = 0.3, alpha = 0.5, color = line_color) +
  scale_y_continuous(breaks = c(10^seq(1, 7, 2)),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     trans = scales::pseudo_log_trans(base = 10)) +
  xlab("Tweet number") +
  ylab("Total users exposed") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right")
gg_exposuretweet
ggsave(gg_exposuretweet, filename = paste0(outpath, "total_exposed_tweetnumber.pdf"), width = 90, height = 45, units = "mm", dpi = 400)


####################
# Relative tweet time tweet number vs exposure
####################
# Plot
gg_expVnum <- exposure_timeseries %>% 
  filter(tweet_number > -1) %>% 
  ggplot(., aes(x = relative_tweet_count, y = relative_cumulative_exposed, group = total_article_number, color = total_article_number)) +
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linetype = "dashed", size = 0.3)+
  geom_line(size = 0.3, alpha = 0.5, color = line_color) +
  xlab("Proportion of total tweets") +
  ylab("Proportion of total users exposed") +
  theme_ctokita() +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "free_x")
gg_expVnum
ggsave(gg_expVnum, filename = paste0(outpath, "relative_tweet_vs_exposure.pdf"), width = 50, height = 90, units = "mm", dpi = 400)



############################## Plot article exposure by ideology ##############################

# Prep data
exposure_ideol <- exposure_timeseries %>% 
  filter(hour_bin >= 0) %>% 
  select(-source_lean, -relative_cumulative_exposed, -relative_tweet_count) %>% 
  gather(key = "ideology_bin", value = "count", 
         -time, -tweet_number, -tweet_id, -user_id, -user_ideology, -follower_count, -new_exposed_users, -cumulative_exposed, -total_article_number, -hour_bin, -source_type, -article_fc_rating, -article_lean) %>% 
  mutate(ideology_bin = gsub("ideol_", "", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("^\\.", "-", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("_\\.", "_-", ideology_bin)) %>% 
  separate(ideology_bin, c("lower", "upper"), sep = "_", convert = TRUE) %>% 
  mutate(ideology_bin = (lower + upper) / 2) %>% 
  select(-lower, -upper) %>% 
  # count number of articles per grouping (useful for average distributions)
  group_by(!!sym(grouping)) %>% 
  mutate(n_articles_in_grouping = length(unique(total_article_number)))


####################
# Total ideologies exposed to articles by type
####################
gg_ideol_total <- exposure_ideol %>% 
  group_by(!!sym(grouping), ideology_bin) %>% 
  summarise(count = sum(count, na.rm = TRUE),
            avg_count = sum(count, na.rm = TRUE) / length(unique(total_article_number))) %>% 
  ggplot(., aes(x = ideology_bin, y = count, fill = ideology_bin)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(-5.5, 5.5), expand = c(0, 0), breaks = seq(-5, 5, 1)) +
  scale_y_continuous(labels = comma) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-2, 2), oob = scales::squish) +
  xlab("User ideology") +
  ylab("Total users exposed to articles") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "free")
gg_ideol_total
ggsave(gg_ideol_total, filename = paste0(outpath, "ideol_total_exposed.pdf"), width = 90, height = 90, units = "mm", dpi = 400)


####################
# Avg ideologies exposed per article (by type)
####################
gg_ideol_avg <- exposure_ideol %>% 
  group_by(!!sym(grouping), ideology_bin) %>% 
  summarise(avg_count = sum(count, na.rm = TRUE) / length(unique(total_article_number))) %>% 
  ggplot(., aes(x = ideology_bin, y = avg_count, fill = ideology_bin)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(-5.5, 5.5), expand = c(0, 0), breaks = seq(-5, 5, 1)) +
  scale_y_continuous(labels = comma) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-2, 2), oob = scales::squish) +
  xlab("User ideology") +
  ylab("Avg. users exposed to article") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "free")
gg_ideol_avg
ggsave(gg_ideol_avg, filename = paste0(outpath, "ideol_avg_exposed.pdf"), width = 90, height = 90, units = "mm", dpi = 400)


####################
# Ideological distributions of exposure
####################
gg_ideol_dist <- exposure_ideol %>% 
  # filter(total_article_number == 28) %>%
  # For each article, determine proportion exposed by ideology bin
  group_by(!!sym(grouping), ideology_bin, total_article_number, n_articles_in_grouping) %>% 
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  group_by(total_article_number, n_articles_in_grouping) %>% 
  mutate(exposed_prop = count / sum(count),
         exposed_prop = ifelse( is.na(exposed_prop), 0, exposed_prop)) %>% 
  # Now determine average distribution shape by article grouping
  group_by(!!sym(grouping), ideology_bin, n_articles_in_grouping) %>% 
  summarise(avg_exposed_prop = sum(exposed_prop) / unique(n_articles_in_grouping)) %>% 
  # Plot
  ggplot(., aes(x = ideology_bin, y = avg_exposed_prop, fill = ideology_bin)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(-5.5, 5.5), expand = c(0, 0), breaks = seq(-5, 5, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-2, 2), oob = scales::squish) +
  xlab("User ideology") +
  ylab("Avg. proportion of article exposure") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "free_x")
gg_ideol_dist  
  
ggsave(gg_ideol_dist, filename = paste0(outpath, "ideol_avg_exposure_distribution.pdf"), width = 90, height = 90, units = "mm", dpi = 400)


####################
# Exposure time series
####################
gg_ideoltime <- exposure_ideol %>% 
  filter(hour_bin >= 0) %>% 
  group_by(!!sym(grouping), hour_bin, ideology_bin) %>% 
  summarise(count = sum(count)) %>%
  ggplot(., aes(x = hour_bin, y = count, fill = ideology_bin)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  scale_fill_gradientn(colours = ideol_pal, 
                       name = "User\nideology",
                       limits = c(-2, 2), 
                       oob = squish) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  xlab("Time since first article share (hrs)") +
  ylab("New users exposed") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm")) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "free")
gg_ideoltime
ggsave(gg_ideoltime, filename = paste0(outpath, "ideol_exposed_hourbin.pdf"), width = 90, height = 90, units = "mm", dpi = 400)



####################
# Find exposure to fake news
####################
gg_fake_expos <- exposure_ideol %>% 
  filter(article_fc_rating == "Fake news") %>% 
  group_by(total_article_number, ideology_bin) %>% 
  summarise(count = sum(count, na.rm = TRUE)) %>% 
  ggplot(., aes(x = ideology_bin, y = count, fill = ideology_bin)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = seq(-3, 6, 3)) +
  scale_y_continuous(labels = comma) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-2, 2), oob = scales::squish) +
  xlab("User ideology") +
  ylab("Avg. users exposed to article") +
  theme_ctokita() +
  theme(legend.position = "none") +
  facet_wrap(~total_article_number)
gg_fake_expos
ggsave(gg_fake_expos, filename = paste0(outpath, "fake_news_exposure_by_article.pdf"), width = 180, height = 180, units = "mm", dpi = 400)
