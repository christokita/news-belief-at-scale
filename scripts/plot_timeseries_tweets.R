########################################
#
# PLOT: Tweets over time, comparing by article veracity or by news source type
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
# Paramters for analysis: grouping of interest, paths to data, paths for output, and filename
####################
# Choose grouping of interest. Options: 
#     (1) article veracity: "article_fc_rating"
#     (2) source: "source_type"
grouping <- "source_type"

# Paths to files/directories
tweet_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/tweets_labeled.csv' #path to fitness cascade data
if (grouping == "article_fc_rating") {
  outpath <- 'output/tweet_timeseries/veracity/'
} else if(grouping == "source_type") {
  outpath <- 'output/tweet_timeseries/source_type/'
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
# Read in data, Calculate time since first sharing of the story
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

# Add dummy rows of pre-first share for plotting purposes
# (1) Create empty dataframe for dummy rows
n_articles <- length(unique(tweets$total_article_number)) #number of unique articles
dummy_rows <- data.frame(matrix(NA, ncol = ncol(tweets), nrow = 2*n_articles))  #create empty dataframe
names(dummy_rows) <- names(tweets) #give same column names
# (2) Create unique set of article IDs and fact-check rating to add to our dummy rows
unique_article_ratings <- tweets %>% 
  select(source_type, source_lean, total_article_number, article_fc_rating) %>% 
  unique()
# (3) Join together
tweets <- dummy_rows %>% 
  select(-article_fc_rating, -source_type, -source_lean) %>% 
  mutate(relative_tweet_time = rep(c(-2, -0.01), n_articles),
         tweet_number = rep(c(-1, 0), n_articles),
         relative_tweet_count = 0,
         total_article_number = rep(unique(tweets$total_article_number), each = 2)) %>% 
  merge(unique_article_ratings, by = "total_article_number") %>% 
  rbind(tweets, .) %>% 
  mutate(hour_bin = cut(relative_tweet_time, breaks = seq(-2, 50, 1), include.lowest = TRUE, right = FALSE, labels = seq(-2, 49))) %>%  #bin by hour tweet appeared
  mutate(hour_bin = as.numeric(as.character(hour_bin))) #convert from factor to plain number

# If analyzing by veracity, drop out non-True/False articles
if (grouping == "article_fc_rating") {
  tweets <- tweets %>% 
    filter(article_fc_rating %in% c("T", "FM"))
}

# Clean up some labels
tweets <- tweets %>% 
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "Fake news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream", ifelse(source_type == "fringe", "Fringe", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )

# Calculate new tweets per hour
tweet_perhour <- tweets %>% 
  group_by(!!sym(grouping), total_article_number, hour_bin) %>% 
  count(.)
tweet_perhour$n[tweet_perhour$hour_bin %in% c(-2, -1)] <- 0 #zero out the count of dummy rows



############################## Plot simple time series of tweet count ##############################

####################
# New tweets per hour over time
####################
gg_tweettime <- tweet_perhour %>% 
  ggplot(., aes(x = hour_bin, y = n, group = total_article_number)) +
  geom_line(size = 0.2, alpha = 0.5, color = line_color) +
  xlab("Time since first article share (hrs)") +
  ylab("Tweets") +
  scale_y_log10() +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right")
gg_tweettime


####################
# Cumulative tweets
####################
gg_totaltweets <- tweets %>% 
  ggplot(., aes(x = relative_tweet_time, y = tweet_number, group = total_article_number)) +
  geom_vline(aes(xintercept = 24), 
             linetype = "dotted",
             size = 0.3,
             color = "grey60") +
  geom_step(size = 0.2, alpha = 0.5, color = line_color) +
  scale_y_continuous(breaks = c(10^seq(1, 5)),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     trans = scales::pseudo_log_trans(base = 10)) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  xlab("Time since first article share (hrs)") +
  ylab("Log total tweets") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right")
gg_totaltweets
ggsave(gg_totaltweets, filename = paste0(outpath,"total_tweets_time.png"), width = 90, height = 45, units = "mm", dpi = 400)

# Color by other factor
if (grouping == "article_fc_rating") {
  other_grouping <- "source_type"
  plot_tag <- "sourcetype"
} else if (grouping == "source_type") {
  other_grouping <- "article_fc_rating"
  plot_tag <- "veracity"
}
gg_totaltweets_color <- tweets %>% 
  ggplot(., aes(x = relative_tweet_time, y = tweet_number, group = total_article_number, color = !!sym(other_grouping))) +
  geom_vline(aes(xintercept = 24), 
             linetype = "dotted",
             size = 0.3,
             color = "grey60") +
  geom_step(size = 0.2, alpha = 0.5) +
  scale_y_continuous(breaks = c(10^seq(1, 5)),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     trans = scales::pseudo_log_trans(base = 10)) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_color_manual(values = c("#EF8354", "#b80d48", "#EF8354", "#404040"), 
                     name = "") +
  xlab("Time since first article share (hrs)") +
  ylab("Log total tweets") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right")
gg_totaltweets_color
ggsave(gg_totaltweets_color, filename = paste0(outpath,"total_tweets_time_", plot_tag, ".png"), width = 90, height = 45, units = "mm", dpi = 400)

####################
# Percentiage of tweets per story
####################
gg_perctweets <- tweets %>% 
  ggplot(., aes(x = relative_tweet_time, y = relative_tweet_count, group = total_article_number)) +
  geom_vline(aes(xintercept = 24), 
             linetype = "dotted",
             size = 0.3,
             color = "grey60") +
  geom_step(size = 0.2, alpha = 0.5, color = line_color) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  xlab("Time since first article share (hrs)") +
  ylab("Proportion of story stweets") +
  theme_ctokita() +
  theme(aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right")
gg_perctweets
ggsave(gg_perctweets, filename = paste0(outpath, "relative_tweets_time.png"), width = 90, height = 45, units = "mm", dpi = 400)


####################
# Saturation time of stories
####################
percentiles <- c(0.5, 0.8)
for (percentile in percentiles) {
  
  percentile_data <- tweets %>% 
    filter(relative_tweet_count >= percentile) %>% 
    group_by(total_article_number) %>% 
    filter(relative_tweet_count == min(relative_tweet_count)) 
  
  gg_saturationcount <- ggplot(percentile_data, aes(x = relative_tweet_time)) +
    geom_histogram(aes(y = stat(density)), binwidth = 2, color = 'white', fill = line_color) +
    xlab(paste0("Time to ", percentile*100, "% sharing saturation (hrs.)")) +
    theme_ctokita() +
    facet_wrap(as.formula(paste("~", grouping)), 
               ncol = 1,
               strip.position = "right")
  gg_saturationcount
  ggsave(gg_saturationcount, filename = paste0(outpath, "story_saturation", percentile*100, ".png"), width = 55, height = 90, units = "mm", dpi = 400)
  
}  

# t.test(relative_tweet_time~article_fc_rating, data = percentile_data)

####################
# New tweets vs retweets over time
####################
gg_retweets <- tweets %>% 
  # data processing
  filter(hour_bin >= 0) %>% 
  group_by(!!sym(grouping), total_article_number, hour_bin) %>% 
  summarize(RTs = sum(as.logical(is_retweet)),
            total_tweets = length(is_retweet)) %>% 
  mutate(perc_RTs = RTs / total_tweets) %>% 
  group_by(!!sym(grouping), hour_bin) %>% 
  summarise(mean_perc_RTs = mean(perc_RTs),
            sd_perc_RTs = sd(perc_RTs)) %>% 
  # graph
  ggplot(., aes(x = hour_bin, y = mean_perc_RTs)) +
  geom_ribbon(aes(ymin = ifelse((mean_perc_RTs - sd_perc_RTs) < 0, 0, mean_perc_RTs - sd_perc_RTs),
                  ymax = ifelse((mean_perc_RTs + sd_perc_RTs) > 1, 1, mean_perc_RTs + sd_perc_RTs)),
              fill = line_color, alpha = 0.2) +
  geom_line(color = line_color) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  # scale_y_continuous(limits = c(0, 1)) +
  xlab("Time since first article share (hrs)") +
  ylab("Prop. retweets") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm")) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right")
gg_retweets
ggsave(gg_retweets, filename = paste0(outpath, "RTpercentage.png"), width = 55, height = 45, units = "mm", dpi = 400)
  


############################## Plot time series of tweeter ideology ##############################

####################
# Ideological distance of tweeters relative to article content
####################
# Calculate distance between ideological category of user and article
# Filter out users without ideological scores, bin by hour
gg_ideoldisttime <- tweets %>% 
  #data processing
  filter(!is.na(user_ideology)) %>% 
  mutate(article_ideol_category = ifelse(article_lean == "L", -1, ifelse(article_lean == "C", 1, 0))) %>% 
  mutate(ideol_diff = user_ideol_category - article_ideol_category) %>% 
  mutate(ideol_distance = sign(user_ideol_category - article_ideol_category)) %>%
  group_by(!!sym(grouping), total_article_number, hour_bin) %>% 
  count(ideol_distance) %>% 
  mutate(freq_ideol_distance = n / sum(n)) %>% 
  group_by(!!sym(grouping), hour_bin, ideol_distance) %>% 
  summarise(freq_ideol_distance = mean(freq_ideol_distance)) %>% 
  filter(hour_bin >= 0) %>%
  #graph
  ggplot(., aes(x = hour_bin, y = freq_ideol_distance, fill = factor(ideol_distance, levels = c(1, 0, -1)))) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_fill_manual(values = rev(ideol_dist_pal[c(1,3,5)]),
                    name = "User ideology relative to\narticle content",
                    labels = c("User more conservative", 
                               "Same ideology", 
                               "User more liberal")) +
  xlab("Time since first article share (hrs)") +
  ylab("Prop. of tweeters") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm")) +
  facet_wrap(as.formula(paste("~", grouping)), 
             strip.position = "right",
             ncol = 1)
gg_ideoldisttime
ggsave(gg_ideoldisttime, filename = paste0(outpath, "relative_ideology_dist.png"), width = 100, height = 45, units = "mm", dpi = 400)


####################
# Ideology distribution of tweeters over time
####################
# Calculate average ideological distribution of tweeters over time
# Bin ideologies, filter out users without ideological scores, bin by hour
gg_ideoltime <- tweets %>% 
  #data processing
  filter(!is.na(user_ideology)) %>% 
  mutate(ideol_bin = cut(user_ideology, breaks = c(-6, -1, 1, 6), labels = c(-1, 0, 1), right = FALSE, include.lowest = TRUE)) %>%
  group_by(!!sym(grouping), total_article_number, hour_bin) %>% 
  count(ideol_bin) %>% 
  mutate(freq_ideol_bin = n / sum(n)) %>% 
  group_by(!!sym(grouping), hour_bin, ideol_bin) %>% 
  summarise(freq_ideol_bin = mean(freq_ideol_bin)) %>% 
  # graph
  ggplot(., aes(x = hour_bin, y = freq_ideol_bin, fill = factor(ideol_bin, levels = c(1, 0, -1)))) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_fill_manual(values = rev(ideol_pal[c(1,3,5)]),
                    name = "Tweeter ideology",
                    labels = c("Conservative", "Moderate", "Liberal")) +
  xlab("Time since first article share (hrs)") +
  ylab("Prop. of tweeters") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm")) +
  facet_wrap(as.formula(paste("~", grouping)),
             ncol = 1,
             strip.position = "right")
gg_ideoltime
ggsave(gg_ideoltime, filename = paste0(outpath,"ideology_dist.png"), width = 90, height = 45, units = "mm", dpi = 400)

# Fine scale breakdown
gg_ideoltime_fine <- tweets %>% 
  #data processing
  filter(!is.na(user_ideology)) %>% 
  mutate(ideol_bin = cut(user_ideology, breaks = c(-6, -1, 1, 6), labels = c(-1, 0, 1), right = FALSE, include.lowest = TRUE)) %>%
  group_by(!!sym(grouping), article_lean, total_article_number, hour_bin) %>% 
  count(ideol_bin) %>% 
  mutate(freq_ideol_bin = n / sum(n)) %>% 
  group_by(!!sym(grouping), article_lean, hour_bin, ideol_bin) %>% 
  summarise(freq_ideol_bin = mean(freq_ideol_bin)) %>% 
  # graph
  ggplot(., aes(x = hour_bin, y = freq_ideol_bin, fill = factor(ideol_bin, levels = c(1, 0, -1)))) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_fill_manual(values = rev(ideol_pal[c(1,3,5)]),
                    name = "Tweeter ideology",
                    labels = c("Conservative", "Moderate", "Liberal")) +
  xlab("Time since first article share (hrs)") +
  ylab("Prop. of tweeters") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm")) +
  facet_grid(as.formula(paste(grouping, "~", "article_lean")))
gg_ideoltime_fine
ggsave(gg_ideoltime_fine, filename = paste0(outpath,"ideology_dist_articlelean.png"), width = 120, height = 45, units = "mm", dpi = 400)


####################
# Raw plot of ideology of tweeters over time
####################
gg_ideoltime_raw <- tweets %>% 
  filter(!is.na(user_ideology)) %>%
  ggplot(., aes(x = relative_tweet_time, y = user_ideology, color = user_ideology)) +
  geom_point(size = 0.3, stroke = 0, position = position_jitter(width = 0.1, height = 0.1)) +
  scale_color_gradientn(colors = ideol_pal, limits = c(-2, 2), oob = scales::squish) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  xlab("Time since first article share (hrs)") +
  ylab("User ideology") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.position = "none") +
  facet_wrap(as.formula(paste("~", grouping)),
             ncol = 1,
             strip.position = "right")
gg_ideoltime_raw
ggsave(gg_ideoltime_raw, filename = paste0(outpath, "ideology_raw.png"), width = 90, height = 45, units = "mm", dpi = 400)

# Break out additionally by article lean
gg_ideoltime_raw <- tweets %>% 
  filter(!is.na(user_ideology)) %>%
  ggplot(., aes(x = relative_tweet_time, y = user_ideology, color = user_ideology)) +
  geom_point(size = 0.3, stroke = 0, position = position_jitter(width = 0.1, height = 0.1)) +
  scale_color_gradientn(colors = ideol_pal, limits = c(-2, 2), oob = scales::squish) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  xlab("Time since first article share (hrs)") +
  ylab("User ideology") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.position = "none") +
  facet_grid(as.formula(paste(grouping, "~", "article_lean")))
gg_ideoltime_raw
ggsave(gg_ideoltime_raw, filename = paste0(outpath, "ideology_raw_articlelean.png"), width = 120, height = 45, units = "mm", dpi = 400)


####################
# Ideology of tweeters: source type by article veracity
####################
# Calculate average ideological distribution of tweeters over time, broken out by source and article veracity
# Bin ideologies, filter out users without ideological scores, bin by hour
gg_ideoltimesource <- tweets %>% 
  filter(!is.na(user_ideology),
         article_fc_rating %in% c("Fake news", "True news")) %>% 
  mutate(ideol_bin = cut(user_ideology, breaks = c(-6, -1, 1, 6), labels = c(-1, 0, 1), right = FALSE, include.lowest = TRUE)) %>%
  group_by(article_fc_rating, source_type, article_fc_rating, total_article_number, hour_bin) %>% 
  count(ideol_bin) %>% 
  mutate(freq_ideol_bin = n / sum(n)) %>% 
  group_by(article_fc_rating, source_type, article_fc_rating, hour_bin, ideol_bin) %>% 
  summarise(freq_ideol_bin = mean(freq_ideol_bin)) %>% 
  ggplot(., aes(x = hour_bin, y = freq_ideol_bin, fill = factor(ideol_bin, levels = c(1, 0, -1)))) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_fill_manual(values = rev(ideol_pal[c(1,3,5)]),
                    name = "User ideology",
                    labels = c("Conservative", "Moderate", "Liberal")) +
  xlab("Time since first article share (hrs)") +
  ylab("Prop. of tweeters") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm")) +
  facet_grid(article_fc_rating~source_type)
gg_ideoltimesource
ggsave(gg_ideoltimesource, filename = paste0(outpath, "ideology_dist_bysourceandveracity.png"), width = 120, height = 45, units = "mm", dpi = 400)



