#########################################
# Name: `plot_belief_insight_against_baseline.R`
# Author: Chris Tokita
# Purpose: Plots to demonstrate reviewer's request that we show how receptivity gives us more insight over baseline exposure.
# Details:
#   (These R scripts assume the use of the `.Rproj` at top of the news-belief-at-scale/ repo. Otherwise, set the working directory to one level above this script.)
#
#   The Variables at the beginning of the script that are in all caps need to be set by the user:
#     `DATA_DIRECTORY`: path to the data directory. (Copies of data are currently stored on external hard drive and high-performance cluster.)
#     `GROUPING`:       determines whether the plots will break out tweet belief according to article veracity ("article_fc_rating") or the source of the article ("source_type").
# 
# Data In:
# `<data storage location>/data_derived/tweets/tweets_labeled.csv`: article tweets with article and tweeter metadata.
# `<data storage location>/data_derived/belief/estimated_belief_over_time.csv`: estimated belief in each article tweet.
# `<data storage location>/data/articles/evaluations.csv`: article fact-check ratings for each unique news article.
# `<data storage location>/data/articles/daily_articles.csv`: article source information for each unique news article.
# 
# Data Out: Plots written to output sub-folder depending on if we are comparing article veracity or news source type. 
# `output/belief/veracity/`
# `output/belief/source_type/`
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
library(brms)
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
belief_path <- paste0(DATA_DIRECTORY, "data_derived/belief/estimated_belief_over_time.csv") #estimated belief per tweet
article_evaluation_path <- paste0(DATA_DIRECTORY, "data/articles/evaluations.csv") #article fact-check rating
article_source_path <- paste0(DATA_DIRECTORY, "data/articles/daily_articles.csv") #article source type

# Set path for plots
outpath <- 'output/belief/'
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

plot_color <- "#495867"
grouping_pal <- c("#F18805", plot_color)


####################
# Load and prepare data 
####################
# Read in tweet data for article info
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

# Get article metadata and count number of articles per GROUPING (useful for average distributions)
article_data <- tweets %>% 
  select(tweet_id, total_article_number, source_type, source_lean, article_fc_rating, article_lean, user_ideology) 

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

# Load belief data 
belief_data <- read.csv(belief_path, header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
  mutate(tweet_number = tweet_number+1) %>%  #python zero index
  rename(time = relative_time) %>% 
  arrange(total_article_number, tweet_number)

# Fix rounding error for new belief
# Some cases have one more new_believing_user than new_exposed_user due to rounding
belief_data$new_believing_users[belief_data$new_believing_users > belief_data$new_exposed_users] = belief_data$new_exposed_users[belief_data$new_believing_users > belief_data$new_exposed_users]

# Merge in relevant article level data
# NOTE: adds one extra row, check after double checking with new data
belief_timeseries <- merge(belief_data, article_data, by = c("tweet_id", "total_article_number"), all = TRUE) 

# Add dummy rows of pre-first share for plotting purposes
# NOTE: This will add rows to the dataframe---two extra rows per article---but it will not effect user counts
# (1) Create empty dataframe for dummy rows
n_articles <- length(unique(belief_timeseries$total_article_number)) #number of unique articles
dummy_rows <- data.frame(matrix(NA, ncol = ncol(belief_timeseries), nrow = 2*n_articles))  #create empty dataframe
names(dummy_rows) <- names(belief_timeseries) #give same column names
# (2) Create unique set of article IDs and fact-check rating to add to our dummy rows
unique_article_ratings <- article_data %>% 
  select(source_type, source_lean, total_article_number, article_fc_rating) %>% 
  unique()
# (3) Join together
belief_timeseries <- dummy_rows %>% 
  select(-article_fc_rating, -source_type, -source_lean) %>% 
  mutate(time = rep(c(-2, -0.01), n_articles),
         tweet_number = rep(c(-1, 0), n_articles),
         new_exposed_users = 0, 
         cumulative_exposed = 0, 
         new_believing_users = 0,
         cumulative_believing = 0,
         total_article_number = rep(unique(article_data$total_article_number), each = 2)) %>% 
  merge(unique_article_ratings, by = "total_article_number") %>% 
  rbind(belief_timeseries, .) %>% 
  mutate(hour_bin = cut(time, breaks = seq(-2, 24*14, 1), include.lowest = TRUE, right = FALSE, labels = seq(-2, 24*14-1))) %>%  #bin by hour tweet appeared
  mutate(hour_bin = as.numeric(as.character(hour_bin))) %>%  #convert from factor to plain number
  group_by(total_article_number) %>% 
  mutate(relative_cumulative_exposed = cumulative_exposed / max(cumulative_exposed),
         relative_cumulative_belief = cumulative_believing / max(cumulative_believing),
         relative_tweet_count = tweet_number / max(tweet_number)) %>% 
  arrange(total_article_number, tweet_number) %>% 
  ungroup()

rm(dummy_rows, article_data, belief_data, article_group_counts)

# If analyzing by veracity, drop out non-True/False articles
if (GROUPING == "article_fc_rating") {
  belief_timeseries <- belief_timeseries %>% 
    filter(article_fc_rating %in% c("T", "FM"))
}

# Clean up some labels
belief_timeseries <- belief_timeseries %>% 
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "False/Misleading news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream outlet", ifelse(source_type == "fringe", "Fringe outlet", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )

# Melt data to make one ideological bin per tweet per row
# This is helpful for plotting belief by ideological bin
belief_ideol <- belief_timeseries %>% 
  filter(hour_bin >= 0) %>% 
  select(-source_lean, -relative_cumulative_exposed, -relative_cumulative_belief, -relative_tweet_count, -follower_count) %>% 
  gather(key = "ideology_bin", value = "count", 
         -time, -tweet_number, -tweet_id, -user_id, -user_ideology, -new_exposed_users, -cumulative_exposed, -new_believing_users, -cumulative_believing, -total_article_number, -hour_bin, -source_type, -article_fc_rating, -article_lean, -n_articles_in_grouping) %>% 
  mutate(ideology_bin = gsub("ideol_", "", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("^\\.", "-", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("_\\.", "_-", ideology_bin)) %>% 
  separate(ideology_bin, c("lower", "upper"), sep = "_", convert = TRUE) %>% 
  mutate(ideology_bin = (lower + upper) / 2) %>% 
  select(-lower, -upper)


####################
# PLOT: Relative cumulative exposure vs belief over first 48 hours
####################
# Prep time bins
missing_time_bins <- expand_grid(total_article_number = unique(belief_timeseries$total_article_number), 
                                 time_bin = seq(0, 24, 0.1))

# Plot
gg_baseline <- belief_timeseries %>% 
  # Remove articles that didn't expose anyone
  filter(!is.na(relative_cumulative_belief)) %>% 
  # Bin into 6 min increments
  mutate(time_bin = (time %/% 0.1) * 0.1) %>% 
  # Add in missing hour bins
  filter(time <= 24) %>% 
  merge(missing_time_bins, by = c("total_article_number", "time_bin"), all = TRUE) %>% 
  fill(cumulative_exposed, cumulative_believing, !!sym(GROUPING)) %>% 
  select(total_article_number, !!sym(GROUPING), time_bin, cumulative_exposed, cumulative_believing) %>%
  # Prep data
  group_by(article_fc_rating, total_article_number, time_bin) %>% 
  summarise(cumulative_exposed = max(cumulative_exposed, na.rm = TRUE),
         cumulative_believing = max(cumulative_believing, na.rm = TRUE)) %>% 
  group_by(!!sym(GROUPING), time_bin) %>% 
  summarise(cumulative_exposed = mean(cumulative_exposed, na.rm = TRUE),
            cumulative_believing = mean(cumulative_believing, na.rm = TRUE)) %>% 
  gather(key = "metric", value = "count", -!!sym(GROUPING), -time_bin) %>% 
  # Plot
  ggplot(., aes(x = time_bin, y = count, color = !!sym(GROUPING), fill = !!sym(GROUPING))) +
  geom_area(aes(alpha = metric), 
            position = "identity", 
            linewidth = 0,
            color = NA) +
  xlab("Time since first article share (hrs.)") +
  ylab("Cumulative users") +
  scale_x_continuous(breaks = seq(0, 24, 6),
                     limits = c(-2, 24), 
                     expand = c(0, 0)) +
  scale_y_continuous(labels = comma,
                     expand = c(0, 0)) +
  scale_color_manual(values = grouping_pal) +
  scale_fill_manual(values = grouping_pal) +
  scale_alpha_manual(values = c(0.65, 0.35)) +
  theme_ctokita() +
  theme(aspect.ratio = NULL,
        legend.position = "none") +
  facet_wrap(as.formula(paste("~", GROUPING)),
             ncol = 1,
             strip.position = "top",
             scales = "free")

gg_baseline
ggsave(gg_baseline, filename = paste0(outpath, subdir_out, "cumulative_users_exposed_vs_receptive.pdf"), width = 45, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Relative cumulative exposure vs belief over first 48 hours
####################
# Prep legend labels for plot
if (GROUPING == "article_fc_rating") {
  legend_name <- "Article rating"
  legend_labels <- c("False/Misleading", "True")
} else if (GROUPING == "source_type") {
  legend_name <- "News Source Type"
  legend_labels <- c("Fringe", "Mainstream")
} 

# Plot
gg_baseline <- belief_timeseries %>% 
  # Remove articles that didn't expose anyone
  filter(!is.na(relative_cumulative_belief)) %>% 
  filter(time >= 0 & time <= 72) %>% 
  # Prep data
  mutate(belief_above_exposure = relative_cumulative_belief - relative_cumulative_exposed) %>% 
  # Plot
  ggplot(., aes(x = time, y = belief_above_exposure, color = !!sym(GROUPING), fill = !!sym(GROUPING), group = total_article_number)) +
  geom_line(linewidth = 0.6) +
  xlab("Time since first tweet (hrs)") +
  ylab("Belief accumulation relative to exposure") +
  scale_color_manual(values = grouping_pal, name = legend_name, labels = legend_labels) +
  scale_fill_manual(values = grouping_pal, name = legend_name, labels = legend_labels) +
  theme_ctokita() +
  theme(legend.position = c(0.75, 0.2))

gg_baseline
