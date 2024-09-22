#########################################
# Name: `plot_exposurebelief_normalied.R`
# Author: Chris Tokita
# Purpose: Plot estimated exposure to and belief in news article normalized across article virality.
# Details:
#   (These R scripts assume the use of the `.Rproj` at top of the news-belief-at-scale/ repo. Otherwise, set the working directory to one level above this script.)
#
#   The Variables at the beginning of the script that are in all caps need to be set by the user:
#     `DATA_DIRECTORY`: path to the data directory. (Copies of data are currently stored on external hard drive and high-performance cluster.)
# 
# Data In:
# `<data storage location>/data_derived/tweets/tweets_labeled.csv`: article tweets with article and tweeter metadata.
# `<data storage location>/data_derived/exposure/estimated_users_exposed_over_time.csv: estimated exposure to each article tweet.
# `<data storage location>/data_derived/belief/estimated_belief_over_time.csv`: estimated belief in each article tweet.
# `<data storage location>/data/articles/evaluations.csv`: article fact-check ratings for each unique news article.
# `<data storage location>/data/articles/daily_articles.csv`: article source information for each unique news article.
#
# Data Out: Plots written to output subfolder focusing on belief (since we are analyzing both belief and exposure combined here) broken out by article veracity. 
# `output/belief/veracity/`
# 
# Machine: Chris' laptop
########################################


####################
# Load packages and set paths
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


####################
# Prepare for analysis: set paths to data, paths for output, and color palettes for plotting
####################
# Set paths for data
tweet_path <- paste0(DATA_DIRECTORY, "data_derived/tweets/tweets_labeled.csv") #tweets
exposure_path <- paste0(DATA_DIRECTORY, "data_derived/exposure/estimated_users_exposed_over_time.csv") #estimated exposure per tweet
belief_path <- paste0(DATA_DIRECTORY, "data_derived/belief/estimated_belief_over_time.csv") #estimated belief per tweet
article_evaluation_path <- paste0(DATA_DIRECTORY, "data/articles/evaluations.csv") #article fact-check rating
article_source_path <- paste0(DATA_DIRECTORY, "data/articles/daily_articles.csv") #article source type

# Set path for plots
outpath <- 'output/belief/veracity/'

# Set plotting palettes
ideol_pal <- rev(brewer.pal(5, "RdBu"))
ideol_pal[3] <- "#e0e0e0"
ideol_limit <- 3 #limit beyond which we squish the color palette

plot_color <- "#495867"
grouping_pal <- c("#F18805", plot_color)


####################
# Load article data
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

article_data <- tweets %>% 
  select(tweet_id, total_article_number, source_type, source_lean, article_fc_rating, article_lean, user_ideology) 


####################
# Load belief data
####################
# Load belief data 
belief_data <- read.csv(belief_path, header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
  mutate(tweet_number = tweet_number+1) %>%  #python zero index
  rename(time = relative_time) %>% 
  arrange(total_article_number, tweet_number)

# Merge in relevant article level data
belief_timeseries <- merge(belief_data, article_data, by = c("tweet_id", "total_article_number"), all = TRUE) %>% 
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "False/Misleading news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream outlet", ifelse(source_type == "fringe", "Fringe outlet", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )

# Add ideological data (data is melted to make one ideological bin per tweet per row)
belief_ideol <- belief_timeseries %>% 
  select(-follower_count) %>% 
  gather(key = "ideology_bin", value = "count", 
         -time, -tweet_number, -tweet_id, -user_id, -user_ideology, -new_exposed_users, -new_believing_users, -cumulative_exposed, -cumulative_believing, -total_article_number, -source_type, -source_lean, -article_fc_rating, -article_lean) %>% 
  mutate(ideology_bin = gsub("ideol_", "", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("^\\.", "-", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("_\\.", "_-", ideology_bin)) %>% 
  separate(ideology_bin, c("lower", "upper"), sep = "_", convert = TRUE) %>% 
  mutate(ideology_bin = (lower + upper) / 2) %>% 
  select(-lower, -upper)
rm(belief_data)


####################
# Load exposure data
####################
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
exposure_timeseries <- merge(exposure_data, article_data, by = c("tweet_id", "total_article_number"), all = TRUE) %>% 
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "False/Misleading news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream outlet", ifelse(source_type == "fringe", "Fringe outlet", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )

# Add ideology data (data is melted to make one ideological bin per tweet per row)
exposure_ideol <- exposure_timeseries %>% 
  gather(key = "ideology_bin", value = "count", 
         -time, -tweet_number, -tweet_id, -user_id, -user_ideology, -follower_count, -new_exposed_users, -cumulative_exposed, -total_article_number, -source_type, -source_lean, -article_fc_rating, -article_lean) %>% 
  mutate(ideology_bin = gsub("ideol_", "", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("^\\.", "-", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("_\\.", "_-", ideology_bin)) %>% 
  separate(ideology_bin, c("lower", "upper"), sep = "_", convert = TRUE) %>% 
  mutate(ideology_bin = (lower + upper) / 2) %>% 
  select(-lower, -upper)
rm(exposure_data, exposure_timeseries, article_data, tweets)


####################
# Calculate normalized exposure
####################
norm_exposure_ideol <- exposure_ideol %>% 
  # For each article, determine proportion exposed by ideology bin
  group_by(article_fc_rating, ideology_bin, total_article_number) %>% 
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  group_by(total_article_number) %>% 
  mutate(exposed_prop = count / sum(count),
         exposed_prop = ifelse( is.na(exposed_prop), 0, exposed_prop)) %>% 
  # Now determine average distribution shape by article GROUPING
  group_by(article_fc_rating, ideology_bin) %>% 
  summarise(avg_exposed_prop = sum(exposed_prop) / length(unique(total_article_number))) #note this doesn't always sum to one within a grouping


####################
# Calculate normalized belief
####################
# Count total exposed and total receptive by article and ideological bin
exposure_by_article_and_ideol <- exposure_ideol %>% 
  group_by(total_article_number, article_fc_rating, ideology_bin) %>% 
  summarise(exposed = sum(count, na.rm = TRUE))

belief_by_article_and_ideol <- belief_ideol %>% 
  group_by(total_article_number, article_fc_rating, ideology_bin) %>% 
  summarise(belief = sum(count, na.rm = TRUE))

# Calculate fraction exposed that are receptive by article and ideological bin
belief_fraction_by_article_and_ideol <- exposure_by_article_and_ideol %>% 
  merge(belief_by_article_and_ideol, by = c("total_article_number", "article_fc_rating", "ideology_bin")) %>% 
  mutate(belief_rate = belief / exposed)
rm(exposure_by_article_and_ideol, belief_by_article_and_ideol)

# Calculate average belief rate per article veracity
belief_rate_avg <- belief_fraction_by_article_and_ideol %>% 
  group_by(article_fc_rating, ideology_bin) %>% 
  summarise(avg_belief_rate = mean(belief_rate, na.rm = TRUE))

# Calculate normalized belief (and exposure)
norm_belief_ideol <- norm_exposure_ideol %>% 
  merge(belief_rate_avg, by = c("article_fc_rating", "ideology_bin")) %>% 
  mutate(avg_belief_prop = avg_exposed_prop * avg_belief_rate) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news")) #only focus on T and F articles
rm(belief_fraction_by_article_and_ideol, belief_rate_avg)


######################################## Normalized Exposure/Belief ########################################

####################
# PLOT: Total exposure and belief, broken out by article veracity
####################
# Plot
# NOTE: labels for "Exposed" and "Believing" are added later in vector art program during figure creation for paper
gg_exposebelief_normalized <- ggplot(norm_belief_ideol, aes(x = ideology_bin, fill = ideology_bin)) +
  # Data
  geom_step(aes(x = ideology_bin - 0.25, y = avg_exposed_prop,  color = ideology_bin - 0.25),
            linewidth = 0.3,
            alpha = 0.8) +
  geom_bar(aes(y = avg_exposed_prop),
           stat = "identity",
           alpha = 0.5,
           width = 0.5) +
  geom_step(aes(x = ideology_bin - 0.25, y = avg_belief_prop),
            linewidth = 0.6,
            color = "white",
            alpha = 0.7) +
  geom_bar(aes(y = avg_belief_prop),
           stat = "identity",
           alpha = 1,
           width = 0.5) +
  #Plot params
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 2)) +
  scale_y_continuous(expand = c(0, 0), 
                     labels = comma) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("User ideology") +
  ylab("Avg. proportion of users") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "top",
             scales = "free")

gg_exposebelief_normalized
ggsave(gg_exposebelief_normalized, filename = paste0(outpath, "combined_normalized_beliefANDexposure.pdf"), width = 45, height = 90, units = "mm", dpi = 400)



######################################## Ideological Extremity of Normalized Exposure/Belief ########################################

####################
# Recalculate normalized exposure and belief by ideological extremity
####################
# Calculate average proportion of exposure by absolute ideological bin
norm_exposure_ideol_extr <- exposure_ideol %>% 
  mutate(absolute_ideology_bin = abs(ideology_bin + 0.25)) %>%  # Shift bins labels from left edge to center to allow for transformation to absolute value
  # For each article, determine proportion exposed by ideology bin
  group_by(article_fc_rating, absolute_ideology_bin, total_article_number) %>% 
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  group_by(total_article_number) %>% 
  mutate(exposed_prop = count / sum(count),
         exposed_prop = ifelse( is.na(exposed_prop), 0, exposed_prop)) %>% 
  # Now determine average distribution shape by article GROUPING
  group_by(article_fc_rating, absolute_ideology_bin) %>% 
  summarise(avg_exposed_prop = sum(exposed_prop) / length(unique(total_article_number))) #note this doesn't always sum to one within a grouping

# Count total exposed and total receptive by article and ideological bin
exposure_by_article_and_ideol_extr <- exposure_ideol %>% 
  mutate(absolute_ideology_bin = abs(ideology_bin + 0.25)) %>%  # Shift bins labels from left edge to center to allow for transformation to absolute value
  group_by(total_article_number, article_fc_rating, absolute_ideology_bin) %>% 
  summarise(exposed = sum(count, na.rm = TRUE))

belief_by_article_and_ideol_extr <- belief_ideol %>% 
  mutate(absolute_ideology_bin = abs(ideology_bin + 0.25)) %>%  # Shift bins labels from left edge to center to allow for transformation to absolute value
  group_by(total_article_number, article_fc_rating, absolute_ideology_bin) %>% 
  summarise(belief = sum(count, na.rm = TRUE))

# Calculate fraction exposed that are receptive by article and ideological bin
belief_fraction_by_article_and_ideol_extr <- exposure_by_article_and_ideol_extr %>% 
  merge(belief_by_article_and_ideol_extr, by = c("total_article_number", "article_fc_rating", "absolute_ideology_bin")) %>% 
  mutate(belief_rate = belief / exposed)
rm(exposure_by_article_and_ideol_extr, belief_by_article_and_ideol_extr)

# Calculate average belief rate per article veracity
belief_rate_avg_extr <- belief_fraction_by_article_and_ideol_extr %>% 
  group_by(article_fc_rating, absolute_ideology_bin) %>% 
  summarise(avg_belief_rate = mean(belief_rate, na.rm = TRUE))

# Calculate normalized belief (and exposure)
norm_belief_ideol_extr <- norm_exposure_ideol_extr %>% 
  merge(belief_rate_avg_extr, by = c("article_fc_rating", "absolute_ideology_bin")) %>% 
  mutate(belief_prop = avg_exposed_prop * avg_belief_rate,
         nonbelief_prop = avg_exposed_prop - belief_prop) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news")) #only focus on T and F articles
rm(belief_fraction_by_article_and_ideol_extr, belief_rate_avg_extr)


####################
# PLOT: Comparing normalized distribution of true and false news by receptivity
####################
# Prep data
long_norm_belief_ideol <- norm_belief_ideol_extr %>% 
  select(article_fc_rating, absolute_ideology_bin, belief_prop, nonbelief_prop) %>% 
  rename("Receptive users" = belief_prop, "Non-receptive users" = nonbelief_prop) %>% 
  pivot_longer(-c(article_fc_rating, absolute_ideology_bin), names_to = "receptivity", values_to = "avg_prop")

# Plot
gg_veracity_ideological_extreme_compare <- 
  ggplot(long_norm_belief_ideol, aes(x = absolute_ideology_bin, y = avg_prop, fill = article_fc_rating)) +
  geom_step(aes(x = absolute_ideology_bin - 0.25, linetype = receptivity, color = article_fc_rating),
            linewidth = 0.3,
            alpha = 0.8) +
  geom_bar(aes(alpha = receptivity),
           position = "identity",
           stat = "identity",
           width = 0.5)  +
  scale_x_continuous(limits = c(-0.25, 6), 
                     expand = c(0, 0), 
                     breaks = seq(0, 6, 1)) +
  scale_y_continuous(expand = c(0, 0),
                     labels = comma) +
  scale_alpha_manual(values = c(0.2, 0.4)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_color_manual(values = grouping_pal) +
  scale_fill_manual(values = grouping_pal) +
  xlab("Relative user ideology") +
  ylab("Avg. proportion of users") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "top",
             scales = "free")

gg_veracity_ideological_extreme_compare
ggsave(gg_veracity_ideological_extreme_compare, filename = paste0("output/belief/veracity/ideol_extremity_normalized.pdf"), width = 45, height = 90, units = "mm", dpi = 400)
