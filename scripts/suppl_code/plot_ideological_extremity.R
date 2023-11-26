#########################################
# Name: `plot_ideological_extremity.R`
# Author: Chris Tokita
# Purpose: Plot ideological extremity of receptive and non-receptive users
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
# Load belief data
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
rm(exposure_data, exposure_timeseries, tweets, article_data)



######################################## Ideological Extremity of Receptive and Non-Receptive Users ########################################

####################
# Calculate non-receptive users among those exposed
####################
# Sort dataframes to make sure they lineup
exposure_ideol <- exposure_ideol %>% 
  arrange(article_fc_rating, total_article_number, time, tweet_number, tweet_id, ideology_bin)

belief_ideol <- belief_ideol %>% 
  arrange(article_fc_rating, total_article_number, time, tweet_number, tweet_id, ideology_bin)

# Calculate nonreceptive users
receptivity_ideol <- exposure_ideol 
receptivity_ideol$count_belief <- belief_ideol$count

receptivity_ideol <- receptivity_ideol %>% 
  select(article_fc_rating, tweet_id, total_article_number, time, tweet_number, ideology_bin, count, count_belief) %>% 
  rename(count_exposed = count) %>% 
  mutate(count_nonbelief = count_exposed - count_belief) # nonreceptive = exposed - receptive

# Create absolute ideology bin, i.e., ideological extremity
receptivity_ideol <- receptivity_ideol %>% 
  mutate(absolute_ideology_bin = abs(ideology_bin + 0.25)) # Shift bins labels from left edge to center to allow for transformation to absolute value


####################
# PLOT: Comparing normalized distribution of receptive vs nonreceptive users exposed for true and false news
####################
# Prep data
receptivity_ideol_dist <- receptivity_ideol %>% 
  group_by(article_fc_rating, absolute_ideology_bin) %>% 
  summarise("Receptive users" = sum(count_belief, na.rm = TRUE),
            "Non-receptive users" = sum(count_nonbelief, na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_longer(-c(article_fc_rating, absolute_ideology_bin), names_to = "receptivity", values_to = "count") %>% 
  group_by(article_fc_rating, receptivity) %>% 
  mutate(proportion = count / sum(count)) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news")) %>% 
  arrange(article_fc_rating, receptivity, absolute_ideology_bin)

# Plot
gg_receptivity_ideological_extreme_compare <- ggplot(receptivity_ideol_dist, aes(x = absolute_ideology_bin, y = proportion, fill = article_fc_rating)) +
  geom_step(aes(color = article_fc_rating, absolute_ideology_bin - 0.25),
            linewidth = 0.3,
            alpha = 0.8) +
  geom_bar(position = "identity",
           stat = "identity",
           alpha = 0.5,
           width = 0.5)  +
  scale_x_continuous(limits = c(-0.25, 6), 
                     expand = c(0, 0), 
                     breaks = seq(0, 6, 1)) +
  scale_y_continuous(limits = c(0, 0.4),
                     expand = c(0, 0), 
                     labels = comma) +
  scale_fill_manual(values = grouping_pal) +
  scale_color_manual(values = grouping_pal) +
  xlab("Relative user ideology") +
  ylab("Proportion of users exposed") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(~receptivity, 
             ncol = 1,
             strip.position = "top",
             scales = "free")

gg_receptivity_ideological_extreme_compare
# ggsave(gg_receptivity_ideological_extreme_compare, filename = "output/belief/veracity/ideol_extremity.pdf", width = 45, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Comparing normalized distribution of true and false news by receptivity
####################
# Prep data
receptivity_ideol_dist <- receptivity_ideol %>% 
  group_by(article_fc_rating, absolute_ideology_bin) %>% 
  summarise("Receptive users" = sum(count_belief, na.rm = TRUE),
            "Non-receptive users" = sum(count_nonbelief, na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_longer(-c(article_fc_rating, absolute_ideology_bin), names_to = "receptivity", values_to = "count") %>% 
  group_by(article_fc_rating) %>% 
  mutate(proportion = count / sum(count)) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news")) %>% 
  arrange(article_fc_rating, receptivity, absolute_ideology_bin)

# Plot
gg_veracity_ideological_extreme_compare <- ggplot(receptivity_ideol_dist, aes(x = absolute_ideology_bin, y = count, fill = article_fc_rating)) +
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
  ylab("Users exposed") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "top",
             scales = "free")

gg_veracity_ideological_extreme_compare
ggsave(gg_veracity_ideological_extreme_compare, filename = "output/belief/veracity/ideol_extremity.pdf", width = 45, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Comparing normalized distribution of exposed users for true and false news
####################
# Prep data
exposure_ideol_dist <- receptivity_ideol %>% 
  group_by(article_fc_rating, absolute_ideology_bin) %>% 
  summarise(exposed_count = sum(count_exposed, na.rm = TRUE)) %>% 
  mutate(proportion = exposed_count / sum(exposed_count)) %>% 
  ungroup() %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news")) 

# Plot
gg_exposure_ideological_extreme_compare <- ggplot(exposure_ideol_dist, aes(x = absolute_ideology_bin, y = proportion, fill = article_fc_rating)) +
  # Data
  geom_step(aes(color = article_fc_rating, absolute_ideology_bin - 0.25),
            linewidth = 0.3,
            alpha = 0.8) +
  geom_bar(position = "identity",
           stat = "identity",
           alpha = 0.5,
           width = 0.5)  +
  scale_x_continuous(limits = c(-0.25, 6), 
                     expand = c(0, 0), 
                     breaks = seq(0, 6, 1)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, 0.4),
                     labels = comma) +
  scale_fill_manual(values = grouping_pal) +
  scale_color_manual(values = grouping_pal) +
  xlab("Relative user ideology") +
  ylab("Proportion of users exposed") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL)

gg_exposure_ideological_extreme_compare
ggsave(gg_exposure_ideological_extreme_compare, filename = "output/exposure/veracity/ideol_extremity.pdf", width = 45, height = 84.5, units = "mm", dpi = 400)






####################
# STATS: Test whether ideological extremity is significantly different for exposure
####################
# Is total user exposure different by news veracity?
exposure_t <- exposure_ideol_dist %>% 
  filter(article_fc_rating == "True news") %>% 
  arrange(absolute_ideology_bin)

exposure_fm <- exposure_ideol_dist %>% 
  filter(article_fc_rating == "False/Misleading news") %>% 
  arrange(absolute_ideology_bin)

chisq.test(
  data.frame(nonreceptive = exposure_t$exposed_count, receptive = exposure_fm$exposed_count) #p<0.0001
)

# Is exposure different between non-receptive and receptive users for a given news veracity
receptive_fm <- receptivity_ideol_dist %>% 
  filter(article_fc_rating == "False/Misleading news", receptivity == "Receptive users") %>% 
  arrange(absolute_ideology_bin)
receptive_t <- receptivity_ideol_dist %>% 
  filter(article_fc_rating == "True news", receptivity == "Receptive users") %>% 
  arrange(absolute_ideology_bin)

chisq.test(
  data.frame(true = receptive_t$count, false = receptive_fm$count) #p<0.0001
)


nonreceptive_t <- receptivity_ideol_dist %>% 
  filter(article_fc_rating == "True news", receptivity == "Non-receptive users") %>% 
  arrange(absolute_ideology_bin)
nonreceptive_fm <- receptivity_ideol_dist %>% 
  filter(article_fc_rating == "False/Misleading news", receptivity == "Non-receptive users") %>% 
  arrange(absolute_ideology_bin)

chisq.test(
  data.frame(true = nonreceptive_t$count, false = nonreceptive_fm$count) #p<0.0001
)


# Do receptive and non-receptive users look different for true news and for false news?
chisq.test(
  data.frame(nonreceptive = nonreceptive_t$count, receptive = receptive_t$count) #p<0.0001
)

chisq.test(
  data.frame(nonreceptive = nonreceptive_fm$count, receptive = receptive_fm$count) #p<0.0001
)

# Do receptive users look different from overall exposure pattern
chisq.test(
  data.frame(exposure = exposure_t$exposed_count, receptive = receptive_t$count) #p<0.0001
)

chisq.test(
  data.frame(exposure = exposure_fm$exposed_count, receptive = receptive_fm$count) #p<0.0001
)

# Do non-receptive users look different from overall exposure pattern
chisq.test(
  data.frame(exposure = exposure_t$exposed_count, nonreceptive = nonreceptive_t$count) #p<0.0001
)

chisq.test(
  data.frame(exposure = exposure_fm$exposed_count, nonreceptive = nonreceptive_fm$count) #p<0.0001
)
