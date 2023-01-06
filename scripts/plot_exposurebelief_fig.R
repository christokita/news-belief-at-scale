#########################################
# Name: `plot_exposurebelief_fig.R`
# Author: Chris Tokita
# Purpose: Plot estimated exposure to and belief in news article for main figures for paper.
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
# 
# Data Out: Plots written to output subfolder focusing on belief (since we are analyzing both belief and exposure combined here) broken out by article veracity. 
# `<data storage location>/output/belief/veracity/`
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
         source_type = ifelse(source_type == "mainstream", "Mainstream", ifelse(source_type == "fringe", "Fringe", source_type)),
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
  select(-lower, -upper) %>% 
  # count number of articles per grouping (useful for average distributions)
  group_by(article_fc_rating) %>% 
  mutate(n_articles_in_grouping = length(unique(total_article_number)))
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
         source_type = ifelse(source_type == "mainstream", "Mainstream", ifelse(source_type == "fringe", "Fringe", source_type)),
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
  select(-lower, -upper) %>% 
  # count number of articles per grouping (useful for average distributions)
  group_by(article_fc_rating) %>% 
  mutate(n_articles_in_grouping = length(unique(total_article_number)))
rm(exposure_data, exposure_timeseries, tweets, article_data)



######################################## Total Exposure/Belief ########################################

####################
# PLOT: Total exposure and belief, broken out by article veracity
####################
# Prep data
exposure_ideol_sum <- exposure_ideol %>% 
  group_by(article_fc_rating, ideology_bin) %>% 
  summarise(exposure_count = sum(count, na.rm = TRUE))

belief_ideol_sum <- belief_ideol %>% 
  group_by(article_fc_rating, ideology_bin) %>% 
  summarise(belief_count = sum(count, na.rm = TRUE))

exposure_belief_data <- merge(exposure_ideol_sum, belief_ideol_sum, by = c("article_fc_rating", "ideology_bin")) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news"))

# Plot
# NOTE: labels for "Exposed" and "Believing" are added later in vector art program during figure creation for paper
gg_exposebelief_total <- ggplot(exposure_belief_data, aes(x = ideology_bin, fill = ideology_bin)) +
  # Data
  geom_step(aes(x = ideology_bin - 0.25, y = exposure_count,  color = ideology_bin - 0.25),
           size = 0.3,
           alpha = 0.8) +
  geom_bar(aes(y = exposure_count),
           stat = "identity",
           alpha = 0.5,
           width = 0.5) +
  geom_step(aes(x = ideology_bin - 0.25, y = belief_count),
            size = 0.6,
            color = "white",
            alpha = 0.7) +
  geom_bar(aes(y = belief_count),
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
  ylab("Total number of users") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(~article_fc_rating, 
             ncol = 1,
             strip.position = "top",
             scales = "free")

gg_exposebelief_total
ggsave(gg_exposebelief_total, filename = "output/belief/veracity/combined_total_beliefANDexposure.pdf", width = 45, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Total exposure and belief, broken out by both (a) article veracity and (b) news source lean
####################
# Prep data
exposure_ideol_sum <- exposure_ideol %>% 
  group_by(article_fc_rating, source_lean, ideology_bin) %>% 
  summarise(exposure_count = sum(count, na.rm = TRUE))

belief_ideol_sum <- belief_ideol %>% 
  group_by(article_fc_rating, source_lean, ideology_bin) %>% 
  summarise(belief_count = sum(count, na.rm = TRUE))

exposure_belief_data <- merge(exposure_ideol_sum, belief_ideol_sum, by = c("article_fc_rating", "source_lean", "ideology_bin")) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news")) %>% 
  mutate(source_lean = ifelse(source_lean == "C", "Conservative", 
                              ifelse(source_lean == "L", "Liberal", 
                                     ifelse(source_lean == "U", "Unclear", "")))) %>% 
  mutate(source_lean = factor(source_lean, levels = c("Liberal", "Unclear", "Conservative")))

# Plot
# NOTE:
# - axis/facet label for "News Source Lean" and "Exposed"/"Believing" are added later in vector art program during figure creation for paper
# - the "Exposed" count is the lighter histogram, while the "Believing" count is the darker histogram
gg_exposebelief_total_source <- ggplot(exposure_belief_data, aes(x = ideology_bin, fill = ideology_bin)) +
  # Data
  geom_step(aes(x = ideology_bin - 0.25, y = exposure_count,  color = ideology_bin - 0.25),
            size = 0.3,
            alpha = 0.8) +
  geom_bar(aes(y = exposure_count),
           stat = "identity",
           alpha = 0.5,
           width = 0.5) +
  geom_step(aes(x = ideology_bin - 0.25, y = belief_count),
            size = 0.6,
            color = "white",
            alpha = 0.7) +
  geom_bar(aes(y = belief_count),
           stat = "identity",
           alpha = 1,
           width = 0.5) +
  #Plot params
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 3)) +
  scale_y_continuous(expand = c(0, 0), 
                     labels = comma) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("User ideology") +
  ylab("Total number of users") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL,
        panel.spacing = unit(5, "mm")) +
  facet_grid(article_fc_rating~source_lean, scale = 'free')

gg_exposebelief_total_source
ggsave(gg_exposebelief_total_source, filename = "output/belief/veracity/combined_total_beliefANDexposure_bysourcelean.pdf", width = 90, height = 80, units = "mm", dpi = 400)


####################
# PLOT: Total exposure and belief, broken out by both (a) article veracity and (b) article lean
####################
# Prep data
exposure_ideol_sum <- exposure_ideol %>% 
  mutate(article_lean = gsub("Neutral|Unclear", "Neutral/Unclear", article_lean)) %>% 
  mutate(article_lean = factor(article_lean, levels = c("Liberal", "Neutral/Unclear", "Conservative"))) %>% 
  group_by(article_fc_rating, article_lean, ideology_bin) %>% 
  summarise(exposure_count = sum(count, na.rm = TRUE))

belief_ideol_sum <- belief_ideol %>% 
  mutate(article_lean = gsub("Neutral|Unclear", "Neutral/Unclear", article_lean)) %>% 
  mutate(article_lean = factor(article_lean, levels = c("Liberal", "Neutral/Unclear", "Conservative"))) %>% 
  group_by(article_fc_rating, article_lean, ideology_bin) %>% 
  summarise(belief_count = sum(count, na.rm = TRUE))

exposure_belief_data <- merge(exposure_ideol_sum, belief_ideol_sum, by = c("article_fc_rating", "article_lean", "ideology_bin")) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news"))

# Plot
# NOTE: 
# - axis/facet label for "Article Lean" and "Exposed"/"Believing" are added later in vector art program during figure creation for paper
# - the "Exposed" count is the lighter histogram, while the "Believing" count is the darker histogram
gg_exposebelief_total_article <- ggplot(exposure_belief_data, aes(x = ideology_bin, fill = ideology_bin)) +
  # Data
  geom_step(aes(x = ideology_bin - 0.25, y = exposure_count,  color = ideology_bin - 0.25),
            size = 0.3,
            alpha = 0.8) +
  geom_bar(aes(y = exposure_count),
           stat = "identity",
           alpha = 0.5,
           width = 0.5) +
  geom_step(aes(x = ideology_bin - 0.25, y = belief_count),
            size = 0.6,
            color = "white",
            alpha = 0.7) +
  geom_bar(aes(y = belief_count),
           stat = "identity",
           alpha = 1,
           width = 0.5) +
  #Plot params
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 3)) +
  scale_y_continuous(expand = c(0, 0),
                     labels = comma) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("User ideology") +
  ylab("Total number of users") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL,
        panel.spacing = unit(5, "mm")) +
  facet_grid(article_fc_rating~article_lean, scale = 'free_y')

gg_exposebelief_total_article
ggsave(gg_exposebelief_total_article, filename = "output/belief/veracity/combined_total_beliefANDexposure_byarticlelean.pdf", width = 90, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Difference in ideological composition of exposure and belief, broken out by article veracity
####################
# Prep data
exposure_ideol_sum <- exposure_ideol %>% 
  group_by(article_fc_rating, ideology_bin) %>% 
  summarise(exposure_count = sum(count, na.rm = TRUE))

belief_ideol_sum <- belief_ideol %>% 
  group_by(article_fc_rating, ideology_bin) %>% 
  summarise(belief_count = sum(count, na.rm = TRUE))

exposure_belief_data <- merge(exposure_ideol_sum, belief_ideol_sum, by = c("article_fc_rating", "ideology_bin")) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news"))


# Calculate cumulative distribution of exposure and belief
exposure_belief_diff <- exposure_belief_data %>% 
  mutate(ideology_bin = ideology_bin + 0.25) %>%  #bins in original data represent left edge
  group_by(article_fc_rating) %>% 
  arrange(article_fc_rating, ideology_bin) %>% 
  mutate(exposure_count = replace_na(exposure_count, 0),
         belief_count = replace_na(belief_count, 0), 
         exposure_pct = exposure_count / sum(exposure_count),
         belief_pct = belief_count / sum(belief_count),
         diff_pct = belief_pct - exposure_pct,
         ratio_pct = belief_pct / exposure_pct - 1,
         exposure_cumulative = cumsum(exposure_count),
         belief_cumulative = cumsum(belief_count),
         exposure_cdf = exposure_cumulative / max(exposure_cumulative),
         belief_cdf = belief_cumulative / max(belief_cumulative),
         diff_cdf = belief_cdf - exposure_cdf)

# Plot ideology bin % of exposure vs % of belief
gg_pct_exposurebelief <- ggplot(exposure_belief_diff, aes(x = exposure_pct, y = belief_pct, color = ideology_bin)) +
  geom_abline(aes(slope = 1, intercept = 0), size = 0.3, linetype = "dotted") +
  geom_point(size = 3, stroke = 0) +
  # scale_color_manual(values = grouping_pal) +
  scale_x_continuous(breaks = seq(0, 0.2, 0.1), limits = c(0, 0.21)) +
  scale_y_continuous(breaks = seq(0, 0.2, 0.1), limits = c(0, 0.21)) +
  scale_color_gradientn(colors = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish, name = "User ideology") +
  scale_shape_manual(values = c(19, 0)) +
  xlab("Proportion of total exposed users") +
  ylab("Propotion of total believing users") +
  facet_wrap(~article_fc_rating, ncol = 1, scales = "free") +
  theme_ctokita() 

gg_pct_exposurebelief
ggsave(gg_pct_exposurebelief, filename = "output/belief/veracity/belief_vs_exposure/belief_vs_exposure.pdf", height = 90, width = 70, units = "mm", dpi = 400)

# % of right-leaning users exposed to true news vs % believing
exposure_belief_diff %>% 
  filter(ideology_bin > 0) %>% 
  group_by(article_fc_rating) %>% 
  summarise(pct_exposed = sum(exposure_pct),
            pct_belief = sum(belief_pct))

####################
# PLOT: Difference in ideological composition of exposure and belief, broken out by (a) article veracity and (b) article lean
####################
# Prep data
exposure_ideol_sum <- exposure_ideol %>% 
  mutate(article_lean = gsub("Neutral|Unclear", "Neutral/Unclear", article_lean)) %>% 
  mutate(article_lean = factor(article_lean, levels = c("Liberal", "Neutral/Unclear", "Conservative"))) %>% 
  group_by(article_fc_rating, article_lean, ideology_bin) %>% 
  summarise(exposure_count = sum(count, na.rm = TRUE))

belief_ideol_sum <- belief_ideol %>% 
  mutate(article_lean = gsub("Neutral|Unclear", "Neutral/Unclear", article_lean)) %>% 
  mutate(article_lean = factor(article_lean, levels = c("Liberal", "Neutral/Unclear", "Conservative"))) %>% 
  group_by(article_fc_rating, article_lean, ideology_bin) %>% 
  summarise(belief_count = sum(count, na.rm = TRUE))

exposure_belief_data <- merge(exposure_ideol_sum, belief_ideol_sum, by = c("article_fc_rating", "article_lean", "ideology_bin")) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news"))

# Calculate cumulative distribution of exposure and belief
exposure_belief_diff <- exposure_belief_data %>% 
  mutate(ideology_bin = ideology_bin + 0.25) %>%  #bins in original data represent left edge
  group_by(article_fc_rating, article_lean) %>% 
  arrange(article_fc_rating, ideology_bin) %>% 
  mutate(exposure_count = replace_na(exposure_count, 0),
         belief_count = replace_na(belief_count, 0), 
         exposure_pct = exposure_count / sum(exposure_count),
         belief_pct = belief_count / sum(belief_count),
         diff_pct = belief_pct - exposure_pct,
         ratio_pct = belief_pct / exposure_pct - 1,
         exposure_cumulative = cumsum(exposure_count),
         belief_cumulative = cumsum(belief_count),
         exposure_cdf = exposure_cumulative / max(exposure_cumulative),
         belief_cdf = belief_cumulative / max(belief_cumulative),
         diff_cdf = belief_cdf - exposure_cdf)

# Plot ideology bin % of exposure vs % of belief
gg_pct_exposurebelief_articlelean <- ggplot(exposure_belief_diff, aes(x = exposure_pct, y = belief_pct, color = ideology_bin)) +
  geom_abline(aes(slope = 1, intercept = 0), size = 0.3, linetype = "dotted") +
  geom_point(size = 3, stroke = 0) +
  scale_x_continuous(breaks = seq(0, 0.3, 0.1), limits = c(0, 0.31), trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(breaks = seq(0, 0.3, 0.1), limits = c(0, 0.31), trans = scales::pseudo_log_trans(base = 10)) +
  scale_color_gradientn(colors = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish, name = "User ideology") +
  scale_shape_manual(values = c(19, 0)) +
  xlab("Proportion of total exposed users") +
  ylab("Propotion of total believing users") +
  facet_grid(article_fc_rating~article_lean, scales = "free") +
  theme_ctokita() +
  theme(panel.spacing = unit(5, "mm"))

gg_pct_exposurebelief_articlelean
ggsave(gg_pct_exposurebelief_articlelean, filename = "output/belief/veracity/belief_vs_exposure/belief_vs_exposure_articlelean.pdf", height = 90, width = 180, units = "mm", dpi = 400)


  
######################################## Exposure/Belief Over Time ########################################

####################
# PLOT: Total exposure and belief, broken out by hours 0-6, 6-12, 12-18, 18-24, and 24+
####################
# Prep data
exposure_ideol_sum <- exposure_ideol %>% 
  mutate(time_hr_bin = floor(time / 2)) %>% 
  mutate(time_window = ifelse(time_hr_bin == 0, "First 2 hrs.", 
                              ifelse(time_hr_bin == 1, "Hr. 2 to 4", 
                                     ifelse(time_hr_bin == 2, "Hr. 4 to 6", 
                                            ifelse(time_hr_bin == 3, "Hr. 6 to 8", 
                                                   "Hr. 8+"))))) %>% 
  group_by(article_fc_rating, time_window, ideology_bin) %>% 
  summarise(exposure_count = sum(count, na.rm = TRUE))

belief_ideol_sum <- belief_ideol %>% 
  mutate(time_hr_bin = floor(time / 2)) %>% 
  mutate(time_window = ifelse(time_hr_bin == 0, "First 2 hrs.", 
                              ifelse(time_hr_bin == 1, "Hr. 2 to 4", 
                                     ifelse(time_hr_bin == 2, "Hr. 4 to 6", 
                                            ifelse(time_hr_bin == 3, "Hr. 6 to 8", 
                                                   "Hr. 8+"))))) %>% 
  group_by(article_fc_rating, time_window, ideology_bin) %>% 
  summarise(belief_count = sum(count, na.rm = TRUE))

exposure_belief_data <- merge(exposure_ideol_sum, belief_ideol_sum, by = c("article_fc_rating", "ideology_bin", "time_window")) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news")) %>% 
  mutate(time_window = factor(time_window, levels = c("First 2 hrs.", "Hr. 2 to 4", "Hr. 4 to 6", "Hr. 6 to 8", "Hr. 8+")))

# Plot
gg_exposebelief_total_timewindow <- ggplot(exposure_belief_data, aes(x = ideology_bin, fill = ideology_bin)) +
  # Data
  geom_step(aes(x = ideology_bin - 0.25, y = exposure_count,  color = ideology_bin - 0.25),
            size = 0.3,
            alpha = 0.8) +
  geom_bar(aes(y = exposure_count),
           stat = "identity",
           alpha = 0.5,
           width = 0.5) +
  geom_step(aes(x = ideology_bin - 0.25, y = belief_count),
            size = 0.6,
            color = "white",
            alpha = 0.7) +
  geom_bar(aes(y = belief_count),
           stat = "identity",
           alpha = 1,
           width = 0.5) +
  #Plot params
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 2),
                     labels = c("-6", "", "", "0", "", "", "6")) +
  scale_y_continuous(expand = c(0, 0), 
                     labels = comma) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("User ideology") +
  ylab("Total number of users") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL,
        panel.spacing.x = unit(0.4, "lines")) +
  facet_grid(article_fc_rating~time_window, scale = "free")

gg_exposebelief_total_timewindow
ggsave(gg_exposebelief_total_timewindow, filename = "output/belief/veracity/combined_total_beliefANDexposure_timewindow.pdf", width = 90, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Time to 50 percent exposure/belief
####################
# While this section of plots/analysis is not directly used in the paper,
# it was an analysis to see if there was any difference with the rate at which True and False/Misleading news articles
# accumulate user exposure and user belief. 
# We find there is not a difference. We may opt to delete this after the review if it is not relevant.

# Prep data
majority_exposure <- belief_timeseries %>% 
  group_by(total_article_number) %>% 
  mutate(relative_exposed = cumulative_exposed / max(cumulative_exposed),
         majority_exposed = relative_exposed >= 0.5) %>% 
  arrange(total_article_number, time, tweet_id) %>% 
  filter(majority_exposed == TRUE) %>% 
  filter(time == min(time),
         relative_exposed == min(relative_exposed)) %>% 
  distinct(time, relative_exposed, .keep_all = TRUE) %>% 
  rename(time_majority_exposure = time) %>% 
  select(total_article_number, time_majority_exposure, article_fc_rating, relative_exposed)
  

majority_belief <- belief_timeseries %>% 
  group_by(total_article_number) %>% 
  mutate(relative_believing = cumulative_believing / max(cumulative_believing),
         majority_believing = relative_believing >= 0.5) %>% 
  arrange(total_article_number, time, tweet_id) %>% 
  filter(majority_believing == TRUE) %>% 
  filter(time == min(time),
         relative_believing == min(relative_believing)) %>% 
  distinct(time, relative_believing, .keep_all = TRUE) %>% 
  rename(time_majority_believing = time) %>% 
  select(total_article_number, time_majority_believing, relative_believing)


majority_point <- merge(majority_exposure, majority_belief, by = c("total_article_number")) %>% 
  gather(key = "type", value = "time", 
         -total_article_number, -article_fc_rating, -relative_exposed, -relative_believing) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news")) %>% 
  mutate(article_fc_rating = factor(article_fc_rating, levels = c("True news", "False/Misleading news")),
         type = ifelse(type == "time_majority_exposure", "Exposure",
                       ifelse(type == "time_majority_believing", "Belief", NA)),
         type = factor(type, levels = c("Exposure", "Belief")))

# Calculate mean
mean_majority_point <- majority_point %>% 
  group_by(article_fc_rating, type) %>% 
  summarise(mean_time = mean(time))

# Calculate statistics
F_belief <- majority_point$time[majority_point$type == "Belief" & majority_point$article_fc_rating == "False/Misleading news"]
T_belief <- majority_point$time[majority_point$type == "Belief" & majority_point$article_fc_rating == "True news"]
t.test(T_belief, F_belief)

F_exposure <- majority_point$time[majority_point$type == "Exposure" & majority_point$article_fc_rating == "False/Misleading news"]
T_exposure <- majority_point$time[majority_point$type == "Exposure" & majority_point$article_fc_rating == "True news"]
t.test(T_exposure, F_exposure)

# Plot
gg_majority_point <- ggplot(data = majority_point, aes(x = time, y = ..density.., group = article_fc_rating, fill = article_fc_rating)) + 
  geom_histogram(bins = 10,
                 alpha = 0.7,
                 position = "identity") +
  geom_vline(data = mean_majority_point, aes(xintercept = mean_time, color = article_fc_rating),
             size = 0.75, linetype = "dotted") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(1, 10, 100, 1000),
                     labels = scales::comma_format(accuracy = 1)) +
  scale_fill_manual(name = "Article rating", labels = c("True", "False/Misleading"),
                    values = grouping_pal[c(2, 1)]) +
  scale_color_manual(name = "Article rating", labels = c("True", "False/Misleading"),
                    values = grouping_pal[c(2, 1)]) +
  xlab("Time to majority article exposure/belief (hrs.)") +
  ylab("Density") +
  facet_wrap(type~.,
             ncol = 2,
             scales = "free_x") +
  theme_ctokita() +
  theme(aspect.ratio = NULL,
        legend.position = "top")

gg_majority_point
ggsave(gg_majority_point, filename = "output/belief/veracity/combined_time_to_majority.pdf", width = 90, height = 100, units = "mm", dpi = 400)
