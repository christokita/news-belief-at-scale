########################################
#
# PLOT: Estimated belief and exposure for paper figure
#
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

tweet_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/tweets_labeled.csv'
ideol_pal <- rev(brewer.pal(5, "RdBu"))
ideol_pal[3] <- "#e0e0e0"
ideol_limit <- 3 #limit beyond which we squish the color palette


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
belief_data <- read.csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/belief/estimated_belief_over_time.csv', 
                        header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
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


# Add ideological data
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
rm(belief_data, belief_timeseries)


####################
# Load expsoure data
####################
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
exposure_timeseries <- merge(exposure_data, article_data, by = c("tweet_id", "total_article_number"), all = TRUE) %>% 
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "False/Misleading news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream", ifelse(source_type == "fringe", "Fringe", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )

# Add ideology data
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


####################
# PLOT: Total exposure and belief
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
# PLOT: Total exposure and belief, broken out by source lean
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

# Plot: True news
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
  facet_grid(article_fc_rating~source_lean, scale = 'free_y')
gg_exposebelief_total_source

ggsave(gg_exposebelief_total_source, filename = "output/belief/veracity/combined_total_beliefANDexposure_bysourcelean.pdf", width = 135, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Total exposure and belief, broken out by article lean
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
  facet_grid(article_fc_rating~article_lean, scale = 'free_y')
gg_exposebelief_total_article

ggsave(gg_exposebelief_total_article, filename = "output/belief/veracity/combined_total_beliefANDexposure_byarticlelean.pdf", width = 135, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Total exposure and belief, broken out by hours 0-6 vs. hours 18-24
####################
# Prep data
exposure_ideol_sum <- exposure_ideol %>% 
  filter( time < 24) %>% #grab first 24 hours
  mutate(time_hr_6 = floor(time / 6)) %>% 
  mutate(time_window = ifelse(time_hr_6 == 0, "First 6 hrs.", 
                              ifelse(time_hr_6 == 1, "Hr. 6 to 12", 
                                     ifelse(time_hr_6 == 2, "Hr. 12 to 18", 
                                            ifelse(time_hr_6 == 3, "Hr. 18 to 24", NA))))) %>% 
  group_by(article_fc_rating, time_window, ideology_bin) %>% 
  summarise(exposure_count = sum(count, na.rm = TRUE))

belief_ideol_sum <- belief_ideol %>% 
  filter( time < 24) %>% #grab first 24 hours
  mutate(time_hr_6 = floor(time / 6)) %>% 
  mutate(time_window = ifelse(time_hr_6 == 0, "First 6 hrs.", 
                              ifelse(time_hr_6 == 1, "Hr. 6 to 12", 
                                     ifelse(time_hr_6 == 2, "Hr. 12 to 18", 
                                            ifelse(time_hr_6 == 3, "Hr. 18 to 24", NA))))) %>% 
  group_by(article_fc_rating, time_window, ideology_bin) %>% 
  summarise(belief_count = sum(count, na.rm = TRUE))

exposure_belief_data <- merge(exposure_ideol_sum, belief_ideol_sum, by = c("article_fc_rating", "ideology_bin", "time_window")) %>% 
  filter(article_fc_rating %in% c("False/Misleading news", "True news")) %>% 
  mutate(time_window = factor(time_window, levels = c("First 6 hrs.", "Hr. 6 to 12", "Hr. 12 to 18", "Hr. 18 to 24")))

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
        panel.spacing.x = unit(0.7, "lines")) +
  facet_grid(article_fc_rating~time_window, scale = 'free')
gg_exposebelief_total_timewindow

ggsave(gg_exposebelief_total_timewindow, filename = "output/belief/veracity/combined_total_beliefANDexposure_timewindow.pdf", width = 90, height = 90, units = "mm", dpi = 400)
