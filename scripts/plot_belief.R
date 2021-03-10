########################################
#
# PLOT: Estimated belief to news articles over time, comparing by article veracity or by news source type
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
grouping <- "article_fc_rating"

# Paths to files/directories
tweet_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/tweets_labeled.csv' #path to fitness cascade data
if (grouping == "article_fc_rating") {
  outpath <- 'output/belief/veracity/'
} else if(grouping == "source_type") {
  outpath <- 'output/belief/source_type/'
}

# Color palette
line_color <- "#495867"
ideol_pal <- rev(brewer.pal(5, "RdBu"))
ideol_pal[3] <- "#e0e0e0"
ideol_dist_pal <- rev(brewer.pal(5, "PuOr"))
ideol_dist_pal[3] <- "#e0e0e0"
ideol_limit <- 3 #limit beyond which we squish the color palette


####################
# Load and prep data 
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
# NOTE: adds one extra row, check after double checking with new data
belief_timeseries <- merge(belief_data, article_data, by = c("tweet_id", "total_article_number"), all = TRUE) 

# Add dummy rows of pre-first share for plotting purposes
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
         total_article_number = rep(unique(article_data$total_article_number), each = 2)) %>% 
  merge(unique_article_ratings, by = "total_article_number") %>% 
  rbind(belief_timeseries, .) %>% 
  mutate(hour_bin = cut(time, breaks = seq(-2, 50, 1), include.lowest = TRUE, right = FALSE, labels = seq(-2, 49))) %>%  #bin by hour tweet appeared
  mutate(hour_bin = as.numeric(as.character(hour_bin))) %>%  #convert from factor to plain number
  group_by(total_article_number) %>% 
  mutate(relative_cumulative_exposed = cumulative_exposed / max(cumulative_exposed),
         relative_tweet_count = tweet_number / max(tweet_number)) %>% 
  arrange(total_article_number, tweet_number)

rm(dummy_rows, article_data, belief_data)

# If analyzing by veracity, drop out non-True/False articles
if (grouping == "article_fc_rating") {
  belief_timeseries <- belief_timeseries %>% 
    filter(article_fc_rating %in% c("T", "FM"))
}

# Clean up some labels
belief_timeseries <- belief_timeseries %>% 
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "False/Misleading news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream", ifelse(source_type == "fringe", "Fringe", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )



############################## Plot belief by ideology ##############################

# Prep data
belief_ideol <- belief_timeseries %>% 
  filter(hour_bin >= 0) %>% 
  select(-source_lean, -relative_cumulative_exposed, -relative_tweet_count, -follower_count) %>% 
  gather(key = "ideology_bin", value = "count", 
         -time, -tweet_number, -tweet_id, -user_id, -user_ideology, -new_exposed_users, -cumulative_exposed, -total_article_number, -hour_bin, -source_type, -article_fc_rating, -article_lean) %>% 
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
# Total ideologies believing articles (by veracity)
####################
gg_ideol_total <- belief_ideol %>% 
  # Calculate totals
  group_by(!!sym(grouping), ideology_bin) %>% 
  summarise(count = sum(count, na.rm = TRUE),
            avg_count = sum(count, na.rm = TRUE) / unique(n_articles_in_grouping)) %>% 
  # Plot
  ggplot(., aes(x = ideology_bin, y = count, fill = ideology_bin, color = ideology_bin)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 2)) +
  scale_y_continuous(labels = comma,
                     expand = c(0, 0)) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("Believing user ideology") +
  ylab("Total users") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "free")
gg_ideol_total
ggsave(gg_ideol_total, filename = paste0(outpath, "ideol_total_belief.pdf"), width = 45, height = 90, units = "mm", dpi = 400)


####################
# Avg number of believing article across ideologies (by type)
####################
gg_ideol_avg <- belief_ideol %>% 
  # Calculate average exposed
  group_by(!!sym(grouping), ideology_bin) %>% 
  summarise(avg_count = sum(count, na.rm = TRUE) / unique(n_articles_in_grouping)) %>% 
  # Plot
  ggplot(., aes(x = ideology_bin, y = avg_count, fill = ideology_bin, color = ideology_bin)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 1)) +
  scale_y_continuous(labels = comma,
                     expand = c(0, 0)) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +  
  xlab("Believing user ideology") +
  ylab("Avg. number of users") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "free")
gg_ideol_avg
ggsave(gg_ideol_avg, filename = paste0(outpath, "ideol_avg_belief.pdf"), width = 90, height = 90, units = "mm", dpi = 400)


####################
# Ideological distributions of belief
####################
gg_ideol_dist <- belief_ideol %>% 
  # Calculate average distribution of belief
  group_by(!!sym(grouping), ideology_bin, total_article_number, n_articles_in_grouping) %>% 
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  group_by(total_article_number) %>% 
  mutate(belief_prop = count / sum(count),
         belief_prop = ifelse( is.na(belief_prop), 0, belief_prop)) %>% 
  group_by(!!sym(grouping), ideology_bin) %>% 
  summarise(avg_belief_prop = sum(belief_prop) / unique(n_articles_in_grouping)) %>% 
  # Plot
  ggplot(., aes(x = ideology_bin, y = avg_belief_prop, fill = ideology_bin, color = ideology_bin)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.05), 
                     limits = c(0, 0.25),
                     expand = c(0, 0)) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_fill_gradientn(colours = ideol_pal, limit = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  xlab("Believing user ideology") +
  ylab("Avg. proportion of article beliefs") +
  theme_ctokita() +
  theme(legend.position = "none",
        aspect.ratio = NULL) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "free_x")
gg_ideol_dist

ggsave(gg_ideol_dist, filename = paste0(outpath, "ideol_avg_belief_distribution.pdf"), width = 45, height = 90, units = "mm", dpi = 400)


####################
# Belief time series
####################
gg_ideoltime <- belief_ideol %>% 
  group_by(!!sym(grouping), hour_bin, ideology_bin) %>% 
  summarise(count = sum(count)) %>%
  ggplot(., aes(x = hour_bin, y = count, fill = ideology_bin, color = ideology_bin)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  scale_fill_gradientn(colours = ideol_pal, 
                       name = "User\nideology",
                       limits = c(-2, 2), 
                       oob = squish) +
  scale_color_gradientn(colours = ideol_pal, 
                       name = "User\nideology",
                       limits = c(-2, 2), 
                       oob = squish) +
  scale_x_continuous(breaks = seq(0, 48, 6),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Time since first article share (hrs)") +
  ylab("New users believing") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_rect(size = 0.6, fill = NA),
        axis.line = element_blank()) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "free")
gg_ideoltime
ggsave(gg_ideoltime, filename = paste0(outpath, "ideol_belief_hourbin.pdf"), width = 90, height = 90, units = "mm", dpi = 400)


# gg_ideoltime_binned <-  exposure_timeseries %>% 
#   mutate(Liberal = ideol_.3.0_.2.5 + ideol_.2.5_.2.0 + ideol_.2.0_.1.5 + ideol_.1.5_.1.0,
#          Moderate = ideol_.1.0_.0.5 + ideol_.0.5_0.0 + ideol_0.0_0.5 + ideol_0.5_1.0,
#          Conservative = ideol_1.0_1.5 + ideol_1.5_2.0 + ideol_2.0_2.5 + ideol_2.5_3.0+ ideol_3.0_3.5 + ideol_3.5_4.0 + ideol_4.0_4.5 + ideol_4.5_5.0 + ideol_5.0_5.5) +
#   select()
