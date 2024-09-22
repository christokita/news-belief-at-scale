#########################################
# Name: `plot_networkmetrics.R`
# Author: Chris Tokita
# Purpose: Plot network metrics of all article retweet networks.
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
# `<data storage location>/data_derived/networks/article_network_metrics.csv: network-level metrics of each article's retweet network.
# `<data storage location>/data_derived/networks/`: directory containing individual retweet network edge and node lists.
# 
# Data Out: Plots written to output sub-folder depending on if we are comparing article veracity or news source type. 
# `output/networks/veracity/`
# `output/networks/source_type/`
# 
# Machine: Chris' laptop
########################################


####################
# Load packages
####################
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
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
exposure_path <- paste0(DATA_DIRECTORY, "data_derived/exposure/estimated_users_exposed_over_time.csv") #estimated exposure per tweet
network_metric_path <- paste0(DATA_DIRECTORY, "data_derived/networks/article_network_metrics.csv") #article retweet network metrics
retweet_network_path <- paste0(DATA_DIRECTORY, "data_derived/networks/")

# Set path for plots
outpath <- 'output/networks/'
if (GROUPING == "article_fc_rating") {
  subdir_out <- 'veracity/'
} else if(GROUPING == "source_type") {
  subdir_out <- 'source_type/'
}

# Set color palette
plot_color <- "#495867"
grouping_pal <- c("#F18805", plot_color)


####################
# Load data 
####################
# Read in network metric data
network_metrics <- read.csv(network_metric_path, header = TRUE)

# Read in network data
rt_edges <- read.csv(paste0(retweet_network_path, 'rtnetwork_edges.csv'), colClasses = c("Source"="character", "Target"="character"))

# Load tweets and articles for article metadata
tweets <- read.csv(tweet_path, header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>%  #discard first 10 articles from analysis
  # clean up metadata
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "False/Misleading news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream outlet", ifelse(source_type == "fringe", "Fringe outlet", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )
article_data <- tweets %>% 
  select(user_id, tweet_id, total_article_number, source_type, source_lean, article_fc_rating, article_lean, user_ideology) 

# Read in exposure data
exposure_data <- read.csv(exposure_path, 
                          header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
  mutate(tweet_number = tweet_number+1) %>%  #python zero index
  rename(time = relative_time) %>% 
  arrange(total_article_number, tweet_number) %>% 
  merge(article_data, by = c("user_id", "tweet_id")) %>% 
  # clean up metadata
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "False/Misleading news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream outlet", ifelse(source_type == "fringe", "Fringe outlet", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )
# If analyzing by veracity, drop out non-True/False articles
if (GROUPING == "article_fc_rating") {
  exposure_data <- exposure_data %>% 
    filter(article_fc_rating %in% c("True news", "False/Misleading news"))
}



############################## Basic Summary Statistics ##############################

####################
# PLOT: Tweet types
####################
# Breakdown of tweet types
gg_tweettypes <- 
  # Prep data
  tweets %>% 
  mutate(is_retweet = as.logical(is_retweet),
         is_quote = as.logical(is_quote),
         tweet_type = ifelse(is_retweet|is_quote, "Retweet", "Original Share"),
         tweet_type = factor(tweet_type, levels = c("Original Share", "Retweet"))) %>% 
  # Plot
  ggplot(., aes(x = tweet_type, fill = tweet_type)) +
  geom_bar(stat = "count") +
  scale_y_continuous(breaks = seq(0, 160000, 20000), 
                     limits = c(0, 120000),
                     expand = c(0, 0),
                     labels = scales::comma) +
  scale_fill_manual(values =  c("#081d58", "#225ea8")) +
  xlab("News article tweet type") +
  ylab("Count") +
  theme_ctokita() +
  theme(legend.position = "none")

gg_tweettypes
ggsave(gg_tweettypes, filename = paste0(outpath, "tweet_type.pdf"), width = 45, heigh = 45, units = "mm", dpi = 400)
  

# Breakdown of retweet types
# NOTE: users can (a) retweet someone directly, (b) retweet someone indirectly, by retweeting a retweet, (c) retweet themself, or (d) retweet by quote tweeting a tweet
gg_RTtypes <- 
  # Prep data
  rt_edges %>% 
  mutate(RT_type = ifelse(grepl("Phantom RT", RT_type), "Direct RT", RT_type),
         RT_type = gsub(" RT", "", RT_type),
         RT_type = factor(RT_type, levels = c("Direct", "Indirect", "Self", "Quote"))) %>% 
  # Plot
  ggplot(., aes(x = RT_type, fill = RT_type)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("#1d91c0", "#41b6c4", "#c7e9b4", "#edf8b1")) +
  scale_y_continuous(breaks = seq(0, 100000, 20000), 
                     limits = c(0, 100000),
                     expand = c(0, 0),
                     labels = scales::comma) +
  xlab("Retweet type") +
  ylab("Count") +
  theme_ctokita() +
  theme(legend.position = "none")

gg_RTtypes
ggsave(gg_RTtypes, filename = paste0(outpath, "RT_type.pdf"), width = 45, heigh = 45, units = "mm", dpi = 400)


####################
# PLOT: Degree distribution, RAW
####################
# Degree frequency data
degree_data <- exposure_data %>% 
  select(user_id, !!sym(GROUPING), follower_count) %>% 
  distinct() %>%
  group_by(!!sym(GROUPING), follower_count) %>% 
  summarise(count = n()) %>% 
  arrange(desc(follower_count)) %>% 
  mutate(prob = count / sum(count),
         cdf = cumsum(count) / sum(count)) %>%
  rename(degree = follower_count) %>% 
  arrange(!!sym(GROUPING), degree)

# Plot
gg_degree_dist <- degree_data %>% 
  filter(degree > 0) %>% 
  ggplot(aes(x = degree, y = prob, color = !!sym(GROUPING))) +
  geom_point(size = 0.8, alpha = 0.3, stroke = 0) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = c(10^(-1:-8)),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^-5, 10^-2)) +
  scale_color_manual(values = grouping_pal) +
  ylab("P(N. followers)") +
  xlab("Number of followers") +
  theme_ctokita() +
  theme(legend.title = element_blank(),
        legend.position = c(0.6, 0.95),
        legend.key.size = unit(2.5, "mm"),
        legend.key.height = unit(0, 'mm'),
        legend.background = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 1.25)))

gg_degree_dist
ggsave(gg_degree_dist, filename = paste0(outpath, subdir_out, "degreedistribution.pdf"), width = 45, height = 45, units = "mm")

                
####################
# PLOT:Retweet distribution, RAW
####################
# Retweet frequency data
if (GROUPING == "article_fc_rating") {
  filtered_tweets <- tweets %>% 
    filter(article_fc_rating %in% c("True news", "False/Misleading news"))
} else {
  filtered_tweets <- tweets
}

RT_freq_data <- filtered_tweets %>% 
  select(total_article_number, !!sym(GROUPING)) %>% 
  distinct() %>% 
  merge(rt_edges, by = "total_article_number") %>% 
  group_by(!!sym(GROUPING), total_article_number, Source) %>% 
  summarise(n_retweets = n()) %>%  #count number of retweets each user has
  ungroup() %>% 
  group_by(!!sym(GROUPING), n_retweets) %>% 
  summarise(count = n()) %>% 
  arrange(desc(n_retweets)) %>% 
  mutate(prob = count / sum(count),
         cdf = cumsum(count) / sum(count)) %>%
  arrange(!!sym(GROUPING), n_retweets)

# Plot
gg_retweetfreq_dist <- RT_freq_data %>% 
  ggplot(aes(x = n_retweets, y = prob, color = !!sym(GROUPING))) +
  geom_point(size = 0.8, alpha = 0.5, stroke = 0) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = c(10^(-1:-8)),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^-4.1, 10^-0)) +
  scale_color_manual(values = c("#F18805", plot_color)) +
  ylab("P(N. rewteets)") +
  xlab("Number of retweets") +
  theme_ctokita() +
  theme(legend.title = element_blank(),
        legend.position = c(0.6, 0.95),
        legend.key.size = unit(2.5, "mm"),
        legend.key.height = unit(0, 'mm'),
        legend.background = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 1.25)))

gg_retweetfreq_dist
ggsave(gg_retweetfreq_dist, filename = paste0(outpath, subdir_out, "retweetdistribution.pdf"), width = 45, height = 45, units = "mm")



############################## Ideology in networks ##############################

####################
# PLOT: Ideological diversity by article veracity
####################
# Bayesian point estimate
# Note: we assume equal variances. Both an eye check and running this with unequal variances confirm they cannot be distinguished.
blm_ideoldiversity <- brm(bf(ideology_sd ~ 0 + article_fc_rating), 
                          data = network_metrics %>% filter(article_fc_rating %in% c("T", "FM")),
                          family = gaussian(), 
                          chains = 4, 
                          warmup = 1000, 
                          iter = 3500)

ideoldiversity_estimates <- as.data.frame( posterior_summary(blm_ideoldiversity, 
                                                           variables = c("article_fc_ratingFM", "article_fc_ratingT")) ) %>% 
  tibble::rownames_to_column() %>% 
  rename(article_fc_rating = rowname) %>%
  mutate(article_fc_rating = gsub("b_article_fc_rating", "", article_fc_rating)) %>% 
  select(-Q2.5, -Q97.5)

# Merge in HDI-based CI
ideoldiversity_estimates <- posterior_samples(blm_ideoldiversity) %>% 
  select(b_article_fc_ratingFM, b_article_fc_ratingT) %>% 
  bayestestR::hdi(., ci = 0.99) %>% 
  as.data.frame() %>% 
  rename(article_fc_rating = Parameter) %>% 
  mutate(article_fc_rating = gsub("b_article_fc_rating", "", article_fc_rating)) %>% 
  merge(ideoldiversity_estimates, ., by = "article_fc_rating")

# Hypothesis test that they aren't the same
hypothesis_ideoldiversity <- hypothesis(blm_ideoldiversity, "article_fc_ratingFM = article_fc_ratingT", alpha = 0.05)
hypothesis_ideoldiversity #BF >> 100, P = 1

t.test(ideology_sd ~ article_fc_rating, data = network_metrics %>% filter(article_fc_rating %in% c("T", "FM"))) #p = 0.03955

# Plot
gg_ideodiversity <- network_metrics %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
           ggplot(., aes(x = article_fc_rating, color = article_fc_rating)) +
  # geom_point(aes(y = ideology_sd),
  #            size = 1.2, stroke = 0, alpha = 0.2,
  #            position = position_jitter(width = 0.05)) +
  ggbeeswarm::geom_quasirandom(aes(y = ideology_sd),
                               size = 1.5, stroke = 0, alpha = 0.3,
                               width = 0.2) +
  geom_errorbar(data = ideoldiversity_estimates, aes(ymin = CI_low, ymax = CI_high), 
                size = 0.6, width = 0) +
  geom_point(data = ideoldiversity_estimates, aes(y = Estimate),
             size = 2.5, stroke = 0) +
  scale_x_discrete(labels = c("False/Misleading", "True")) +
  scale_y_continuous(breaks = seq(0, 2, 0.5),
                     limits = c(0, 2),
                     expand = c(0, 0)) +
  scale_color_manual(values = grouping_pal) +
  ylab("Ideol. diversity of tweeters") +
  xlab("Article veracity") +
  theme_ctokita() +
  theme(legend.position = "none")

gg_ideodiversity
ggsave(gg_ideodiversity, filename = paste0(outpath, "ideodiversity_byveracity.pdf"), width = 45, height = 45, units = "mm")


############################## Network structure ##############################

####################
# PLOT: Network density by article veracity
####################
gg_veracitydensity <- network_metrics %>% 
  filter(article_fc_rating %in% c("FM", "T")) %>% 
  # filter(total_tweets > 20) %>% 
  ggplot(., aes(x = article_fc_rating, y = network_density)) +
  geom_point(size = 1, stroke = 0, alpha = 0.5,
             position = position_jitter(width = 0.03)) +
  scale_x_discrete(labels = c("False/Misleading", "True")) +
  ylab("Network density") +
  xlab("Article veracity") +
  theme_ctokita()

gg_veracitydensity
ggsave(gg_veracitydensity, filename = paste0(outpath, "networkdensity_byveracity.pdf"), width = 45, height = 45, units = "mm", dpi = 400)
