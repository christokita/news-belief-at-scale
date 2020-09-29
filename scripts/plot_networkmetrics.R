########################################
#
# PLOT: Network metrics of articles
#
########################################

####################
# Load packages
####################
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
source("scripts/_plot_themes/theme_ctokita.R")

####################
# Paramters for analysis: paths to data, paths for output, and filename
####################
# Choose grouping of interest. Options: 
#     (1) article veracity: "article_fc_rating"
#     (2) source: "source_type"
grouping <- "article_fc_rating"

# Paths
tweet_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/tweets_labeled.csv'
exposure_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/exposure/estimated_users_exposed_over_time.csv'
network_metric_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/networks/article_network_metrics.csv'
retweet_network_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/networks/'
outpath <- 'output/networks/'

# Path to specific subdirectory for grouping specific resultes
if (grouping == "article_fc_rating") {
  subdir_out <- 'veracity/'
} else if(grouping == "source_type") {
  subdir_out <- 'source_type/'
}


# Color palette
plot_color <- "#495867"


####################
# Load data 
####################
# Read in network metric data
network_metrics <- read.csv(network_metric_path, header = TRUE)

# Read in network data
rt_edges <- read.csv(paste0(retweet_network_path, 'rtnetwork_edges.csv'), colClasses = c("Source"="character", "Target"="character"))

# Load tweets for article metadata
article_data <- read.csv(tweet_path, header = TRUE, colClasses = c("user_id"="character", "tweet_id"="character")) %>% 
  filter(total_article_number > 10) %>% #discard first 10 articles from analysis
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
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "Fake news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream", ifelse(source_type == "fringe", "Fringe", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )
# If analyzing by veracity, drop out non-True/False articles
if (grouping == "article_fc_rating") {
  exposure_data <- exposure_data %>% 
    filter(article_fc_rating %in% c("True news", "Fake news"))
}

############################## Basics ##############################

####################
# Plot: Degree distribution
####################
degree_data <- exposure_data %>% 
  select(user_id, !!sym(grouping), follower_count) %>% 
  distinct() %>%
  group_by(!!sym(grouping), follower_count) %>% 
  summarise(count = n()) %>% 
  mutate(prob = count / sum(count))

gg_degree_dist <- degree_data %>% 
  filter(follower_count > 0) %>%
  ggplot(., aes(x = follower_count, y = prob, color = !!sym(grouping))) +
  geom_point(size = 0.8, alpha = 0.7, stroke = 0) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = c(0.01, 0.001, 0.0001, 0.00001),
                labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                limits = c(0.00001, 0.01)) +
  scale_color_manual(values = c("#F18805", plot_color)) +
  ylab("P(N. followers)") +
  xlab("Number of followers") +
  theme_ctokita() +
  theme(legend.title = element_blank(),
        legend.position = c(0.75, 0.9),
        legend.key.size = unit(2.5, "mm"),
        legend.key.height = unit(0, 'mm')) +
  guides(color = guide_legend(override.aes = list(size = 1.25)))
gg_degree_dist

ggsave(gg_degree_dist, filename = paste0(outpath, subdir_out, "degreedistribution.png"), width = 45, height = 45, units = "mm")

####################
# Plot: Tweet types
####################
gg_tweettypes <- rt_edges %>% 
  mutate(RT_type = ifelse(RT_type == "Presumed Phantom RT", "Phantom RT", RT_type)) %>% 
  mutate(RT_type = gsub(" RT", "", RT_type)) %>% 
  ggplot(., aes(x = RT_type, fill = RT_type)) +
  geom_histogram(stat = "count") +
  scale_fill_manual(values = c("#1B264F", "#274690", "#576CA8", "#302B27", "Green")) +
  xlab("Retweet type") +
  ylab("Count") +
  theme_ctokita() +
  theme(legend.position = "none")
gg_tweettypes
ggsave(gg_tweettypes, filename = paste0(outpath, "RT_type.png"), width = 60, heigh = 60, units = "mm", dpi = 400)


############################## Ideology in networks ##############################

####################
# Plot: Ideological diversity by article veracity
####################
gg_ideodiversity <- network_metrics %>% 
  filter(article_fc_rating %in% c("T", "FM")) %>% 
           ggplot(., aes(x = article_fc_rating, y = ideology_sd)) +
  geom_point(size = 1, stroke = 0, alpha = 0.5,
             position = position_jitter(width = 0.03)) +
  scale_x_discrete(labels = c("Fake news", "Real news")) +
  ylab("Ideol. diversity of article tweeters") +
  xlab("") +
  theme_ctokita()
gg_ideodiversity
ggsave(gg_ideodiversity, filename = paste0(outpath, "ideodiversity_byveracity.png"), width = 45, height = 45, units = "mm")



############################## Network structure ##############################

####################
# Plot: Network density by article veracity
####################
gg_veracitydensity <- network_metrics %>% 
  filter(article_fc_rating %in% c("FM", "T")) %>% 
  # filter(total_tweets > 20) %>% 
  ggplot(., aes(x = article_fc_rating, y = network_density)) +
  geom_point(size = 1, stroke = 0, alpha = 0.5,
             position = position_jitter(width = 0.03)) +
  scale_x_discrete(labels = c("Fake news", "Real news")) +
  ylab("Network density") +
  xlab("") +
  theme_ctokita()
gg_veracitydensity
ggsave(gg_veracitydensity, filename = paste0(outpath, "networkdensity_byveracity.png"), width = 45, height = 45, units = "mm", dpi = 400)
