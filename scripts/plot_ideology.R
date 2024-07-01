#########################################
# Name: `plot_ideology.R`
# Author: Chris Tokita
# Purpose: Plot the distribution of ideologies among the Tweeters of tracked news articles.
# Details:
#   (These R scripts assume the use of the `.Rproj` at top of the news-belief-at-scale/ repo. Otherwise, set the working directory to one level above this script.)
#
#   The Variables at the beginning of the script that are in all caps need to be set by the user:
#     `DATA_DIRECTORY`: path to the data directory. (Copies of data are currently stored on external hard drive and high-performance cluster.)
#     `GROUPING`:       determines whether the plots will break out tweet belief according to article veracity ("article_fc_rating") or the source of the article ("source_type").
# 
# Data In:
# `<data storage location>/data_derived/tweets/tweets_labeled.csv`: article tweets with article and tweeter metadata.
# `<data storage location>/data_derived/ideological_scores/estimated_ideol_distributions/follower_ideology_distribution_shapes.csv`: estimated ideological distribution of each tweeter's followers.
# 
# Data Out: Plots written to output sub-folder depending on if we are comparing article veracity or news source type. 
# `output/tweeter_ideology/veracity/`
# `output/tweeter_ideology/source_type/`
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
library(ggExtra)
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
follower_ideology_path <- paste0(DATA_DIRECTORY, "data_derived/ideological_scores/estimated_ideol_distributions/follower_ideology_distribution_shapes.csv") #estimated ideological distribution of each tweeter's followers

# Set path for plots
outpath <- 'output/tweeter_ideology/'
if (GROUPING == "article_fc_rating") {
  subdir_out <- 'veracity/'
} else if(GROUPING == "source_type") {
  subdir_out <- 'source_type/'
}


# Set color palette
plot_color <- "#495867"
ideol_pal <- rev(brewer.pal(5, "RdBu"))
ideol_pal[3] <- "#e0e0e0"
source_pal <- c(ideol_pal[c(5,1)], "#9D69A3", "#C5CBD3")
ideol_limit <- 3 #limit beyond which we squish the color palette


####################
# Load data 
####################
# Read in data
tweets <- read.csv(tweet_path, header = TRUE, colClasses = c('user_id'='character')) %>% 
  mutate(article_ideology = article_con_feel - article_lib_feel) %>% 
  filter(total_article_number > 10) %>%  #discard first 10 articles from analysis
  # clean up labels
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "False/Misleading news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream outlet", ifelse(source_type == "fringe", "Fringe outlet", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) ) %>% 
  # bin ideologies
  mutate(bin = cut(user_ideology, breaks = seq(-5.5, 5.5, 0.5))) %>% 
  mutate(lower_edge = as.numeric( gsub("^[^0-9]([-\\.0-9]+),.*", "\\1", bin, perl = TRUE) ),
         upper_edge = as.numeric( gsub(".*,([-\\.0-9]+)[^0-9]$", "\\1", bin, perl = TRUE) )) %>% 
  mutate(ideology_bin = (lower_edge + upper_edge) / 2) %>% 
  select(-bin, -lower_edge, -upper_edge) %>% 
  # count number of articles per GROUPING (useful for average distributions) %>% 
  group_by(!!sym(GROUPING)) %>% 
  mutate(n_articles_in_GROUPING = length(unique(total_article_number)))

# If analyzing by veracity, drop out non-True/False articles
if (GROUPING == "article_fc_rating") {
  tweets <- tweets %>% 
    filter(article_fc_rating %in% c("True news", "False/Misleading news"))
}


# Load follower distribution shape data
follower_ideology_distributions <- read.csv(follower_ideology_path) %>% 
  mutate(user_id = gsub("\"" , "", user_id_str)) %>% 
  merge(tweets, by = 'user_id')

# Confirm we have proper user set
stopifnot( length(unique(tweets$user_id)) == length(unique(follower_ideology_distributions$user_id)) )



############################## Summary of tweeter ideology ##############################

####################
# PLOT: Total tweets by ideology of tweeter
####################
gg_ideoltweets <- tweets %>% 
  # Calculate tweets per ideology bin
  filter(!is.na(user_ideology)) %>% 
  group_by(!!sym(GROUPING), ideology_bin) %>% 
  summarise(count = n()) %>% 
  # plot
  ggplot(., aes(x = ideology_bin, y = count, fill = ideology_bin)) +
  geom_bar(stat = "identity", width = 0.5, color = NA) +
  xlab("Tweeter ideology") +
  ylab("Number of tweets") +
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 2)) +
  scale_y_continuous(expand = c(0, 0),
                     label = scales::comma) +
  scale_fill_gradientn(colors = ideol_pal, limits = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_color_gradientn(colors = ideol_pal, limits = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.position = "none") +
  facet_wrap(as.formula(paste("~", GROUPING)), 
             scales = 'free', 
             ncol = 1, 
             strip.position = "right")

gg_ideoltweets
ggsave(gg_ideoltweets, file = paste0(outpath, subdir_out, "total_tweets_by_ideology.pdf"), width = 45, height = 90, units = "mm", dpi = 400)
ggsave(gg_ideoltweets + theme(strip.text = element_blank()), file = paste0(outpath, subdir_out, "total_tweets_by_ideology_FIG.pdf"), width = 45, height = 45, units = "mm", dpi = 400) #facet labels removed and manually added back during figure creation in vector art program


####################
# PLOT: Average ideology distribution of tweeters
####################
gg_ideoldist_tweeters <- tweets %>% 
  # Calculate distribution of tweeter ideologies per article
  filter(!is.na(user_ideology)) %>% 
  group_by(n_articles_in_GROUPING, !!sym(GROUPING), total_article_number, ideology_bin) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  group_by(!!sym(GROUPING), total_article_number) %>% 
  mutate(tweeter_prop = count / sum(count)) %>% 
  # Now determine average distribution shape by article GROUPING
  group_by(!!sym(GROUPING), ideology_bin) %>% 
  summarise(avg_tweeter_prop = sum(tweeter_prop) / unique(n_articles_in_GROUPING)) %>% 
  # Plot
  ggplot(., aes(x = ideology_bin, y = avg_tweeter_prop, fill = ideology_bin)) +
  geom_bar(stat = "identity", width = 0.5, color = NA) +
  xlab("Tweeter ideology") +
  ylab("Avg. proportion of article tweeters") +
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 2)) +
  scale_y_continuous(limits = c(0, 0.25),
                     expand = c(0, 0),
                     breaks = seq(0, 0.3, 0.05),
                     label = scales::comma) +
  scale_fill_gradientn(colors = ideol_pal, limits = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_color_gradientn(colors = ideol_pal, limits = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.position = "none") +
  facet_wrap(as.formula(paste("~", GROUPING)), 
             scales = 'free_x', 
             ncol = 1, 
             strip.position = "right")

gg_ideoldist_tweeters
ggsave(gg_ideoldist_tweeters, file = paste0(outpath, subdir_out, "average_tweeter_ideoldist.pdf"), width = 45, height = 90, units = "mm", dpi = 400)


####################
# PLOT: Average ideology distribution of retweets only
####################
gg_ideoldist_retweet <- tweets %>% 
  # Calculate distribution of tweeter ideologies per article
  filter(!is.na(user_ideology),
         is_retweet == "True") %>% 
  group_by(n_articles_in_GROUPING, !!sym(GROUPING), total_article_number, ideology_bin) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  group_by(!!sym(GROUPING), total_article_number) %>% 
  mutate(tweeter_prop = count / sum(count)) %>% 
  # Now determine average distribution shape by article GROUPING
  group_by(!!sym(GROUPING), ideology_bin) %>% 
  summarise(avg_tweeter_prop = sum(tweeter_prop) / unique(n_articles_in_GROUPING)) %>% 
  # Plot
  ggplot(., aes(x = ideology_bin, y = avg_tweeter_prop, fill = ideology_bin, color = ideology_bin)) +
  geom_bar(stat = "identity", width = 0.5, color = NA) +
  xlab("Retweeter ideology") +
  ylab("Avg. proportion of article retweeters") +
  scale_x_continuous(limits = c(-6, 6), 
                     expand = c(0, 0), 
                     breaks = seq(-6, 6, 2)) +
  scale_y_continuous(limits = c(0, 0.25),
                     expand = c(0, 0),
                     breaks = seq(0, 0.3, 0.05)) +
  scale_fill_gradientn(colors = ideol_pal, limits = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  scale_color_gradientn(colors = ideol_pal, limits = c(-ideol_limit, ideol_limit), oob = scales::squish) +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.position = "none") +
  facet_wrap(as.formula(paste("~", GROUPING)), 
             scales = 'free_x', 
             ncol = 1, 
             strip.position = "right")

gg_ideoldist_retweet
ggsave(gg_ideoldist_retweet, file = paste0(outpath, subdir_out, "average_retweeter_ideoldist.pdf"), width = 45, height = 90, units = "mm", dpi = 400)



############################## Plot follower distribution follower_ideology_distributions ##############################

####################
# Plot shape of follower distributions
####################
gg_dist_follower_ideology_distributions <- 
  # Calculate mean SD for each skew normal distribution (see: https://en.wikipedia.org/wiki/Skew_normal_distribution)
  follower_ideology_distributions %>% 
  filter(basis == "followers") %>% 
  select(user_id, mu, sigma, alpha, user_ideology) %>% 
  mutate(delta = alpha / sqrt(1 + alpha**2)) %>% 
  mutate(mean = mu + sigma * delta * sqrt(2/pi)) %>% 
  mutate(variance = sigma**2 * (1 - ((2*delta**2) / pi)) ) %>% 
  mutate(sd = sqrt(variance)) %>% 
  distinct() %>% 
  # Plot
  ggplot(., aes(x = mean, y = sd, color = user_ideology)) +
  geom_point(alpha = 0.1, size = 1, stroke = 0) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-2, 2), oob = scales::squish, name = "Tweeter\nideology") +
  xlab("Estimated mean ideology\nof tweeter's followers") +
  ylab("Estimated s.d. ideology\nof tweeter's followers") +
  theme_ctokita() +
  theme(plot.background = element_blank(),
        legend.position = "right")

gg_dist_follower_ideology_distributions
ggsave(gg_dist_follower_ideology_distributions, filename = paste0(outpath, "follower_dist_follower_ideology_distributions.png"), width = 65, height = 50, units = "mm", dpi = 400, bg = "transparent")


####################
# Plot follower distribution mu vs. user_ideology
####################
gg_ideology_comp <- follower_ideology_distributions %>% 
  select(user_id, mu, sigma, user_ideology, !!(GROUPING)) %>% 
  distinct() %>% 
  filter(!is.na(user_ideology)) %>% 
  ggplot(., aes(x = user_ideology, y = mu, color = user_ideology)) +
  geom_point(alpha = 0.1, size = 1, stroke = 0) +
  scale_y_continuous(breaks = seq(-6, 6, 1)) +
  scale_x_continuous(breaks = seq(-6, 6, 1)) +
  scale_color_gradientn(name = "User\nideology",
                        colours = ideol_pal, 
                        limit = c(-2, 2), 
                        oob = scales::squish) +
  xlab("User ideology") +
  ylab("Est. mean follower ideology") +
  theme_ctokita() +
  theme(legend.position = "none", 
        plot.background = element_blank())

gg_ideology_comp
ggsave(gg_ideology_comp, filename = paste0(outpath, "ideology_comparison_user-follower.png"), width = 45, height = 45, units = "mm", dpi = 400, bg = "transparent")



