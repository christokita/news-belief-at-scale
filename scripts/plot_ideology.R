########################################
#
# PLOT: Distribution of ideologies of tweeters
#
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
# Paramters for analysis: paths to data, paths for output, and filename
####################
# Choose grouping of interest. Options: 
#     (1) article veracity: "article_fc_rating"
#     (2) source: "source_type"
grouping <- "article_fc_rating"

# Paths to data
tweeter_score_path <- '/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/tweets_labeled.csv' 
if (grouping == "article_fc_rating") {
  outpath <- 'output/tweeter_ideology/veracity/'
} else if(grouping == "source_type") {
  outpath <- 'output/tweeter_ideology/source_type/'
}

# Color palette
plot_color <- "#495867"
ideol_pal <- rev(brewer.pal(5, "RdBu"))
ideol_pal[3] <- "#e0e0e0"
source_pal <- c(ideol_pal[c(5,1)], "#9D69A3", "#C5CBD3")
ideol_limit <- 3 #limit beyond which we squish the color palette


####################
# Load data 
####################
# Read in data
tweets <- read.csv(tweeter_score_path, header = TRUE, colClasses = c('user_id'='character')) %>% 
  mutate(article_ideology = article_con_feel - article_lib_feel) %>% 
  filter(total_article_number > 10) %>%  #discard first 10 articles from analysis
  # clean up labels
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "False/Misleading news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream", ifelse(source_type == "fringe", "Fringe", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) ) %>% 
  # bin ideologies
  mutate(bin = cut(user_ideology, breaks = seq(-5.5, 5.5, 0.5))) %>% 
  mutate(lower_edge = as.numeric( gsub("^[^0-9]([-\\.0-9]+),.*", "\\1", bin, perl = TRUE) ),
         upper_edge = as.numeric( gsub(".*,([-\\.0-9]+)[^0-9]$", "\\1", bin, perl = TRUE) )) %>% 
  mutate(ideology_bin = (lower_edge + upper_edge) / 2) %>% 
  select(-bin, -lower_edge, -upper_edge) %>% 
  # count number of articles per grouping (useful for average distributions) %>% 
  group_by(!!sym(grouping)) %>% 
  mutate(n_articles_in_grouping = length(unique(total_article_number)))

# If analyzing by veracity, drop out non-True/False articles
if (grouping == "article_fc_rating") {
  tweets <- tweets %>% 
    filter(article_fc_rating %in% c("True news", "False/Misleading news"))
}


# Load follower distribution shape data
shapes <- read.csv("/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/estimated_ideol_distributions/follower_ideology_distribution_shapes.csv") %>% 
  mutate(user_id = gsub("\"" , "", user_id_str)) %>% 
  merge(tweets, by = 'user_id')

# Confirm we have proper user set
length(unique(tweets$user_id)) == length(unique(shapes$user_id))



############################## Summary of tweeter ideology ##############################

####################
# Plot: Total tweets by ideology
####################
gg_ideoltweets <- tweets %>% 
  # Calculate tweets per ideology bin
  filter(!is.na(user_ideology)) %>% 
  group_by(!!sym(grouping), ideology_bin) %>% 
  summarise(count = n()) %>% 
  # plot
  ggplot(., aes(x = ideology_bin, y = count, fill = ideology_bin, color = ideology_bin)) +
  geom_bar(stat = "identity") +
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
  facet_wrap(as.formula(paste("~", grouping)), 
             scales = 'free', 
             ncol = 1, 
             strip.position = "right")
gg_ideoltweets
ggsave(gg_ideoltweets, file = paste0(outpath, "total_tweets_by_ideology.pdf"), width = 45, height = 90, units = "mm", dpi = 400)
ggsave(gg_ideoltweets + theme(strip.text = element_blank()), file = paste0(outpath, "total_tweets_by_ideology_FIG.pdf"), width = 45, height = 45, units = "mm", dpi = 400)


####################
# Plot: Average ideology distribution of tweeters
####################
gg_ideoldist_tweeters <- tweets %>% 
  # Calculate distribution of tweeter ideologies per article
  filter(!is.na(user_ideology)) %>% 
  group_by(n_articles_in_grouping, !!sym(grouping), total_article_number, ideology_bin) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  group_by(!!sym(grouping), total_article_number) %>% 
  mutate(tweeter_prop = count / sum(count)) %>% 
  # Now determine average distribution shape by article grouping
  group_by(!!sym(grouping), ideology_bin) %>% 
  summarise(avg_tweeter_prop = sum(tweeter_prop) / unique(n_articles_in_grouping)) %>% 
  # Plot
  ggplot(., aes(x = ideology_bin, y = avg_tweeter_prop, fill = ideology_bin, color = ideology_bin)) +
  geom_bar(stat = "identity") +
  xlab("Tweeter ideology") +
  ylab("Avg. proportion of article tweeters") +
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
  facet_wrap(as.formula(paste("~", grouping)), 
             scales = 'free_x', 
             ncol = 1, 
             strip.position = "right")
gg_ideoldist_tweeters
ggsave(gg_ideoldist_tweeters, file = paste0(outpath, "average_tweeter_ideoldist.pdf"), width = 45, height = 90, units = "mm", dpi = 400)


####################
# Plot: Average ideology distribution of retweets only
####################
gg_ideoldist_retweet <- tweets %>% 
  # Calculate distribution of tweeter ideologies per article
  filter(!is.na(user_ideology),
         is_retweet == "True") %>% 
  group_by(n_articles_in_grouping, !!sym(grouping), total_article_number, ideology_bin) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  group_by(!!sym(grouping), total_article_number) %>% 
  mutate(tweeter_prop = count / sum(count)) %>% 
  # Now determine average distribution shape by article grouping
  group_by(!!sym(grouping), ideology_bin) %>% 
  summarise(avg_tweeter_prop = sum(tweeter_prop) / unique(n_articles_in_grouping)) %>% 
  # Plot
  ggplot(., aes(x = ideology_bin, y = avg_tweeter_prop, fill = ideology_bin, color = ideology_bin)) +
  geom_bar(stat = "identity") +
  xlab("Tweeter ideology") +
  ylab("Avg. proportion of article re-tweeters") +
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
  facet_wrap(as.formula(paste("~", grouping)), 
             scales = 'free_x', 
             ncol = 1, 
             strip.position = "right")
gg_ideoldist_retweet
ggsave(gg_ideoldist_retweet, file = paste0(outpath, "average_retweeter_ideoldist.pdf"), width = 45, height = 90, units = "mm", dpi = 400)



############################## Plot follower distribution shapes ##############################

####################
# Plot shape of follower distributions
####################
gg_dist_shapes <- shapes %>% 
  filter(basis == "followers") %>% 
  select(user_id, mu, sigma, user_ideology) %>% 
  distinct() %>% 
  ggplot(., aes(x = mu, y = sigma, color = user_ideology)) +
  geom_point(alpha = 0.1, size = 1, stroke = 0) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-2, 2), oob = scales::squish, name = "User\nideology") +
  xlab("Estimated mean ideology of user's followers") +
  ylab("Estimated s.d. ideology of user's followers") +
  theme_ctokita() +
  theme(plot.background = element_blank(),
        legend.position = "right")
gg_dist_shapes

ggsave(gg_dist_shapes, filename = paste0(outpath, "follower_dist_shapes.png"), width = 90, height = 90, units = "mm", dpi = 400, bg = "transparent")


####################
# Plot follower distribution mu vs. user_ideology
####################
gg_ideology_comp <- shapes %>% 
  select(user_id, mu, sigma, user_ideology, !!(grouping)) %>% 
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

# gg_ideology_comp <- ggMarginal(gg_ideology_comp, type = "histogram", fill = plot_color, alpha = 0.6, size = 3, colour = NA,
#                                xparams = list( bins = 30 ),
#                                yparams = list( bins = 30 ))
gg_ideology_comp

ggsave(gg_ideology_comp, filename = paste0(outpath, "ideology_comparison_user-follower.png"), width = 45, height = 45, units = "mm", dpi = 400, bg = "transparent")



