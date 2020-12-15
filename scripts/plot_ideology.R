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
  outpath <- 'output/ideology/veracity/'
} else if(grouping == "source_type") {
  outpath <- 'output/ideology/source_type/'
}

# Color palette
plot_color <- "#495867"
ideol_pal <- rev(brewer.pal(5, "RdBu"))
ideol_pal[3] <- "#e0e0e0"
source_pal <- c(ideol_pal[c(5,1)], "#9D69A3", "#C5CBD3")


####################
# Load data 
####################
# Read in data
tweets <- read.csv(tweeter_score_path, header = TRUE, colClasses = c('user_id'='character')) %>% 
  mutate(article_ideology = article_con_feel - article_lib_feel) %>% 
  filter(total_article_number > 10) %>%  #discard first 10 articles from analysis
  # clean up labels
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "Fake news", 
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
    filter(article_fc_rating %in% c("True news", "Fake news"))
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
  ggplot(., aes(x = ideology_bin, y = count, fill = ideology_bin)) +
  geom_bar(stat = "identity") +
  xlab("Tweeter ideology") +
  ylab("Number of tweets") +
  scale_x_continuous(limits = c(-5.5, 5.5), expand = c(0, 0), breaks = seq(-5, 5, 1)) +
  scale_fill_gradientn(colors = ideol_pal, limits = c(-2, 2), oob = scales::squish) +
  scale_color_gradientn(colors = ideol_pal, limits = c(-2, 2), oob = scales::squish) +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.position = "none") +
  facet_wrap(as.formula(paste("~", grouping)), 
             scales = 'free', 
             ncol = 1, 
             strip.position = "right")
gg_ideoltweets
ggsave(gg_ideoltweets, file = paste0(outpath, "total_tweets_by_ideology.pdf"), width = 90, height = 90, units = "mm", dpi = 400)


# Plot ideological extremity
gg_veracityextr <- tweets %>% 
  filter(!is.na(user_ideology)) %>% 
  ggplot(., aes(x = user_ideol_extremity, fill = ..x..)) +
  geom_histogram(position = 'identity', binwidth = 0.25) +
  theme_ctokita() +
  facet_wrap(as.formula(paste("~", grouping)), 
             scales = 'free', 
             ncol = 1, 
             strip.position = "right")
gg_veracityextr


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
  ggplot(., aes(x = ideology_bin, y = avg_tweeter_prop, fill = ideology_bin)) +
  geom_bar(stat = "identity") +
  xlab("Tweeter ideology") +
  ylab("Avg. proportion of article tweeters") +
  scale_x_continuous(limits = c(-5.5, 5.5), expand = c(0, 0), breaks = seq(-5, 5, 1)) +
  scale_fill_gradientn(colors = ideol_pal, limits = c(-2, 2), oob = scales::squish) +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.position = "none") +
  facet_wrap(as.formula(paste("~", grouping)), 
             scales = 'free', 
             ncol = 1, 
             strip.position = "right")
gg_ideoldist_tweeters
ggsave(gg_ideoldist_tweeters, file = paste0(outpath, "average_tweeter_ideoldist.pdf"), width = 90, height = 90, units = "mm", dpi = 400)


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
  ggplot(., aes(x = ideology_bin, y = avg_tweeter_prop, fill = ideology_bin)) +
  geom_bar(stat = "identity") +
  xlab("Tweeter ideology") +
  ylab("Avg. proportion of article re-tweeters") +
  scale_x_continuous(limits = c(-5.5, 5.5), expand = c(0, 0), breaks = seq(-5, 5, 1)) +
  scale_fill_gradientn(colors = ideol_pal, limits = c(-2, 2), oob = scales::squish) +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.position = "none") +
  facet_wrap(as.formula(paste("~", grouping)), 
             scales = 'free', 
             ncol = 1, 
             strip.position = "right")
gg_ideoldist_retweet
ggsave(gg_ideoldist_retweet, file = paste0(outpath, "average_retweeter_ideoldist.pdf"), width = 90, height = 90, units = "mm", dpi = 400)

# Plot by total article number
gg_ideoldist_article <- tweets %>% 
  # Calculate distribution of tweeter ideologies per article
  filter(!is.na(user_ideology),
         # is_retweet == "True",
         article_fc_rating == "True news") %>% 
  group_by(!!sym(grouping), total_article_number, ideology_bin) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  group_by(!!sym(grouping), total_article_number) %>% 
  mutate(tweeter_prop = count / sum(count)) %>% 
  # Plot
  ggplot(., aes(x = ideology_bin, y = tweeter_prop, fill = ideology_bin)) +
  geom_bar(stat = "identity") +
  xlab("Tweeter ideology") +
  ylab("Avg. proportion of article re-tweeters") +
  scale_x_continuous(limits = c(-5.5, 5,5), expand = c(0, 0), breaks = seq(-5, 5, 1)) +
  scale_fill_gradientn(colors = ideol_pal, limits = c(-2, 2), oob = scales::squish) +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.position = "none") +
  facet_wrap(~total_article_number, 
             scales = 'free')
gg_ideoldist_article

############################## Plot follower distribution shapes ##############################

####################
# Plot shape of follower distributions
####################
gg_dist_shapes <- shapes %>% 
  filter(basis == "followers") %>% 
  select(user_id, mu, sigma, user_ideology) %>% 
  distinct() %>% 
  ggplot(., aes(x = mu, y = sigma, color = user_ideology)) +
  geom_point(alpha = 0.15, size = 0.3, stroke = 0) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-2, 2), oob = scales::squish, name = "User\nideology") +
  xlab("Est. mean follower ideology") +
  ylab("Est. s.d. follower ideology") +
  theme_ctokita()
gg_dist_shapes

ggsave(gg_dist_shapes, filename = paste0(outpath, "follower_dist_shapes.pdf"), width = 60, height = 45, units = "mm", dpi = 400)


####################
# Plot shape of follower distributions by grouping
####################
gg_dist_grouping <- shapes %>% 
  filter(basis == "followers") %>% 
  select(user_id, mu, !!(grouping), sigma, user_ideology) %>% 
  distinct() %>% 
  ggplot(., aes(x = mu, y = sigma, color = user_ideology)) +
  geom_point(alpha = 0.15, size = 0.3, stroke = 0) +
  scale_color_gradientn(colours = ideol_pal, limit = c(-2, 2), oob = scales::squish, name = "User\nideology") +
  xlab("Est. mean follower ideology") +
  ylab("Est. s.d. follower ideology") +
  theme_ctokita() +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "free")
gg_dist_grouping



####################
# Plot follower distribution mu vs. user_ideology
####################
gg_ideology_comp <- shapes %>% 
  select(user_id, mu, sigma, user_ideology, !!(grouping)) %>% 
  distinct() %>% 
  filter(!is.na(user_ideology)) %>% 
  ggplot(., aes(x = user_ideology, y = mu, color = user_ideology)) +
  geom_point(alpha = 0.2, size = 0.3, stroke = 0) +
  scale_y_continuous(breaks = seq(-3, 5, 1)) +
  scale_x_continuous(breaks = seq(-3, 5, 1)) +
  scale_color_gradientn(name = "User\nideology",
                        colours = ideol_pal, 
                        limit = c(-2, 2), 
                        oob = scales::squish) +
  xlab("User ideology") +
  ylab("Est. mean follower ideology") +
  theme_ctokita() +
  theme(legend.position = "none")

gg_ideology_comp <- ggMarginal(gg_ideology_comp, type = "histogram", fill = plot_color, alpha = 0.6, size = 3, colour = NA,
                               xparams = list( bins = 30 ),
                               yparams = list( bins = 30 ))
gg_ideology_comp

ggsave(gg_ideology_comp, filename = paste0(outpath, "ideology_comparison_user-follower.pdf"), width = 60, height = 60, units = "mm", dpi = 400)




############################## Ideological distribution of tweeters ##############################

####################
# Plot: Tweeter ideology distribution by article source lean
####################
# All source types
gg_fmtweeters <- tweets %>% 
  filter(article_fc_rating == "FM",
         !is.na(user_ideology)) %>% 
  ggplot(., aes(x = user_ideology, group = source_lean, fill = source_lean)) +
  geom_histogram(alpha = 0.6, position = 'identity', binwidth = 0.25) +
  xlab("Tweeter ideology") +
  ylab("Log number of tweets") +
  scale_fill_manual(values = source_pal[c(1,2,4)], name = "Article type", labels = c("False, Conservative source", "False, Liberal source", "False, Unclear source")) +
  scale_x_continuous(limits = c(-5.5, 5.5), expand = c(0, 0)) +
  scale_y_continuous(trans = "log10", limits = c(1, 12000), expand = c(0,0)) +
  theme_ctokita() +
  theme(aspect.ratio = 0.5)
gg_fmtweeters
ggsave(gg_fmtweeters, file = paste0(outpath, "FMtweeters_ideologies_bysource.pdf"), width = 120, height = 45, units = "mm", dpi = 400)

# Just liberal and conservative sources
gg_fmtweeters_libcon <- tweets %>% 
  filter(source_lean %in% c("C", "L"),
         article_fc_rating == "FM") %>% 
  ggplot(., aes(x = user_ideology, group = source_lean, fill = source_lean)) +
  geom_histogram(alpha = 0.6, position = 'identity', binwidth = 0.25) +
  xlab("Tweeter ideology") +
  ylab("Number of tweets") +
  scale_fill_manual(values = source_pal[c(1,2,4)], name = "Article type", labels = c("False, Conservative source", "False, Liberal source")) +
  scale_x_continuous(limits = c(-5.5, 5.5)) +
  theme_ctokita() +
  theme(aspect.ratio = 0.5)
gg_fmtweeters_libcon
ggsave(gg_fmtweeters_libcon, file = paste0(outpath, "FMtweeters_ideologies_lib-con-sources.pdf"), width = 120, height = 45, units = "mm", dpi = 400)


####################
# Plot: Tweeter ideology by article lean
####################

###### All articles #####
# Distribution of user ideology by article lean
gg_articlelean <- ggplot(data = tweets, aes(x = user_ideology, group = article_lean, fill = article_lean)) +
  geom_histogram(position = 'identity', binwidth = 0.25) +
  xlab("Tweeter ideology") +
  ylab("Number of tweets") +
  scale_x_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  scale_fill_manual(values = source_pal, name = "Article lean", labels = c("Conservative", "Liberal", "Neutral", "Unclear")) +
  facet_wrap(~article_lean, scales = 'free', ncol = 1) +
  theme_ctokita() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.2)
gg_articlelean
ggsave(gg_articlelean, file = paste0(outpath, "ideologies_byarticlelean.pdf"), width = 90, height = 70, units = "mm", dpi = 400)

# Tweeter ideology vs article ideology
gg_articleideol <- tweets %>% 
  filter(article_fc_rating %in% c("Fake news", "True news"),
         !is.na(user_ideology)) %>% 
  ggplot(., aes(x = user_ideology, y = article_ideology, color = article_lean)) +
  geom_hline(aes(yintercept = 0), size = 0.3, linetype = "dotted") +
  geom_vline(aes(xintercept = 0), size = 0.3, linetype = "dotted") +
  geom_point(size = 0.5,
             alpha = 0.4, 
             stroke = 0, 
             position = position_jitter(0.1)) +
  xlab("Tweeter ideology") +
  ylab("Article ideology") +
  scale_color_manual(values = source_pal, name = "Article lean", labels = c("Conservative", "Liberal", "Neutral", "Unclear")) +
  scale_x_continuous(breaks = seq(-6, 6, 2), limits = c(-4.7, 4.7)) +
  scale_y_continuous(breaks = seq(-6, 6, 2), limits = c(-4.6, 4.6)) +
  theme_ctokita() +
  facet_wrap(~article_fc_rating, scales = "free")
gg_articleideol
ggsave(gg_articleideol, file = paste0(outpath, "ideology_vs_articleideology.pdf"), width = 100, height = 45, units = "mm", dpi = 400)


###### Just FM articles #####
gg_articlelean <- tweets %>% 
  filter(article_fc_rating == "Fake news") %>% 
  ggplot(., aes(x = user_ideology, group = article_lean, fill = article_lean)) +
  geom_histogram(position = 'identity', binwidth = 0.25) +
  xlab("Tweeter ideology") +
  ylab("Number of tweets") +
  scale_x_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  scale_fill_manual(values = source_pal, name = "Article lean", labels = c("Conservative", "Liberal", "Neutral", "Unclear")) +
  facet_wrap(~article_lean, scales = 'free', ncol = 1) +
  theme_ctokita() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.2)
gg_articlelean
ggsave(gg_articlelean, file = paste0(outpath, "FMtweeters_ideologies_byarticlelean.pdf"), width = 90, height = 70, units = "mm", dpi = 400)



