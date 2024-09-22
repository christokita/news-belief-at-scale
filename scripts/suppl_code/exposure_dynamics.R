########################################
#
# PLOT: Show how choice in exposure dynamic affects competing tweets
#
########################################

####################
# Load packages
####################
library(truncnorm)
library(ggplot2)
source("scripts/_plot_themes/theme_ctokita.R")

# Important paths
outpath <- "output/interventions/time_to_exposure/"

####################
# Parameters
####################
n_followers = 100000
mean_exposure_time <- 1
sd_exposure_time <- 1
response_delay <- 2

####################
# Paramters
####################
fake_story <- data.frame(story_type = "fake_news",
                         exposure = rtruncnorm(n = n_followers, a = 0, b = Inf, mean = mean_exposure_time, sd = sd_exposure_time)) 

response <- data.frame(story_type = "correction_response",
                       exposure = rtruncnorm(n = n_followers, a = 0, b = Inf, mean = mean_exposure_time, sd = sd_exposure_time) + response_delay)

tweet_exposure <- rbind(fake_story, response) %>% 
  mutate(story_type = factor(story_type, levels = c("fake_news", "correction_response")))


####################
# Plot just fake news exposure distribution
####################
gg_exposure_dist <- ggplot(tweet_exposure %>% filter(story_type == "fake_news"), aes(x = exposure, fill = story_type)) +
  geom_histogram(breaks = seq(0, 16, 0.25), 
                 position = "identity",
                 alpha = 0.5) +
  scale_x_continuous(breaks = seq(0,24, 2)) +
  scale_fill_manual(values = c("#F64740", "#8F6593"),
                    name = "") +
  xlab("Time (hrs.)") +
  ylab("N. users exposed") +
  theme_ctokita()
gg_exposure_dist

ggsave(gg_exposure_dist, filename = paste0(outpath, "exposuredsit_mean", mean_exposure_time, "_sd", sd_exposure_time, ".png"),
       width = 90, height = 45, units = "mm", dpi = 400)

####################
# Plot competition between fake news and correction
####################
gg_competition <- ggplot(tweet_exposure, aes(x = exposure, fill = story_type)) +
  geom_histogram(breaks = seq(0, 16, 0.25), 
                 position = "identity",
                 alpha = 0.5) +
  scale_x_continuous(breaks = seq(0,24, 2)) +
  scale_fill_manual(values = c("#F64740", "#8F6593"),
                    name = "") +
  xlab("Time (hrs.)") +
  ylab("N. users exposed") +
  theme_ctokita()
gg_competition

ggsave(gg_competition, filename = paste0(outpath, "tweetcompetition_mean", mean_exposure_time, "_sd", sd_exposure_time, ".png"),
       width = 90, height = 45, units = "mm", dpi = 400)
