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
library(brms)
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

plot_color <- "#495867"
grouping_pal <- c("#F18805", plot_color)


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

# Fix rounding error for new belief
# Some cases have one more new_believing_user than new_exposed_user due to rounding
belief_data$new_believing_users[belief_data$new_believing_users > belief_data$new_exposed_users] = belief_data$new_exposed_users[belief_data$new_believing_users > belief_data$new_exposed_users]

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
         new_believing_users = 0,
         cumulative_believing = 0,
         total_article_number = rep(unique(article_data$total_article_number), each = 2)) %>% 
  merge(unique_article_ratings, by = "total_article_number") %>% 
  rbind(belief_timeseries, .) %>% 
  mutate(hour_bin = cut(time, breaks = seq(-2, 24*14, 1), include.lowest = TRUE, right = FALSE, labels = seq(-2, 24*14-1))) %>%  #bin by hour tweet appeared
  mutate(hour_bin = as.numeric(as.character(hour_bin))) %>%  #convert from factor to plain number
  group_by(total_article_number) %>% 
  mutate(relative_cumulative_exposed = cumulative_exposed / max(cumulative_exposed),
         relative_cumulative_belief = cumulative_believing / max(cumulative_believing),
         relative_tweet_count = tweet_number / max(tweet_number)) %>% 
  arrange(total_article_number, tweet_number) %>% 
  ungroup()

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

# Prep data (data is melted to make one ideological bin per tweet per row)
belief_ideol <- belief_timeseries %>% 
  filter(hour_bin >= 0) %>% 
  select(-source_lean, -relative_cumulative_exposed, -relative_cumulative_belief, -relative_tweet_count, -follower_count) %>% 
  gather(key = "ideology_bin", value = "count", 
         -time, -tweet_number, -tweet_id, -user_id, -user_ideology, -new_exposed_users, -cumulative_exposed, -new_believing_users, -cumulative_believing, -total_article_number, -hour_bin, -source_type, -article_fc_rating, -article_lean) %>% 
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
             strip.position = "top",
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
             strip.position = "top",
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
                     limits = c(0, 0.2),
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
             strip.position = "top",
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
  geom_bar(position = "fill", stat = "identity", width = 1, size = 0.01) +
  scale_fill_gradientn(colours = ideol_pal, 
                       name = "User\nideology",
                       limits = c(-2, 2), 
                       oob = squish) +
  scale_color_gradientn(colours = ideol_pal, 
                       name = "User\nideology",
                       limits = c(-2, 2), 
                       oob = squish) +
  scale_x_continuous(breaks = seq(0, 72, 12),
                     expand = c(0, 0),
                     limits = c(-0.5, 72.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Time since first article share (hrs)") +
  ylab("Proportion of newly believing users") +
  theme_ctokita() +
  theme(aspect.ratio = NULL, 
        legend.box.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_rect(size = 0.6, fill = NA),
        axis.line = element_blank()) +
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "top",
             scales = "free")
gg_ideoltime
ggsave(gg_ideoltime, filename = paste0(outpath, "ideol_belief_hourbin.pdf"), width = 90, height = 90, units = "mm", dpi = 400)



####################
# Belief as function of exposure over time
####################
# Prep data
belief_per_exposure <- belief_timeseries %>% 
  mutate(belief_per_exposure = new_believing_users / new_exposed_users) %>% 
  filter(time >= 0 & time <= 48)

# Fit bayesian trend line of belief-per-exposure over time
if (grouping == "article_fc_rating") {
  belief_split <- belief_per_exposure %>% 
    filter(!is.na(belief_per_exposure)) %>% 
    split(.$article_fc_rating)
  group_names <- unique(belief_timeseries$article_fc_rating)
} else if (grouping == "source_type") {
  belief_split <- belief_per_exposure %>% 
    filter(!is.na(belief_per_exposure)) %>% 
    split(.$source_type)
  group_names <- unique(belief_timeseries$source_type)
}

regression_belief <- brm_multiple(data = belief_split,
                                 formula = belief_per_exposure ~ 1 + time + I(time^2),
                                 # formula = belief_per_exposure ~ 1 + time,
                                 family = gaussian(), #assume normally distributed error
                                 prior = c(prior(uniform(-100, 100), class = Intercept),
                                           prior(uniform(-100, 100), class = b),
                                           prior(uniform(-100, 100), class = sigma)),
                                 iter = 3500,
                                 warmup = 1000,
                                 chains = 4,
                                 seed = 323,
                                 combine = FALSE)

# Get stats on slope
slope_fm <- posterior_samples(regression_belief[[1]], pars = "time", as.array = TRUE)
slope_t <- posterior_samples(regression_belief[[2]], pars = "time", as.array = TRUE)
hypothesis_t_less_fm <- sum(slope_fm > slope_t) / length(slope_fm) #what set of posterior samples. P = 1
bayes_factor <- hypothesis_t_less_fm / (1 - hypothesis_t_less_fm) #Inf

# Get fitted values from model to data range/space
x_values <- data.frame(time = seq(0, 48, 0.1))
fit_belief <- lapply(seq(1:length(group_names)), function(i) {
  group <- group_names[i]
  fit_line <- fitted(regression_belief[[i]], newdata = x_values) %>% 
    as.data.frame() %>% 
    mutate(group = group,
           time = x_values$time)
})
fit_belief <- do.call("rbind", fit_belief)



# Plot trend line and binned mean of belief-per-exposure over time
minutes_per_bin <- 5

gg_belief_rate_exposure <- belief_per_exposure %>% 
  # Prep data
  mutate(binned_time = floor(time / (minutes_per_bin / 60)  ), #assign to bin
         binned_time = binned_time / (60 / minutes_per_bin)) %>% #translate bin into real time
  group_by(!!sym(grouping), binned_time) %>% 
  summarise(belief_per_exposure = mean(belief_per_exposure, na.rm = TRUE)) %>% 
  rename(group = !!sym(grouping)) %>% 
  # Plot
  ggplot(., aes(x = binned_time, y = belief_per_exposure, color = group)) +
  geom_point(stroke = 0, alpha = 0.15, size = 1) +
  geom_line(data = fit_belief, aes(x = time, y = Estimate),
            size = 0.6) +
  scale_x_continuous(breaks = seq(0, 48, 12),
                     limits = c(0, 48)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1),
                     limits = c(0.2, 0.9),
                     expand = c(0, 0)) +
  scale_color_manual(values = grouping_pal, 
                     name = "Article rating", 
                     labels = c("False/Misleading", "True")) +
  xlab("Time since first article share (hrs.)") +
  ylab("Beliefs per exposure") +
  # facet_wrap(~group,
  #            ncol = 1,
  #            strip.position = "top",
  #            scales = "free_y") +
  theme_ctokita() +
  theme(legend.position = "none")
gg_belief_rate_exposure 
ggsave(gg_belief_rate_exposure, filename = paste0(outpath, "beliefs_per_exposure.pdf"), width = 45, height = 45, units = "mm", dpi = 400)


# Slopes of bayesian regression
regression_lines_summary <- data.frame(matrix(ncol = 5, nrow = 0))
names(regression_lines_summary) <- c("group", "Estimate", "Est. Error", "Q2.5", "Q97.5")
for (i in seq(1, length(regression_belief) ) ) {
  estimate_vals <- data.frame( fixef(regression_belief[[i]])) #grab fixed effect values
  which_row <- row.names(estimate_vals) == "time" #which row has our slope estimate
  estimate_vals <- estimate_vals[which_row, ] #grab slope row
  estimate_vals$group <- group_names[i] #label with grouoping
  rownames(regression_lines_summary) <- NULL #remove row names
  regression_lines_summary <- rbind(regression_lines_summary, estimate_vals) #bind together
}


####################
# Relative cumulative belief over first 24 hours
####################
gg_24hr_belief <- belief_timeseries %>% 
  # Prep data
  filter(time >= 0 & time <= 48) %>% 
  group_by(total_article_number) %>% 
  arrange(time) %>% 
  mutate(relative_cumulative_exposed = cumulative_exposed / max(cumulative_exposed),
         relative_cumulative_belief = cumulative_believing / max(cumulative_believing)) %>% 
  group_by(!!sym(grouping), hour_bin) %>% 
  summarise(mean_relative_cumulative_belief = mean(relative_cumulative_belief, na.rm = TRUE),
            sd_belief = sd(relative_cumulative_belief, na.rm = TRUE)) %>% 
  mutate(lower = pmax(0, mean_relative_cumulative_belief - sd_belief),
         upper = pmin(1, mean_relative_cumulative_belief + sd_belief)) %>% #don't allow to go above 1 or below 0
  # Plot
  ggplot(., aes(x = hour_bin, y = mean_relative_cumulative_belief, color = !!sym(grouping), fill = !!sym(grouping))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2) +
  geom_line(size = 0.6) +
  xlab("Hours since first article share") +
  ylab("Relative cumulative belief") +
  scale_x_continuous(breaks = seq(0, 24, 6),
                     limits = c(0, 24)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), 
                     labels = c("0.0", "", "0.5", "", "1.0"),
                     expand = c(0, 0), 
                     limits = c(0, 1)) +
  scale_color_manual(values = grouping_pal, name = "Article rating", labels = c("False/Misleading", "True")) +
  scale_fill_manual(values = grouping_pal, name = "Article rating", labels = c("False/Misleading", "True")) +
  # facet_wrap(as.formula(paste("~", grouping)),
  #            ncol = 1,
  #            strip.position = "top") +
  theme_ctokita() +
  theme(legend.position = c(0.75, 0.2))
  
gg_24hr_belief
ggsave(gg_24hr_belief, filename = paste0(outpath, "relative_cumulative_belief.pdf"), width = 45, height = 45, units = "mm", dpi = 400)


# Just raw article data
gg_24hr_belief <- belief_timeseries %>% 
  ggplot(., aes(x = time, y = relative_cumulative_belief, color = !!sym(grouping), group = total_article_number)) +
  geom_line(size = 0.6, alpha = 0.15) +
  xlab("Hours since first article share") +
  ylab("Relative cumulative belief over first 24 hrs.") +
  scale_x_continuous(breaks = seq(0, 24, 6),
                     limits = c(0, 24)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), 
                     labels = c("0.0", "", "0.5", "", "1.0"),
                     expand = c(0, 0), 
                     limits = c(0, 1)) +
  scale_color_manual(values = grouping_pal, name = "Article rating", labels = c("False/Misleading", "True")) +
  facet_wrap(as.formula(paste("~", grouping)),
             ncol = 1,
             strip.position = "top") +
  theme_ctokita() +
  theme(legend.position = c(0.75, 0.2))

gg_24hr_belief


####################
# Follower count of users sharing news
####################
gg_follower_count <- ggplot(belief_timeseries, aes(x = time, y = new_exposed_users, color = !!sym(grouping))) +
  geom_point(stroke = 0, size = 2, alpha = 0.5) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(limits = c(0, 7*24)) +
  scale_color_manual(values = grouping_pal, name = "Article rating", labels = c("False/Misleading", "True")) +
  theme_ctokita() +
  facet_wrap(as.formula(paste("~", grouping)),
             ncol = 1,
             strip.position = "top",
             scales = "free_y") 
gg_follower_count
