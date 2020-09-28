########################################
#
# PLOT: Diversity of exposure to news
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
outpath <- 'output/ideology/'

# Color palette
line_color <- "#495867"
ideol_pal <- rev(brewer.pal(5, "RdBu"))
ideol_pal[3] <- "#e0e0e0"
ideol_dist_pal <- rev(brewer.pal(5, "PuOr"))
ideol_dist_pal[3] <- "#e0e0e0"


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

# Load exposure data 
#
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
exposure_timeseries <- merge(exposure_data, article_data, by = c("tweet_id", "total_article_number"), all = TRUE) 

# Add dummy rows of pre-first share for plotting purposes
# (1) Create empty dataframe for dummy rows
n_articles <- length(unique(exposure_timeseries$total_article_number)) #number of unique articles
dummy_rows <- data.frame(matrix(NA, ncol = ncol(exposure_timeseries), nrow = 2*n_articles))  #create empty dataframe
names(dummy_rows) <- names(exposure_timeseries) #give same column names
# (2) Create unique set of article IDs and fact-check rating to add to our dummy rows
unique_article_ratings <- article_data %>% 
  select(source_type, source_lean, total_article_number, article_fc_rating) %>% 
  unique()
# (3) Join together
exposure_timeseries <- dummy_rows %>% 
  select(-article_fc_rating, -source_type, -source_lean) %>% 
  mutate(time = rep(c(-2, -0.01), n_articles),
         tweet_number = rep(c(-1, 0), n_articles),
         new_exposed_users = 0, 
         cumulative_exposed = 0, 
         total_article_number = rep(unique(article_data$total_article_number), each = 2)) %>% 
  merge(unique_article_ratings, by = "total_article_number") %>% 
  rbind(exposure_timeseries, .) %>% 
  mutate(hour_bin = cut(time, breaks = seq(-2, 50, 1), include.lowest = TRUE, right = FALSE, labels = seq(-2, 49))) %>%  #bin by hour tweet appeared
  mutate(hour_bin = as.numeric(as.character(hour_bin))) %>%  #convert from factor to plain number
  group_by(total_article_number) %>% 
  mutate(relative_cumulative_exposed = cumulative_exposed / max(cumulative_exposed),
         relative_tweet_count = tweet_number / max(tweet_number)) %>% 
  arrange(total_article_number, tweet_number)

rm(dummy_rows, article_data, exposure_data)

# If analyzing by veracity, drop out non-True/False articles
if (grouping == "article_fc_rating") {
  exposure_timeseries <- exposure_timeseries %>% 
    filter(article_fc_rating %in% c("T", "FM"))
}

# Clean up some labels
exposure_timeseries <- exposure_timeseries %>% 
  mutate(article_fc_rating = ifelse(article_fc_rating == "T", "True news", ifelse(article_fc_rating == "FM", "Fake news", 
                                                                                  ifelse(article_fc_rating == "CND", "Borderline", 
                                                                                         ifelse(article_fc_rating == "No Mode!", "No mode", article_fc_rating)))),
         source_type = ifelse(source_type == "mainstream", "Mainstream", ifelse(source_type == "fringe", "Fringe", source_type)),
         article_lean = ifelse(article_lean == "C", "Conservative", ifelse(article_lean == "L", "Liberal",
                                                                           ifelse(article_lean == "N", "Neutral", 
                                                                                  ifelse(article_lean == "U", "Unclear", source_type)))) )



############################## Plot article ideological diversity of exposure ##############################
# Prep data
exposure_ideol <- exposure_timeseries %>% 
  select(-source_lean, -relative_cumulative_exposed, -relative_tweet_count) %>% 
  gather(key = "ideology_bin", value = "count", 
         -time, -tweet_number, -tweet_id, -user_id, -user_ideology, -follower_count, -new_exposed_users, -cumulative_exposed, -total_article_number, -hour_bin, -source_type, -article_fc_rating, -article_lean) %>% 
  mutate(ideology_bin = gsub("ideol_", "", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("^\\.", "-", ideology_bin)) %>% 
  mutate(ideology_bin = gsub("_\\.", "_-", ideology_bin)) %>% 
  separate(ideology_bin, c("lower", "upper"), sep = "_", convert = TRUE) %>% 
  mutate(ideology_bin = (lower + upper) / 2) %>% 
  select(-lower, -upper)


# Function to calcualte diversity as expected ideological difference between two randomly chosen exposed users
ideological_diversity <- function(distance_matrix, ideological_counts) {
  
  diversity_data <- c()
  for (i in 1:nrow(ideological_counts)) {
    # Normalize counts into proportions and compute ideology X ideology encounter probability
    count_row <- ideological_counts[i , ]
    normalized_counts <- as.matrix( count_row / sum(count_row) )
    encounter_probs <- t(normalized_counts) %*% normalized_counts
    
    # Compute ideological distance X probability of encounter
    weighted_difference <- encounter_probs *  distance_matrix 
    expected_difference <- sum(weighted_difference)
    diversity_data <- c(diversity_data, expected_difference)
  }
  return(diversity_data)
}

# Distance matrix between ideology bins
bin_centers <- unique(exposure_ideol$ideology_bin)
ideol_distance_matrix <- as.matrix( dist(bin_centers, diag = TRUE, upper = TRUE, method = "euclidean") )

####################
# Diversity of exposure at tweet level
####################
# Calculate 
diversity_tweet <- exposure_ideol %>% 
  filter(hour_bin >= 0) %>% 
  select(total_article_number, !!sym(grouping), tweet_id, user_id, tweet_number, time, hour_bin, new_exposed_users, ideology_bin, count) %>% 
  pivot_wider(names_from = ideology_bin, values_from = count) 
ideol_bins_cols <- grep("[-.0-9]+", names(diversity_tweet)) #list of ideology bin columns
diversity_index <- data.frame(ideol_diversity = ideological_diversity(distance_matrix = ideol_distance_matrix, ideological_counts = diversity_tweet[ , ideol_bins_cols]))
diversity_tweet <- diversity_tweet %>% 
  select(!!sym(grouping), total_article_number, user_id, tweet_id, time, hour_bin, new_exposed_users) %>% 
  cbind(diversity_index)

# Raw plot of tweets
gg_expos_diversity_raw <- diversity_tweet %>% 
  filter(hour_bin < 31, 
         !is.na(ideol_diversity)) %>% 
  ggplot(., aes(x = time, y = ideol_diversity, group = total_article_number)) +
  geom_point(size = 0.1, alpha = 0.2, color = line_color) +
  theme_ctokita() + 
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "fixed")
gg_expos_diversity_raw

# Fit bayesian regression to data
library(brms)
diversity_split <- diversity_tweet %>%
  ungroup() %>%
  select(time, ideol_diversity, !!sym(grouping)) %>% 
  filter(!is.na(ideol_diversity))
if (grouping == "article_fc_rating") {
  diversity_split <- diversity_split %>% 
    split(.$article_fc_rating)
  group_names <- unique(diversity_tweet$article_fc_rating)
} else if (grouping == "source_type") {
  diversity_split <- diversity_split %>% 
    split(.$source_type)
  group_names <- unique(diversity_tweet$source_type)
}
regression_diversity <- brm_multiple(data = diversity_split,
                                     # formula = ideol_diversity ~ 1 + time + I(time^2),
                                     formula = ideol_diversity ~ 1 + time, 
                                     prior = c(prior(uniform(-10, 10), class = Intercept),
                                               prior(normal(0, 1), class = b),
                                               prior(normal(0, 1), class = sigma)),
                                     iter = 3000,
                                     warmup = 1000,
                                     chains = 4,
                                     seed = 323,
                                     combine = FALSE)

# Get fitted values from model to data range/space
x_values <- data.frame(time = seq(0, 32, 0.1))
fit_diversity <- lapply(seq(1:length(group_names)), function(i) {
  group <- group_names[i]
  fit_line <- fitted(regression_diversity[[i]], newdata = x_values) %>% 
    as.data.frame() %>% 
    mutate(group = group,
           time = seq(0, 32, 0.1))
})
fit_diversity <- do.call("rbind", fit_diversity)

# Plot fitted model
gg_diversity_fit <- ggplot(fit_diversity, aes(x = time, y = Estimate, group = group, color = group)) +
  geom_ribbon(aes(ymax = Q97.5, ymin = Q2.5, fill = group),
              color = NA, alpha = 0.2) +
  geom_line(size = 0.3) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_y_continuous(breaks = seq(0, 3, 0.05)) +
  scale_color_manual(values = c("#F18805", line_color), name = "") +
  scale_fill_manual(values = c("#F18805", line_color), name = "") +
  xlab("Time since first article share (hrs)") +
  ylab("Ideological diversity of exposed users") +
  theme_ctokita()
gg_diversity_fit


####################
# Diversity of exposure at article level (per hour)
####################
# Calculate diversity of exposure at article level (per hour)
diveristy_article <- exposure_ideol %>% 
  filter(hour_bin >= 0) %>% 
  group_by(total_article_number, !!sym(grouping), hour_bin, ideology_bin) %>% 
  summarise(count = sum(count)) %>% 
  pivot_wider(names_from = ideology_bin, values_from = count) 
ideol_bins_cols <- grep("[-.0-9]+", names(diveristy_article)) #list of ideology bin columns
diversity_index <- data.frame(ideol_diversity = ideological_diversity(distance_matrix = ideol_distance_matrix, ideological_counts = diveristy_article[ , ideol_bins_cols]))
diveristy_article <- diveristy_article %>% 
  select(!!sym(grouping), total_article_number, hour_bin) %>% 
  cbind(diversity_index)

# Raw plot of tweets
gg_expos_diversity_raw <- diveristy_article %>% 
  filter(hour_bin < 31, 
         !is.na(ideol_diversity)) %>% 
  ggplot(., aes(x = hour_bin, y = ideol_diversity, group = total_article_number)) +
  geom_point(size = 0.3, alpha = 0.2, color = line_color) +
  theme_ctokita() + 
  facet_wrap(as.formula(paste("~", grouping)), 
             ncol = 1,
             strip.position = "right",
             scales = "fixed")
gg_expos_diversity_raw


# Fit bayesian regression to data
library(brms)
diversity_split <- diveristy_article %>%
  ungroup() %>%
  select(hour_bin, ideol_diversity, !!sym(grouping)) %>% 
  filter(!is.na(ideol_diversity))
if (grouping == "article_fc_rating") {
  diversity_split <- diversity_split %>% 
    split(.$article_fc_rating)
  group_names <- unique(diveristy_article$article_fc_rating)
} else if (grouping == "source_type") {
  diversity_split <- diversity_split %>% 
    split(.$source_type)
  group_names <- unique(diveristy_article$source_type)
}
regression_diversity <- brm_multiple(data = diversity_split,
                                     # formula = ideol_diversity ~ 1 + time + I(time^2),
                                     formula = ideol_diversity ~ 1 + hour_bin, 
                                     prior = c(prior(uniform(-10, 10), class = Intercept),
                                               prior(normal(0, 1), class = b),
                                               prior(normal(0, 1), class = sigma)),
                                     iter = 3000,
                                     warmup = 1000,
                                     chains = 4,
                                     seed = 323,
                                     combine = FALSE)

# Get fitted values from model to data range/space
x_values <- data.frame(hour_bin = seq(0, 32, 0.1))
fit_diversity <- lapply(seq(1:length(group_names)), function(i) {
  group <- group_names[i]
  fit_line <- fitted(regression_diversity[[i]], newdata = x_values) %>% 
    as.data.frame() %>% 
    mutate(group = group,
           hour_bin = seq(0, 32, 0.1))
})
fit_diversity <- do.call("rbind", fit_diversity)

# Plot fitted model
gg_diversity_fit <- ggplot(fit_diversity, aes(x = hour_bin, y = Estimate, group = group, color = group)) +
  geom_ribbon(aes(ymax = Q97.5, ymin = Q2.5, fill = group),
              color = NA, alpha = 0.2) +
  geom_line(size = 0.3) +
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  scale_y_continuous(breaks = seq(0, 3, 0.05)) +
  scale_color_manual(values = c("#F18805", line_color), name = "") +
  scale_fill_manual(values = c("#F18805", line_color), name = "") +
  xlab("Time since first article share (hrs)") +
  ylab("Ideological diversity of exposed users") +
  theme_ctokita()
gg_diversity_fit

