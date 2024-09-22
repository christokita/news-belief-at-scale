########################################
#
# PLOT: Ideology check - self-reported ideology (categorical) vs estimated Twitter ideology (continuous)
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
ideology_scores_path <- '/Volumes/CKT-DATA/news-belief-at-scale/data/ideology_check/selfreported_vs_pabloscore.csv' 
outpath <- 'output/ideology_basis/'

# Color palette
plot_color <- "#495867"
ideol_pal_small <- rev(brewer.pal(5, "RdBu"))
ideol_pal_small[3] <- "#e0e0e0"

ideol_pal_large <- rev(brewer.pal(7, "RdBu"))
ideol_pal_large[4] <- "#e0e0e0"



####################
# Load data 
####################
# Read in self-reported vs. pablo score data
ideologies <- read.csv(ideology_scores_path, header = TRUE) %>% 
  mutate(ideology_category = ideo5)


####################
# PLOT: scatter Twitter vs. self-reported ideology
####################
# Calculate mean by ideology category
ideology_means <- ideologies %>% 
  group_by(ideology_category) %>% 
  summarise(mean_pablo_score = mean(pablo_score),
            sd_pablo_score = sd(pablo_score),
            se_pablo_score = sd(pablo_score) / sqrt(length(pablo_score)))

# Create composite mean for "somewhat X" categories
somewhat_liberal <- ideology_means %>% 
  filter(ideology_category %in% c("Liberal", "Moderate")) %>% 
  summarise(mean_pablo_score = mean(mean_pablo_score)) %>% 
  mutate(ideology_category = "Somewhat liberal",
         sd_pablo_score = NA,
         se_pablo_score = NA)

somewhat_conservative <- ideology_means %>% 
  filter(ideology_category %in% c("Conservative", "Moderate")) %>% 
  summarise(mean_pablo_score = mean(mean_pablo_score)) %>% 
  mutate(ideology_category = "Somewhat conservative",
         sd_pablo_score = NA,
         se_pablo_score = NA)

ideology_means <- rbind(ideology_means, somewhat_liberal)
ideology_means <- rbind(ideology_means, somewhat_conservative) %>% 
  mutate(ideology_category = factor(ideology_category, levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative", "Not sure"))) %>% 
  arrange(ideology_category)

ideologies <- ideologies %>% 
  mutate(ideology_category = factor(ideology_category, levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative", "Not sure")))

# Plot
gg_ideology_comparison <- ggplot() +
  geom_hline(yintercept = 0, linetype = 'dotted', size = 0.3) +
  # Have to double plot group means to force factors to plot in order
  geom_point(data = ideology_means, aes(x = ideology_category, y = mean_pablo_score, fill = ideology_category),
             size = 3, 
             shape = 21,
             color = 'white') +  
  geom_point(data = ideologies, aes(x = ideology_category, y = pablo_score, color = ideology_category),
             size = 1.5, 
             stroke = 0, 
             alpha = 0.25,
             position = position_jitter(width = 0.1, height = 0)) +
  geom_point(data = ideology_means, aes(x = ideology_category, y = mean_pablo_score, fill = ideology_category),
             size = 3,
             shape = 21,
             color = 'white') +
  xlab('Self-reported ideology') +
  ylab('Inferred ideology from Twitter') +
  scale_x_discrete(labels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative", "Not sure")) +
  scale_y_continuous(breaks = seq(-4, 4, 1), limits = c(-4, 4), expand = c(0, 0)) +
  scale_color_manual(values = c(ideol_pal_small, 'black')) +
  scale_fill_manual(values = c(ideol_pal_large, 'black')) +
  theme_ctokita() +
  theme(legend.position = 'none',
        aspect.ratio = 0.6,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.03))

gg_ideology_comparison
ggsave(gg_ideology_comparison, filename = paste0(outpath, 'empirical_ideology_means.pdf'), width = 90, height = 70, units = 'mm', dpi = 400)


####################
# PLOT: scatter Twitter vs. self-reported ideology with resampled data for "Somewhat" category
####################
# Get IDs of users to reassign to "somewhat category"
somewhat_liberal_samples <- ideologies %>% 
  filter(ideology_category %in% c("Liberal", "Moderate")) %>% 
  group_by(ideology_category) %>% 
  sample_frac(0.33)

somewhat_conservative_samples <- ideologies %>% 
  filter(ideology_category %in% c("Conservative", "Moderate")) %>% 
  filter(!smappid %in% somewhat_liberal_samples$smappid) %>% 
  group_by(ideology_category) %>% 
  sample_frac(0.33)

# Reassign selected users to "Somewhat" categories
ideologies_synthetic <- ideologies %>% 
  mutate(data_type = "Actual category")
ideologies_synthetic$ideology_category[ideologies_synthetic$smappid %in% somewhat_liberal_samples$smappid] <- "Somewhat liberal"
ideologies_synthetic$data_type[ideologies_synthetic$smappid %in% somewhat_liberal_samples$smappid] <- "Reassigned category"

ideologies_synthetic$ideology_category[ideologies_synthetic$smappid %in% somewhat_conservative_samples$smappid] <- "Somewhat conservative"
ideologies_synthetic$data_type[ideologies_synthetic$smappid %in% somewhat_conservative_samples$smappid] <- "Reassigned category"


ideologies_synthetic <- ideologies_synthetic %>% 
  mutate(ideology_category = factor(ideology_category, levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative", "Not sure")))

# Calculate mean by ideology category
ideology_means_synthetic <- ideologies_synthetic %>% 
  group_by(ideology_category) %>% 
  summarise(mean_pablo_score = mean(pablo_score),
            sd_pablo_score = sd(pablo_score),
            se_pablo_score = sd(pablo_score) / sqrt(length(pablo_score))) %>% 
  mutate(ideology_category = factor(ideology_category, levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative", "Not sure")))

# Plot
gg_ideology_comparison_synthetic <- ggplot() +
  geom_hline(yintercept = 0, linetype = 'dotted', size = 0.3) +
  # Have to double plot group means to force factors to plot in order
  geom_point(data = ideology_means_synthetic, aes(x = ideology_category, y = mean_pablo_score, fill = ideology_category),
             size = 3, 
             shape = 21,
             color = 'white') +  
  geom_point(data = ideologies_synthetic, aes(x = ideology_category, y = pablo_score, color = ideology_category),
             size = 1.5, 
             stroke = 0, 
             alpha = 0.25,
             position = position_jitter(width = 0.1, height = 0)) +
  geom_point(data = ideology_means_synthetic, aes(x = ideology_category, y = mean_pablo_score, fill = ideology_category),
             size = 3,
             shape = 21,
             color = 'white') +
  xlab('Self-reported ideology') +
  ylab('Inferred ideology from Twitter') +
  scale_x_discrete(labels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative", "Not sure")) +
  scale_y_continuous(breaks = seq(-4, 4, 1), limits = c(-4, 4), expand = c(0, 0)) +
  scale_color_manual(values = c(ideol_pal_large, 'black')) +
  scale_fill_manual(values = c(ideol_pal_large, 'black')) +
  theme_ctokita() +
  theme(legend.position = 'none',
        aspect.ratio = 0.6,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.03))

gg_ideology_comparison_synthetic
ggsave(gg_ideology_comparison_synthetic, filename = paste0(outpath, 'synthetic_ideology_means.pdf'), width = 90, height = 70, units = 'mm', dpi = 400)


####################
# Calculate bins using these empirical categories
####################
empirical_bins <- ideology_means %>% 
  filter(ideology_category != "Not sure") %>% 
  mutate(bin_edge_upper = (mean_pablo_score + lead(mean_pablo_score)) / 2,
         bin_edge_lower = (mean_pablo_score + lag(mean_pablo_score)) / 2,
         bin_size = bin_edge_upper - bin_edge_lower)

empirical_bins_synethic <- ideology_means_synthetic %>% 
  filter(ideology_category != "Not sure") %>% 
  mutate(bin_edge_upper = (mean_pablo_score + lead(mean_pablo_score)) / 2,
         bin_edge_lower = (mean_pablo_score + lag(mean_pablo_score)) / 2,
         bin_size = bin_edge_upper - bin_edge_lower)


####################
# PLOT: confusion matrix self-reported ideology vs. Twitter ideology with resampled data for "Somewhat" category
####################
# Bin twitter ideology
# Categorical bins will be both those used in the paper and empirical bins
binned_ideologies_synthetic <- ideologies_synthetic %>% 
  filter(ideology_category != "Not sure") %>% 
  mutate(twitter = ifelse(pablo_score <= -2.5, "Very liberal",
                                            ifelse(-2.5 < pablo_score & pablo_score <= -1.5, "Liberal",
                                                   ifelse(-1.5 < pablo_score & pablo_score <= -0.5, "Somewhat liberal",
                                                          ifelse(-0.5 < pablo_score & pablo_score <= 0.5, "Moderate",
                                                                 ifelse(0.5 < pablo_score & pablo_score <= 1.5, "Somewhat conservative",
                                                                        ifelse(1.5 < pablo_score & pablo_score <= 2.5, "Conservative",
                                                                               ifelse(2.5 < pablo_score, "Very conservative", NA)))))))
         ) %>%
  # mutate(twitter_empirical = ifelse(pablo_score <= empirical_bins_synethic$bin_edge_upper[empirical_bins_synethic$ideology_category == "Very liberal"], "Very liberal",
  #                                         ifelse(empirical_bins_synethic$bin_edge_lower[empirical_bins_synethic$ideology_category == "Liberal"] < pablo_score & pablo_score <= empirical_bins_synethic$bin_edge_upper[empirical_bins_synethic$ideology_category == "Liberal"], "Liberal",
  #                                                ifelse(empirical_bins_synethic$bin_edge_lower[empirical_bins_synethic$ideology_category == "Somewhat liberal"] < pablo_score & pablo_score <= empirical_bins_synethic$bin_edge_upper[empirical_bins_synethic$ideology_category == "Somewhat liberal"], "Somewhat liberal",
  #                                                       ifelse(empirical_bins_synethic$bin_edge_lower[empirical_bins_synethic$ideology_category == "Moderate"] < pablo_score & pablo_score <= empirical_bins_synethic$bin_edge_upper[empirical_bins_synethic$ideology_category == "Moderate"], "Moderate",
  #                                                              ifelse(empirical_bins_synethic$bin_edge_lower[empirical_bins_synethic$ideology_category == "Somewhat conservative"] < pablo_score & pablo_score <= empirical_bins_synethic$bin_edge_upper[empirical_bins_synethic$ideology_category == "Somewhat conservative"], "Somewhat conservative",
  #                                                                     ifelse(empirical_bins_synethic$bin_edge_lower[empirical_bins_synethic$ideology_category == "Conservative"] < pablo_score & pablo_score <= empirical_bins_synethic$bin_edge_upper[empirical_bins_synethic$ideology_category == "Conservative"], "Conservative",
  #                                                                            ifelse(empirical_bins_synethic$bin_edge_lower[empirical_bins_synethic$ideology_category == "Very conservative"] < pablo_score, "Very conservative", NA)))))))
  #        ) %>%
  rename(selfreported=ideology_category) %>% 
  mutate(twitter = factor(twitter, levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative")),
         selfreported_ = factor(selfreported, levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative")))

# Calculate confusion matrix
ideology_confusion_matrix <- binned_ideologies_synthetic %>% 
  group_by(selfreported, twitter) %>% 
  count() %>% 
  ungroup() %>% 
  group_by(selfreported) %>% 
  mutate(proportion = n / sum(n))

# Fill in missing values of confusion matrix
ideology_levels = factor(c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative"), levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative"))

ideology_confusion_matrix <- ideology_confusion_matrix %>% 
  merge(expand.grid(selfreported = ideology_levels, twitter = ideology_levels), by = c("selfreported", "twitter"), all = TRUE)

# Plot confusion matrix
gg_ideology_confusion_matrix <- ggplot(data = ideology_confusion_matrix, aes(x = twitter, y = selfreported, fill = proportion)) +
  geom_tile() +
  # scale_fill_viridis_c(name = "Proportion within\nself-reported category", option = "inferno", na.value = "black") +
  scale_fill_gradientn(name = "Proportion within\nself-reported category", colours = brewer.pal(5, "Purples"), na.value = "#f2f0f7") +
  scale_x_discrete(expand = c(0, 0), 
                   position = 'top',
                   labels = c("L+", "L", "L-", "M", "C-", "C", "C+")) +
  scale_y_discrete(expand = c(0, 0), 
                   limits = rev(ideology_levels),
                   labels = c("Very conservative ( C+ )",  "Conservative ( C   )", "Somewhat conservative ( C- )", "Moderate ( M  )", "Somewhat Liberal ( L- )", "Liberal ( L   )", "Very liberal ( L+ )")) +
  xlab("Ideology inferred from twitter") +
  ylab("Self-reported ideology") +
  theme_ctokita() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0),
        panel.border = element_rect(linewidth = 0.6, fill = NA),
        axis.line = element_blank())

gg_ideology_confusion_matrix
ggsave(gg_ideology_confusion_matrix, filename = paste0(outpath, 'ideology_confusion_matrix.pdf'), width = 90, height = 45, units = 'mm', dpi = 400)

# CALCULATE: left-right divide accuracy
ideology_mixup <- ideology_confusion_matrix %>% 
  mutate(selfreported = tolower( gsub("Very |Somewhat ", "", selfreported) ),
         twitter = tolower( gsub("Very |Somewhat ", "", twitter) )) %>%
  group_by(selfreported, twitter) %>%
  summarise(n = sum(n, na.rm = TRUE)) %>% 
  mutate(proportion = n / sum(n))


####################
# PLOT: belief agreement matrix among ideology categories in surveys
####################
# We want to compare this to the confusion matrix to get a sense of whether miscategorization
# within left- and right-leaning categories is hugely impact

# Load belief data
article_belief <- read.csv('/Volumes/CKT-DATA/news-belief-at-scale/data/article_belief/response_distribution.csv') %>% 
  filter(total_article_number >= 11,
         fc_rating %in% c("T", "FM")) %>% 
  select(total_article_number, fc_rating, contains("_t")) %>% 
  select(-art_type, -Haven.t.thought.much.about.it_t) %>% #drop article type and "unsure" category
  gather(key = "ideology_category", value = "belief_rate", -total_article_number, -fc_rating) %>% 
  mutate(ideology_category = gsub("_t", "", ideology_category),
         ideology_category = gsub("\\.", " ", ideology_category),
         ideology_category = gsub("Moderate  Middle of the road", "Moderate", ideology_category),
         ideology_category = gsub("Extremely Liberal", "Very liberal", ideology_category),
         ideology_category = gsub("Extremely Conservative", "Very conservative", ideology_category),
         ideology_category = gsub("Slightly Liberal", "Somewhat liberal", ideology_category),
         ideology_category = gsub("Slightly Conservative", "Somewhat conservative", ideology_category)) %>% 
  mutate(ideology_category = factor(ideology_category, levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative")),
         fc_rating = ifelse(fc_rating == "T", "True news", ifelse(fc_rating == "FM", "False/Misleading news", 
                                                                  ifelse(fc_rating == "CND", "Borderline", 
                                                                         ifelse(fc_rating == "No Mode!", "No mode", fc_rating))))) %>% 
  arrange(total_article_number, ideology_category)

# Calculate pairwise distance
belief_distance <- article_belief %>% 
  group_by(total_article_number, fc_rating) %>%
  do({
    total_article_number <- unique(.$total_article_number)
    fc_rating <- unique(.$fc_rating)
    ideologies <- .$ideology_category
    values <- .$belief_rate
    dist_matrix <- as.matrix(dist(values, diag = TRUE, upper = TRUE))
    diag(dist_matrix) <- 0  # Add self-distances (zeros)
    pairwise_df <- as.data.frame(dist_matrix)
    colnames(pairwise_df) <- ideologies
    pairwise_df$ideology1 <- ideologies
    pairwise_df <- pairwise_df %>%
      pivot_longer(-ideology1, names_to = 'ideology2', values_to = 'distance')
    pairwise_df$total_article_number <- total_article_number
    pairwise_df$fc_rating <- fc_rating
    pairwise_df
  }) %>%
  ungroup() %>%
  mutate(ideology1 = factor(ideology1, levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative")),
         ideology2 = factor(ideology2, levels = c("Very liberal", "Liberal", "Somewhat liberal", "Moderate", "Somewhat conservative", "Conservative", "Very conservative"))) %>% 
  select(total_article_number, fc_rating, ideology1, ideology2, distance) %>% 
  arrange(total_article_number, ideology1, ideology2)

avg_belief_distance <- belief_distance %>% 
  group_by(fc_rating, ideology1, ideology2) %>% 
  summarise(avg_distance = mean(distance))

# Plot distance matrix
dist_colors <- rev(brewer.pal(5, "PuOr"))

gg_ideology_belief_distance <- ggplot(data = avg_belief_distance, aes(x = ideology1, y = ideology2, fill = avg_distance)) +
  geom_tile() +
  scale_fill_gradientn(name = "Avg. difference in\narticle belief rate", colors = dist_colors) +
  scale_x_discrete(expand = c(0, 0), 
                   position = "top",
                   labels = c("L+", "L", "L-", "M", "C-", "C", "C+")) +
  scale_y_discrete(expand = c(0, 0), 
                   limits = rev(ideology_levels),
                   labels = rev(c("L+", "L", "L-", "M", "C-", "C", "C+"))) +
  xlab("Self-reported ideology") +
  ylab("Self-reported ideology") +
  theme_ctokita() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 0.6, fill = NA),
        axis.line = element_blank(),
        strip.placement = 'outside') +
  facet_grid(~fc_rating,
             switch = 'both')
         
gg_ideology_belief_distance
ggsave(gg_ideology_belief_distance, filename = paste0(outpath, 'ideology_belief_agreement_matrix.pdf'), width = 90, height = 45, units = 'mm', dpi = 400)
