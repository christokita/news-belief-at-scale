#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 15:37:26 2020

@author: ChrisTokita

SCRIPT
Grab a few hundred moderates according to the ideology score to see whether ~0 corresponds to true moderates or not.
"""
import pandas as pd

 # high level directory (external HD or cluster storage)
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD


# Load tweet data, esnure in proper format
labeled_tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                             dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                                      'user_id': object, 'tweet_id': 'int64'})

# Output samples - raw scores
labeled_tweets['user_id_str'] = "\"" + labeled_tweets['user_id'] + "\""
moderate_sample = labeled_tweets[['user_id', 'user_id_str', 'user_name', 'user_ideology', 'user_ideol_extremity']].drop_duplicates()
moderate_sample = moderate_sample.nsmallest(200, 'user_ideol_extremity')
moderate_sample.to_csv(data_directory + "data_derived/ideological_scores/sample_moderate_users.csv", index = False)

# Output samples of normalized scores
# Normalize
upper_bound = 3.5
lower_bound= -3.5
max_val = max(labeled_tweets.user_ideology)
min_val = min(labeled_tweets.user_ideology)
labeled_tweets['normalized_ideology'] = (labeled_tweets['user_ideology'] - min_val) / (max_val - min_val)
labeled_tweets['normalized_ideology'] = ( (upper_bound - lower_bound) *  labeled_tweets['normalized_ideology'] ) + lower_bound

# Select examples
labeled_tweets['normalized_ideol_extremity'] = abs(labeled_tweets['normalized_ideology'])
norm_moderate_sample = labeled_tweets[['user_id', 'user_id_str', 'user_name', 'user_ideology', 'normalized_ideology', 'normalized_ideol_extremity']].drop_duplicates()
norm_moderate_sample = norm_moderate_sample.nsmallest(200, 'normalized_ideol_extremity')
norm_moderate_sample.to_csv(data_directory + "data_derived/ideological_scores/sample_moderate_users_normalized.csv", index = False)



import matplotlib.pyplot as plt
plt.scatter(labeled_tweets['user_ideology'], labeled_tweets['normalized_ideology'])
plt.hist(labeled_tweets['normalized_ideology'], bins = np.arange(-3.5, 4, 1))