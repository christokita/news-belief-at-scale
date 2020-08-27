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

# Output samples
labeled_tweets['user_id_str'] = "\"" + labeled_tweets['user_id'] + "\""
moderate_sample = labeled_tweets[['user_id', 'user_id_str', 'user_name', 'user_ideology', 'user_ideol_extremity']].drop_duplicates()
moderate_sample = moderate_sample.nsmallest(200, 'user_ideol_extremity')
moderate_sample.to_csv(data_directory + "data_derived/ideological_scores/sample_moderate_users.csv", index = False)