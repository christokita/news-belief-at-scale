#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 16:15:11 2020

@author: ChrisTokita

SCRIPT:
Add user ideology and article-level data to tweet data.
"""

####################
# Load libraries and packages, set important parameters, load data
####################
import pandas as pd
import numpy as np
import re
import os
    
# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD


####################
# Add article-level information
####################
# Load tweets
tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_parsed.csv")

# Get source lean ratings
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")
articles = articles.rename(columns = {'total article number': 'total_article_number'})
articles['source_lean'] = articles['source'].str.replace('rss_|ct_', '')
lean_dict = {'con': 'C', 'lib': 'L', 'unclear': 'U'}
articles['source_lean'] = articles['source_lean'].map(lean_dict)

# Get source veracity (mainstream or fringe)
articles['source_type'] = articles['source'].str.replace('_con|_lib|_unclear', '')
articles['source_type'] = articles['source_type'].str.replace('rss', 'fringe')
articles['source_type'] = articles['source_type'].str.replace('ct', 'mainstream')

# Filter down article metadata to columns of interest
articles = articles[['total_article_number', 'source_lean', 'source_type']]

# Load and prepare article veracity evaluations
news_evaluations = pd.read_csv(data_directory + "data/articles/evaluations.csv")
news_evaluations = news_evaluations.rename(columns = {'article_num': 'total_article_number',
                                                      'mode of FC': 'article_fc_rating'})
news_evaluations = news_evaluations[['total_article_number', 'article_fc_rating']]
news_evaluations = news_evaluations.dropna()

# Load and prepare article lean ratings
article_ratings = pd.read_csv(data_directory + "data/articles/article_level.csv")
article_ratings = article_ratings.rename(columns = {'article_num': 'total_article_number', 
                                                    'main_tag': 'article_main_tag',
                                                    'lean': 'article_lean',
                                                    'other_tags': 'article_other_tags',
                                                    'con_feel': 'article_con_feel',
                                                    'lib_feel': 'article_lib_feel'})
article_ratings = article_ratings.drop(columns = ['daily article number', 'partisan_diff', 
                                                  'Unnamed: 8', 'Unnamed: 9', 
                                                  'Unnamed: 10', 'Unnamed: 11'])

# Merge in article and source info
labeled_tweeters = tweets.merge(articles, how = 'left', on = 'total_article_number')
labeled_tweeters = labeled_tweeters.merge(news_evaluations, how = 'left', on = 'total_article_number')
labeled_tweeters = labeled_tweeters.merge(article_ratings, how = 'left', on = 'total_article_number')

####################
# Add user-level information
####################
# Load and prepare ideological scores
ideological_scores = pd.read_csv(data_directory + "data_derived/ideological_scores/unique_tweeters_ideology_scores.csv")
ideological_scores = ideological_scores.rename(columns = {'id_str': 'user_id', 'pablo_score': 'user_ideology'})
ideological_scores['user_ideol_extremity'] = abs(ideological_scores['user_ideology'])

# Create broad ideological categories: Liberal, Moderate, Conservative
ideo_threshold = 1 #how far from zero is the cutoff for non-Moderates?
ideological_scores.loc[ideological_scores.user_ideology < -ideo_threshold, 'user_ideol_category'] = -1 #liberal
ideological_scores.loc[ideological_scores.user_ideology > ideo_threshold, 'user_ideol_category'] = 1 #conservative
ideological_scores.loc[ideological_scores['user_ideology'].between(-ideo_threshold, ideo_threshold, inclusive = True), 'user_ideol_category'] = 0 #moderate


# Merge in ideological scores
labeled_tweeters = labeled_tweeters.merge(ideological_scores, how = 'left', on = 'user_id')

####################
# Calculate missing ideologies for tweeters 
####################
# These were users who SMAPP did not already have calculated
# We will calculate their ideology using a simple mean of their friends

# Determine missing ideologies
missing_ideologies = labeled_tweeters[['user_id', 'user_ideology']].drop_duplicates()
missing_ideologies = missing_ideologies[pd.isna(missing_ideologies.user_ideology)]

# Load friend ideology scores
friend_ideologies = pd.read_csv(data_directory + "data_derived/ideological_scores/cleaned_friends_ideology_scores.csv")

# For each tweeter, calculate mean ideology of friends
friend_files = os.listdir(data_directory + "data/friends/")
#for i in np.arange(missing_ideologies.shape[0]):
for i in np.arange(300):

    # Get user ID and load followers
    user = missing_ideologies['user_id'].iloc[i]
    regex = re.compile(r"[0-9].*_%s.csv" % user)
    file = list(filter(regex.match, friend_files))
    
    # Parse friends and calculate mean. If no followers, next person
    if len(file) == 0: #no followers, no follower file
            next
    else:
        try:
            follower_list = np.genfromtxt(data_directory + "data/followers/" + file[0], dtype = int)
            follower_list = follower_list[1:len(follower_list)] #remove header, will raise error if empty
            friend_scores = friend_ideologies[friend_ideologies['user_id'].isin(follower_list)]
            avg_score = np.mean(friend_scores['pablo_score'])
            missing_ideologies['user_ideology'].iloc[i] = avg_score
        except:
            next
    
# Merge in new scores
labeled_tweeters = labeled_tweeters.set_index('user_id').user_ideology.fillna(missing_ideologies.set_index('user_id').user_ideology).reset_index()
labeled_tweeters['manual_ideology_score'] = labeled_tweeters['user_id'].isin(missing_ideologies['user_id'])

####################
# Write to file
####################
# Save
labeled_tweeters.to_csv(data_directory + "data_derived/tweets/tweets_labeled.csv", index = False)