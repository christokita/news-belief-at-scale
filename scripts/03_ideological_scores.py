#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 16:15:11 2020

@author: ChrisTokita

SCRIPT:
Analyze ideological scores of tweeters
"""

####################
# Load libraries and packages, set important parameters, load data
####################
import pandas as pd
import numpy as np
import re
import os
import multiprocessing as mp
import datetime as dt
    
# high level directory (external HD or cluster storage)
#data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD


####################
# Add article-level information
####################
# Load tweets
tweets = pd.read_csv(data_directory + "data_derived/tweets/parsed_tweets.csv")

# Get source lean ratings
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")
articles = articles.rename(columns = {'total article number': 'total_article_number'})
articles['source_lean'] = articles['source'].str.replace('rss_|ct_', '')
lean_dict = {'con': 'C', 'lib': 'L', 'unclear': 'U'}
articles['source_lean'] = articles['source_lean'].map(lean_dict)
articles = articles[['total_article_number', 'source_lean']]

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

# Merge in ideological scores
labeled_tweeters = labeled_tweeters.merge(ideological_scores, how = 'left', on = 'user_id')


####################
# Write to file
####################
# Save
labeled_tweeters.to_csv(data_directory + "data_derived/tweets/all_tweets_labeled.csv", index = False)