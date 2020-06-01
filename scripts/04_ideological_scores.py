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
# Determine article sets: Liberal FM, Conservative FM, Liberal True, Conservative True
####################
# Load articles
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")

# Get article IDs of fake and true articles (as evaluated by fact checkers)
news_evaluations = pd.read_csv(data_directory + "data/articles/evaluations.csv")
fakenews_ids = news_evaluations["article_num"][news_evaluations['mode of FC'] == "FM"]
truenews_ids = news_evaluations["article_num"][news_evaluations['mode of FC'] == "T"]

# Grab articles by lean and veracity
true_con_articles = articles[articles['total article number'].isin(truenews_ids) & articles['source'].str.match('.*con')]
true_lib_articles = articles[articles['total article number'].isin(truenews_ids) & articles['source'].str.match('.*lib')]
fm_con_articles = articles[articles['total article number'].isin(fakenews_ids) & articles['source'].str.match('.*con')]
fm_lib_articles = articles[articles['total article number'].isin(fakenews_ids) & articles['source'].str.match('.*lib')]


####################
# Determine ideological distribution 
####################
# Load tweets and ideological scores
tweets = pd.read_csv(data_directory + "data_derived/tweets/parsed_tweets.csv")
ideological_scores = pd.read_csv(data_directory + "data_derived/ideological_scores/unique_fm_tweeters-scored.csv")
ideological_scores = ideological_scores.rename(columns = {'id_str': 'user_id'})
ideological_scores = ideological_scores.drop(columns = ['Unnamed: 0'])

# Subset out by article type
fm_con_tweeters = tweets[tweets['total_article_number'].isin(fm_con_articles['total article number'])]
fm_lib_tweeters = tweets[tweets['total_article_number'].isin(fm_lib_articles['total article number'])]

# Merge in ideological scores
fm_con_tweeters = fm_con_tweeters.merge(ideological_scores, how = 'left', on = 'user_id')
fm_lib_tweeters = fm_lib_tweeters.merge(ideological_scores, how = 'left', on = 'user_id')

# Label and join into master data set
fm_con_tweeters['article_type'] = 'false-conservative'
fm_lib_tweeters['article_type'] = 'false-liberal'
labeled_tweeters = fm_con_tweeters.append(fm_lib_tweeters)

# Plot
import matplotlib.pyplot as plt

plt.figure(figsize = (8,6))
bin_list = np.arange(-4.1, 4.1, 0.2)
plt.hist(fm_con_tweeters['ideology_score'], alpha = 0.5, label = "Con FM Articles", color = '#d54c54', bins = bin_list)
plt.hist(fm_lib_tweeters['ideology_score'], alpha = 0.5, label = "Lib FM Articles", color = '#006195', bins = bin_list)