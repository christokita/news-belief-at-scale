#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 15:14:12 2020

@author: ChrisTokita

SCRIPT:
Determine who shared fake news articles and who was potentially exposed to it. 
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import re
import os

# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/"

####################
# Parse articles to get false news articles
####################
# Get article IDs of false/misleading articles (as evaluated by fact checkers)
news_evaluations = pd.read_csv(data_directory + "data/articles/evaluations.csv")
fakenews_ids = news_evaluations["article_num"][news_evaluations['mode of FC'] == "FM"]

# Get URL set for FM articles
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")
fm_articles = articles[articles['total article number'].isin(fakenews_ids)]

# Load and parse tweets according to full and shortened URLs
tweets = pd.read_csv(data_directory + "data_derived/tweets/parsed_tweets.csv")

fm_tweets = pd.DataFrame(columns = tweets.columns)
for i, l in enumerate(fm_articles):
    link = fm_articles['link'].iloc[i]
    shortlink = fm_articles['short link'].iloc[i]
    try:
        has_full_link = tweets['urls_expanded'].str.contains(link)
        has_short_link = tweets['urls_expanded'].str.contains(shortlink)
    except:
        has_short_link = pd.Series(np.repeat(False, tweets.shape[0]))
    has_link = has_full_link | has_short_link #boolean operator to find which indices have one of the two possible links
    fm_tweets = fm_tweets.append(tweets[has_link])
    
    
####################
# Determine which users tweeted false articles
####################
# Get user IDs of tweeters of FM articles, 
fm_tweeters = fm_tweets['user_id']
fm_tweeters, fake_freq = np.unique(fm_tweeters, return_counts = True)

# plot frequency of fake news article sharing
import matplotlib.pyplot as plt  
plt.hist(fake_freq, log = True)


####################
# Determine how many users were potentially exposed to these articles
####################
# loop through files and get set of unique followers who were exposed to the article
follower_files = os.listdir(data_directory + "data/followers/")
follower_files = [file for file in follower_files if re.match('[0-9]', file)] #filter out hidden copies of same files
exposed_followers = np.array([])
no_follower_list = np.array([])
for user_id in fm_tweeters:
    regex = re.compile(r"[0-9].*_%s.csv" % str(user_id))
    file = list(filter(regex.match, follower_files))
    try:
        if len(file) > 1:
            print("WARNING: user_id = %d matches multiple follower list files." % user_id)
        follower_list = pd.read_csv(data_directory + "data/followers/" + file[0])
        exposed_followers = np.append(exposed_followers, follower_list)
        exposed_followers = np.unique(exposed_followers)
    except:
        no_follower_list = np.append(no_follower_list, user_id)
    
# Get unique IDs of tweeters
tweet_data = pd.read_csv(data_directory + "data_derived/tweets/parsed_tweets.csv")
tweeters = np.unique(tweet_data['user_id'])
exposed_followers = np.setdiff1d(exposed_followers, tweeters) #remove followers who are in tweeter set

# Save
exposed_followers = pd.DataFrame(exposed_followers, columns = ['user_id'])
exposed_followers.to_csv(follower_files + "data_derived/followers/followers_exposed_fakenews.csv")

no_follower_list = pd.DataFrame(no_follower_list, columns = ['user_id'])
no_follower_list.to_csv(follower_files + "data_derived/followers/nofollowersfound_exposed_fakenews.csv")
