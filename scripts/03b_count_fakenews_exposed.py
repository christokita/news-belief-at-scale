#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 10:44:52 2020

@author: ChrisTokita

SCRIPT:
Parse user data to determine the unique number and IDs of users (and user types).
This script uses a high performance computing cluster, assuming a slurm-based scheduling system.
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import os

# path to data
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/"
tweet_file = data_directory + "data_derived/tweets/parsed_tweets.csv"
path_to_exposed = data_directory + "data_derived/followers/processed_exposed_followers/"
path_nofollowers_fm_tweeters = data_directory + "data_derived/followers/nofollowers_fm_tweeters/"

# Desired outpaths for writing new data files
outpath = data_directory + "data_derived/"
outpath_follower_data = data_directory + "data_derived/followers/"
outpath_tweet_data = data_directory + "data_derived/tweets/"


####################
# Count up unique followers
####################
# Get unique IDs of tweeters
tweet_data = pd.read_csv(tweet_file)
tweeters = np.unique(tweet_data['user_id'])
tweeters = tweeters.astype(int)
count_tweeters = len(tweeters)

# Count up followers
exposed_files = sorted( os.listdir(path_to_exposed) )
num_files = len(exposed_files)
all_exposed = np.array([], dtype = int)
for file in exposed_files:
    data = np.loadtxt(path_to_exposed + file, skiprows = 1, dtype = int) #first row is header
    all_exposed = np.append(all_exposed, data)
    all_exposed = np.unique(all_exposed)
    del(data)
all_exposed_minus_tweeters = np.setdiff1d(all_exposed, tweeters)
count_all_exposed = all_exposed.shape[0]
count_all_exposed_minus_tweeters = all_exposed_minus_tweeters.shape[0]
all_exposed = pd.DataFrame(all_exposed, columns = ["user_id"])
all_exposed_minus_tweeters = pd.DataFrame(all_exposed_minus_tweeters, columns = ["user_id"])
    

# Count FM tweeters who had no followers
exposed_nofollower_files = sorted( os.listdir(path_nofollowers_fm_tweeters) )
no_followers = np.array([], dtype = int)
for file in exposed_nofollower_files:
    data = np.loadtxt(path_nofollowers_fm_tweeters + file, skiprows = 1, dtype = int) #first row is header
    no_followers = np.append(no_followers, data)
    no_followers = np.unique(no_followers)
    del(data)
count_no_followers = no_followers.shape[0]
no_followers = pd.DataFrame(no_followers, columns = ["user_id"])


####################
# Count up FM tweeters
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
for j in range(len(fm_articles)):
    link = fm_articles['link'].iloc[j]
    shortlink = fm_articles['short link'].iloc[j]
    try:
        has_full_link = tweets['urls_expanded'].str.contains(link)
        has_short_link = tweets['urls_expanded'].str.contains(shortlink)
    except:
        has_short_link = pd.Series(np.repeat(False, tweets.shape[0]))
    has_link = has_full_link | has_short_link #boolean operator to find which indices have one of the two possible links
    fm_tweets = fm_tweets.append(tweets[has_link])
    
# Get user IDs of tweeters of FM articles, 
fm_tweeters = fm_tweets['user_id']
fm_tweeters, fake_freq = np.unique(fm_tweeters, return_counts = True)
fm_tweeters =pd.DataFrame(fm_tweeters, columns = ['user_id'])
count_fm_tweeters = fm_tweeters.shape[0]


####################
# Save
####################
# Measure number of users and write out to file
unique_users = pd.DataFrame({'user_type': ["Tweeters", "FM tweeters", "Followers exposed to FM Articles", "Exposed followers excluding all weeters", "FM Tweeters w/o followers"], 
                             'count': [count_tweeters, count_fm_tweeters, count_all_exposed, count_all_exposed_minus_tweeters, count_no_followers]})
unique_users.to_csv(outpath + "exposed_user_summary.csv", index = False)
all_exposed.to_csv(outpath_follower_data + "exposed_fakenews_followers.csv", index = False)
no_followers.to_csv(outpath_follower_data + "fm_tweeter_nofollowers.csv", index = False)
fm_tweeters.to_csv(outpath_tweet_data + "unique_fm_tweeters.csv", index = False)