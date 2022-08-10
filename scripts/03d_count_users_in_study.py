#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 11:56:51 2020

@author: ChrisTokita

SCRIPT
Count just the users in our study. We exclude the first 10 articles from analysis.

NOTE
Again, for counting purposes, int64 is faster than handling IDs as str/object.
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import os
import re

# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD
outpath = data_directory + "data_derived/"
   

####################
# Load tweets, count unique tweeters
#################### 
# Load tweet data, filter out unused articles
tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                     dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                              'user_id': 'int64', 'tweet_id': 'int64'})
tweets = tweets[tweets.total_article_number > 10].copy()

# Count up tweets and unique tweeters
n_tweets = tweets.shape[0]
n_tweeters = tweets['user_id'].unique().shape[0]

####################
# Count up followers
#################### 
follower_files = os.listdir(data_directory + "data/followers/")
follower_files = [file for file in follower_files if re.match('^[0-9]', file)] #filter out hidden copies of same files
followers = np.array([], dtype = 'int64')

# Loop through user IDs and add to list followers
for user_id in tweets['user_id'].unique():
    regex = re.compile(r"[0-9].*_%s.csv" % user_id)
    file = list(filter(regex.match, follower_files))
    if len(file) > 1:
        print("WARNING: user_id = %d matches multiple follower list files." % user_id)
    try:
        follower_list = np.genfromtxt(data_directory + "data/followers/" + file[0], dtype = 'int64')
        follower_list = follower_list[1:len(follower_list)] #remove header, will raise error if empty
        followers = np.append(followers, follower_list)
        del follower_list
    except:
        continue
        
# Get unique of lists
followers = np.unique(followers)        

# Count up unique followers
n_followers = len(followers)
del followers

####################
# Count up friends
#################### 
friend_files = os.listdir(data_directory + "data/friends/")
friend_files = [file for file in friend_files if re.match('^[0-9]', file)] #filter out hidden copies of same files
friends = np.array([], dtype = 'int64')

# Loop through user IDs and add to list followers
for user_id in tweets['user_id'].unique():
    regex = re.compile(r"[0-9].*_%s.csv" % user_id)
    file = list(filter(regex.match, friend_files))
    if len(file) > 1:
        print("WARNING: user_id = %d matches multiple follower list files." % user_id)
    try:
        friend_list = np.genfromtxt(data_directory + "data/friends/" + file[0], dtype = 'int64')
        friend_list = friend_list[1:len(friend_list)] #remove header, will raise error if empty
        friends = np.append(friends, friend_list)
        del friend_list
    except:
        continue
        
# Get unique of lists
friends = np.unique(friends)        

# Count up unique followers
n_friends = len(friends)
del friends


####################
# Ouatput summary
#################### 
unique_users = pd.DataFrame({'study_object': ["Tweets", "Tweeters", "Followers", "Friends"], 
                             'unique_count': [n_tweets,
                                              n_tweeters, 
                                              n_followers,
                                              n_friends]})
unique_users.to_csv(outpath + "summary_users_in_study.csv", index = False)