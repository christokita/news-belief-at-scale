#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 10:44:52 2020

@author: ChrisTokita
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import os
import re

# path to data
tweet_file = "/Volumes/CKT-DATA/data_derived/tweets/parsed_tweets.csv"
path_to_followers = "/Volumes/CKT-DATA/data/followers/"

# Desired outpaths for writing new data files
outpath_highlevel_data = "/Volumes/CKT-DATA/data_derived/"
outpath_tweeters = "/Volumes/CKT-DATA/data_derived/tweets/"
outpath_followers = "/Volumes/CKT-DATA/data_derived/followers/"


####################
# Unique user count
####################
# Get unique IDs of tweeters
tweet_data = pd.read_csv(tweet_file)
tweeters = np.unique(tweet_data['user_id'])

# Get unique IDs of followers
follower_files = os.listdir(path_to_followers)
follower_files = [file for file in follower_files if re.match('[0-9]', file)] #filter out hidden copies of same files

followers = np.array([])
for file in follower_files:
    follower_list = pd.read_csv(path_to_followers + file)
    followers = np.append(followers, follower_list)
followers = np.unique(followers)
followers = followers.astype(int)
followers_excl_tweeters = np.setdiff1d(followers, tweeters) #remove followers who are in tweeter set



# Measure number of users and write out to file
unique_users = pd.DataFrame({'user_type': ["Tweeters", "Followers", "Followers excluding tweeters", "Total"], 
                             'count': [len(tweeters), len(followers), len(followers_excl_tweeters), len(tweeters) + len(followers_excl_tweeters)]})
tweeters = pd.DataFrame(tweeters, columns = ['tweeter_id'])
followers = pd.DataFrame(followers, columns = ['follower_id'])
followers_excl_tweeters = pd.DataFrame(followers_excl_tweeters, columns = ['follower_id'])

unique_users.to_csv(outpath_highlevel_data + "unique_user_summary.csv", index = False)
tweeters.to_csv(outpath_tweeters + "unique_tweeters.csv", index = False)
followers.to_csv(outpath_followers + "unique_followers.csv", index = False)
followers_excl_tweeters.to_csv(outpath_followers + "unique_followers_excludingtweeters.csv", index = False)
