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

all_followers = np.array([])
for file in follower_files:
    followers = pd.read_csv(path_to_followers + file)
    all_followers = np.append(all_followers, followers)
all_followers = np.unique(all_followers)

# Measure number of users and write out to file
unique_users = pd.DataFrame({'user_type': ["Tweeters", "Followers", "Total"], 'count': [len(tweeters), len(all_followers), len(tweeters) + len(all_followers)]})
all_followers = pd.DataFrame(all_followers, columns = ['follwer_id'])
tweeters = pd.DataFrame(tweeters, columns = ['tweeter_id'])

unique_users.to_csv(outpath_highlevel_data + "unique_user_summary.csv", index = False)
tweeters.to_csv(outpath_tweeters + "unique_tweeters.csv", index = False)
all_followers.to_csv(outpath_followers + "unique_followers.csv", index = False)
