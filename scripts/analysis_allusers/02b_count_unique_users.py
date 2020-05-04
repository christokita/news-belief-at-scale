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
path_to_followers = data_directory + "data_derived/followers/processed_followers/"
path_to_followers_excl_tweeters = data_directory + "data_derived/followers/processed_followers_excludingtweeters/"

# Desired outpaths for writing new data files
outpath = data_directory + "data_derived/"
outpath_follower_data = data_directory + "data_derived/followers/"



####################
# Count up unique followers
####################
# Get unique IDs of tweeters
tweet_data = pd.read_csv(tweet_file)
tweeters = np.unique(tweet_data['user_id'])
count_tweeters = len(tweeters)

# Count up followers
count_all_followers = 0
follower_files = sorted( os.listdir(path_to_followers) )
num_files = len(follower_files)
all_followers = np.array([], dtype = int)
for file in follower_files:
    data = np.loadtxt(path_to_followers + file, skiprows = 1, dtype = int) #first row is header
    all_followers = np.append(all_followers, data)
    all_followers = np.unique(all_followers)
    del(data)
count_all_followers = all_followers.shape[0]
all_followers = pd.DataFrame(all_followers, columns = ["follower_id"])
    

# Count up followers, excluding tweeters
follower_excl_tweeters_files = sorted( os.listdir(path_to_followers_excl_tweeters) )
just_followers = np.array([], dtype = int)
for file in follower_excl_tweeters_files:
    data = np.loadtxt(path_to_followers_excl_tweeters + file, skiprows = 1, dtype = int) #first row is header
    just_followers = np.append(just_followers, data)
    just_followers = np.unique(just_followers)
    del(data)
count_just_followers = just_followers.shape[0]
just_followers = pd.DataFrame(just_followers, columns = ["follower_id"])

# Measure number of users and write out to file
unique_users = pd.DataFrame({'user_type': ["Tweeters", "Followers", "Followers excluding tweeters", "Total"], 
                             'count': [count_tweeters, count_all_followers, count_just_followers, count_tweeters+count_just_followers]})
unique_users.to_csv(outpath + "unique_user_summary.csv", index = False)
all_followers.to_csv(outpath_follower_data + "unique_followers.csv", index = False)
just_followers.to_csv(outpath_follower_data + "unique_followers_excl_tweeters.csv", index = False)