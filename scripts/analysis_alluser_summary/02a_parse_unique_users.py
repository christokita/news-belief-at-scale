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
import re
import sys

# Get file numbers
i = int(sys.argv[1]) # get which chunk of the 

# path to data
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/"
tweet_file = data_directory + "data_derived/tweets/parsed_tweets.csv"
path_to_followers = data_directory + "data/followers/"

# Desired outpaths for writing new data files
outpath_highlevel_data = data_directory + "data_derived/"
outpath_tweeters = data_directory + "data_derived/tweets/"
outpath_followers = data_directory + "data_derived/followers/"


####################
# Unique user count
####################
# Get unique IDs of tweeters
tweet_data = pd.read_csv(tweet_file)
tweeters = np.unique(tweet_data['user_id'])

# Get unique IDs of followers
follower_files = os.listdir(path_to_followers)
follower_files = [file for file in follower_files if re.match('[0-9]', file)] #filter out hidden copies of same files
follower_files = follower_files[ (865*i) : (865*(i+1)) ] #there are 86,370 files and we want 100 chunks

followers = np.array([])
for file in follower_files:
    follower_list = pd.read_csv(path_to_followers + file)
    followers = np.append(followers, follower_list)
followers = np.unique(followers)
followers = followers.astype(int)
followers_excl_tweeters = np.setdiff1d(followers, tweeters) #remove followers who are in tweeter set

# Measure number of users and write out to file
tweeters = pd.DataFrame(tweeters, columns = ['tweeter_id'])
followers = pd.DataFrame(followers, columns = ['follower_id'])
followers_excl_tweeters = pd.DataFrame(followers_excl_tweeters, columns = ['follower_id'])

chunk_label = str(i).zfill(2)
tweeters.to_csv(outpath_tweeters + "unique_tweeters.csv", index = False)
followers.to_csv(outpath_followers + "processed_followers/unique_followers_" + chunk_label + ".csv", index = False)
followers_excl_tweeters.to_csv(outpath_followers + "processed_followers_excludingtweeters/unique_followers_excludingtweeters_" + chunk_label + ".csv", index = False)
