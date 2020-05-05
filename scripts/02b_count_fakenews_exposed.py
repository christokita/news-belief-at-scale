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
path_to_tweets = data_directory + "data_derived/tweets/"
path_to_exposed = data_directory + "data_derived/followers/processed_exposed_followers/"
path_nofollowers_fm_tweeters = data_directory + "data_derived/followers/nofollowers_fm_tweeters/"
path_to_fm_friends = data_directory + "data_derived/friends/processed_fm_tweeter_friends/"
path_to_nofriends_fm_tweeters = data_directory + "data_derived/friends/nofriends_fm_tweeters/"

# Desired outpaths for writing new data files
outpath = data_directory + "data_derived/"
outpath_follower_data = data_directory + "data_derived/followers/"
outpath_tweet_data = data_directory + "data_derived/tweets/"
outpath_friend_data = data_directory + "data_derived/friends/"


####################
# Count up FM tweeters
####################
# Load FM tweets
fm_tweets = pd.read_csv(path_to_tweets + "FM_tweets.csv")
    
# Get user IDs of tweeters of FM articles, 
fm_tweeters = fm_tweets['user_id'].astype(int)
fm_tweeters, fake_freq = np.unique(fm_tweeters, return_counts = True)
fm_tweeters = pd.DataFrame(fm_tweeters, columns = ['user_id'])
count_fm_tweeters = fm_tweeters.shape[0]

# Save FM tweeters
fm_tweeters.to_csv(outpath_tweet_data + "unique_fm_tweeters.csv", index = False)


####################
# Count up unique followers exoised to FM news
####################
# Get unique IDs of tweeters
tweet_data = pd.read_csv(path_to_tweets + "parsed_tweets.csv")
tweeters = tweet_data['user_id'].astype(int)
tweeters = np.unique(tweeters)
count_tweeters = len(tweeters)

# Count up followers
exposed_files = sorted( os.listdir(path_to_exposed) )
num_files = len(exposed_files)
all_exposed = np.array([], dtype = int)
for file in exposed_files:
    data = np.genfromtxt(path_to_exposed + file, skip_header = 1, dtype = int) #first row is header
    all_exposed = np.append(all_exposed, data)
    all_exposed = np.unique(all_exposed)
    del data
all_exposed_minus_tweeters = np.setdiff1d(all_exposed, fm_tweeters)
count_all_exposed = all_exposed.shape[0]
count_all_exposed_minus_tweeters = all_exposed_minus_tweeters.shape[0]
all_exposed = pd.DataFrame(all_exposed, columns = ["user_id"], dtype = int)
all_exposed_minus_tweeters = pd.DataFrame(all_exposed_minus_tweeters, columns = ["user_id"], dtype = int)
    

# Count FM tweeters who had no followers
exposed_nofollower_files = sorted( os.listdir(path_nofollowers_fm_tweeters) )
no_followers = np.array([], dtype = int)
for file in exposed_nofollower_files:
    data = np.genfromtxt(path_nofollowers_fm_tweeters + file, dtype = int) #first row is header
    try:
        data = data[1:len(data)] #if empty, will return error
        no_followers = np.append(no_followers, data)
        no_followers = np.unique(no_followers)
        del data
    except:
        next
count_no_followers = no_followers.shape[0]
no_followers = pd.DataFrame(no_followers, columns = ["user_id"], dtype = int)

# Save
all_exposed.to_csv(outpath_follower_data + "exposed_fakenews_followers.csv", index = False)
all_exposed_minus_tweeters.to_csv(outpath_follower_data + "exposed_fakenews_followers_excl_fmtweeters.csv", index = False)
no_followers.to_csv(outpath_follower_data + "fm_tweeter_nofollowers.csv", index = False)
del all_exposed, all_exposed_minus_tweeters, no_followers

####################
# Count up friends of FM tweeters
####################
# Count up followers
friend_files = sorted( os.listdir(path_to_fm_friends) )
num_files = len(friend_files)
all_friends = np.array([], dtype = int)
for file in friend_files:
    data = np.genfromtxt(path_to_fm_friends + file, skip_header = 1, dtype = int) #first row is header
    all_friends = np.append(all_friends, data)
    all_friends = np.unique(all_friends)
    del data
all_friends_minus_tweeters = np.setdiff1d(all_friends, fm_tweeters)
count_all_friends = all_friends.shape[0]
count_all_friends_minus_tweeters = all_friends_minus_tweeters.shape[0]
all_friends = pd.DataFrame(all_friends, columns = ["user_id"], dtype = int)
all_friends_minus_tweeters = pd.DataFrame(all_friends_minus_tweeters, columns = ["user_id"], dtype = int)
    

# Count FM tweeters who have no friends
fm_tweeters_nofriends_files = sorted( os.listdir(path_to_nofriends_fm_tweeters) )
no_friends = np.array([], dtype = int)
for file in fm_tweeters_nofriends_files:
    data = np.genfromtxt(path_to_nofriends_fm_tweeters + file, dtype = int) #first row is header
    try:
        data = data[1:len(data)] #if empty, will return error
        no_friends = np.append(no_friends, data)
        no_friends = np.unique(no_friends)
        del data
    except:
        next
count_no_friends = no_friends.shape[0]
no_friends = pd.DataFrame(no_friends, columns = ["user_id"], dtype = int)

# Save
all_friends.to_csv(outpath_friend_data + "friends_fm_tweeters.csv", index = False)
all_friends_minus_tweeters.to_csv(outpath_friend_data + "friends_fm_tweeters_excl_fmtweeters.csv", index = False)
no_friends.to_csv(outpath_friend_data + "fm_tweeter_nofriends.csv", index = False)
del all_friends, all_friends_minus_tweeters, no_friends

####################
# Summarize
####################
# Measure number of users and write out to file
unique_users = pd.DataFrame({'user_type': ["Tweeters", "FM tweeters", 
                                           "Followers exposed to FM Articles", "Exposed followers excluding FM tweeters", "FM tweeters w/o followers",
                                           "Friends of FM tweeters", "Friends of FM tweeters excluding FM tweeters", "FM tweeters w/o friends"], 
                             'count': [count_tweeters, count_fm_tweeters, 
                                       count_all_exposed, count_all_exposed_minus_tweeters, count_no_followers,
                                       count_all_friends, count_all_friends_minus_tweeters, count_no_friends]})
unique_users.to_csv(outpath + "exposed_user_summary.csv", index = False)
