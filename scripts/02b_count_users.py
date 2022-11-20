#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: `02b_count_users.py`
Date: April 14, 2020
Author: Chris Tokita
Purpose: Count ALL (even those not used in study) user lists in our dataset---i.e., tweeters, followers, friends.

Data In: 
    `<data storage location>/data_derived/tweets/tweets_parsed.csv`: Consolidated tweet dataframe from `01_parse_tweets.py`.
    `<data storage location>/data_derived/followers/`: consolidated follower lists from `02a_parse_users.py` containing Twitter user IDs of followers.
    `<data storage location>/data_derived/friends/`: consolidated friend lists from `02a_parse_users.py` containing Twitter user IDs of friends.

Data Out: `<external harddrive>/data_derived/summary_unique_users.csv`
    CSV file summarizing the number of unique tweeters, friends, and followers in our complete dataset, 
    before filtering out data that won't be used in study.

Machine: High-performance computing cluster
    This script is batched to the cluster using `slurm_scripts/count_users.cmd`
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import os

# path to data
data_directory = "/scratch/gpfs/ctokita/news-belief-at-scale/"
path_to_tweets = data_directory + "data_derived/tweets/"
path_to_followers = data_directory + "data_derived/followers/processed_followers/"
path_tweeters_nofollowers = data_directory + "data_derived/followers/tweeters_nofollowers/"
path_to_friends = data_directory + "data_derived/friends/processed_friends/"
path_to_tweeters_nofriends = data_directory + "data_derived/friends/tweeters_nofriends/"

# Desired outpaths for writing new data files
outpath = data_directory + "data_derived/"
outpath_follower_data = data_directory + "data_derived/followers/"
outpath_tweet_data = data_directory + "data_derived/tweets/"
outpath_friend_data = data_directory + "data_derived/friends/"


####################
# Count up tweeters
####################
# Load tweets
tweets = pd.read_csv(path_to_tweets + "tweets_parsed.csv",
                     dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                              'user_id': 'int64', 'tweet_id': 'int64'})
tweeters = tweets['user_id']
    
# Get user IDs of tweeters
tweeters, tweet_freq = np.unique(tweeters, return_counts = True)
tweeters = pd.DataFrame(tweeters, columns = ['user_id'], dtype = 'int64')
count_tweeters = tweeters.shape[0]

# Save tweeters
tweeters.to_csv(outpath_tweet_data + "unique_tweeters.csv", index = False)


####################
# Count up unique followers of tweeters
####################
# Count up followers
follower_files = sorted( os.listdir(path_to_followers) )
num_files = len(follower_files)
all_followers = np.array([], dtype = 'int64')
for file in follower_files:
    data = np.genfromtxt(path_to_followers + file, skip_header = 1, dtype = 'int64') #first row is header
    all_followers = np.append(all_followers, data)
    all_followers = np.unique(all_followers)
    del data
followers_minus_tweeters = np.setdiff1d(all_followers, tweeters)
count_followers = all_followers.shape[0]
count_followers_minus_tweeters = followers_minus_tweeters.shape[0]
all_followers = pd.DataFrame(all_followers, columns = ["user_id"], dtype = 'int64')
followers_minus_tweeters = pd.DataFrame(followers_minus_tweeters, columns = ["user_id"], dtype = 'int64')
    

# Count tweeters who had no followers
nofollower_files = sorted( os.listdir(path_tweeters_nofollowers) )
no_followers = np.array([], dtype = 'int64')
for file in nofollower_files:
    data = np.genfromtxt(path_tweeters_nofollowers + file, dtype = 'int64') #first row is header
    try:
        data = data[1:len(data)] #if empty, will return error
        no_followers = np.append(no_followers, data)
        no_followers = np.unique(no_followers)
        del data
    except:
        next
count_no_followers = no_followers.shape[0]
no_followers = pd.DataFrame(no_followers, columns = ["user_id"], dtype = 'int64')

# Save
all_followers.to_csv(outpath_follower_data + "unique_followers.csv", index = False)
followers_minus_tweeters.to_csv(outpath_follower_data + "unique_followers_excl_tweeters.csv", index = False)
no_followers.to_csv(outpath_follower_data + "tweeters_nofollowers.csv", index = False)
del all_followers, followers_minus_tweeters, no_followers

####################
# Count up friends of tweeters
####################
# Count up friends
friend_files = sorted( os.listdir(path_to_friends) )
num_files = len(friend_files)
all_friends = np.array([], dtype = 'int64')
for file in friend_files:
    data = np.genfromtxt(path_to_friends + file, skip_header = 1, dtype = 'int64') #first row is header
    all_friends = np.append(all_friends, data)
    all_friends = np.unique(all_friends)
    del data
all_friends_minus_tweeters = np.setdiff1d(all_friends, tweeters)
count_all_friends = all_friends.shape[0]
count_all_friends_minus_tweeters = all_friends_minus_tweeters.shape[0]
all_friends = pd.DataFrame(all_friends, columns = ["user_id"], dtype = 'int64')
all_friends_minus_tweeters = pd.DataFrame(all_friends_minus_tweeters, columns = ["user_id"], dtype = 'int64')
    

# Count tweeters who have no friends
tweeters_nofriends_files = sorted( os.listdir(path_to_tweeters_nofriends) )
no_friends = np.array([], dtype = 'int64')
for file in tweeters_nofriends_files:
    data = np.genfromtxt(path_to_tweeters_nofriends + file, dtype = 'int64') #first row is header
    try:
        data = data[1:len(data)] #if empty, will return error
        no_friends = np.append(no_friends, data)
        no_friends = np.unique(no_friends)
        del data
    except:
        next
count_no_friends = no_friends.shape[0]
no_friends = pd.DataFrame(no_friends, columns = ["user_id"], dtype = 'int64')

# Save
all_friends.to_csv(outpath_friend_data + "unique_friends.csv", index = False)
all_friends_minus_tweeters.to_csv(outpath_friend_data + "unique_friends_excl_tweeters.csv", index = False)
no_friends.to_csv(outpath_friend_data + "tweeters_nofriends.csv", index = False)
del all_friends, all_friends_minus_tweeters, no_friends

####################
# Summarize
####################
# Measure number of users and write out to file
unique_users = pd.DataFrame({'user_type': ["Tweeters", 
                                           "Followers", "Followers excluding tweeters", "Tweeters w/o followers",
                                           "Friends of tweeters", "Friends excludingtweeters", "Tweeters w/o friends"], 
                             'count': [count_tweeters, 
                                       count_followers, count_followers_minus_tweeters, count_no_followers,
                                       count_all_friends, count_all_friends_minus_tweeters, count_no_friends]})
unique_users.to_csv(outpath + "summary_unique_users.csv", index = False)
