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
                              'user_id': object, 'tweet_id': object, 
                              'retweeted_user_id': object, 'retweet_id': object,
                              'quoted_user_id': object, 'quoted_id': object})
tweeters = tweets['user_id']
    
# Get user IDs of tweeters
tweeters, tweet_freq = np.unique(tweeters, return_counts = True)
tweeters = pd.DataFrame(tweeters, columns = ['user_id'])
count_tweeters = tweeters.shape[0]

# Save tweeters
tweeters.to_csv(outpath_tweet_data + "unique_tweeters.csv", index = False)


####################
# Count up unique followers of tweeters
####################
# Count up followers
follower_files = sorted( os.listdir(path_to_followers) )
num_files = len(follower_files)
all_followers = np.array([], dtype = object)
for file in follower_files:
    data = np.genfromtxt(path_to_followers + file, skip_header = 1, dtype = str) #first row is header
    all_followers = np.append(all_followers, data)
    all_followers = np.unique(all_followers)
    del data
followers_minus_tweeters = np.setdiff1d(all_followers, tweeters)
count_followers = all_followers.shape[0]
count_followers_minus_tweeters = followers_minus_tweeters.shape[0]
all_followers = pd.DataFrame(all_followers, columns = ["user_id"], dtype = object)
followers_minus_tweeters = pd.DataFrame(followers_minus_tweeters, columns = ["user_id"], dtype = object)
    

# Count tweeters who had no followers
nofollower_files = sorted( os.listdir(path_tweeters_nofollowers) )
no_followers = np.array([], dtype = object)
for file in nofollower_files:
    data = np.genfromtxt(path_tweeters_nofollowers + file, dtype = str) #first row is header
    try:
        data = data[1:len(data)] #if empty, will return error
        no_followers = np.append(no_followers, data)
        no_followers = np.unique(no_followers)
        del data
    except:
        next
count_no_followers = no_followers.shape[0]
no_followers = pd.DataFrame(no_followers, columns = ["user_id"], dtype = object)

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
all_friends = np.array([], dtype = object)
for file in friend_files:
    data = np.genfromtxt(path_to_friends + file, skip_header = 1, dtype = str) #first row is header
    all_friends = np.append(all_friends, data)
    all_friends = np.unique(all_friends)
    del data
all_friends_minus_tweeters = np.setdiff1d(all_friends, tweeters)
count_all_friends = all_friends.shape[0]
count_all_friends_minus_tweeters = all_friends_minus_tweeters.shape[0]
all_friends = pd.DataFrame(all_friends, columns = ["user_id"], dtype = object)
all_friends_minus_tweeters = pd.DataFrame(all_friends_minus_tweeters, columns = ["user_id"], dtype = object)
    

# Count tweeters who have no friends
tweeters_nofriends_files = sorted( os.listdir(path_to_tweeters_nofriends) )
no_friends = np.array([], dtype = object)
for file in tweeters_nofriends_files:
    data = np.genfromtxt(path_to_tweeters_nofriends + file, dtype = str) #first row is header
    try:
        data = data[1:len(data)] #if empty, will return error
        no_friends = np.append(no_friends, data)
        no_friends = np.unique(no_friends)
        del data
    except:
        next
count_no_friends = no_friends.shape[0]
no_friends = pd.DataFrame(no_friends, columns = ["user_id"], dtype = object)

# Save
all_friends.to_csv(outpath_friend_data + "unique_friends.csv", index = False)
all_friends_minus_tweeters.to_csv(outpath_friend_data + "unique_followers_excl_tweeters.csv", index = False)
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
unique_users.to_csv(outpath + "unique_user_summary.csv", index = False)
