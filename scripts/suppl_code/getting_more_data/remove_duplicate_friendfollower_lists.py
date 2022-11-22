#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 18:23:11 2022

@author: ChrisTokita

SCRIPT
Remove duplicate follower/friend lists for users. 
We recently got new friend and follower lists for our new tweeters as well as users who we only had a friend or follower list but not both.
Thus, we have duplicates. We will go through our data, find cases where we have more than one friend/follower list for a user, and will delete the shorter one.
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import re
import os

# high level directory (external HD or cluster storage)
data_directory = '/Volumes/CKT-DATA/fake-news-diffusion/'


####################
# load data
####################
tweets = pd.read_csv(data_directory + 'data_derived/tweets/tweets_parsed.csv',
                     dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                              'user_id': object, 'tweet_id': object, 
                              'retweeted_user_id': object, 'retweet_id': object,
                              'quoted_user_id': object, 'quoted_id': object})
tweeters = tweets['user_id']
tweeters = np.unique(tweeters)

# Get friend and follower lists in our dataset
follower_files = os.listdir(data_directory + 'data/followers/')
friend_files = os.listdir(data_directory + 'data/friends/')


####################
# Function to determine longer friend/follower list and remove all others for a user
####################
def delete_shorter_duplicate_lists(file_list: list, path_to_files: str):
    """
    :param file_list:       friend/follower lists to check.
    :param path_to_files:   path to where these files are stored.
    :return:                none, but deletes shorter list(s) from file.
    """
    # Compare list lengths
    list_lengths = np.array([])
    for file in file_list:
        loaded_list = np.genfromtxt(path_to_files + file, dtype = str)
        n_items = loaded_list.size - 1 #first row is header
        list_lengths = np.append(list_lengths, n_items)

    # Find longest file
    longest_list = np.argmax(list_lengths)
    
    # Delete longest file from list of files
    del file_list[longest_list]
    
    # Delete file
    for file_to_remove in file_list:
        os.remove(path_to_files + file_to_remove)
        print('deleted: {}{}'.format(path_to_files,file_to_remove))



####################
# Go through and remove duplicate files
####################
for user_id in tweeters:
    
    # Find friend/follower lists
    regex = re.compile(r'^[0-9]{4}__[0-9]{2}__[0-9]{2}__%s\.csv' % user_id)
    tweeter_friend_lists = list(filter(regex.match, friend_files))
    tweeter_follower_lists = list(filter(regex.match, follower_files))
    
    # Determine if more than one friend file
    if len(tweeter_friend_lists) > 1:
        
        delete_shorter_duplicate_lists(file_list = tweeter_friend_lists,
                                       path_to_files = data_directory + 'data/friends/')
        
        
    # Determine if more than one follow file, and if so, delete shorter list(s)
    if len(tweeter_follower_lists) > 1:
        
        delete_shorter_duplicate_lists(file_list = tweeter_follower_lists,
                                       path_to_files = data_directory + 'data/followers/')
