#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 15:14:12 2020

@author: ChrisTokita

SCRIPT:
Determine who shared news articles and who their firnds and followers are 

Note: 
The overall logic of this script is to 
(1) load in user IDs as int64 to prevent truncation. 
(2) convert to str?
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import re
import os
import sys
import math

# Get file numbers
i = int(sys.argv[1]) # get which chunk of the 

# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/"

####################
# load data
####################
tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_parsed.csv",
                     dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                              'user_id': object, 'tweet_id': object, 
                              'retweeted_user_id': object, 'retweet_id': object,
                              'quoted_user_id': object, 'quoted_id': object})
tweeters = tweets['user_id']

####################
# Determine how many users were potentially exposed to these articles
####################
# Take chunk of FM article tweeters and combine followers
step_size = math.ceil(len(tweeters) / 300)
tweeters = tweeters[ (step_size*i) : (step_size*(i+1)) ]
tweeters = np.unique(tweeters) 

# loop through files and get set of unique followers who were exposed to the article
follower_files = os.listdir(data_directory + "data/followers/")
follower_files = [file for file in follower_files if re.match('^[0-9]', file)] #filter out hidden copies of same files
followers = np.array([], dtype = object)
no_follower_list = np.array([], dtype = object)
for user_id in tweeters:
    regex = re.compile(r"^[0-9]{4}__[0-9]{2}__[0-9]{2}__%s\.csv" % user_id)
    file = list(filter(regex.match, follower_files))
    if len(file) > 1:
        print("WARNING: user_id = %s matches multiple follower list files." % user_id)
    try:
        follower_list = np.genfromtxt(data_directory + "data/followers/" + file[0], dtype = str)
        follower_list = follower_list[1:len(follower_list)] #remove header, will raise error if empty
        followers = np.append(followers, follower_list)
        followers = np.unique(followers)
    except:
        no_follower_list = np.append(no_follower_list, user_id)
        
# Save
chunk_label = str(i).zfill(2)
followers = pd.DataFrame(followers, columns = ['user_id'], dtype = object)
no_follower_list = pd.DataFrame(no_follower_list, columns = ['user_id'], dtype = object)
followers.to_csv(data_directory + "data_derived/followers/processed_followers/followers_" + chunk_label + ".csv", index = False)
no_follower_list.to_csv(data_directory + "data_derived/followers/tweeters_nofollowers/tweeters_nofollowers_" + chunk_label + ".csv", index = False)
del followers, no_follower_list


####################
# Determine the friends of tweeters (i.e., who they are following)
####################
# loop through friend files and get set of unique followers who were exposed to the article
friend_files = os.listdir(data_directory + "data/friends/")
friend_files = [file for file in friend_files if re.match('^[0-9]', file)] #filter out hidden copies of same files
friends = np.array([], dtype = object)
no_friends_list = np.array([], dtype = object)
for user_id in tweeters:
    regex = re.compile(r"[0-9]{4}__[0-9]{2}__[0-9]{2}__%s\.csv" % user_id)
    file = list(filter(regex.match, friend_files))
    if len(file) > 1:
        print("WARNING: user_id = %s matches multiple friend list files." % user_id)
    try:
        friend_list = np.genfromtxt(data_directory + "data/friends/" + file[0], dtype = str)
        friend_list = friend_list[1:len(friend_list)] #remove header, will raise error if empty
        friends = np.append(friends, friend_list)
        friends = np.unique(friends)
    except:
        no_friends_list = np.append(no_friends_list, user_id)
        
# Save
friends = pd.DataFrame(friends, columns = ['user_id'], dtype = object)
no_friends_list = pd.DataFrame(no_friends_list, columns = ['user_id'], dtype = object)
friends.to_csv(data_directory + "data_derived/friends/processed_friends/friends_" + chunk_label + ".csv", index = False)
no_friends_list.to_csv(data_directory + "data_derived/friends/tweeters_nofriends/tweeters_nofriends_" + chunk_label + ".csv", index = False)
del friends, no_friends_list