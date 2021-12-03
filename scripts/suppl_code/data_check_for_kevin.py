#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 16:55:00 2021

@author: ChrisTokita

SCRIPT:
    
Sept 13, 2021
Check two parts of the data for Kevin Aslett:
    (1) Confirm users who appear not to have friends and/or followers
    (2) Grab the first and last tweets of each story
    
Nov 28, 2021
Check two parts again. We now have more tweets for our main tracked articles after an expanded search. I also found a mistake in the gathering of no-friend-no-follower tweeters on Sept 13
    (1) Determine which new tweeters in our dataset need to be searched for friends/followers
    (2) Determine which users still need to be searched for friend/followers from our previous set of tweeters
"""


####################
# Load libraries and packages, set important parameter
####################
import pandas as pd
import numpy as np
import os
import json

# Paths
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD




######################################## September 13, 2021 ########################################

####################
# Compile list of users with no followers and no friends
####################
# NOTE: I made a mistake here and I should've parsed data_derived/friends/tweeters_nofriends and data_derived/followers/tweeters/no_followers
no_followers = pd.DataFrame(columns = ['user_id'],  dtype = object)
no_friends = pd.DataFrame(columns = ['user_id'],  dtype = object)

for file in os.listdir(data_directory + 'data_derived/followers/nofollowers_fm_tweeters/'):
    user_list = pd.read_csv(data_directory + 'data_derived/followers/nofollowers_fm_tweeters/' + file, dtype = object)
    no_followers = no_followers.append(user_list)
    del user_list
no_followers['user_id_str'] = "\"" + no_followers['user_id'] + "\""
    
for file in os.listdir(data_directory + 'data_derived/friends/nofriends_fm_tweeters/'):
    user_list = pd.read_csv(data_directory + 'data_derived/friends/nofriends_fm_tweeters/' + file, dtype = object)
    no_friends = no_friends.append(user_list)
    del user_list
no_friends['user_id_str'] = "\"" + no_friends['user_id'] + "\""

    
####################
# Get first and last tweet from each story
####################
labeled_tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                             dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                                      'user_id': 'int64', 'tweet_id': 'int64'})
labeled_tweets['tweet_time'] = pd.to_datetime(labeled_tweets['tweet_time'], format = '%a %b %d %H:%M:%S %z %Y')


minmax_tweets = labeled_tweets.groupby('total_article_number').agg({'tweet_time': [np.min, np.max, 'count']})

tweet_times = pd.DataFrame({'total_article_number': minmax_tweets.index,
                            'n_tweets_in_data': minmax_tweets['tweet_time']['count'],
                            'first_tweet_time': minmax_tweets['tweet_time']['amin'], 
                            'last_tweet_time': minmax_tweets['tweet_time']['amax']})
tweet_times['time_window_in_data'] = tweet_times['last_tweet_time'] - tweet_times['first_tweet_time']


####################
# Write to file
####################
today = pd.to_datetime("today")
today = str(today.year) + "-" + str(today.month) + "-" + str(today.day)
no_followers.to_csv(data_directory + 'data_derived/_data_checks/users_no_followers_{}.csv'.format(today), index = False)
no_friends.to_csv(data_directory + 'data_derived/_data_checks/users_no_friends_{}.csv'.format(today), index = False)
tweet_times.to_csv(data_directory + 'data_derived/_data_checks/tweets_firstlast_{}.csv'.format(today), index = False)




######################################## November 28, 2021 ########################################

####################
# Determine which new tweeters need to be searched for friends/followers
####################
new_tweets_file = data_directory + 'data/tweets/expanded_window_tweets.json'
f = open(new_tweets_file)
new_tweets = json.loads(f)
