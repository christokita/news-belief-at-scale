#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  8 22:13:14 2022

@author: ChrisTokita

SCRIPT:
Prep new tweets, tweeters, and tweeter friend/follower lists from expanded scope of search for article tweets

We did an extensive search to make sure we got all the missing tweets of our tracked news articles, as well as the friend/follower lists from new tweeters in this extended search.
The collection of these data came in a format that is different from the original data, so here we will format them so that they can be combined with our older data in our analysis pipeline.

"""

####################
# Load packages and data, set important parameters
####################
import os
import json
import pandas as pd
import numpy as np

# Paths
data_directory = '/Volumes/CKT-DATA/fake-news-diffusion/' #external HD


####################
# Combine individual json tweet files into one json file
####################
# Get list of individual tweets
tweet_list = os.listdir(data_directory + 'data/tweets/expanded_window_tweets/')
tweet_list = [x for x in tweet_list if not x.startswith('.')]

# Loop through and collect tweets into one json file
expanded_window_tweets = []
for tweet_json in tweet_list:
    
    with open(data_directory + 'data/tweets/expanded_window_tweets/' + tweet_json) as f:
        tweet_json_str = json.load(f) #read in json string
        tweet = json.loads(tweet_json_str) #convert json string to proper json dictionary
        expanded_window_tweets.append(tweet)
    

# Write to file
with open(data_directory + 'data/tweets/expanded_window_tweets.json', 'w') as out_file:
    
    for json_obj in expanded_window_tweets:
        json.dump(json_obj,  out_file)
        out_file.write(',\n')
        
# Free up memory
del expanded_window_tweets
        
        
####################
# Break up combined follower lists into indiviudal follower lists like the rest of the data
####################
# Read in follower list(s) for new tweeters
new_tweeter_followers = pd.read_csv(data_directory + 'data/followers/_new_tweeters_followers_jan2022/FULL_New_Tweeters_Followers_IDs.csv', dtype = object)

# Loop through each unique tweeter, get their individual list of followers, and write to file
unique_accounts_for_followers = np.unique(new_tweeter_followers['account_id'])
for user in unique_accounts_for_followers:
    followers = new_tweeter_followers[new_tweeter_followers.account_id == user]
    followers = followers[['follower_id']]
    followers = followers.rename(columns = {'follower_id': 'user_id_followers'}) #this is how our original data is formatted
    followers.to_csv(data_directory + 'data/followers/2022__01__09__' + user + '.csv', index = False) #standard file format 
    
    
####################
# Break up combined friend lists into indiviudal friend lists like the rest of the data
####################
# Read in friend list(s) for new tweeters
new_tweeter_friends = pd.read_csv(data_directory + 'data/friends/_new_tweeters_friends_jan2022/FULL_New_Tweeters_Friends_IDs.csv', dtype = object)

# Loop through each unique tweeter, get their individual list of friends, and write to file
unique_accounts_for_friends = np.unique(new_tweeter_friends['account_id'])
for user in unique_accounts_for_friends:
    friends = new_tweeter_friends[new_tweeter_friends.account_id == user]
    friends = friends[['friend_id']]
    friends = friends.rename(columns = {'friend_id': 'user_id_friends'}) #this is how our original data is formatted
    friends.to_csv(data_directory + 'data/friends/2022__01__09__' + user + '.csv', index = False) #standard file format 