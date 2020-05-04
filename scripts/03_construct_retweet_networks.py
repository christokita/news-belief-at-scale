#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 16:15:11 2020

@author: ChrisTokita

SCRIPT:
Construct the rewteet network of fake news articles
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import re
import os
import multiprocessing as mp

# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/"
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/"


####################
# Load fake news tweets
####################
fm_tweets = pd.read_csv(data_directory + "data_derived/tweets/FM_tweets.csv")
fm_tweets_sorted = fm_tweets.sort_values(by = ['retweet_count'], ascending = False)
top_tweet = fm_tweets_sorted['tweet_id'][0]
top_tweeter = fm_tweets_sorted['user_id'][0]

all_tweets = pd.read_csv(data_directory + "data_derived/tweets/parsed_tweets.csv")

np.sum(fm_tweets['is_retweet'])
np.sum(fm_tweets['is_quote_nominal'])
np.sum(fm_tweets['is_quote'])


####################
# Construct network
####################
friend_files = os.listdir(data_directory + "data/friends/")
friend_files = [file for file in friend_files if re.match('[0-9]', file)] #filter out hidden copies of same files
follower_files = os.listdir(data_directory + "data/followers/")
follower_files = [file for file in follower_files if re.match('[0-9]', file)] #filter out hidden copies of same files


def determine_retweet_edges(i, data, friend_files, follower_files):

    # Grab specific tweet and preliinary information    
    tweet = data.iloc[i,:]
    user_id = tweet['user_id']
    rt_edge = pd.DataFrame({'from': user_id, 'to': None}, dtype = int, index = [0])
    
    
    # If retweet determine who they were actually retweeting
    if tweet['is_retweet']:
        rt_edge = parse_retweet(tweet, user_id, rt_edge, friend_files, follower_files)
    elif tweet['is_quote'] or tweet['is_quote_nominal']:
        pass
    
    
def parse_retweet(tweet, user_id, rt_edge, friend_files, follower_files, fm_tweets):
    # Function that will take a retweet and determine who the RT came from: original tweeter or friend that tweeted
    #
    # OUTPUT
    # - rt_edge: returns the rt_edge dataframe row with appropriate user_id filled in the 'to' column
    
    # Get friends of focal user
    regex_friend = re.compile(r"[0-9].*_%s.csv" % str(user_id))
    fr_file = list(filter(regex_friend.match, friend_files))
    friends = np.loadtxt(data_directory + "data/friends/" + fr_file[0], skiprows = 1, dtype = int)
    
    # Get followers of retweeted user
    rt_user_id = int(tweet['retweeted_user_id'])
    regex_rt_followers = re.compile(r"[0-9].*_%s.csv" % str(rt_user_id))
    fol_file = list(filter(regex_rt_followers.match, follower_files))
    rt_user_followers = np.loadtxt(data_directory + "data/followers/" + fol_file[0], skiprows = 1, dtype = int)
    
    # If retweeted user is followed by focal user, count that as flow of tweet
    if rt_user_id in friends:
        rt_edge['to'] = rt_user_id
        
    # Otherwise determine which friends follow the retweeted user
    else:
        friends_following_rt_user = np.intersect1d(friends, rt_user_followers)
        fm_tweeters = fm_tweets['user_id']
        candidate_friends = np.intersect1d(friends_following_rt_user, fm_tweeters)
        # Left off here!