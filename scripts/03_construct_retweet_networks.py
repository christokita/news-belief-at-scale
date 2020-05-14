#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 16:15:11 2020

@author: ChrisTokita

SCRIPT:
Construct the rewteet network of fake news articles
"""

####################
# Load libraries and packages, set important parameters, load data
####################
import pandas as pd
import numpy as np
import re
import os
import multiprocessing as mp
import datetime as dt

# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/"
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/"

# Load data
fm_tweets = pd.read_csv(data_directory + "data_derived/tweets/FM_tweets.csv")
all_tweets = pd.read_csv(data_directory + "data_derived/tweets/parsed_tweets.csv")

#fm_tweets_sorted = fm_tweets.sort_values(by = ['retweet_count'], ascending = False)
#top_tweet = fm_tweets_sorted['tweet_id'][0]
#top_tweeter = fm_tweets_sorted['user_id'][0]
#
#np.sum(fm_tweets['is_retweet'])
#np.sum(fm_tweets['is_quote_nominal'])
#np.sum(fm_tweets['is_quote'])


####################
# Construct network
####################
# Note i = 8 is example self-retweet
# Note i = 9,10 is example of non-direct retweetn of non-existent friend??
# Note i = 10

friend_files = [file for file in os.listdir(data_directory + "data/friends/") if re.match('^[0-9]', file)] #filter out hidden copies of same files
follower_files = [file for file in os.listdir(data_directory + "data/followers/") if re.match('^[0-9]', file)] #filter out hidden copies of same files


def determine_retweet_edges(i, data, friend_files, follower_files):

    # Grab specific tweet and preliinary information    
    tweet = data.iloc[i,:]
    user_id = tweet['user_id'].astype(int)
    rt_edge = pd.DataFrame({'Source': user_id, 'Target': None, 'Type': None}, index = [0])
    
    
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
    
    # Filter FM tweets to only those talking about this news article
    article_id = tweet['total_article_number']
    article_tweets = fm_tweets[fm_tweets['total_article_number'] == article_id]
    
    # Get friends of focal user and id of original tweeter who was RT'd
    regex_friend = re.compile(r"[0-9].*_%s.csv" % user_id)
    fr_file = list(filter(regex_friend.match, friend_files))
    friends = np.genfromtxt(data_directory + "data/friends/" + fr_file[0], skip_header = 1, dtype = int)
    rt_user_id = tweet['retweeted_user_id'].astype(int)
    
    # If retweeted user is followed by focal user, count that as flow of tweet
    if rt_user_id in friends:
        rt_edge['Target'] = rt_user_id
        rt_edge['Type'] = "Direct RT"
        
    # Check if self-retweet
    elif rt_user_id == user_id:
        rt_edge['Target'] = rt_user_id
        rt_edge['Type'] = "Self RT"
        
    # Otherwise determine if indirect RT or phantom RT
    else:
        
        try:
            # Grab FM tweeters who are both i's friends and tweeted the news article in question
            article_tweeters = np.unique(article_tweets['user_id'].astype(int))
            candidate_friends = np.intersect1d(friends, article_tweeters)
            rt_time = dt.datetime.strptime(tweet['tweet_time'], '%a %b %d %H:%M:%S %z %Y')
            candidate_tweets = fm_tweets[(fm_tweets['user_id'].isin(candidate_friends))].copy()
            
            # Filter candidate tweets to only those talking about the same story and that occured before i's retweet
            candidate_tweets.loc[:,'tweet_time'] = pd.to_datetime(candidate_tweets['tweet_time'], format = '%a %b %d %H:%M:%S %z %Y')
            candidate_tweets = candidate_tweets[candidate_tweets['tweet_time'] < rt_time] #tweets before focal individual
            retweeted = candidate_tweets.loc[ candidate_tweets['tweet_time'].idxmax() ] #most recent tweet before focal tweet
            
            # Set indirect RT edge
            rt_edge['Target'] = retweeted['user_id'].astype(int)
            rt_edge['Type'] = "Direct RT"
            
        except:
            # set phantom RT edge
            rt_edge['Target'] = rt_user_id
            rt_edge['Type'] = "Phantom RT"
    
    return rt_edge


def parse_quotedtweet(tweet, user_id, rt_edge, friend_files, follower_files, fm_tweets):
    pass