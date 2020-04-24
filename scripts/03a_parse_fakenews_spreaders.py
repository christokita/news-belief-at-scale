#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 15:14:12 2020

@author: ChrisTokita

SCRIPT:
Determine who shared fake news articles and who was potentially exposed to it. 
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import re
import os
import sys

# Get file numbers
i = int(sys.argv[1]) # get which chunk of the 

# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/"

####################
# Parse articles to get false news articles
####################
# Get article IDs of false/misleading articles (as evaluated by fact checkers)
news_evaluations = pd.read_csv(data_directory + "data/articles/evaluations.csv")
fakenews_ids = news_evaluations["article_num"][news_evaluations['mode of FC'] == "FM"]

# Get URL set for FM articles
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")
fm_articles = articles[articles['total article number'].isin(fakenews_ids)]

# Load and parse tweets according to full and shortened URLs
tweets = pd.read_csv(data_directory + "data_derived/tweets/parsed_tweets.csv")
fm_tweets = pd.DataFrame(columns = tweets.columns)
for j in range(len(fm_articles)):
    link = fm_articles['link'].iloc[j]
    shortlink = fm_articles['short link'].iloc[j]
    try:
        has_full_link = tweets['urls_expanded'].str.contains(link)
        has_short_link = tweets['urls_expanded'].str.contains(shortlink)
    except:
        has_short_link = pd.Series(np.repeat(False, tweets.shape[0]))
    has_link = has_full_link | has_short_link #boolean operator to find which indices have one of the two possible links
    fm_tweets = fm_tweets.append(tweets[has_link])
    
    
####################
# Determine which users tweeted false articles
####################
# Get user IDs of tweeters of FM articles, 
fm_tweeters = fm_tweets['user_id']
fm_tweeters, fake_freq = np.unique(fm_tweeters, return_counts = True)

# Write FM tweets and tweeters to file
fm_tweets.to_csv(data_directory + "data_derived/tweets/FM_tweets.csv", index = False)

# plot frequency of fake news article sharing
import matplotlib.pyplot as plt  
plt.hist(fake_freq, log = True)


####################
# Determine how many users were potentially exposed to these articles
####################
# Take chunk of FM article tweeters and combine followers
fm_tweeters = fm_tweeters[ (260*i) : (260*(i+1)) ] #There are 25,546 users who tweeted a FM article

# loop through files and get set of unique followers who were exposed to the article
follower_files = os.listdir(data_directory + "data/followers/")
follower_files = [file for file in follower_files if re.match('[0-9]', file)] #filter out hidden copies of same files
exposed_followers = np.array([])
no_follower_list = np.array([])
for user_id in fm_tweeters:
    regex = re.compile(r"[0-9].*_%s.csv" % str(user_id))
    file = list(filter(regex.match, follower_files))
    try:
        if len(file) > 1:
            print("WARNING: user_id = %d matches multiple follower list files." % user_id)
        follower_list = pd.read_csv(data_directory + "data/followers/" + file[0])
        exposed_followers = np.append(exposed_followers, follower_list)
        exposed_followers = np.unique(exposed_followers)
    except:
        no_follower_list = np.append(no_follower_list, user_id)
exposed_followers = exposed_followers.astype(int)
no_follower_list = no_follower_list.astype(int)


####################
# Determine the friends of FM tweeters (i.e., who they are following)
####################
# loop through friend files and get set of unique followers who were exposed to the article
friend_files = os.listdir(data_directory + "data/friends/")
friend_files = [file for file in friend_files if re.match('[0-9]', file)] #filter out hidden copies of same files
fm_friends = np.array([])
no_friends_list = np.array([])
for user_id in fm_tweeters:
    regex = re.compile(r"[0-9].*_%s.csv" % str(user_id))
    file = list(filter(regex.match, friend_files))
    try:
        if len(file) > 1:
            print("WARNING: user_id = %d matches multiple follower list files." % user_id)
        friend_list = pd.read_csv(data_directory + "data/friends/" + file[0])
        fm_friends = np.append(fm_friends, friend_list)
        fm_friends = np.unique(fm_friends)
    except:
        no_friends_list = np.append(no_friends_list, user_id)
fm_friends = fm_friends.astype(int)
no_friends_list = no_friends_list.astype(int)



####################
# Save
####################
chunk_label = str(i).zfill(2)

exposed_followers = pd.DataFrame(exposed_followers, columns = ['user_id'])
no_follower_list = pd.DataFrame(no_follower_list, columns = ['user_id'])
fm_friends = pd.DataFrame(fm_friends, columns = ['user_id'])
no_friends_list = pd.DataFrame(no_friends_list, columns = ['user_id'])

exposed_followers.to_csv(data_directory + "data_derived/followers/processed_exposed_followers/followers_exposed_fakenews_" + chunk_label + ".csv", index = False)
no_follower_list.to_csv(data_directory + "data_derived/followers/nofollowers_fm_tweeters/nofollowers_fm_tweeters_" + chunk_label + ".csv", index = False)
fm_friends.to_csv(data_directory + "data_derived/friends/processed_fm_tweeter_friends/friends_fm_" + chunk_label + ".csv", index = False)
no_friends_list.to_csv(data_directory + "data_derived/friends/nofriends_fm_tweeters/nofriends_fm_tweeters_" + chunk_label + ".csv", index = False)