#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 15:14:12 2020

@author: ChrisTokita

SCRIPT:
Determine who shared fake news articles and who was potentially exposed to it. 

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
# Parse articles to assign article number to tweet and filter to fake news tweets
####################
# Get article IDs of false/misleading articles (as evaluated by fact checkers)
news_evaluations = pd.read_csv(data_directory + "data/articles/evaluations.csv")
fakenews_ids = news_evaluations["article_num"][news_evaluations['mode of FC'] == "FM"]

# Load articles dataset to get their URLs
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")

# Function to clean up links for better matching
def simplify_link(link):
    if not pd.isnull(link):
        link = re.sub('http.*//', '', link)
        link = re.sub('^www\.', '', link)
        link = re.sub('\?.*$', '', link)
        link = re.sub('/$', '', link)
    return link

# Load and parse tweets according to full and shortened URLs, assigning article number ID
tweets = pd.read_csv(data_directory + "data_derived/tweets/parsed_tweets.csv")
tweets = tweets.drop(['total_article_number'], axis=1, errors='ignore') #drop article ID column if it had previously been assigned
tweets = tweets.join(pd.DataFrame(np.repeat(np.nan, tweets.shape[0]), columns = ['total_article_number']))
for j in range(articles.shape[0]):
    # Prep links for pattern matching
    link = simplify_link( articles['link'].iloc[j] )
    shortlink = simplify_link( articles['short link'].iloc[j] )
    # Search through URLS
    has_full_link = tweets['urls_expanded'].str.contains(link) | tweets['quoted_urls_expanded'].str.contains(link)
    if not pd.isna(shortlink):
        has_short_link_main = tweets['urls_expanded'].str.contains(shortlink) | tweets['urls'].str.contains(shortlink)
        has_short_link_quoted = tweets['quoted_urls_expanded'].str.contains(shortlink) | tweets['quoted_urls'].str.contains(shortlink)
        has_short_link = has_short_link_main | has_short_link_quoted
    elif pd.isna(shortlink):
        has_short_link = pd.Series(np.repeat(False, tweets.shape[0])) #if shortlink is nan
    has_link = has_full_link | has_short_link #boolean operator to find which indices have one of the two possible links
    # Assign article number ID
    tweets.loc[has_link, 'total_article_number'] = articles['total article number'].iloc[j]
    
# Write tweets and filtered FM tweets and tweeters to file
tweets.to_csv(data_directory + "data_derived/tweets/parsed_tweets.csv", index = False)
fm_tweets = tweets[tweets['total_article_number'].isin(fakenews_ids)]
fm_tweets.to_csv(data_directory + "data_derived/tweets/FM_tweets.csv", index = False)

# ID and save tweets without IDs
missing_ids = tweets[pd.isnull(tweets['total_article_number'])]
if missing_ids.shape[0] > 0:
    missing_article_ids = [x  for x in  np.arange(1, 166) if x not in np.unique(tweets['total_article_number'])]
    missing_articles = articles[articles['total article number'].isin(missing_article_ids)]
    missing_articles.to_csv(data_directory + "data_derived/articles/articles_notfoundintweets.csv", index = False)
    missing_ids.to_csv(data_directory + "data_derived/tweets/noarticleID_tweets.csv", index = False)
     
# Get user IDs of tweeters of FM articles, 
fm_tweeters = fm_tweets['user_id'].astype(int)
fm_tweeters = np.unique(fm_tweeters)

####################
# Determine how many users were potentially exposed to these articles
####################
# Take chunk of FM article tweeters and combine followers
step_size = math.ceil(len(fm_tweeters) / 100)
fm_tweeters = fm_tweeters[ (step_size*i) : (step_size*(i+1)) ] #There are 26,493 users who tweeted a FM article

# loop through files and get set of unique followers who were exposed to the article
follower_files = os.listdir(data_directory + "data/followers/")
follower_files = [file for file in follower_files if re.match('^[0-9]', file)] #filter out hidden copies of same files
exposed_followers = np.array([], dtype = int)
no_follower_list = np.array([], dtype = int)
for user_id in fm_tweeters:
    regex = re.compile(r"[0-9].*_%s.csv" % user_id)
    file = list(filter(regex.match, follower_files))
    if len(file) > 1:
        print("WARNING: user_id = %d matches multiple follower list files." % user_id)
    try:
        follower_list = np.genfromtxt(data_directory + "data/followers/" + file[0], dtype = int)
        follower_list = follower_list[1:len(follower_list)] #remove header, will raise error if empty
        exposed_followers = np.append(exposed_followers, follower_list)
        exposed_followers = np.unique(exposed_followers)
    except:
        no_follower_list = np.append(no_follower_list, user_id)


####################
# Determine the friends of FM tweeters (i.e., who they are following)
####################
# loop through friend files and get set of unique followers who were exposed to the article
friend_files = os.listdir(data_directory + "data/friends/")
friend_files = [file for file in friend_files if re.match('^[0-9]', file)] #filter out hidden copies of same files
fm_friends = np.array([], dtype = int)
no_friends_list = np.array([], dtype = int)
for user_id in fm_tweeters:
    regex = re.compile(r"[0-9].*_%s.csv" % user_id)
    file = list(filter(regex.match, friend_files))
    if len(file) > 1:
        print("WARNING: user_id = %d matches multiple follower list files." % user_id)
    try:
        friend_list = np.genfromtxt(data_directory + "data/friends/" + file[0], dtype = int)
        friend_list = friend_list[1:len(friend_list)] #remove header, will raise error if empty
        fm_friends = np.append(fm_friends, friend_list)
        fm_friends = np.unique(fm_friends)
    except:
        no_friends_list = np.append(no_friends_list, user_id)


####################
# Save
####################
chunk_label = str(i).zfill(2)

exposed_followers = pd.DataFrame(exposed_followers, columns = ['user_id'], dtype = int)
no_follower_list = pd.DataFrame(no_follower_list, columns = ['user_id'], dtype = int)
fm_friends = pd.DataFrame(fm_friends, columns = ['user_id'], dtype = int)
no_friends_list = pd.DataFrame(no_friends_list, columns = ['user_id'], dtype = int)

exposed_followers.to_csv(data_directory + "data_derived/followers/processed_exposed_followers/followers_exposed_fakenews_" + chunk_label + ".csv", index = False)
no_follower_list.to_csv(data_directory + "data_derived/followers/nofollowers_fm_tweeters/nofollowers_fm_tweeters_" + chunk_label + ".csv", index = False)
fm_friends.to_csv(data_directory + "data_derived/friends/processed_fm_tweeter_friends/friends_fm_" + chunk_label + ".csv", index = False)
no_friends_list.to_csv(data_directory + "data_derived/friends/nofriends_fm_tweeters/nofriends_fm_tweeters_" + chunk_label + ".csv", index = False)
