#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 16:55:00 2021

@author: ChrisTokita

SCRIPT:
    
Oct 19, 2020
Check what new users were found when we searched for tweets sharing URLs of articles that otherwise appeared to not have any associated tweets.
Note: This part of the script was actually done on Jan 16, 2022 because I realized I hadn't done this data check yet.
    
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
import re

# Paths
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD



######################################## October 19, 2020 ########################################

####################
# Determine which new tweeters need to be searched for friends/followers
####################
# Read in tweets and get list of new tweeters
missing_tweeters = np.array([])
missing_tweets_file = data_directory + 'data/tweets/crowdsource_factchecking_missing_article_tweets.json'
with open(missing_tweets_file, 'r') as f:
    for tweet in f:
        missing_tweet = json.loads(tweet)
        missing_tweeters = np.append(missing_tweeters, missing_tweet['user']['id_str'])

missing_tweeters = np.unique(missing_tweeters) #drops 96 duplicate user IDs

# Compile list of users for which we (nominally) have their list of friend/followers
follower_lists = os.listdir(data_directory + 'data/followers/')
follower_lists = [re.sub('[0-9]{4}__[0-9]{2}__[0-9]{2}__', '', x) for x in follower_lists]
follower_lists = [re.sub('\.csv', '', x) for x in follower_lists]
follower_lists = np.array(follower_lists)
follower_lists = np.sort(follower_lists)[1:] #drop .DS_store

friend_lists = os.listdir(data_directory + 'data/friends/')
friend_lists = [re.sub('[0-9]{4}__[0-9]{2}__[0-9]{2}__', '', x) for x in friend_lists]
friend_lists = [re.sub('\.csv', '', x) for x in friend_lists]
friend_lists = np.array(friend_lists)
friend_lists = np.sort(friend_lists)[1:] #drop .DS_store

have_lists = np.intersect1d(follower_lists, friend_lists) #only count those we have both friend and follower list for


# Determine which new tweeters aren't in our set that is already accounted for
missing_tweeters_need_to_check = np.setdiff1d(missing_tweeters, have_lists)
missing_tweeters_need_to_check = pd.DataFrame(data = missing_tweeters_need_to_check, 
                                          columns = ['user_id'],
                                          dtype = object)
missing_tweeters_need_to_check['user_id_str'] = "\"" + missing_tweeters_need_to_check['user_id'] + "\""


####################
# Write to file
####################
datestamp = '2020-10-19'
missing_tweeters_need_to_check.to_csv(data_directory + 'data_derived/_data_checks/missing_tweeters_needing_friendsfollowers_{}.csv'.format(datestamp), index = False)




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
# Read in tweets and get list of new tweeters
new_tweeters = np.array([])
new_tweets_file = data_directory + 'data/tweets/expanded_window_tweets.json'
with open(new_tweets_file, 'r') as f:
    for tweet in f:
        new_tweet = json.loads(tweet)
        new_tweeters = np.append(new_tweeters, new_tweet['user']['id_str'])

new_tweeters = np.unique(new_tweeters) #drops 2,287 duplicate user IDs

# Compile list of users for which we (nominally) have their list of friend/followers
follower_lists = os.listdir(data_directory + 'data/followers/')
follower_lists = [re.sub('[0-9]{4}__[0-9]{2}__[0-9]{2}__', '', x) for x in follower_lists]
follower_lists = [re.sub('\.csv', '', x) for x in follower_lists]
follower_lists = np.array(follower_lists)
follower_lists = np.sort(follower_lists)[1:] #drop .DS_store

friend_lists = os.listdir(data_directory + 'data/friends/')
friend_lists = [re.sub('[0-9]{4}__[0-9]{2}__[0-9]{2}__', '', x) for x in friend_lists]
friend_lists = [re.sub('\.csv', '', x) for x in friend_lists]
friend_lists = np.array(friend_lists)
friend_lists = np.sort(friend_lists)[1:] #drop .DS_store

have_lists = np.intersect1d(follower_lists, friend_lists) #only count those we have both friend and follower list for
missing_friend_list = np.setdiff1d(follower_lists, have_lists) #if in follower_lists but not have_list (i.e., both lists present), then missing friend list by logic
missing_follower_list = np.setdiff1d(friend_lists, have_lists)  #if in friend_lists but not have_list (i.e., both lists present), then missing follower list by logic


# Determine which new tweeters aren't in our set that is already accounted for
new_tweeters_need_to_check = np.setdiff1d(new_tweeters, have_lists)
new_tweeters_need_to_check = pd.DataFrame(data = new_tweeters_need_to_check, 
                                          columns = ['user_id'],
                                          dtype = object)
new_tweeters_need_to_check['user_id_str'] = "\"" + new_tweeters_need_to_check['user_id'] + "\""



####################
# Re-do from above (Sept 13, 2021): Compile list of users with no followers and no friends
####################
no_followers = pd.DataFrame(columns = ['user_id'],  dtype = object)
no_friends = pd.DataFrame(columns = ['user_id'],  dtype = object)

# Compile our list of no-friend/no-follower users that we flag during data processing
for file in os.listdir(data_directory + 'data_derived/followers/tweeters_nofollowers/'):
    user_list = pd.read_csv(data_directory + 'data_derived/followers/tweeters_nofollowers/' + file, dtype = object)
    no_followers = no_followers.append(user_list)
    del user_list
no_followers['user_id_str'] = "\"" + no_followers['user_id'] + "\""
no_followers = no_followers.drop_duplicates()
    
for file in os.listdir(data_directory + 'data_derived/friends/tweeters_nofriends/'):
    user_list = pd.read_csv(data_directory + 'data_derived/friends/tweeters_nofriends/' + file, dtype = object)
    no_friends = no_friends.append(user_list)
    del user_list
no_friends['user_id_str'] = "\"" + no_friends['user_id'] + "\""
no_friends = no_friends.drop_duplicates()


# Filter out users we already ran
already_checked_no_followers = pd.read_csv(data_directory + 'data_derived/_data_checks/users_no_followers_2021-9-13.csv')
already_checked_no_friends = pd.read_csv(data_directory + 'data_derived/_data_checks/users_no_friends_2021-9-13.csv')

no_followers_filtered = no_followers[~no_followers['user_id_str'].isin(already_checked_no_followers['user_id_str'])]
no_friends_filtered = no_friends[~no_friends['user_id_str'].isin(already_checked_no_friends['user_id_str'])]


# Add in users who appear to be missing one list
still_missing_follower_list = np.setdiff1d(missing_follower_list, no_followers_filtered['user_id']) #make sure they aren't in our no_followers_filtered list
still_missing_follower_list = pd.DataFrame(still_missing_follower_list, columns = ['user_id'])
still_missing_follower_list['user_id_str'] = "\"" + still_missing_follower_list['user_id'] + "\""

still_missing_friend_list = np.setdiff1d(missing_friend_list, no_friends_filtered['user_id']) #make sure they aren't already in our no_friends_filtered list
still_missing_friend_list = pd.DataFrame(still_missing_friend_list, columns = ['user_id'])
still_missing_friend_list['user_id_str'] = "\"" + still_missing_friend_list['user_id'] + "\""


# Combine
no_followers_complete = no_followers_filtered.append(still_missing_follower_list, ignore_index = True)
no_followers_complete = no_followers_complete.drop_duplicates()

no_friends_complete = no_friends_filtered.append(still_missing_friend_list, ignore_index = True)
no_friends_complete = no_friends_complete.drop_duplicates()


####################
# Write to file
####################
today = pd.to_datetime("today")
today = str(today.year) + "-" + str(today.month) + "-" + str(today.day)
no_followers_complete.to_csv(data_directory + 'data_derived/_data_checks/allusers_no_followers_{}.csv'.format(today), index = False)
no_friends_complete.to_csv(data_directory + 'data_derived/_data_checks/allusers_no_friends_{}.csv'.format(today), index = False)
new_tweeters_need_to_check.to_csv(data_directory + 'data_derived/_data_checks/new_tweeters_needing_friendsfollowers_{}.csv'.format(today), index = False)




######################################## April 23, 2022 ########################################

####################
# Get last tweet from each story
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
tweet_times = tweet_times.sort_values(['time_window_in_data'])

####################
# Write to file
####################
today = pd.to_datetime("today")
today = str(today.year) + "-" + str(today.month) + "-" + str(today.day)
tweet_times.to_csv(data_directory + 'data_derived/_data_checks/tweets_firstlast_{}.csv'.format(today), index = False)



