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
import multiprocessing as mp

# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/"
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/"



############################## Parse ##############################

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
fm_tweeters = np.unique(fm_tweeters)

# Write FM tweets and tweeters to file
fm_tweets.to_csv(data_directory + "data_derived/tweets/FM_tweets.csv", index = False)


####################
# Determine which followers were potentially exposed to fake news
####################
# Function to parse a portion of exposed followers
def parse_exposed_followers(i, ids, directory):
    
    # Take chunk of FM article tweeters and combine followers. There are 25,546 users who tweeted a FM article
    ids = ids[ (260*i) : (260*(i+1)) ]

    # loop through files and get set of unique followers who were exposed to the article
    follower_files = os.listdir(directory + "data/followers/")
    follower_files = [file for file in follower_files if re.match('[0-9]', file)] #filter out hidden copies of same files
    exposed_followers = pd.DataFrame(columns = ['user_id'], dtype = int)
    no_followers = pd.DataFrame(columns = ['user_id'], dtype = int)
    for user_id in ids:
        regex = re.compile(r"[0-9].*_%s.csv" % str(user_id))
        file = list(filter(regex.match, follower_files))
        try:
            if len(file) > 1:
                print("WARNING: user_id = %d matches multiple follower list files." % user_id)
            follower_list = pd.read_csv(directory + "data/followers/" + file[0], names = ['user_id'], header = 0)
            exposed_followers = exposed_followers.append(follower_list, ignore_index = True)
            del follower_list
            exposed_followers = exposed_followers.drop_duplicates()
        except:
            no_followers = no_followers.append(user_id, ignore_index = True)
            
    # Write to file        
    exposed_followers.to_csv(directory + "data_derived/followers/processed_exposed_followers/followers_exposed_fakenews_" + i + ".csv", index = False)
    no_followers.to_csv(directory + "data_derived/followers/nofollowers_fm_tweeters/nofollowers_fm_tweeters_" + i + ".csv", index = False)

# Parse exposed followers
pool = mp.Pool(mp.cpu_count())
for i in range(100):
    pool.apply(parse_exposed_followers, args = (i, fm_tweeters, data_directory))
pool.close()


####################
# Determine the friends of FM tweeters (i.e., who they are following)
####################
# Function to parse relevant friends
def parse_fm_friends(i, ids, directory):
    
    # Take chunk of FM article tweeters and combine followers. There are 25,546 users who tweeted a FM article
    ids = ids[ (260*i) : (260*(i+1)) ]
    
    # loop through friend files and get set of unique followers who were exposed to the article
    friend_files = os.listdir(directory + "data/friends/")
    friend_files = [file for file in friend_files if re.match('[0-9]', file)] #filter out hidden copies of same files
    fm_friends = pd.DataFrame(columns = ['user_id'], dtype = int)
    no_friends_list = pd.DataFrame(columns = ['user_id'], dtype = int)
    for user_id in ids:
        regex = re.compile(r"[0-9].*_%s.csv" % str(user_id))
        file = list(filter(regex.match, friend_files))
        try:
            if len(file) > 1:
                print("WARNING: user_id = %d matches multiple follower list files." % user_id)
            friend_list = pd.read_csv(directory + "data/friends/" + file[0])
            fm_friends = fm_friends.append(friend_list, ignore_index = True)
            del friend_list
            fm_friends = fm_friends.drop_duplicates()
        except:
            no_friends_list = no_friends_list.append(user_id, ignore_index = True)
            
    # Write to file
    fm_friends.to_csv(directory + "data_derived/friends/processed_fm_tweeter_friends/friends_fm_" + i + ".csv", index = False)
    no_friends_list.to_csv(directory + "data_derived/friends/nofriends_fm_tweeters/nofriends_fm_tweeters_" + i + ".csv", index = False)

# Parse exposed followers
pool = mp.Pool(mp.cpu_count())
for i in range(100):
    pool.apply(parse_exposed_followers, args = (i, fm_tweeters, data_directory))
pool.close()



############################## Count ##############################

# path to data
path_to_tweets = data_directory + "data_derived/tweets/"
path_to_exposed = data_directory + "data_derived/followers/processed_exposed_followers/"
path_nofollowers_fm_tweeters = data_directory + "data_derived/followers/nofollowers_fm_tweeters/"
path_to_fm_friends = data_directory + "data_derived/friends/processed_fm_tweeter_friends/"
path_to_nofriends_fm_tweeters = data_directory + "data_derived/friends/nofriends_fm_tweeters/"

# Desired outpaths for writing new data files
outpath = data_directory + "data_derived/"
outpath_follower_data = data_directory + "data_derived/followers/"
outpath_tweet_data = data_directory + "data_derived/tweets/"
outpath_friend_data = data_directory + "data_derived/friends/"


####################
# Count up FM tweeters
####################
count_fm_tweeters = fm_tweeters.shape[0]


####################
# Count up unique followers
####################
# Get unique IDs of tweeters
tweeters = np.unique(tweets['user_id'])
tweeters = tweeters.astype(int)
count_tweeters = len(tweeters)

# Count up followers
exposed_files = sorted( os.listdir(path_to_exposed) )
num_files = len(exposed_files)
all_exposed = np.array([], dtype = int)
for file in exposed_files:
    data = np.loadtxt(path_to_exposed + file, skiprows = 1, dtype = int) #first row is header
    all_exposed = np.append(all_exposed, data)
    all_exposed = np.unique(all_exposed)
    del(data)
all_exposed_minus_tweeters = np.setdiff1d(all_exposed, fm_tweeters)
count_all_exposed = all_exposed.shape[0]
count_all_exposed_minus_tweeters = all_exposed_minus_tweeters.shape[0]
all_exposed = pd.DataFrame(all_exposed, columns = ["user_id"])
all_exposed_minus_tweeters = pd.DataFrame(all_exposed_minus_tweeters, columns = ["user_id"])  

# Count FM tweeters who had no followers
exposed_nofollower_files = sorted( os.listdir(path_nofollowers_fm_tweeters) )
no_followers = np.array([], dtype = int)
for file in exposed_nofollower_files:
    data = np.loadtxt(path_nofollowers_fm_tweeters + file, skiprows = 1, dtype = int) #first row is header
    no_followers = np.append(no_followers, data)
    no_followers = np.unique(no_followers)
    del(data)
count_no_followers = no_followers.shape[0]
no_followers = pd.DataFrame(no_followers, columns = ["user_id"])


####################
# Count up friends of FM tweeters
####################
# Count up followers
friend_files = sorted( os.listdir(path_to_fm_friends) )
num_files = len(friend_files)
all_friends = np.array([], dtype = int)
for file in friend_files:
    data = np.loadtxt(path_to_fm_friends + file, skiprows = 1, dtype = int) #first row is header
    all_friends = np.append(all_friends, data)
    all_friends = np.unique(all_friends)
    del(data)
all_friends_minus_tweeters = np.setdiff1d(all_friends, fm_tweeters)
count_all_friends = all_friends.shape[0]
count_all_friends_minus_tweeters = all_friends_minus_tweeters.shape[0]
all_friends = pd.DataFrame(all_friends, columns = ["user_id"])
all_friends_minus_tweeters = pd.DataFrame(all_friends_minus_tweeters, columns = ["user_id"])
    
# Count FM tweeters who have no friends
fm_tweeters_nofriends_files = sorted( os.listdir(path_to_nofriends_fm_tweeters) )
no_friends = np.array([], dtype = int)
for file in fm_tweeters_nofriends_files:
    data = np.loadtxt(path_to_nofriends_fm_tweeters + file, skiprows = 1, dtype = int) #first row is header
    no_friends = np.append(no_friends, data)
    no_friends = np.unique(no_friends)
    del(data)
count_no_friends = no_friends.shape[0]
no_friends = pd.DataFrame(no_friends, columns = ["user_id"])


####################
# Save
####################
# Measure number of users and write out to file
unique_users = pd.DataFrame({'user_type': ["Tweeters", "FM tweeters", 
                                           "Followers exposed to FM Articles", "Exposed followers excluding FM tweeters", "FM tweeters w/o followers",
                                           "Friends of FM tweeters", "Friends of FM tweeters excluding FM tweeters", "FM tweeters w/o friends"], 
                             'count': [count_tweeters, count_fm_tweeters, 
                                       count_all_exposed, count_all_exposed_minus_tweeters, count_no_followers,
                                       count_all_friends, count_all_friends_minus_tweeters, count_no_friends]})
unique_users.to_csv(outpath + "exposed_user_summary.csv", index = False)

all_exposed.to_csv(outpath_follower_data + "exposed_fakenews_followers.csv", index = False)
all_exposed_minus_tweeters.to_csv(outpath_follower_data + "exposed_fakenews_followers_excl_fmtweeters.csv", index = False)
no_followers.to_csv(outpath_follower_data + "fm_tweeter_nofollowers.csv", index = False)
fm_tweeters.to_csv(outpath_tweet_data + "unique_fm_tweeters.csv", index = False)
all_friends.to_csv(outpath_friend_data + "friends_fm_tweeters.csv", index = False)
all_friends_minus_tweeters.to_csv(outpath_friend_data + "friends_fm_tweeters_excl_fmtweeters.csv", index = False)
no_friends.to_csv(outpath_friend_data + "fm_tweeter_nofriends.csv", index = False)