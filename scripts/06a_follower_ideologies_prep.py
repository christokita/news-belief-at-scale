#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 14:28:50 2020

@author: ChrisTokita

SCRIPT:
Prepare data for analysis of all follower ideologies
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import os
import re
import math
import multiprocessing as mp


# Function to load followers
def load_followers(file, data_directory):
    """
    Function that will load the follower files or return an empty array if there are no followers for this user

    OUTPUT
    - followers:   array of follower user IDs (numpy array, str)
    """
    if len(file) == 0: #no followers, no follower file
        followers = np.array([])
    else:
        followers = np.genfromtxt(data_directory + "data/followers/" + file[0], dtype = str)
        try:
            followers = followers[1:len(followers)] #remove header, will raise error if empty
        except:
            followers = np.array([]) #no followers, empty file
    return followers

# Create dataframe of tweeters and their followers
def match_followers_to_ideologies(user_ids, ideologies, data_directory, outpath, batch, n_batches):
    """
    Function that will...
    
    OUTPUT
    
    """
    
    # Grab proper batch of user_ids
    batch_size =  math.ceil(len(user_ids) / n_batches)
    start = batch*batch_size
    end = (batch+1)*batch_size
    user_id_batch = user_ids[start:end]
    
    # Loop through tweet user IDs and grab all follower IDs
    follower_files = os.listdir(data_directory + "data/followers/")
    followers_matched = pd.DataFrame(columns = ['tweeter_id', 'follower_id'], dtype = object)
    for user in user_id_batch:
    
        # Load followers
        regex = re.compile(r"[0-9].*_%s.csv" % user)
        file = list(filter(regex.match, follower_files))
        followers = load_followers(file, data_directory) 
        
        # Append to dataset
        new_data = pd.DataFrame({'tweeter_id': user, 'follower_id': followers}, dtype = object)
        followers_matched = followers_matched.append(new_data, ignore_index = True)
        
    # Merge in ideologies and save
    del new_data, followers, follower_files 
    followers_matched = followers_matched.merge(ideologies, how = "inner", on = 'follower_id')
    followers_matched = followers_matched.drop(columns = ['follower_id'])
    followers_matched = followers_matched.rename(columns = {'tweeter_id': 'user_id', 'pablo_score': 'follower_ideology'})
    followers_matched = followers_matched.sort_values(by = 'user_id')
    followers_matched.to_csv(outpath + 'temp_matched_follower_ideology_' + str(batch) + '.csv', index = False)
    return None

####################
# Create user x follower ideology dataset
####################
if __name__ == '__main__':
    
    # high level directory (external HD or cluster storage)
    data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#    data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD
    
    # Load tweets and follower ideology data 
    tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                         dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                                  'user_id': object, 'tweet_id': object, 
                                  'retweeted_user_id': object, 'retweet_id': object,
                                  'quoted_user_id': object, 'quoted_id': object})
        
    follower_ideologies = pd.read_csv(data_directory + "data_derived/ideological_scores/cleaned_followers_ideology_scores.csv",
                                      dtype = {'user_id': object, 'pablo_score': float})
    follower_ideologies = follower_ideologies.drop(columns = ['accounts_followed'])
    follower_ideologies = follower_ideologies.rename(columns = {'user_id': 'follower_id'})

    # Get unique users
    unique_users = tweets.user_id[~pd.isna(tweets['total_article_number'])].unique()
    del tweets
    
    # Create temporary folder to hold partial results
    temp_dir = data_directory + 'data_derived/ideological_scores/TEMP_matched_follower_ideologies/'
    os.makedirs(temp_dir, exist_ok = True)

    # Process tweeters in batches and write temporary files out
    pool = mp.Pool(mp.cpu_count())
    n_batches = 20
    for i in np.arange(n_batches):
        pool.apply_async(match_followers_to_ideologies, args = (unique_users, follower_ideologies, data_directory, temp_dir, i, n_batches))
    pool.close()
    pool.join()
    
    # Bind together temporary out files into final dataset
    temp_files = os.listdir(temp_dir)
    for file in temp_files:
        temp_data = pd.read_csv(temp_dir + file, dtype = {'user_id': object})
        if 'followers_data' not in globals():
            followers_data = temp_data
        else:
            followers_data = followers_data.append(temp_data, ignore_index = True, sort = False)
        del temp_data
    
    # save
    followers_data.to_csv(data_directory + 'data_derived/ideological_scores/paired_tweeter-follower_ideology.csv', index = False)

    # Delete temp data
    for file in temp_files:
        os.remove(temp_dir + file)
    os.rmdir(temp_dir)