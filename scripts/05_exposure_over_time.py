#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 10:05:11 2020

@author: ChrisTokita

SCRIPT:
Count the users exposed to a given article over time
"""

####################
# Load libraries and packages, set important parameters, load data
####################
import pandas as pd
import numpy as np
import re
import os
import multiprocessing as mp

####################
# Functions to determine unique exposed users over time
####################

def unique_exposed_over_time(story_id, tweets, data_directory):
    # Function that will determine how many unique users had been exposed to a story over time.
    #
    # OUTPUT:
    # - exposed_over_time: a dataframe that contains the time series of tweets and exposure over time for that story.

    # Filter tweets to those sharing story and sort by time
    selected_tweets = tweets[tweets.total_article_number == story_id].copy()
    selected_tweets['tweet_time'] = pd.to_datetime(selected_tweets['tweet_time'], format = '%a %b %d %H:%M:%S %z %Y')
    selected_tweets = selected_tweets.sort_values(by = 'tweet_time')
    
    # Calculate time since first share
    selected_tweets['article_first_time'] = min(selected_tweets['tweet_time'])
    selected_tweets['relative_tweet_time'] = selected_tweets['tweet_time'] - selected_tweets['article_first_time']
    selected_tweets['relative_tweet_time'] = selected_tweets['relative_tweet_time'] / np.timedelta64(1, 'h') #convert to hours
    
    # Create data frame to collect exposure data
    exposed_over_time = pd.DataFrame(columns = ['time', 'tweet_number', 'user_id', 'new_exposed_users', 'cumulative_exposed'])
    
    # Go through each tweet in order and determine unique exposed users
    exposed_already = np.array([], dtype = int)
    follower_files = os.listdir(data_directory + "data/followers/")
    for i in np.arange(selected_tweets.shape[0]):
        
        # Get user ID and load followers
        user_id = selected_tweets['user_id'].iloc[i]
        regex = re.compile(r"[0-9].*_%s.csv" % user_id)
        file = list(filter(regex.match, follower_files))
        if len(file) == 0: #no followers, no follower file
            follower_list = np.array([])
        else:
            follower_list = np.genfromtxt(data_directory + "data/followers/" + file[0], dtype = int)
            try:
                follower_list = follower_list[1:len(follower_list)] #remove header, will raise error if empty
            except:
                follower_list = np.array([]) #no followers, empty file
            
        
        # Determine new unique exposed users
        new_exposed = np.setdiff1d(follower_list, exposed_already)
        exposed_already = np.append(exposed_already, new_exposed)
        
        # Summarise and append to data set
        exposed_over_time = exposed_over_time.append({'time': selected_tweets['relative_tweet_time'].iloc[i], 
                                                      'tweet_number': i, 
                                                      'user_id': selected_tweets['user_id'].iloc[i],
                                                      'new_exposed_users': len(new_exposed),
                                                      'cumulative_exposed': len(exposed_already)}, 
                                                    ignore_index = True)
    # Add additional info for data frame and return
    exposed_over_time['total_article_number'] = story_id
    return(exposed_over_time)
    
    
####################
# Main script: Count exposed users over time
####################
if __name__ == '__main__':
    
    # high level directory (external HD or cluster storage)
    data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#    data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD
    
    # Load tweet data, esnure in proper format
    labeled_tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                                 dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                                          'user_id': 'Int64', 'tweet_id': 'Int64', 
                                          'retweeted_user_id': 'Int64', 'retweet_id': 'Int64',
                                          'quoted_user_id': 'Int64', 'quoted_id': 'Int64'})
    
    # Get unique articles
    unique_articles = labeled_tweets['total_article_number'].unique()
    unique_articles.sort()
    unique_articles = unique_articles[~np.isnan(unique_articles)] #drop nan
    unique_articles = unique_articles.astype(int)
    
    # Check if some articles have already been processed
    # If so, remove from list of articles to process
    processed_articles = os.listdir(data_directory + "data_derived/timeseries/individual_articles/")
    processed_articles = [file for file in processed_articles if re.match('^article', file)] #filter out hidden copies of same files
    processed_articles = [re.search('([0-9]+)', file).group(1) for file in processed_articles]
    processed_articles = np.array(processed_articles, dtype = int)
    unique_articles = np.setdiff1d(unique_articles, processed_articles)
    
    # Process articles in parallel, saving the data for each article individually
    pool = mp.Pool(mp.cpu_count())
    for story in unique_articles:
        story_exposed = pool.apply(unique_exposed_over_time, args = (story, labeled_tweets, data_directory))
        story_exposed.to_csv(data_directory + "data_derived/timeseries/individual_articles/article_" + str(story) + ".csv", index = False)
    pool.close()
    pool.join()

    # Compile into main dataframe and write out
    exposed_timeseries = pd.DataFrame(columns = ['time', 'tweet_number', 'user_id', 'new_exposed_users', 'cumulative_exposed', 'total_article_number'])
    article_files = os.listdir(data_directory + "data_derived/timeseries/individual_articles/")
    article_files = [file for file in article_files if re.match('^article', file)] #filter out hidden copies of same files
    for file in article_files:
        story_exposed = pd.read_csv(data_directory + "data_derived/timeseries/individual_articles/" + file,
                                    dtype = {'user_id': int})
        exposed_timeseries = exposed_timeseries.append(story_exposed, sort = False)
    exposed_timeseries = exposed_timeseries.to_csv(data_directory + "data_derived/timeseries/users_exposed_over_time.csv", index = False)