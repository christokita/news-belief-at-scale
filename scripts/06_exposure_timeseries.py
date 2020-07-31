#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 12:00:20 2020

@author: ChrisTokita

SCRIPT:
Determine exposure by political type over time.
"""


####################
# Load libraries and packages, set important parameters, load data
####################
import pandas as pd
import numpy as np
import re
import os
import math
import sys
#import multiprocessing as mp

# Get chunk number from SLURM script
i = int(sys.argv[1])

####################
# Functions to determine unique exposed users over time
####################

def unique_exposed_over_time(story_id, tweets, data_directory, ideol_bin_size):
    """
    Function that will determine how many unique users had been exposed to a story over time.
    
    OUTPUT:
    - exposed_over_time: a dataframe that contains the time series of tweets and exposure over time for that story.
    """
    
    # Filter tweets to those sharing story and sort by time
    selected_tweets = tweets[tweets.total_article_number == story_id].copy()
    selected_tweets['tweet_time'] = pd.to_datetime(selected_tweets['tweet_time'], format = '%a %b %d %H:%M:%S %z %Y')
    selected_tweets = selected_tweets.sort_values(by = 'tweet_time')
    
    # Calculate time since first share
    selected_tweets['article_first_time'] = min(selected_tweets['tweet_time'])
    selected_tweets['relative_tweet_time'] = selected_tweets['tweet_time'] - selected_tweets['article_first_time']
    selected_tweets['relative_tweet_time'] = selected_tweets['relative_tweet_time'] / np.timedelta64(1, 'h') #convert to hours
    
    # Load follower ideology data
    follower_ideologies = pd.read_csv(data_directory + "data_derived/ideological_scores/cleaned_followers_ideology_scores.csv",
                                      dtype = {'user_id': object, 'pablo_score': float})
    
    # Determine bin edges from size. Ideologies are in range [-2.654, 5.001]
    ideol_bins, bin_labels = create_ideology_bins(follower_ideologies, ideol_bin_size)
    
    # Create data frame to collect exposure data
    cols = ['time', 'tweet_number', 'user_id', 'new_exposed_users', 'cumulative_exposed'] + list(bin_labels)
    exposed_over_time = pd.DataFrame(columns = cols)
    del cols
    
    # Go through each tweet in order and determine unique exposed users, and their ideology
    exposed_already = np.array([], dtype = object)
    follower_files = os.listdir(data_directory + "data/followers/")
    for j in np.arange(selected_tweets.shape[0]):
        
        # Get user ID and appropriate file
        user_id = selected_tweets['user_id'].iloc[j]
        regex = re.compile(r"[0-9].*_%s.csv" % user_id)
        file = list(filter(regex.match, follower_files))
        
        # Load followers
        followers = load_followers(file, data_directory)
        
        # Determine new unique exposed users
        new_exposed, exposed_already = determine_new_exposed(followers, exposed_already)
        
        # Get ideologies of these followers     
        ideologies = follower_ideologies.loc[follower_ideologies['user_id'].isin(new_exposed), 'pablo_score']
        ideol_counts = np.histogram(ideologies, bins = ideol_bins)
        ideol_counts = pd.DataFrame([ideol_counts[0]], columns = bin_labels)
        
        # Summarise and append to data set
        new_row = pd.DataFrame({'time': selected_tweets['relative_tweet_time'].iloc[j], 
                                                      'tweet_number': j, 
                                                      'user_id': selected_tweets['user_id'].iloc[j],
                                                      'new_exposed_users': len(new_exposed),
                                                      'cumulative_exposed': len(exposed_already)}, index = [0])
        new_row = pd.concat([new_row, ideol_counts], axis = 1, sort = False)
        exposed_over_time = exposed_over_time.append(new_row, ignore_index = True, sort = False)
        
    # Add additional info for data frame and return
    exposed_over_time['total_article_number'] = story_id
    return(exposed_over_time)
    
    
def create_ideology_bins(ideologes, bin_size):
    """
    Function that computes bins for ideology scores
    
    OUTPUT:
    - bins:     bin edges (numpy array).
    - labels:   labels for bins (list of str).
    """
    lower_edge = math.floor(ideologes['pablo_score'].min() / bin_size) * bin_size
    upper_edge = math.ceil(ideologes['pablo_score'].max() / bin_size) * bin_size
    bins = np.arange(lower_edge, upper_edge + 0.001, bin_size) #stop number in range is non-inclusive, so need to add small amount
    labels = np.delete(bins, -1) #delete upper edge of last bin
    labels = ["ideol_" + str(s) + "_" + str(s+0.5) for s in labels]
    return bins, labels
    

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

def determine_new_exposed(followers, exposed_already):
    """
    Function to determine what new unique followers have been exposed to the news article
    
    OUTPUT
    - new_exposed:  array of user IDs of users that were exposed for the first time. 
    """
    if len(followers) > 0: #if array is empty, np.setdiff1d will raise annoying error.
        new_exposed = np.setdiff1d(followers, exposed_already)
        exposed_already = np.append(exposed_already, new_exposed)
    else: #bypass annoying error with np.setdiff1d if followers is an empty array
        new_exposed = np.array([])
    return new_exposed, exposed_already

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
                                          'user_id': object, 'tweet_id': object, 
                                          'retweeted_user_id': object, 'retweet_id': object,
                                          'quoted_user_id': object, 'quoted_id': object})
    
    # Get unique articles
    unique_articles = labeled_tweets['total_article_number'].unique()
    unique_articles.sort()
    unique_articles = unique_articles[~np.isnan(unique_articles)] #drop nan
    unique_articles = unique_articles.astype(int)
    
    # Select specific article and process exposure over time
    story = unique_articles[i]
    story_exposed = unique_exposed_over_time(story_id = story, 
                                             tweets = labeled_tweets, 
                                             data_directory = data_directory, 
                                             ideol_bin_size = 0.5)
    story_exposed.to_csv(data_directory + "data_derived/timeseries/individual_articles/article_" + str(story) + ".csv", index = False)

    # Check if all stories are processed at this point.
    processed_articles = os.listdir(data_directory + "data_derived/timeseries/individual_articles/")
    processed_articles = [file for file in processed_articles if re.match('^article', file)] #filter out hidden copies of same files
    processed_articles = [re.search('([0-9]+)', file).group(1) for file in processed_articles]
    processed_articles = np.array(processed_articles, dtype = int)
    articles_left_to_do = np.setdiff1d(unique_articles, processed_articles)
    
    # If all articles processed, compile into main dataframe and write out
    if len(articles_left_to_do) == 0:   
        col_names = story_exposed.columns
        exposed_timeseries = pd.DataFrame(columns = col_names)
        article_files = os.listdir(data_directory + "data_derived/timeseries/individual_articles/")
        article_files = [file for file in article_files if re.match('^article', file)] #filter out hidden copies of same files
        for file in article_files:
            story_exposed = pd.read_csv(data_directory + "data_derived/timeseries/individual_articles/" + file,
                                        dtype = {'user_id': int})
            exposed_timeseries = exposed_timeseries.append(story_exposed, sort = False)
        exposed_timeseries = exposed_timeseries.to_csv(data_directory + "data_derived/timeseries/users_exposed_over_time.csv", index = False)