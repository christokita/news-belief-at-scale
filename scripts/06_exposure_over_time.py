#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 12:00:20 2020

@author: ChrisTokita

SCRIPT:
Determine exposure by political type over time.
"""


####################
# Load libraries and packages, set important parameter
####################
import pandas as pd
import numpy as np
import re
import os
import math
import sys

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
    selected_tweets = selected_tweets.sort_values('relative_tweet_time')
    
    # Load follower IDs for each tweeter
    followers_per_tweeter = gather_followers(selected_tweets)
    
    # Count total followers for each user
    follower_count = pd.DataFrame(followers_per_tweeter['user_id'].value_counts())
    follower_count = follower_count.rename(columns = {'user_id': 'followers'})
    follower_count['user_id'] = follower_count.index.astype('int64')
    follower_count = follower_count.reset_index(drop = True)
    
    # Drop duplicates so we only keep first instance of exposure per follower
    followers_per_tweeter = followers_per_tweeter.drop_duplicates(subset = "follower_id")
    
    # Merge in follower ideology data
    follower_ideologies = pd.read_csv(data_directory + "data_derived/ideological_scores/cleaned_followers_ideology_scores.csv",
                                      dtype = {'user_id': 'int64', 'pablo_score': float})
    follower_ideologies = follower_ideologies.rename(columns = {'user_id': 'follower_id'})
    follower_ideologies = follower_ideologies.drop(columns = ['accounts_followed'])
    followers_per_tweeter = followers_per_tweeter.merge(follower_ideologies, how = "left", on = "follower_id")
    
    # Determine bin edges from size. Ideologies are in range [-2.654, 5.001]
    ideol_bins, bin_labels = create_ideology_bins(follower_ideologies, ideol_bin_size)
    
    # Create data frame to collect exposure data
    cols = ['relative_time', 'tweet_number', 'tweet_id', 'user_id', 'follower_count', 'new_exposed_users', 'cumulative_exposed'] + list(bin_labels)
    exposed_over_time = pd.DataFrame(columns = cols)
    exposed_over_time['tweet_id'] = exposed_over_time['tweet_id'].astype('int64')
    exposed_over_time['user_id'] = exposed_over_time['user_id'].astype('int64')
    del cols
    
    # Construct our data set
    exposed_over_time = create_exposure_dataset(collection_df = exposed_over_time, 
                                                follower_count = follower_count,
                                                article_tweets = selected_tweets, 
                                                follower_data = followers_per_tweeter,
                                                bins = ideol_bins,
                                                bin_labels = bin_labels)
    
    
    # Add additional info for data frame and return
    exposed_over_time['total_article_number'] = story_id
    return(exposed_over_time)
    
    
def gather_followers(article_tweets):
    """
    Function that will go through each tweet and will load the follower IDs for each tweeter.
    
    OUTPUT
    - followers_exposed: dataframe listing user ID of each tweeter and all user IDs of their respective followers.
    """
    # Load follower IDs for each tweeter
    follower_files = os.listdir(data_directory + "data/followers/")
    followers_exposed = pd.DataFrame(columns = ['user_id', 'follower_id'], dtype = 'int64')
    for j in np.arange(article_tweets.shape[0]):
        
        # Get user ID and appropriate file
        user_id = article_tweets['user_id'].iloc[j]
        regex = re.compile(r"[0-9].*_%s.csv" % user_id)
        file = list(filter(regex.match, follower_files))
        
        # Load followers, bind to dataframe
        followers = load_followers(file, data_directory)
        if len(followers) > 0:
            new_row = pd.DataFrame({'user_id': user_id, 'follower_id': followers}, dtype = 'int64')
            followers_exposed = followers_exposed.append(new_row, ignore_index = True)
    return followers_exposed
    
    
def load_followers(file, data_directory):
    """
    Function that will load the follower files or return an empty array if there are no followers for this user

    OUTPUT
    - followers:   array of follower user IDs (numpy array, str)
    """
    if len(file) == 0: #no followers, no follower file
        followers = np.array([])
    else:
        followers = np.genfromtxt(data_directory + "data/followers/" + file[0], dtype = 'int64')
        try:
            followers = followers[1:len(followers)] #remove header, will raise error if empty
        except:
            followers = np.array([]) #no followers, empty file
    return followers


def create_exposure_dataset(collection_df, follower_count, article_tweets, follower_data, bins, bin_labels):
    """
    Function that will go construct our dataset
    """
    cumulative_exposed = 0
    for j in np.arange(article_tweets.shape[0]):
        
        # Get tweet identifying items
        user_id = article_tweets['user_id'].iloc[j]
        tweet_id = article_tweets['tweet_id'].iloc[j]
        rel_time = article_tweets['relative_tweet_time'].iloc[j]
        tweet_number = j
        
        # Get follower count
        if user_id in set(follower_count.user_id):
            total_followers = follower_count.followers[follower_count.user_id == user_id].iloc[0]
        else:
            total_followers = 0
            
        # Subset out unique, newly exposed followers of this user
        their_followers = follower_data[follower_data['user_id'] == user_id]
        
        # Count total new exposed, update cumulative exposed
        new_exposed_count = their_followers.shape[0]
        cumulative_exposed += new_exposed_count
        
        # Get ideologies of these followers     
        ideol_counts = np.histogram(their_followers['pablo_score'], bins = bins)
        ideol_counts = pd.DataFrame([ideol_counts[0]], columns = bin_labels)
        
        # Summarise and create new row
        new_row = pd.DataFrame({'relative_time': rel_time, 
                                'tweet_number': tweet_number, 
                                'tweet_id': tweet_id,
                                'user_id': user_id,
                                'follower_count': total_followers,
                                'new_exposed_users': new_exposed_count,
                                'cumulative_exposed': cumulative_exposed}, index = [0])
        new_row = pd.concat([new_row, ideol_counts], axis = 1, sort = False)
        collection_df = collection_df.append(new_row, ignore_index = True, sort = False)
    
    # return    
    return collection_df
        

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
                                          'user_id': 'int64', 'tweet_id': 'int64'})
    
    # Get unique articles
    unique_articles = labeled_tweets['total_article_number'].unique()
    unique_articles.sort()
    unique_articles = unique_articles[~np.isnan(unique_articles)] #drop nan
    unique_articles = unique_articles.astype(int)
    
    # Check specific article. If it is already processed, skip proccessing and finish script
    story = unique_articles[i]
    processed_articles = os.listdir(data_directory + "data_derived/exposure/individual_articles/")
    processed_articles = [file for file in processed_articles if re.match('^article', file)] #filter out hidden copies of same files
    processed_articles = [re.search('([0-9]+)', file).group(1) for file in processed_articles]
    processed_articles = np.array(processed_articles, dtype = int)
    if story in processed_articles:
        sys.exit(0)
    
    # Process exposure over time
    story_exposed = unique_exposed_over_time(story_id = story, 
                                             tweets = labeled_tweets, 
                                             data_directory = data_directory, 
                                             ideol_bin_size = 0.5)
    story_exposed.to_csv(data_directory + "data_derived/exposure/individual_articles/article_" + str(story) + ".csv", index = False)

    # Check if all stories are processed at this point.
    processed_articles = os.listdir(data_directory + "data_derived/exposure/individual_articles/")
    processed_articles = [file for file in processed_articles if re.match('^article', file)] #filter out hidden copies of same files
    processed_articles = [re.search('([0-9]+)', file).group(1) for file in processed_articles]
    processed_articles = np.array(processed_articles, dtype = int)
    articles_left_to_do = np.setdiff1d(unique_articles, processed_articles)
    
    # If all articles processed, compile into main dataframe and write out
    if len(articles_left_to_do) == 0:   
        col_names = story_exposed.columns
        exposed_timeseries = pd.DataFrame(columns = col_names)
        article_files = os.listdir(data_directory + "data_derived/exposure/individual_articles/")
        article_files = [file for file in article_files if re.match('^article', file)] #filter out hidden copies of same files
        for file in article_files:
            story_exposed = pd.read_csv(data_directory + "data_derived/exposure/individual_articles/" + file,
                                        dtype = {'user_id': 'int64'})
            exposed_timeseries = exposed_timeseries.append(story_exposed, sort = False)
        exposed_timeseries = exposed_timeseries.to_csv(data_directory + "data_derived/exposure/users_exposed_over_time.csv", index = False)