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
import time
from scipy import stats

# Get chunk number from SLURM script
i = int(sys.argv[1])

####################
# Functions to determine unique exposed users over time
####################

def unique_exposed_over_time(story_id, tweets, data_directory, ideol_bin_size, mean_exposure_time, sd_exposure_time):
    """
    Function that will determine how many unique users had been exposed to a story over time.
    
    Parameters:
        story_id (int):             id of article of interest.
        tweets (dataframe):         dataframe containing all tweets.
        data_directory (str):       path to high-level data directory.
        ideol_bin_size (float):     width of bins for ideology, with break at zero.
        mean_expsoure_time (float): mean time for a user to be exposed to a tweet from someone they follow.
        sd_exposure_time (float):   standard devistion of time for a user to be exposed to a tweet from someone they follow.
    
    Returns:
        exposed_over_time (dataframe): contains the time series of tweets and exposure over time for that story.
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
    selected_tweets['tweet_number'] = np.arange(selected_tweets.shape[0])
    
    # We will now create a list of uniquely exposed followers, matched to ideology.
    # Becuase we want to also count users that we do not have ideology scores for, 
    # we will not use our paired tweeter-follower ideology dataset crated in script 5a
    
    # Load follower IDs for each tweeter and include important tweet metadata
    tweets_with_followers = match_followers_to_tweet(article_tweets = selected_tweets, 
                                                     story_id = story_id, 
                                                     data_directory = data_directory)
    
    # Count total followers for each user
    follower_count = pd.DataFrame(tweets_with_followers[['user_id', 'follower_id']].drop_duplicates()['user_id'].value_counts()) #make sure to get filter first to unique tweeter-follower pairs before counting
    follower_count = follower_count.rename(columns = {'user_id': 'followers'})
    follower_count['user_id'] = follower_count.index.astype('int64')
    follower_count = follower_count.reset_index(drop = True)
    
    # Account for delay in exposure to each follower
    exposure_time_offset = time_distribute_exposure(n_followers = tweets_with_followers.shape[0],
                                                    mean_exposure_time = mean_exposure_time,
                                                    sd_exposure_time = sd_exposure_time)
    tweets_with_followers['relative_exposure_time'] = tweets_with_followers['relative_tweet_time'] + exposure_time_offset
    
    # Drop duplicates so we only keep first instance of exposure per follower
    tweets_with_followers = tweets_with_followers.sort_values('relative_exposure_time')
    tweets_with_followers = tweets_with_followers.drop_duplicates(subset = "follower_id", keep = "first")
    
    # Load ideology scores--of both followers and tweeters (since tweeters can also be followers)
    follower_ideologies = pd.read_csv(data_directory + "data_derived/ideological_scores/cleaned_followers_ideology_scores.csv",
                                      dtype = {'user_id': 'int64', 'pablo_score': float})
    follower_ideologies = follower_ideologies.rename(columns = {'user_id': 'follower_id'})
    follower_ideologies = follower_ideologies.drop(columns = ['accounts_followed']) 
    tweeter_ideologies = tweets[['user_id', 'user_ideology']]
    tweeter_ideologies = tweeter_ideologies.rename(columns = {'user_id': 'follower_id', 'user_ideology': 'pablo_score'})
    tweeter_ideologies = tweeter_ideologies[~pd.isna(tweeter_ideologies['pablo_score'])]
    follower_ideologies = follower_ideologies.append(tweeter_ideologies, ignore_index = True)
    follower_ideologies = follower_ideologies.drop_duplicates()
    del tweeter_ideologies
    
    # Merge in follower ideology data
    tweets_with_followers = tweets_with_followers.merge(follower_ideologies, how = "left", on = "follower_id")

    # Determine bin edges from size. Ideologies are in range [-2.654, 5.009]
    ideol_bins, bin_labels = create_ideology_bins(follower_ideologies.pablo_score, ideol_bin_size)
    
    # Create data frame to collect exposure data
    cols = ['relative_time', 'tweet_number', 'tweet_id', 'user_id', 'follower_count', 'new_exposed_users', 'cumulative_exposed'] + list(bin_labels)
    exposed_over_time = pd.DataFrame(columns = cols)
    exposed_over_time['tweet_id'] = exposed_over_time['tweet_id'].astype('int64')
    exposed_over_time['user_id'] = exposed_over_time['user_id'].astype('int64')
    del cols
    
    # Construct our data set
    selected_tweets = selected_tweets.sort_values('relative_tweet_time') #paranoid double check that tweets are properly ordered
    exposed_over_time = create_exposure_dataset(collection_df = exposed_over_time, 
                                                follower_count = follower_count,
                                                article_tweets = selected_tweets, 
                                                follower_data = tweets_with_followers,
                                                bins = ideol_bins,
                                                bin_labels = bin_labels)
    
    
    # Add additional info for data frame and return
    exposed_over_time['total_article_number'] = story_id
    return(exposed_over_time)
    
    
def match_followers_to_tweet(article_tweets, story_id, data_directory):
   """
   Function that will load followers that would have seen a given tweet.
   
   Parameters:
       article_tweets (dataframe): all tweets for this story/article.
       story_id (int):             id of story/article of interest.
       data_directory (str):       path to high-level data directory.
   
   Returns:
       matched_followers (dataframe):   dataframe listing tweet ID with the set of follower IDs that could have potentially seen it (i.e., the followers of the user who tweeted it).
   """
   
   # Set place to store compiled lists
   data_subdir = data_directory + 'data_derived/exposure/compiled_exposed_follower_lists/'
   os.makedirs(data_subdir, exist_ok = True)
   
   # If file already exists, load; otherwise, compile and save.
   file_name = data_subdir + 'exposedfollowers_story' + str(story_id) + '.csv'
   col_names = ['user_id', 'follower_id', 'tweet_id', 'tweet_time', 'tweet_number', 'relative_tweet_time']
   if not os.path.exists(file_name):
   
       # Get list of unique users and tweets in this set 
       unique_users = article_tweets['user_id'].unique()
       tweet_list = article_tweets[['user_id', 'tweet_id', 'tweet_time', 'tweet_number', 'relative_tweet_time']].copy()
       del article_tweets
   
       # Create paired list of tweeters and followers, batching out to a temp csv   
       follower_files = os.listdir(data_directory + "data/followers/")
       for user_id in unique_users:
           regex = re.compile(r"[0-9].*_%s.csv" % user_id)
           file = list(filter(regex.match, follower_files))
           followers = load_followers(file, data_directory)
           matched_followers = pd.DataFrame({'user_id': user_id, 'follower_id': followers}, dtype = 'int64')
           matched_followers = matched_followers.merge(tweet_list, how = "left", on = "user_id")
           matched_followers.to_csv(file_name, mode = "a", index = False, header = False)
           col_names = matched_followers.columns
           del matched_followers, followers
       
   # Load in compiled list, sort, and return
   matched_followers = pd.read_csv(file_name, header = None, names = col_names)
   matched_followers['tweet_time'] = pd.to_datetime(matched_followers['tweet_time'], format = '%Y-%m-%d %H:%M:%S%z')
   matched_followers = matched_followers.sort_values(by = ['tweet_number'])
   return matched_followers
    
    
def load_followers(file, data_directory):
    """
    Function that will load the follower files or return an empty array if there are no followers for this user.

    Parameters:
        file (str):           file name of this user's follower data file.
        data_directory (str): path to high-level data directory.
        
    Returns:
        followers (int64): array of follower user IDs.
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


def time_distribute_exposure(n_followers, mean_exposure_time, sd_exposure_time):
    """
    Distribute the exposure to a tweet over the time following the tweets, so that exposure for all followers isn't instantaneous.
    For now, we assume a truncated normal distribution cut off at time +0 and time +48 (the upper bound shouldnbe sufficiently large).
    
    Parameters:
        n_followers (int):          number of followers for which to calcualte the exposure time delay.
        mean_expsoure_time (float): mean of truncated normal distribution.
        sd_exposure_time (float):   standard devistion of truncated normal distribution.
        
    Returns:
        added_exposure_time (float): vector that describes at what point after the tweet that the follower was exposed.
    """
    
    lower, upper = 0, 48
    a = (lower - mean_exposure_time) / sd_exposure_time
    b = (upper - mean_exposure_time) / sd_exposure_time
    exposure_dist = stats.truncnorm(a, b, loc = mean_exposure_time, scale = sd_exposure_time)
    added_exposure_time = exposure_dist.rvs(n_followers)
    return added_exposure_time


def create_exposure_dataset(collection_df, follower_count, article_tweets, follower_data, bins, bin_labels):
    """
    Function that will go construct our dataset.
    
    Parameters:
        collection_df (dataframe):  dataframe that will collect our data.
        follower_count (dataframe): data for how many followers each tweeter has.
        article_tweets (dataframe): all tweets for this story/article.
        follower_data (dataframe):  list of all uniquely exposed followers--i.e., each unique follower only apears once at the point they were exposed--and their ideologies, for those that we have that data.
        bins (numpy array):         list of all bin edges for binning follower ideologies.
        bin_labels (numpy array):   list of labels for each bin.
        
    Returns:
        collection_df (dataframe): completed dataset of exposure.
    """
    cumulative_exposed = 0
    already_tweeted = [] #list to keep track of which users have already shown up sequentially
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
        
        # Bin ideologies of exposed followers     
        ideol_counts = np.histogram(their_followers['pablo_score'], bins = bins)
        ideol_counts = pd.DataFrame([ideol_counts[0]], columns = bin_labels)
        
        # Count total new exposed and update cumulative exposed
        # If user had already shared tweeted story before, then we can't double count their followers as exposed.
        if user_id in already_tweeted:
            new_exposed_count = 0
            ideol_counts[:] = 0
        else:
            new_exposed_count = their_followers.shape[0]
        cumulative_exposed += new_exposed_count
        already_tweeted.append(user_id)
        
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
        

def create_ideology_bins(ideologies, bin_size):
    """
    Function that computes and creates bins for ideology scores. "Center" bin break is at 0.
    
    Parameters:
        ideologies (pandas series): list of all ideology scores that we have for followers.
        bin_size (float):           desired width of bins for ideology.
    
    Returns:
        bins (numpy array):     bin edges .
        labels (list of str):   labels for bins.
    """
    lower_edge = math.floor(ideologies.min() / bin_size) * bin_size
    upper_edge = math.ceil(ideologies.max() / bin_size) * bin_size
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
                                             ideol_bin_size = 0.5,
                                             mean_exposure_time = 1,
                                             sd_exposure_time = 2)
    story_exposed.to_csv(data_directory + "data_derived/exposure/individual_articles/article_" + str(story) + ".csv", index = False)

    # Check if all stories are processed at this point.
    time.sleep(10)
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