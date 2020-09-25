#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 14:28:50 2020

@author: ChrisTokita

SCRIPT:
Calculate the distributions of ideologies of the population and of each individual's followers
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

# high level directory (external HD or cluster storage)
#data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD

# Paths for ideology distrubtion data
outpath = data_directory + "data_derived/ideological_scores/estimated_ideol_distributions/"
pop_dist_file = outpath + "_population_distribution.csv"
ind_dist_file = outpath + "follower_ideology_distributions.csv"


####################
# Functions for analysis
####################
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
    Function that will take a set of user IDs and subset a batch of them,
    then match these users to their followers that have ideology scores
    
    OUTPUT
    - nothing. It writes the save this dataset to file.
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
    followers_matched = followers_matched.rename(columns = {'tweeter_id': 'user_id', 'pablo_score': 'follower_ideology'})
    followers_matched = followers_matched.sort_values(by = 'user_id')
    followers_matched.to_csv(outpath + 'temp_matched_follower_ideology_' + str(batch) + '.csv', index = False)
    return None

def create_ideology_bins(ideologies = None, bin_ends = None, bin_size = None):
    """
    Function that computes and creates bins for ideology scores. "Center" bin break is at 0.
    
    OUTPUT:
    - bins:     bin edges (numpy array).
    - labels:   labels for bins (list of str).
    """
    if ideologies is not None:
        lower_edge = math.floor(ideologies.min() / bin_size) * bin_size
        upper_edge = math.ceil(ideologies.max() / bin_size) * bin_size
        bins = np.arange(lower_edge, upper_edge + 0.001, bin_size) #stop number in range is non-inclusive, so need to add small amount
    elif bin_ends is not None:
        bins = np.arange(bin_ends[0], bin_ends[1] + 0.001, bin_size)
    labels = np.delete(bins, -1) #delete upper edge of last bin
    labels = ["ideol_" + str(s) + "_" + str(s+bin_size) for s in labels]
    return bins, labels


# Main script
if __name__ == '__main__':
    
    ####################
    # Load (or create) data set of paired tweeters-follower ideology scores
    ####################
    if os.path.exists(data_directory + 'data_derived/ideological_scores/paired_tweeter-follower_ideology.csv'):
        followers_data = pd.read_csv(data_directory + 'data_derived/ideological_scores/paired_tweeter-follower_ideology.csv',
                                     dtype = {'user_id': 'int64', 'follower_id': 'int64'})
    
    else:
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
        
        # Add tweeter ideologies (since they can also be followers) 
        tweeter_ideologies = tweets[['user_id', 'user_ideology']]
        tweeter_ideologies = tweeter_ideologies.rename(columns = {'user_id': 'follower_id', 'user_ideology': 'pablo_score'})
        tweeter_ideologies = tweeter_ideologies[~pd.isna(tweeter_ideologies['pablo_score'])]
        tweeter_ideologies = tweeter_ideologies.drop_duplicates()
        follower_ideologies = follower_ideologies.append(tweeter_ideologies, ignore_index = True)
        follower_ideologies = follower_ideologies.drop_duplicates()
        del tweeter_ideologies
        
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
            
    
    ####################
    # Calculate population-level distribution of ideologies
    ####################
    # Unique followers
    unique_followers_data = followers_data[['follower_id', 'follower_ideology']].drop_duplicates()
    
    # Determine bin edges from size. Ideologies are in range [-2.654, 5.009]
    ideol_bin_size = 0.5
    ideol_bins, bin_labels = create_ideology_bins(ideologies = unique_followers_data.follower_ideology, 
                                                  bin_size = ideol_bin_size)
    
    # Bin ideologies   
    bin_count, bin_edges = np.histogram(unique_followers_data.follower_ideology, bins = ideol_bins)
    ideol_dist = pd.DataFrame({'bin_label': bin_labels, 'count': bin_count, 'frequency': bin_count/sum(bin_count)})
    ideol_dist['bin_left_edge'] =  ideol_dist['bin_label'].str.extract(r"ideol_([-.0-9]+)_").astype(float)
    ideol_dist['bin_right_edge'] =  ideol_dist['bin_label'].str.extract(r"ideol_[-.0-9]+_([-.0-9]+)").astype(float)
    ideol_dist['bin_center'] = (ideol_dist['bin_left_edge'] + ideol_dist['bin_right_edge']) / 2
    ideol_dist = ideol_dist[['bin_label', 'bin_left_edge', 'bin_center', 'bin_right_edge', 'count', 'frequency']]

    # Plot
    import matplotlib.pyplot as plt
    plt.bar(x = ideol_dist['bin_center'], 
            height = ideol_dist['frequency'], 
            width = ideol_bin_size, 
            edgecolor = "white", align = "center")

    # Save
    ideol_dist.to_csv(pop_dist_file, index = False)
    del unique_followers_data
    
    ####################
    # Calculate distribution of each tweeter's followers
    ####################
    # Load tweets
    tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                     dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                              'user_id': 'int64', 'tweet_id': 'int64'})
    tweeters = tweets['user_id'].unique()
    del tweets
    
    # Create bins
    ideol_bin_size = 0.5
    ideol_bins, bin_labels = create_ideology_bins(bin_ends = [-3, 5.5], bin_size = ideol_bin_size)

    # Count up each users followers
    follower_dists = followers_data[['user_id', 'follower_ideology']].groupby(['user_id']).apply(lambda x: np.histogram(x, bins = ideol_bins)[0])
    follower_dists = pd.DataFrame(follower_dists.values.tolist(), columns = bin_labels, index = follower_dists.index)
    
    # Normalize into frequency
    norm_follower_dists = follower_dists.div(follower_dists.sum(axis = 1), axis = 0)
    norm_follower_dists = norm_follower_dists.reset_index()
    norm_follower_dists['user_id_str'] = "\"" + norm_follower_dists['user_id'].astype(str) + "\""
    norm_follower_dists['basis'] = "followers"
    norm_follower_dists = norm_follower_dists[ ['user_id', 'user_id_str', 'basis'] + bin_labels ]
    
    # Add in missing tweeters
    missing_tweeters = tweeters[~np.isin(tweeters, norm_follower_dists.user_id)]
    missing_dists = pd.DataFrame([ideol_dist['frequency'].values], columns = bin_labels)
    missing_dists = missing_dists.loc[missing_dists.index.repeat(len(missing_tweeters))]
    missing_dists['user_id'] = missing_tweeters
    missing_dists['basis'] = "population"
    norm_follower_dists = norm_follower_dists.append(missing_dists, ignore_index = True, sort = False)
    
    # Save
    norm_follower_dists.to_csv(ind_dist_file, index = False)
