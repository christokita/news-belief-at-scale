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

# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD

# Paths for ideology distrubtion data
outpath = data_directory + "data_derived/ideological_scores/estimated_ideol_distributions/"
prior_file = outpath + "_population_estimate.csv"


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

# Likelihood function for estimating s.d. of ideology
def logL_sd(mu, sigma, samples):
    """
    Function that calculates the likelihood of the n observations, given a normal distribution.
    
    INPUTS
    - n: number of observations
    - V: sample variance of observations
    - sigma: s.d. of normal distrubtion we assume the observations are coming from
    """
    n = len(samples)
    first_part = -n/2 * np.log(2*math.pi)
    second_part = -n/2 * np.log(sigma**2)
    third_part = -1/(2 * sigma**2) * np.sum( (samples-mu)**2 )
    logL = first_part + second_part + third_part
    return logL

# Main script
if __name__ == '__main__':
    

    ####################
    # Load tweeters, followers, and ideology scores
    ####################
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


    ####################
    # Create data set matching follower ideology to each tweeter
    ####################
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
    followers_data = followers_data[['follower_id', 'follower_ideology']].drop_duplicates()
    
    # Get population sample mean
    mean_population_ideol = followers_data['follower_ideology'].mean() #we assume the distribution is centered on this
    
    # Set up uninformative prior
    sigmas = np.arange(0.1, 4.1, 0.1)
    log_prior = np.log(1 / sigmas**2)    
    
    # Calculate log likelihood
    log_L = np.array([])
    for sigma in sigmas:
        val = logL_sd(mu = mean_population_ideol, 
                      sigma = sigma, 
                      samples = followers_data['follower_ideology'])
        log_L = np.append(log_L, val)
        
    # Calculate posterior and save for future use
    log_post = log_prior + log_L
    log_pop_posterior = pd.DataFrame({'mu': mean_population_ideol, 'sigma': sigmas, 'log_pr': log_post})
    log_pop_posterior.to_csv(prior_file, index = False)