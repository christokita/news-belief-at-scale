#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:13:26 2020

@author: ChrisTokita

SCRIPT
Conduct bayesian inference on each user's distrubtion of follower ideology.

NOTE
This script is intended to run on an HPC cluster with SLURM scheduling that will use a job array.
It will grab the job number in the array of 200 jobs, will divide the dataset up accordingly, and will grab the corresponding chunk.
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import os
import math
import re
import sys
import time

# Get batch number
batch = int(sys.argv[1]) 
n_batches = 200

# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD
#
# Path to necessary data
outpath = data_directory + "data_derived/ideological_scores/estimated_ideol_distributions/"
prior_file = outpath + "_population_estimate.csv"


####################
# Functions for estimating follower ideology 
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
def match_followers_to_ideologies(user_ids, ideologies, data_directory, batch, n_batches):
    """
    Function that will take a set of user IDs and subset a batch of them,
    then match these users to their followers that have ideology scores
    
    OUTPUT
    - followers_matched: data set containing all users in batch and their followers with ideology scores.
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
    return followers_matched

# Likelihood function for s.d.
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

# Function to estimate ideology of followers
def estimate_ideology_distribution(user_id, ideology_samples, sigmas, log_prior):
    """
    This function will estimate the distribution of ideologies of a user.
    It will assume the distrubtion follows the mean of the samples and will esetimate the s.d.
    
    INPUT:
    - user_id: id of user who's followers we are estimating.
    - ideology_samples: ideology scores for set of followers.
    - sigmas: list of sigma (s.d.) values to calculate likelihood.
    - prior: assumption of log likelihood over same sigmas given to function.
    """
    # Calculate log likelihood
    est_mu = np.mean(ideology_samples)
    log_L = np.array([])
    for sigma in sigmas:
        val = logL_sd(mu = est_mu, 
                      sigma = sigma, 
                      samples = ideology_samples)
        log_L = np.append(log_L, val)
        
    # Calculate posterior, our population level posterior is our prior here
    log_posterior = log_prior + log_L
    
    # Grab maximum a posteriori estimation, return data row
    est_sig = sigmas[log_posterior == max(log_posterior)]
    follower_estimate = pd.DataFrame({'user_id': user_id, 
                                      'user_id_str': "\"" + user_id + "\"", 
                                      'mu': est_mu, 
                                      'sigma': est_sig})
    return follower_estimate


####################
# Load data
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
follower_ideologies = follower_ideologies.rename(columns = {'user_id': 'follower_id', 'pablo_score': 'follower_ideology'})

# Unique tweeters to loop through
tweeters = tweets['user_id'].unique()

# Add tweeter ideologies (since they can also be followers) 
# This dataset represents all followers for which we have ideology scores.
tweeter_ideologies = tweets[['user_id', 'user_ideology']].drop_duplicates()
tweeter_ideologies = tweeter_ideologies.rename(columns = {'user_id': 'follower_id', 'user_ideology': 'follower_ideology'})
tweeter_ideologies = tweeter_ideologies[~pd.isna(tweeter_ideologies['follower_ideology'])]
follower_ideologies = follower_ideologies.append(tweeter_ideologies, ignore_index = True)
follower_ideologies = follower_ideologies.drop_duplicates()
del tweeter_ideologies, tweets

# We will use ideological category bins of size 1, with moderate being [-0.5, 0.5]
# But per SMAPP instructions, we will not use the normalized scores.
# Normalize ideologies between -3.5 and 3.5.
min_ideol = follower_ideologies['follower_ideology'].min()
max_ideol = follower_ideologies['follower_ideology'].max()
follower_ideologies['norm_ideology'] = (follower_ideologies['follower_ideology'] - min_ideol) / (max_ideol - min_ideol)
follower_ideologies['norm_ideology'] = 3.5 * (2*follower_ideologies['norm_ideology'] - 1)
del min_ideol, max_ideol

# Load our prior, which is the population-level posterior
log_pop_posterior = pd.read_csv(prior_file)
    

####################
# Estimate distribution of followers for each unique tweeter in this batch
####################
# Take our batch of user IDs and match them to their followers that have ideology scores
ideology_data = match_followers_to_ideologies(user_ids = tweeters, 
                                              ideologies = follower_ideologies, 
                                              data_directory = data_directory, 
                                              batch = batch, 
                                              n_batches = n_batches)
del follower_ideologies, tweeters

# Grab our unique tweeters in this batch
user_ids = ideology_data['user_id'].unique()

# Set the sigmas we will be checking 
sigmas = np.arange(0.1, 4.1, 0.1)

# Loop through our tweeters and estimate ideology distribution of followers
estimated_ideology_batch = pd.DataFrame(columns = ['user_id', 'user_id_str', 'mu', 'sigma'])
for user_id in user_ids:
    followers = ideology_data.follower_ideology[ideology_data.user_id == user_id]
    follower_estimate = estimate_ideology_distribution(user_id = user_id, 
                                                       ideology_samples = followers, 
                                                       sigmas = sigmas, 
                                                       log_prior = log_pop_posterior['log_pr'])
    estimated_ideology_batch = estimated_ideology_batch.append(follower_estimate, ignore_index = True)


# Create temporary folder to hold partial results
temp_dir = outpath + 'TEMP_follower_ideology_dist_shapes/'
os.makedirs(temp_dir, exist_ok = True)
estimated_ideology_batch.to_csv(temp_dir + 'batch_' + str(batch).zfill(3) + '.csv', index = False)

# If this is the last batch to finish, bind all temoprary files together and save
time.sleep(10)
temp_files = os.listdir(temp_dir)
if len(temp_files) == n_batches:
    estimated_ideology_distributions = pd.DataFrame(columns = estimated_ideology_batch.columns)
    for file in temp_files:
        batch_data = pd.read_csv(temp_dir + file, dtype = {'user_id': str})
        estimated_ideology_distributions = estimated_ideology_distributions.append(batch_data)
        os.remove(temp_dir + file)
    estimated_ideology_distributions.to_csv(outpath + 'follower_ideology_distribution_shapes.csv', index = False)
    os.rmdir(temp_dir)