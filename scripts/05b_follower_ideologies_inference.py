#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 10:41:23 2020

@author: ChrisTokita
"""


####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import scipy.stats as stats
import os
import math
import re
import sys
import time
import pymc3 as pm


# Get batch number
batch = int(sys.argv[1]) 
n_batches = 200

# high level directory (external HD or cluster storage)
#data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD

# Path to necessary data
outpath = data_directory + "data_derived/ideological_scores/estimated_ideol_distributions/"
prior_file = outpath + "_population_posterior_samples.csv"


####################
# Functions for estimating follower ideology 
####################
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
    return user_id_batch, followers_matched

def from_posterior(param, samples, lower_bound, upper_bound):
    """
    Function to sample from the population-level posterior to form an empirical prior for future use at the individual level.
    
    Modified from: https://docs.pymc.io/notebooks/updating_priors.html
    """
    smin, smax = np.min(samples), np.max(samples)
    x = np.linspace(smin, smax, 1000)
    y = stats.gaussian_kde(samples.ravel())(x)

    # what was never sampled should have a small probability but not 0,
    # so we'll extend the domain and use linear approximation of density on it
    x = np.concatenate([[lower_bound], x, [upper_bound]])
    y = np.concatenate([[0], y, [0]])
    return Interpolated(param, x, y)



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


####################
# Try inference
####################
# Get our user and the sample of their followers with ideology scores
user_ids, ideology_data = match_followers_to_ideologies(tweeters, follower_ideologies, data_directory, 
                                                        batch = batch, n_batches = n_batches)

# Grab our unique tweeters in this batch that do have follower ideology scores
users_with_scored_followers = ideology_data['user_id'].unique()

# Load our samples from the population-level posterior. This will form our prior.
prior = pd.read_csv(prior_file)
prior_mu = np.array([ prior['mu_samples'] ]).T #needs to be column format for from_posterior function
prior_sigma = np.array([ prior['sigma_samples'] ]).T #needs to be column format for from_posterior function

# Loop through our tweeters and estimate ideology distribution of followers
estimated_ideology_batch = pd.DataFrame(columns = ['user_id', 'user_id_str', 'mu', 'sigma'])
for user_id in users_with_scored_followers[60:61]:
    
    # Grab follower ideologies we have for this user
    followers = ideology_data.follower_ideology[ideology_data.user_id == user_id]
    
    # Determine MAP estimate 
    with pm.Model() as model:
        # Create empirical prior from posterior
        mu = from_posterior('mu', prior_mu, lower_bound = -5, upper_bound = 5)
        sigma = from_posterior('sigma', prior_sigma, lower_bound = 0, upper_bound = 5)
        alpha = from_posterior('sigma', prior_alpha, lower_bound = -10, upper_bound = 10)
        
        # Likelihood function
        observed_data = pm.Normal('observed_data', mu = mu, sigma = sigma, observed = followers)
        
        # Best estimate of distribution shape
        est_parameters = pm.find_MAP()
        
    follower_estimate = pd.DataFrame({'user_id': user_id, 
                                      'user_id_str': "\"" + user_id + "\"", 
                                      'mu': est_parameters['mu'], 
                                      'sigma': est_parameters['sigma']}, index = [0])
    estimated_ideology_batch = estimated_ideology_batch.append(follower_estimate, ignore_index = True)

estimated_ideology_batch['mu_basis'] = "followers"

    

    

x = np.linspace(-4, 4, 100)
plt.plot(x, stats.norm(loc = est_parameters['mu'], scale = est_parameters['sigma']).pdf(x))    
plt.hist(followers, density = True)