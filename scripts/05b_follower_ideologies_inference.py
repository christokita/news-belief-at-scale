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
n_batches = 400

# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD

# Path to necessary data
outpath = data_directory + "data_derived/ideological_scores/estimated_ideol_distributions/"
prior_file = outpath + "_population_posterior_samples.csv"
map_estimate_file = outpath + "_population_MAP_estimate.csv"

# Supress theano/pymc3 logging
import logging
logger = logging.getLogger("pymc3")
logger.propagate = False


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
    return pm.distributions.Interpolated(param, x, y)



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

# NOTE: We build our ideology dataset this way because it is faster than loading and filter the pre-made paired tweeter-follower ideology data

####################
# Try inference
####################
# Get our user and the sample of their followers with ideology scores
user_ids, ideology_data = match_followers_to_ideologies(tweeters, follower_ideologies, data_directory, 
                                                        batch = batch, n_batches = n_batches)
del follower_ideologies

# Grab our unique tweeters in this batch that do have follower ideology scores
users_with_scored_followers = ideology_data['user_id'].unique()

# Load our samples from the population-level posterior. This will form our prior.
prior = pd.read_csv(prior_file)
prior_mu = np.array([ prior['mu_samples'] ]).T #needs to be column format for from_posterior function
prior_sigma = np.array([ prior['sigma_samples'] ]).T #needs to be column format for from_posterior function

# Load MAP estimate of population-level ideology distribution
population_MAP_estimate = pd.read_csv(map_estimate_file)

# Unpooled model
n_samples = 500
burn_in = int(n_samples/2)
estimated_ideology_batch = pd.DataFrame(columns = ['user_id', 'user_id_str', 'n_follower_samples', 'mu', 'sigma', 'basis'])   
for user in users_with_scored_followers:
    
    # Grab user's followers. 
    followers = ideology_data[ideology_data.user_id == user]
    n_follower_samples = len(followers.follower_ideology)
    
    # For now, if more than 500k samples, down sample to 500k
    if n_follower_samples > 500000:
        followers = followers.sample(n = 500000, replace = False, random_state = 323)
    follower_samples = np.array([ followers['follower_ideology'] ]).T
    
    with pm.Model() as unpooled_model:
        # Create empirical prior from posterior
        mu = from_posterior('mu', prior_mu, lower_bound = -6, upper_bound = 6)
        sigma = from_posterior('sigma', prior_sigma, lower_bound = 0, upper_bound = 10)
         
        # Likelihood function
        observed_data = pm.Normal('observed_data', mu = mu, sigma = sigma, observed = follower_samples)
        
        # Sample from posterior
        step = pm.NUTS(target_accept = 0.8)
        trace = pm.sample(draws = n_samples, start = None, init = 'advi_map', step = step, random_seed = 323, cores = 2, chains = 2)
 
    # Get estimate of parameters
    posterior_mu = trace.get_values('mu', burn = burn_in, combine = True)
    posterior_sigma = trace.get_values('sigma', burn = burn_in, combine = True)
    
    estimated_ideology_batch = estimated_ideology_batch.append({'user_id': user,
                                                                'user_id_str': "\"" + user + "\"",
                                                                'n_follower_samples': n_follower_samples,
                                                                'mu': np.mean(posterior_mu), 
                                                                'sigma': np.mean(posterior_sigma),
                                                                'basis': "followers"}, ignore_index = True)
        
    
####################
# Use population-level estimates for users we don't have follower samples for
####################
# For tweeters that didn't have followers, use population-level estimate
missing_users = [x for x in user_ids if x not in users_with_scored_followers]
missing_user_estimates = pd.DataFrame({'user_id': missing_users,
                                      'user_id_str': ["\"" + x + "\"" for x in missing_users]}) 
missing_user_estimates['mu'] = population_MAP_estimate['mu'].iloc[0]
missing_user_estimates['sigma'] = population_MAP_estimate['sigma'].iloc[0]
missing_user_estimates['n_follower_samples'] = 0
missing_user_estimates['basis'] = "population"

# Append to all estiamtes and save to temporary folder to hold partial results
estimated_ideology_batch = estimated_ideology_batch.append(missing_user_estimates, ignore_index = True)
estimated_ideology_batch['batch'] = batch
temp_dir = outpath + 'TEMP_follower_ideology_dist_shapes/'
os.makedirs(temp_dir, exist_ok = True)
estimated_ideology_batch.to_csv(temp_dir + 'batch_' + str(batch).zfill(3) + '.csv', index = False)


####################
# Combine all data that was processed in parallel
####################
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

    
####################
# Uncomment to see individual fits
####################
#select_est_parameters = select_user = 0
#
#import matplotlib.pyplot as plt
#x = np.linspace(-4, 4, 100)
#plt.plot(x, stats.norm(loc = est_parameters['mu'], scale = est_parameters['sigma']).pdf(x))    
#plt.hist(followers.follower_ideology, density = True, bins = np.arange(-5, 5, 0.25))