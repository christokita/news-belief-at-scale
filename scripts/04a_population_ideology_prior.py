#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: `04a_population_ideology_prior.py`
Date: August 28, 2020
Author: Chris Tokita
Purpose: Conduct MCMC bayesian inference of population-level distribution of follower ideology.
Details:
    (Copies of data are currently stored on external harddrive and high-performance cluster.)
    The results of this inference will be used as the prior for the inference of tweeter's follower ideology.

Data In: CSV files of tweets and follower ideologies.
    `<data storage location>/data_derived/tweets/tweets_labeled.csv`
    `<data storage location>/data_derived/ideological_scores/cleaned_followers_ideology_scores.csv`

Data Out: CSV files of posterior samples and of the point estimate for Twitter population ideology.
    `<data storage location>/data_derived/ideological_scores/estimated_ideol_distributions/`

        `_population_posterior_samples.csv`
            `mu_samples`:    posterior samples for population ideology mean.
            `sigma_samples`: posterior samples for pouplation ideology standard deviation.

        `_population_MAP_estimate.csv`
            `mu`:    MAP estimate of mean ideology of population of Twitter users.
            `sigma`: MAP estimate of standard deviation of ideology of population of Twitter users.

Machine: High-performance computing cluster
    This script is batched to the cluster using `slurm_scripts/compute_population_ideology_dist.cmd`
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import pymc3 as pm


# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD

# Paths for ideology distrubtion data
outpath = data_directory + "data_derived/ideological_scores/estimated_ideol_distributions/"
posterior_samples_file = outpath + "_population_posterior_samples.csv"
population_map_file = outpath + "_population_MAP_estimate.csv"


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


# Main script
if __name__ == '__main__':
    
    ####################
    # Create list of unique follower/user ideology scores
    ####################
    # Load follower ideology data 
    follower_ideologies = pd.read_csv(data_directory + "data_derived/ideological_scores/cleaned_followers_ideology_scores.csv",
                                      dtype = {'user_id': object, 'pablo_score': float})
    follower_ideologies = follower_ideologies.drop(columns = ['accounts_followed'])
    follower_ideologies = follower_ideologies.rename(columns = {'user_id': 'follower_id'})
    
    # Add tweeter ideologies (since they can also be followers) 
    tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                 dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                          'user_id': object, 'tweet_id': object, 
                          'retweeted_user_id': object, 'retweet_id': object,
                          'quoted_user_id': object, 'quoted_id': object})
    tweeter_ideologies = tweets[['user_id', 'user_ideology']]
    tweeter_ideologies = tweeter_ideologies.rename(columns = {'user_id': 'follower_id', 'user_ideology': 'pablo_score'})
    tweeter_ideologies = tweeter_ideologies[~pd.isna(tweeter_ideologies['pablo_score'])]
    tweeter_ideologies = tweeter_ideologies.drop_duplicates()
    follower_ideologies = follower_ideologies.append(tweeter_ideologies, ignore_index = True)
    follower_ideologies = follower_ideologies.drop_duplicates()
    del tweeter_ideologies
    
    # Rename columns
    follower_ideologies = follower_ideologies.rename(columns={"pablo_score": "follower_ideology"})
        
    
    ####################
    # Infer population level ideology distrubiton shape
    ####################
    # Prepare data for inference
    # samples = follower_ideologies['follower_ideology'].sample(frac = 0.1, replace = False, random_state = 323) #downsample due to sheer size of data
    samples = follower_ideologies['follower_ideology'].sample(n = 500000, replace = False, random_state = 323) #downsample due to sheer size of data
    del follower_ideologies
    sample_mu = np.array([np.mean(samples)])
    sample_sigma = np.array([np.std(samples)])
    
    # Run inference using NUTS
    n_cores = 2
    total_samples = 3000
    burn_in = 1000
    niter = int( (total_samples + n_cores*burn_in) / n_cores ) #account for removal of burn in at the beginning of each chain
    with pm.Model() as model:
        # define priors
        # Based on P. Barbera's work we know ideology scores tend to be normally distributed around 0
        # For standard deviation, typically an exponential distribution is used
        mu = pm.Normal('mu', mu = 0, sigma = 2, shape = sample_mu.shape)
        sigma = pm.Exponential('sigma', lam = 2, shape = sample_sigma.shape)
        
        # define likelihood
        observed_data = pm.Normal('observed_data', mu = mu, sigma = sigma, observed = samples)
        
        # inference
        map_estimate = pm.find_MAP()
        step = pm.NUTS(target_accept = 0.90)
        trace = pm.sample(draws = niter, start = None, init = 'advi_map', step = step, random_seed = 323, cores = n_cores, chains = n_cores)
     
    print("Done with inference!")    
        
    # Get samples of population-level posterior for use as prior in individual-level inference later
    del samples
    posterior_mu = trace.get_values('mu', burn = burn_in, combine = True)
    posterior_sigma = trace.get_values('sigma', burn = burn_in, combine = True)
    posterior_samples = pd.DataFrame({'mu_samples': posterior_mu.flatten(), 
                                      'sigma_samples': posterior_sigma.flatten()})
    posterior_samples.to_csv(posterior_samples_file, index = False)
    
    # Get MAP estimate
    population_MAPestimate = pd.DataFrame({'mu': map_estimate['mu'], 'sigma': map_estimate['sigma']})
    population_MAPestimate.to_csv(population_map_file, index = False)
    
#    import matplotlib.pyplot as plt
#    import scipy.stats as stats
#    plt.hist(samples, density = True)
#    x = np.linspace(-5, 5, 100)
#    plt.plot(x, stats.norm(loc = map_estimate['mu'], scale = map_estimate['sigma']).pdf(x))
#    plt.plot(x, stats.norm(loc = np.mean(posterior_mu), scale = np.mean(posterior_sigma)).pdf(x))