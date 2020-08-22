#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:13:26 2020

@author: ChrisTokita

SCRIPT
Conduct bayesian inference on each user's distrubtion of follower ideology
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import os
#import scipy.stats as stats
#import matplotlib.pyplot as plt
import math
import multiprocessing as mp


# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD
outpath = data_directory + "data_derived/ideological_scores/estimated_ideol_distributions/"




####################
# Functions for estimating follower ideology 
####################
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

## Likelihood function for s.d.
#def logL_sd(mu, sigma, samples):
#    """
#    Function that calculates the likelihood of the n observations, given sigma of distribution.
#    
#    INPUTS
#    - n: number of observations
#    - V: sample variance of observations
#    - sigma: s.d. of normal distrubtion we assume the observations are coming from
#    """
#    a, b = (-3.5 - mu) / sigma, (3.5 - mu) / sigma
#    probs = stats.truncnorm.pdf(samples, a, b, loc = mu, scale = sigma)
#    probs[probs == 0] = 10**-30
#    logL = np.sum( np.log(probs) )
#    return logL

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
    return_data = pd.DataFrame({'user_id': user_id, 
                                'user_id_str': "\"" + user_id + "\"", 
                                'mu': est_mu, 
                                'sigma': est_sig})
    return return_data
    


####################
# Calculate distribution of ideologies, partially in parallel
####################
if __name__ == "__main__":
    
    ####################
    # Load data
    ####################
    # Combine follower and tweeter ideology data
    ideology_data = pd.read_csv(data_directory + "data_derived/ideological_scores/paired_tweeter-follower_ideology.csv",
                                dtype = {'user_id': object})
    
    # We will use ideological category bins of size 1, with moderate being [-0.5, 0.5]
    # But per SMAPP instructions, we will not use the normalized scores.
    # Normalize ideologies between -3.5 and 3.5.
    min_ideol = ideology_data['follower_ideology'].min()
    max_ideol = ideology_data['follower_ideology'].max()
    ideology_data['norm_ideology'] = (ideology_data['follower_ideology'] - min_ideol) / (max_ideol - min_ideol)
    ideology_data['norm_ideology'] = 3.5 * (2*ideology_data['norm_ideology'] - 1)
    del min_ideol, max_ideol
    
    ####################
    # Create prior for standard deviation based on population of scores we have
    ####################
    # Check if this estiamte already exists, if so load it... 
    prior_file = outpath + "_population_estimate.csv"
    if os.path.isfile(prior_file):
        log_pop_posterior = pd.read_csv(prior_file)
        
    # ...otherwise calculate the population estimate of follower ideology  
    if not os.path.isfile(prior_file):
        
        # Get population sample mean
        mean_population_ideol = ideology_data['follower_ideology'].mean() #we assume the distribution is centered on this
        
        # Set up uninformative prior
        sigmas = np.arange(0.1, 4.1, 0.1)
        log_prior = np.log(1 / sigmas**2)    
        
        # Calculate log likelihood
        log_L = np.array([])
        for sigma in sigmas:
            val = logL_sd(mu = mean_population_ideol, 
                          sigma = sigma, 
                          samples = ideology_data['follower_ideology'])
            log_L = np.append(log_L, val)
            print("Done: calculated log L for sigma = " + str(sigma))
            
        # Calculate posterior and save for future use
        log_post = log_prior + log_L
        log_pop_posterior = pd.DataFrame({'sigma': sigmas, 'log_pr': log_post})
        log_pop_posterior.to_csv(prior_file, index = False)
    
    
    #MAP_sigma = log_pop_posterior.sigma[log_pop_posterior.log_pr == max(log_pop_posterior['log_pr'])]
    #x = np.linspace(-5, 5, 100)
    #plt.plot(x, stats.norm.pdf(x, mean_population_ideol, MAP_sigma))
    #plt.hist(ideology_data.follower_ideology, density = True, bins = np.arange(-5.5, 5.5, 1))

    
    ####################
    # Estimate distribution of followers for each tweeter
    ####################
    # Function to collect results from multi-threading 
    results = []
    def collect_results(result):
        """Uses apply_async's callback to setup up a separate Queue for each process"""
        results.append(result)
    
    # Grab our unique tweeters
    tweeters = ideology_data['user_id'].unique()
    
    # Set the sigmas we will be checking 
    sigmas = np.arange(0.1, 4.1, 0.1)
    
    # Loop through our tweeters and estimate ideology distribution of followers
    pool = mp.Pool(mp.cpu_count())
    for user_id in tweeters:
        followers = ideology_data.follower_ideology[ideology_data.user_id == user_id]
        pool.apply_async(estimate_ideology_distribution, args = (user_id, followers, sigmas, log_pop_posterior['log_pr'], ), callback = collect_results) 
    pool.close()
    pool.join()

    # Get results and save
    estimated_ideologies = pd.concat(results).reset_index(drop = True)
    estimated_ideologies.to_csv(outpath + 'follower_ideology_dist_shapes.csv', index = False)