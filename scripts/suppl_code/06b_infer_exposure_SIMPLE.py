#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 09:36:18 2020

@author: ChrisTokita

SCRIPTS
Exposure over time using inferred distribution of followers. We are determining total exposure.
"""

####################
# Load libraries and packages, set important parameters, load data
####################
import pandas as pd
import numpy as np
import scipy.stats as stats
import re
import multiprocessing as mp

# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD

####################
# Functions for analysis
####################
def update_exposure_data(exposure_df):
    # Get ideological bins
    bin_labels = [col for col in exposure_df if re.search("ideol_[-0-9.]+_[-0-9.]+", col)]
    
    # Loop through tweeters in exposure dataset
    for i in np.arange(exposure_df.shape[0]):
        
        # If user has no followers/new followers, skip
        if exposure_df['new_exposed_users'].iloc[i] == 0:
            continue
        
        # Get user info
        user_id = exposure_df['user_id'].iloc[i]
        follower_dist = follower_distributions[follower_distributions.user_id == user_id]
        follower_dist = follower_dist[bin_labels]
        
        # Determine how many users we don't have ideological scores for
        new_exposed_users = exposure_df['new_exposed_users'].iloc[i]
        accounted_new_exposed = sum(exposure_df[bin_labels].iloc[i])
        unaccounted_new_exposed = new_exposed_users - accounted_new_exposed
        
        # Infer 
        inferred_exposed = np.round(unaccounted_new_exposed * follower_dist.iloc[0])
        
        # Add to the concrete counts of binned ideology
        total_estimated_exposure = exposure_df[bin_labels].iloc[i] + inferred_exposed
        exposure_df.loc[exposure_df.index[i], bin_labels] = total_estimated_exposure
    return(exposure_df)



####################
# Loop through our exposure data and infer additional users that were exposed
####################
if __name__ == '__main__':
    
    # Load exposure timeseries data
    exposure_data = pd.read_csv(data_directory + "data_derived/exposure/users_exposed_over_time.csv", 
                                dtype = {'user_id': object, 'tweet_id': object})
    exposure_data = exposure_data.fillna(0)
    
    # Load follower distributions of each tweeter
    follower_distributions = pd.read_csv(data_directory + "data_derived/ideological_scores/estimated_ideol_distributions/follower_binned_ideology_distributions.csv", 
                                         dtype = {'user_id': object})

    # Prepare for parallel computation
    n_cores = mp.cpu_count()
    batched_data = np.array_split(exposure_data, n_cores)
    pool = mp.Pool(mp.cpu_count())
    result_list = pool.map(update_exposure_data, batched_data)
    pool.close()
    pool.join()

    # Write to csv
    est_exposure_data = pd.concat(result_list)
    est_exposure_data = est_exposure_data.sort_values(by = ['total_article_number', 'tweet_number'])
    est_exposure_data.to_csv(data_directory + "data_derived/exposure/estimated_users_exposed_over_time_SIMPLE.csv", index = False)
