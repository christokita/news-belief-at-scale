#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: `05b_infer_exposure_ideology.py`
Date: Aug 24, 2020
Author: Chris Tokita
Purpose: Estimate what user ideologies were exposed over time to the news articles.
Details:
    (Copies of data are currently stored on external harddrive and high-performance cluster.)
 
Data In: CSV files the compiled exposure over time dataframe and the known Twitter user ideologies.
    `<data storage location>/data_derived/exposure/users_exposed_over_time.csv`
    `<data storage location>/data_derived/ideological_scores/cleaned_followers_ideology_scores.csv`

Data Out: CSV file of exposure over time dataframe, including information about the user ideologies exposed.
    `<data storage location>/data_derived/exposure/estimated_users_exposed_over_time.csv`

Machine: High-performance computing cluster
    This script is batched to the cluster using `slurm_scripts/estimate_exposure_ideology.cmd`
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
def get_bin_edges(df_columns):
    """
    This funciton will take the array of the exposure dataset's column names and extract the ideological bin edges
    """
    
    bin_labels = [col for col in df_columns if re.search("ideol_[-0-9.]+_[-0-9.]+", col)]
    bin_edges = [label.replace("ideol_", "") for label in bin_labels]
    bin_edges = np.array([string.split("_") for string in bin_edges])
    bin_edges = bin_edges.flatten().astype(float)
    bin_edges = np.unique(bin_edges)
    bin_edges = np.delete(bin_edges, 0)
    return pd.DataFrame({"bin_label": bin_labels, "upper_bin_edge": bin_edges})

def calculate_extra_exposed(total_users, mu, sigma, upper_bin_edges):
    """
    This function takes the follower distriubtion shape and calculates how many users were exposed (binning into the ideological bins)
    """
    
    bin_size = upper_bin_edges[1] - upper_bin_edges[0]  #assuming regularly sized bins
    follower_dist = stats.norm(loc = mu, scale = sigma)
    exposed_probs = np.array([])
    for upper_edge in upper_bin_edges:
        # If it's the lowest bin, calculate space of -inf to upper edge
        if upper_edge == min(upper_bin_edges):
            prob = follower_dist.cdf(upper_edge)
        # If it's the upper most bin, calculate lower edge to inf
        elif upper_edge == max(upper_bin_edges):
            prob = 1 - follower_dist.cdf(upper_edge - bin_size)
        # Otherwise integrate between upper and lower edges
        else:
            prob = follower_dist.cdf(upper_edge) - follower_dist.cdf(upper_edge - bin_size)
        exposed_probs = np.append(exposed_probs, prob)
    # Round to nearest full individual and return
    exposed_counts = total_users * exposed_probs
    exposed_counts = np.round(exposed_counts)
    return exposed_counts

def update_exposure_data(exposure_df):
    # Get ideological bins
    ideol_bins = get_bin_edges(exposure_df.columns.values)
    
    # Loop through tweeters in exposure dataset
    for i in np.arange(exposure_df.shape[0]):
        
        # If user has no followers/new followers, skip
        if exposure_df['new_exposed_users'].iloc[i] == 0:
            continue
        
        # Get user info
        user_id = exposure_df['user_id'].iloc[i]
        
        # Determine how many users we don't have ideological scores for
        new_exposed_users = exposure_df['new_exposed_users'].iloc[i]
        accounted_new_exposed = sum(exposure_df[ideol_bins['bin_label']].iloc[i])
        unaccounted_new_exposed = new_exposed_users - accounted_new_exposed
        
        # Get distribution shape
        follower_mu = follower_distributions.mu[follower_distributions.user_id == user_id].iloc[0]
        follower_sigma = follower_distributions.sigma[follower_distributions.user_id == user_id].iloc[0]
        
        # Infer 
        inferred_exposed = calculate_extra_exposed(total_users = unaccounted_new_exposed,
                                                   mu = follower_mu,
                                                   sigma = follower_sigma,
                                                   upper_bin_edges = ideol_bins['upper_bin_edge'])
        
        # Add to the concrete counts of binned ideology
        total_estimated_exposure = exposure_df[ideol_bins['bin_label']].iloc[i] + inferred_exposed
        exposure_df.loc[exposure_df.index[i], ideol_bins['bin_label']] = total_estimated_exposure
    return(exposure_df)



####################
# Loop through our exposure data and infer additional users that were exposed
####################
if __name__ == '__main__':
    
    # Load exposure timeseries data
    exposure_data = pd.read_csv(data_directory + "data_derived/exposure/users_exposed_over_time.csv", 
                                dtype = {'user_id': object, 'tweet_id': object})
    
    # Load inferred follower distribution shape for each tweeter
    follower_distributions = pd.read_csv(data_directory + "data_derived/ideological_scores/estimated_ideol_distributions/follower_ideology_distribution_shapes.csv", dtype = {'user_id': object})

    # Prepare for parallel computation
    n_cores = mp.cpu_count()
    batched_data = np.array_split(exposure_data, n_cores)
    pool = mp.Pool(mp.cpu_count())
    result_list = pool.map(update_exposure_data, batched_data)
    pool.close()
    pool.join()

    # Write to csv
    updated_exposure_data = pd.concat(result_list)
    updated_exposure_data = updated_exposure_data.sort_values(by = ['total_article_number', 'tweet_number'])
    updated_exposure_data.to_csv(data_directory + "data_derived/exposure/estimated_users_exposed_over_time.csv", index = False)
