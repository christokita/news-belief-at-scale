#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 14:36:31 2020

@author: ChrisTokita

SCRIPT
Estimate the number of people who believed the articles they were exposed to on Twitter.
"""

####################
# Load packages and data, set important parameters
####################
import pickle
import pandas as pd
import numpy as np
import copy
import re
import multiprocessing as mp

# high level directory (external HD or cluster storage)
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD
#data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
outpath = data_directory + "data_derived/belief/"

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

def format_belief_data(ideological_bins, story_id, belief_data):
    """
    This will take the estimated exposed to an tweet and will estimate how many believed it (by ideology)
    """
    # Grab belief rates for this story
    story_belief = belief_data[belief_data.total_article_number == story_id]
    
    # Create dataframe to story belief rate data
    belief_rate_df = pd.DataFrame(columns = ["bin", "upper_bin_edge", "ideology", "True", "Unsure", "False"])
    
    for j in np.arange(ideological_bins.shape[0]):
        
        if ideological_bins['upper_bin_edge'].iloc[j] <= -2.5:
            individual_type = story_belief[story_belief.ideology_score == -3]
        if -2.5 < ideological_bins['upper_bin_edge'].iloc[j] <= -1.5:
            individual_type = story_belief[story_belief.ideology_score == -2]
        if -1.5 < ideological_bins['upper_bin_edge'].iloc[j] <= -0.5:
            individual_type = story_belief[story_belief.ideology_score == -1]
        if -0.5 < ideological_bins['upper_bin_edge'].iloc[j] <= 0.5:
            individual_type = story_belief[story_belief.ideology_score == 0]
        if 0.5 < ideological_bins['upper_bin_edge'].iloc[j] <= 1.5:
            individual_type = story_belief[story_belief.ideology_score == 1]
        if 1.5 < ideological_bins['upper_bin_edge'].iloc[j] <= 2.5:
            individual_type = story_belief[story_belief.ideology_score == 2]
        if ideological_bins['upper_bin_edge'].iloc[j] > 2.5:
            individual_type = story_belief[story_belief.ideology_score == 3]
            
        new_row = pd.DataFrame({"bin": ideological_bins['bin_label'].iloc[j], 
                                "upper_bin_edge": ideological_bins['upper_bin_edge'].iloc[j],
                                "ideology": individual_type['ideology'].iloc[0], 
                                "True": individual_type.belief_freq[individual_type.belief == "True"].iloc[0], 
                                "Unsure": individual_type.belief_freq[individual_type.belief == "Unsure"].iloc[0], 
                                "False": individual_type.belief_freq[individual_type.belief == "False"].iloc[0]}, 
                            index = [0])
        belief_rate_df = belief_rate_df.append(new_row, ignore_index = True)
    
    return belief_rate_df


def estimate_belief(estimated_exposure_df):
    """
    This function will estimate how many people potentially believed the article.
    """
    # Get ideological bins
    ideol_bins = get_bin_edges(estimated_exposure_df.columns.values)
    
    # Loop through tweeters in exposure dataset
    estimated_belief_df = copy.deepcopy(estimated_exposure_df)
    for i in np.arange(estimated_exposure_df.shape[0]):
        
        # Format belief data
        story = estimated_exposure_df.total_article_number.iloc[i]
        belief_rates = format_belief_data(ideological_bins = ideol_bins, 
                                          story_id = story, 
                                          belief_data = belief_data)
        
        # Estimate belief among those exposed to this tweet/article 
        # NOTE: For now just doing true
        est_belief = belief_rates['True'].values * estimated_exposure_df.loc[estimated_exposure_df.index[i], ideol_bins['bin_label']].values
        estimated_belief_df.loc[estimated_belief_df.index[i], ideol_bins['bin_label']] = np.round(est_belief).astype(int)
        
    return estimated_belief_df


if __name__ == '__main__':
    
    ####################
    # Load and shape belief data
    ####################
    # Load belief data, turn into dataframe
    belief_file = data_directory + '/data/article_belief/response_distribution.p'
    belief_data = pickle.load( open(belief_file, "rb") )
    belief_data = pd.DataFrame.from_dict(belief_data, orient = 'index')
    belief_data['total_article_number'] = belief_data.index
    
    # Drop unneccessary columns 
    belief_data = belief_data.drop(columns = ['date', 'art_type'])
    belief_data = belief_data.loc[:, ~belief_data.columns.str.match('.*response_count')] #drop response count columns
    
    # Melt and form new columns for ideology and belief
    belief_data = pd.melt(belief_data, id_vars = ['total_article_number', 'fc_rating', 'mean_fc_likert', 'total_responses'])
    belief_data['belief'] = belief_data['variable'].str.extract(r'_([a-z])$')
    belief_data['ideology'] = belief_data['variable'].str.extract(r'^(.*)_[a-z]$')
    belief_data['ideology'] = belief_data['ideology'].str.replace(': Middle of the road', '') #remove extra text from moderates
    belief_data = belief_data[belief_data.ideology != "Haven't thought much about it"]
    belief_data['ideology_score'] = belief_data['ideology'].map({'Extremely Liberal': -3, 'Liberal': -2, 'Slightly Liberal': -1,
                                                                 'Moderate':0, 
                                                                 'Slightly Conservative': 1, 'Conservative': 2, 'Extremely Conservative': 3})
    
    # Clean up dataframe
    belief_data = belief_data.drop(columns = ['variable'])
    belief_data = belief_data.rename(columns = {'value': 'belief_freq'})
    belief_data['belief'] = belief_data['belief'].map({'t': 'True', 'f': 'False', 'c': 'Unsure'})
    
    
    ####################
    # Load estimated exposure data
    ####################
    estimated_exposure = pd.read_csv(data_directory + "data_derived/exposure/estimated_users_exposed_over_time.csv", dtype = {'user_id': 'int64', 'tweet_id': 'int64'})
    estimated_exposure = estimated_exposure[estimated_exposure.total_article_number > 10] #remove first 10 articles
    
    
    ####################
    # Compute belief in parallel
    ####################
    # Prepare for parallel computation
    n_cores = mp.cpu_count()
    batched_data = np.array_split(estimated_exposure, n_cores)
    pool = mp.Pool(mp.cpu_count())
    result_list = pool.map(estimate_belief, batched_data)
    pool.close()
    pool.join()

    # Write to csv
    estimated_belief_data = pd.concat(result_list)
    estimated_belief_data = estimated_belief_data.sort_values(by = ['total_article_number', 'tweet_number'])
    estimated_belief_data.to_csv(outpath + "estimated_belief_over_time.csv", index = False)
    
    
    