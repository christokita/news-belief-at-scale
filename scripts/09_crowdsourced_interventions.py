#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 09:28:23 2021

@author: ChrisTokita

SCRIPT:
Simulate crowd-sourced fact-checking and intervention
"""

####################
# Load packages 
####################
import numpy as np
import pandas as pd
import re
import os
import sys

# Load custom intervention functions from 
model_interventions = __import__("08a_model_interventions")
simulate_intervention = model_interventions.simulate_intervention
match_followers_to_tweet = model_interventions.match_followers_to_tweet


####################
# Parameters for simulation
####################
# Important paths
# data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
crowd_files = ['crowdsdata_fringe.csv', 'crowdsdata_mainstream.csv'] #add in more simulated crowd data as we wish
outpath = data_directory + "data_derived/interventions/"

# Crowd-sourcing paramters
which_story = int(sys.argv[1]) #get from command line
n_replicates = 30
crowd_size = 30 #number of individuals asked to evaluate article
crowd_type = "random" #type of crowd composition ("random" or "balanced", i.e., ideologically balanced)

# Intervention parameters
visibility_reduction = 0.75
sharing_reduction = 0
intervention_time = 1 #assumed time of crowd-sourced intervention


####################
# Load data 
####################
# Bring in crowd-data
for file in crowd_files:
    data = pd.read_csv(data_directory + 'data/crowdsourced_factcheck/' + file, dtype = {'in_robust_mode': object})
    if 'crowd_data' not in globals():
        crowd_data = data.copy()
    else:
        crowd_data = crowd_data.append(data)
    del data
    
    
# Load tweet data, esnure in proper format
labeled_tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                             dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                                      'user_id': 'int64', 'tweet_id': 'int64', 
                                      'retweeted_user_id': object, 'retweet_id': object})
labeled_tweets[['retweeted_user_id', 'retweet_id']] = labeled_tweets[['retweeted_user_id', 'retweet_id']].fillna(-1) #convert RT IDs to int (can't declare them upfront due to NaNs)
labeled_tweets[['retweeted_user_id', 'retweet_id']] = labeled_tweets[['retweeted_user_id', 'retweet_id']].astype('int64')

# Filter to article and tweets of interest
stories = labeled_tweets['total_article_number'].unique()
story = int(stories[which_story])
story_tweets = labeled_tweets[labeled_tweets.total_article_number == story].copy()
story_tweets['tweet_time'] = pd.to_datetime(story_tweets['tweet_time'], format = '%a %b %d %H:%M:%S %z %Y')
del labeled_tweets

# Load retweet network for this article
RT_network = pd.read_csv(data_directory + "data_derived/networks/specific_article_networks/article" + str(story)+ "_edges.csv",
                          dtype = {'Source': 'int64', 'Target': 'int64', 'source_tweet_id': 'int64', 'target_tweet_id': 'int64'})


####################
# Reshape data 
####################
# Drop professional fact-checker data
drop_cols = [col for col in crowd_data.columns if re.search("fc_", col)]
drop_cols = drop_cols + ['mode', 'in_robust_mode', 'no_false', 'all_true_set', 'any_true_set', 'no_true', 'all_false_set'	, 'any_false_set']
crowd_data = crowd_data.drop(columns = drop_cols)
del drop_cols

# Rename crowd columns
new_cols = [re.sub("resp_cat_", "veracity_categ", col) for col in crowd_data.columns]
new_cols = [re.sub("resp_veracity_", "veracity_rating", col) for col in new_cols]
new_cols = [re.sub("ideology_resp_", "ideology", col) for col in new_cols]
new_cols = [re.sub("pol_know_", "political_knowledge", col) for col in new_cols]
new_cols = [re.sub("crt_", "crtical_reasoning", col) for col in new_cols]
crowd_data.columns = new_cols
del new_cols

# Identify each row as unique simulated crowd for that article
crowd_data['crowd_id'] = crowd_data.groupby('article_num').cumcount()

# Format wide to long
crowd_data = pd.wide_to_long(crowd_data,
                             stubnames = ['veracity_categ', 'veracity_rating', 'ideology', 'political_knowledge', 'crtical_reasoning'],
                             i = ['article_num', 'crowd_id'],
                             j = 'individual')
crowd_data = crowd_data.reset_index()

# Label crowd memebers as liberal, conservative, or moderate
crowd_data['ideology_category'] = ''
crowd_data.loc[crowd_data.ideology < 0, 'ideology_category'] = 'Liberal'
crowd_data.loc[crowd_data.ideology > 0, 'ideology_category'] = 'Conservative'
crowd_data.loc[crowd_data.ideology == 0, 'ideology_category'] = 'Moderate'


####################
# Function for generating crowd verdict
####################
def generate_factcheck_crowds(crowd_data, crowd_size, article_id, n_crowds, crowd_type = "random"):
    """
    Function to create a list of 600 fact-checking crowds for a given article.
    
    INPUT:
        - crowd_size (int):   number of individuals in fact-checking crowd.
        - crowd_type (str):   type of crowd to generate, random or ideologically-balanced ("balanced")    

    OUTPUT:
        - 

    """
    # Loop through individually simulated crowds and analyze their fact-checking
    article_crowds = crowd_data[crowd_data.article_num == article_id].copy()
    crowd_ratings = pd.DataFrame(columns = ['article_id', 'crowd_id', 'crowd_size', 'true', 'false', 'cnd', 'rating_mode', 'veracity_mean', 'veracity_median'])
    if crowd_type == "random":
        for i in np.arange(n_crowds):
            
            # Basic veracity evaluation
            crowd_set = article_crowds[article_crowds.crowd_id == i].copy()
            crowd_set = crowd_set.sort_values(by = ['individual'])
            selected_crowd = crowd_set.iloc[0:crowd_size,] #grab frist ::crowd_size:: individuals
            rating_true = sum(selected_crowd.veracity_categ == 't') / crowd_size
            rating_false = sum(selected_crowd.veracity_categ == 'f') / crowd_size
            rating_cnd = sum(selected_crowd.veracity_categ == 'c') / crowd_size
            mean_veracity = selected_crowd.veracity_rating.mean()
            median_veracity = selected_crowd.veracity_rating.median()
            
            # Determine modal rating
            rating_mode = selected_crowd.veracity_categ.mode()
            crowd_rating_mode = np.nan
            if len(rating_mode) > 1:
                crowd_rating_mode = True
            elif (rating_mode.iloc[0] == 't') or (rating_mode.iloc[0] == 'c'):
                crowd_rating_mode = True
            elif rating_mode.iloc[0] == 'f':
                crowd_rating_mode = False
                
            # Add this crowd's rating
            crowd_ratings = crowd_ratings.append({'article_id': article_id,
                                                  'crowd_id': i, 
                                                  'crowd_size': crowd_size, 
                                                  'true': rating_true, 
                                                  'false': rating_false, 
                                                  'cnd': rating_cnd, 
                                                  'rating_mode': crowd_rating_mode,
                                                  'veracity_mean': mean_veracity,
                                                  'veracity_median': median_veracity},
                                                 ignore_index = True)
            
    # Return
    return crowd_ratings
   

####################
# Generate crowd ratings
#################### 
for article in crowd_data.article_num.unique():
    crowd_ratings = generate_factcheck_crowds(crowd_data, 
                                              crowd_size = crowd_size,
                                              article_id = article, 
                                              n_crowds = n_replicates,
                                              crowd_type = crowd_type)
    if 'crowd_evaluations' not in globals():
        crowd_evaluations = crowd_ratings.copy()
    else: 
        crowd_evaluations = crowd_evaluations.append(crowd_ratings, ignore_index = True)
    del crowd_ratings
crowd_evaluations.to_csv(outpath + 'simulated_crowd_factchecking.csv', index = False)

####################
# Prep data for simulated interventions
#################### 
# Grab only crowd evluation for this article
article_crowd_evals = crowd_evaluations[crowd_evaluations.article_id == story].copy()
del crowd_evaluations

# Organize tweets by time since first share of article
story_tweets['article_first_time'] = min(story_tweets['tweet_time'])
story_tweets['relative_tweet_time'] = story_tweets['tweet_time'] - story_tweets['article_first_time']
story_tweets['relative_tweet_time'] = story_tweets['relative_tweet_time'] / np.timedelta64(1, 'h') #convert to hours
story_tweets = story_tweets.sort_values('relative_tweet_time')   
story_tweets['tweet_number'] = np.arange(story_tweets.shape[0])
story_tweets = story_tweets.reset_index(drop = True)

# Only keep columns needed
story_tweets = story_tweets[['user_id', 'tweet_id', 'is_retweet', 'is_quote',
                             'tweet_time', 'article_first_time', 'relative_tweet_time', 'tweet_number', 
                             'total_article_number', 'article_fc_rating', 'source_type']]

# Create list of tweets that includes list of followers per tweets, and simulated
tweets_with_followers = match_followers_to_tweet(tweets = story_tweets, 
                                                 story_id = story, 
                                                 data_directory = data_directory)

    

####################
# Simulate intervention
#################### 
# Prep directory for output
sub_dir = "{}crowdsourced_reduce_sharing{}_visibility{}/".format(outpath, str(sharing_reduction), str(visibility_reduction))
os.makedirs(sub_dir, exist_ok = True)

# Loop through replicate simulations
all_intervention_tweets = []
all_exposure_timeseries = []

# First, no intervention
noint_tweets, noint_exposure_time = simulate_intervention(tweets = story_tweets, 
                                                          paired_tweets_followers = tweets_with_followers, 
                                                          RT_network = RT_network,
                                                          sharing_reduction = 0, #NO INTERVENTION
                                                          visibility_reduction = 0, #NO INTERVENTION
                                                          intervention_time = 0, 
                                                          replicate_number = -1,
                                                          mean_time_to_exposure = 1,
                                                          sd_time_to_exposure = 2)
noint_tweets['simulation_type'] = noint_exposure_time['simulation_type'] = "baseline"
all_intervention_tweets.append(noint_tweets)
all_exposure_timeseries.append(noint_exposure_time)

# Next, interventions based on crowd-sourced fact checking
for i in np.arange(n_replicates):
    
    # Establish if T/F by rule
    article_crowd_evals['labeled_false'] = article_crowd_evals['veracity_mean'] < 4
    
    # Simulate intervention if marked by crowd, else leave untouched (i.e., use baseline data)
    if article_crowd_evals['labeled_false'].iloc[i] == True:
        replicate_tweets, replicate_exposure_time = simulate_intervention(tweets = story_tweets, 
                                                                          paired_tweets_followers = tweets_with_followers, 
                                                                          RT_network = RT_network,
                                                                          sharing_reduction = sharing_reduction,
                                                                          visibility_reduction = visibility_reduction, 
                                                                          intervention_time = intervention_time, 
                                                                          replicate_number = i,
                                                                          mean_time_to_exposure = 1,
                                                                          sd_time_to_exposure = 2)
        replicate_tweets['simulation_type'] = replicate_exposure_time['simulation_type'] = "intervention"
    elif article_crowd_evals['labeled_false'].iloc[i] == False:
        replicate_tweets = noint_tweets.copy()
        replicate_exposure_time = noint_exposure_time.copy()
        replicate_tweets['simulation_type'] = replicate_exposure_time['simulation_type'] = "no intervention"
        
    # Append to data set
    all_intervention_tweets.append(replicate_tweets)
    all_exposure_timeseries.append(replicate_exposure_time)
    del replicate_tweets, replicate_exposure_time
    
# Bind together and save
all_intervention_tweets = pd.concat(all_intervention_tweets)
all_exposure_timeseries = pd.concat(all_exposure_timeseries)
all_intervention_tweets.to_csv(sub_dir + 'article' + str(story) + "_crowdsourced_tweets.csv", index = False)
all_exposure_timeseries.to_csv(sub_dir + 'article' + str(story) + "_crowdsourced_exposetime.csv", index = False)
