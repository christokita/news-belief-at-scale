#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: `03a_clean_ideology_scores.py`
Date: July 24, 2020
Author: Chris Tokita
Purpose: Compile and format all the ideology scores of followers/friends calcualted by SMaPP.
Details:
    (Copies of data are currently stored on external hard drive and high-performance cluster.)

Data In: CSV files containing list of Twitter user IDs and corresponding ideology (Pablo) score.
    `<data storage location>/data_derived/ideological_scores/`
        `unique_followers_excl_tweeters--with_scores.csv`
        `unique_followers_zerohedge--with_scores.csv`
        `unique_friends_excl_tweeters--with_scores.csv`
        `unique_tweeters_ideology_scores--with_scores.csv`

Data Out: CSV files containing compiled and cleaned ideology scores for tweeters, followers, and friends.
    `<data storage location>/data_derived/ideological_scores/`
        `cleaned_followers_ideology_scores.csv`
        `cleaned_friends_ideology_scores.csv`
        `cleaned_tweeter_ideology_scores.csv`

Machine: Chris' laptop
"""

import pandas as pd

# Define function to clean the dataframes
def clean_ideology_scores(scores):
    # This function will intake the raw files as handed to me and will:
    # (1) Remove the index column
    # (2) Rename columns for consistency--user_id, accounts_followed, pablo_score
    # (3) Drop rows that have NA values for the score, as to shrink the size of the data
    
    scores = scores.drop(columns = ['Unnamed: 0'], errors = 'ignore')
    scores = scores.rename(columns = {'0_norm': 'pablo_score',
                                      'id_str': 'user_id'})
    scores = scores.dropna(subset = ['pablo_score'])
    return scores
    
####################
# Load, process, write follower data
####################
# Original list
followers_original = pd.read_csv('/Volumes/CKT-DATA/news-belief-at-scale/data_derived/ideological_scores/unique_followers_excl_tweeters--with_scores.csv',
                        dtype = {'user_id': object, 'id_str': object})
followers_original_cleaned = clean_ideology_scores(followers_original)
del followers_original

# Missing account list: we were previously missing scores for this account's followers, so we calculated them separately
followers_zerohedge = pd.read_csv('/Volumes/CKT-DATA/news-belief-at-scale/data_derived/ideological_scores/unique_followers_zerohedge--with_scores.csv',
                                  dtype = {'user_id': object})
followers_zerohedge_cleaned = clean_ideology_scores(followers_zerohedge)
del followers_zerohedge

# Bind together make unique, write to file
followers_cleaned = followers_original_cleaned.append(followers_zerohedge_cleaned)
followers_cleaned.shape[0] - followers_original_cleaned.shape[0] #check how many extra followers we got
del followers_original_cleaned, followers_zerohedge_cleaned

followers_cleaned = followers_cleaned.drop_duplicates(subset = ['user_id']) 
followers_cleaned.to_csv('/Volumes/CKT-DATA/news-belief-at-scale/data_derived/ideological_scores/cleaned_followers_ideology_scores.csv', index = False)
del followers_cleaned


####################
# Load and clean friend lists
####################
friends = pd.read_csv('/Volumes/CKT-DATA/news-belief-at-scale/data_derived/ideological_scores/unique_friends_excl_tweeters--with_scores.csv',
                      dtype = {'user_id': object, 'id_str': object})
friends_cleaned = clean_ideology_scores(friends)
friends_cleaned.to_csv('/Volumes/CKT-DATA/news-belief-at-scale/data_derived/ideological_scores/cleaned_friends_ideology_scores.csv', index = False)
del friends, friends_cleaned


####################
# Load and clean tweeter lists
####################
tweeters = pd.read_csv('/Volumes/CKT-DATA/news-belief-at-scale/data_derived/ideological_scores/unique_tweeters_ideology_scores.csv',
                      dtype = {'user_id': object, 'id_str': object})
tweeters_cleaned = clean_ideology_scores(tweeters)
tweeters_cleaned.to_csv('/Volumes/CKT-DATA/news-belief-at-scale/data_derived/ideological_scores/cleaned_tweeter_ideology_scores.csv', index = False)
