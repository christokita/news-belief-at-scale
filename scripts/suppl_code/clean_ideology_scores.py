#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 12:36:49 2020

@author: ChrisTokita

SCRIPT:
Take in the ideology scores of followers/friends from SMAPP and format for use in this project.

We had to piece together the ideology lists from disparate sources:
    (1) the original set of unique friend, followers, tweeters with ideology scores
    (2) missing accounts that we realized we didn't have friend/follower lists for (these will be tagged with the account name or id)
    (3) the new set of friend/follower ideology scores from new batch of Covid-19 articles
"""

import pandas as pd

# Define function to clean the dataframes
def clean_ideology_scores(scores):
    # This function will intake the raw files as handed to me and will:
    # (1) Remove the index column
    # (2) Rename columns for consistency--user_id, accounts_followed, pablo_score
    # (3) Drop rows that have NA values for the score, as to shrink the size of the data
    
    scores = scores.drop(columns = ['Unnamed: 0'], errors = 'ignore')
    scores = scores.rename(columns = {'0_norm': 'pablo_score'})
    scores = scores.dropna(subset = ['pablo_score'])
    return scores
    
####################
# Load, process, write follower data
####################
# Original list
followers_original = pd.read_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/unique_followers_excl_tweeters--with_scores.csv',
                        dtype = {'user_id': object})
followers_original_cleaned = clean_ideology_scores(followers_original)
del followers_original

# Missing account list
followers_zerohedge = pd.read_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/unique_followers_zerohedge--with_scores.csv',
                                  dtype = {'user_id': object})
followers_zerohedge_cleaned = clean_ideology_scores(followers_zerohedge)
del followers_zerohedge

# Bind together make unique, write to file
followers_cleaned = followers_original_cleaned.append(followers_zerohedge_cleaned)
followers_original_cleaned.shape[0] + followers_zerohedge_cleaned.shape[0] == followers_cleaned.shape[0] #sanity check
del followers_original_cleaned, followers_zerohedge_cleaned

followers_cleaned = followers_cleaned.drop_duplicates(subset = ['user_id']) 
followers_cleaned.to_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/cleaned_followers_ideology_scores.csv', index = False)
del followers_cleaned



####################
# Load and clean friend lists
####################
friends = pd.read_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/unique_friends_excl_tweeters--with_scores.csv',
                      dtype = {'user_id': object})
friends_cleaned = clean_ideology_scores(friends)
friends_cleaned.to_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/cleaned_friends_ideology_scores.csv', index = False)
del friends, friends_cleaned


####################
# Load and clean tweeter lists
####################
tweeters = pd.read_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/unique_tweeters_ideology_scores.csv',
                      dtype = {'user_id': object})
tweeters_cleaned = tweeters.drop(columns = ['Unnamed: 0'])
tweeters_cleaned = tweeters_cleaned.rename(columns = {'0_norm': 'pablo_score'})
tweeters_cleaned.to_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/cleaned_tweeter_ideology_scores.csv', index = False)


tweeters_old = pd.read_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/unique_tweeters_ideology_scores_OLD.csv',
                      dtype = {'user_id': object})
