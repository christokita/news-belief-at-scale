#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 12:36:49 2020

@author: ChrisTokita

SCRIPT:
Take in the ideology scores of followers/friends from SMAPP and format for use in this project.
"""

import pandas as pd

# Define function to clean the dataframes
def clean_ideology_scores(scores):
    # This function will intake the raw files as handed to me and will:
    # (1) Remove the index column
    # (2) Rename columns for consistency--user_id, accounts_followed, pablo_score
    # (3) Drop rows that have NA values for the score, as to shrink the size of the data
    
    scores = scores.drop(columns = ['Unnamed: 0'])
    scores = scores.rename(columns = {'0_norm': 'pablo_score'})
    scores = scores.dropna(subset = ['pablo_score'])
    return scores
    
# Load, process, write data
followers = pd.read_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/unique_followers_excl_tweeters--with_scores.csv',
                        dtype = {'user_id': object})
followers_cleaned = clean_ideology_scores(followers)
followers_cleaned.to_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/cleaned_followers_ideology_scores.csv', index = False)
del followers, followers_cleaned

friends = pd.read_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/unique_friends_excl_tweeters--with_scores.csv',
                      dtype = {'user_id': object})
friends_cleaned = clean_ideology_scores(friends)
friends_cleaned.to_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/cleaned_friends_ideology_scores.csv', index = False)
del friends, friends_cleaned

tweeters = pd.read_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/unique_tweeters_ideology_scores.csv',
                      dtype = {'user_id': object})
tweeters_cleaned = tweeters.drop(columns = ['Unnamed: 0'])
tweeters_cleaned = tweeters_cleaned.rename(columns = {'0_norm': 'pablo_score'})
tweeters_cleaned.to_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/cleaned_tweeter_ideology_scores.csv', index = False)


tweeters_old = pd.read_csv('/Volumes/CKT-DATA/fake-news-diffusion/data_derived/ideological_scores/unique_tweeters_ideology_scores_OLD.csv',
                      dtype = {'user_id': object})
