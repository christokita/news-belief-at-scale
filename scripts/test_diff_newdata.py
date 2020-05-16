#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 16:59:48 2020

@author: ChrisTokita
"""
import pandas as pd
import numpy as np


########## New tweet parsing ##########
#
# Check differencign that better parsing URLs in RTs and quoted tweets bring
#
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/"

# Followers exposed to FM news
followers_old = pd.read_csv(data_directory + "data_derived/followers/exposed_fakenews_followers_excl_fmtweeters_OLD.csv", dtype = int)
followers_new = pd.read_csv(data_directory + "data_derived/followers/exposed_fakenews_followers_excl_fmtweeters.csv", dtype = int)
new_followers = sorted(np.setdiff1d(followers_new['user_id'], followers_old['user_id']))
nolonger_followers = sorted(np.setdiff1d(followers_old['user_id'], followers_new['user_id']))

# Friends of FM tweeters
friends_old = pd.read_csv(data_directory + "data_derived/friends/friends_fm_tweeters_excl_fmtweeters_OLD.csv", dtype = int)
friends_new = pd.read_csv(data_directory + "data_derived/friends/friends_fm_tweeters_excl_fmtweeters.csv", dtype = int)
new_friends = sorted(np.setdiff1d(friends_new['user_id'], friends_old['user_id']))
nolonger_friends = sorted(np.setdiff1d(friends_old['user_id'], friends_new['user_id']))


# FM tweeters
tweeters_old = pd.read_csv(data_directory + "data_derived/tweets/unique_fm_tweeters_OLD.csv", dtype = int)
tweeters_new = pd.read_csv(data_directory + "data_derived/tweets/unique_fm_tweeters.csv", dtype = int)
new_tweeters = sorted(np.setdiff1d(tweeters_new['user_id'], tweeters_old['user_id']))
nolonger_tweeters = sorted(np.setdiff1d(tweeters_old['user_id'], tweeters_new['user_id']))


# Print summaries
print(")
