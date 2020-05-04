#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 16:59:48 2020

@author: ChrisTokita
"""
import pandas as pd
import numpy as np

new_data = pd.read_csv("../exposed_fakenews_followers_excl_fmtweeters_NEW.csv")
old_data = pd.read_csv("../exposed_fakenews_followers_excl_fmtweeters_OLD.csv")
old_data_prev = pd.read_csv("/Volumes/CKT-DATA/fake-news-diffusion/data_derived/followers/exposed_fakenews_followers_excl_fmtweeters.csv")
fm_tweeters = pd.read_csv("/Volumes/CKT-DATA/fake-news-diffusion/data_derived/tweets/unique_fm_tweeters.csv")

followers_old = pd.read_csv("/Volumes/CKT-DATA/fake-news-diffusion/data/followers/2020__03__03__2654.csv") 
followers_new = pd.read_csv("/Volumes/CKT-DATA/fake-news-diffusion/data/followers/2020__03__03__2654.csv", dtype = str) 

diff1 = sorted(np.setdiff1d(new_data['user_id'], old_data['user_id']))
diff2 = sorted(np.setdiff1d(old_data['user_id'], new_data['user_id']))
diff3 = sorted(np.setdiff1d(old_data_prev['user_id'], old_data['user_id']))
diff4 = sorted(np.setdiff1d(old_data_prev['user_id'], new_data['user_id']))