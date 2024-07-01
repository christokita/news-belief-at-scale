#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 16:58:50 2020

@author: ChrisTokita

SCRIPT:
Get list of unchecked followers who we do not already have ideology scores for.
This is in the case we get handed new follower lists and want to check them against SMAPP's ideology score database.
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import os
import re

# high level directory (external HD or cluster storage)
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD



####################
# Function to load followers
####################
def load_followers(file, data_directory):
    """
    Function that will load the follower files or return an empty array if there are no followers for this user

    OUTPUT
    - followers:   array of follower user IDs (numpy array, str)
    """
    if len(file) == 0: #no followers, no follower file
        followers = np.array([])
    else:
        followers = np.genfromtxt(data_directory + "data/followers/" + file[0], dtype = str)
        try:
            followers = followers[1:len(followers)] #remove header, will raise error if empty
        except:
            followers = np.array([]) #no followers, empty file
    return followers


####################
# Load in this user's followers
####################
follower_files = os.listdir(data_directory + "data/followers/")
user_id = '18856867'

regex = re.compile(r"[0-9].*_%s.csv" % user_id)
file = list(filter(regex.match, follower_files))
followers = load_followers(file, data_directory) 


####################
# Get set not already included in unique followers for which we do have ideology scores
####################
follower_ideologies = pd.read_csv(data_directory + "data_derived/ideological_scores/cleaned_followers_ideology_scores.csv",
                                  dtype = {'user_id': object, 'pablo_score': float})

followers_to_check = np.setdiff1d(followers, follower_ideologies.user_id)