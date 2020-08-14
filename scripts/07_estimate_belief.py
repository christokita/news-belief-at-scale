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

# high level directory (external HD or cluster storage)
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD


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
# preliminary plot
####################
import matplotlib.pyplot as plt

belief_type = 'True'
data_subset = belief_data[belief_data.belief == belief_type]
#data_subset = data_subset[data_subset.total_article_number == 28]
data_subset = data_subset[data_subset.fc_rating == 'FM']
data_subset = data_subset[['ideology_score', 'belief_freq']]
sum_data = data_subset.groupby(['ideology_score']).mean()

plt.bar(sum_data.index, sum_data['belief_freq'])
plt.xlabel("Ideology")
plt.ylabel("\% believing %s" % belief_type)