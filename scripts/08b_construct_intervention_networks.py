#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 10:45:24 2020

@author: ChrisTokita

SCRIPT:
Grab a specific example RT network from an intervention of interest for plotting purposes
"""

####################
# Load packages and data, set important parameters
####################
import pandas as pd
import numpy as np
import re

# high level directory (external HD or cluster storage)
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD
outpath = data_directory + "data_derived/interventions/"


# Load story and intervention of interest
story_id = 28
intervention_time = 2
visibility_reduction = 0.5
intervention_tweets_file = data_directory + "data_derived/interventions/reduceviz{}_t{}/article{}_intervention_tweets.csv".format(visibility_reduction, intervention_time, story_id)
intervention_tweets = pd.read_csv(intervention_tweets_file, dtype = {'user_id': object, 'tweet_id': object})

# Load full RT network
node_file = data_directory + "data_derived/networks/specific_article_networks/article{}_nodes.csv".format(story_id)
edge_file = data_directory + "data_derived/networks/specific_article_networks/article{}_edges.csv".format(story_id)
RT_nodes = pd.read_csv(node_file, dtype = {'user_id': object})
RT_edges = pd.read_csv(edge_file, dtype = {'Source': object, 'Target': object, 'source_tweet_id': object, 'target_tweet_id': object})


####################
# Get intervention RT network
####################
# Filter to intervention simulation of interest
intervention = intervention_tweets[intervention_tweets.replicate == 1]

# Filter RT network to tweets of in intervention
RT_edges_intervention = RT_edges[RT_edges['target_tweet_id'].isin(intervention['tweet_id'])]
RT_edges_intervention = RT_edges_intervention[RT_edges_intervention['source_tweet_id'].isin(intervention['tweet_id'])]
RT_nodes_intervention = RT_nodes[RT_nodes['user_id'].isin(intervention['user_id'])]

# Add total exposed by user 
total_exposed = intervention.groupby(['user_id'])['new_exposed_users'].sum()
total_exposed = total_exposed.reset_index()
RT_nodes_intervention = RT_nodes_intervention.merge(total_exposed, on = 'user_id')

# Save
RT_nodes_intervention.to_csv(data_directory + "data_derived/interventions/intervention_networks/article{}_reduceviz{}_t{}_nodes.csv".format(story_id, visibility_reduction, intervention_time), index = False)
RT_edges_intervention.to_csv(data_directory + "data_derived/interventions/intervention_networks/article{}_reduceviz{}_t{}_edges.csv".format(story_id, visibility_reduction, intervention_time), index = False)
