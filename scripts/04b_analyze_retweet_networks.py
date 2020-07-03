#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 16:15:11 2020

@author: ChrisTokita

SCRIPT:
Analyze the rewteet network of fake news articles
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import re
import os
import networkx as nx
import scipy.stats as stats

    
# high level directory (external HD or cluster storage)
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD

####################
# Load and bind together rt networks (temporary)
####################
# Load data
#rt_edges = pd.read_csv(data_directory + "/data_derived/networks/rtnetwork_edges.csv")
#rt_nodes = pd.read_csv(data_directory + "/data_derived/networks/rtnetwork_nodes.csv")

# Load and bind together rt networks (temporary)
article_files = os.listdir(data_directory + "data_derived/networks/specific_article_networks/")
article_edge_files = [file for file in article_files if re.match('^article[0-9]+_edges', file)] #filter out hidden copies of same files
article_node_files = [file for file in article_files if re.match('^article[0-9]+_nodes', file)] #filter out hidden copies of same files
for file in article_node_files:
    nodes = pd.read_csv(data_directory + "data_derived/networks/specific_article_networks/" + file)
    if 'rt_nodes' not in globals():
        rt_nodes = nodes
    else:
        rt_nodes = rt_nodes.append(nodes, sort = False)
    del nodes
for file in article_edge_files:
    edges = pd.read_csv(data_directory + "data_derived/networks/specific_article_networks/" + file)
    if 'rt_edges' not in globals():
        rt_edges = edges
    else:
        rt_edges = rt_edges.append(edges, sort = False)
    del edges

# Add in article information
    
####################
# Add article-level information
####################
# Get source lean ratings
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")
articles = articles.rename(columns = {'total article number': 'total_article_number'})
articles['source_lean'] = articles['source'].str.replace('rss_|ct_', '')
lean_dict = {'con': 'C', 'lib': 'L', 'unclear': 'U'}
articles['source_lean'] = articles['source_lean'].map(lean_dict)
articles = articles[['total_article_number', 'source_lean']]

# Load and prepare article veracity evaluations
news_evaluations = pd.read_csv(data_directory + "data/articles/evaluations.csv")
news_evaluations = news_evaluations.rename(columns = {'article_num': 'total_article_number',
                                                      'mode of FC': 'article_fc_rating'})
news_evaluations = news_evaluations[['total_article_number', 'article_fc_rating']]
news_evaluations = news_evaluations.dropna()

# Load and prepare article lean ratings
article_ratings = pd.read_csv(data_directory + "data/articles/article_level.csv")
article_ratings = article_ratings.rename(columns = {'article_num': 'total_article_number', 
                                                    'main_tag': 'article_main_tag',
                                                    'lean': 'article_lean',
                                                    'other_tags': 'article_other_tags',
                                                    'con_feel': 'article_con_feel',
                                                    'lib_feel': 'article_lib_feel'})
article_ratings = article_ratings.drop(columns = ['daily article number', 'partisan_diff', 
                                                  'Unnamed: 8', 'Unnamed: 9', 
                                                  'Unnamed: 10', 'Unnamed: 11'])
    
# Merge in article and source info
rt_nodes = rt_nodes.merge(articles, how = 'left', on = 'total_article_number')
rt_nodes = rt_nodes.merge(news_evaluations, how = 'left', on = 'total_article_number')
rt_nodes = rt_nodes.merge(article_ratings, how = 'left', on = 'total_article_number')


####################
# Calculate network metrics by article
####################
# Get unique articles
unique_articles = rt_nodes['total_article_number'].unique().astype(int)
unique_articles.sort()

# Loop over articles and calculate basic metrics
article_network_metrics = pd.DataFrame(columns = ['total_article_number', 'total_tweets', 'total_RTs',
                                                  'network_density',
                                                  'ideology_mean', 'ideology_sd', 'ideology_se',
                                                  'article_fc_rating', 'article_lean'])
for article_id in unique_articles:
    
    # Filter to nodes and edges of interest
    article_nodes = rt_nodes[rt_nodes.total_article_number == article_id]
    article_edges = rt_edges[rt_edges.total_article_number == article_id]
    
    # Create graph
    g = nx.DiGraph() #directed
    g.add_nodes_from(article_nodes['user_id'])
    g.add_edges_from( list(zip(article_edges['Source'], article_edges['Target'])) )
    
    # Calculate metrics of interest
    network_density = nx.density(g)
    total_tweets = sum(article_nodes['tweet_count'])
    total_RTs = article_edges.shape[0]
    ideology_mean = np.mean(article_nodes['user_ideology'])
    ideology_sd = np.std(article_nodes['user_ideology'])
    ideology_se = stats.sem(article_nodes['user_ideology'], nan_policy = 'omit')
    
    # Append
    article_data = pd.DataFrame([[article_id, total_tweets, total_RTs,
                                 network_density,
                                 ideology_mean, ideology_sd, ideology_se,
                                 article_nodes['article_fc_rating'].unique()[0], article_nodes['article_lean'].unique()[0]]], 
                                columns = article_network_metrics.columns, index = [0])
    article_network_metrics = article_network_metrics.append(article_data, sort = False)