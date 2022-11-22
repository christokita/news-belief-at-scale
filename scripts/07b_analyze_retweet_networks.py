#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: `07b_analyze_retweet_networks.py`
Date: April 30, 2020
Author: Chris Tokita
Purpose: Analyze article retweet networks.
Details:
    (Copies of data are currently stored on external harddrive and high-performance cluster.)
 
Data In: CSV files of (a) retweet network edges, (b) retweet network nodes, (c) article list, (d) article fact-check evaluations, and (e) article ideological lean ratings.
    (a) `<data storage location>/data_derived/networks/rtnetwork_edges.csv`
    (b) `<data storage location>/data_derived/networks/rtnetwork_nodes.csv`
    (c) `<data storage location>/"data/articles/daily_articles.csv`
    (d) `<data storage location>/data/articles/evaluations.csv`
    (e) `<data storage location>/data/articles/article_level.csv`

Data Out: CSV file of network metrics for each article.
    `<data storage location>/data_derived/networks/article_network_metrics.csv`

Machine: High-performance computing cluster
    This script is batched to the cluster using `slurm_scripts/simulate_intervention.cmd`
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import networkx as nx
import scipy.stats as stats

    
# high level directory (external HD or cluster storage)
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD

####################
# Load rt networks (temporary)
####################
# Load data
rt_edges = pd.read_csv(data_directory + "/data_derived/networks/rtnetwork_edges.csv",
                       dtype = {'Source': object, 'Target': object, 'tweet_id': object})
rt_nodes = pd.read_csv(data_directory + "/data_derived/networks/rtnetwork_nodes.csv",
                       dtype = {'user_id': object, 'ID': object})

# Add article-level information
# Get source lean ratings
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")
articles = articles.rename(columns = {'total article number': 'total_article_number'})
articles['source_lean'] = articles['source'].str.replace('rss_|ct_', '', regex = True)
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
unique_articles = unique_articles[unique_articles > 10] #excluding first 10 articles

# Loop over articles and calculate basic metrics
article_network_metrics = pd.DataFrame(columns = ['total_article_number', 'total_tweets', 'total_RTs',
                                                  'network_density',
                                                  'ideology_mean', 'ideology_sd', 'ideology_se',
                                                  'article_fc_rating', 'article_lean'])

# igraph uses C and it can't handle long numbers, so for each graph convert user_id to a within-graph id
unique_ids = rt_nodes['user_id'].unique()
    
    
for article_id in unique_articles:
    
    # Filter to nodes and edges of interest
    article_nodes = rt_nodes[rt_nodes.total_article_number == article_id]
    article_edges = rt_edges[rt_edges.total_article_number == article_id]
    
    # Create graph
    g = nx.DiGraph() #directed
    g.add_nodes_from(article_nodes['user_id'])
    g.add_edges_from( list(zip(article_edges['Source'], article_edges['Target'])) )
    
    # TEST: import ideology scores for assortativity calculation
#    ideologies = article_nodes['user_ideology'].replace(np.nan, 0)
#    ideologies = round(ideologies*1000, 0)
#    ideologies = ideologies.astype(int)
#    nx.set_node_attributes(g, ideologies, name = 'user_ideology')
#    assort = nx.numeric_assortativity_coefficient(g, 'user_ideology')
    
    # Try igraph assortativity
    # igraph uses C and it can't handle long numbers, so for each graph convert user_id to a within-graph id
    unique_ids = article_nodes['user_id'].unique()
    node_ids = np.array([ np.where(unique_ids == x) for x in article_nodes['user_id'] ]).flatten()
    sources = np.array([ np.where(unique_ids == x) for x in article_edges['Source'] ]).flatten()
    targets = np.array([ np.where(unique_ids == x) for x in article_edges['Target'] ]).flatten()
#    g_assort = igraph.Graph(vertex_attrs = {"label": node_ids, 
#                                            edges = list(zip(sources, targets)), 
#                                            directed = True)
    
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
    
# Write to file
article_network_metrics.to_csv(data_directory + "data_derived/networks/article_network_metrics.csv", index = False)
