#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 16:15:11 2020

@author: ChrisTokita

SCRIPT:
Construct the rewteet network of fake news articles

Note: normally I prefer columns to all have lowercase names, but Gephi needs "Source" and "Target" to be capitalized.
"""

####################
# Load libraries and packages, set important parameters, load data
####################
import pandas as pd
import numpy as np
import re
import os
import multiprocessing as mp    

# Helpful for debugging/test cases
# Note i = 8 is example self-retweet
# Note i = 9,10 is example of phantom retweet
# Note i = 2014, 2603, 3777 are quotes
# Note tweet_id = '1197665879598624769' should be an indirect RT

####################
# Functions to parse retweets and quoted tweets
####################


def article_tweet_edges(article_id, tweets, articles, fuzzy_matched_URLs, friend_files):
    """
    Function that will take a set of tweets, contruct an edge+node list, and save to file.
    Most inputs are self-explanatory, but fuzzy_matched_URLs is the list of URLs generated for previously unmatched tweets (see script 03b).

    OUTPUT:
    - rt_edge: dataframe row with user ids in the 'Source' and 'Target' column, denoting the flow of the article, along with a label for the type of flow.
               This will return None if it encounters an original sharing of the article.
    """
    
    # Filter tweets to those sharing story and sort by time
    selected_tweets = tweets[tweets.total_article_number == article_id].copy()
    n_tweets = len(selected_tweets)
    
    # Loop over tweets sharing this specific article and create edgelist
    retweet_edges = pd.DataFrame(columns = ['Source', 'Target', 'RT_type', 'total_article_number', 'source_tweet_id', 'target_tweet_id'])
    for i in range(n_tweets):
        edge = determine_retweet_edges(i, selected_tweets, articles, fuzzy_matched_URLs, friend_files)
        if edge is not None:
            edge['target_tweet_id'] = selected_tweets.tweet_id.iloc[i]
            retweet_edges = retweet_edges.append(edge, sort = False)
            
    # Create total nodelist
    nodes = selected_tweets[['user_id', 'user_name', 'user_ideology', 'total_article_number']]
    node_tweet_count = pd.value_counts(nodes['user_id']).to_frame().reset_index()
    node_tweet_count = node_tweet_count.rename(columns = {'index': 'user_id', 'user_id': 'tweet_count'}) #count up multiple tweets by same user
    nodes = nodes.drop_duplicates()
    nodes = nodes.merge(node_tweet_count, on = 'user_id') #add column of tweet counts
    nodes['ID'] = nodes['user_id'] #for Gephi
    nodes['Label'] = nodes['user_name'] #for Gephi
    
    # Write to file
    retweet_edges.to_csv(data_directory + "data_derived/networks/specific_article_networks/article" + str(article_id) + "_edges.csv", index = False)
    nodes.to_csv(data_directory + "data_derived/networks/specific_article_networks/article" + str(article_id) + "_nodes.csv", index = False) 
    return None

def determine_retweet_edges(i, tweets, articles, manually_matched_URLs, friend_files):
    """
    Function that will take a single tweet and determine if it was a RT or quoted tweet, and if so, determine who the RT came from.
    
    OUTPUT:
    - rt_edge: dataframe row with user ids in the 'Source' and 'Target' column, denoting the flow of the article, along with a label for the type of flow.
               This will return None if it encounters an original sharing of the article.
    """

    # Grab specific tweet and preliinary information    
    tweet = tweets.iloc[i,:] #i is passed to this function, so grab specific tweet of interest
    user_id = tweet['user_id']
    rt_edge = pd.DataFrame({'Source': None, 'Target': user_id, 'RT_type': None, 'total_article_number': tweet['total_article_number'],
                            'source_tweet_id': None, 'target_tweet_id': tweet['tweet_id']}, index = [0])
    
    # If retweet determine who they were actually retweeting
    if tweet['is_retweet']:
        rt_edge = parse_retweet(tweet, user_id, rt_edge, friend_files, tweets)
    elif tweet['is_quote']:
        rt_edge = parse_quotedtweet(tweet, rt_edge, articles, manually_matched_URLs)
    else:
        rt_edge = None
    return rt_edge
    
    
def parse_retweet(tweet, user_id, rt_edge, friend_files, all_tweets):
    """
    Function that will take a retweet and determine who the RT came from: original tweeter or friend that tweeted
    
    OUTPUT
    - rt_edge: returns the rt_edge dataframe row with appropriate user_id filled in the 'Source' column
    """
    
    # Ensure IDs and times are in proper format
    all_tweets = all_tweets.astype({'retweeted_user_id': object, 'retweet_id': object, 'user_id': object, 'tweet_id': object}) 
    
    # Filter tweets to only those talking about this news article
    article_id = tweet['total_article_number']
    article_tweets = all_tweets[all_tweets['total_article_number'] == article_id]
    
    #  Get RT info
    rt_user_id = tweet['retweeted_user_id']
    rt_id = tweet['retweet_id']
    
    # Get friends of focal user and id of original tweeter who was RT'd
    regex_friend = re.compile(r"[0-9].*_%s.csv" % user_id)
    fr_file = list(filter(regex_friend.match, friend_files))
    if len(fr_file) == 0: #don't have friend file for user, could be because it is private
        rt_edge['Source'] = rt_user_id
        rt_edge['RT_type'] = "Presumed Phantom RT"
        return rt_edge
    else: 
        friends = np.genfromtxt(data_directory + "data/friends/" + fr_file[0], dtype = str)
        friends = np.delete(friends, 0) #drop header
    
    # If retweeted user is followed by focal user, count that as flow of tweet
    if rt_user_id in friends:
        rt_edge['Source'] = rt_user_id
        rt_edge['RT_type'] = "Direct RT"
        rt_edge['source_tweet_id'] = rt_id
        
    # Check if self-retweet
    elif rt_user_id == user_id:
        rt_edge['Source'] = rt_user_id
        rt_edge['RT_type'] = "Self RT"
        rt_edge['source_tweet_id'] = tweet['tweet_id']
        
    # Otherwise determine if indirect RT or phantom RT
    else:
        
        # Grab tweeters who are both i's friends and tweeted the news article in question
        article_tweeters = np.unique(article_tweets['user_id'])
        candidate_friends = np.intersect1d(friends, article_tweeters)
        rt_time = tweet['tweet_time']
        candidate_tweets = article_tweets[(article_tweets['user_id'].isin(candidate_friends))].copy()
        
        # Filter candidate tweets to only those that are RTing the same tweet and that occured before i's retweet
        candidate_tweets = candidate_tweets[candidate_tweets['retweet_id'] == rt_id]
        candidate_tweets = candidate_tweets[candidate_tweets['tweet_time'] < rt_time] #tweets before focal individual
        
        try:
            # Grab most recent tweet before focal tweet, Set indirect RT edge
            retweeted = candidate_tweets.loc[ candidate_tweets['tweet_time'].idxmax() ]
            rt_edge['Source'] = retweeted['user_id']
            rt_edge['RT_type'] = "Indirect RT"
            rt_edge['source_tweet_id'] = retweeted['tweet_id']
            
        except:
            # set phantom RT edge
            rt_edge['Source'] = rt_user_id
            rt_edge['RT_type'] = "Phantom RT"
            rt_edge['source_tweet_id'] = rt_id
    
    return rt_edge


def parse_quotedtweet(tweet, rt_edge, articles, fuzzy_matched_URLs):
    """
    Function that will take a quoted tweet and determine if this counts as the original FM news share or a RT
    
    OUTPUT
    - rt_edge: returns the rt_edge dataframe row with appropriate user_id filled in the 'Source' column
    """
    
    # Filter to appropriate article
    article = articles[articles['total article number'] == tweet['total_article_number']]
    
    # Prep links for pattern matching
    link = simplify_link( article.iloc[0]['link'] )
    shortlink = simplify_link( article.iloc[0]['short link'] )
    
    # If links are not present (np.nan) in tweet, replace with '' for string matching purposes
    # You can't test if a string is in x, if x is not a string (e.g., x = np.nan)
    for col in ['urls', 'urls_expanded', 'quoted_urls', 'quoted_urls_expanded']:
        if pd.isnull(tweet[col]):
            tweet.loc[col] = ''
    
    # Search through URLS to determine if in main (i.e., extra text) tweet or the quoted tweet
    if pd.isna(shortlink):
        in_main = link in tweet['urls_expanded']
        in_quoted = link in tweet['quoted_urls_expanded']
    else: 
        in_main = (link in tweet['urls_expanded']) | (shortlink in tweet['urls'])
        in_quoted = (link in tweet['quoted_urls_expanded']) | (shortlink in tweet['quoted_urls'])
        
    # Determine what to do
    if in_quoted:
        rt_edge['Source'] = tweet['quoted_user_id']
        rt_edge['RT_type'] = "Quote"
        rt_edge['source_tweet_id'] = tweet['quoted_id']
    elif in_main and not in_quoted:
        rt_edge = None #this is treated as an original tweet of the article
    elif not in_main and not in_quoted:
        result = find_unmatched_url(tweet, fuzzy_matched_URLs)
        if result is not None:
            rt_edge['Source'] = tweet['quoted_user_id']
            rt_edge['RT_type'] = "Quote"
            rt_edge['source_tweet_id'] = tweet['quoted_id']
        else:
            rt_edge = None
    
    return rt_edge

# Function to clean up links for better matching
def simplify_link(link):
    if not pd.isna(link):
        link = re.sub('http.*//', '', link)
        link = re.sub('^www\.', '', link)
        link = re.sub('\?.*$', '', link)
        link = re.sub('/$', '', link)
    return link


# Function to check for the hard to match tweets (see script: 03b)
def find_unmatched_url(tweet, fuzzy_matched_URLs):
    """
    In the case of a quoted tweet that we can't find the URL from our original list, look for a match in the set of URLs we fuzzy matched.
    
    OUTPUT
    - rt_edge: returns the rt_edge dataframe row with appropriate user_id filled in the 'Source' column
    """
    # Get matched URLs for this article in question
    URLs_of_interest = fuzzy_matched_URLs[fuzzy_matched_URLs.total_article_number == tweet.total_article_number]
    article_URLs = URLs_of_interest['urls_expanded']
    
    # Loop through and determine where it shows up
    for url in article_URLs:
        in_body = url in tweet['urls_expanded']
        in_quote = url in tweet['quoted_urls_expanded']
        
        if in_body:
            return None
        elif in_quote:
            return "retweet"
        
    # If still not found, raise note of caution
    print("For user %s in tweet i=%d, a matching link was not found in either the main or quoted tweet text! Check this." % (tweet['user_id'], tweet.name))

####################
# Main script: Construct network edges
####################
if __name__ == '__main__':
    
    # high level directory (external HD or cluster storage)
    data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#    data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD
    
    # Load tweet data, esnure in proper format
    tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                         dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                                  'user_id': object, 'tweet_id': object, 
                                  'retweeted_user_id': object, 'retweet_id': object,
                                  'quoted_user_id': object, 'quoted_id': object})
    number_of_tweets = tweets.shape[0]
    
    # Format tweet time
    tweets.loc[:,'tweet_time'] = pd.to_datetime(tweets['tweet_time'], format = '%a %b %d %H:%M:%S %z %Y')

    # Get unique articles
    unique_articles = tweets['total_article_number'].unique()
    unique_articles.sort()
    unique_articles = unique_articles[~np.isnan(unique_articles)] #drop nan
    unique_articles = unique_articles.astype(int)
    
    # Check if some articles have already been processed
    # If so, remove from list of articles to process
    processed_articles = os.listdir(data_directory + "data_derived/networks/specific_article_networks/")
    processed_articles = [file for file in processed_articles if re.match('^article', file)] #filter out hidden copies of same files
    processed_articles = [re.search('([0-9]+)', file).group(1) for file in processed_articles]
    processed_articles = np.array(processed_articles, dtype = int)
    unique_articles = np.setdiff1d(unique_articles, processed_articles)
    
    # Load articles and fuzzy/manually-matched URLs
    articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")
    checked_unmatched_urls = pd.read_csv(data_directory + 'data_derived/articles/manuallyconfirmed_unmatched_urls.csv') 
    checked_unmatched_urls.loc[~pd.isna(checked_unmatched_urls['manual_matched_ID']), 'total_article_number'] = checked_unmatched_urls.loc[~pd.isna(checked_unmatched_urls['manual_matched_ID']), 'manual_matched_ID'] #plugs in manually fixed article number
    confirmed_matches =  checked_unmatched_urls[checked_unmatched_urls['good_match'] == True]

    # Get list of friend files for reference
    friend_files = [file for file in os.listdir(data_directory + "data/friends/") if re.match('^[0-9]', file)] #filter out hidden copies of same files

    # Process article tweets in parallel, writing individual article networks to file
    pool = mp.Pool(mp.cpu_count())
    for article_id in unique_articles:
        pool.apply_async(article_tweet_edges, args = (article_id, tweets, articles, confirmed_matches, friend_files))
    pool.close()
    pool.join()
    
    # Compile individual article networks into main edge and node list
    retweet_edges = pd.DataFrame(columns = ['Source', 'Target', 'RT_type', 'total_article_number', 'tweet_id'])
    retweet_nodes = pd.DataFrame(columns = ['user_id', 'user_name', 'user_ideology', 'ID', 'Label'])
    article_files = os.listdir(data_directory + "data_derived/networks/specific_article_networks/")
    article_edge_files = [file for file in article_files if re.match('^article[0-9]+_edges', file)] #filter out hidden copies of same files
    article_node_files = [file for file in article_files if re.match('^article[0-9]+_nodes', file)] #filter out hidden copies of same files
    for j in range(len(article_node_files)):
        edges = pd.read_csv(data_directory + "data_derived/networks/specific_article_networks/" + article_edge_files[j])
        nodes = pd.read_csv(data_directory + "data_derived/networks/specific_article_networks/" + article_node_files[j])
        retweet_edges = retweet_edges.append(edges, sort = False)
        retweet_nodes = retweet_nodes.append(nodes, sort = False)
    
    # Write to file
    retweet_edges.to_csv(data_directory + "data_derived/networks/rtnetwork_edges.csv", index = False)
    retweet_nodes.to_csv(data_directory + "data_derived/networks/rtnetwork_nodes.csv", index = False)
    