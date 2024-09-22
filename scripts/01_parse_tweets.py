#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: `01_parse_tweets.py`
Date: April 3, 2020
Author: Chris Tokita
Purpose: Parse the JSON files containing tweet data, turn it into a pandas dataframe, and write to csv.
Details:
    (Copies of data are currently stored on external hard drive and high-performance cluster.)

Data In: JSON files containing tweet data pulled directly from Twitter API.
    `<data storage location>/data/tweets/`    

Data Out: The CSV file contains the extracted tweet information in a flat, dataframe format.
    `<data storage location>/data_derived/tweets/`

Machine: Chris' laptop
"""

####################
# Load libraries and packages
####################
import json
import pandas as pd
import numpy as np
import copy


####################
# Define custom function to load and parse JSON tweet data
####################

def tweet_parser(filename):
    """
    Parse JSON tweet data and return a compiled dataset
    
    INPUT:
    - filename: path to specific JSON file containing tweet data (string).
    
    OUTPUT:
    - tweet_data: dataframe of all tweet data in the specified file, with a subset of the relevant data headings (pandas dataframe).
    """
    
    with open(filename,) as json_file:
        # Parse each individual tweet in data set 
        for tweet in json_file:
    
            # Load individual tweet
            tweet_obj = json.loads(tweet)

            # Get tweet attributes of interest
            user_id = tweet_obj['user']['id_str']
            user_name = tweet_obj['user']['screen_name']
            tweet_time = tweet_obj['created_at']
            tweet_id = tweet_obj['id_str']
            favorite_count = tweet_obj['favorite_count']
            retweet_count = tweet_obj['retweet_count']
            quote_count = tweet_obj['quote_count']
            reply_count = tweet_obj['reply_count']
            is_quote = 'quoted_status' in tweet_obj
            is_retweet = 'retweeted_status' in tweet_obj
            tweet_url = "https://twitter.com/" + user_name + "/status/" + str(tweet_id)
            RTofselfRT = False #flag that finds  RT of self-RT by finding 'quoted_status' in 'retweeted_status'
            
            # Parse other tweet information based on tweet type
            # If the object is a retweet, we will treat that as the main tweet body for text and urls.
            if is_retweet:
                retweet_obj = tweet_obj['retweeted_status']
                retweet_id = retweet_obj['id_str']
                retweet_user_id = retweet_obj['user']['id_str']
                retweet_user_name = retweet_obj['user']['screen_name']
                tweet_text, url, url_expanded, url_count, is_extended = get_tweet_text(retweet_obj)
    
                # This catches where someone RT a self-RT.
                if 'quoted_status' in retweet_obj:
                    RTofselfRT = True
                
            else:
                tweet_text, url, url_expanded, url_count, is_extended = get_tweet_text(tweet_obj)
                retweet_id = retweet_user_id = retweet_user_name = np.nan
                
            # If object is quoted tweet, the quoted tweet text and urls will be handled separately.
            if is_quote:
                quoted_obj = tweet_obj['quoted_status']
                quoted_id = quoted_obj['id_str']
                quoted_user_id =  quoted_obj['user']['id_str']
                quoted_user_name = quoted_obj['user']['screen_name']
                quoted_text, quoted_url, quoted_url_expanded, quoted_url_count, quote_extended = get_tweet_text(quoted_obj)
            else: 
                quoted_id = quoted_user_id = quoted_user_name = np.nan
                quoted_text = quoted_url = quoted_url_expanded = quoted_url_count=quote_extended = np.nan
                
            # Catch for quoted tweets that quoted since-deleted tweet (would be missing quoted_status object)
            quoted_tweet_deleted = False
            if tweet_obj['is_quote_status']:
                is_quote = True
                quoted_tweet_deleted = True
            
            # Append to tweet dataframe
            this_tweet = pd.DataFrame({'user_id': user_id, 
                                       'user_name': user_name,
                                       'tweet_time': tweet_time,
                                       'tweet_text': tweet_text,
                                       'tweet_id': tweet_id,
                                       'urls': url,
                                       'urls_expanded': url_expanded,
                                       'url_count': url_count,
                                       'favorite_count': favorite_count,
                                       'retweet_count': retweet_count,
                                       'quote_count': quote_count,
                                       'reply_count': reply_count,
                                       'is_quote': is_quote,
                                       'is_retweet': is_retweet, 
                                       'is_extended_tweet': is_extended,
                                       'retweeted_user_id': retweet_user_id,
                                       'retweeted_user_name': retweet_user_name,
                                       'retweet_id': retweet_id,
                                       'quoted_user_id': quoted_user_id,
                                       'quoted_user_name': quoted_user_name,
                                       'quoted_id': quoted_id,
                                       'quoted_text': quoted_text,
                                       'quoted_urls': quoted_url,
                                       'quoted_urls_expanded': quoted_url_expanded,
                                       'quoted_url_count': quoted_url_count,
                                       'quoted_is_extended': quote_extended,
                                       'quoted_tweet_deleted': quoted_tweet_deleted,
                                       'RTofselfRT': RTofselfRT,
                                       'tweet_url': tweet_url}, index = [0])
            try:
                tweet_data = tweet_data.append(this_tweet, ignore_index = True, sort = False)
            except:
                tweet_data = this_tweet
            
    return tweet_data


def get_tweet_text(tweet_object):
    """
    Function to get tweet text and URLs from tweet data.
    
    INPUT
    - tweet_object: JSON-parsed data. Can be main tweet object, retweeted_status object, or quoted_status object
    """
    
    # If extended tweet (>140 char), go to proper sub-object.
    is_extended_tweet = 'extended_tweet' in tweet_object
    if is_extended_tweet:
        extended_tweet_object = tweet_object['extended_tweet']
        tweet_text = extended_tweet_object['full_text']
        url_collection = extended_tweet_object['entities']['urls']
    else: 
        url_collection = tweet_object['entities']['urls']
        tweet_text = tweet_object['text']
        
    # Parse urls (can contain multiple and are in nested list)
    url, url_expanded, url_count = parse_urls(url_collection)
    return tweet_text, url, url_expanded, url_count, is_extended_tweet


def parse_urls(url_collection):
    """
    Function to parse a collection of URLs
    
    INPUT
    - url_collection: the URL object from the tweet object or sub-object (retweet_status or quoted_status)
    """
    
    url_count = len(url_collection)
    if url_count > 0:
        urls_list = []
        urls_expanded_list = []
        for url_set in url_collection:
            urls_list.append(url_set['url'])
            urls_expanded_list.append(url_set['expanded_url'])
        url = ",".join(urls_list)
        url_expanded = ",".join(urls_expanded_list)
    else:
        url = ""
        url_expanded = ""
    return url, url_expanded, url_count
    


####################
# Load data
####################
# Speficy files to be read
save_file_name = "tweets_parsed"
# source_directory = "/Volumes/CKT-DATA/news-belief-at-scale/" #external hard drive
source_directory = "/scratch/gpfs/ctokita/news-belief-at-scale/" #storage on HPC cluster
json1 = "data/tweets/crowdsource_factchecking_prem2.json" #original set of article tweets
json2 = "data/tweets/crowdsource_factchecking.json" #original set of article tweets
json3 = "data/tweets/crowdsource_factchecking_enterprise_trial.json" #original set of article tweets
json4 = "data/tweets/crowdsource_factchecking_missing_article_tweets.json" #found tweets for original articles without any associated tweets
json5 = "data/tweets/expanded_window_tweets.json" #tweets from expanded time window search
json_files = [json1, json2, json3, json4, json5]
json_files = [source_directory + file for file in json_files]

# Parse data
if 'tweet_data' in globals():
    del tweet_data
for file in json_files:
    data = tweet_parser(file)
    try:
        tweet_data = tweet_data.append(data, ignore_index = True, sort = False)
    except:
        tweet_data = copy.deepcopy(data)
    del data

# Drop duplicate tweets (~5,800 originally discovered)
tweet_data = tweet_data.drop_duplicates()

# Save
out_path = source_directory + "data_derived/tweets/" + save_file_name + ".csv"
tweet_data.to_csv(out_path, index = False)