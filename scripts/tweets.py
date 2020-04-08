#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 16:14:18 2020

@author: ChrisTokita

SCRIPT:
Looking at tweet data
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
    # Parse JSON tweet data and return a compiled dataset
    #
    # INPUT:
    # - filename: path to specific JSON file containing tweet data (string).
    #
    # OUTPUT:
    # - tweet_data: dataframe of all tweet data in the specified file, with a subset of the relevant data headings (pandas dataframe).
    #
    tweet_data = pd.DataFrame(columns = ['user_id', 'user_name',
                                         'tweet_time', 'tweet_text', 'tweet_id',
                                         'urls', 'urls_expanded',
                                         'is_quote', 'is_retweet', 'tweet_url'])
    with open(filename,) as json_file:
        # Parse each individual tweet in data set 
        for tweet in json_file:
    
            # Load individual tweet
            data = json.loads(tweet)

            # Get tweet attributes of interest
            user_id = data['user']['id']
            user_name = data['user']['screen_name']
            tweet_time = data['created_at']
            tweet_text = data['text']
            tweet_id = data['id']
            quote_status = data['is_quote_status']
            tweet_url = "https://twitter.com/" + user_name + "/status/" + str(tweet_id)
            
            # Parse urls (cane contain multiple and are in nested list)
            url_collection = data['entities']['urls']
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
                urls_list = ""
                urls_expanded_list = ""
            
            # Check if rewteeted and if so, get info
            if 'retweeted_status' in data:
                retweet_status = True
                rt_user = data['retweeted_status']['user']
                retweet_user_id = rt_user['id']
            else:
                retweet_status = False
                retweet_user_id = np.nan
            
            # Append to tweet dataframe
            this_tweet = pd.DataFrame({'user_id': user_id, 
                                       'user_name': user_name,
                                       'tweet_time': tweet_time,
                                       'tweet_text': tweet_text,
                                       'tweet_id': tweet_id,
                                       'urls': url,
                                       'urls_expanded': url_expanded,
                                       'url_count': url_count,
                                       'is_quote': quote_status,
                                       'is_retweet': retweet_status, 
                                       'retweeted_user_id': retweet_user_id,
                                       'tweet_url': tweet_url}, index = [0])
            tweet_data = tweet_data.append(this_tweet, ignore_index = True, sort = False)
            
    return tweet_data


####################
# Load data
####################
# Speficy files to be read
json1 = "../data/tweets/crowdsource_factchecking_prem2.json"
json2 = "../data/tweets/crowdsource_factchecking.json"
json3 = "../data/tweets/crowdsource_factchecking_enterprise_trial.json"
json_files = [json1, json2, json3]

# Parse data
if 'tweet_data' in globals():
    del tweet_data
for file in json_files:
    data = tweet_parser(file)
    if 'tweet_data' not in globals():
        tweet_data = copy.deepcopy(data)
    else:
        tweet_data = tweet_data.append(data, ignore_index = True, sort = False)

