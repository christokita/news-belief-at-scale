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
    
    with open(filename,) as json_file:
        # Parse each individual tweet in data set 
        for tweet in json_file:
    
            # Load individual tweet
            data = json.loads(tweet)

            # Get tweet attributes of interest
            user_id = data['user']['id_str']
            user_name = data['user']['screen_name']
            tweet_time = data['created_at']
            tweet_id = data['id_str']
            favorite_count = data['favorite_count']
            retweet_count = data['retweet_count']
            quote_count = data['quote_count']
            reply_count = data['reply_count']
            quote_status_nominal = data['is_quote_status']
            quote_status = 'quoted_status' in data
            retweet_status = 'retweeted_status' in data
            tweet_url = "https://twitter.com/" + user_name + "/status/" + str(tweet_id)
            
            # If extended tweet (>140 char) get full data
            extended_tweet = 'extended_tweet' in data
            if extended_tweet:
                extended_data = data['extended_tweet']
                tweet_text = extended_data['full_text']
                url_collection = extended_data['entities']['urls']
            else: 
                url_collection = data['entities']['urls']
                tweet_text = data['text']
            # Parse urls (can contain multiple and are in nested list)
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
            
            # Check if rewteeted or quoted tweet, and if so, get info
            if quote_status:
                retweet_id = data['quoted_status']['id_str']
                rt_user = data['quoted_status']['user']
                retweet_user_id = rt_user['id_str']
                retweet_user_name = rt_user['screen_name']
            elif retweet_status:
                retweet_id = data['retweeted_status']['id_str']
                rt_user = data['retweeted_status']['user']
                retweet_user_id = rt_user['id_str']
                retweet_user_name = rt_user['screen_name']
            else:
                retweet_id = np.nan
                retweet_user_id = np.nan
                retweet_user_name = np.nan
            
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
                                       'is_quote': quote_status,
                                       'is_quote_nominal': quote_status_nominal,
                                       'is_retweet': retweet_status, 
                                       'is_extended_tweet': extended_tweet,
                                       'retweeted_user_id': retweet_user_id,
                                       'rewteeted_user_name': retweet_user_name,
                                       'retweet_id': retweet_id,
                                       'tweet_url': tweet_url}, index = [0])
            try:
                tweet_data = tweet_data.append(this_tweet, ignore_index = True, sort = False)
            except:
                tweet_data = this_tweet
            
    return tweet_data


####################
# Load data
####################
# Speficy files to be read
save_file_name = "parsed_tweets"
#source_directory = "../" #local files   
source_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external hard drive
json1 = "data/tweets/crowdsource_factchecking_prem2.json"
json2 = "data/tweets/crowdsource_factchecking.json"
json3 = "data/tweets/crowdsource_factchecking_enterprise_trial.json"
json_files = [json1, json2, json3]
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

# Save
out_path = source_directory + "data_derived/tweets/" + save_file_name + ".csv"
tweet_data.to_csv(out_path, index = False)