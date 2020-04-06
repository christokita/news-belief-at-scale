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


####################
# Load data
####################
filename = "../data/tweets/sample_tweets.json"
#filename = "../data/tweets/crowdsource_factchecking_prem2.json"

tweet_data = pd.DataFrame(columns = ['user_id', 'user_name',
                                     'tweet_time', 'tweet_text', 'tweet_id',
                                     'urls', 'urls_expanded',
                                     'is_quote', 'is_retweet', 'tweet_url'])
with open(filename,) as json_file:
    
    # Parse each tweet in data set 
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
        retweet_status = data['retweeted']
        tweet_url = "https://twitter.com/" + user_name + "/status/" + str(tweet_id)
        
        # Parse urls (cane contain multiple and are in nested list)
        url_collection = data['entities']['urls']
        url_count = 0
        urls_list = []
        urls_expanded_list = []
        for url_set in url_collection:
            urls_list.append(url_set['url'])
            urls_expanded_list.append(url_set['expanded_url'])
            url_count += 1
        url = ",".join(urls_list)
        url_expanded = ",".join(urls_expanded_list)
        
        # Append to tweet dataframe
        print(url)
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
                                   'tweet_url': tweet_url}, index = [0])
        tweet_data = tweet_data.append(this_tweet, ignore_index = True)
