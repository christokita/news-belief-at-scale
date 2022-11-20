#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: `03c_finalize_tweet_data.py`
Date: July 26, 2020
Author: Chris Tokita
Purpose: Finalize our tweet dataset for the study.
Details:
    (Copies of data are currently stored on external harddrive and high-performance cluster.)
    This script will do two main things:
        (1) Manually confirm/label tweets that still are not associated with a particular article ID
        (2) Filter our tweets to only those that fall within a one-week window after the first article share.

Data In: CSV file of tweets with added tweeter ideology and article information, but where some tweets still need to be matched to an article.
    `<data storage location>/data_derived/tweets/tweets_labeled_pre_manual_match.csv`

Data Out: CSV containing fully labeled (all matched to article and having tweeter ideology and article metadata) tweet dataset for the study.
    `<data storage location>/data_derived/tweets/tweets_labeled.csv`

Machine: Chris' laptop
"""

import pandas as pd
import re
import numpy as np
import requests
from bs4 import BeautifulSoup
from thefuzz import fuzz
import time
import math

# high level directory (external HD or cluster storage)
#data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD


####################
# Step 1:  Get list of unique URLs and try to match headline from fetched URL to headline in dataset
####################
# Load labeled tweets, which already have metadata attached
tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled_pre_manual_match.csv",
                     dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                              'user_id': object, 'tweet_id': object, 
                              'retweeted_user_id': object, 'retweet_id': object,
                              'quoted_user_id': object, 'quoted_id': object}) 

# Function to clean up links for better matching
def simplify_link(link):
    if not pd.isnull(link):
        link = re.sub('http.*//', '', link)
        link = re.sub('^www\.', '', link)
        link = re.sub('\?.*$', '', link)
        link = re.sub('/$', '', link)
    return link

# Load originally unmatched articles, merge with good_matches to get proper tweet_id
unmatched = pd.read_csv(data_directory + '/data_derived/tweets/noarticleID_tweets.csv',
                            dtype = {'user_id': object, 'tweet_id': object})
unmatched = unmatched.drop_duplicates() #duplicate tweets from previous version of tweets_parsed that had duplicates

# Get list of URLs,
unique_unmatched_urls = unmatched['urls_expanded'].str.split(',').tolist()
unique_unmatched_urls = pd.DataFrame(unique_unmatched_urls).stack()
umatched_quoted_urls = unmatched['quoted_urls_expanded'].str.split(',').tolist()
umatched_quoted_urls = [x for x in umatched_quoted_urls if str(x) != 'nan']
umatched_quoted_urls = pd.DataFrame(umatched_quoted_urls).stack()
unique_unmatched_urls = unique_unmatched_urls.append(umatched_quoted_urls)

# Clean, count unique URLs
for i in range(len(unique_unmatched_urls)):
    unique_unmatched_urls.iloc[i] = simplify_link(unique_unmatched_urls.iloc[i])
unique_unmatched_urls = pd.DataFrame(unique_unmatched_urls.value_counts().reset_index())
unique_unmatched_urls.columns = ['urls_expanded', 'count']



####################
# Step 2:  Fetch headlines from URL
####################
# Function to fetch headline from URL
def fetch_headline(url, wait_time = 0):
    # Header to allow get request to go onto webpages that need User-Agent
    headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_5) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/50.0.2661.102 Safari/537.36'}
    # Get content
    try:
        url_request = requests.get("http://" + url,  timeout = 15, headers = headers)
        contents = BeautifulSoup(url_request.text, features = "html.parser")
        headline = contents.find('title')
        if (headline != None) & (headline != '403 - Forbidden'):
            headline = headline.text
        else:
         headline = ''
    except:
        headline = '~ERORR~'
    # Sleep and return
    time.sleep(wait_time)
    return headline
    
# Get headlines
headlines = np.array([])  
for url in unique_unmatched_urls['urls_expanded']:
    headline = fetch_headline(url)
    headlines = np.append(headlines, headline)
unique_unmatched_urls['fetched_headline'] = headlines
    
# Retry on URLs that gave us errors
for i, row in unique_unmatched_urls.iterrows():
    if row['fetched_headline'] == '~ERORR~':
        headline = fetch_headline(url, wait_time = 2)
        unique_unmatched_urls.loc[i, 'fetched_headline'] = headline
        del headline


####################
# Step 3:  Match fetched headline to article headlines in dataset
####################
# Function to compute best match headline in fuzzy match
def best_fuzzy_match(headline, articles):
    best_match_score = 0
    best_match_headline = None
    best_match_articleID = None
    for j in range(articles.shape[0]):
        text = articles.headline[j]
        match_score = fuzz.ratio(headline, text)
        if match_score > best_match_score:
            best_match_score = match_score
            best_match_headline = articles.headline[j]
            best_match_articleID = articles.total_article_number[j]
    return best_match_score, best_match_headline, best_match_articleID

# Match up to headlines with fuzzy match
articles = pd.read_csv(data_directory + 'data/articles/daily_articles.csv')
articles = articles.rename(columns = {'total article number': 'total_article_number'})
unique_unmatched_urls['fuzzy_match_score'] = np.nan
unique_unmatched_urls['matched_headline'] = np.nan
unique_unmatched_urls['total_article_number'] = np.nan

for headline in unique_unmatched_urls['fetched_headline']:
    if pd.isna(headline):
        continue
    which_row = unique_unmatched_urls['fetched_headline'] == headline
    score, matched, articledID = best_fuzzy_match(headline = headline, articles = articles)
    unique_unmatched_urls.loc[which_row, 'fuzzy_match_score'] = score
    unique_unmatched_urls.loc[which_row, 'matched_headline'] = matched
    unique_unmatched_urls.loc[which_row, 'total_article_number'] = articledID
    
# Output for manual confirmation
unique_unmatched_urls['good_match'] = False
unique_unmatched_urls['manual_matched_ID'] = None 
unique_unmatched_urls.loc[unique_unmatched_urls['fuzzy_match_score'] >= 52, 'good_match'] = True #these appear to be good matches    
unique_unmatched_urls.to_csv(data_directory + 'data_derived/articles/unique_unmatched_urls.csv', index = False)


####################
# Step 4: Assign article ID of manually confirmed URLs
####################
# Take in mannually confirmed matches, account for manual edits to article assignment
checked_unmatched_urls = pd.read_csv(data_directory + 'data_derived/articles/manuallyconfirmed_unmatched_urls.csv') 
checked_unmatched_urls.loc[~pd.isna(checked_unmatched_urls['manual_matched_ID']), 'total_article_number'] = checked_unmatched_urls.loc[~pd.isna(checked_unmatched_urls['manual_matched_ID']), 'manual_matched_ID']
checked_unmatched_urls.loc[~pd.isna(checked_unmatched_urls['manual_matched_ID']), 'good_match'] = True #manually edited article assignment is now a good match
        
# Ensure there are no duplicate tweets in our dataset (we got an error trying to merge on tweet ID due to duplicates)
tweets = tweets.sort_values('total_article_number') \
               .drop_duplicates(subset = ['tweet_id'], keep = 'first')    
               
# Ensure that we aren't going to accidentally fill in a tweet that is already labeled
unmatched = unmatched.drop_duplicates(subset = ['tweet_id']) #same problem with duplicates as above
has_ID = tweets[~tweets.total_article_number.isna()]['tweet_id']
unmatched = unmatched[~unmatched.tweet_id.isin(has_ID)]
        

# Filter to URLs with good matches to known articles and assign article ID to our unmatched tweets
# NOTE: We will loop from least shared URL to most in the event that a tweet shared multiple URLs--it'll get assigned the most popular URL article ID.
confirmed_matches =  checked_unmatched_urls[checked_unmatched_urls['good_match'] == True]
confirmed_matches = confirmed_matches.sort_values(['count'], ascending = True) \
        .reset_index(drop = True)
        
for j in range(confirmed_matches.shape[0]):
    # Prep links for pattern matching
    link = confirmed_matches['urls_expanded'].iloc[j]
    # Search through URLS
    has_link = unmatched['urls_expanded'].str.contains(link, na = False) | unmatched['quoted_urls_expanded'].str.contains(link, na = False)
    # Assign article number ID
    article_id = confirmed_matches['total_article_number'].iloc[j].copy()
    unmatched.loc[has_link, 'total_article_number'] = article_id
    
# Update our labeled tweets
matched = unmatched[~pd.isna(unmatched['total_article_number'])]
five_percent_progress = math.ceil(matched.shape[0]/20)
for i in range(matched.shape[0]):
    
    # Print prgoress
    if (i % five_percent_progress) == 0:
        progress = i / five_percent_progress * 5
        print(f'{int(progress)}%...')
        
    # Grab tweet ID and corersponding assigned article ID
    tweet_id = matched['tweet_id'].iloc[i]
    matched_article_number = matched['total_article_number'].iloc[i]
    
    # Update total article number
    tweets.loc[tweets.tweet_id == tweet_id, 'total_article_number'] = matched_article_number
    del tweet_id, matched_article_number

# Flag manually confirmed/matched tweets
tweets['manual_article_assignment'] = False
tweets.loc[tweets['tweet_id'].isin(matched.tweet_id), 'manual_article_assignment'] = True



####################
# Step 5: Add missing metadata to our newly identified articles
####################
# Function to add metadata
def update_tweet_article_metadata(manually_matched_tweets, all_tweets, article_data):
    cols_of_interest = ['source_lean', 'source_type', 'article_fc_rating', 'article_main_tag', 'article_other_tags', 'article_lean', 'article_con_feel', 'article_lib_feel']
    articles_to_add_metadata = manually_matched_tweets['total_article_number'].unique()
    for article_num in articles_to_add_metadata:
        which_rows = all_tweets['tweet_id'].isin(manually_matched_tweets['tweet_id']) & (all_tweets['total_article_number'] == article_num)
        metadata = article_data[article_data['total_article_number'] == article_num]
        metadata = metadata[cols_of_interest]
        all_tweets.loc[which_rows, cols_of_interest] = metadata.values
    return all_tweets

# Load metadata for articles
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")
articles = articles.rename(columns = {'total article number': 'total_article_number'})
articles['source_lean'] = articles['source'].str.replace('rss_|ct_', '', regex = True)
lean_dict = {'con': 'C', 'lib': 'L', 'unclear': 'U'}
articles['source_lean'] = articles['source_lean'].map(lean_dict)
articles['source_type'] = articles['source'].str.replace('_con|_lib|_unclear', '', regex = True)
articles['source_type'] = articles['source_type'].str.replace('rss', 'fringe', regex = True)
articles['source_type'] = articles['source_type'].str.replace('ct', 'mainstream', regex = True)
articles = articles[['total_article_number', 'source_lean', 'source_type']]

news_evaluations = pd.read_csv(data_directory + "data/articles/evaluations.csv")
news_evaluations = news_evaluations.rename(columns = {'article_num': 'total_article_number',
                                                      'mode of FC': 'article_fc_rating'})
news_evaluations = news_evaluations[['total_article_number', 'article_fc_rating']]
news_evaluations = news_evaluations.dropna()

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

# Update metadata for newly identified tweets
article_metadata = articles.merge(news_evaluations, on = 'total_article_number')
article_metadata = article_metadata.merge(article_ratings, on = 'total_article_number')


tweets =  update_tweet_article_metadata(manually_matched_tweets = matched, 
                                        all_tweets = tweets, 
                                        article_data = article_metadata)

####################
# Step 6: Filter to tweets that are only up to one week after first link share
####################
# Determine first time each article was shared
tweets['tweet_time_dt'] = pd.to_datetime(tweets['tweet_time'], format = '%a %b %d %H:%M:%S %z %Y')
first_tweets = tweets.groupby('total_article_number')['tweet_time_dt'] \
    .min() \
    .reset_index() \
    .rename(columns = {'tweet_time_dt': 'article_first_time'})
        
# Determine time since first share
tweets = tweets.merge(first_tweets, on = 'total_article_number', how = 'left')
tweets['relative_tweet_time'] = tweets['tweet_time_dt'] - tweets['article_first_time']
tweets['relative_tweet_time'] = tweets['relative_tweet_time'] / np.timedelta64(1, 'h') #convert to hours

# Number sequence of tweets within article
tweets['tweet_number'] = tweets.sort_values('relative_tweet_time').groupby('total_article_number').cumcount()

# Filter to tweets up to one week after the first link share and finalize dataset
tweets = tweets[tweets.relative_tweet_time <= 7*24] #one-week limit
tweets = tweets.sort_values(by = ['total_article_number', 'relative_tweet_time'])
tweets = tweets.drop(columns = ['tweet_time_dt'])

# Save!
tweets.to_csv(data_directory + "data_derived/tweets/tweets_labeled.csv", index = False)
