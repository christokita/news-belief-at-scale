#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 12:11:34 2020

@author: ChrisTokita

Script:
Add in missing article IDs to tweets that were not matched in the first pass.
"""

import pandas as pd
import re
import numpy as np
import requests
from bs4 import BeautifulSoup
from fuzzywuzzy import fuzz

# high level directory (external HD or cluster storage)
#data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD


####################
# Step 1: Add in tweets with good fuzzy match score to headline in our article set
####################
"""
This fuzzy matching was carried out by Will Godel
The text of a tweet was compared with the headline of the article, as listed in our dataset.
After manually checking these matches, it appears that all headlines with scores >=87 are good matches.
"""

# Load fuzzy match scores, filter to good match scores
fuzzy_matches = pd.read_csv(data_directory + '/data_derived/tweets/noarticleID_tweets_with_score.csv',
                            dtype ={'user_id': object, 'tweet_id': object})
fuzzy_matches = fuzzy_matches.drop_duplicates()
fuzzy_matches = fuzzy_matches.rename(columns = {'total article number': 'total_article_number'})
good_matches = fuzzy_matches[fuzzy_matches.fuzzy_score >= 87]
good_matches = good_matches[['tweet_id', 'total_article_number']].drop_duplicates() #duplicate tweets from previous version of tweets_parsed that had duplicates

# Load labeled tweets, which already have metadata attached, and add in article IDs for good fuzzy-matched articles
tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                     dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                              'user_id': object, 'tweet_id': object, 
                              'retweeted_user_id': object, 'retweet_id': object,
                              'quoted_user_id': object, 'quoted_id': object})
tweets[['tweet_id', 'total_article_number']] = tweets.set_index('tweet_id').total_article_number.fillna(good_matches.set_index('tweet_id').total_article_number).reset_index()

# Load metadata for articles
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")
articles = articles.rename(columns = {'total article number': 'total_article_number'})
articles['source_lean'] = articles['source'].str.replace('rss_|ct_', '')
lean_dict = {'con': 'C', 'lib': 'L', 'unclear': 'U'}
articles['source_lean'] = articles['source_lean'].map(lean_dict)
articles['source_type'] = articles['source'].str.replace('_con|_lib|_unclear', '')
articles['source_type'] = articles['source_type'].str.replace('rss', 'fringe')
articles['source_type'] = articles['source_type'].str.replace('ct', 'mainstream')
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
    
# Add metadata for newly matched articles
article_metadata = articles.merge(news_evaluations, on = 'total_article_number')
article_metadata = article_metadata.merge(article_ratings, on = 'total_article_number')
tweets =  update_tweet_article_metadata(manually_matched_tweets = good_matches, 
                                        all_tweets = tweets, 
                                        article_data = article_metadata)

# Save
tweets['manual_article_assignment'] = False
tweets.loc[tweets['tweet_id'].isin(good_matches.tweet_id), 'manual_article_assignment'] = True
tweets.to_csv(data_directory + "data_derived/tweets/tweets_labeled.csv", index = False)

####################
# Step 2: Get list of unique URLs and try to match headline from fetched URL to headline in dataset
####################
# Load labeled tweets, which already have metadata attached
tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                     dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                              'user_id': object, 'tweet_id': object, 
                              'retweeted_user_id': object, 'retweet_id': object,
                              'quoted_user_id': object, 'quoted_id': object}) 

########## Pre-manual confirmation ##########

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
still_unmatched = unmatched[~unmatched['tweet_id'].isin(good_matches['tweet_id'])]
unique_unmatched_urls = still_unmatched['urls_expanded'].str.split(',').tolist()
unique_unmatched_urls = pd.DataFrame(unique_unmatched_urls).stack()
umatched_quoted_urls = still_unmatched['quoted_urls_expanded'].str.split(',').tolist()
umatched_quoted_urls = [x for x in umatched_quoted_urls if str(x) != 'nan']
umatched_quoted_urls = pd.DataFrame(umatched_quoted_urls).stack()
unique_unmatched_urls = unique_unmatched_urls.append(umatched_quoted_urls)

# Clean, count unique URLs
for i in range(len(unique_unmatched_urls)):
    unique_unmatched_urls.iloc[i] = simplify_link(unique_unmatched_urls.iloc[i])
unique_unmatched_urls = pd.DataFrame(unique_unmatched_urls.value_counts().reset_index())
unique_unmatched_urls.columns = ['urls_expanded', 'count']

# Function to fetch headline from URL
def fetch_headline(url):
    # Header to allow get request to go onto webpages that need User-Agent
    headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_5) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/50.0.2661.102 Safari/537.36'}
    # Get content
    try:
        url_request = requests.get("http://" + url,  timeout = 10, headers = headers)
        contents = BeautifulSoup(url_request.text, features = "html.parser")
        headline = contents.find('title')
        headline = headline.text
        if headline == '403 - Forbidden':
           headline = '' 
    except:
        headline = ''
    return headline
    
# Get headlines
headlines = np.array([])  
for url in unique_unmatched_urls['urls_expanded']:
    headline = fetch_headline(url)
    headlines = np.append(headlines, headline)
unique_unmatched_urls['fetched_headline'] = headlines
    
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
unique_unmatched_urls.loc[unique_unmatched_urls['fuzzy_match_score'] >= 49, 'good_match'] = True #these appear to be good matches    
unique_unmatched_urls.to_csv(data_directory + 'data_derived/articles/unique_unmatched_urls.csv', index = False)


########## Post-manual confirmation ##########

# Take in mannually confirmed matches, account for manual edits to article assignment
# T/F good match was changed for one article about Caroll Spinney
checked_unmatched_urls = pd.read_csv(data_directory + 'data_derived/articles/manuallyconfirmed_unmatched_urls.csv') 
checked_unmatched_urls.loc[~pd.isna(checked_unmatched_urls['manual_matched_ID']), 'total_article_number'] = checked_unmatched_urls.loc[~pd.isna(checked_unmatched_urls['manual_matched_ID']), 'manual_matched_ID']

# Filter to good matches and assign article ID
confirmed_matches =  checked_unmatched_urls[checked_unmatched_urls['good_match'] == True]
for j in range(confirmed_matches.shape[0]):
    # Prep links for pattern matching
    link = confirmed_matches['urls_expanded'].iloc[j]
    # Search through URLS
    has_link = still_unmatched['urls_expanded'].str.contains(link, na = False) | still_unmatched['quoted_urls_expanded'].str.contains(link, na = False)
    # Assign article number ID
    article_id = confirmed_matches['total_article_number'].iloc[j].copy()
    still_unmatched.loc[has_link, 'total_article_number'] = article_id
    
# Merge in and add metadata
round2_matched = still_unmatched[~pd.isna(still_unmatched['total_article_number'])]
tweets[['tweet_id', 'total_article_number']] = tweets.set_index('tweet_id').total_article_number.fillna(round2_matched.set_index('tweet_id').total_article_number).reset_index()
tweets.loc[tweets['tweet_id'].isin(round2_matched.tweet_id), 'manual_article_assignment'] = True
tweets =  update_tweet_article_metadata(manually_matched_tweets = round2_matched, 
                                        all_tweets = tweets, 
                                        article_data = article_metadata)

# Save!
tweets.to_csv(data_directory + "data_derived/tweets/tweets_labeled.csv", index = False)


####################
# Step 3: Output list final unmatched tweets and manually tag
####################
# Grab unidentified URLs in remaining unidentified tweets
last_unmatched_tweets = tweets[pd.isna(tweets['total_article_number'])]
last_unmatched_tweets = last_unmatched_tweets[['tweet_id', 'user_id', 'user_name', 'tweet_time', 'tweet_text',
                                               'urls', 'urls_expanded', 'url_count', 'quoted_text', 'quoted_urls',
                                               'quoted_urls_expanded', 'quoted_url_count', 'tweet_url', 'total_article_number']]
last_unmatched_tweets.update('"' + last_unmatched_tweets['tweet_id'] + '"') #add quotes so that excel doesn't mess up the long numbers when the csv file is opened
last_unmatched_tweets.to_csv(data_directory + 'data_derived/articles/tweets_to_manually_tag.csv', index = False)

# Load after manually tagging
manually_tagged_tweets = pd.read_csv(data_directory + 'data_derived/articles/tweets_to_manually_tag.csv',
                                     dtype = {'tweet_id': object, 'user_id': object})
manually_tagged_tweets = manually_tagged_tweets[['tweet_id', 'total_article_number']]
manually_tagged_tweets['tweet_id'] = manually_tagged_tweets['tweet_id'].str.replace('"', '')

# Load labeled tweets, which already have metadata attached
tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                     dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                              'user_id': object, 'tweet_id': object, 
                              'retweeted_user_id': object, 'retweet_id': object,
                              'quoted_user_id': object, 'quoted_id': object}) 

# Merge in and add metadata
tweets[['tweet_id', 'total_article_number']] = tweets.set_index('tweet_id').total_article_number.fillna(manually_tagged_tweets.set_index('tweet_id').total_article_number).reset_index()
tweets.loc[tweets['tweet_id'].isin(manually_tagged_tweets.tweet_id), 'manual_article_assignment'] = True
tweets =  update_tweet_article_metadata(manually_matched_tweets = manually_tagged_tweets, 
                                        all_tweets = tweets, 
                                        article_data = article_metadata)
# Save!
tweets.to_csv(data_directory + "data_derived/tweets/tweets_labeled.csv", index = False)