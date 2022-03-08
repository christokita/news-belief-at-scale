#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 16:15:11 2020

@author: ChrisTokita

SCRIPT:
Add user ideology and article-level data to tweet data.
"""

####################
# Load libraries and packages, set important parameters, load data
####################
import pandas as pd
import numpy as np
import re
import os
    
# high level directory (external HD or cluster storage)
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
#data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD

####################
# Parse articles to assign article number to tweet
####################
# Load articles dataset to get their URLs
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")

# Function to clean up article links in our dataset for better matching
def simplify_link(link):
    if not pd.isnull(link):
        link = re.sub('http.*//', '', link)
        link = re.sub('^www\.', '', link)
        link = re.sub('\?.*$', '', link)
        link = re.sub('/$', '', link)
    return link

# Load and parse tweets according to full and shortened URLs, assigning article number ID
tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_parsed.csv",
                     dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                              'user_id': object, 'tweet_id': object, 
                              'retweeted_user_id': object, 'retweet_id': object,
                              'quoted_user_id': object, 'quoted_id': object})
tweets = tweets.drop(['total_article_number'], axis = 1, errors = 'ignore') #drop article ID column if it had previously been assigned
tweets = tweets.join(pd.DataFrame(np.repeat(np.nan, tweets.shape[0]), columns = ['total_article_number']))
for j in range(articles.shape[0]):
    # Prep links for pattern matching
    link = simplify_link( articles['link'].iloc[j] )
    shortlink = simplify_link( articles['short link'].iloc[j] )
    # Search through URLS
    has_full_link = tweets['urls_expanded'].str.contains(link, na = False) | tweets['quoted_urls_expanded'].str.contains(link, na = False)
    if not pd.isna(shortlink):
        has_short_link_main = tweets['urls_expanded'].str.contains(shortlink, na = False) | tweets['urls'].str.contains(shortlink, na = False)
        has_short_link_quoted = tweets['quoted_urls_expanded'].str.contains(shortlink, na = False) | tweets['quoted_urls'].str.contains(shortlink, na = False)
        has_short_link = has_short_link_main | has_short_link_quoted
    elif pd.isna(shortlink):
        has_short_link = pd.Series(np.repeat(False, tweets.shape[0])) #if shortlink is nan
    has_link = has_full_link | has_short_link #boolean operator to find which indices have one of the two possible links
    # Assign article number ID
    tweets.loc[has_link, 'total_article_number'] = articles['total article number'].iloc[j]
    
# find and save tweets without article IDs
missing_ids = tweets[pd.isnull(tweets['total_article_number'])]
if missing_ids.shape[0] > 0:
    missing_article_ids = [x  for x in  np.arange(1, 166) if x not in np.unique(tweets['total_article_number'])]
    missing_articles = articles[articles['total article number'].isin(missing_article_ids)]
    missing_articles.to_csv(data_directory + "data_derived/articles/articles_notfoundintweets.csv", index = False)
    missing_ids.to_csv(data_directory + "data_derived/tweets/noarticleID_tweets.csv", index = False)
    unique_urls = missing_ids[['urls', 'urls_expanded', 'tweet_text']].drop_duplicates()
    unique_urls.to_csv(data_directory + "data_derived/tweets/noarticleID_uniqueURLs.csv", index = False) #create a list of only unique URLs


####################
# Add article-level information
####################
# Get source lean ratings
articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")
articles = articles.rename(columns = {'total article number': 'total_article_number'})
articles['source_lean'] = articles['source'].str.replace('rss_|ct_', '')
lean_dict = {'con': 'C', 'lib': 'L', 'unclear': 'U'}
articles['source_lean'] = articles['source_lean'].map(lean_dict)

# Get source veracity (mainstream or fringe)
articles['source_type'] = articles['source'].str.replace('_con|_lib|_unclear', '')
articles['source_type'] = articles['source_type'].str.replace('rss', 'fringe')
articles['source_type'] = articles['source_type'].str.replace('ct', 'mainstream')

# Filter down article metadata to columns of interest
articles = articles[['total_article_number', 'source_lean', 'source_type']]

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
labeled_tweets = tweets.merge(articles, how = 'left', on = 'total_article_number')
labeled_tweets = labeled_tweets.merge(news_evaluations, how = 'left', on = 'total_article_number')
labeled_tweets = labeled_tweets.merge(article_ratings, how = 'left', on = 'total_article_number')


####################
# Add user-level information
####################
# Load and prepare ideological scores
ideological_scores = pd.read_csv(data_directory + "data_derived/ideological_scores/cleaned_tweeter_ideology_scores.csv",
                                 dtype = {'user_id': object})
ideological_scores = ideological_scores.rename(columns = {'pablo_score': 'user_ideology'})

# Merge in ideological scores
labeled_tweets = labeled_tweets.merge(ideological_scores, how = 'left', on = 'user_id')


####################
# Calculate missing ideologies for tweeters 
####################
# These were users who SMAPP did not already have calculated
# We will calculate their ideology using a simple mean of their friends

# Determine missing ideologies
missing_ideologies = labeled_tweets[['user_id', 'user_ideology']].drop_duplicates()
missing_ideologies = missing_ideologies[pd.isna(missing_ideologies.user_ideology)]

# Load friend ideology scores
friend_ideologies = pd.read_csv(data_directory + "data_derived/ideological_scores/cleaned_friends_ideology_scores.csv",
                                dtype = {'user_id': object})

# For each tweeter, calculate mean ideology of friends
friend_files = os.listdir(data_directory + "data/friends/")
for user in missing_ideologies['user_id']:

    # Get user ID and load followers
    regex = re.compile(r"[0-9].*_%s.csv" % user)
    file = list(filter(regex.match, friend_files))
    
    # Parse friends and calculate mean. If no followers, next person
    if len(file) == 0: #no followers, no follower file
            next
    else:
        try:
            friend_list = np.genfromtxt(data_directory + "data/followers/" + file[0], dtype = str)
            friend_list = friend_list[1:len(friend_list)] #remove header, will raise error if empty
            friend_scores = friend_ideologies[friend_ideologies['user_id'].isin(friend_list)]
            missing_ideologies.loc[missing_ideologies.user_id == user, 'user_ideology'] =  np.mean(friend_scores['pablo_score'])
        except:
            next
    
# Merge in new scores
labeled_tweets[['user_id', 'user_ideology']] = labeled_tweets.set_index('user_id').user_ideology.fillna(missing_ideologies.set_index('user_id').user_ideology).reset_index()
labeled_tweets['manual_ideology_score'] = labeled_tweets['user_id'].isin(missing_ideologies['user_id'])


####################
# Calculate additional ideological metrics
####################
# Calculate ideological extremity
labeled_tweets['user_ideol_extremity'] = abs(labeled_tweets['user_ideology'])

# Create broad ideological categories: Liberal, Moderate, Conservative
ideo_threshold = 1 #how far from zero is the cutoff for non-Moderates?
labeled_tweets.loc[labeled_tweets.user_ideology < -ideo_threshold, 'user_ideol_category'] = -1 #liberal
labeled_tweets.loc[labeled_tweets.user_ideology > ideo_threshold, 'user_ideol_category'] = 1 #conservative
labeled_tweets.loc[labeled_tweets['user_ideology'].between(-ideo_threshold, ideo_threshold, inclusive = True), 'user_ideol_category'] = 0 #moderate


####################
# Write to file
####################
# Save
labeled_tweets.to_csv(data_directory + "data_derived/tweets/tweets_labeled_pre_manual_match.csv", index = False)