#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 12:23:28 2020

@author: ChrisTokita

SCRIPT
Simulate simple interventions to prevent the spread of fake news
"""

####################
# Load packages and data, set important parameters
####################
import pandas as pd
import numpy as np
import os
import re
import sys

# high level directory (external HD or cluster storage)
#data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD
data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
outpath = data_directory + "data_derived/interventions/"

# Parameters for simulation
n_replicates = 10
intervention_time = 6
visibility_reduction = 0.5
which_story = int(sys.argv[1]) #get from command line


####################
# Functions for modeling interventions
####################
# Simualte the proposed intervention in which users are x% less likely to see a tweet after time point t
def simulate_intervention(tweets, paired_tweets_followers, RT_network, visibility_reduction, intervention_time, replicate_number):
    """
    Function that will simulate an intervention in which users are ::visibility_reduction::% less likely to see a tweet after time point ::intervention_time::
        
    OUTPUT:
        - 
    """
    # Run intervention on dataset
    # We assume people are ::odds_reduction:: less likely to see tweet sharing fake news, 
    # therefore RTs are ::odds_reduction:: less likely to happen after intervention kicks in at ::intervention_time::
    print("Calculate RT removals")
    story_RTs = tweets[tweets.is_retweet == True].copy()
    affected_RTs = story_RTs[story_RTs.relative_tweet_time >= intervention_time]
    removed = np.random.choice([True, False], 
                                size = affected_RTs.shape[0], 
                                replace = True, 
                                p = [visibility_reduction, 1-visibility_reduction])
    removed_RTs = affected_RTs.iloc[removed].copy()
    
    # Remove RTs of RTs that were already removed by intervention (i.e., can't RT a tweet that didn't happen)
    print("Calculate indirect RT removals")
    removed_indirect_RTs = RT_network.target_tweet_id[RT_network['source_tweet_id'].isin(removed_RTs.tweet_id)]
    removed_RTs = removed_RTs.append( affected_RTs[affected_RTs['tweet_id'].isin(removed_indirect_RTs)] ) #remove RTs of RTs that now did not happen
        
    # Create set of tweets that actually would've occured under this intervention
    print("Filter tweets")
    intervention_tweets = tweets[~tweets['tweet_id'].isin(removed_RTs.tweet_id)].copy()
    intervention_tweets['replicate'] = replicate_number
    intervention_tweets['tweet_number'] = np.arange(intervention_tweets.shape[0])
    
    # Estimate who would've seen a given tweet under intervention and determine first exposure of each unique individual
    print("Estimate exposure")
    exposed_followers = estimate_exposure_under_intervention(tweets = intervention_tweets, 
                                                             paired_tweets_followers = paired_tweets_followers, 
                                                             visibility_reduction = visibility_reduction, 
                                                             intervention_time = intervention_time)
    exposed_followers = exposed_followers.sort_values(by = ['relative_tweet_time', 'tweet_number'])
    exposed_followers = exposed_followers.drop_duplicates(subset = ['follower_id'], keep = 'first')
    exposed_per_tweet = exposed_followers['tweet_id'].value_counts()
    exposed_per_tweet = pd.DataFrame({'tweet_id': exposed_per_tweet.index, 'new_exposed_users': exposed_per_tweet.values})
    
    # Merge tweets and exposure count, and return
    print("Merge tweets and exposure")
    intervention_tweets = intervention_tweets.merge(exposed_per_tweet, on = 'tweet_id', how = 'left')
    intervention_tweets['new_exposed_users'] = intervention_tweets['new_exposed_users'].fillna(0)
    intervention_tweets = intervention_tweets.sort_values(by = ['relative_tweet_time', 'tweet_number'])
    intervention_tweets['cumulative_exposed'] = intervention_tweets['new_exposed_users'].cumsum()
    return intervention_tweets


# Estimate exposure of unique users to the story, assuming the story, assuming users are x% less likely to see a tweet after time point t
def estimate_exposure_under_intervention(tweets, paired_tweets_followers, visibility_reduction, intervention_time):
    """
    Function that will determine estimate many unique users would be exposed to a story that is under an intevention in viewership.
    We assume people are ::odds_reduction:: less likely to see tweet sharing article
    
    OUTPUT:
    - exposed_over_time: a dataframe that contains the time series of tweets and exposure over time for that story.
    """
    
    exposed_followers = []
    for index, tweet in tweets.iterrows():
        tweet_id = tweet['tweet_id']
        followers = paired_tweets_followers[paired_tweets_followers.tweet_id == tweet_id]
        if tweet['relative_tweet_time'] >= intervention_time:
            tweet_visible = np.random.choice([True, False], 
                                             size = followers.shape[0], 
                                             replace = True, 
                                             p = [1-visibility_reduction, visibility_reduction])
            followers = followers[tweet_visible] #only grab followers who would've seen the tweet
        exposed_followers.append(followers)
    exposed_followers = pd.concat(exposed_followers)
    return exposed_followers
      
    
def match_followers_to_tweet(tweets, data_directory):
    """
    Function that will load followers that would have seen a given tweet.
    
    OUTPUT
    - matched_followers:   dataframe listing tweet ID with the set of follower IDs that could have potentially seen it (i.e., the followers of the user who tweeted it).
    """
    
    # Get list of unique users and tweets in this set 
    unique_users = tweets['user_id'].unique()
    tweet_list = tweets[['user_id', 'tweet_id', 'tweet_time', 'tweet_number', 'relative_tweet_time']]
    
    # Create paired list of tweeters and followers
    paired_followers = []
    follower_files = os.listdir(data_directory + "data/followers/")
    for user_id in unique_users:
        regex = re.compile(r"[0-9].*_%s.csv" % user_id)
        file = list(filter(regex.match, follower_files))
        followers = load_followers(file, data_directory)
        new_row = pd.DataFrame({'user_id': user_id, 'follower_id': followers}, dtype = 'int64')
        paired_followers.append(new_row)
    paired_followers = pd.concat(paired_followers)
    
    # Merge into our list of tweets
    matched_followers = tweet_list.merge(paired_followers, how = "left", on = "user_id")
    matched_followers = matched_followers.sort_values('relative_tweet_time')
    matched_followers = matched_followers.reset_index(drop = True)
    return matched_followers
        
    
def load_followers(file, data_directory):
    """
    Function that will load the follower files or return an empty array if there are no followers for this user

    OUTPUT
    - followers:   array of follower user IDs (numpy array, str)
    """
    if len(file) == 0: #no followers, no follower file
        followers = np.array([])
    else:
        followers = np.genfromtxt(data_directory + "data/followers/" + file[0], dtype = 'int64')
        try:
            followers = followers[1:len(followers)] #remove header, will raise error if empty
        except:
            followers = np.array([]) #no followers, empty file
    return followers


####################
# Load fake news tweets
####################
# Load tweet data, esnure in proper format
labeled_tweets = pd.read_csv(data_directory + "data_derived/tweets/tweets_labeled.csv",
                             dtype = {'quoted_urls': object, 'quoted_urls_expanded': object, #these two columns cause memory issues if not pre-specified dtype
                                      'user_id': 'int64', 'tweet_id': 'int64', 
                                      'retweeted_user_id': object, 'retweet_id': object})
# Convert RT IDs to int (can't declare them upfront due to NaNs)
labeled_tweets[['retweeted_user_id', 'retweet_id']] = labeled_tweets[['retweeted_user_id', 'retweet_id']].fillna(-1)
labeled_tweets[['retweeted_user_id', 'retweet_id']] = labeled_tweets[['retweeted_user_id', 'retweet_id']].astype('int64')

# Filter to fake news tweets
fm_tweets = labeled_tweets[labeled_tweets.article_fc_rating == "FM"].copy()
fm_stories = fm_tweets['total_article_number'].unique()
story = int(fm_stories[which_story])
del labeled_tweets


####################
# Prep data for this specific story
####################
# Grab tweets of this specific article
story_tweets = fm_tweets[fm_tweets.total_article_number == story].copy()
story_tweets['tweet_time'] = pd.to_datetime(story_tweets['tweet_time'], format = '%a %b %d %H:%M:%S %z %Y')

# Organize tweets by time since first share of article
story_tweets['article_first_time'] = min(story_tweets['tweet_time'])
story_tweets['relative_tweet_time'] = story_tweets['tweet_time'] - story_tweets['article_first_time']
story_tweets['relative_tweet_time'] = story_tweets['relative_tweet_time'] / np.timedelta64(1, 'h') #convert to hours
story_tweets = story_tweets.sort_values('relative_tweet_time')   
story_tweets['tweet_number'] = np.arange(story_tweets.shape[0])
story_tweets = story_tweets.reset_index(drop = True)

# Only keep columns needed
story_tweets = story_tweets[['user_id', 'tweet_id', 'is_retweet', 'is_quote',
                             'tweet_time', 'article_first_time', 'relative_tweet_time', 'tweet_number', 
                             'total_article_number', 'article_fc_rating', 'source_type']]

# Create list of tweets that includes list of followers per tweets
print("Match tweets to followers")
tweets_with_followers = match_followers_to_tweet(tweets = story_tweets, data_directory = data_directory)

# Load retweet network
RT_network = pd.read_csv(data_directory + "data_derived/networks/specific_article_networks/article" + str(story)+ "_edges.csv",
                         dtype = {'Source': 'int64', 'Target': 'int64', 'source_tweet_id': 'int64', 'target_tweet_id': 'int64'})


####################
# Loop through replicate simulations
####################
# Prep directory for output
sub_dir = outpath + 'reduceviz' + str(visibility_reduction) + '_t' + str(intervention_time) + '/'
os.makedirs(sub_dir, exist_ok = True)

# Loop through replicate simulations
for i in np.arange(n_replicates):
    replicate_sim = simulate_intervention(tweets = story_tweets, 
                                          paired_tweets_followers = tweets_with_followers, 
                                          RT_network = RT_network,
                                          visibility_reduction = visibility_reduction, 
                                          intervention_time = intervention_time, 
                                          replicate_number = i)
    replicate_sim.to_csv(sub_dir + 'article' + str(story) + "_intervention_rep" + str(i) + ".csv")
