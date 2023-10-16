#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: `07b_analyze_retweet_networks.py`
Date: October 16, 2020
Author: Chris Tokita
Purpose: Simulate simple platofmr-level interventions to prevent the spread of misinformation on Twitter.
Details:
    (Copies of data are currently stored on external hard drive and high-performance cluster.)
    To fully simulate interventions, we need to essentially rerun the entire exposure/belief estimation pipeline each time.
    Thus, this is a long script since it needs to figure out which retweets did not occur 
    and then rerun estimations of exposure/belief under the intervention parameters.
 
Data In: CSV files of (a) tweets, (b) known ideologies of followers, (c) inferred distributions for ideology of tweeters' followers, (d) article belief survey data, (e) article retweet networks.
    (a) `<data storage location>/data_derived/tweets/tweets_labeled.csv`
    (b) `<data storage location>/data_derived/ideological_scores/cleaned_followers_ideology_scores.csv`
    (c) `<data storage location>/data_derived/ideological_scores/estimated_ideol_distributions/follower_ideology_distribution_shapes.csv`
    (d) `<data storage location>/data/article_belief/response_distribution.p`
    (e) `<data storage location>/data_derived/networks/specific_article_networks/`

Data Out: CSV files of article tweets and exposure over time under specific intervention parameters (e.g., sharing reduction, visibility reduction, belief reduction)
    `<data storage location>/data_derived/interventions/{}reduce_sharing{}_visibility{}_belief{}/`

Machine: High-performance computing cluster
    This script is batched to the cluster using `slurm_scripts/simulate_intervention.cmd`
"""

####################
# Load packages 
####################
import pandas as pd
import numpy as np
import scipy.stats as stats 
import os
import re
import sys
import pickle



####################
# Paths and parameters for simulation
####################
# high level directory (external HD or cluster storage)
# data_directory = "/Volumes/CKT-DATA/news-belief-at-scale/" #external HD
data_directory = "/scratch/gpfs/ctokita/news-belief-at-scale/" #HPC cluster storage
outpath = data_directory + "data_derived/interventions/"

# Parameters for simulation
n_replicates = 50
visibility_reduction = 0.75
sharing_reduction = 0
belief_reduction = 0



####################
# Functions for modeling interventions
####################
# Simualte the proposed intervention in which users are x% less likely to see a tweet after time point t
def simulate_intervention(tweets, paired_tweets_followers, ideologies, follower_ideol_distributions, RT_network, article_belief_data, sharing_reduction, visibility_reduction, belief_reduction, intervention_time, replicate_number, mean_time_to_exposure, sd_time_to_exposure):
    """
    Function that will simulate an intervention in which users are 
    ::sharing_reduction::% less likely to share and 
    ::visibility_reduction::% less likely to see
    a tweet after time point ::intervention_time::.
        
    OUTPUT:
        - intervention_tweets (dataframe):   list of tweets that occured under this simulation. 
        - exposure_over_time (dataframe):    binned count of newly and cumulatively exposed users over time.
    """
    # Run intervention on dataset
    # We assume people are ::odds_reduction:: less likely to see tweet sharing fake news, 
    # therefore RTs are ::odds_reduction:: less likely to happen after intervention kicks in at ::intervention_time::
    story_RTs = tweets[tweets.is_retweet == True].copy()
    affected_RTs = story_RTs[story_RTs.relative_tweet_time >= intervention_time]
    removed = np.random.choice([True, False], 
                                size = affected_RTs.shape[0], 
                                replace = True, 
                                p = [sharing_reduction, 1-sharing_reduction])
    removed_RTs = affected_RTs.iloc[removed].copy()
    
    # Remove RTs of RTs that were already removed by intervention (i.e., can't RT a tweet that didn't happen)
    removed_indirect_RTs = RT_network.target_tweet_id[RT_network['source_tweet_id'].isin(removed_RTs.tweet_id)]
    removed_RTs = removed_RTs.append( affected_RTs[affected_RTs['tweet_id'].isin(removed_indirect_RTs)] ) #remove RTs of RTs that now did not happen
        
    # Create set of tweets that actually would've occured under this intervention
    intervention_tweets = tweets[~tweets['tweet_id'].isin(removed_RTs.tweet_id)].copy()
    intervention_tweets['replicate'] = replicate_number
    intervention_tweets['tweet_number'] = np.arange(intervention_tweets.shape[0])
    intervention_tweets.insert(0, 'sharing_reduction', sharing_reduction) #insert column labeling intervention info at beginning
    intervention_tweets.insert(0, 'visibility_reduction', visibility_reduction) #insert column labeling intervention info at beginning
    intervention_tweets.insert(0, 'intervention_time', intervention_time) #insert column labeling intervention info at beginning
    
    # Estimate who would've seen a given tweet under intervention and determine first exposure of each unique individual
    exposed_followers = estimate_exposure_under_intervention(tweets = intervention_tweets, 
                                                             paired_tweets_followers = paired_tweets_followers, 
                                                             visibility_reduction = visibility_reduction, 
                                                             intervention_time = intervention_time,
                                                             mean_exposure_time = mean_time_to_exposure,
                                                             sd_exposure_time = sd_time_to_exposure)
    exposed_followers = exposed_followers.sort_values(by = ['relative_exposure_time'])
    exposed_followers = exposed_followers.drop_duplicates(subset = ['follower_id'], keep = 'first')
    exposed_per_tweet = exposed_followers['tweet_id'].value_counts()
    exposed_per_tweet = pd.DataFrame({'tweet_id': exposed_per_tweet.index, 'new_exposed_users': exposed_per_tweet.values})
    
    # Estimate belief for a tweet under intervention
    believing_followers = estimate_belief(exposed_users = exposed_followers, 
                                          ideologies = ideologies, 
                                          follower_ideol_distributions = follower_ideol_distributions, 
                                          article_id = tweets['total_article_number'].unique(), 
                                          article_belief_data = article_belief_data,
                                          intervention_time = intervention_time,
                                          belief_reduction = belief_reduction)
    
    # Bin exposure and belief by time
    exposure_bin_counts, bin_edges = np.histogram(exposed_followers.relative_exposure_time, bins = 14*24*4, range = (0, 14*24)) #bin over first 2 weeks
    belief_bin_counts, _ = np.histogram(believing_followers.relative_exposure_time, bins = 14*24*4, range = (0, 14*24)) #bin over first 2 weeks
    exposure_over_time = pd.DataFrame({'intervention_time': intervention_time,
                                       'visibility_reduction': visibility_reduction,
                                       'sharing_reduction': sharing_reduction,
                                       'belief_reduction': belief_reduction,
                                       'total_article_number': intervention_tweets.total_article_number.unique().item(), 
                                       'replicate': replicate_number,
                                       'time': bin_edges[1:], 
                                       'new_exposed_users': exposure_bin_counts,
                                       'new_believing_users': belief_bin_counts})
    exposure_over_time['cumulative_exposed'] = exposure_over_time['new_exposed_users'].cumsum()
    exposure_over_time['cumulative_believing'] = exposure_over_time['new_believing_users'].cumsum()
    
    # Merge tweets and exposure count, and return
    intervention_tweets = intervention_tweets.merge(exposed_per_tweet, on = 'tweet_id', how = 'left')
    intervention_tweets['new_exposed_users'] = intervention_tweets['new_exposed_users'].fillna(0)
    intervention_tweets = intervention_tweets.sort_values(by = ['relative_tweet_time', 'tweet_number'])
    intervention_tweets['cumulative_exposed'] = intervention_tweets['new_exposed_users'].cumsum()
    return intervention_tweets, exposure_over_time

# Estimate exposure of unique users to the story, assuming the story, assuming users are x% less likely to see a tweet after time point t
def estimate_belief(exposed_users, ideologies, follower_ideol_distributions, article_id, article_belief_data, intervention_time, belief_reduction = 0.0):
    """
    Function that calculates the estimated belief among users exposued to tweets
    
    PARAMTERS:
    - exposed_users (dataframe):                 dataframe of all paired tweeters and exposed users
    - ideologies (dataframe):                    dataframe of all known user ideologies
    - follower_ideol_distributions(dataframe):   dataframe of the estimated distribution of ideology each tweeter's followers
    - article_id (int):                          article number
    - article_belief_data (dataframe):           survey data of frequency of individual belief by demographics (e.g., ideology)
    - interention_time (float):                  relative time at which an intervention kicks in.
    - belief_reduction (float):                  reduction in belief, accounting for intervention. Must be between 0.0 and 1.0, with 0.0 = 0% reduction and 1.0 = 100% reduction.
    
    OUTPUT:
    - 
    """
    # Add in ideology scores
    belief_users = exposed_users.merge(ideologies, on = 'follower_id', how = 'left')
    
    # Loop through users and draw missing ideologies from user's follower ideology distribution
    unique_tweeters = belief_users['user_id'].unique()
    for tweeter in unique_tweeters:
        # Get user's follower ideology distribution shape
        follower_mu = follower_ideol_distributions.mu[follower_ideol_distributions.user_id == tweeter].iloc[0]
        follower_sigma = follower_ideol_distributions.sigma[follower_ideol_distributions.user_id == tweeter].iloc[0]
        follower_alpha = follower_distributions.alpha[follower_ideol_distributions.user_id == tweeter].iloc[0]
        
        # Draw ideology samples for missing users from fit skew normal distribution
        n_samples = belief_users.loc[(belief_users.user_id == tweeter) & pd.isna(belief_users.pablo_score)].shape[0] #number of unknown follower ideologies for this user
        ideol_samples = stats.skewnorm(
            loc=follower_mu, 
            scale=follower_sigma, 
            a=follower_alpha
            ).rvs(n_samples)
        belief_users.loc[(belief_users.user_id == tweeter) & pd.isna(belief_users.pablo_score), 'pablo_score'] = ideol_samples #fill in scores
        
    # Bin user ideologies
    belief_users['ideology_bin'] = pd.cut(belief_users['pablo_score'],
                                          bins = [-50, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 50],
                                          labels = [-3, -2, -1, 0, 1, 2, 3])
    
    # Calculate belief of each individual probablistically 
    article_belief_data = article_belief_data[article_belief_data.belief == "True"] #only focus on true belief
    article_belief_data = article_belief_data[['belief_freq', 'belief', 'ideology_score']]
    article_belief_data = article_belief_data.rename(columns = {'ideology_score': 'ideology_bin'})
    belief_users = belief_users.merge(article_belief_data, on = 'ideology_bin')
    belief_users.loc[belief_users.relative_exposure_time >= intervention_time, 'belief_freq'] *= (1 - belief_reduction) #reduce belief freq for those after intervention kicks in
    belief_users['follower_belief'] = belief_users['belief_freq'].apply(lambda x: np.random.binomial(1, x)) #draw belief from binomial distribution
    
    # Filter to just believeing users and return
    belief_users = belief_users[belief_users.follower_belief == 1].copy()
    belief_users = belief_users.drop(columns = ['pablo_score', 'ideology_bin', 'belief_freq',])
    return belief_users
    

# Estimate exposure of unique users to the story, assuming the story, assuming users are x% less likely to see a tweet after time point t
def estimate_exposure_under_intervention(tweets, paired_tweets_followers, visibility_reduction, intervention_time, mean_exposure_time, sd_exposure_time):
    """
    Function that will determine estimate many unique users would be exposed to a story that is under an intevention in viewership.
    We assume people are ::odds_reduction:: less likely to see tweet sharing article
    
    OUTPUT:
    - exposed_over_time: a dataframe that contains the time series of tweets and exposure over time for that story.
    """
    
    # Filter out paired tweet-follower sets of removed tweets
    exposed_followers = paired_tweets_followers[paired_tweets_followers['tweet_id'].isin(tweets.tweet_id)].copy()
    del paired_tweets_followers
    
    # Estimate exposure time for each user
    exposure_time_offset = time_distribute_exposure(n_followers = exposed_followers.shape[0],
                                                    mean_exposure_time = mean_exposure_time,
                                                    sd_exposure_time = sd_exposure_time)
    exposed_followers['relative_exposure_time'] = exposed_followers['relative_tweet_time'] + exposure_time_offset
    exposed_followers = exposed_followers.sort_values(by = ['relative_exposure_time'])
    
    # Simulate intervention of reduced visibility
    under_intervention = exposed_followers.relative_exposure_time > intervention_time
    tweet_visible_pre = np.repeat(True, sum(~under_intervention))
    tweet_visible_post = np.random.choice([True, False], 
                                          size = sum(under_intervention), 
                                          replace = True, 
                                          p = [1-visibility_reduction, visibility_reduction])
    tweet_visible = np.append(tweet_visible_pre, tweet_visible_post)
    exposed_followers = exposed_followers[tweet_visible]    
    return exposed_followers


def time_distribute_exposure(n_followers, mean_exposure_time, sd_exposure_time):
    """
    Distribute the exposure to a tweet over the time following the tweets, so that exposure for all followers isn't instantaneous
    For now, we assume a truncated normal distribution cut off at time +0 and time +48 (the upper bound should matter)
    
    Parameters:
        n_followers (int):          number of followers for which to calcualte the exposure time delay.
        mean_expsoure_time (float): mean of truncated normal distribution.
        sd_exposure_time (float):   standard devistion of truncated normal distribution.
        
    Returns:
        added_exposure_time (float): vector that describes at what point after the tweet that the follower was exposed.
    """
    
    lower, upper = 0, 48
    a = (lower - mean_exposure_time) / sd_exposure_time
    b = (upper - mean_exposure_time) / sd_exposure_time
    exposure_dist = stats.truncnorm(a, b, loc = mean_exposure_time, scale = sd_exposure_time)
    added_exposure_time = exposure_dist.rvs(n_followers)
    return added_exposure_time

    
def match_followers_to_tweet(tweets, story_id, data_directory):
    """
    Function that will load followers that would have seen a given tweet.
    
    OUTPUT
    - matched_followers:   dataframe listing tweet ID with the set of follower IDs that could have potentially seen it (i.e., the followers of the user who tweeted it).
    """
    
    # Set place to store compiled lists
    data_subdir = data_directory + 'data_derived/interventions/compiled_exposed_follower_lists/'
    os.makedirs(data_subdir, exist_ok = True)
    
    # If file already exists, load; otherwise, compile and save.
    file_name = data_subdir + 'exposedfollowers_story' + str(story_id) + '.csv'
    col_names = ['user_id', 'follower_id', 'tweet_id', 'tweet_time', 'tweet_number', 'relative_tweet_time']
    if not os.path.exists(file_name):
    
        # Get list of unique users and tweets in this set 
        unique_users = tweets['user_id'].unique()
        tweet_list = tweets[['user_id', 'tweet_id', 'tweet_time', 'tweet_number', 'relative_tweet_time']].copy()
        del tweets
    
        # Create paired list of tweeters and followers, batching out to a temp csv   
        follower_files = os.listdir(data_directory + "data/followers/")
        for user_id in unique_users:
            regex = re.compile(r"[0-9].*_%s.csv" % user_id)
            file = list(filter(regex.match, follower_files))
            followers = load_followers(file, data_directory)
            matched_followers = pd.DataFrame({'user_id': user_id, 'follower_id': followers}, dtype = 'int64')
            matched_followers = matched_followers.merge(tweet_list, how = "left", on = "user_id")
            matched_followers = matched_followers[col_names] #make sure columns are in right order
            matched_followers.to_csv(file_name, mode = "a", index = False, header = False)
            del matched_followers, followers
        
    # Load in compiled list, sort, and return
    matched_followers = pd.read_csv(file_name, header = None, names = col_names)
    matched_followers['tweet_time'] = pd.to_datetime(matched_followers['tweet_time'], format = '%Y-%m-%d %H:%M:%S%z')
    matched_followers = matched_followers.sort_values(by = ['tweet_number'])
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



# Execution of actual intervention simulation
if __name__ == "__main__":
    
    ####################
    # Load fake news tweets
    ####################
    # Determine which article we are simulating
    which_story = int(sys.argv[1]) #get from command line

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
    fm_stories = fm_stories[fm_stories > 10] #we don't use articles 1-10 in our study
    story = int(fm_stories[which_story])
    
    
    ####################
    # Load user ideologies
    ####################
    # Load ideology scores--of both followers and tweeters (since tweeters can also be followers)
    follower_ideologies = pd.read_csv(data_directory + "data_derived/ideological_scores/cleaned_followers_ideology_scores.csv",
                                      dtype = {'user_id': 'int64', 'pablo_score': float})
    follower_ideologies = follower_ideologies.rename(columns = {'user_id': 'follower_id'})
    follower_ideologies = follower_ideologies.drop(columns = ['accounts_followed']) 
    tweeter_ideologies = labeled_tweets[['user_id', 'user_ideology']]
    tweeter_ideologies = tweeter_ideologies.rename(columns = {'user_id': 'follower_id', 'user_ideology': 'pablo_score'})
    tweeter_ideologies = tweeter_ideologies[~pd.isna(tweeter_ideologies['pablo_score'])]
    follower_ideologies = follower_ideologies.append(tweeter_ideologies, ignore_index = True)
    follower_ideologies = follower_ideologies.drop_duplicates()
    del tweeter_ideologies, labeled_tweets
    
    # Load inferred follower distribution shape for each tweeter
    follower_distributions = pd.read_csv(data_directory + "data_derived/ideological_scores/estimated_ideol_distributions/follower_ideology_distribution_shapes.csv", dtype = {'user_id': 'int64'})
    
    ####################
    # Load and shape belief data
    ####################
    # Load belief data, turn into dataframe
    belief_file = data_directory + '/data/article_belief/response_distribution.p'
    belief_data = pickle.load( open(belief_file, "rb") )
    belief_data = pd.DataFrame.from_dict(belief_data, orient = 'index')
    belief_data['total_article_number'] = belief_data.index
    
    # Drop unneccessary columns 
    belief_data = belief_data.drop(columns = ['date', 'art_type'])
    belief_data = belief_data.loc[:, ~belief_data.columns.str.match('.*response_count')] #drop response count columns
    
    # Melt and form new columns for ideology and belief
    belief_data = pd.melt(belief_data, id_vars = ['total_article_number', 'fc_rating', 'mean_fc_likert', 'total_responses'])
    belief_data['belief'] = belief_data['variable'].str.extract(r'_([a-z])$')
    belief_data['ideology'] = belief_data['variable'].str.extract(r'^(.*)_[a-z]$')
    belief_data['ideology'] = belief_data['ideology'].str.replace(': Middle of the road', '') #remove extra text from moderates
    belief_data = belief_data[belief_data.ideology != "Haven't thought much about it"]
    belief_data['ideology_score'] = belief_data['ideology'].map({'Extremely Liberal': -3, 'Liberal': -2, 'Slightly Liberal': -1,
                                                                 'Moderate':0, 
                                                                 'Slightly Conservative': 1, 'Conservative': 2, 'Extremely Conservative': 3})
    
    # Clean up dataframe
    belief_data = belief_data.drop(columns = ['variable'])
    belief_data = belief_data.rename(columns = {'value': 'belief_freq'})
    belief_data['belief'] = belief_data['belief'].map({'t': 'True', 'f': 'False', 'c': 'Unsure'})
    
    # Grab only belief for this article
    story_belief_data = belief_data[belief_data.total_article_number == story].copy()
    
    
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
    
    # Create list of tweets that includes list of followers per tweets, and simulated
    tweets_with_followers = match_followers_to_tweet(tweets = story_tweets, 
                                                     story_id = story, 
                                                     data_directory = data_directory)
    
    # Load retweet network
    RT_network = pd.read_csv(data_directory + "data_derived/networks/specific_article_networks/article" + str(story)+ "_edges.csv",
                              dtype = {'Source': 'int64', 'Target': 'int64', 'source_tweet_id': 'int64', 'target_tweet_id': 'int64'})
    
        
    ####################
    # Loop through replicate simulations
    ####################
    # Prep directory for output
    sub_dir = "{}reduce_sharing{}_visibility{}_belief{}/".format(outpath, str(sharing_reduction), str(visibility_reduction), str(belief_reduction))
    os.makedirs(sub_dir, exist_ok = True)
    
    # Loop through replicate simulations
    all_intervention_tweets = []
    all_exposure_timeseries = []
    
    # First, no intervention
    noint_tweets, noint_exposure_time = simulate_intervention(tweets = story_tweets, 
                                                              paired_tweets_followers = tweets_with_followers, 
                                                              ideologies = follower_ideologies,
                                                              follower_ideol_distributions = follower_distributions,
                                                              RT_network = RT_network,
                                                              article_belief_data = story_belief_data,
                                                              sharing_reduction = 0, #NO INTERVENTION
                                                              visibility_reduction = 0, #NO INTERVENTION
                                                              belief_reduction = 0, #NO INTERVENTION
                                                              intervention_time = 0, 
                                                              replicate_number = -1,
                                                              mean_time_to_exposure = 1,
                                                              sd_time_to_exposure = 2)
    noint_tweets['simulation_type'] = noint_exposure_time['simulation_type'] = "no intervention"
    noint_tweets['visibility_reduction'] = noint_exposure_time['visibility_reduction'] = visibility_reduction #label which simulation run this baseline is part of
    noint_tweets['sharing_reduction'] = noint_exposure_time['sharing_reduction'] = sharing_reduction #label which simulation run this baseline is part of
    all_intervention_tweets.append(noint_tweets)
    all_exposure_timeseries.append(noint_exposure_time)
    del noint_tweets, noint_exposure_time
    
    
    # Next, intervention in effect
    for intervention_time in range(0, 13):
        for i in np.arange(n_replicates):
            replicate_tweets, replicate_exposure_time = simulate_intervention(tweets = story_tweets, 
                                                                              paired_tweets_followers = tweets_with_followers, 
                                                                              ideologies = follower_ideologies,
                                                                              follower_ideol_distributions = follower_distributions,
                                                                              RT_network = RT_network,
                                                                              article_belief_data = story_belief_data,
                                                                              sharing_reduction = sharing_reduction,
                                                                              visibility_reduction = visibility_reduction, 
                                                                              belief_reduction = belief_reduction,
                                                                              intervention_time = intervention_time, 
                                                                              replicate_number = i,
                                                                              mean_time_to_exposure = 1,
                                                                              sd_time_to_exposure = 2)
            replicate_tweets['simulation_type'] = replicate_exposure_time['simulation_type'] = "intervention"
            all_intervention_tweets.append(replicate_tweets)
            all_exposure_timeseries.append(replicate_exposure_time)
            del replicate_tweets, replicate_exposure_time
        
    # Bind together and save
    all_intervention_tweets = pd.concat(all_intervention_tweets)
    all_exposure_timeseries = pd.concat(all_exposure_timeseries)
    all_intervention_tweets.to_csv(sub_dir + 'article' + str(story) + "_intervention_tweets.csv", index = False)
    all_exposure_timeseries.to_csv(sub_dir + 'article' + str(story) + "_intervention_exposetime.csv", index = False)