#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 15:14:12 2020

@author: ChrisTokita

SCRIPT:
Determine who shared fake news articles, who was potentially exposed to it, and who fake news tweeters are following. 
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np
import re
import os
import math
import multiprocessing as mp


####################
# Functions for parsing in parallel
####################
# Function to parse a portion of exposed followers
def parse_users(i, ids, directory, step_size, data_type):
    
    # Take chunk of FM article tweeters and combine followers. There are 26,163 users who tweeted a FM article
    ids = ids[ (step_size*i) : (step_size*(i+1)) ]

    # Get appropriate list of files for friend/follower data
    if data_type == "followers":
        user_files = os.listdir(directory + "data/followers/")
    elif data_type == "frineds":
        user_files = os.listdir(directory + "data/friends/")
    user_files = [file for file in user_files if re.match('[0-9]', file)] #filter out hidden copies of same files
    
    # Loop through user IDs and get friends/followers of that user
    users = np.array([], dtype = str)
    no_user_data = np.array([], dtype = str)
    for user_id in ids:
        regex = re.compile(r"[0-9].*_%s.csv" % user_id)
        file = list(filter(regex.match, user_files))
        try:
            if len(file) > 1:
                print("WARNING: user_id = %d matches multiple follower list files." % user_id)
            user_list = np.genfromtxt(data_directory + "data/followers/" + file[0], dtype = str)
            user_list = user_list[1:len(user_list)] #remove header, will raise error if empty
            users = np.append(users, user_list)
            users = np.unique(users)
        except:
            no_user_data = np.append(no_user_data, user_id)
            
    # Write to file  
    chunk_label = str(i).zfill(2)      
    users = pd.DataFrame(users, columns = ['user_id'], dtype = str)
    no_user_data = pd.DataFrame(no_user_data, columns = ['user_id'], dtype = str)
    if data_type == "followers":
        users.to_csv(directory + "data_derived/followers/processed_exposed_followers/followers_exposed_fakenews_" + chunk_label + ".csv", index = False)
        no_user_data.to_csv(directory + "data_derived/followers/nofollowers_fm_tweeters/nofollowers_fm_tweeters_" + chunk_label + ".csv", index = False)
    elif data_type == "friends":
        users.to_csv(directory + "data_derived/followers/processed_exposed_followers/followers_exposed_fakenews_" + chunk_label + ".csv", index = False)
        no_user_data.to_csv(directory + "data_derived/followers/nofollowers_fm_tweeters/nofollowers_fm_tweeters_" + chunk_label + ".csv", index = False)
    
   
####################
# Main script
####################
if __name__ == "__main__":
    
    # high level directory (external HD or cluster storage)
    data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/"
    
    ############################## Parse ##############################
    
    ####################
    # Parse articles to get false news articles
    ####################
    # Get article IDs of false/misleading articles (as evaluated by fact checkers)
    news_evaluations = pd.read_csv(data_directory + "data/articles/evaluations.csv")
    fakenews_ids = news_evaluations["article_num"][news_evaluations['mode of FC'] == "FM"]
    
    # Get URL set for FM articles
    articles = pd.read_csv(data_directory + "data/articles/daily_articles.csv")
    fm_articles = articles[articles['total article number'].isin(fakenews_ids)]
    
    # Load and parse tweets according to full and shortened URLs
    tweets = pd.read_csv(data_directory + "data_derived/tweets/parsed_tweets.csv")
    fm_tweets = pd.DataFrame(columns = tweets.columns)
    for j in range(len(fm_articles)):
        link = fm_articles['link'].iloc[j]
        shortlink = fm_articles['short link'].iloc[j]
        try:
            has_full_link = tweets['urls_expanded'].str.contains(link)
            has_short_link = tweets['urls_expanded'].str.contains(shortlink)
        except:
            has_short_link = pd.Series(np.repeat(False, tweets.shape[0]))
        has_link = has_full_link | has_short_link #boolean operator to find which indices have one of the two possible links
        fm_tweets = fm_tweets.append(tweets[has_link])  
    
    
    ####################
    # Determine which users tweeted false articles
    ####################
    # Get user IDs of tweeters of FM articles, 
    fm_tweeters = pd.DataFrame(fm_tweets['user_id'], columns = ['user_id'], dtype = str)
    fm_tweeters = fm_tweeters.drop_duplicates()
    fm_tweeters = fm_tweeters.sort_values(by = ['user_id'])
    
    # Write FM tweets and tweeters to file
    fm_tweets.to_csv(data_directory + "data_derived/tweets/FM_tweets.csv", index = False)
    
    
    ####################
    # Parse followers to FM tweets and friends of FM tweeters
    ####################
    step_size = math.ceil( len(fm_tweeters)/100 ) #determine the proper step size to get through all FM tweeters
    pool = mp.Pool(mp.cpu_count())
    for i in range(100):
        pool.apply(parse_users, args = (i, fm_tweeters['user_id'], data_directory, step_size, "followers"))
        pool.apply(parse_users, args = (i, fm_tweeters['user_id'], data_directory, step_size, "friends"))
    pool.close()
    pool.join()



    ############################## Count ##############################
    
    # Function to compile processed user data files
    def compile_userdata(path_to_data, processed_files, tweeter_set = None):
        compiled_users = pd.DataFrame(columns = ['user_id'], dtype = str)
        for file in processed_files:
            data = pd.read_csv(path_to_data + file, dtype = str)
            compiled_users = compiled_users.append(data, ignore_index = True, sort = False)
            compiled_users = compiled_users.drop_duplicates()
            compiled_users = compiled_users.sort_values(by = ['user_id'])
            del(data)
        if tweeter_set is None:
            return compiled_users
        else:
            compiled_users_minus_tweeters = pd.DataFrame(np.setdiff1d(compiled_users, tweeter_set), columns = ['user_id'], dtype = int)
            return compiled_users, compiled_users_minus_tweeters
    
    
    ####################
    # Count up FM tweeters
    ####################
    count_fm_tweeters = fm_tweeters.shape[0]
    
    
    ####################
    # Count up unique followers exposed to FM news
    ####################
    # Get unique IDs of tweeters
    tweeters = tweets['user_id'].astype(str)
    tweeters = np.unique(tweeters)
    count_tweeters = len(tweeters)
    
    # Count up followers
    path_to_exposed = data_directory + "data_derived/followers/processed_exposed_followers/"
    exposed_files = sorted( os.listdir(path_to_exposed) )
    all_exposed, all_exposed_minus_tweeters = compile_userdata(path_to_data = path_to_exposed,
                                                               processed_files = exposed_files, 
                                                               tweeter_set = fm_tweeters)
    count_all_exposed = all_exposed.shape[0]
    count_all_exposed_minus_tweeters = all_exposed_minus_tweeters.shape[0]
    
    
    # Count FM tweeters who had no followers
    path_nofollowers_fm_tweeters = data_directory + "data_derived/followers/nofollowers_fm_tweeters/"
    exposed_nofollower_files = sorted( os.listdir(path_nofollowers_fm_tweeters) )
    no_followers = compile_userdata(path_to_data = path_nofollowers_fm_tweeters,
                                    processed_files = exposed_nofollower_files)
    count_no_followers = no_followers.shape[0]
    
    
    ####################
    # Count up friends of FM tweeters
    ####################
    # Count up friends
    path_to_fm_friends = data_directory + "data_derived/friends/processed_fm_tweeter_friends/"
    friend_files = sorted( os.listdir(path_to_fm_friends) )
    all_friends, all_friends_minus_tweeters = compile_userdata(path_to_data = path_to_fm_friends,
                                                               processed_files = friend_files, 
                                                               tweeter_set = fm_tweeters)
    count_all_friends = all_friends.shape[0]
    count_all_friends_minus_tweeters = all_friends_minus_tweeters.shape[0]
        
    # Count FM tweeters who have no friends
    path_to_nofriends_fm_tweeters = data_directory + "data_derived/friends/nofriends_fm_tweeters/"
    fm_tweeters_nofriends_files = sorted( os.listdir(path_to_nofriends_fm_tweeters) )
    no_friends = compile_userdata(path_to_data = path_to_nofriends_fm_tweeters,
                                  processed_files = fm_tweeters_nofriends_files)
    count_no_friends = no_friends.shape[0]
    
    
    ####################
    # Save
    ####################
    # Measure number of users and write out to file
    unique_users = pd.DataFrame({'user_type': ["Tweeters", "FM tweeters", 
                                               "Followers exposed to FM Articles", "Exposed followers excluding FM tweeters", "FM tweeters w/o followers",
                                               "Friends of FM tweeters", "Friends of FM tweeters excluding FM tweeters", "FM tweeters w/o friends"], 
                                 'count': [count_tweeters, count_fm_tweeters, 
                                           count_all_exposed, count_all_exposed_minus_tweeters, count_no_followers,
                                           count_all_friends, count_all_friends_minus_tweeters, count_no_friends]})
    unique_users.to_csv(data_directory + "data_derived/exposed_user_summary.csv", index = False)
    
    # Save all compiled user lists
    fm_tweeters.to_csv(data_directory + "data_derived/tweets/unique_fm_tweeters.csv", index = False)
    all_exposed.to_csv(data_directory + "data_derived/followers/exposed_fakenews_followers.csv", index = False)
    all_exposed_minus_tweeters.to_csv(data_directory + "data_derived/followers/exposed_fakenews_followers_excl_fmtweeters.csv", index = False)
    no_followers.to_csv(data_directory + "data_derived/followers/fm_tweeter_nofollowers.csv", index = False)
    all_friends.to_csv(data_directory + "data_derived/friends/friends_fm_tweeters.csv", index = False)
    all_friends_minus_tweeters.to_csv(data_directory + "data_derived/friends/friends_fm_tweeters_excl_fmtweeters.csv", index = False)
    no_friends.to_csv(data_directory + "data_derived/friends/fm_tweeter_nofriends.csv", index = False)
