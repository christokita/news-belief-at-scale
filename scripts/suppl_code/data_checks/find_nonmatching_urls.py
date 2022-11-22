#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 13:36:56 2020

@author: ChrisTokita
"""
import pandas as pd


data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/"

tweets = pd.read_csv(data_directory + "data_derived/tweets/parsed_tweets.csv")
unassigned_tweets = tweets[pd.isna(tweets['total_article_number'])] #get unassigned tweets
unassigned_tweets = unassigned_tweets[unassigned_tweets['url_count'] == 1]