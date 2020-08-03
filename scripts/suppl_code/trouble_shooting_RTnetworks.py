#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 13:53:12 2020

@author: ChrisTokita
"""
article_id = 28
test = tweets[tweets.total_article_number == article_id]
test = test[test.is_retweet]

x = 1

tweet = test.iloc[x]
user_id = tweet.user_id 
