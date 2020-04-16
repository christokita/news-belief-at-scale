# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 13:54:07 2019

@author: 15037
"""

import csv

csv_file = './Original_Data/familiesbelongtogether_tweets.csv'


#Pull Out the tweet_id, user_id, type_tweet (Retweet,Quoted,Original),
#Original_text , Retweet_text

with open(csv_file, 'r') as csvfile:
    reader = csv.DictReader(csvfile, delimiter = '|')
    headers = reader.fieldnames
    
print(headers)

tweet_data = []

with open(csv_file, 'r',encoding="UTF-8") as csvfile:
    reader = csv.DictReader(csvfile, delimiter = '|')
    for row in reader:
        tweet_data.append({'tweet_id': row['tweet_id'],
                           'user_id': row['user_id'],
                           'original_tweet_text': row['tweet_text'],
                           'retweet_status': row['retweet_status'],
                           'original_follower_count' : row['num_followers'],
                           'full_json': row['full_json']})
type(tweet_data[0])


Quote_tweet_diff_1 = []
Quote_tweet_diff_2 = []
for i in range(0,len(tweet_data)):
    json = eval(tweet_data[i]['full_json'])
    if tweet_data[i]['retweet_status'] == 'true':
        r_stats = json['retweeted_status']
        dic = {'type_tweet': 'Retweet', 'retweet_quote_text': r_stats['text']}
        tweet_data[i] = {**tweet_data[i], **dic}
        tweet_data[i]['original_tweet_text'] = 'NA' 
    else:
        if json['is_quote_status'] == True:
            try:
                q_stats = json['quoted_status']
                try:
                    dic = {'type_tweet': 'Quoted', 'retweet_quote_text': q_stats['text']}
                    tweet_data[i] = {**tweet_data[i], **dic}
                except:
                    dic = {'type_tweet': 'Quoted', 'retweet_quote_text': q_stats['full_text']}
                    tweet_data[i] = {**tweet_data[i], **dic}
            except:
                dic = {'type_tweet': 'Quoted', 'retweet_quote_text': 'NA'}
                tweet_data[i] = {**tweet_data[i], **dic}
        else:
            dic = {'type_tweet': 'Original', 'retweet_quote_text': 'NA'}
            tweet_data[i] = {**tweet_data[i], **dic}
        

for i in range(0,len(tweet_data)):
    del(tweet_data[i]['full_json'])
    del(tweet_data[i]['retweet_status'])

with open('Generated_Data\Edited_tweets_families_belong_together.csv', 'w',encoding='UTF-8') as csvfile:
    fieldnames = ['tweet_id',
                  'user_id',
                  'type_tweet',
                  'original_tweet_text',
                  'original_follower_count', 
                  'retweet_quote_text']
    writer = csv.DictWriter(csvfile, lineterminator='\n', fieldnames=fieldnames)
    writer.writeheader()
    for i in range(0,len(tweet_data)):
        writer.writerow(tweet_data[i])





        