#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:13:26 2020

@author: ChrisTokita

SCRIPT
Conduct bayesian inference on each user's distrubtion of follower ideology
"""

####################
# Load libraries and packages, set important parameters
####################
import pandas as pd
import numpy as np

import arviz as az
import matplotlib.pyplot as plt
import pymc3 as pm


# high level directory (external HD or cluster storage)
#data_directory = "/scratch/gpfs/ctokita/fake-news-diffusion/" #HPC cluster storage
data_directory = "/Volumes/CKT-DATA/fake-news-diffusion/" #external HD


####################
# Load data
####################
ideology_data = pd.read_csv(data_directory + "data_derived/ideological_scores/paired_tweeter-follower_ideology_example.csv",
                            dtype = {'user_id': object})



####################
# Build Bayesian model
####################
basic_model = pm.Model()