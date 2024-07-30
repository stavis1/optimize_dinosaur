#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:25:11 2024

@author: 4vt
"""

def run_trial(sarray_i):
    import pandas as pd

    from optimize_dinosaur.processing_tools import run_job
    
    jobs = pd.read_csv('initial_trials.tsv', sep = '\t')
    job = {col:param for col,param in zip(jobs.columns, jobs.iloc[sarray_i,:])}
    run_job(job)
    
 