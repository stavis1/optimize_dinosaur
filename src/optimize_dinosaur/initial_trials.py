#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:25:11 2024

@author: 4vt
"""

def run_initial_job(sarray_i, pipeline):
    import pandas as pd

    jobs = pd.read_csv('initial_trials.tsv', sep = '\t', dtype = str, keep_default_na=False)
    job = {col:param for col,param in zip(jobs.columns, jobs.iloc[sarray_i,:])}
    pipeline.run_job(job)
    
 