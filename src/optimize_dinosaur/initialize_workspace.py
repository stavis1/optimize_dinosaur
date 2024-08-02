#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:06:37 2024

@author: 4vt
"""

def make_workspace(target, pipeline):
    import os
    import pandas as pd

    os.chdir(target)
    params = pipeline.get_params()
    
    #set up results files
    empty_params = pd.DataFrame({k:[] for k in params.keys()})
    empty_params.to_csv('attempted_solutions.tsv', sep = '\t', index = False)
    
    outcomes = pd.DataFrame({m:[] for m in list(pipeline.get_metrics().keys()) + ['runtime']})
    empty_params = pd.concat([empty_params, outcomes])
    empty_params.to_csv('outcomes.tsv', sep = '\t', index = False)
    
    #add initial trials jobs file
    def make_job(p,i):
        return [v[0] if param != p else v[i] for param,v in params.items()]
    
    jobs = pd.DataFrame([make_job(p, i) for p in params.keys() for i in range(len(params[p]))])
    jobs.columns = params.keys()
    jobs = jobs.drop_duplicates()   
    jobs.to_csv('initial_trials.tsv', sep = '\t', index = False)

    #run pipeline specific setup
    pipeline.setup_workspace()


