#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:06:37 2024

@author: 4vt
"""

def make_workspace(target):
    import os
    import subprocess
    import pandas as pd
    from optimize_dinosaur.shared_data import params    
    os.chdir(target)
    
    if not os.path.exists('Dinosaur.jar'):
        subprocess.run('wget https://github.com/fickludd/dinosaur/releases/download/1.2.0/Dinosaur-1.2.0.free.jar -O Dinosaur.jar',
                       shell = True)

    #set up results files
    empty_params = pd.DataFrame({k:[] for k in params.keys()})
    empty_params.to_csv('attempted_solutions.tsv', sep = '\t', index = False)
    
    outcomes = pd.DataFrame({'mean_relative_error':[],
                             'quant_overlap_jaccard':[],
                             'quant_depth':[]})
    empty_params = pd.concat([empty_params, outcomes])
    empty_params.to_csv('outcomes.tsv', sep = '\t', index = False)
    
    #add initial trials jobs file
    def make_job(p,i):
        return [v[0] if param != p else v[i] for param,v in params.items()]
    
    jobs = pd.DataFrame([make_job(p, i) for p in params.keys() for i in range(len(params[p]))])
    jobs.columns = params.keys()
    jobs = jobs.drop_duplicates()
    
    jobs.to_csv('initial_trials.tsv', sep = '\t', index = False)
