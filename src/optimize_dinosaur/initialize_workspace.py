#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:06:37 2024

@author: 4vt
"""

def make_workspace(target):
    import os
    import pandas as pd
    from optimize_dinosaur.shared_data import params
    
    os.chdir(target)
    empty_params = pd.DataFrame({k:[] for k in params.keys()})
    empty_params.to_csv('attempted_solutions.tsv', sep = '\t', index = False)
    outcomes = pd.DataFrame({'mean_relative_error':[],
                             'quant_overlap_jaccard':[],
                             'quant_depth':[]})
    empty_params = pd.concat([empty_params, outcomes])
    empty_params.to_csv('outcomes.tsv', sep = '\t', index = False)

