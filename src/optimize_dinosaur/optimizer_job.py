#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 14:24:48 2024

@author: 4vt
"""

def breeding_population(outcomes):
    from sortedcontainers import SortedList
    
    depth_index = SortedList(zip(outcomes['quant_depth'], outcomes.index))
    mre_index = SortedList(zip(outcomes['mean_relative_error'], outcomes.index), key = lambda x: -x[0])
    
    breeding_pop = []
    for i, depth, mre in zip(outcomes.index, 
                             outcomes['quant_depth'], 
                             outcomes['mean_relative_error']):
        depth_dominated = depth_index[depth_index.bisect_right((depth,)):]
        if depth_dominated:
            depth_dominated = set(o[1] for o in depth_dominated) if type(depth_dominated[0]) == tuple else set([depth_dominated[1]])
        else:
            depth_dominated = set([])
        mre_dominated = mre_index[mre_index.bisect_right((mre,)):]
        if mre_dominated:
            mre_dominated = set(o[1] for o in mre_dominated) if type(mre_dominated[0]) == tuple else set([mre_dominated[1]])
        else:
            mre_dominated = set([])
        if len(depth_dominated.intersection(mre_dominated)) <= 1:
            breeding_pop.append(i)
    
    return breeding_pop

def run_job(sarry_i):
    import os
    
    import pandas as pd
    import numpy as np
    
    from optimize_dinosaur.shared_data import params
    from optimize_dinosaur.processing_tools import run_job
    
    rng = np.random.default_rng(os.getpid())
    
    #read finished run data and attempted solutions
    attempts = pd.read_csv('attempted_solutions.tsv', sep = '\t')
    attempts = set(zip(*[attempts[c] for c in attempts.columns]))
    
    outcomes = pd.read_csv('outcomes.tsv', sep = '\t')
        
    #find breeding populaiton
    breeding_pop = breeding_population(outcomes)
    
    #select parents
    p1, p2 = rng.choice(breeding_pop, 2, replace = False)
    p1 = outcomes.loc[p1, params.keys()]
    p2 = outcomes.loc[p2, params.keys()]
    
    parent_params = {k:p for k,p in zip(params.keys(), zip(p1, p2))}
    
    #make child
    choices = rng.choice((0,1), len(params.keys()), replace = True)
    new_params = {k:parent_params[k][c] for k,c in zip(params.keys(), choices)}
    for k in new_params.keys():
        if type(new_params[k]) != str:
            new_params[k] = new_params[k].item()

    while tuple(new_params.values()) in attempts:
        mutate_param = str(rng.choice(list(params.keys()), 1)[0])
        new_params[mutate_param] = rng.choice(params[mutate_param], 1)[0].item()

    #run job
    run_job(new_params)


