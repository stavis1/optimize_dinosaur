#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 14:24:48 2024

@author: 4vt
"""

def breeding_population(outcomes, pipeline):
    from sortedcontainers import SortedList
    metric1, metric2 = pipeline.get_metrics().items()
    
    m1_index = SortedList(zip(outcomes[metric1[0]], outcomes.index), key = lambda x: x[0]*metric1[1])
    m2_index = SortedList(zip(outcomes[metric2[0]], outcomes.index), key = lambda x: x[0]*metric2[1])
    
    breeding_pop = []
    for i, m1, m2 in zip(outcomes.index, 
                             outcomes[metric1[0]], 
                             outcomes[metric2[0]]):
        m1_dominated = m1_index[m1_index.bisect_right((m1,)):]
        if m1_dominated:
            m1_dominated = set(o[1] for o in m1_dominated) if type(m1_dominated[0]) == tuple else set([m1_dominated[1]])
        else:
            m1_dominated = set([])
        m2_dominated = m2_index[m2_index.bisect_right((m2,)):]
        if m2_dominated:
            m2_dominated = set(o[1] for o in m2_dominated) if type(m2_dominated[0]) == tuple else set([m2_dominated[1]])
        else:
            m2_dominated = set([])
        if len(m1_dominated.intersection(m2_dominated)) <= 1:
            breeding_pop.append(i)
    
    return breeding_pop

def run_optimizer_job(sarry_i, pipeline):
    import os
    
    import pandas as pd
    import numpy as np
    
    rng = np.random.default_rng(os.getpid())
    params = pipeline.get_params()
    
    #read finished run data and attempted solutions
    attempts = pd.read_csv('attempted_solutions.tsv', sep = '\t')
    attempts = set(zip(*[attempts[c] for c in attempts.columns]))
    
    outcomes = pd.read_csv('outcomes.tsv', sep = '\t')
        
    #find breeding population
    breeding_pop = breeding_population(outcomes, pipeline)
    
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
            
    #muatate to ensure solution uniqueness
    while tuple(new_params.values()) in attempts:
        mutate_param = str(rng.choice(list(params.keys()), 1)[0])
        new_params[mutate_param] = rng.choice(params[mutate_param], 1)[0].item()

    #run job
    pipeline.run_job(new_params)


def genetic_slurm_array_submission(pipeline, target, N):
    import subprocess

    with open('init_run_script.sbatch', 'w') as sbatch:
        sbatch.write('\n'.join(['#!/bin/bash',
                                '#SBATCH -A ACF-UTK0011',
                                '#SBATCH -p campus',
                                '#SBATCH --qos=campus',
                                f'#SBATCH -t {pipeline.timeout}',
                                '#SBATCH --nodes=1',
                                f'#SBATCH -c {pipeline.cores}',
                                f'#SBATCH --mem={pipeline.memory}g',
                                f'#SBATCH -J {pipeline.name}',
                                f'#SBATCH --output=out_{pipeline.name}_%j_%a.log',
                                f'#SBATCH --error=err_{pipeline.name}_%j_%a.log',
                                '#SBATCH --mail-type=ALL',
                                '#SBATCH --mail-user=stavis@vols.utk.edu',
                                f'#SBATCH --array=0-{N-1}\n',]))
        sbatch.write(' '.join(['python -m optimize_dinosaur',
                               '-t genetic_job',
                               f'-d {target}',
                               f'-p {pipeline.name}',
                               '-i $SLURM_ARRAY_TASK_ID']))
    
    subprocess.run('sbatch init_run_script.sbatch', shell = True)

