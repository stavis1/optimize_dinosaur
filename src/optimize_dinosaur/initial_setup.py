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
    temp_dir = f'{pipeline.name}_optimization'
    files = [f for f in os.listdir() if os.path.isfile(f)]
    os.mkdir(temp_dir)
    os.chdir(temp_dir)
    for file in files:
        os.link(f'../{file}', file)
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

def initial_slurm_array_submission(target, pipeline):
    import os
    import subprocess
    import pandas as pd

    trials = pd.read_csv('initial_trials.tsv', sep = '\t')    
    
    with open('init_run_script.sbatch', 'w') as sbatch:
        sbatch.write('#!/bin/bash\n')
        sbatch.write(' '.join(['python -m optimize_dinosaur',
                               '-t initial_job',
                               f'-d {target}',
                               f'-p {pipeline.name}',
                               '-i $SLURM_ARRAY_TASK_ID']))
    
    slurm_params_file = os.path.abspath(__file__)
    slurm_params_file = os.path.split(os.path.split(slurm_params_file)[0])[0]
    slurm_params_file = os.path.join(slurm_params_file, 'slurm_params.txt')
    with open(slurm_params_file, 'r') as params_file:
        local_params = [l.strip() for l in params_file]
    
    command = ['sbatch',
               f'-t {pipeline.timeout}',
               '--nodes=1',
               f'-c {pipeline.cores}',
               f'--mem={pipeline.memory}g',
               f'-J {pipeline.name}',
               f'--output=out_{pipeline.name}_%j_%a.log',
               f'--error=err_{pipeline.name}_%j_%a.log',
               f'--array=0-{trials.shape[0]-1}']
    command += local_params + ['init_run_script.sbatch']
    print(command, flush = True)
    command = ' '.join(command)
    
    subprocess.run(command, shell = True)
