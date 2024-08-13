#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 16:42:32 2024

@author: 4vt
"""

from optimize_dinosaur import pipelines
pipemethods = ('get_params', 'set_params', 'get_metrics', 'setup_workspace', 'run_job')
pipeline_objects = [p() for p in pipelines.__dict__.values() if all(hasattr(p, n) for n in pipemethods)]

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--task', action = 'store', choices = ['initialize',
                                                                 'initialize_workspace',
                                                                 'genetic_optimization',
                                                                 'initial_job',
                                                                 'genetic_job'],
                    help = 'Which task to perform')
parser.add_argument('-i', '--index', action = 'store', type = int, required = False, default = -1,
                    help = 'The slurm array index')
parser.add_argument('-d', '--directory', action = 'store', required = True,
                    help = 'The directory containing input data')
parser.add_argument('-p', '--pipeline', action = 'store', required = True, choices = [p.name for p in pipeline_objects],
                    help = 'The pipeline to work on')
parser.add_argument('-n', '--N_jobs', action = 'store', type = int, default = 95,
                    help = 'The number of genetic optimization jobs to run')
args = parser.parse_args()
import os
args.directory = os.path.abspath(args.directory)

from optimize_dinosaur.initial_setup import make_workspace, initial_slurm_array_submission
from optimize_dinosaur.initial_trials import run_initial_job
from optimize_dinosaur.optimizer_job import run_optimizer_job, genetic_slurm_array_submission

pipeline = next(p for p in pipeline_objects if p.name == args.pipeline)

if args.task == 'initialize_workspace':
    make_workspace(args.directory, pipeline)
    
elif args.task == 'initialize':
    make_workspace(args.directory, pipeline)
    initial_slurm_array_submission(args.directory, pipeline)

else:
    if os.path.split(args.directory)[-1] == f'{pipeline.name}_optimization':
        print('It looks like you passed a results directory to --directory, please always use the same argument to --directory',
              flush = True)
    os.chdir(os.path.join(args.directory, f'{pipeline.name}_optimization'))
    
    if args.task == 'genetic_optimization':
        genetic_slurm_array_submission(pipeline, args.directory, args.N_jobs)
    
    elif args.task == 'initial_job':
        run_initial_job(args.index, pipeline)
    
    elif args.task == 'genetic_job':
        run_optimizer_job(args.index, pipeline)

