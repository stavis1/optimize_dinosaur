#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 16:42:32 2024

@author: 4vt
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--task', action = 'store', choices = ['initialize_workspace', 
                                                                 'initial_trials',
                                                                 'genetic_optimization'],
                    help = 'Which task to perform')
parser.add_argument('-i', '--index', action = 'store', type = int, required = False, default = -1,
                    help = 'The slurm array index')
parser.add_argument('-d', '--directory', action = 'store', required = False, default = False,
                    help = 'The directory to initialize')
args = parser.parse_args()

from optimize_dinosaur.initialize_workspace import make_workspace
from optimize_dinosaur.initial_trials import run_trial
from optimize_dinosaur.optimizer_job import run_job

if args.task == 'initialize_workspace':
    make_workspace(args.directory)

elif args.task == 'initial_trials':
    run_trial(args.index)

elif args.task == 'genetic_optimization':
    run_job(args.index)

