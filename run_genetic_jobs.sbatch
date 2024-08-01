#!/bin/bash

#SBATCH -A ACF-UTK0011
#SBATCH -p campus
#SBATCH --qos=campus
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -J genetic_runs
#SBATCH --output=initial_trials_out_%j_%a.log
#SBATCH --error=initial_trials_err_%j_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stavis@vols.utk.edu

conda run -p ~/optimize_dinosaur/env/od_env python -m optimize_dinosaur -t genetic_optimization -i $SLURM_ARRAY_TASK_ID


