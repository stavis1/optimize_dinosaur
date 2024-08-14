#!/bin/bash

#SBATCH -A ACF-UTK0011
#SBATCH -p campus
#SBATCH --qos=campus
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -J setup
#SBATCH --output=setup_out_%j_%a.log
#SBATCH --error=setup_err_%j_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stavis@vols.utk.edu

#all parameters should be passed to the sbatch call wrapped in quotes
python -m optimize_dinosaur $1
