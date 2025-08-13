#!/bin/bash

#SBATCH -p cpu
#SBATCH --mem=50g
#SBATCH -c 8
#SBATCH --time=60:00:00
#SBATCH --job-name=exitron_pipe
#SBATCH --output=%A.out
#SBATCH --error=%A.err

python optimized_exitron_pipe.py
