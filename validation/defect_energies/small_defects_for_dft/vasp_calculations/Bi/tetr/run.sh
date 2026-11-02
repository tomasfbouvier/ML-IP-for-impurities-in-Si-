#!/bin/bash
#SBATCH --job-name=example
#SBATCH --account=djurabek
#SBATCH --partition=medium
#SBATCH --time=01:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=384 --cpus-per-task=1  # The product should be 384

# ASE vasp environment variables

srun /projappl/project_2019925/vasp-env/view/bin/vasp_std
