#!/bin/bash
#SBATCH --job-name=vasp
#SBATCH --account=djurabek
#SBATCH --time=36:00:00
##SBATCH --mem-per-cpu=4G
#SBATCH --partition=medium
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --output=job-%J.out
#SBATCH --error=job-%J.err

#module load python-data
module load vasp/5.4.4.pl2

# ASE vasp environment variables
export ASE_VASP_COMMAND="srun vasp_std"
export VASP_PP_PATH=/appl/soft/phys/vasp

python3 run_vasp_db_loop.py
