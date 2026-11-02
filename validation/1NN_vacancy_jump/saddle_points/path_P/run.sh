#!/bin/bash
#SBATCH --job-name=vasp
#SBATCH --account=djurabek
#SBATCH --time=36:00:00
##SBATCH --mem-per-cpu=4G
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --output=job-%J.out
#SBATCH --error=job-%J.err
#SBATCH --ntasks-per-node=384 --cpus-per-task=1  # The product should be 384

module load python-data

# ASE vasp environment variables
export ASE_VASP_COMMAND="srun /projappl/project_2019925/vasp-env/view/bin/vasp_std"
export VASP_PP_PATH=/scratch/project_2019925

python3 run_vasp_db_loop.py
