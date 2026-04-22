#!/bin/bash
#SBATCH --account=djurabek
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=125

nreplica=5
dump_file=results_neb.dump
srun --oversubscribe lmp_g++_openmpi 	-partition ${nreplica}x25 -v dump_file ${dump_file} -v num_replicas ${nreplica} -in in.neb >> log.neb 

rm screen*
rm log.lammps*

tail -1 log.neb >> log2.neb
python3 parse_energy_barriers.py
python3 neb_final.py -o neb_movie -r results_neb.dump.*
#rm log*
#rm results_neb*



