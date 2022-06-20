#!/bin/sh
#PBS -N NAMENAME
#PBS -l nodes=1:ppn=64
#PBS -l walltime=3:00:00

ulimit -s unlimited
cd ${PBS_O_WORKDIR}
nodes_num=$(cat ${PBS_NODEFILE} | wc -l)
module load intel/18u2

date > output.$PBS_JOBID
# Change it to your own MPI program
mpirun -np ${nodes_num} /home/liyang1/Software/CalcProg/VASP/platform/w001/vasp-544-patched/vasp_intel18u2/bin/vasp_ncl >> output.$PBS_JOBID
# Log stop timestamp
date >> output.$PBS_JOBID