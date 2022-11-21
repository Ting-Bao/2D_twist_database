#!/bin/bash
#PBS -N opmxband
#PBS -l nodes=1:ppn=64
#PBS -l walltime=96:00:00

cd ${PBS_O_WORKDIR}
nodes_num=$(cat ${PBS_NODEFILE} | wc -l)
#module load intel/18u2
#module load intel/19u5
#sleep 20h


source /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/bin/mpivars.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin

#for i in {0..29}; do
#for j in {0..9}; do

#cd /home/xyz/runzhang/2d-bi2se3_train_set_vdw12/bi2se3-2l-soc-300-xy0.1-z0.2/config/$i"_"$j/

#mpirun -np ${nodes_num} /home/liyang1/Software/CalcProg/OpenMX/platform/w001/openmx3.9-patched_18u2/source/openmx openmx_in.dat > openmx.std
#mpirun -np 16 /home/liyang1/Software/CalcProg/OpenMX/platform/w001/openmx3.9-patched_18u2/source/openmx openmx_in.dat -nt 4 > openmx.std

mpirun -np 64 /home/xurz/bin/openmx_399_19u5 openmx_in.dat  > openmx.std
cat openmx.out >> openmx.scfout

#done
#done
