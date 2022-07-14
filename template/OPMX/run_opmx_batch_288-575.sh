#!/bin/bash
#PBS -N NAMENAME
#PBS -l nodes=1:ppn=64
#PBS -l walltime=800:00:00

#cd ${PBS_O_WORKDIR}
#cp $PBS_NODEFILE node
nodes_num=$(cat ${PBS_NODEFILE} | wc -l)
#module load intel/18u2
module load intel/19u5
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lihe/apps/hdf5-1.12.1-parallel/lib64  #openmx_eij

# 72 144 216 288 360 432 504 576
for i in {288..575}; do
for j in 0; do

cd CONFIGFILEPATH/config/$i"_"$j/
#mkdir output

mpirun -np ${nodes_num} /home/xurz/bin/openmx_399_19u5 openmx_in.dat > openmx.std
#mpirun -np 32 /home/lihe/apps/openmx/3.9.7_19u5_Eij/openmx openmx_in.dat -nt 2 > openmx.std

cat openmx.out >> openmx.scfout
rm -r openmx_rst *.cube

done
done

# use the python on the calculation server
usePYTHON="/home/baot/bin/miniconda3/envs/deeph20220712/bin/python"
# do preprocess and delete unnecessary file 
# cd CONFIGFILEPATH
