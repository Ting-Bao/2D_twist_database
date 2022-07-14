#!/bin/bash
#PBS -N CONTENT1
#PBS -l nodes=1:ppn=64
#PBS -l walltime=800:00:00

#cd ${PBS_O_WORKDIR}
#cp $PBS_NODEFILE node
nodes_num=$(cat ${PBS_NODEFILE} | wc -l)
#module load intel/18u2
module load intel/19u5
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lihe/apps/hdf5-1.12.1-parallel/lib64  #openmx_eij

CalcPath=CONTENT2
StorePath=CONTENT3
Name=CONTENT1

# 72 144 216 288 360 432 504 576
for i in {0..575}; do
for j in 0; do

cd ${CalcPath}/config/$i"_"$j/
#mkdir output

mpirun -np ${nodes_num} /home/xurz/bin/openmx_399_19u5 openmx_in.dat > openmx.std
#mpirun -np 32 /home/lihe/apps/openmx/3.9.7_19u5_Eij/openmx openmx_in.dat -nt 2 > openmx.std
cat openmx.out >> openmx.scfout
rm -r openmx_rst *.cube

done
done
# use the python on the calculation server
usePYTHON="/home/baot/bin/miniconda3/envs/deeph20220712/bin/python"

source /home/xyz/.bashrc_baot

deeph-preprocess --config ${CalcPath}/preprocess.ini
deeph-train --config ${CalcPath}/gen_graph.ini

mkdir ${StorePath}/processed/${Name}
mkdir ${StorePath}/graph/${Name}
mkdir ${StorePath}/config_reserved/${Name}

mv ${CalcPath}/processed ${StorePath}/${Name}/

# make sure file stored correctly
# then remove the config file
if [ -f  /home/xyz/nas_disk/baot/processed_dataset/${Name}/processed/575_0/rh.h5]
then
    rm -rf /home/xyz/baot/dataset/${Name}/config
fi

echo "${Name}" >> ${CalcPath}/../finish.txt
echo "${Name}" >> ${StorePath}/../finish.txt