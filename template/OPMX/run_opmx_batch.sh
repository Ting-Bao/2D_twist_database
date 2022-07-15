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

#CalcPath have nameid while StorePath not
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

cd ${CalcPath}
deeph-preprocess --config ${CalcPath}/preprocess.ini >> output
deeph-train --config ${CalcPath}/gen_graph.ini >> output

mkdir ${StorePath}/processed/${Name}
mkdir ${StorePath}/graph/${Name}
mkdir ${StorePath}/config_reserved/${Name}
mkdir ${StorePath}/config_reserved/${Name}/config

mv ${CalcPath}/processed ${StorePath}/processed/${Name}/processed
mv ${CalcPath}/graph ${StorePath}/processed/${Name}/graph
mv ${CalcPath}/output ${StorePath}/config_reserved/${Name}/output

for i in {0..575}; do
for j in 0; do
mkdir ${StorePath}/config_reserved/${Name}/config/$i"_"$j
mv ${CalcPath}/config/$i"_"$j/openmx.out ${StorePath}/config_reserved/${Name}/config/$i"_"$j/openmx.out
done
done

# make sure file stored correctly
# then remove the config file
if [ -f  /home/xyz/nas_disk/baot/processed_dataset/${Name}/processed/575_0/rh.h5]
then
    rm -rf ${CalcPath}/config
fi

echo "${Name}" >> ${CalcPath}/../finish.txt
echo "${Name}" >> ${StorePath}/../finish.txt