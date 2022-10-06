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

source /home/liyang1/Software/EnvsProg/IntelLib/oneapi2022.w001/setvars.sh
cd ${CalcPath}/config/$i"_"$j/
#mkdir output
target="The calculation was normally finished."
findfile=${CalcPath}/config/$i"_"$j/openmx.std
if [`grep -c "$target" $findfile` -ne '0'];then
    echo "jump $findfile"
else
    mpirun -np 64 /home/liyang1/Software/CalcProg/OpenMX/platform/w001/openmx3.9-patched_oneapi/work/openmx openmx_in.dat  > openmx.std
fi


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

date >> ${CalcPath}/../finish.txt
echo "${Name}" >> ${CalcPath}/../finish.txt