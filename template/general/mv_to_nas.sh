CalcPath=CONTENT2
StorePath=CONTENT3
Name=CONTENT1

mkdir ${StorePath}/processed/${Name}
mkdir ${StorePath}/graph/${Name}
mkdir ${StorePath}/config_reserved/${Name}
mkdir ${StorePath}/config_reserved/${Name}/config

mv ${CalcPath}/processed ${StorePath}/processed/${Name}/
mv ${CalcPath}/graph ${StorePath}/graph/${Name}/
mv ${CalcPath}/result ${StorePath}/graph/${Name}/
mv ${CalcPath}/output ${StorePath}/config_reserved/${Name}/output
mv ${CalcPath}/config ${StorePath}/config_reserved/${Name}/


#for i in {0..575}; do
#for j in 0; do
#mkdir ${StorePath}/config_reserved/${Name}/config/$i"_"$j
#mv ${CalcPath}/config/$i"_"$j/openmx.out ${StorePath}/config_reserved/${Name}/config/$i"_"$j/openmx.out
#done
#done

# make sure file stored correctly
# then remove the config file
#if [ -f  /home/xyz/nas_disk/baot/dataset/processed/${Name}/processed/575_0/rh.h5 ]
#then
#    rm -rf ${CalcPath}/config
#fi

date >> ${StorePath}/../finish.txt
echo "${Name}" >> ${StorePath}/../finish.txt
