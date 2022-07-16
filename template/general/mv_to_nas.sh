CalcPath=CONTENT2
StorePath=CONTENT3
Name=CONTENT1

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

echo "${Name}" >> ${StorePath}/../finish.txt
