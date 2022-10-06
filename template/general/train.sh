echo '[begin gen graph]' >> sysinfo.txt
date >> sysinfo.txt
nvidia-smi >> sysinfo.txt
conda activate 20220907&& which python >> sysinfo.txt
pythonpath=/root/miniconda3/envs/20220907/bin/python
e3path=/root/baot/DeepH-E3-220907/
${pythonpath} ${e3path}/deephe3-train.py ./train.ini -n 8 >> sysinfo.txt
date >> sysinfo.txt
echo '[end gen graph]' >> sysinfo.txt
