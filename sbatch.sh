#!/bin/bash                                                                       
#SBATCH -c 2
#SBATCH -N 1
#SBATCH -A scope
#SBATCH -p scope-shared
#SBATCH -t 15-0:00
#SBATCH --mem 5G
#SBATCH -J MSt_MAGs
#SBATCH -o slurm/mst_mags.%A.%a.out
#SBATCH -e slurm/mst_mags.%A.%a.out

#export PATH=/home/jmeppley/opt/conda/bin:$PATH
#source activate ./env
. /home/jmeppley/.bash_conda
conda activate ./env

snakemake \
    --config max_threads=9 \
    -j 999 -p -k --nt --local-cores 2 \
    --latency-wait 30 \
    --use-conda \
    --cluster-config cluster.yaml \
    --cluster "sbatch --mem {cluster.mem} -c {threads} -t {cluster.t} -A scope -p scope-shared" \


