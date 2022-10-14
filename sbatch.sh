#!/bin/bash
#SBATCH -c 2
#SBATCH -N 1
#SBATCH -A scope
#SBATCH -p scope-shared
#SBATCH -t 15-0:00
#SBATCH --mem 5G
#SBATCH -J PRGN_MAGs
#SBATCH -o slurm/prgn_mags.%A.out
#SBATCH -e slurm/prgn_mags.%A.out

# default is to run everything ("all" target)
# some versions of snakemake don't submit all jobs leading
# up to a checkpoint, so it may make sense to run
# in stages using the "up_to_checkpoints" target
TARGET=${1:all}

REPO_DIR=/home/jmeppley/repos/mags_from_separate_assemblies
. /home/jmeppley/.bash_conda
conda activate $REPO_DIR/snakemake.env

SNAKEFILE=${REPO_DIR}/Snakefile

mkdir -p slurm/snake

snakemake \
    -s ${SNAKEFILE} \
    --config \
        max_threads=9 \
        fasta_template="../assembly/{sample}/contigs.all.fasta" \
        fastq_template="../assembly/{sample}/{seq_run}.clean.fastq" \
    -j 999 -p -k --nt --local-cores 2 \
    --latency-wait 30 \
    --use-conda \
    --rerun-incomplete \
    --cluster-config cluster.yaml \
    --cluster "sbatch --mem {cluster.mem} -c {threads} -t {cluster.t} \
               -A scope -p scope-shared \
               -o slurm/snake/task.%A.out -e slurm/snake/task.%A.out" \
    ${TARGET}

