#!/bin/bash

# m_matschiner Mon Dec 20 13:26:27 CET 2021

# Job name:
#SBATCH --job-name=msprime
#
# Wall clock limit:
#SBATCH --time=10:00:00
#
# Processor and memory usage:
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G

## Set up job environment:
set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variables as an error
module --quiet purge # Reset the modules to the system default

# Load modules.
module load Python/3.8.6-GCCcore-10.2.0 #add module name

# Get the command-line arguments.
num_inds=${1}
pop_size=${2}
ploidy=${3}
chr_length=${4}
gene_flow_period=${5}
species_tree=${6}
migration_matrix=${7}
species_table=${8}
prefix=${9}

# Simulate tree sequences with msprime.
python3 simulate.py --num_indiv ${num_inds} \
    --Ne ${pop_size} \
    --chr_length ${chr_length} \
    --gene_flow_period ${gene_flow_period} \
    --mut_rate 0 \
    --species_list ${species_table} \
    --ploidy ${ploidy} \
    ${species_tree} \
    ${migration_matrix} \
    ${prefix}
