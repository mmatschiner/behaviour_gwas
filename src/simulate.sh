# m_matschiner Tue Mar 10 10:48:07 CET 2020

# Set the account.
acct=nn9883k

# Make the results directory.
mkdir -p ../res/msprime

# Make the log directory.
mkdir -p ../log/msprime

# Set the tree file.
species_tree=../data/trees/b1.tre

# Set the table with a subset of species.
species_table=../data/tables/species.txt

# Set simulation parameters.
n_inds=1
pop_size=80000
ploidy=1
chr_length=100000
gene_flow_period=1000000

# Perform simulations with different rate matrices.
for migration_matrix in ../data/migration_matrices/*.txt
do
    # Set the prefix.
    matrix_id=`basename ${migration_matrix%.txt}`
    prefix=../res/msprime/${matrix_id}

    # Set the log file.
    log=../log/msprime/${matrix_id}.txt
    rm -f ${log}

    # Simulate tree sequences with msprime.
    sbatch -A ${acct} -o ${log} simulate.slurm ${n_inds} ${pop_size} ${ploidy} ${chr_length} ${gene_flow_period} ${species_tree} ${migration_matrix} ${species_table} ${prefix}
    exit
done
