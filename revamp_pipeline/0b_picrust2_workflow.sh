#!/bin/bash
#SBATCH --job-name=picrust2         # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --cpus-per-task=6
#SBATCH --ntasks=1
#SBATCH --mem=120gb                     # Job memory request
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=testserial.%j.out    # Standard output log
#SBATCH --error=testserial.%j.err     # Standard error log

cd $SLURM_SUBMIT_DIR

module load PICRUSt2/2.3.0

picrust2_pipeline.py -s dna-sequences.fasta -i rarefied_no-mitochondria-no-chloroplast-blank-filtered-YS-only-Sub.biom -o picrust2_out_rarefied_table_no_sub_pipeline -p 6
