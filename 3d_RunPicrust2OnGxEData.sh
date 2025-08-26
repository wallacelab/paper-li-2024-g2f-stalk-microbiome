#! /bin/bash

# Run PICRUST2 on the GxE stalk endophyte data

# Define needed files & parameters
rep_seqs=3_GxE/3d_rep_seqs.fasta   # Representative sequences; from QIIME
abundance_table=3_GxE/3a_asv_table.rarefied.filt_for_gxe.tsv   # Tab-separated abundance table; exported in step 3a_
outdir=3_GxE/3d_picrust2_predictions  # Output directory
ncores=7  # Parallel processing cores

# Activate conda environment with PICRUST2 v2.5.3; (After 2.6, the database changed and we got very different results, so sticking with the original)
. $(conda info --root)/etc/profile.d/conda.sh 
conda activate picrust2-2.5.3

# Unzip rep seqs from original DADA2 pipeline
unzip -p 0_data_files/dada2_rep-seqs.qza 4e5ce913-2a8c-4916-8016-f8ed3b04d2c1/data/dna-sequences.fasta > $rep_seqs
 
# Run PICRUST2 Pipeline
picrust2_pipeline.py -s $rep_seqs -i $abundance_table -o $outdir -p $ncores

# Add pathway descriptions
add_descriptions.py -i $outdir/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o $outdir/pathways_out/path_abun_unstrat_descrip.tsv.gz
