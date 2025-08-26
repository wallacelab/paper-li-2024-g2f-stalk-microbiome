#!/bin/bash

#import the raw data. 
#They are already demutiplexed Pairedend data
#Use According import format


#!/bin/bash

ls /home/hl46161/G2F_data_dada2/raw_data_fastq_results/ | awk -F_ '{print $1}' >> G2F_sample_id.txt

##

trimmed_fastq_dir=/home/hl46161/new_G2F_dada2/trimmed_g2f_data

remove_location=(GAH1 TXH1 TXH2)

##removed fastq samples from GAH1 TXH1 TXH2 location 
for file in $(ls $trimmed_fastq_dir| grep TXH2)
do
   echo $file
   rm $trimmed_fastq_dir/$file

done



#########

for file in $(ls $trimmed_fastq_dir| grep .*.fastq.gz)
do
   echo $file
   fastqc $trimmed_fastq_dir/$file  -o $trimmed_fastq_dir/trimmed_data_fastq_results/

done

#########

multiqc ./trimmed_data_fastq_results/ -o trimmed_multiqc_report/

####################


for sample in $(ls $trimmed_fastq_dir); do
       
       echo $sample
       if [[ $sample == *_"R1"_* ]];
       then
          revised_sample=`echo $sample | cut -d'_' -f 1`
          echo $revised_sample
          echo $revised_sample,$trimmed_fastq_dir/$sample,forward >> manifest_file.txt
       else
          revised_sample=`echo $sample | cut -d'_' -f 1`
          echo $revised_sample
          echo $revised_sample,$trimmed_fastq_dir/$sample,reverse >> manifest_file.txt
       fi
done


##############


#import G2F data in paired reads format
 qiime tools import \
     --type 'SampleData[PairedEndSequencesWithQuality]' \
     --input-path manifest_file.txt \
     --output-path paired-end-demux.qza \
     --input-format PairedEndFastqManifestPhred33
     

    
# # #summarize the raw sequences 

 qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv
  
 
############ import the sequences into dada for denoise  do not have to join reads into pair 

   qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 240 \
  --o-table dada2_table.qza \
  --o-representative-sequences dada2_rep-seqs.qza \
  --o-denoising-stats dada2_denoising-stats.qza
  
########### summarize the table 

  qiime feature-table summarize \
  --i-table dada2_table.qza \
  --o-visualization dada2_table.qzv \
  --m-sample-metadata-file 2018_G2F_metadata_enviromental.tsv

########## use the silva 138 classifier to classify taxa 

qiime feature-table tabulate-seqs \
  --i-data dada2_rep-seqs.qza \
  --o-visualization dada2_rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file dada2_denoising-stats.qza \
  --o-visualization dada2_denoising-stats.qzv
  
      qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --i-reads dada2_rep-seqs.qza \
  --o-classification silva-138-99_nb_DADA2_taxonomy.qza
  
  qiime metadata tabulate \
  --m-input-file silva-138-99_nb_DADA2_taxonomy.qza \
  --o-visualization silva-138-99_nb_DADA2_taxonomy.qzv
  
    qiime taxa barplot \
  --i-table dada2_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --m-metadata-file 2018_G2F_metadata_enviromental.tsv \
  --o-visualization dada2-pretrained-taxa-bar-plots.qzv
  
## filter out mitochondria and chloroplast

    qiime taxa filter-table \
  --i-table dada2_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast.qza
  
   qiime feature-table summarize \
     --i-table dada2_table-no-mitochondria-no-chloroplast.qza \
     --o-visualization dada2_table-no-mitochondria-no-chloroplast.qzv 

# visualize the table after filtration 

     qiime taxa barplot \
  --i-table dada2_table-no-mitochondria-no-chloroplast.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --m-metadata-file 2018_G2F_metadata_enviromental.tsv \
  --o-visualization dada2_table-no-mitochondria-no-chloroplast-pretrained-taxa-bar-plots.qzv 
  
  
## generate phylogentic tree 

   qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences dada2_rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree dada2-unrooted-tree.qza \
  --o-rooted-tree dada2-rooted-tree.qza

## generate the rarefaction graph

  qiime diversity alpha-rarefaction \
  --i-table dada2_table-no-mitochondria-no-chloroplast.qza \
  --p-max-depth 35000 \
  --o-visualization dada2_table-no-mitochondria-no-chloroplast-alpha-rarefaction-35000.qzv



  qiime diversity core-metrics-phylogenetic \
   --i-phylogeny dada2-rooted-tree.qza \
   --i-table dada2_table-no-mitochondria-no-chloroplast.qza  \
   --p-sampling-depth 1500 \
   --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
   --output-dir core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-corrected_metadata-1500
   


       qiime feature-table filter-samples\
  --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered.qza \
  --m-metadata-file G2F_metadata_2019_duplicate_pedigree.txt \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree.qza
  
  
  
  qiime diversity core-metrics-phylogenetic \
   --i-phylogeny dada2-rooted-tree.qza \
   --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree.qza  \
   --p-sampling-depth 1500 \
   --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
   --output-dir core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500
   
  
     qiime feature-table summarize \
     --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree.qza \
     --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree.qzv
  


while read line ; do
    echo $line
    path=core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500
    
    qiime diversity alpha-group-significance \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file G2F_metadata_2019_duplicate_pedigree.txt \
  --o-visualization $path/${line}-group-sginificance.qzv
  
       qiime diversity alpha-correlation \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file G2F_metadata_2019_duplicate_pedigree.txt \
  --p-intersect-ids TRUE \
  --o-visualization $path/${line}_correlation.qzv
    
    
done <alpha_diversity_metrics.txt



while read line ; do
    echo $line
    path=core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500
    
    qiime diversity beta-group-significance \
  --i-distance-matrix $path/${line}.qza \
  --m-metadata-file G2F_metadata_2019_duplicate_pedigree.txt \
  --o-visualization $path/${line}_corrected_pedigree.qzv \
  --m-metadata-column Corrected_pedigree \
  --p-pairwise
 
   qiime diversity beta-group-significance \
  --i-distance-matrix $path/${line}.qza \
  --m-metadata-file G2F_metadata_2019_duplicate_pedigree.txt \
  --o-visualization $path/${line}_location.qzv \
  --m-metadata-column location \
  --p-pairwise
  
  
done <beta_diversity_metrics.txt

  
  
##############################################################################3


 qiime feature-table filter-samples\
  --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree.qza \
  --m-metadata-file G2f_2019_sub_corrected_metadata_no_blank_with_sum_weather.txt \
  --p-where "[location] IN ('NCH1')" \
  --p-exclude-ids TRUE \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1.qza
  
  
  qiime feature-table group \
  --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1.qza \
  --p-axis 'sample' \
  --m-metadata-file G2f_2019_sub_corrected_metadata_no_blank_with_sum_weather.txt \
  --m-metadata-column location \
  --p-mode 'sum' \
  --o-grouped-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location.qza
  

 qiime feature-table summarize \
 --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location.qza \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location.qzv 

     qiime diversity core-metrics-phylogenetic \
   --i-phylogeny dada2-rooted-tree.qza \
   --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location.qza  \
   --p-sampling-depth 44000 \
   --m-metadata-file G2F_2019_environmental_data_revised_for_mantel.csv \
   --output-dir core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-44000

############################################################
   
results_path=core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-44000
    

correlation_method=("pearson")

for file in $results_path; \
do \
  echo $file
  
  
while read line ; do
    echo $line
    
       qiime metadata distance-matrix \
  --m-metadata-file G2F_2019_environmental_data_revised_for_mantel.csv \
  --m-metadata-column $line \
  --o-distance-matrix $file/${line}_distance_matrix.qza 

  
for method in $correlation_method; \
  do 
      echo $method

  qiime diversity mantel \
  --i-dm1 $file/${line}_distance_matrix.qza \
  --i-dm2 $file/weighted_unifrac_distance_matrix.qza \
  --p-method $method \
  --p-permutations 999 \
  --p-label1 ${line}_concentration \
  --p-label2 weighted_unifrac_distance \
  --p-intersect-ids TRUE \
  --o-visualization $file/${line}_weighted_unifrac_pearson_association.qzv
  
    qiime diversity mantel \
  --i-dm1 $file/${line}_distance_matrix.qza  \
  --i-dm2 $file/unweighted_unifrac_distance_matrix.qza \
  --p-method $method \
  --p-permutations 999 \
  --p-label1 $line \
  --p-label2 unweighted_unifrac_distance \
  --p-intersect-ids TRUE \
  --o-visualization $file/${line}_unweighted_unifrac_pearson_association.qzv

  
done
  
done <column_name_list.txt

done


############################################################

     qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/1500_depth_rarefied_table_species.qza 
  

     qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/1500_depth_rarefied_table_genus.qza 

       qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/1500_depth_rarefied_table_family.qza 
  
       qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/1500_depth_rarefied_table_class.qza 
  
       qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/1500_depth_rarefied_table_order.qza 
  
       qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/1500_depth_rarefied_table_phylumn.qza 
  

 qiime tools export \
    --input-path  dada2_table-no-mitochondria-no-chloroplast_genus.qza \
    --output-path ./exported_table/
     
    biom convert -i ./exported_table/feature-table.biom \
       -o ./exported_table/dada2_table-no-mitochondria-no-chloroplast_genus.txt --to-tsv
       

           qiime tools export \
    --input-path  dada2_table-no-mitochondria-no-chloroplast_family.qza \
    --output-path ./exported_table/
     
    biom convert -i ./exported_table/feature-table.biom \
       -o ./exported_table/dada2_table-no-mitochondria-no-chloroplast_family.txt --to-tsv
       
           qiime tools export \
    --input-path  dada2_table-no-mitochondria-no-chloroplast_order.qza  \
    --output-path ./exported_table/
     
    biom convert -i ./exported_table/feature-table.biom \
       -o ./exported_table/dada2_table-no-mitochondria-no-chloroplast_order.txt --to-tsv
       
           qiime tools export \
    --input-path  dada2_table-no-mitochondria-no-chloroplast_class.qza \
    --output-path ./exported_table/
     
    biom convert -i ./exported_table/feature-table.biom \
       -o ./exported_table/dada2_table-no-mitochondria-no-chloroplast_class.txt --to-tsv
       
           qiime tools export \
    --input-path  dada2_table-no-mitochondria-no-chloroplast_phylumn.qza \
    --output-path ./exported_table/
     
    biom convert -i ./exported_table/feature-table.biom \
       -o ./exported_table/dada2_table-no-mitochondria-no-chloroplast_phylumn.txt --to-tsv
       
       
       qiime feature-table core-features 
       

qiime feature-table core-features \
 --i-table 1500_depth_rarefied_table_phylumn.qza \
 --p-min-fraction 0.2 \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-core-microbiome-phylumn-1500.qzv

qiime feature-table core-features \
 --i-table 1500_depth_rarefied_table_class.qza \
 --p-min-fraction 0.2 \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-core-microbiome-class-1500.qzv
 
 qiime feature-table core-features \
 --i-table 1500_depth_rarefied_table_order.qza \
 --p-min-fraction 0.2 \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-core-microbiome-order-1500.qzv
 
 qiime feature-table core-features \
 --i-table 1500_depth_rarefied_table_family.qza \
 --p-min-fraction 0.2 \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-core-microbiome-family-1500.qzv
 
 qiime feature-table core-features \
 --i-table 1500_depth_rarefied_table_genus.qza \
 --p-min-fraction 0.2 \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-core-microbiome-genus-1500.qzv


 qiime feature-table core-features \
 --i-table 1500_depth_rarefied_table_species.qza \
 --p-min-fraction 0.2 \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-core-microbiome-species-1500.qzv
 
 
#############################################################


 qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_species.qza \
  --p-min-samples 70 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-species_20_filtered.qza 
  
    
   qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_genus.qza \
  --p-min-samples 70 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-genus_20_filtered.qza 
  
          qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_family.qza \
  --p-min-samples 70 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-family_20_filtered.qza
  
        qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_order.qza \
  --p-min-samples 70 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-order_20_filtered.qza 
  
     qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_class.qza \
  --p-min-samples 70 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-class_20_filtered.qza 
  
  
          qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_phylumn.qza \
  --p-min-samples 70 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-phylumn_20_filtered.qza 

##########################################################################3


   qiime taxa barplot \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --m-metadata-file G2F_metadata_2019_duplicate_pedigree.txt \
  --o-visualization dada2_table-no-mitochondria-no-chloroplast-duplicate-pedigree-pretrained-taxa-bar-plots.qzv 
  
  
  
 qiime feature-table filter-samples\
  --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree.qza \
  --m-metadata-file G2f_2019_sub_corrected_metadata_no_blank_with_sum_weather.txt \
  --p-where "[location] IN ('GAH2')" \
  --p-exclude-ids FALSE \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-GAH2.qza


 qiime feature-table summarize \
 --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-GAH2.qza \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-GAH2.qzv 

    qiime diversity core-metrics-phylogenetic \
   --i-phylogeny dada2-rooted-tree.qza \
   --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-GAH2.qza  \
   --p-sampling-depth 1500 \
   --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
   --output-dir core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-GAH2-1500
   
   
#################################################################################



 qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_species.qza \
  --p-min-samples 208 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-species_60_filtered.qza 
  
    
   qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_genus.qza \
  --p-min-samples 208 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-genus_60_filtered.qza 
  
          qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_family.qza \
  --p-min-samples 208 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-family_60_filtered.qza
  
        qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_order.qza \
  --p-min-samples 208 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-order_60_filtered.qza 
  
     qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_class.qza \
  --p-min-samples 208 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-class_60_filtered.qza 
  
  
          qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_phylumn.qza \
  --p-min-samples 208 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-phylumn_60_filtered.qza 

  
###############################################

location_string="GAH2 IAH2 IAH4 INH1 MNH1 MOH1 NCH1 NEH1 NEH2 OHH1 SCH1"

for experiment_field in $location_string
do
 echo $experiment_field
 
  qiime feature-table filter-samples\
  --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree.qza \
  --m-metadata-file G2f_2019_sub_corrected_metadata_no_blank_with_sum_weather.txt \
  --p-where "[location] IN ('${experiment_field}')" \
  --p-exclude-ids FALSE \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${experiment_field}.qza

  qiime feature-table summarize \
 --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${experiment_field}.qza \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${experiment_field}.qzv 

    qiime diversity core-metrics-phylogenetic \
   --i-phylogeny dada2-rooted-tree.qza \
   --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${experiment_field}.qza  \
   --p-sampling-depth 1500 \
   --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
   --output-dir core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${experiment_field}-1500
   
done


location_string="GAH2 IAH2 IAH4 INH1 MNH1 MOH1 NCH1 NEH1 NEH2 OHH1 SCH1"

for experiment_field in $location_string
do
 echo $experiment_field
 
 directory_experiemnt_field=core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${experiment_field}-1500
 
 echo $directory_experiemnt_field
 
 fraction_list="0.6 0.7"
 taxanomy_levels="2 3 4 5 6 7"
 
for fraction_ratio in $fraction_list
do

for level in $taxanomy_levels
do
    

      qiime taxa collapse \
  --i-table ${directory_experiemnt_field}/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level $level \
  --o-collapsed-table ${directory_experiemnt_field}/rarefied_table_${level}_level_taxonomy.qza 
  
      qiime feature-table filter-features-conditionally \
  --i-table ${directory_experiemnt_field}/rarefied_table_${level}_level_taxonomy.qza  \
  --p-abundance 0.0000000001 \
  --p-prevalence $fraction_ratio \
  --o-filtered-table ${directory_experiemnt_field}/rarefied_table_${level}_level_taxonomy_${fraction_ratio}_filtered.qza
  
   qiime feature-table summarize \
 --i-table ${directory_experiemnt_field}/rarefied_table_${level}_level_taxonomy_${fraction_ratio}_filtered.qza \
 --o-visualization ${directory_experiemnt_field}/rarefied_table_${level}_level_taxonomy_${fraction_ratio}_filtered.qzv


done

done

done

 
 
 
 
 for experiment_field in $location_string
do
 echo $experiment_field
 
 directory_experiemnt_field=core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${experiment_field}-1500
 
 echo $directory_experiemnt_field
 
 fraction_list="0.2 0.3 0.4 0.5"
 taxanomy_levels="2 3 5 6 7"
 

for level in $taxanomy_levels
do
         qiime taxa collapse \
  --i-table ${directory_experiemnt_field}/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level $level \
  --o-collapsed-table ${directory_experiemnt_field}/rarefied_table_${level}_level_taxonomy.qza  

   qiime feature-table core-features \
 --i-table ${directory_experiemnt_field}/rarefied_table_${level}_level_taxonomy.qza \
 --p-min-fraction 0.2 \
 --o-visualization ${directory_experiemnt_field}/rarefied_table_${level}_level_taxonomy_core_microbiome.qzv
  


done

done
 

    qiime taxa collapse \
  --i-table ${directory_experiemnt_field}/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/1500_depth_rarefied_table_species.qza 



############################################################

Correct_pedigree_string="2369XLH123HT B14AXH95 B14AXMO17 B14AXOH43 B37XH95 B37XMO17 B37XOH43 B73XMO17 B73XPHM49 B73XPHN82 CG119XCG108 CG44XCGR01 F42XH95 F42XMO17 F42XOH43 LH74XPHN82 OH43XB37 PHG39XPHN82 PHW52XPHM49 PHW52XPHN82"

for Correct_pedigree in $Correct_pedigree_string
do
 echo $Correct_pedigree
 
  qiime feature-table filter-samples\
  --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree.qza \
  --m-metadata-file G2f_2019_sub_corrected_metadata_no_blank_change_X.csv \
  --p-where "[Corrected_pedigree] IN ('${Correct_pedigree}')" \
  --p-exclude-ids FALSE \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}.qza

  qiime feature-table summarize \
 --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}.qza \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}.qzv 

    qiime diversity core-metrics-phylogenetic \
   --i-phylogeny dada2-rooted-tree.qza \
   --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}.qza  \
   --p-sampling-depth 1500 \
   --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
   --output-dir core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}-1500
   
done

Correct_pedigree_string="2369XLH123HT B14AXH95 B14AXMO17 B14AXOH43 B37XH95 B37XMO17 B37XOH43 B73XMO17 B73XPHM49 B73XPHN82 CG119XCG108 CG44XCGR01 F42XH95 F42XMO17 F42XOH43 LH74XPHN82 OH43XB37 PHG39XPHN82 PHW52XPHM49 PHW52XPHN82"

for Correct_pedigree in $Correct_pedigree_string
do
 echo $Correct_pedigree
 
 directory_Correct_pedigree=core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}-1500
 
 echo $directory_Correct_pedigree
 
 fraction_list="0.5 0.6 0.7 0.8"
 taxanomy_levels="2 3 4 5 6 7"
 
for fraction_ratio in $fraction_list
do

for level in $taxanomy_levels
do

      qiime taxa collapse \
  --i-table ${directory_Correct_pedigree}/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level $level \
  --o-collapsed-table ${directory_Correct_pedigree}/rarefied_table_${level}_level_taxonomy.qza 
  
      qiime feature-table filter-features-conditionally \
  --i-table ${directory_Correct_pedigree}/rarefied_table_${level}_level_taxonomy.qza  \
  --p-abundance 0.0000000001 \
  --p-prevalence $fraction_ratio \
  --o-filtered-table ${directory_Correct_pedigree}/rarefied_table_${level}_level_taxonomy_${fraction_ratio}_filtered.qza
  
   qiime feature-table summarize \
 --i-table ${directory_Correct_pedigree}/rarefied_table_${level}_level_taxonomy_${fraction_ratio}_filtered.qza \
 --o-visualization ${directory_Correct_pedigree}/rarefied_table_${level}_level_taxonomy_${fraction_ratio}_filtered.qzv


done

done

done


############################################################

 qiime feature-table filter-samples \
  --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree.qza \
  --m-metadata-file G2f_2019_sub_corrected_metadata_no_blank_with_sum_weather.txt \
  --p-where "[location] IN ('IAH2')" \
  --p-exclude-ids FALSE \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-IAH2.qza
  

qiime deicode rpca \
    --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-rarefraction-ordination.qza \
    --o-distance-matrix g2f_aitchison_distance.qza

    
     qiime diversity beta-group-significance \
    --i-distance-matrix g2f_aitchison_distance.qza \
    --m-metadata-file G2f_2019_sub_corrected_metadata_no_blank_with_sum_weather.txt \
    --p-method permanova \
    --p-pairwise TRUE \
    --m-metadata-column location \
    --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-aitchison-location-significance.qzv

    
         qiime diversity beta-group-significance \
    --i-distance-matrix g2f_aitchison_distance.qza \
    --m-metadata-file G2f_2019_sub_corrected_metadata_no_blank_with_sum_weather.txt \
    --p-method permanova \
    --p-pairwise TRUE \
    --m-metadata-column pedigree \
    --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-aitchison-pedigree-significance.qzv

    
     qiime emperor biplot \
    --i-biplot dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-rarefraction-ordination.qza \
    --m-sample-metadata-file G2f_2019_sub_corrected_metadata_no_blank_with_sum_weather_tab.csv \
    --m-feature-metadata-file silva-138-99_nb_DADA2_taxonomy.qza \
    --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-rarefraction-ordination.qzv \
    --p-number-of-features 8

    
    qiime deicode rpca \
    --i-table August-table-no-mitochondria-no-chloroplast.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot August_aitchison-no-rarefraction-ordination.qza \
    --o-distance-matrix August_aitchison_distance.qza

     qiime diversity beta-group-significance \
    --i-distance-matrix August_aitchison_distance.qza \
    --m-metadata-file living_mulch_data_metafile_massalin_no_blank.tsv \
    --p-method permanova \
    --p-pairwise TRUE \
    --m-metadata-column treatment \
    --o-visualization August_aitchison_treatment_significance.qzv

    
         qiime emperor biplot \
    --i-biplot August_aitchison-no-rarefraction-ordination.qza \
    --m-sample-metadata-file living_mulch_data_metafile_massalin_no_blank.tsv \
    --m-feature-metadata-file living_mulch_no_blank_filtered_taxonomy.qza \
    --o-visualization August_aitchison_distance-biplot.qzv \
    --p-number-of-features 8
        
    
    
    qiime deicode rpca \
    --i-table June-table-no-mitochondria-no-chloroplast.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot June_aitchison-no-rarefraction-ordination.qza \
    --o-distance-matrix June_aitchison_distance.qza

     qiime diversity beta-group-significance \
    --i-distance-matrix June_aitchison_distance.qza \
    --m-metadata-file living_mulch_data_metafile_massalin_no_blank.tsv \
    --p-method permanova \
    --p-pairwise TRUE \
    --m-metadata-column treatment \
    --o-visualization June_aitchison_treatment_significance.qzv

    
         qiime emperor biplot \
    --i-biplot June_aitchison-no-rarefraction-ordination.qza \
    --m-sample-metadata-file living_mulch_data_metafile_massalin_no_blank.tsv \
    --m-feature-metadata-file living_mulch_no_blank_filtered_taxonomy.qza \
    --o-visualization June_aitchison_distance-biplot.qzv \
    --p-number-of-features 8
       
    
    
    
    qiime deicode rpca \
    --i-table May-table-no-mitochondria-no-chloroplast.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot May_aitchison-no-rarefraction-ordination.qza \
    --o-distance-matrix May_aitchison_distance.qza

     qiime diversity beta-group-significance \
    --i-distance-matrix May_aitchison_distance.qza \
    --m-metadata-file living_mulch_data_metafile_massalin_no_blank.tsv \
    --p-method permanova \
    --p-pairwise TRUE \
    --m-metadata-column treatment \
    --o-visualization May_aitchison_treatment_significance.qzv

    
         qiime emperor biplot \
    --i-biplot May_aitchison-no-rarefraction-ordination.qza \
    --m-sample-metadata-file living_mulch_data_metafile_massalin_no_blank.tsv \
    --m-feature-metadata-file living_mulch_no_blank_filtered_taxonomy.qza \
    --o-visualization May_aitchison_distance-biplot.qzv \
    --p-number-of-features 8

    
    
results_path=core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000
    

correlation_method=("spearman")

for file in $results_path; \
do \
  echo $file
  
while read line ; do
    echo $line

for method in $correlation_method; \
  do 
      echo $method

  qiime diversity mantel \
  --i-dm1 $file/August_{$line}_distance_matrix.qza \
  --i-dm2 $file/August_aitchison_distance.qza \
  --p-method $method \
  --p-permutations 999 \
  --p-label1 ${line}_concentration \
  --p-label2 aitchison_distance \
  --p-intersect-ids TRUE \
  --o-visualization $file/${line}_aitchison_pearson_association.qzv
 
  
done

done <column_name_list.txt
  
done

############################################################################################


  qiime diversity core-metrics-phylogenetic \
   --i-phylogeny dada2-rooted-tree.qza \
   --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1.qza  \
   --p-sampling-depth 1500 \
   --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
   --output-dir core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-1500
   
  
     qiime feature-table summarize \
     --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1.qza \
     --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1.qzv
  


while read line ; do
    echo $line
    path=core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-1500
    
    qiime diversity alpha-group-significance \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file G2f_2019_sub_corrected_metadata_for_bioenv_numeric.tsv \
  --o-visualization $path/${line}-group-sginificance.qzv
  
       qiime diversity alpha-correlation \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file G2f_2019_sub_corrected_metadata_for_bioenv_numeric.tsv \
  --o-visualization $path/${line}_correlation.qzv
    
    
done <alpha_diversity_metrics.txt



while read line ; do
    echo $line
    path=core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-1500
    
    qiime diversity beta-group-significance \
  --i-distance-matrix $path/${line}.qza \
  --m-metadata-file G2F_metadata_2019_duplicate_pedigree.txt \
  --o-visualization $path/${line}_corrected_pedigree.qzv \
  --m-metadata-column Corrected_pedigree \
  --p-pairwise
 
   qiime diversity beta-group-significance \
  --i-distance-matrix $path/${line}.qza \
  --m-metadata-file G2F_metadata_2019_duplicate_pedigree.txt \
  --o-visualization $path/${line}_location.qzv \
  --m-metadata-column location \
  --p-pairwise
  
  
done <beta_diversity_metrics.txt


qiime diversity bioenv \
  --i-distance-matrix core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-1500/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file G2f_2019_sub_corrected_metadata_for_bioenv_numeric.tsv \
  --o-visualization g2f_bioenviroment_duplicate_dataset.qzv 

                                                             
#dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1.qza

###################################################################################


while read line ; do
    echo $line
    path=core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-1500
    
    qiime diversity alpha-group-significance \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file G2F_metadata_2019_duplicate_pedigree.txt \
  --o-visualization $path/${line}-group-sginificance.qzv
  
       qiime diversity alpha-correlation \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file G2F_metadata_2019_duplicate_pedigree.txt \
  --o-visualization $path/${line}_correlation.qzv
    
    
done <alpha_diversity_metrics.txt

########################################################################################


qiime feature-table filter-samples\
  --i-table rarefied_table.qza \
  --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
  --p-where "[location] IN ('NCH1')" \
  --p-exclude-ids TRUE \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-1500.qza
  
  
  qiime feature-table group \
  --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-1500.qza \
  --p-axis 'sample' \
  --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
  --m-metadata-column location \
  --p-mode 'sum' \
  --o-grouped-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table.qza
  
  qiime feature-table summarize \
 --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table.qza \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table.qzv

     qiime diversity core-metrics-phylogenetic \
   --i-phylogeny dada2-rooted-tree.qza \
   --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table.qza  \
   --p-sampling-depth 18000 \
   --m-metadata-file G2F_2019_environmental_data_revised_for_mantel.csv \
   --output-dir core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000

##############################################################################


qiime feature-table group \
  --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-1500.qza \
  --p-axis 'sample' \
  --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
  --m-metadata-column location \
  --p-mode 'sum' \
  --o-grouped-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table.qza
  
  

  
Correct_pedigree_string="2369XLH123HT B14AXH95 B14AXMO17 B14AXOH43 B37XH95 B37XMO17 B37XOH43 B73XMO17 B73XPHM49 B73XPHN82 CG119XCG108 CG44XCGR01 F42XH95 F42XMO17 F42XOH43 LH74XPHN82 OH43XB37 PHG39XPHN82 PHW52XPHM49 PHW52XPHN82"

for Correct_pedigree in $Correct_pedigree_string
do
 echo $Correct_pedigree
 
 directory_Correct_pedigree=core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}-1500
 
 echo $directory_Correct_pedigree
 
 qiime feature-table group \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}-1500/rarefied_table.qza \
  --p-axis 'sample' \
  --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
  --m-metadata-column location \
  --p-mode 'sum' \
  --o-grouped-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}-by-location.qza
  
   qiime feature-table summarize \
 --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}-by-location.qza \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}-by-location.qzv 

    qiime diversity core-metrics-phylogenetic \
   --i-phylogeny dada2-rooted-tree.qza \
   --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}-by-location.qza  \
   --p-sampling-depth 3000 \
   --m-metadata-file G2F_2019_environmental_data_revised_for_mantel.csv \
   --output-dir core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-${Correct_pedigree}-by-location-3000
  

done

#######################################################################################

qiime feature-table group \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-1500/rarefied_table.qza \
  --p-axis 'sample' \
  --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
  --m-metadata-column location \
  --p-mode 'sum' \
  --o-grouped-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-by-location.qza
  
  
  qiime feature-table summarize \
 --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-by-location.qza \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-by-location.qzv
 
  
    qiime diversity core-metrics-phylogenetic \
   --i-phylogeny dada2-rooted-tree.qza \
   --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-by-location.qza  \
   --p-sampling-depth 18000 \
   --m-metadata-file G2F_2019_environmental_data_revised_for_mantel.csv \
   --output-dir core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-by-location-18000

##########################################################################################3

 qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_species.qza \
  --p-min-samples 33 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-species_10_filtered.qza 
  
    
   qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_genus.qza \
  --p-min-samples 33 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-genus_10_filtered.qza 
  
          qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_family.qza \
  --p-min-samples 33 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-family_10_filtered.qza
  
        qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_order.qza \
  --p-min-samples 33 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-order_10_filtered.qza 
  
     qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_class.qza \
  --p-min-samples 33 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-class_10_filtered.qza 
  
  
    qiime feature-table filter-features \
  --i-table 1500_depth_rarefied_table_phylumn.qza \
  --p-min-samples 33 \
  --o-filtered-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-duplicate-phylumn_10_filtered.qza 

  
  ##########################################################################################
  
  
    qiime feature-table group \
  --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-1500.qza \
  --p-axis 'sample' \
  --m-metadata-file G2f_2019_sub_corrected_metadata.tsv \
  --m-metadata-column location \
  --p-mode 'sum' \
  --o-grouped-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table.qza
  
  qiime feature-table summarize \
 --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table.qza \
 --o-visualization dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table.qzv

     qiime diversity core-metrics-phylogenetic \
   --i-phylogeny dada2-rooted-tree.qza \
   --i-table dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table.qza  \
   --p-sampling-depth 18000 \
   --m-metadata-file G2F_2019_environmental_data_revised_for_mantel.csv \
   --output-dir core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000
   
   
     qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/18000_location_depth_rarefied_table_species.qza 
  
   
    qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/18000_location_depth_rarefied_table_genus.qza 
  
   qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/18000_location_depth_rarefied_table_family.qza 
  
   qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/18000_location_depth_rarefied_table_class.qza 
  
   qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/18000_location_depth_rarefied_table_order.qza 
  
   qiime taxa collapse \
  --i-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/rarefied_table.qza \
  --i-taxonomy silva-138-99_nb_DADA2_taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/18000_location_depth_rarefied_table_phylumn.qza 
  

   
