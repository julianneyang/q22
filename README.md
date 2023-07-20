# q22
Scripts and data files accompanying publication

# PICRUST2
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/make_tree_for_picrust2.sh rep-seqs.qza 8 picrust2_default_sepp_ref.qza 
Call bash make_tree_for_picrust2.sh rep-seqs.qza number_of_threads reference_database.qza
Saved Phylogeny[Rooted] to: picrust_tree/tree.qza
Saved Placements to: picrust_tree/placements.qza
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/q2-picrust2.sh Q22_ASV.qza picrust_tree/tree.qza picrust2_Q22
Call bash q2-picrust2.sh table.qza tree.qza outputdirname
Running the below commands:
hsp.py -i 16S -t /tmp/tmpb6wwhqr6/placed_seqs.tre -p 1 -n -o /tmp/tmpb6wwhqr6/picrust2_out/16S_predicted.tsv.gz -m mp -e 0.5

hsp.py -i EC -t /tmp/tmpb6wwhqr6/placed_seqs.tre -p 20 -n -o /tmp/tmpb6wwhqr6/picrust2_out/EC_predicted.tsv.gz -m mp -e 0.5

hsp.py -i KO -t /tmp/tmpb6wwhqr6/placed_seqs.tre -p 20 -n -o /tmp/tmpb6wwhqr6/picrust2_out/KO_predicted.tsv.gz -m mp -e 0.5

metagenome_pipeline.py -i /tmp/tmpb6wwhqr6/intable.biom -f /tmp/tmpb6wwhqr6/picrust2_out/EC_predicted.tsv.gz -o /tmp/tmpb6wwhqr6/picrust2_out/EC_metagenome_out --max_nsti 2 -m /tmp/tmpb6wwhqr6/picrust2_out/16S_predicted.tsv.gz
26 of 2234 ASVs were above the max NSTI cut-off of 2.0 and were removed from the downstream analyses.
26 of 2234 ASVs were above the max NSTI cut-off of 2.0 and were removed from the downstream analyses.

metagenome_pipeline.py -i /tmp/tmpb6wwhqr6/intable.biom -f /tmp/tmpb6wwhqr6/picrust2_out/KO_predicted.tsv.gz -o /tmp/tmpb6wwhqr6/picrust2_out/KO_metagenome_out --max_nsti 2 -m /tmp/tmpb6wwhqr6/picrust2_out/16S_predicted.tsv.gz
26 of 2234 ASVs were above the max NSTI cut-off of 2.0 and were removed from the downstream analyses.
26 of 2234 ASVs were above the max NSTI cut-off of 2.0 and were removed from the downstream analyses.

pathway_pipeline.py -i /tmp/tmpb6wwhqr6/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -o /tmp/tmpb6wwhqr6/picrust2_out/pathways_out -p 20

Saved FeatureTable[Frequency] to: picrust2_output_Q22_ASV.qza/ko_metagenome.qza
Saved FeatureTable[Frequency] to: picrust2_output_Q22_ASV.qza/ec_metagenome.qza
Saved FeatureTable[Frequency] to: picrust2_output_Q22_ASV.qza/pathway_abundance.qza

(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files/picrust2_output_Q22_ASV.qza$ for file in *; do bash ../../../../fast-16s-analysis/shell_scripts/4-qiime-tools-export.sh $file;done
Takes qza input file as input and cranks out tsv and summary.txt file ec_metagenome.qza
Exported ec_metagenome.qza as BIOMV210DirFmt to directory export_ec_metagenome
Takes qza input file as input and cranks out tsv and summary.txt file ko_metagenome.qza
Exported ko_metagenome.qza as BIOMV210DirFmt to directory export_ko_metagenome
Takes qza input file as input and cranks out tsv and summary.txt file pathway_abundance.qza
Exported pathway_abundance.qza as BIOMV210DirFmt to directory export_pathway_abundance
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files/picrust2_output_Q22_ASV.qza$ add_descriptions.py -i export_pathway_abundance/feature-table.tsv -m METACYC -o export_pathway_abundance/annotated_pathway.tsv

# Collapse Taxa
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/collapse-taxa.sh Q22_ASV.qza taxonomy.qza
# Taxa Barplots 
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/import-taxonomy.sh Q22_Taxonomy.tsv 
Imported Q22_Taxonomy.tsv as TSVTaxonomyDirectoryFormat to taxonomy.qza

(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/site_subsets$ for file in *; do bash ../../../fast-16s-analysis/shell_scripts/taxabarplot.sh $file ../starting_files/taxonomy.qza ../starting_files/Q22_Metadata.tsv; done

(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/site_subsets$ bash ../../../fast-16s-analysis/shell_scripts/groupsamples.sh ILE_Q22_ASV.qza ../starting_files/Q22_Metadata.tsv Q22 
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
Saved FeatureTable[Frequency] to: groupby_Q22_ILE_Q22_ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/site_subsets$ bash ../../../fast-16s-analysis/shell_scripts/groupsamples.sh COL_Q22_ASV.qza ../starting_files/Q22_Metadata.tsv Q22 
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
Saved FeatureTable[Frequency] to: groupby_Q22_COL_Q22_ASV.qza

(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/site_subsets$ bash ../../../fast-16s-analysis/shell_scripts/taxabarplot.sh groupby_Q22_COL_Q22_ASV.qza ../starting_files/taxonomy.qza site_metadata.tsv 
Call bash taxabarplot.sh table taxonomy metadata
Saved Visualization to: barplot_groupby_Q22_COL_Q22_ASV.qza_dir.qzv
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/site_subsets$ bash ../../../fast-16s-analysis/shell_scripts/taxabarplot.sh groupby_Q22_ILE_Q22_ASV.qza ../starting_files/taxonomy.qza site_metadata.tsv 
Call bash taxabarplot.sh table taxonomy metadata
Saved Visualization to: barplot_groupby_Q22_ILE_Q22_ASV.qza_dir.qzv

# Alpha Diversity 
qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/1-qiime-tools-import.sh Q22_ASV.tsv 
Enter filepath to tsv containing ASV count data Q22_ASV.tsv
Imported Q22_ASV.biom as BIOMV210Format to Q22_ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ ls
'DADA2 for Julianne.R'   errR.Rdata     Q22_ASV.qza           Q22_ASV.tsv        Q22_Taxonomy.tsv
 errF.Rdata              Q22_ASV.biom   Q22_ASV_summary.txt   Q22_Metadata.tsv   rep-seqs.fna
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ nano Q22_ASV_summary.txt 
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh Q22_ASV.qza Q22_Metadata.tsv Site ILE
Enter filepath to the .qza table Q22_ASV.qza
Enter filepath to the metadata file in tsv format Q22_Metadata.tsv
Enter the column name by which to subset the data Site
Enter the value in the column by which to subset the data ILE
Saved FeatureTable[Frequency] to: ILE_Q22_ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ nano Q22_Metadata.tsv 
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh Q22_ASV.qza Q22_Metadata.tsv Site COL
Enter filepath to the .qza table Q22_ASV.qza
Enter filepath to the metadata file in tsv format Q22_Metadata.tsv
Enter the column name by which to subset the data Site
Enter the value in the column by which to subset the data COL
Saved FeatureTable[Frequency] to: COL_Q22_ASV.qza

qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ mkdir site_subsets
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ mv ILE_Q22_ASV.qza COL_Q22_ASV.qza site_subsets/
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files$ cd site_subsets/
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files/site_subsets$ for file in *; do bash ../../../../fast-16s-analysis/shell_scripts/4-qiime-tools-export.sh $file;done
Takes qza input file as input and cranks out tsv and summary.txt file COL_Q22_ASV.qza
Exported COL_Q22_ASV.qza as BIOMV210DirFmt to directory export_COL_Q22_ASV
Takes qza input file as input and cranks out tsv and summary.txt file ILE_Q22_ASV.qza
Exported ILE_Q22_ASV.qza as BIOMV210DirFmt to directory export_ILE_Q22_ASV
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files/site_subsets$ nano export_COL_Q22_ASV/biom-summary.txt 
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files/site_subsets$ nano export_COL_Q22_ASV/biom-summary.txt 
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files/site_subsets$ bash ../../../../fast-16s-analysis/shell_scripts/rarefy_alpha_diversity.sh COL_Q22_ASV.qza 3371
Call bash rarefy_alpha_diversity.sh table.qza sampling_depth_integer
Saved FeatureTable[Frequency] to: d3371_COL_Q22_ASV.qza
Saved SampleData[AlphaDiversity] to: alpha_COL_Q22_ASV/shannon.qza
Saved SampleData[AlphaDiversity] to: alpha_COL_Q22_ASV/chao1.qza
Saved SampleData[AlphaDiversity] to: alpha_COL_Q22_ASV/otus.qza
Saved SampleData[AlphaDiversity] to: alpha_COL_Q22_ASV/pielou_e.qza
Exported chao1.qza as AlphaDiversityDirectoryFormat to directory chao1_dir
Exported otus.qza as AlphaDiversityDirectoryFormat to directory otus_dir
Exported pielou_e.qza as AlphaDiversityDirectoryFormat to directory pielou_e_dir
Exported shannon.qza as AlphaDiversityDirectoryFormat to directory shannon_dir


(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files/site_subsets$ nano export_ILE_Q22_ASV/biom-summary.txt 
(qiime2-2022.2) julianne@laptop:~/Documents/q22/Q22_Microbiome/starting_files/site_subsets$ bash ../../../../fast-16s-analysis/shell_scripts/rarefy_alpha_diversity.sh ILE_Q22_ASV.qza 15940
Call bash rarefy_alpha_diversity.sh table.qza sampling_depth_integer
Saved FeatureTable[Frequency] to: d15940_ILE_Q22_ASV.qza
Saved SampleData[AlphaDiversity] to: alpha_ILE_Q22_ASV/shannon.qza
Saved SampleData[AlphaDiversity] to: alpha_ILE_Q22_ASV/chao1.qza
Saved SampleData[AlphaDiversity] to: alpha_ILE_Q22_ASV/otus.qza
Saved SampleData[AlphaDiversity] to: alpha_ILE_Q22_ASV/pielou_e.qza
Exported chao1.qza as AlphaDiversityDirectoryFormat to directory chao1_dir
Exported otus.qza as AlphaDiversityDirectoryFormat to directory otus_dir
Exported pielou_e.qza as AlphaDiversityDirectoryFormat to directory pielou_e_dir
Exported shannon.qza as AlphaDiversityDirectoryFormat to directory shannon_dir
