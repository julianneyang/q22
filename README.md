# q22
Scripts and data files accompanying publication

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
