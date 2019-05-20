## Outliers analysis  
##### Dependencies:
python3
pandas  
numpy  
scipy.stats  
argparse  

Example code in make_outliers_table.sh:

```
#!/bin/bash
updown="up"
location_of_py_file="make_outliers_table.py"
location_of_data_file="test_input_for_outliers.tsv"
iqrs_over_median=0.5 #Note 1.5 IQRs is suggested, this is just for test data.
gene_column_name="Genes"
output_prefix="test_outliers_output_${updown}"
sample_names_file="test_samples.txt"
aggregate=True
write_frac_table="False"

python3.5 ${location_of_py_file} \
--input_df ${location_of_data_file} \
--iqrs_over_median ${iqrs_over_median} \
--gene_column_name ${gene_column_name} \
--output_prefix ${output_prefix} \
--sample_names_file ${sample_names_file} \
--aggregate ${aggregate} \
--up_or_down ${updown} \
--write_frac_table ${write_frac_table}



```

##### Argument explanation:
*location_of_py_file:* Path to make_outliers_table.py file.

*location_of_data_file:* Path to input file. Sample names must match column names. Needs gene names column as well. Can also include other annotation columns, but those will not be propagated to output.  

*iqrs_over_median:* Threshold for a value to be considered an "outlier". Number of inter-quartile ranges over the median, calculated for each gene, to be used as a threshold.  1.5 is suggested.

*gene_column_name:* The name of the column to be used for aggregating values per gene. For instance, if phospho-peptide data is used as input, it may be helpful to aggregate outliers at the gene or isoform level.  

*output_file:* Name/location to put outliers output.   

*sample_names_file:* Column/sample names to be used for analysis. Columns not included in this list will be ignored.   

*aggregate:* Whether to combine rows based on an annotation column. For instance, add up outliers and non-outliers for multiple phospho sites based on a gene column. Choices are True or False, default True.  

*up_or_down*: Whether analysis should mark outliers that are **greater than** median + IQR threshold, or **less than** median - IQR threshold.  

*write_frac_table:* Optional output. Writes table where values are the fraction of sites per sample per gene that were outliers. Not great for stats, but good for visualization. Choices are True or False, default False.    


##### Output explanation:
The output will have a column corresponding to the gene_column_name from the input. For each sample listed in sample_names_file, there will be 2 output columns. For sample X, one column will be called X_outliers which counts the number of outliers for that gene/row for that sample. The second columns will be called X_notOutliers which counts the number of sites per gene/row that are not outliers in that sample. For input data where each each value in gene_column_name is unique, and doesn't have missing values, the values from X_outliers and X_notOutliers will always add up to 1. In other words, each sample will be either an outlier or not an outlier (or be a missing value). This format is helpful for calculating comparisons in downstream analysis.  

Fraction table writes a tab separated file with the fraction of sites per sample per gene that were outliers. This file is good for visualizing custom groups of genes.


## Group comparisons using outliers

##### Dependencies:
python3
pandas  
numpy  
matplotlib  
seaborn  
scipy.stats  
sys  
argparse  
datetime  

Example code in make_comparisons_on_outliers.sh:

```
#!/bin/bash

location_of_py_file="compare_groups_outliers.py"
location_of_outliers_file="test_outliers_output_up.tsv"
gene_column_name="Genes"
fdr_cut_off=0.05
genes_to_highlight="genes_of_interest.txt"
blue_or_red="red"
output_qvals="True"
frac_filter="0.3"

output_prefix="test_group1_comparison"
group1_label="group_of_interest"
group1_list="test_group1.txt"
group2_label="not_in_group_of_interest"
group2_list="test_group2.txt"
group_colors="test_group_colors.txt"

python3.5  ${location_of_py_file} \
--outliers_table  ${location_of_outliers_file} \
--gene_column_name ${gene_column_name} \
--fdr_cut_off ${fdr_cut_off} \
--output_prefix ${output_prefix} \
--group1_label ${group1_label} \
--group1_list ${group1_list} \
--group2_label ${group2_label} \
--group2_list ${group2_list} \
--genes_to_highlight ${genes_to_highlight} \
--blue_or_red ${blue_or_red} \
--group_colors ${group_colors} \
--output_qvals ${output_qvals} \
--frac_filter ${frac_filter}

```

##### Argument explanation:
*location_of_py_file:* Path to compare_groups_outliers.py file.  

*location_of_outliers_file:* Path to outliers output from make_outliers_table.py.  

*gene_column_name:* Name of column used for labeling rows, usually gene names.

*fdr_cut_off:* Threshold for signficance, used for selecting genes that are significantly different between two groups, after multiple testing correction.  

*genes_to_highlight:* Optional. Can input a list of genes that will be marked with "\*" in the output heatmap.  

*blue_or_red:* Color map colors for output heatmap. Useful for differentiating between up and down outliers. Choices are red or blue, default is red.

*output_prefix:* Prefix used for output files.

*group1_label:* Label used for group of interest. If inputting colors key for groups, must match color key label.  

*group1_list:* Path to list of samples belonging to primary group of interest. Samples not included in this list or group2_list will be ignored in analysis.

*group2_label:* Label used for samples **not** in group of interest. If inputting colors key for groups, must match color key label.  

*group2_list:* Path to list of samples **not** in group of interest. Samples not included in this list or group1_list will be ignored in analysis.  

*group_colors:* Tab separated list of group labels with paired hex color for label. If it's not supplies, default colors will be used.  

*output_qvals:* Choose whether to output adjusted p-values for all genes into a table. Choices are True and False, default False.   

*frac_filter:* None or a fraction between 0 and 1. Filter genes for *frac_filter* of samples in group1 having an outlier. Useful for making sure results are not driven by a small subset of samples with high outlier values. 0.3 default.    

##### Output explanation:
This tool has 2 default outputs, and 1 optional output. The first output is a heatmap, visualizing genes that are significantly enriched in the group of interest, compared to the second group. The second output is a list of genes/labels for rows that were found to be significantly enriched in the group of interest. Optionally, you can also output a table of gene names and adjusted p-values.
