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
