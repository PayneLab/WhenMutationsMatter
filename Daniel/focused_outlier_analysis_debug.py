#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import cptac
import binarization_functions as bf
import blackSheepCPTACmoduleCopy as blsh
import gseapy as gp
from gseapy.plot import barplot, heatmap, dotplot
import json
import requests
import random

cptac.download(dataset='renalccrcc')

re = cptac.RenalCcrcc()
proteomics = re.get_proteomics()
clinical = re.get_clinical()

#Data cleaning. We want to get rid of non_clear cell carcinoma, 
#as it skews our analysis, as it is not what we are studying
non_clear_cell_filter = clinical['histologic_type'] == 'non-Clear Cell renal cell carcinoma'
clinical = clinical[~non_clear_cell_filter]

#We also need to get rid of those samples in the datasets we plan to use. In this case, Proteomics.
rows_to_drop = ['S011', 'S097', 'S104', 'S107', 'S115', 'S119', 'S123']

proteomics.drop(rows_to_drop)

print("Rows Dropped:")
print(rows_to_drop)

columns_to_explore = ['tumor_size_in_cm', 
                      'tumor_focality', 
                      'history_of_cancer', 
                      'vital_status_at_12months_follow_up', 
                      'tumor_status_at_12months_follow_up', 
                      'tumor_necrosis', 
                      'margin_status']

#Create a subsetted copy of the original Clinical DataFrame
annotations = pd.DataFrame(clinical[columns_to_explore].copy())

#Then we will make any non-binary columns binary
col = 'tumor_size_in_cm'
col_mean = annotations[col].mean()
annotations[col]= bf.binarizeCutOff(annotations, col, col_mean, 
                                    "Above_Mean("+str(round(col_mean, 2))+")", 
                                    "Below_Mean("+str(round(col_mean, 2))+")")

col2 = 'tumor_status_at_12months_follow_up'
annotations[col2] = annotations[col2].replace('Unknown', np.nan)

outliers_up = blsh.make_outliers_table(proteomics, iqrs=1.5, 
                                       up_or_down='up', 
                                       aggregate=False, 
                                       frac_table=False)

length = int(len(outliers_up.index) / 2)
only_outliers_up = outliers_up[:length]
outliers_up_dict = {}
for i in range(length):
    key = proteomics.index[i]
    my_filter = only_outliers_up.iloc[i] == 1.0
    value = list(only_outliers_up.iloc[i][my_filter].index)
    outliers_up_dict[key] = value

bf.renameDuplicateColumns(outliers_up)

results_up = blsh.compare_groups_outliers(outliers_up, 
                                          annotations)

results_up = results_up.dropna(axis=0, how='all')

#Drop Columns with less than 4 significant up-regulated enrichments
sig_cols_up = []
for col in results_up.columns:
    sig_col = bf.significantEnrichments(results_up, col)
    if sig_col is not None and len(sig_col) >= 4:
        sig_cols_up.append(sig_col)
    else:
        results_up = results_up.drop(col, axis=1)

#Here we will link clinical attributes with significantly up-regulated genes
sig_genes_up = {}
for i, col in enumerate(sig_cols_up):
    list_of_genes = list(col.index)
    sig_genes_up[sig_cols_up[i].columns[0][:-9]] = list_of_genes   

# To perform a request specifically for inhibitors, you may opt to use a loop
# While this option is slower than the map, it is more specific, and may be worth
# waiting a few more seconds to reduce manual sifting through interaction types
inhibitors = []
for genes in sig_genes_up.values():
    inhibitors.append(bf.dgidb_get_request(genes, interaction_types=['inhibitor']))

print('UP-REGULATED INHIBITOR REQUEST:\n')
inhibitor_dict = {}
for i, request in enumerate(inhibitors):
    clinical_attribute = sig_cols_up[i].columns[0][:-9]
    inhibitor_dict[clinical_attribute] = bf.dgidb_json_parse(request, genes=True)
print(json.dumps(inhibitor_dict, indent=4))

patient_drugs_genes_up = bf.compare_enrichments_with_drugs(outliers_up_dict, clinical)

patients_to_check_up = []
for i in range(3):
    index = random.randrange(0, len(only_outliers_up))
    patients_to_check_up.append(clinical.index[index])

#%%
personalized_up = {}
for patient in patients_to_check_up:
    import pdb; pdb.set_trace()
    patients_up_dict = bf.dgidb_get_request(outliers_up_dict[patient], 
                                            interaction_types = ['inhibitor'])
    parsed = bf.dgidb_json_parse(patients_up_dict, genes=True)
    personalized_up[patient] = parsed
print(json.dumps(personalized_up, indent = 4))
