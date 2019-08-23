import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Create boxplot and stripplot with pval annotation
def cis_plot(df, gene, omics_name, pval, mutation_type="Mutated"):
    omics_col = gene+"_"+omics_name
    
    # get right order for boxplots
    if mutation_type == "Mutated":
        comparison_list = ['Wildtype', 'Mutated']
    elif mutation_type == "Missense":
        comparison_list = ['Wildtype', 'Missense']
    elif mutation_type == "Truncation":
        comparison_list = ['Wildtype', 'Truncation']
        
    # get pval from dataframe or float 
    if isinstance(pval, pd.DataFrame):
        pval_series = pval['P_Value']
        num_pval = float(pval_series[0])
        str_pval = str(pval_series[0])
    elif isinstance(pval, float):
        num_pval = pval
        str_pval = str(pval)
        
    # Boxplot and Stripplot
    plt.rcParams['figure.figsize']=(8,5)
    sns.set(font_scale = 1.3)
    cis_boxplot = sns.boxplot(data = df, x = 'binary_mutations',
                              y = omics_col, order = comparison_list, showfliers = False)  
    cis_boxplot.set_title(
        gene + " Effect on " + gene +" "+omics_name.capitalize()+" in Kidney Tumors\n P-Value = "+str_pval[:6]+"\n")
    cis_boxplot = sns.stripplot(data= df, x = 'binary_mutations',
                                y = omics_col,jitter = True, color = ".3", order = comparison_list)
    cis_boxplot.set(xlabel = "\n"+gene + " Mutation Status in Tumors", ylabel = omics_name.capitalize())
    cis_boxplot.set_xticklabels(cis_boxplot.get_xticklabels())
    
    # pval annotation
    bonferroni_cutoff = .05/6
    if num_pval <= bonferroni_cutoff:
        pval_symbol = "*"
    else:
        pval_symbol = "ns"
    
    x1, x2 = 0, 1   # columns (first column: 0, see plt.xticks())
    y, h = df[omics_col].max() + .05, .05  
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, color= '.3')
    plt.text((x1+x2)*.5,y+h, pval_symbol, horizontalalignment='center', verticalalignment='bottom', color = "black")

    plt.show()
    plt.clf()
    plt.close()


#statistical annotation
def pval_annotation(df, pval_symbol_text, col_A=0, col_B=1, below=False):
    
    x1, x2 = col_A, col_B   # columns (first column: 0, see plt.xticks())
    if below == True:
        y, h = df[col_A].max() + .05, .05  
    else:
        y, h = df[col_A].max() + .2, .05 
    plt.plot([x1, x1, x2, x2], #draw horizontal line
             [y, y+h, y+h, y], #vertical line
             lw=1.5, color= '.3')
    plt.text((x1+x2)*.5, # half between x coord
             y+h, pval_symbol_text, horizontalalignment='center', verticalalignment='bottom', color = ".3")


def add_significance_col(results_df, num_comparisons):
    "bonferroni multiple hypothesis"""
    alpha = .05
    bonferroni_cutoff = alpha / num_comparisons
    
    pval = results_df['P_Value']
    if pval[0] <= bonferroni_cutoff:
        results_df['Significant'] = True
    else: 
        results_df['Significant'] = False
    return results_df

def format_cis_comparison_data(cancer_object, omics_name, gene):
    import numpy as np
    # Step 1 - Create dataframe in order to do comparisons with wrap_ttest - drop nan values
    omics_and_mutations = cancer_object.join_omics_to_mutations(
        mutations_genes = gene, omics_df_name = omics_name, omics_genes = gene).dropna()
    #print(omics_and_mutations.head())
    # Check if values in omics data (if not found in proteomics, after na dropped dataframe should be empty)
    if omics_and_mutations[gene+"_"+omics_name].empty:
        print('Not possible to do T-test. No data for', gene, 'in', omics_name)
        return None
    else:
        
        # Step 2 - Create the binary column needed to do the comparison
        omics_and_mutations['binary_mutations'] = np.where(
            omics_and_mutations[gene+'_Mutation_Status'] == 'Wildtype_Tumor', 'Wildtype', 'Mutated')

        # Step 3 - Format the dataframe correctly for the T-test(just omics and binary columns for tumors)
        tumors = omics_and_mutations.loc[omics_and_mutations['Sample_Status'] == 'Tumor'] #drop Normal samples
        columns_to_drop = [gene+"_Mutation", gene+"_Location", gene+"_Mutation_Status", "Sample_Status"]
        omics_binary_mutations = tumors.drop(columns_to_drop, axis = 1)
        #check if only one column of omics data (total 2 columns)
        if len(omics_binary_mutations.columns) != 2:
            print('exeption with columns. check omics data')
            return None
            
        else:
            return omics_binary_mutations

def get_missense_truncation_comparison(cancer_object, omics_name, gene):
    import numpy as np
    #get omics data and tumors
    omics_and_mutations = cancer_object.join_omics_to_mutations(
                mutations_genes = gene, omics_df_name = omics_name, omics_genes = gene)
    tumors = omics_and_mutations.loc[omics_and_mutations['Sample_Status'] == 'Tumor'] #drop Normal samples


    somatic_mutations = cancer_object.get_mutations().reset_index()

    if cancer_object.get_cancer_type() == 'colon':
        missence_truncation_groups = {'frameshift substitution': 'Truncation', 
            'frameshift deletion': 'Truncation', 'frameshift insertion': 'Truncation', 
            'stopgain': 'Truncation', 'stoploss': 'Truncation', 'nonsynonymous SNV': 'Missense',
            'nonframeshift insertion': 'Missense','nonframeshift deletion': 'Missense', 
            'nonframeshift substitution': 'Missense'}
    else: 
        missence_truncation_groups = {'In_Frame_Del': 'Missense', 'In_Frame_Ins': 'Missense',
            'Missense_Mutation': 'Missense', 'Frame_Shift_Del': 'Truncation','Nonsense_Mutation': 'Truncation', 
            'Splice_Site': 'Truncation', 'Frame_Shift_Ins': 'Truncation','Nonstop_Mutation':'Truncation'}

    mutations_replaced_M_T = somatic_mutations.replace(missence_truncation_groups)
    mutations_replaced_M_T = mutations_replaced_M_T.loc[mutations_replaced_M_T['Gene'] == gene]

    # group mutation categories
    miss = mutations_replaced_M_T.loc[mutations_replaced_M_T['Mutation'] == 'Missense']
    trunc = mutations_replaced_M_T.loc[mutations_replaced_M_T['Mutation'] == 'Truncation']

    #get lists of unique samples for missence and trucation categories
    miss_unique_samples = list(miss['Sample_ID'].unique())
    trunc_unique_samples = list(trunc['Sample_ID'].unique())
    
    #check if there is only one type of mutation for the specific gene
    if miss_unique_samples == []:
        print('Only truncation type mutations found for', gene+'.', 
             'Not possible to compare missense with wildtype.')
        truncation_omics = tumors.loc[tumors.index.isin(trunc_unique_samples)]
        truncation_omics = truncation_omics.assign(binary_mutations = 'Truncation')
        columns_to_drop = [gene+"_Mutation", gene+"_Location", gene+"_Mutation_Status", "Sample_Status"]
        binary_mut_omics = truncation_omics.drop(columns_to_drop, axis = 1)
        return binary_mut_omics
    elif trunc_unique_samples == []:
        print('Only missence type mutations found for', gene+'.', 
             'Not possible to compare truncation with wildtype.')
        missence_omics = tumors.loc[tumors.index.isin(miss_unique_samples)]
        missence_omics = missence_omics.assign(binary_mutations = 'Missense')
        columns_to_drop = [gene+"_Mutation", gene+"_Location", gene+"_Mutation_Status", "Sample_Status"]
        binary_mut_omics = missence_omics.drop(columns_to_drop, axis = 1)
        return binary_mut_omics

    # Step 2 - Create the binary column needed to do the comparison
    # Get mutation catagories with omics data
    missence_omics = tumors.loc[tumors.index.isin(miss_unique_samples)]
    missence_omics = missence_omics.assign(binary_mutations = 'Missense')
    truncation_omics = tumors.loc[tumors.index.isin(trunc_unique_samples)]
    truncation_omics = truncation_omics.assign(binary_mutations = 'Truncation')
    binary_mut_omics = missence_omics.append(truncation_omics)

    # Step 3 - Format the dataframe correctly for the T-test(just omics and binary columns for tumors)
    columns_to_drop = [gene+"_Mutation", gene+"_Location", gene+"_Mutation_Status", "Sample_Status"]
    binary_mut_omics = binary_mut_omics.drop(columns_to_drop, axis = 1)

    return binary_mut_omics


# PHOSPHOPROTEOMICS

def format_phospho_cis_comparison(cancer_object, omics_name, gene, specific_phospho):
    import numpy as np
    # Step 1 - Create dataframe in order to do comparisons with wrap_ttest
    omics_and_mut = cancer_object.join_omics_to_mutations(
        mutations_genes = gene, omics_df_name = omics_name, omics_genes = gene)

    # Step 2 - Create the binary column needed to do the comparison
    omics_and_mut['binary_mutations'] = omics_and_mut[gene+'_Mutation_Status'].apply(
        lambda x: 'Wildtype' if x == 'Wildtype_Tumor' else 'Mutated')

    # Step 3 - Format the dataframe correctly for the T-test(just omics and binary columns for tumors)
    tumors = omics_and_mut.loc[omics_and_mut['Sample_Status'] == 'Tumor'] #drop Normal samples
    binary_phospho = tumors[[specific_phospho, 'binary_mutations']].dropna(axis = 0)

    return binary_phospho


def get_missence_truncation_phospho(cancer_object, omics_name, gene, specific_phospho):
    import numpy as np
    #get omics data and tumors
    omics_and_mutations = cancer_object.join_omics_to_mutations(
                mutations_genes = gene, omics_df_name = omics_name, omics_genes = gene)
    tumors = omics_and_mutations.loc[omics_and_mutations['Sample_Status'] == 'Tumor'] #drop Normal samples

    somatic_mutations = cancer_object.get_mutations().reset_index()

    if cancer_object.get_cancer_type() == 'colon':
        missence_truncation_groups = {'frameshift substitution': 'Truncation', 
            'frameshift deletion': 'Truncation', 'frameshift insertion': 'Truncation', 
            'stopgain': 'Truncation', 'stoploss': 'Truncation', 'nonsynonymous SNV': 'Missence',
            'nonframeshift insertion': 'Missence','nonframeshift deletion': 'Missence', 
            'nonframeshift substitution': 'Missence'}
    else: 
        missence_truncation_groups = {'In_Frame_Del': 'Missence', 'In_Frame_Ins': 'Missence',
            'Missense_Mutation': 'Missence', 'Frame_Shift_Del': 'Truncation','Nonsense_Mutation': 'Truncation', 
            'Splice_Site': 'Truncation', 'Frame_Shift_Ins': 'Truncation','Nonstop_Mutation':'Truncation'}

    mutations_replaced_M_T = somatic_mutations.replace(missence_truncation_groups)
    mutations_replaced_M_T = mutations_replaced_M_T.loc[mutations_replaced_M_T['Gene'] == gene]

    # group mutation categories
    miss = mutations_replaced_M_T.loc[mutations_replaced_M_T['Mutation'] == 'Missence']
    trunc = mutations_replaced_M_T.loc[mutations_replaced_M_T['Mutation'] == 'Truncation']

    #get lists of unique samples for missence and trucation categories
    miss_unique_samples = list(miss['Sample_ID'].unique())
    trunc_unique_samples = list(trunc['Sample_ID'].unique())
    
    #check if there is only one type of mutation for the specific gene
    if miss_unique_samples == []:
        print('Only truncation type mutations found for', gene+'.', 
             'Not possible to compare mutation types.')
        truncation_omics = tumors.loc[tumors.index.isin(trunc_unique_samples)]
        truncation_omics = truncation_omics.assign(binary_mutations = 'Truncation')
        columns_to_drop = [gene+"_Mutation", gene+"_Location", gene+"_Mutation_Status", "Sample_Status"]
        binary_mut_omics = truncation_omics.drop(columns_to_drop, axis = 1)
        return binary_mut_omics
    elif trunc_unique_samples == []:
        print('Only missence type mutations found for', gene+'.', 
             'Not possible to compare mutation types.')
        missence_omics = tumors.loc[tumors.index.isin(miss_unique_samples)]
        missence_omics = missence_omics.assign(binary_mutations = 'Missence')
        columns_to_drop = [gene+"_Mutation", gene+"_Location", gene+"_Mutation_Status", "Sample_Status"]
        binary_mut_omics = missence_omics.drop(columns_to_drop, axis = 1)
        return binary_mut_omics

    # Step 2 - Create the binary column needed to do the comparison
    # Get mutation catagories with omics data
    missence_omics = tumors.loc[tumors.index.isin(miss_unique_samples)]
    missence_omics = missence_omics.assign(binary_mutations = 'Missence')
    truncation_omics = tumors.loc[tumors.index.isin(trunc_unique_samples)]
    truncation_omics = truncation_omics.assign(binary_mutations = 'Truncation')
    binary_mut_omics = missence_omics.append(truncation_omics)

    # Step 3 - Format the dataframe correctly for the T-test(just omics and binary columns for tumors)
    columns_to_drop = [gene+"_Mutation", gene+"_Location", gene+"_Mutation_Status", "Sample_Status"]
    binary_mut_omics = binary_mut_omics.drop(columns_to_drop, axis = 1)
    
    # select specific col and drop nan rows
    binary_phospho = binary_mut_omics[[specific_phospho, 'binary_mutations']].dropna(axis = 0)

    return binary_phospho