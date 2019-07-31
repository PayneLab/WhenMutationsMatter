#Import dependencies
import numpy as np
import pandas as pd
import requests
import json

#Binarization Functions
def binarizeCutOff(df, column, cut_off, replace_high, replace_low):
    """Input parameters:
           df:
               The Clinical DataFrame
               
           column: 
               A column in df to binarize
           
           cut_off:
               The benchmark at which you compare the original
               samples in the dataframe. For example, you 
               may set a cut_off of 65 for the clinical['Age'] 
               column to separate Retired and Working individuals. 
           
           replace_low: 
               The value that will replace original values
               lower than cut_off in the chosen column.
           
           replace_high: 
               The value that will replace original values
               greater than cut_off in the chosen column.
       
       Method description:
           This should replace the original value in the 
           clinical dataframe with replace_low if the original
           value is lower than cut_off. Otherwise, it will
           replace the original value with replace_high if the
           original value is greater than or equal to cut_off.
           Because np.where() overwrites NaN values, they will
           be saved by Sample index, and then put back in at
           the end using pd.DataFrame.where() 
       
       Return Value:
           A one column dataframe with the specified, replaced values.
    """
    
    new_df = df.copy()
    final_df = df.copy()
    nulls = new_df[column].isnull()
    
    new_df[column] = np.where(new_df[column] >= cut_off, replace_high, replace_low)
    final_df[column] = new_df[column].where(~nulls, np.nan)
    
    return(final_df[column])


def binarizeRange(df, column, lower_bound, upper_bound, in_range='In_Range', out_of_range='Out_Of_Range'):
    """Input parameters:
           df:
               The Clinical DataFrame
           
           column: 
               A column in df to binarize
           
           lower_bound:
               The lowest value of a researchers's specified
               range. For example, a value of 18.5 for BMI
               could be the function's lower_bound.
           
           upper_bound: 
               The highest value of a researchers's specified
               range. For example, a value of 24.9 for BMI
               could be the function's upper_bound.
                      
           in_range (default value = 'In_Range'): 
               The value that will replace original values
               within the specified range between lower_bound
               and upper_bound.
               
           out_of_range (default value = 'Out_Of_Range'):
               The value that will replace original values
               outside of the specified range between
               lower_bound and upper_bound.
       
       Method description:
           This should replace the original value in the 
           clinical dataframe with in_range if the original
           value is within the range from lower_bound to 
           upper_bound. Otherwise, it will replace the 
           original value with out_of_range. Because np.where() 
           overwrites NaN values, they will be saved by Sample 
           index, and then put back in at the end using
           pandas.DataFrame.where().
       
       Return Value:
           The return value should be a list of values with
           two options, either in_range or out_of_range. This
           function is for columns with continuous variables, 
           such as integers or floats
    """
    
    new_df = df.copy()
    final_df = df.copy()
    nulls = new_df[column].isnull()

    new_df[column] = np.where((new_df[column] >= lower_bound) & 
                              (new_df[column] < upper_bound), 
                              in_range, out_of_range)
    final_df[column] = new_df[column].where(~nulls, np.nan)
    
    return(final_df[column])


def binarizeCategorical(df, column, dictionary):
    """
    Input Parameters:
        
        df:
            The Clinical DataFrame
            Example: df = en.get_clinical()
        
        column:
            A column in df to binarize
            Example: column = df['Race']
        
        dictionary: 
            The dictionary to use in mapping df[column], where
            keys are the responses to be replaced, and values 
            that will replace the keys.
            
            Example: 
            my_dictionary = {'White':'European', 
                             'Black or African American':'Not_European', 
                             'Not Reported':'Not_European', 
                             'Asian':'Not_European'}
            
    Method Description:
        This function takes a dataframe of clinical observations, and 
        a column containing categorical variables, and maps a dictionary 
        to replace values within that column.

    Return Value:
        A one column dataframe with the specified, replaced values.
    """
    new_df = df.copy()
    new_df[column] = df[column].map(dictionary).fillna(np.nan)
    
    return(new_df[column])


def significantEnrichments(df, column, p_value=0.05):
    """
    Input Parameters
        df:
            The DataFrame returned from the compare_groups_outliers
            function in cptac.Algorithms().
            
        column:
            The column within df that you want to filter for 
            significance. For example, 'Histologic_type_Serous_enrichment_FDR'.
        
        p_value:
            This is the cut off you choose to filter for significant
            enrichments. The default value of 0.05, for instance, will
            print out the proteins and their associated p-values ONLY
            for proteins whose p-values are less than 0.05.
            
    Methods Description:
        This function will drop NaN values from the column you choose to 
        filter, and then will filter that column based on the chosen p_value.
        It will then print out the number of significant proteins compared to
        the original number of proteins, and the p-values associated with each
        significant protein that are less than 0.05.
    
    Return Value:
        The DataFrame of significant proteins and their p-values will be returned.
    """
    #Get rid of '_enrichment_FDR' for nicer print statement
    words = column.split('_')
    attribute = []
    for item in words:
        if item == 'enrichment' or item == 'FDR':
            break
        else:
            attribute.append(item)
    attribute = "_".join(attribute)
    
    #Drop NaN values and filter by p_value
    results = df[column].dropna().copy()
    sig_results = results[results < p_value]
    sig_results = pd.DataFrame(sig_results)
    sig_results.columns = [attribute+'_P_values']
    
    if len(sig_results) < 1:
        
        return
    
    elif len(sig_results) == 1:
        print('1 significant protein enrichment in '
              +attribute+':\n')
        
        return(sig_results)
    
    else:
        print(str(len(sig_results))
              +' significant protein enrichments in '
              +attribute+'\n')
        
        return(sig_results)

def renameDuplicateColumns(outliers_df, dict_of_counts=False):
    #Dealing with the duplicate genes problem
    outliers_list = list(outliers_df.columns)
    if len(outliers_list) == len(set(outliers_list)):
        print("There are no duplicates")

        return

    else:
        duplicates = {} #Keys = gene name, Values = counts
        duplicates_list = []

        #Create duplicate list and populate duplicate dictionary
        for i in range(len(outliers_list) - 1):
            if outliers_list[i] == outliers_list[i + 1]:
                duplicates_list.append(outliers_list[i])
                duplicates[outliers_list[i]] = 1

        #Go through and check for duplicates, and then rename duplicates accordingly
        #Can this be made more efficient?
        for i in range(len(outliers_list) - 1):
            #Check if the current == the next
            if outliers_list[i] == outliers_list[i + 1]:
                count = 1
                outliers_list[i + count] += "({})".format(str(count))
            #Check if the next == the current before the (#)
            elif outliers_list[i + 1] == outliers_list[i][0:-3]:
                count += 1
                outliers_list[i + 1] += "({})".format(str(count))
            #Check if the next == the current before the (#) if there are greater than 9 duplicates
            elif outliers_list[i + 1] == outliers_list[i][0:-4]:
                count += 1
                outliers_list[i + 1] += "({})".format(str(count))
        
        #Renaming Columns to avoid errors
        outliers_df.columns = outliers_list

        if dict_of_counts == True:
            
            #Count how many duplicates of each gene there are and update counts in dict        
            for i in range(len(duplicates_list) - 1):
                if duplicates_list[i] == duplicates_list[i + 1]:
                    duplicates[duplicates_list[i]] += 1

            return duplicates
        
        return

def dgidb_list_parameters():
    """
    This function takes no parameters, and will list possible options
    for querying the DGIdb with the dgidb_get_request(). dgidb_list_parameters()
    will be called in dgidb_get_request() if incorrect parameter values are given.
    """
    
    param_list = ['genes_or_drugs_list, genes, drugs, interaction_sources, interaction_types, fda_approved_drug, immunotherapy, anti_neoplastic, clinically_actionable, druggable_genome, drug_resistance, gene_categories, source_trust_levels']
    
    print('\nFull list of parameters:')
    for param in param_list:
        print(param)
    
    print('\nRequired Parameters:')
    print('genes_or_drugs_list (list)')
    print('genes (string: "true" or "false")')
    print('drugs (string: "true" or "false"\n')
    
    print('Optional Parameters')
    print('interaction_sources (string: "TTD", "DrugBank", "TALC", "something")')
    print('interaction_types (string: "inhibitor", "activator")')
    print('fda_approved_drug (string: "true" or "false")')
    print('immunotherapy (string: "true" or "false")')
    print('anti_neoplastic (string: "true" or "false")')
    print('clinically_actionable (string: "true" or "false")')
    print('druggable_genome (string: "true" or "false")')
    print('drug_resistance (string: "true" or "false")')
    print('gene_categories (string: see www.dgidb.org/api for options)')
    print('source_trust_levels (string: see www.dgidb.org/api for options)')
    
    print('\nFor more information on function parameters and DGIdb, see:')
    print('http://www.dgidb.org/api')
    
    return
    
#GET request maker for the DGIdb    
def dgidb_get_request(genes_or_drugs_list,
                      genes='N/A', drugs='N/A',
                      interaction_sources=[], 
                      interaction_types=[], 
                      fda_approved_drug='N/A', 
                      immunotherapy='N/A', 
                      anti_neoplastic='N/A', 
                      clinically_actionable='N/A', 
                      druggable_genome='N/A', 
                      drug_resistance='N/A', 
                      gene_categories=[], 
                      source_trust_levels=[]):
    """
    Input Parameters (parameter descriptions are from www.dgidb.org/api)
        genes_or_drugs_list (list; Required):
            A comma delimited list of gene or drug names/symbols. 
            
        genes='N/A' (NOTE, either drugs or genes must be included):
            If this value is set to 'yes', the corresponding genes_list
            will be created into a string of values and connected to the
            HTTP request.
        
        drugs='N/A' (NOTE, either drugs or genes must be included):
            If this value is set to 'yes', the corresponding drugs_list
            will be created into a string of values and connected to the
            HTTP request.
            
        interaction_sources=[] (list; optional):
            A comma delimited list of interaction types to include in the result set. 
            If this field is omitted, all interaction types will be included.
            Valid Inputs: ["DrugBank","PharmGKB","TALC","TEND","TTD"]
            
        interaction_types=[] (list; optional):
            A comma delimited list of interaction types to include in the result set. If 
            this field is omitted, all interaction types will be included.
            Valid inputs: ["activator", "inhibitor", "unknown"]
        
        fda_approved_drug='N/A' (string: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            fda-approved drugs. If this field is omitted, interactions for all types of
            drugs will be included.
            Valid inputs: 'true', 'false'
            
        immunotherapy='N/A' (string: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            immunotherapeutic drugs. If this field is omitted, interactions for all types
            of drugs will be included.
            Valid inputs: 'true', 'false'
            
        anti_neoplastic='N/A' (string: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            antineoplastic drugs. If this field is omitted, interactions for all types of 
            drugs will be included.
            Valid inputs: 'true', 'false'
            
        clinically_actionable='N/A' (string: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            clinically actionable genes. If this field is omitted, interactions for all 
            types of genes will be included.
            Valid inputs: 'true', 'false'
            
        druggable_genome='N/A' (string: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            the durggable genome. If this field is omitted, interactions for all types of 
            genes will be included.
            Valid inputs: 'true', 'false'
            
        drug_resistance='N/A' (string: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            drug-resistant genes. If this field is omitted, interactions for all types of 
            genes will be included.
            Valid inputs: 'true', 'false'
            
        gene_categories=[] (list: optional):
            A comma delimited list of gene categories to include in the result set. 
            If this field is omitted, all gene categories will be included.
            Valid inputs (not case-sensitive): ["KINASE", "DNA REPAIR", "TUMOR SUPPRESSOR"]
        
        source_trust_levels=[] (list: optional):
            A comma delimited list of source trust levels to include in the result set. 
            If this field is omitted, all trust levels will be included.
            Valid inputs (case-sensitive): ["Expert curated", "Non-curated"]
            
    Methods Description:
        This function will take a list of either drugs or genes, and any other
        parameters to filter and create a specific HTTP GET request to the 
        Drug Gene Interaction Database (DGIdb). It is advised to use this function 
        after an Outlier analysis is performed and significant gene enrichments are 
        found. This will allow you to quickly  find drugs suitable for cancer
        treatment based on patients' gene enrichments for higher precision.
    
    Return Value:
        A Python requests object containing JSON with the information from your GET
        request to the DGIdb.
    """
    url = 'http://www.dgidb.org/api/v2/interactions.json?'
    valid_interaction_sources = ["DrugBank","PharmGKB","TALC","TEND","TTD"]
    valid_interaction_types = ["activator", "inhibitor", "unknown"]
    valid_gene_categories = ["KINASE", "DNA REPAIR", "TUMOR SUPPRESSOR"] #check api if case sensitive
    valid_source_trust_levels = ["Expert curated", "Non-curated"]
    try:
    
        if genes == 'yes':
            genes_string = 'genes='
            for gene in genes_or_drugs_list:
                if gene == genes_or_drugs_list[0]:
                    genes_string += gene
                else:
                    genes_string += ',' + gene

            url += genes_string

        elif drugs == 'yes':
            drugs_string = 'drugs='
            for drug in genes_or_drugs_list:
                if drug == genes_or_drugs_list[0]:
                    drugs_string += drug
                else:
                    drugs_string += ',' + drug
            url += drugs_string
        
        if drugs == "N/A" and genes == "N/A":
            raise Exception("genes_or_drugs is ambiguous. Please specify if it is a list of genes or drugs.")

        if interaction_sources == 'true' or interaction_sources == 'false':
            url += '&interaction_sources=' + interaction_sources
           
        elif interaction_sources != "N/A":
            raise Exception("Invalid input interaction_sources: {}".format(interaction_sources))

        if interaction_types == 'true' or interaction_types == 'false':
            url += '&interaction_types=' + interaction_types
            
        elif interaction_types != "N/A":
            raise Exception("Invalid input interaction_types: {}".format(interaction_types))

        if fda_approved_drug == 'true' or fda_approved_drug == 'false':
            url += '&fda_approved_drug=' + fda_approved_drug
        
        elif fda_approved_drug != "N/A":
            raise Exception("Invalid input: fda_approved_drug. Must be 'true' or 'false'. Got {}".format(fda_approved_drug))
            
        if immunotherapy == 'true' or immunotherapy == 'false':
            url += '&immunotherapy=' + immunotherapy
        
        elif immunotherapy != "N/A":
            raise Exception("Invalid input: immunotherapy. Must be 'true' or 'false'. Got {}".format(immunotherapy))
        
        if anti_neoplastic == 'true' or anti_neoplastic == 'false':
            url += '&anti_neoplastic=' + anti_neoplastic
        
        elif anti_neoplastic != "N/A":
            raise Exception("Invalid input: anti_neoplastic. Must be 'true' or 'false'. Got {}".format(anti_neoplastic))
        
        if clinically_actionable == 'true' or clinically_actionable == 'false':
            url += '&clinically_actionable=' + clinically_actionable
            
        elif clinically_actionable != "N/A":
            raise Exception("Invalid input: clinically_actionable. Must be 'true' or 'false'. Got {}".format(clinically_actionable))

        if druggable_genome == 'true' or druggable_genome == 'false':
            url += '&druggable_genome=' + druggable_genome
            
        elif druggable_genome != "N/A":
            raise Exception("Invalid input: druggable_genome. Must be 'true' or 'false'. Got {}".format(druggable_genome))

        if drug_resistance == 'true' or drug_resistance == 'false':
            url += '&drug_resistance=' + drug_resistance
        
        elif drug_resistance != "N/A":
            raise Exception("Invalid input: drug_resistance. Must be 'true' or 'false'. Got  {}".format(drug_resistance))
        
        if len(gene_categories == 1):
            url += '&gene_categories=' + gene_categories[0]
            
        elif len(set(gene_categories) > 1):
            unique_categories = set(gene_categories)
            for cat in unique_categories:
                if cat not in valid_gene_categories:
                    raise Exception("Invalid gene category: {}".format(cat))
                elif cat == unique_categories[0]:
                    url += '&gene_categories=' + cat + ","
                elif cat == unique_categories[-1]:
                    url += cat
                else:
                    url += cat + ","
                    
        if source_trust_levels != 'N/A': #continue refactoring from here
            url += '&source_trust_levels=' + source_trust_levels

       #Cannot feasibly call an error because of the numerous combinations this can take
        #Possible values available at www.dgidb.org/api
            
        print("This is the full URL to your GET request:")
        print(url)
        r = requests.get(url)
        print("If this request failed, it is likely due to invalid inputs for gene_categories or source_trust_levels")
        print("See www.dgidb.org/api for further explanation and resources for valid parameter inputs")
        
        return r
    
    except Exception as e:
        print("Invalid parameter value: ")
        print(e)
        dgidb_list_parameters()
        
        return