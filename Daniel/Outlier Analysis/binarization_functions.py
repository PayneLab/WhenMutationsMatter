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
        #print('1 significant protein enrichment in '
        #      +attribute+':\n')
        
        return(sig_results)
    
    else:
        #print(str(len(sig_results))
        #      +' significant protein enrichments in '
        #      +attribute+'\n')
        
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

def get_dgidb_parameters():
    """
    This function takes no parameters, and will list possible options
    for querying the DGIdb with the dgidb_get_request() function. 
    dgidb_list_parameters() will be called in dgidb_get_request() if 
    incorrect parameter values are given.
    """
        
    print('\nRequired Parameters:')
    print('genes_or_drugs_list (list)')
    print('\nOne of the following:')
    print('genes (bool) OR drugs (bool)\n')
    
    print('Optional Parameters:')
    print('interaction_sources (list (case-sensitive): ["DrugBank","PharmGKB","TALC","TEND","TTD"])')
    print('interaction_types (list (not case-sensitive): ["activator", "inhibitor", "unknown"])')
    print('fda_approved_drug (bool)')
    print('immunotherapy (bool)')
    print('anti_neoplastic (bool)')
    print('clinically_actionable (bool)')
    print('druggable_genome (bool)')
    print('drug_resistance (bool)')
    print('gene_categories (list (not case-sensitive): ["kinase", "dna repair", "tumor suppressor"])')
    print('source_trust_levels (list (not case-sensitive): ["Expert curated", "Non-curated"])')
    
    print('\nFor more information on function parameters and DGIdb, see:')
    print('http://www.dgidb.org/api')
    
    return
      
def dgidb_get_request(genes_or_drugs_list,
                      genes=True, drugs=False,
                      interaction_sources=[], 
                      interaction_types=[], 
                      fda_approved_drug=False, 
                      immunotherapy=False, 
                      anti_neoplastic=False, 
                      clinically_actionable=False, 
                      druggable_genome=False, 
                      drug_resistance=False, 
                      gene_categories=[], 
                      source_trust_levels=[]):
    """
    Input Parameters (parameter descriptions are from www.dgidb.org/api)
        genes_or_drugs_list (list; Required):
            A comma delimited list of gene or drug names/symbols. 
            
        genes=True (bool: NOTE, either drugs or genes must be included):
            If this value is set to True, the corresponding genes_list
            will be created into a string of values and connected to the
            HTTP request.
        
        drugs=False (bool: NOTE, either drugs or genes must be included):
            If this value is set to True, the corresponding drugs_list
            will be created into a string of values and connected to the
            HTTP request. If both drugs and genes are set to True, you will
            run into errors. It would be best to query for drugs and genes
            separately, rather than in the same query.
            
        interaction_sources=[] (list; optional):
            A comma delimited list of interaction types to include in the result set. 
            If this field is omitted, all interaction types will be included.
            Valid Inputs: ["DrugBank","PharmGKB","TALC","TEND","TTD"]
            
        interaction_types=[] (list; optional):
            A comma delimited list of interaction types to include in the result set. If 
            this field is omitted, all interaction types will be included.
            Valid inputs: ["activator", "inhibitor", "unknown"]
        
        fda_approved_drug=False (bool: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            fda-approved drugs. If this field is omitted, interactions for all types of
            drugs will be included.
            
        immunotherapy=False (bool: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            immunotherapeutic drugs. If this field is omitted, interactions for all types
            of drugs will be included.
            
        anti_neoplastic=False (bool: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            antineoplastic drugs. If this field is omitted, interactions for all types of 
            drugs will be included.
            
        clinically_actionable=False (bool: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            clinically actionable genes. If this field is omitted, interactions for all 
            types of genes will be included.
            
        druggable_genome=False (bool: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            the durggable genome. If this field is omitted, interactions for all types of 
            genes will be included.
            
        drug_resistance=False (bool: optional):
            A flag denoting whether or not to limit interactions to only the ones involving 
            drug-resistant genes. If this field is omitted, interactions for all types of 
            genes will be included.
            
        gene_categories=[] (list: optional):
            A comma delimited list of gene categories to include in the result set. 
            If this field is omitted, all gene categories will be included.
            Valid inputs (not case-sensitive): ["kinase", "dna repair", "tumor suppressor"]
        
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
    valid_gene_categories = ["kinase", "dna repair", "tumor suppressor"]
    valid_source_trust_levels = ["expert curated", "non-curated"]
    try:
    
        if genes == True:
            url += 'genes='
            length = len(genes_or_drugs_list)
            for i, gene in enumerate(genes_or_drugs_list):
                if '(' in gene:
                    genes_or_drugs_list[i] = gene[:-3]     

            url += ','.join(genes_or_drugs_list)

        elif drugs == True:
            url += 'drugs='
            url += ','.join(genes_or_drugs_list)
        
        if (drugs == False and genes == False):
            raise Exception("'drugs' and 'genes' are both False. Please specify if genes_or_drugs is a list of genes or drugs.")
            
        elif (drugs == True and genes == True):
            raise Exception("'drugs' and 'genes' are both True. Please specify if genes_or_drugs is a list of genes or drugs.")
          
        if len(interaction_sources) >= 1:
            url += '&interaction_sources='
            for source in interaction_sources:
                if source not in valid_interaction_sources:
                    raise Exception("Invalid input interaction_sources: {}".format(source))
            url += ','.join(interaction_sources)
                    
        if len(interaction_types) >= 1:
            url += '&interaction_types='
            for types in interaction_types:
                if types not in valid_interaction_types:
                    raise Exception("Invalid input interaction_types: {}".format(types))
            url += ','.join(interaction_types)

        if fda_approved_drug == True:
            url += '&fda_approved_drug=true'
            
        if immunotherapy == True:
            url += '&immunotherapy=true'
            
        if anti_neoplastic == True:
            url += '&anti_neoplastic=true'
        
        if clinically_actionable == True:
            url += '&clinically_actionable=true'
            
        if druggable_genome == True:
            url += '&druggable_genome=true'
            
        if drug_resistance == True:
            url += '&drug_resistance=true'
         
        if len(gene_categories) >= 1:
            url += '&gene_categories='
            for i, category in enumerate(gene_categories):
                element = category.lower()
                if element not in valid_gene_categories:
                    raise Exception("Invalid gene category: {}".format(category))
                gene_categories[i] = element
                gene_categories[i] = '%20'.join(gene_categories[i])
            url += ','.join(gene_categories)
                
                    
        length = len(source_trust_levels)
        if length == 1:
            element = source_trust_levels[0].lower()
            if element not in valid_source_trust_levels:
                raise Exception("Invalid input source_trust_levels: {}".format(source_trust_levels[0]))
            if element == "expert curated":
                url += '&source_trust_levels=Expert%20curated'
            else:
                url += '&source_trust_levels=Non-curated'
        
        elif length > 1:
            url += '&source_trust_levels='
            for item in source_trust_levels:
                if item.lower() not in valid_source_trust_levels:
                    raise Exception("Invalid source_trust_level: {}".format(item))
                
            url += 'Expert%20curated,Non-curated'
                    
            
        #print("This is the full URL to your GET request:")
        #print(url)
        r = requests.get(url)
        #print("\nSee www.dgidb.org/api or this function's docstring for further explanation and resources on valid parameter inputs")
        
        return r.json()
    
    except Exception as e:
        print("Invalid parameter value: ")
        print(e)
        get_dgidb_parameters()
        
        return
    
def dgidb_json_parse(json_obj, genes=False, drugs=False):
    if (drugs == False and genes == False) or (drugs == True and genes == True):
        print("'genes' or 'drugs' must be set to True.")
        
        return
    
    elif genes == True:
        json_obj = json_obj['matchedTerms']
        drugs_dict = {}

        for item in json_obj:
            if len(item['interactions']) > 0:
                interactions_dict = {}
                for interaction in item['interactions']:
                    interactions_dict[interaction['drugName']] = interaction['interactionTypes']
                drugs_dict[item['geneName']] = interactions_dict

        if len(drugs_dict) == 0:
            #print('No Gene/Drug interactions')

            return

        return drugs_dict
    
    elif drugs == True:
        json_obj = json_obj['matchedTerms']
        genes_dict = {}

        for item in json_obj:
            if len(item['interactions']) > 0:
                interactions_dict = {}
                for interaction in item['interactions']:
                    interactions_dict[interaction['geneName']] = interaction['interactionTypes']
                genes_dict[item['drugName']] = interactions_dict

        if len(genes_dict) == 0:
            #print('No Drug/Gene interactions')

            return

        return genes_dict
    
def compare_enrichments_with_drugs(outliers_dict, clinical_df):
    patient_drugs_genes = {}
    for key in outliers_dict.keys():
        if key in clinical_df.index:
            comparisons = {}
            if isinstance(clinical_df.loc[key]['patient_medications'], str):
                drug_list = clinical_df.loc[key]['patient_medications'].split(',')
                for i, drug in enumerate(drug_list):
                    paren = drug.find('(')
                    if paren != -1:
                        drug_list[i] = drug[0:paren]
                comparisons['current_medications'] = drug_list
                patient_drugs_genes[key] = comparisons
                query = dgidb_get_request(genes=False, drugs=True,
                                          genes_or_drugs_list = drug_list)
                parsed = dgidb_json_parse(query, drugs=True)
                check = False
                if parsed is None:
                    continue
                for k,v in parsed.items():
                    for i, j in v.items():
                        if i in outliers_dict[key]:
                            comparisons[k] = {i : j}
                            check = True
                if check is False:
                    patient_drugs_genes['Gene/Drug Interactions'] = "No connection between current drugs and significantly enriched genes"
    print(json.dumps(patient_drugs_genes, indent=4))
    return patient_drugs_genes
