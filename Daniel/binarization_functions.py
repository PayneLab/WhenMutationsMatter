#Import dependencies
import numpy as np
import pandas as pd

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
        print('No significant results in '+attribute+'\n')
        
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