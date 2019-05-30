import pandas as pd
import numpy as np
import scipy.stats

# Making outliers tables

agg_col = "agg_col"

def convertToOutliers(df, samples, NUM_IQRs, up_or_down):
    df['row_iqr'] = scipy.stats.iqr(df[samples], axis=1, nan_policy='omit')
    df['row_median'] = np.nanquantile(df[samples], q=0.5, axis=1)

    outlier_df = pd.DataFrame()

    if up_or_down == 'up':
        df['row_medPlus'] = (df['row_median'] + (NUM_IQRs*df['row_iqr']))
        outlier_df = pd.DataFrame()
        outlier_df[samples] = df[samples].gt(df['row_medPlus'], axis=0).astype(int)

    elif up_or_down == 'down':
        df['row_medMinus'] = (df['row_median'] - (NUM_IQRs*df['row_iqr']))
        outlier_df = pd.DataFrame()
        outlier_df[samples] = df[samples].lt(df['row_medMinus'], axis=0).astype(int)

    outlier_df[df[samples].isnull()] = np.nan

    return outlier_df


def convertToCounts(df, samples, aggregate):
    '''
    If aggregate is set to True, adds up outliers in multiple rows,
    if there are repeated values in agg_col.
    If aggregation is set to False, just sets up the correct columns to feed
    into the comparison generator.
    '''
    not_outlier_cols = [x + '_notOutliers' for x in samples]
    outlier_cols = [x + '_outliers' for x in samples]

    if aggregate:
        df[agg_col] = [ind.split('-')[0] for ind in df.index]
        output_df = pd.DataFrame()
        output_df[not_outlier_cols] = df.groupby(by=agg_col)[samples].agg(lambda x: pd.Series(x == 0).sum())
        output_df[outlier_cols] = df.groupby(by=agg_col)[samples].agg(lambda x: pd.Series(x == 1).sum())
    elif not aggregate:
        output_df = pd.DataFrame(index=df.index)
        output_df[outlier_cols] = df[samples]
        output_df[not_outlier_cols] = 1 - df[samples]
    return output_df


def makeFracTable(df, samples):
    cols_outliers = [x + '_outliers' for x in samples]
    cols_notOutliers = [x + '_notOutliers' for x in samples]
    df = df.fillna(0)
    num_total_psites = df[cols_notOutliers].values + df[cols_outliers].values
    with np.errstate(divide='ignore', invalid='ignore'):
        frac_outliers = df[cols_outliers].values / num_total_psites

    frac_outliers = pd.DataFrame(frac_outliers, index=df.index, columns=samples)

    return frac_outliers


def make_outliers_table(df, iqrs=1.5, up_or_down='up', aggregate=True, frac_table=False):
    df = df.transpose()
    samples = df.columns

    df = convertToOutliers(df, samples, iqrs, up_or_down)

    df = convertToCounts(df, samples, aggregate)

    if frac_table:
        fracTable = makeFracTable(df, samples)
        return df.transpose(), fracTable.transpose()

    return df.transpose()


def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):

    from numpy import array, empty
    pvalues = array(pvalues)
    n = sum(~np.isnan(pvalues))
    new_pvalues = empty(len(pvalues))

    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues

    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values = [x for x in values if ~np.isnan(x[0])]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":

        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values = [x for x in values if ~np.isnan(x[0])]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]

    new_pvalues[np.isnan(pvalues)] = np.nan
    for i, val in enumerate(new_pvalues):
        if ~np.isnan(val):
            new_pvalues[i] = min(1, val)

    return new_pvalues


# Comparing groups
def getSampleLists(annotations, col):
    groups = list(pd.Series(annotations[col].value_counts().keys()).dropna())

    group0 = list(annotations.loc[annotations[col]==groups[0], :].index)
    group1 = list(annotations.loc[annotations[col]==groups[1], :].index)
    return groups[0], group0, groups[1], group1


def filterOutliers(df, group0_list, group1_list, frac_filter):

    group0_outliers = [x +'_outliers' for x in group0_list]
    group0_notOutliers = [x +'_notOutliers' for x in group0_list]
    group1_outliers = [x +'_outliers' for x in group1_list]
    group1_notOutliers = [x +'_notOutliers' for x in group1_list]

    if frac_filter!=None:
        min_num_outlier_samps = len(group0_list) * frac_filter
        num_outlier_samps = (df[group0_outliers] > 0).sum(axis=1)
        df = df.loc[num_outlier_samps >= min_num_outlier_samps, :]

    # Filter for higher proportion of outliers in group0 than group1
    group0_outlier_rate = df[group0_outliers].sum(axis=1).divide(df[group0_outliers+group0_notOutliers].sum(axis=1), axis=0)
    group1_outlier_rate = df[group1_outliers].sum(axis=1).divide(df[group0_outliers+group1_notOutliers].sum(axis=1), axis=0)

    df = df.loc[group0_outlier_rate>group1_outlier_rate, :]

    return df


def testDifferentGroupsOutliers(group0_list, group1_list, outlier_table):

    outliers_group0_list = [x +'_outliers' for x in group0_list]
    notOutliers_group0_list = [x +'_notOutliers' for x in group0_list]
    outliers_group1_list = [x +'_outliers' for x in group1_list]
    notOutliers_group1_list = [x +'_notOutliers' for x in group1_list]

    outlier_table['Outlier0'] = outlier_table[outliers_group0_list].sum(axis=1)
    outlier_table['NotOutlier0'] = outlier_table[notOutliers_group0_list].sum(axis=1)
    outlier_table['Outlier1'] = outlier_table[outliers_group1_list].sum(axis=1)
    outlier_table['NotOutlier1'] = outlier_table[notOutliers_group1_list].sum(axis=1)

    outlier_table['fisherp'] = outlier_table.apply((lambda r: scipy.stats.fisher_exact([[r['Outlier0'],
                                                                                     r['Outlier1']],
                                                                                    [r['NotOutlier0'],
                                                                                     r['NotOutlier1']]])[1]),axis=1)

    outlier_table['fisherFDR'] = correct_pvalues_for_multiple_testing(list(outlier_table['fisherp']),
                                     correction_type = "Benjamini-Hochberg")
    return outlier_table['fisherFDR']


def compare_groups_outliers(outliers, annotations, frac_filter=0.3):
    outliers = outliers.transpose()
    results_df = pd.DataFrame(index=outliers.index)
    for comp in annotations.columns:
        group0_label, group0, group1_label, group1 = getSampleLists(annotations, comp)
        label0 = '%s_%s_enrichment_FDR' %(comp, group0_label)
        df = filterOutliers(outliers, group0, group1, frac_filter)
        if len(df) > 0:
            print("Testing %s rows for enrichment in %s %s samples" %(len(df), comp, group0_label))
            col = pd.DataFrame(testDifferentGroupsOutliers(group0, group1, df))
            col.columns = [label0]
            results_df = pd.concat([results_df, col], axis=1, join='outer', sort=False)
        else:
            print("No rows had outliers in at least %s of %s %s samples" %(frac_filter, comp, group0_label))

        label1 = '%s_%s_enrichment_FDR' % (comp, group1_label)
        df = filterOutliers(outliers, group1, group0, frac_filter)
        if len(df) > 0:
            print("Testing %s rows for enrichment in %s %s samples" % (len(df), comp, group1_label))
            col = pd.DataFrame(testDifferentGroupsOutliers(group0, group1, df))
            col.columns = [label1]
            results_df = pd.concat([results_df, col], axis=1, join='outer', sort=False)
        else:
            print("No rows had outliers in at least %s of %s %s samples" % (frac_filter, comp, group1_label))

    return results_df
