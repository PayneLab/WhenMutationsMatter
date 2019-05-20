from __future__ import division
import pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns
import matplotlib.patches as mpatches
import scipy.stats
import sys
import argparse
from datetime import datetime

def fileToList(group_list):
    with open(group_list, 'r') as fh:
        return [line.strip() for line in fh.readlines()]

def fileToDict(tsv_map_file_name):
    with open(tsv_map_file_name, 'r') as fh:
        return {line.split()[0]:line.split()[1] for line in fh.readlines()}

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

def assignColors(group_colors, group1_label, group2_label, group1, group2):
    if group_colors is not None:
        group_color_map = fileToDict(group_colors)

        groups_dict = {sample:group1_label for sample in group1}
        groups_dict2 = {sample:group2_label for sample in group2}
        groups_dict.update(groups_dict2)

        sample_color_map = {sample:group_color_map[groups_dict[sample]] for sample in group1+group2}

    else:
        sample_color_map = {sample:'#571D41' for sample in group1}
        sample_color_map.update({sample:'#F5F5F5' for sample in group2})

        group_color_map = {group1_label:'#571D41', group2_label:'#F5F5F5'}
    return group_color_map, sample_color_map

def filterOutliers(df, group1_list, group2_list, gene_column_name, frac_filter):

    group1_outliers = [x+'_outliers' for x in group1_list]
    group1_notOutliers = [x+'_notOutliers' for x in group1_list]
    group2_outliers = [x+'_outliers' for x in group2_list]
    group2_notOutliers = [x+'_notOutliers' for x in group2_list]

    if frac_filter!=None:
        min_num_outlier_samps = len(group1_list)*frac_filter
        # Filter for 30% of samples having outliers in group1
        num_outlier_samps = (df[group1_outliers] > 0).sum(axis=1)
        df = df.loc[num_outlier_samps >= min_num_outlier_samps, :].reset_index(drop=True)

    # Filter for higher proportion of outliers in group1 than group2
    group1_outlier_rate = df[group1_outliers].sum(axis=1).divide(df[group1_outliers+group1_notOutliers].sum(axis=1), axis=0)
    group2_outlier_rate = df[group2_outliers].sum(axis=1).divide(df[group1_outliers+group2_notOutliers].sum(axis=1), axis=0)

    df = df.loc[group1_outlier_rate>group2_outlier_rate, :].reset_index(drop=True)

    with np.errstate(divide='ignore',invalid='ignore'):
        frac_outliers = pd.concat([pd.DataFrame(df[group1_outliers].values/(df[group1_outliers].values+df[group1_notOutliers].values),
                                                columns=group1_list),
                                   pd.DataFrame(df[group2_outliers].values/(df[group2_outliers].values+df[group2_notOutliers].values),
                                                columns=group2_list)],
                                   axis=1).reset_index(drop=True)

    frac_outliers[gene_column_name] = df[gene_column_name]

    return df, frac_outliers

def testDifferentGroupsOutliers(group1_list, group2_list, outlier_table):

    outliers_group1_list = [x+'_outliers' for x in group1_list]
    notOutliers_group1_list = [x+'_notOutliers' for x in group1_list]
    outliers_group2_list = [x+'_outliers' for x in group2_list]
    notOutliers_group2_list = [x+'_notOutliers' for x in group2_list]

    outlier_table['Outlier1'] = outlier_table[outliers_group1_list].sum(axis=1)
    outlier_table['NotOutlier1'] = outlier_table[notOutliers_group1_list].sum(axis=1)

    outlier_table['Outlier2'] = outlier_table[outliers_group2_list].sum(axis=1)
    outlier_table['NotOutlier2'] = outlier_table[notOutliers_group2_list].sum(axis=1)

    outlier_table['fisherp'] = outlier_table.apply((lambda r: scipy.stats.fisher_exact([[r['Outlier1'],
                                                                                     r['Outlier2']],
                                                                                    [r['NotOutlier1'],
                                                                                     r['NotOutlier2']]])[1]),axis=1)

    outlier_table['fisherFDR'] = correct_pvalues_for_multiple_testing(list(outlier_table['fisherp']),
                                     correction_type = "Benjamini-Hochberg")
    return outlier_table['fisherFDR']

def makeHeatMap(heatmap_table, group_color_map,sample_color_map, group1,
                output_prefix, blue_or_red, genes_to_highlight=None):

    # Set up colors for heatmap
    sns.set(font = 'arial', style = 'white', color_codes=True, font_scale = 1)
    if blue_or_red == 'red':
        cmap = sns.cubehelix_palette(start=0.857, rot=0.00, gamma=1.5, hue=1,
                                    light=1, dark=0.2, reverse=False, as_cmap=True)
    elif blue_or_red=='blue':
        cmap=sns.cubehelix_palette(start=3, rot=0.00, gamma=1.5, hue=1,
                                    light=1, dark=0.2, reverse=False, as_cmap=True)
    cmap.set_bad('#F5F5F5')

    # Set up colors for columns
    heatmap_table.columns = [x.split('_')[0] for x in heatmap_table.columns]
    group = heatmap_table.columns
    column_colors = group.map(sample_color_map)

    # Sort the table
    heatmap_table['mean_group1'] = heatmap_table.loc[:, group1].mean(axis=1)
    heatmap_table = heatmap_table.sort_values('mean_group1')
    heatmap_table = heatmap_table.drop('mean_group1', axis=1)

    # Make heatmap
    g = sns.clustermap(heatmap_table,
                       cmap=cmap,
                       col_cluster = False,
                       row_cluster = False,
                       col_colors = column_colors,
                       xticklabels=False,
                       vmin=0,
                       vmax=1,
                       cbar_kws={'label':'fraction outliers'},
                       # figsize=(0.1*len(list(heatmap_table)), 0.05*len(heatmap_table))
                          )
    g.ax_row_dendrogram.set_visible(False)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    g.cax.set_position((0.15,0.12,0.03,0.6)) #move colorbar to right
    ax = g.ax_heatmap
    ax.set_ylabel('') #change the gene label

    # Make the legend the describes the different sample groups
    handles = [mpatches.Patch(color=color, label=group) for group, color in group_color_map.items()]

    fig = plt.gcf()
    fig.legend(handles=handles, bbox_to_anchor=(0.6, 0.10))

    # Star genes to highlight on heatmap
    if genes_to_highlight != None:
        new_labels = []
        for gene in ax.get_yticklabels():
            gene = gene.get_text()
            if gene in genes_to_highlight:
                new_labels.append(str(gene)+'*')
            else:
                new_labels.append(gene)
        ax.set_yticklabels(new_labels)

    # Save and show the plot
    plt.savefig('%s.pdf' %output_prefix, dpi=500, bbox_inches='tight', pad_inches=0.5)
    return g

def runComparison(outliers, group1, group2,gene_column_name,
                  fdr_cut_off=0.05, frac_filter=0.3):
# Filter for multiple samples with outliers in group1_list
    outliers, frac_outliers = filterOutliers(outliers, group1, group2, gene_column_name, frac_filter)
    if len(outliers) == 0:
        print("No rows have outliers in %s of group1 samples" % (frac_filter))
        sys.exit()
    print('Testing %s genes for enrichment in group1' %(len(outliers)))

# Doing statistical test on different groups
    outliers['FDR'] = testDifferentGroupsOutliers(group1, group2, outliers)

    outliers['significant'] = (outliers['FDR'] < fdr_cut_off)
    sig_diff_count = sum(outliers['significant'])
    print('%s signficantly differential genes at FDR %s' % (sig_diff_count, fdr_cut_off))
    return outliers[[gene_column_name, 'FDR']], frac_outliers, sig_diff_count

def writeSigGenesToFile(outliers, gene_column_name, output_prefix, group1_label):
    sig_genes = outliers.loc[(outliers['significant']==True), gene_column_name]
    with open('%s.outlier_sites_in_%s.txt' %(output_prefix, group1_label), 'w') as fh:
        for gene in sig_genes:
            fh.write('%s\n'%gene)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Input samples in two groups and df")
    parser.add_argument('--outliers_table', type=str,
                        help='tsv with proteins as rows and samples as columns. Phospho data also needs a "counts" column')
    parser.add_argument('--gene_column_name', type=str, default='geneSymbol',
                        help='Specify column name for protein IDs')
    parser.add_argument('--fdr_cut_off', type=float, default=0.05,
                        help='Set FDR cut off for which genes are considered sig diff')
    parser.add_argument('--output_prefix', type=str, default=str(datetime.now().date()),
                        help='Default is current date. Set for prefix for gene lists and heatmap')
    parser.add_argument('--group1_label', type=str, default='group1',
                        help='Label for group1 for heatmap')
    parser.add_argument('--group1_list', type=str,
                        help='List of samples from group 1. Samples separated by new lines, no header')
    parser.add_argument('--group2_label', type=str, default='group2',
                        help='Label for group2 for heatmap')
    parser.add_argument('--group2_list', type=str,
                        help='List of samples from group 2. Samples separated by new lines, no header')
    parser.add_argument('--group_colors', type=str,
                        help='tsv, no header, with group names in 1st column and colors in 2nd column. Alternatively can map colors directly to samples. Group labels must match group labels given above.')
    parser.add_argument('--genes_to_highlight', type=str, default=None,
                        help='List of genes to highlight with * in y tick labels')
    parser.add_argument('--blue_or_red', type=str, default='red',
                        help='Color scale for heatmap')
    parser.add_argument('--output_qvals', type=str, choices=['True', 'False'], default=False,
                        help='Save q-values to a table?')
    parser.add_argument('--frac_filter', type=str, default=0.3,
                        help='None or fraction (bn 0-1) of outlier samples in group1 needed per gene')

    args = parser.parse_args()

    outliers = pd.read_csv(args.outliers_table, sep='\t')
    gene_column_name = args.gene_column_name
    fdr_cut_off = args.fdr_cut_off
    output_prefix = args.output_prefix
    group1_label = args.group1_label
    group2_label = args.group2_label
    group1 = [x for x in fileToList(args.group1_list) if x+'_outliers' in outliers.columns]
    group2 = [x for x in fileToList(args.group2_list) if x+'_outliers' in outliers.columns]
    genes_to_highlight = args.genes_to_highlight
    if genes_to_highlight != None:
        genes_to_highlight = fileToList(genes_to_highlight)
    blue_or_red = args.blue_or_red
    output_qvals = args.output_qvals == 'True'
    frac_filter = args.frac_filter

    if frac_filter=='None':
        frac_filter = None
    else:
        try:
            frac_filter = float(frac_filter)
            if frac_filter > 1:
                print('Non-valid frac_filter value')
                sys.exit()
            if frac_filter < 0:
                print('Non-valid frac_filter value')
                sys.exit()
        except:
            print('Non-valid frac_filter value')
            sys.exit()


# Assigning colors to samples
    group_color_map, sample_color_map = assignColors(args.group_colors, group1_label, group2_label, group1, group2)

# Filter for multiple samples with outliers in group1_list
    outliers, frac_outliers, sig_diff_count = runComparison(outliers, group1, group2, gene_column_name, fdr_cut_off, frac_filter)

    if output_qvals == True:
        outliers.to_csv('%s_comparison_qvals.txt'%output_prefix, sep='\t', index=False)

#If enough genes, make heatmap
    if sig_diff_count >= 1:
        heatmap_table = frac_outliers.loc[(outliers['significant'] == True), [gene_column_name] + group1 + group2]
        heatmap_table = heatmap_table.set_index(gene_column_name)

        makeHeatMap(heatmap_table, group_color_map, sample_color_map, group1,
                    output_prefix, blue_or_red, genes_to_highlight)

#Write significantly different genes to a file
        writeSigGenesToFile(outliers, gene_column_name, output_prefix, group1_label)
