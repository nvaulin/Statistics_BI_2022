import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
from os.path import join
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.multitest import multipletests


def parsing_args():
    parser = argparse.ArgumentParser(
        usage='DGE_statistical_analysis.py [path_to_first_input] [path_to_second_input] [output_prefix] [method]',
        description='''Specification of the location of input and output files and a method of p-value correction
         in the case of multiple testing''')

    method_help = '''Method of p-values adjustment for multiple testing.
    The method can be None for no correction of one of those available in `statsmodels.stats.multitest.multipletests`:
    bonferroni or b, sidak or s, holm-sidak or hs, holm or h,
    simes-hochberg or sh, hommel, fdr_bh, fdr_by, fdr_tsbh, fdr_tsbky
    '''

    parser.add_argument('input_1', type=argparse.FileType('r'),
                        help='Path to the first input file with gene expressions')
    parser.add_argument('input_2', type=argparse.FileType('r'),
                        help='Path to the second input file with gene expressions')
    parser.add_argument('output', type=str, help='Output file name prefix')
    parser.add_argument('-m', type=str, help=method_help, metavar="method", default=None)

    args = parser.parse_args()

    output = args.output
    method = args.method
    expressions_1 = pd.read_csv(args.input_1)
    expressions_2 = pd.read_csv(args.input_2)

    return output, expressions_1, expressions_2, method


def plotter(expressions_1, expressions_2):
    """
    Provides plots generation for each gene asked via terminal
    """
    genes = []
    gene = input("Please type gene you are interested: ")
    genes.append(gene)
    while True:
        gene = input("Any other gene? Or type 'no': ")
        if gene == 'no':
            break
        else:
            genes.append(gene)
    for gene in genes:
        df_plot = pd.concat([expressions_1[[gene]], expressions_2[[gene]]], axis=1)
        df_plot.columns = ['Cell_type_1', 'Cell_type_2']
        hist_plotter(df_plot, gene)
        boxes_plotter(df_plot, gene)


def hist_plotter(df_plot, gene):
    plt.figure()
    exp_hist = sns.histplot(df_plot, stat="density")
    exp_hist.set_title(f'{gene} expression')
    exp_hist.legend(labels=['Cell_type_1', 'Cell_type_2'])
    exp_hist = exp_hist.get_figure()
    path = join("Pictures", f"Expressions_hist_{gene}.png")
    exp_hist.savefig(path)


def boxes_plotter(df_plot, gene, sample_size=250, n_samples=1000):
    mean_expr_1 = expressions_mean(df_plot[['Cell_type_1']], sample_size, n_samples)
    mean_expr_2 = expressions_mean(df_plot[['Cell_type_2']], sample_size, n_samples)
    mean_expr = [mean_expr_1, mean_expr_2]
    plt.figure()
    meanboxes = sns.boxplot(data=mean_expr)
    meanboxes.set_title(f'Mean {gene} expression')
    meanboxes.set_xticklabels(("Cell_type_1", "Cell_type_2"))
    meanboxes = meanboxes.get_figure()
    path = join("Pictures", f"Expressions_boxplots_{gene}.png")
    meanboxes.savefig(path)


def expressions_mean(expressions, sample_size=250, n_samples=1000):
    mean_expressions = []
    for n in range(n_samples):
        mean_expressions.append(expressions.sample(sample_size).mean())
    return mean_expressions


def calculate_ci(sample, std=None):
    if std is None:
        std = sample.std()
    mean = sample.mean()
    se = std / np.sqrt(len(sample))
    left_b = mean - 1.96 * se
    right_b = mean + 1.96 * se
    return (left_b, right_b)


def check_intervals_intersection(intervals):
    first_ci = intervals[0]
    second_ci = intervals[1]
    return (first_ci[1] <= second_ci[0]) or (second_ci[1] <= first_ci[0])


def table_ci_preproc(table, celltype):
    """
    Preprocesses data for CI-test (including CI evaluation)
    :param table:
    :param celltype:
    :return:
    """
    table = table.iloc[:, :-1].apply(calculate_ci, axis=0).transpose()
    table[f'{celltype}_interval'] = list(zip(table[0], table[1]))
    table.drop([0, 1], axis=1, inplace=True)
    return table


def table_z_preproc(table, celltype):
    """
    Preprocesses data for Z-test
    """
    table = table.iloc[:, :-1].transpose()
    table[f'{celltype}'] = table.values.tolist()
    table.drop(list(range(0, 500)), axis=1, inplace=True)
    return table


def check_dge_with_ci(first_table, second_table, celltypes=("Cell_type_1", "Cell_type_2")):
    first_table = table_ci_preproc(first_table, celltypes[0])
    second_table = table_ci_preproc(second_table, celltypes[1])
    intervals = pd.concat([first_table, second_table], axis=1)
    ci_test_results = intervals.apply(check_intervals_intersection, axis=1)
    return ci_test_results


def check_dge_with_ztest(first_table, second_table, method, celltypes=("Cell_type_1", "Cell_type_2")):
    first_table = table_z_preproc(first_table, celltypes[0])
    second_table = table_z_preproc(second_table, celltypes[1])
    dge_ztest = pd.concat([first_table, second_table], axis=1)
    z_test_results = dge_ztest.apply(lambda row: ztest(row.Cell_type_1, row.Cell_type_2), axis=1)
    dge_ztest[['z_score', 'p_value']] = pd.DataFrame(z_test_results.tolist(), index=z_test_results.index)
    if method is None:
        dge_ztest['Significant_differences'] = dge_ztest['p_value'] < 0.05
        return dge_ztest['p_value'], None, dge_ztest['Significant_differences']
    else:
        dge_ztest['p_value_adj'] = list(multipletests(dge_ztest['p_value'], method=method)[1])
        dge_ztest['Significant_differences'] = dge_ztest['p_value_adj'] < 0.05
        return dge_ztest['p_value'], dge_ztest['p_value_adj'], dge_ztest['Significant_differences']


if __name__ == '__main__':
    output, expressions_1, expressions_2, method = parsing_args()
    to_plot = input("Do you want to get some expression plots for any gene? [y/n]: ")
    if to_plot == 'y':
        plotter(expressions_1, expressions_2)

    ci_test_results = check_dge_with_ci(expressions_1, expressions_2)
    z_test_p_values, z_test_p_values_adj, z_test_results = check_dge_with_ztest(expressions_1, expressions_2, method)
    mean_diff = expressions_1.mean(axis=0) - expressions_2.mean(axis=0)

    results = {
        "mean_diff": mean_diff,
        "z_test_results": z_test_results,
        "z_test_p_values": z_test_p_values
    }

    if method is not None:
        results['z_test_p_values_adj'] = z_test_p_values_adj
    else:
        results['ci_test_results'] = ci_test_results

    results = pd.DataFrame(results)
    results.to_csv(f"{output}.csv")
