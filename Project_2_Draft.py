import pandas as pd
import numpy as np
import argparse
#import skbio.diversity.beta import pw_distances

def read_csv_file(in_File):
    input_File = in_File

    #df = pd.read_csv(input_File, skiprows = [0], index_col = 0)
    #return (df)

    col_names = ['s1_p1', 's1_p2', 's1_p3', 's1_p4', 's1_p5', 's1_p6', 's1_p7', \
    's1_p8', 's1_p9', 's1_p10', 's1_p11', 's1_p12']
    #df_data_only = pd.read_csv(input_File, skiprows = [0, 1], index_col = [0], \
    #names = col_names, header = None)
    df_data_only = pd.read_csv(input_File, skiprows = [0, 1], index_col = [0, 13], \
    names = col_names, header = None)
    return (df_data_only)

def set_test(in_df):
    result = set(in_df)
    return (result)

def argparse_practice():
    parser = argparse.ArgumentParser(description = 'Tenative description.')
    parser.add_argument('--foo', help = 'foo of the %(prog)s program')
    args = parser.parse_args()
    #return args.accumulate(args.integers)
    #return args

def bray_curtis_distance(table, sample1_id, sample2_id):
    """
        I took this method from an online tutorial.
        http://readiab.org/book/0.1.3/3/1#4
    """
    numerator = 0
    denominator = 0
    sample1_counts = table[sample1_id]
    sample2_counts = table[sample2_id]
    for sample1_count, sample2_count in zip(sample1_counts, sample2_counts):
        numerator += abs(sample1_count - sample2_count)
        denominator += sample1_count + sample2_count
    return (numerator / denominator)

#newick_tree1 = StringIO('(((((OTU1:0.5, OTU2:0.5):0.5, OTU3:1.0):1.0), (OTU4:0.75, OTU5:0.75):1.25))root;')
#tree1 = TreeNode.read(newick_tree1)

def unweighted_UniFrac(tree, table, sample_id1, sample_id2, verbose = False):
        """
            I took this method from the same tutorial as above.
            http://readiab.org/book/0.1.3/3/1#4
            I have not been able to incorporate this into my Project just yet.
        """
    observed_nodes1 = get_observed_nodes(tree, table, sample_id1, verbose = verbose)
    observed_nodes2 = get_observed_nodes(tree, table, sample_id2, verbose = verbose)
    observed_branch_length = sum(o.length for o in observed_nodes1 | observed_nodes2)
    shared_branch_length = sum(o.length for o in observed_nodes1 & observed_nodes2)
    unique_branch_length = observed_branch_length - shared_branch_length
    unweighted_unifrac = unique_branch_length / observed_length
    return unweighted_unifrac

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d import Axes3D

def scatter_3d(ord_results, df, column, color_map, title = '', axis1 = 0, \
                axis2 = 1, axis3 = 2):
    """
        This was also taken from the tutorial.
        http://readiab.org/book/0.1.3/3/1#4
        I plan on making my Project output a 3D graph.
    """
    coord_matrix = ord_results.site.T
    ids = ord_results.ste_ids
    colors = [color_map[df[column][id_]] for id_ in ord_results.site_ids]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    xs = coord_matrix[axis1]
    ys = coord_matrix[axis2]
    zs = coord_matrix[axis3]
    plot = ax.scatter(xs, ys, zs, c = colors)

    ax.set_xlabel('PC %d' % (axis1 + 1))
    ax.set_ylabel('PC %d' % (axis2 + 1))
    ax.set_zlabel('PC %d' % (axis3 + 1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.set_title(title)
    return fig

def main():
    #otu1_file = read_csv_file(r'scenario1_otus - scenario1_otus.csv')
    #otu1_file = read_csv_file('scenario1_otus.csv')
    otu1_file = read_csv_file('scenario1_otus_with_taxonomy.csv')
    print(otu1_file)

    """set_otu1 = set_test(otu1_file)
    print(set_otu1)"""

    bray_curtis_test = bray_curtis_distance(otu1_file, 's1_p1', 's1_p2')
    print(bray_curtis_test)

    #arg_test = argparse_practice()
    #print(arg_test)



if __name__ == '__main__':
    main()
