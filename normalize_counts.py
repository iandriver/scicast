from rpy2.robjects.packages import STAP


with open('normalize_functions.R', 'r') as f:
    r_file = f.read()
norm_functions = STAP(r_file, "normalize_functions")

#path to fpkm file (usually cuffnorm output)
path_to_file = '/Volumes/Seq_data/count-picard_Pdgfra_ctrl'
#default file name will use genes.fpkm_table from cuffnorm
file_name = 'DESeq__pdgfra_ctrl_all_matrix_norm.txt'
#provide base name for output files
base_name ='Pdgfra_ctrl_counts'
norm_functions.run_deseq2(file_path, save_name, group1_terms, group1_name, group2_terms, group2_name)

def is_valid_file(parser, arg):
    """
    Check if arg is a valid file that already exists on the file system.

    Parameters
    ----------
    parser : argparse object
    arg : str

    Returns
    -------
    arg
    """
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


def get_parser():
    """Get parser object for script cluster.py."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--filepath",
                        dest="filepath",
                        type=lambda x: is_valid_file(parser, x),
                        required =True,
                        help="Filepath to cell-gene matrix file")
    parser.add_argument("-n","--name",
                        dest="base_name",
                        default='',
                        type=str,
                        help="The base name for files that will be generated")
    parser.add_argument("-g","--gene_number",
                        dest="gene_number",
                        default=200,
                        type=int,
                        help="The number of genes that will be selected from top pca genes on each iteration.")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        help="Print more")
    parser.add_argument("-cell_group",
                        dest="cell_list_filename",
                        default=False,
                        help="Provide a filename with two columns: 'SampleID' and 'GroupID'. If GroupID is empty the SampleID list will be used to restrict cells prior to analysis.")
    parser.add_argument("-gene_list",
                        dest="gene_list_filename",
                        default=False,
                        help="Path to a file with a two columns with headers 'GeneID' and 'GroupID' (Singular format). GeneID list will be used to create a new matrix file with only those genes included.")
    parser.add_argument("-group_sig_test",
                        action="store_true",
                        dest="group_sig_test",
                        help="If cell groups are provided, perform significance testing between all groups (independent of any cluster groups).")
    parser.add_argument("-stability",
                        dest="test_clust_stability",
                        default=0,
                        help="Provide a number of iterations to test how stable clustering is as the number of top PCA genes changes from 100-1000. Output will be clustering heatmaps for each iteration a summary of changes as gene number varies.")
    parser.add_argument("--metric",
                        dest="metric",
                        default='euclidean',
                        help="The distance metric to use. The distance function can be: braycurtis, canberra, chebyshev, cityblock, correlation, cosine, dice, euclidean, hamming, jaccard, kulsinski, mahalanobis, matching, minkowski, rogerstanimoto, russellrao, seuclidean, sokalmichener, sokalsneath, sqeuclidean, yule.")
    parser.add_argument("--method",
                        dest="method",
                        default='average',
                        help="The linkage method to use. The linkage algorithm can be: single, complete, average, weighted, centroid median or ward.")
    parser.add_argument("--depth",
                        dest="depth_of_clustering",
                        type=int,
                        default=20,
                        help="The size in cell number at which sub-clustering analysis will stop clustering, pca and correlation analysis.")
    parser.add_argument("--genes_corr",
                        dest="genes_corr",
                        default='',
                        help="If you want to look for correation on a specific gene or set of genes enter them as a comma seperated list i.e. 'Gapdh,Actb'.")
    parser.add_argument("--all_sig",
                        dest="all_sig",
                        action="store_true",
                        help="If you want to look for correation on a specific gene or set of genes enter them as a comma seperated list i.e. 'Gapdh,Actb'.")
    parser.add_argument("--z",
                        dest="z_direction",
                        default=0,
                        help="Either 0 (rows) or 1 (columns) or None. Whether or not to calculate z-scores for the rows or the columns. Z scores are: z = (x - mean)/std, so values in each row (column) will get the mean of the row (column) subtracted, then divided by the standard deviation of the row (column). This ensures that each row (column) has mean of 0 and variance of 1.")
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)
