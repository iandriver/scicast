
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter,SUPPRESS
import sys as _sys
from gettext import gettext as _
import os

def is_valid_file(parser, arg):
    import os
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


class MyParser(ArgumentParser):
    def error(self, message):
        usage = self.usage
        self.usage = None
        self.print_usage(_sys.stderr)
        self.exit(2, _('%s: error: %s\n') % (self.prog, message))
        self.usage = usage
def check_gui_parser():
    p = MyParser(description=__doc__,formatter_class=ArgumentDefaultsHelpFormatter, usage=SUPPRESS)
    p.add_argument("-gui",
                        dest="gui_true",
                        default = False,
                        action = "store_true",
                        help="To run gui instead of command line.")
    gui_args = p.parse_args()
    return(gui_args)

def get_parser():
    """Get parser object for script cluster.py."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "-filepath",
                        dest="filepath",
                        type=lambda x: is_valid_file(parser, x),
                        required =True,
                        help="Filepath to cell-gene matrix file")
    parser.add_argument("-n","-base_name",
                        dest="base_name",
                        default='scicast_analysis',
                        type=str,
                        help="The base name for files that will be generated.")
    parser.add_argument("-g","-gene_number",
                        dest="gene_number",
                        default=200,
                        type=int,
                        help="The number of genes that will be selected from top pca genes on each iteration.")
    parser.add_argument("-v", "-verbose",
                        action="store_true",
                        dest="verbose",
                        default=False,
                        help="Print more")
    parser.add_argument("-cell_group","-cell_list_filename",
                        dest="cell_list_filename",
                        default=False,
                        help="Provide a filename with two columns: 'SampleID' and 'GroupID'. If GroupID is empty the SampleID list will be used to restrict cells prior to analysis.")
    parser.add_argument("-gene_list",
                        dest="gene_list_filename",
                        default=False,
                        type=str,
                        help="Path to a file with a two columns with headers 'GeneID' and 'GroupID' (Singular format). GeneID list will be used to create a new matrix file with only those genes included.")
    parser.add_argument("-group_sig_test",
                        action="store_true",
                        default=False,
                        dest="group_sig_test",
                        help="If cell groups are provided, perform significance testing between all groups (independent of any cluster groups).")
    parser.add_argument("-stability","-test_clust_stability",
                        dest="test_clust_stability",
                        type=int,
                        default=0,
                        help="Provide a number of iterations to test how stable clustering is as the number of top PCA genes changes from 100-1000. Output will be clustering heatmaps for each iteration a summary of changes as gene number varies.")
    parser.add_argument("-metric",
                        dest="metric",
                        default='euclidean',
                        type=str,
                        help="The distance metric to use. The distance function can be: braycurtis, canberra, chebyshev, cityblock, correlation, cosine, dice, euclidean, hamming, jaccard, kulsinski, mahalanobis, matching, minkowski, rogerstanimoto, russellrao, seuclidean, sokalmichener, sokalsneath, sqeuclidean, yule.")
    parser.add_argument("-method",
                        dest="method",
                        default='weighted',
                        type=str,
                        help="The linkage method to use. The linkage algorithm can be: single, complete, average, weighted, centroid, median or ward.")
    parser.add_argument("-depth","-depth_of_clustering",
                        dest="depth_of_clustering",
                        type=int,
                        default=20,
                        help="The size in cell number at which sub-clustering analysis will stop clustering, pca and correlation analysis.")
    parser.add_argument("-genes_corr",
                        dest="genes_corr",
                        default='',
                        type=str,
                        help="If you want to look for correlation on a specific gene or set of genes enter them as a comma seperated list i.e. 'Gapdh,Actb'.")
    parser.add_argument("-all_sig",
                        dest="all_sig",
                        default=False,
                        action="store_true",
                        help="Do significance testing on all hierarchical clustering groups. The minimum number of cells in a group is set by --depth.")
    parser.add_argument("-already_log2",
                        dest="already_log2",
                        default=False,
                        action="store_true",
                        help="If the gene matrix is already log2 transformed i.e. rld or vsd DESeq2, this will disable log2 transform.")
    parser.add_argument("-z","-z_direction",
                        dest="z_direction",
                        default=0,
                        help="Either 0 (rows) or 1 (columns) or None. Whether or not to calculate z-scores for the rows or the columns. Z scores are: z = (x - mean)/std, so values in each row (column) will get the mean of the row (column) subtracted, then divided by the standard deviation of the row (column). This ensures that each row (column) has mean of 0 and variance of 1.")
    parser.add_argument("-no_corr",
                        dest="no_corr",
                        action="store_true",
                        default=False,
                        help="Don't run correlation search. Default is on.")
    parser.add_argument("-no_heatmaps",
                        dest="no_heatmaps",
                        action="store_true",
                        default=False,
                        help="Don't run heatmaps and pca. Default is will generate heatmaps and pca.")
    parser.add_argument("-exclude_genes",
                        dest="exclude_genes",
                        default=False,
                        type=str,
                        help="Filepath to list of genes to be excluded from analysis (at all levels of analysis). Header must be 'GeneID'.")
    parser.add_argument("-limit_cells",
                        dest="limit_cells",
                        action="store_true",
                        default=False,
                        help="With cell group file, will exclude all other cells from analysis (at all levels of analysis). Header must be 'SampleID'.")
    parser.add_argument("-add_ellipse",
                        dest="add_ellipse",
                        action="store_true",
                        default=False,
                        help="When present colored ellipses will be added to cell and gene PCA plots. Must provide gene and/or cell groups.")
    parser.add_argument("-annotate_gene_subset",
                        dest="annotate_gene_subset",
                        default=False,
                        help="Provide path or filename (if in same file) to file with genes to be annotated on gene PCA. Must have 'GeneID' header.")
    parser.add_argument("-annotate_cell_pca",
                        dest="annotate_cell_pca",
                        action="store_true",
                        default=False,
                        help="Option will annotate cell PCA with cell names. Default is off (False).")
    parser.add_argument("-color_cells",
                        dest="color_cells",
                        default=False,
                        type=str,
                        help="Follow with cell group names and pairs (matching cell group names). Forces color and marker style on group. Just color or color and marker can be provided.  Example: Group1,color,marker Group2,red,v")
    parser.add_argument("-color_genes",
                        dest="color_genes",
                        default=False,
                        type=str,
                        help="Can be 'same' if names are the same as cell groups. Follow with gene group names and pairs (matching gene group names). Forces color and marker style on group. Just color or color and marker can be provided.  Example: Group1,color,marker Group2,red,v")
    parser.add_argument("-qgraph_plot",
                        dest="qgraph_plot",
                        default=False,
                        type=str,
                        help="Can be 'cell' 'gene' or 'both' will plot qgraph gene and/or cell network and pca correlation network plot to pdf. Requires both gene and cell groups to be provided.")
    parser.add_argument("-sig_unique",
                        dest="sig_unique",
                        action="store_true",
                        default=False,
                        help="group_sig_test will by default find unique sets of significant genes, so that the whole list has no duplicates. This switches the top significant genes amoung groups to allow repeats.")
    parser.add_argument("-kmeans_cluster_range",
                        dest="kmeans_cluster_range",
                        type=str,
                        default='2,4',
                        help="range of cluster sizes to run kmeans clustering (inclusive)")
    parser.add_argument("-kmeans_sig_test",
                        dest="kmeans_sig_test",
                        default=False,
                        action='store_true',
                        help="Run significance test between groups defined by kmeans clusters.")
    parser.add_argument("-use_TSNE",
                        dest="use_TSNE",
                        default=False,
                        action='store_true',
                        help="Use t-SNE dimensions for k-means clustering rather than SVD.")
    parser.add_argument("-image_format",
                        dest="image_format",
                        default="pdf",
                        help="Specify image format for output (pdf,png,tiff,etc.)")


    args = parser.parse_args()

    yes_or_no_options = ["Don't Run Heatmaps", "Don't Run Correlation", "Verbose", "Test Significance by Groups (User Defined)", "Test Significance by Unbiased Clusters", "Exclude Cells Not in User Cell Groups", "Add Ellipse", "Add Cell Names to PCA", "Display Only Unique Signifcant Genes", "Run Significance Test for kmeans clusters", "Input Matrix is already log2", "use t-SNE (for kmeans clustering)"]
    yes_or_no_answers = [args.no_heatmaps,args.no_corr, args.verbose,args.group_sig_test, args.all_sig, args.limit_cells, args.add_ellipse, args.annotate_cell_pca, args.sig_unique, args.kmeans_sig_test, args.already_log2, args.use_TSNE]

    all_options_dict = {}
    all_options_dict['filepath'] = args.filepath
    all_options_dict['output_name'] = args.base_name
    all_options_dict['method'] = args.method
    all_options_dict['metric'] =args.metric
    all_options_dict['gene_number'] =args.gene_number
    all_options_dict['depth'] = args.depth_of_clustering
    all_options_dict['cell_file'] = args.cell_list_filename
    all_options_dict['gene_file'] = args.gene_list_filename
    all_options_dict['zdir'] = args.z_direction
    all_options_dict['qgraph'] = args.qgraph_plot
    all_options_dict['kmeans_cluster_range'] = args.kmeans_cluster_range
    all_options_dict['exclude_genes'] = args.exclude_genes
    all_options_dict['already_log2'] = args.already_log2
    all_options_dict['color_cells'] = args.color_cells
    all_options_dict['color_genes'] = args.color_genes
    all_options_dict['test_clust_stability'] = args.test_clust_stability
    all_options_dict['genes_corr'] = args.genes_corr
    all_options_dict['annotate_gene_subset'] = args.annotate_gene_subset
    all_options_dict['image_format'] = args.image_format
    for var, flag in zip(yes_or_no_answers, yes_or_no_options):
        all_options_dict[flag] = var
    try:
        from .sci_load import Sci_load
    except (SystemError, ValueError, ImportError):
        from sci_load import Sci_load
    scil = Sci_load()
    opts_all = scil.load_options(all_options_dict = all_options_dict)

    return opts_all

if __name__ == "__main__":
    pass
