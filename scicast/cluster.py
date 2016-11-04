#!/usr/bin/env python

from __future__ import absolute_import, division, print_function
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import pandas as pd
import os
import sys
import collections



def main():
    from .scicast_argparse import check_gui_parser
    try:
        gui = check_gui_parser()
        run_gui = gui.gui_true
    except:
        run_gui = False
    if run_gui:
        from .tkinter_scicast import Window
        scicast_window = Window()
        scicast_window.mainloop()
        from .sci_load import Sci_load
        scil = Sci_load()
        try:
            args = scil.load_options(all_options_dict = scicast_window.all_dict)
        except AttributeError:
            sys.exit('Please provide (at a minimum) a valid path to a file and click Run scicast.')


    else:
        from .scicast_argparse import get_parser
        args = get_parser()



    try:
        if args.base_name:
            new_file = os.path.join(os.path.dirname(args.filepath),args.base_name+'_scicast_analysis')
        else:
            new_file = os.path.join(os.path.dirname(args.filepath),'scicast_analysis')
    except AttributeError:
        sys.exit('Please provide a valid path to a file.')
    if args.verbose:
        print('Making new folder for results of SCICAST clustering: '+new_file)
    try:
        os.mkdir(new_file)
    except OSError:
        if args.verbose:
            print(new_file+' already exists. Files will be overwritten.')

    if args.gene_list_filename:
        if os.path.isfile(args.gene_list_filename):
            gene_file = args.gene_list_filename
        elif os.path.isfile(os.path.join(os.path.dirname(args.filepath),args.gene_list_filename)):
            gene_file = os.path.join(os.path.dirname(args.filepath),args.gene_list_filename)
        else:
            sys.exit('Error: Cannot find gene list file. Please place the gene list file in the same directory or provide a full path.')
    else:
        gene_file = False


    if args.cell_list_filename:
        if os.path.isfile(args.cell_list_filename):
            cell_file = args.cell_list_filename
        elif os.path.isfile(os.path.join(os.path.dirname(args.filepath),args.cell_list_filename)):
            cell_file = os.path.join(os.path.dirname(args.filepath),args.cell_list_filename)
        else:
            sys.exit('Error: Cannot find cell list file. Please place the gene list file in the same directory or provide a full path.')
    else:
        cell_file=False

    try:
        by_cell = pd.read_csv(args.filepath, sep="\t", index_col=0)
    except OSError:
        sys.exit('Please provide a valid path to a file.')

    dup_gene_list = [item for item, count in collections.Counter(by_cell.index).items() if count > 1]
    if len(dup_gene_list) >0:
        by_cell.drop(dup_gene_list, inplace=True)
    #create list of genes
    gene_list_inital = by_cell.index.tolist()
    #create cell list
    cell_list = [x for x in list(by_cell.columns.values)]
    df_by_cell1 = pd.DataFrame(by_cell, columns=cell_list, index=gene_list_inital)

    from .matrix_filter import Matrix_filter
    matrix_data = Matrix_filter(cell_df=df_by_cell1 , args=args, cell_list_filepath=cell_file, gene_list_filepath=gene_file)



    #run heatmap clustering stability function
    stability_iteration_num = int(args.test_clust_stability)
    if  stability_iteration_num != 0:
        from .stability_test import clust_stability
        stability_ratio = clust_stability(args, matrix_data, stability_iteration_num)



    if ',' not in args.kmeans_cluster_range and '-' not in args.kmeans_cluster_range and args.kmeans_cluster_range == '0':
        kmeans_range = []
    elif ',' in args.kmeans_cluster_range:
        kmeans_range = [int(k) for k in args.kmeans_cluster_range.split(',')]
    elif '-' in args.kmeans_cluster_range:
        kmeans_range = [int(k) for k in args.kmeans_cluster_range.split('-')]
    elif len(args.kmeans_cluster_range) == 1:
        kmeans_range = [int(args.kmeans_cluster_range), int(args.kmeans_cluster_range)]

    #run heatmaps and PCA only if no_heatmaps flag is not provided
    if not args.no_heatmaps:
        from .dim_reduction import plot_kmeans,plot_SVD,plot_TSNE,return_top_pca_gene
        from .heatmaps import clust_heatmap,make_tree_json,make_subclusters
        plot_SVD(args, matrix_data, title='all_cells')
        plot_TSNE(args, matrix_data, title='all_cells')

        if kmeans_range:
            plot_kmeans(args, matrix_data, kmeans_range=kmeans_range, title='all_cells')
        if isinstance(matrix_data.log2_df_cell_gene_restricted, pd.DataFrame):
            top_pca_gene_df, top_pca = return_top_pca_gene(args, matrix_data.log2_df_cell_gene_restricted)
        else:
            top_pca_gene_df, top_pca = return_top_pca_gene(args, matrix_data.log2_df_cell)

        cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, matrix_data)
        cc = make_tree_json(cell_linkage, plotted_df_by_gene, new_file)
        make_subclusters(args, cc, matrix_data)

    if not args.no_corr:
        from .correlation import corr_plot
        from .dim_reduction import return_top_pca_gene
        top_pca_gene_df, top_pca = return_top_pca_gene(args, matrix_data.log2_df_cell)

        top_genes_search = top_pca[0:50]
        corr_plot(top_genes_search, top_pca_gene_df, args, matrix_data, title = 'All_cells')

    if args.qgraph_plot == 'both':
        from .R_qgraph import run_qgraph
        run_qgraph(args, matrix_data, gene_or_cell='gene')
        run_qgraph(args, matrix_data, gene_or_cell='cell')
    elif args.qgraph_plot == 'gene':
        from .R_qgraph import run_qgraph
        run_qgraph(args, matrix_data, gene_or_cell='gene')
    elif args.qgraph_plot == 'cell':
        from .R_qgraph import run_qgraph
        run_qgraph(args, matrix_data, gene_or_cell='cell')
    if args.all_sig:
        from .heatmaps import find_twobytwo
        from .heatmaps import make_tree_json, clust_heatmap
        cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, matrix_data, plot=False)
        cc = make_tree_json(cell_linkage, plotted_df_by_gene, new_file)
        find_twobytwo(cc, args, matrix_data)
    if args.group_sig_test and args.cell_list_filename:
        from .significance_testing import multi_group_sig
        multi_group_sig(args, matrix_data)


if __name__ == "__main__":
    main()
