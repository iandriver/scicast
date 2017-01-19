#!/usr/bin/env python

from __future__ import absolute_import, division, print_function
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import pandas as pd
import os
import sys
import collections
import datetime

'''Cluster main function looks for either -gui flag and takes tkinter flags or
expects argparse arguments.Minimum input is a comma,tab, or space delimted
file of cells or experiments as columns and genes as the index.
'''
def make_command_log(args,new_filepath):
    command_str = ''
    now = datetime.datetime.now()
    command_str+= '-f '+args.filepath
    command_str+= ' -n '+args.base_name
    command_str+= ' -method '+args.method
    command_str+= ' -metric '+args.metric
    command_str+= ' -g '+str(args.gene_number)
    command_str+= ' -depth '+str(args.depth_of_clustering)
    if args.cell_list_filename:
        command_str+= ' -cell_group '+args.cell_list_filename
    if args.gene_list_filename:
        command_str+= ' -gene_list '+args.gene_list_filename
    command_str+= ' -z '+str(args.z_direction)
    if args.qgraph_plot:
        command_str+= ' -qgraph_plot '+args.qgraph_plot
    if args.exclude_genes:
        command_str+= ' -exclude_genes '+args.exclude_genes
    command_str+= ' -kmeans_cluster_range '+args.kmeans_cluster_range
    if args.color_cells:
        command_str+= ' -color_cells "'+args.color_cells+'"'
    if args.color_genes:
        command_str+= ' -color_genes "'+args.color_genes+'"'
    if args.genes_corr != '':
        command_str+= ' -genes_corr '+args.genes_corr
    if args.test_clust_stability != 0:
        command_str+= '-stability '+str(args.test_clust_stability)
    if args.annotate_gene_subset:
        command_str+= ' -annotate_gene_subset '+args.annotate_gene_subset
    if args.no_heatmaps:
        command_str+= ' -no_heatmaps'
    if args.no_corr:
        command_str+= ' -no_corr'
    if args.verbose:
        command_str+= ' -v'
    if args.group_sig_test:
        command_str+= ' -group_sig_test'
    if args.all_sig:
        command_str+= ' -all_sig'
    if args.limit_cells:
        command_str+= ' -limit_cells'
    if args.add_ellipse:
        command_str+= ' -add_ellipse'
    if args.annotate_cell_pca:
        command_str+= ' -annotate_cell_pca'
    if args.kmeans_sig_test:
        command_str+= ' -kmeans_sig_test'
    if args.already_log2:
        command_str+= ' -already_log2'
    if args.sig_unique:
        command_str+= ' -sig_unique'
    if args.use_TSNE:
        command_str+= ' -use_TSNE'
    command_str+= ' -image_format '+str(args.image_format)

    new_file = new_filepath+'/scicast_command_log.txt'
    with open(new_file, 'a+') as f:
        f.write(now.strftime("%Y-%m-%d %H:%M")+'\n')
        f.write(command_str+'\n')

def main():
    #check for -gui flag
    try:
        from .scicast_argparse import check_gui_parser
    except (SystemError, ValueError, ImportError):
        from scicast_argparse import check_gui_parser
    try:
        gui = check_gui_parser()
        run_gui = gui.gui_true
    except:
        run_gui = False
    if run_gui:
        try:
            from .tkinter_scicast import Window
        except (SystemError, ValueError, ImportError):
            from tkinter_scicast import Window
        scicast_window = Window()
        scicast_window.mainloop()
        try:
            from .sci_load import Sci_load
        except:
            from sci_load import Sci_load
        scil = Sci_load()
        try:
            args = scil.load_options(all_options_dict = scicast_window.all_dict)
        except AttributeError:
            sys.exit('Please provide (at a minimum) a valid path to a file and click Run scicast.')

    #if -gui flag isn't found check for all other flags
    else:
        try:
            from .scicast_argparse import get_parser
        except (SystemError, ValueError, ImportError):
            from scicast_argparse import get_parser
        args = get_parser()

    #try to make file in the directory of the file provided
    try:
        if args.base_name:
            new_file = os.path.join(os.path.dirname(args.filepath),args.base_name+'_scicast_analysis')
        #if no base name is provided just make the file as 'scicast_analysis'
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

    make_command_log(args,new_file)

    #check for gene list file if full path is provided and in the same directory as the main file
    if args.gene_list_filename:
        if os.path.isfile(args.gene_list_filename):
            gene_file = args.gene_list_filename
        elif os.path.isfile(os.path.join(os.path.dirname(args.filepath),args.gene_list_filename)):
            gene_file = os.path.join(os.path.dirname(args.filepath),args.gene_list_filename)
        else:
            sys.exit('Error: Cannot find gene list file. Please place the gene list file in the same directory or provide a full path.')
    else:
        gene_file = False

    #check for cell list file if full path is provided and in the same directory as the main file
    if args.cell_list_filename:
        if os.path.isfile(args.cell_list_filename):
            cell_file = args.cell_list_filename
        elif os.path.isfile(os.path.join(os.path.dirname(args.filepath),args.cell_list_filename)):
            cell_file = os.path.join(os.path.dirname(args.filepath),args.cell_list_filename)
        else:
            sys.exit('Error: Cannot find cell list file. Please place the gene list file in the same directory or provide a full path.')
    else:
        cell_file=False

    #make main file in to pandas datframe
    try:
        by_cell = pd.read_csv(args.filepath, sep="\t", index_col=0)
    except OSError:
        sys.exit('Please provide a valid path to a file.')

    #check for duplicate gene names and if they exist drop them from the dataframe
    dup_gene_list = [item for item, count in collections.Counter(by_cell.index).items() if count > 1]
    if len(dup_gene_list) >0:
        by_cell.drop(dup_gene_list, inplace=True)
    #create list of genes
    gene_list_inital = by_cell.index.tolist()
    #create cell list
    cell_list = [x for x in list(by_cell.columns.values)]
    df_by_cell1 = pd.DataFrame(by_cell, columns=cell_list, index=gene_list_inital)

    #process matrix data into matrix class object matrix_data
    try:
        from .matrix_filter import Matrix_filter
    except (SystemError, ValueError, ImportError):
        from matrix_filter import Matrix_filter
    matrix_data = Matrix_filter(cell_df=df_by_cell1 , args=args, cell_list_filepath=cell_file, gene_list_filepath=gene_file)



    #if stability iteration number is provided run heatmap clustering stability function
    stability_iteration_num = int(args.test_clust_stability)
    if  stability_iteration_num != 0:
        try:
            from .stability_test import clust_stability
        except (SystemError, ValueError, ImportError):
            from stability_test import clust_stability
        clust_stability(args, matrix_data, stability_iteration_num)


    #parse kmeans cluster range if seperate by a dash or a comma
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
        try:
            from .dim_reduction import plot_kmeans,plot_SVD,plot_TSNE,return_top_pca_gene
            from .heatmaps import clust_heatmap,make_tree_json,make_subclusters
        except (SystemError, ValueError, ImportError):
            from dim_reduction import plot_kmeans,plot_SVD,plot_TSNE,return_top_pca_gene
            from heatmaps import clust_heatmap,make_tree_json,make_subclusters
        plot_SVD(args, matrix_data, title='all_cells')
        plot_TSNE(args, matrix_data, title='all_cells')

        #if kmeans range is provided run kmeans cluster detection
        if kmeans_range:
            plot_kmeans(args, matrix_data, kmeans_range=kmeans_range, title='all_cells')
        #
        if isinstance(matrix_data.log2_df_cell_gene_restricted, pd.DataFrame):
            top_pca_gene_df, top_pca = return_top_pca_gene(args, matrix_data.log2_df_cell_gene_restricted)
        else:
            top_pca_gene_df, top_pca = return_top_pca_gene(args, matrix_data.log2_df_cell)

        cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, matrix_data, matrix_subset=top_pca_gene_df)
        cc = make_tree_json(cell_linkage, plotted_df_by_gene, new_file)
        make_subclusters(args, cc, matrix_data)

    if not args.no_corr:
        try:
            from .correlation import corr_plot
            from .dim_reduction import return_top_pca_gene
        except (SystemError, ValueError, ImportError):
            from correlation import corr_plot
            from dim_reduction import return_top_pca_gene
        top_pca_gene_df, top_pca = return_top_pca_gene(args, matrix_data.log2_df_cell)

        top_genes_search = top_pca[0:50]
        corr_plot(top_genes_search, top_pca_gene_df, args, matrix_data, title = 'All_cells')

    if args.qgraph_plot == 'both':
        try:
            from .R_qgraph import run_qgraph
        except (SystemError, ValueError, ImportError):
            from R_qgraph import run_qgraph
        run_qgraph(args, matrix_data, gene_or_cell='gene')
        run_qgraph(args, matrix_data, gene_or_cell='cell')
    elif args.qgraph_plot == 'gene':
        try:
            from .R_qgraph import run_qgraph
        except (SystemError, ValueError, ImportError):
            from R_qgraph import run_qgraph
        run_qgraph(args, matrix_data, gene_or_cell='gene')
    elif args.qgraph_plot == 'cell':
        try:
            from .R_qgraph import run_qgraph
        except (SystemError, ValueError, ImportError):
            from R_qgraph import run_qgraph
        run_qgraph(args, matrix_data, gene_or_cell='cell')
    if args.all_sig:
        try:
            from .heatmaps import find_twobytwo
            from .heatmaps import make_tree_json, clust_heatmap
        except:
            from heatmaps import find_twobytwo
            from heatmaps import make_tree_json, clust_heatmap
        cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, matrix_data, plot=False)
        cc = make_tree_json(cell_linkage, plotted_df_by_gene, new_file)
        find_twobytwo(cc, args, matrix_data)
    if args.group_sig_test and args.cell_list_filename:
        try:
            from .significance_testing import multi_group_sig
        except (SystemError, ValueError, ImportError):
            from significance_testing import multi_group_sig
        multi_group_sig(args, matrix_data)


if __name__ == "__main__":
    main()
