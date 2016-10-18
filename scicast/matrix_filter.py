import pickle as pickle
import numpy as np
import pandas as pd
import os
import sys
from subprocess import call
import matplotlib
#matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
import scipy
import json
from sklearn.decomposition import PCA as skPCA
from scipy.spatial import distance
from scipy.cluster import hierarchy
import seaborn as sns
from matplotlib.colors import rgb2hex, colorConverter
from pprint import pprint
import difflib
from operator import itemgetter
import itertools
from functools import reduce
import matplotlib.ticker as ticker
import math
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import matplotlib.patches as patches


def _index_to_label(index):
    """Convert a pandas index or multiindex to an axis label."""
    if isinstance(index, pd.MultiIndex):
        return "-".join(map(str, index.names))
    else:
        return index.name

class matrix_filter(object):

    def __init__(self, data, gene_list=None, cell_list=None, gene_exclude=None, cell_exclude=None, num_expressed=3):
        if isinstance(data, pd.DataFrame):
            matrix_data = data.values
        else:
            matrix_data = np.asarray(data)
            data = pd.DataFrame(matrix_data)

            cell_names = _index_to_label(data.columns)
            gene_names = _index_to_label(data.index)

            self.cell_names = cell_names
            self.gene_names = gene_names
            self.data = data
            self.matrix_data = matrix_data

    def make_new_matrix_gene(org_matrix_by_gene, gene_list_source, exclude_list=""):
        if isinstance(gene_list_source,str):
            gene_df = pd.read_csv(gene_list_source, sep=None)
            try:
                gene_list = list(set(gene_df['GeneID'].tolist()))
                if exclude_list != "":
                    gene_list = [g for g in gene_list if g not in exclude_list]
            except KeyError:
                sys.exit("Error: Please provide Gene list file with 'GeneID' as header.")
            try:
                group_list = gene_df['GroupID'].tolist()
            except KeyError:
                sys.exit("Error: Please provide Gene list file with 'GroupID' as header.")
            try:
                gmatrix_df = org_matrix_by_gene[gene_list]
            except KeyError as error_gene:
                cause = error_gene.args[0]
                absent_gene = cause.split('\'')[1]
                print(absent_gene+' not in matrix file.')
                new_list = [x for x in gene_list if x not in [absent_gene]]
                gmatrix_df = org_matrix_by_gene[new_list]
            cmatrix_df = gmatrix_df.transpose()
            cell_list1 = cmatrix_df.columns.values
            new_cmatrix_df = cmatrix_df[cell_list1]
            new_gmatrix_df = new_cmatrix_df.transpose()
            return new_cmatrix_df, new_gmatrix_df
        elif isinstance(gene_list_source,list):
            if exclude_list == "":
                gene_list = gene_list_source
            else:
                gene_list = [g for g in gene_list_source if g not in exclude_list]
            try:
                gmatrix_df = org_matrix_by_gene[gene_list]
            except KeyError as error_gene:
                cause = error_gene.args[0]
                absent_gene = cause.split('\'')[1]
                print(absent_gene+' not in matrix file.')
                new_list = [x for x in gene_list if x not in [absent_gene]]
                gmatrix_df = org_matrix_by_gene[new_list]
            cmatrix_df = gmatrix_df.transpose()
            cell_list1 = cmatrix_df.columns.values
            new_cmatrix_df = cmatrix_df[cell_list1]
            new_gmatrix_df = new_cmatrix_df.transpose()
            return new_cmatrix_df, new_gmatrix_df
        else:
            sys.exit("Error: gene list must be filepath or a list.")

    def make_new_matrix_cell(org_matrix_by_cell, cell_list_file):
        cell_df = pd.read_csv(cell_list_file, sep=None)
        cell_list_new = list(set([cell.strip('\n') for cell in cell_df['SampleID'].tolist()]))
        cell_list_old = org_matrix_by_cell.columns.tolist()
        overlap = [c for c in cell_list_new if c in cell_list_old]
        not_in_matrix = [c for c in cell_list_new if c not in cell_list_old]
        if not_in_matrix != []:
            print('These cells were in the cell list provided by not found in the matrix provided:')
            print(not_in_matrix)
        new_cmatrix_df = org_matrix_by_cell[overlap]
        new_gmatrix_df = new_cmatrix_df.transpose()
        return new_cmatrix_df, new_gmatrix_df


    def preprocess_df(np_by_cell, gen_list, number_expressed=3):
        g_todelete = []
        for g1, gene in enumerate(np_by_cell):
            cells_exp = (gene >= 1.0).sum()
            if cells_exp < number_expressed:
                g_todelete.append(g1)
        g1_todelete = sorted(g_todelete, reverse = True)
        for pos in g1_todelete:
            if type(gen_list[pos]) != float:
                if args.verbose:
                    print('Gene '+gen_list[pos]+' not expressed in '+str(number_expressed)+' cells.')
                pass
            del gen_list[pos]
        n_by_cell = np.delete(np_by_cell, g1_todelete, axis=0)
        return n_by_cell, gen_list

    def find_top_common_genes(log2_df_by_cell, num_common=25):
        top_common_list = []
        count = 0
        done = False
        log2_df_by_gene = log2_df_by_cell.transpose()
        log2_df2_gene = pd.DataFrame(pd.to_numeric(log2_df_by_gene, errors='coerce'))
        log_mean = log2_df2_gene.mean(axis=0).sort_values(ascending=False)
        log2_sorted_gene = log2_df_by_gene.reindex_axis(log2_df_by_gene.mean(axis=0).sort_values(ascending=False).index, axis=1)
        for gene in log2_sorted_gene.columns.tolist():
            if sum(genes < 1 for genes in log2_df_by_gene[gene])<6:
                if count < num_common:
                    count+=1
                    top_common_list.append(gene)
            if count == num_common:
                done = True
                break
        if done:
            return log2_df_by_gene[top_common_list].transpose()
        else:
            return [0]

    def log2_oulierfilter(df_by_cell, plot=False):
        log2_df = np.log2(df_by_cell+1)
        top_log2 = find_top_common_genes(log2_df)
        if all(top_log2) != 0:
            log2_df2= pd.to_numeric(pd.DataFrame(log2_df), errors='coerce')
            log_mean = top_log2.mean(axis=0).sort_values(ascending=False)
            log2_sorted = top_log2.reindex_axis(top_log2.mean(axis=0).sort_values(ascending=False).index, axis=1)
            xticks = []
            keep_col= []
            log2_cutoff = np.average(np.average(log2_sorted))-2*np.average(np.std(log2_sorted))
            for col, m in zip(log2_sorted.columns.tolist(),log2_sorted.mean()):
                if m > log2_cutoff:
                    keep_col.append(col)
                    xticks.append(col+' '+str("%.2f" % m))
            excluded_cells = [x for x in log2_sorted.columns.tolist() if x not in keep_col]
            filtered_df_by_cell = df_by_cell[keep_col]
            filtered_df_by_gene = filtered_df_by_cell.transpose()
            filtered_log2 = np.log2(filtered_df_by_cell[filtered_df_by_cell>0])
            if plot:
                ax = sns.boxplot(data=filtered_log2, whis= .75, notch=True)
                ax = sns.stripplot(x=filtered_log2.columns.values, y=filtered_log2.mean(axis=0), size=4, jitter=True, edgecolor="gray")
                xtickNames = plt.setp(ax, xticklabels=xticks)
                plt.setp(xtickNames, rotation=90, fontsize=9)
                plt.show()
                plt.clf()
                sns.distplot(filtered_log2.mean())
                plt.show()
            log2_expdf_cell = np.log2(filtered_df_by_cell+1)
            log2_expdf_gene = log2_expdf_cell.transpose()
            return log2_expdf_cell, log2_expdf_gene
        else:
            print("no common genes found")
            return log2_df, log2_df.transpose()
