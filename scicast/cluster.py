import numpy as np
import pandas as pd
import os
import sys
from subprocess import call
import matplotlib
matplotlib.use('TkAgg')
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
import matplotlib.patches as patches
from collections import defaultdict
import collections
from sklearn.manifold import TSNE
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
import matplotlib.mlab as mlab





def make_new_matrix_gene(org_matrix_by_gene, gene_list_source, exclude_list=""):
    if isinstance(gene_list_source,str):
        gene_df = pd.read_csv(open(gene_list_source,'rU'), sep=None, engine='python')
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
            cause1 = error_gene.args[0].strip(' not in index')
            cause = [v.strip('\n\' ') for v in cause1.strip('[]').split(' ')]
            absent_gene = cause
            print(' '.join(absent_gene)+' not in matrix file.')
            new_list = [x for x in gene_list if x not in absent_gene]
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
    cell_df = pd.read_csv(open(cell_list_file,'rU'), sep=None, engine='python')
    cell_list_new = list(set([cell.strip('\n') for cell in cell_df['SampleID'].tolist()]))
    cell_list_old = org_matrix_by_cell.columns.tolist()
    overlap = [c for c in cell_list_new if c in cell_list_old]
    not_in_matrix = [c for c in cell_list_new if c not in cell_list_old]
    if not_in_matrix != []:
        print('These cells were in the cell list provided, but not found in the matrix provided:')
        print(not_in_matrix)
    new_cmatrix_df = org_matrix_by_cell[overlap]
    new_gmatrix_df = new_cmatrix_df.transpose()
    return new_cmatrix_df, new_gmatrix_df


def threshold_genes(by_gene, number_expressed=1, gene_express_cutoff=1.0):
    by_gene.apply(lambda column: (column >= 1).sum())
    return


def find_top_common_genes(log2_df_by_cell, num_common=25):
    top_common_list = []
    count = 0
    done = False
    log2_df_by_gene = log2_df_by_cell.transpose()
    log2_df2_gene = log2_df_by_gene.apply(pd.to_numeric,errors='coerce')
    log_mean = log2_df2_gene.mean(axis=0).sort_values(ascending=False)
    try:
        log2_sorted_gene = log2_df_by_gene.reindex_axis(log2_df_by_gene.mean(axis=0).sort_values(ascending=False).index, axis=1)
    except ValueError:
        overlap_list = [item for item, count in collections.Counter(log2_df_by_cell.index).items() if count > 1]
        print(overlap_list, len(overlap_list))
        sys.exit('Error: Duplicate GeneIDs are present.')
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

def log2_oulierfilter(df_by_cell, plot=False, already_log2=False):
    if not already_log2:
        log2_df = np.log2(df_by_cell+1)

    else:
        log2_df = df_by_cell
    top_log2 = find_top_common_genes(log2_df)
    if all(top_log2) != 0:
        log2_df2= log2_df.apply(pd.to_numeric,errors='coerce')
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
        if not already_log2:
            filtered_log2 = np.log2(filtered_df_by_cell[filtered_df_by_cell>0])
        else:
            filtered_log2 = filtered_df_by_cell[filtered_df_by_cell>0]
        if plot:
            ax = sns.boxplot(data=filtered_log2, whis= .75, notch=True)
            ax = sns.stripplot(x=filtered_log2.columns.values, y=filtered_log2.mean(axis=0), size=4, jitter=True, edgecolor="gray")
            xtickNames = plt.setp(ax, xticklabels=xticks)
            plt.setp(xtickNames, rotation=90, fontsize=9)
            plt.show()
            plt.clf()
            sns.distplot(filtered_log2.mean())
            plt.show()
        if not already_log2:
            log2_expdf_cell = np.log2(filtered_df_by_cell+1)
        else:
            log2_expdf_cell = filtered_df_by_cell
        log2_expdf_gene = log2_expdf_cell.transpose()
        return log2_expdf_cell, log2_expdf_gene
    else:
        print("no common genes found")
        return log2_df, log2_df.transpose()


def augmented_dendrogram(*args, **kwargs):
    plt.clf()
    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        for i, d in zip(ddata['icoord'], ddata['dcoord'], ):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y >= 200000:
                plt.plot(x, y, 'ro')
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -8),
                         textcoords='offset points',
                         va='top', ha='center')
    plt.show()
    plt.savefig(os.path.join(new_file,'augmented_dendrogram.png'))

def cluster_indices(cluster_assignments):
    n = cluster_assignments.max()
    indices = []
    for cluster_number in range(1, n + 1):
        indices.append(np.where(cluster_assignments == cluster_number)[0])
    return indices

def clust_members(r_link, cutoff):
    clust = fcluster(r_link,cutoff)
    num_clusters = clust.max()
    indices = cluster_indices(clust)
    return  num_clusters, indices

def print_clust_membs(indices, cell_list):
    for k, ind in enumerate(indices):
        print("cluster", k + 1, "is", [cell_list[x] for x in ind])

def plot_tree(dendr, path_filename, pos=None, save=False):
    icoord = scipy.array(dendr['icoord'])
    dcoord = scipy.array(dendr['dcoord'])
    color_list = scipy.array(dendr['color_list'])
    xmin, xmax = icoord.min(), icoord.max()
    ymin, ymax = dcoord.min(), dcoord.max()
    if pos:
        icoord = icoord[pos]
        dcoord = dcoord[pos]
    for xs, ys, color in zip(icoord, dcoord, color_list):
        plt.plot(xs, ys, color)
    plt.xlim(xmin-10, xmax + 0.1*abs(xmax))
    plt.ylim(ymin, ymax + 0.1*abs(ymax))
    if save:
        plt.savefig(os.path.join(path_filename,'plot_dendrogram.png'))
    plt.show()


# Create a nested dictionary from the ClusterNode's returned by SciPy
def add_node(node, parent):
	# First create the new node and append it to its parent's children
	newNode = dict( node_id=node.id, children=[] )
	parent["children"].append( newNode )

	# Recursively add the current node's children
	if node.left: add_node( node.left, newNode )
	if node.right: add_node( node.right, newNode )


cc = []
# Label each node with the names of each leaf in its subtree
def label_tree(n, id2name):
    # If the node is a leaf, then we have its name
    if len(n["children"]) == 0:
        leafNames = [ id2name[n["node_id"]] ]

    # If not, flatten all the leaves in the node's subtree
    else:
        leafNames = reduce(lambda ls, c: ls + label_tree(c,id2name), n["children"], [])

    cc.append((len(leafNames), [x.strip('\n') for x in leafNames]))
    cc.sort(key=lambda tup: tup[0], reverse = True)

    # Delete the node id since we don't need it anymore and
    # it makes for cleaner JSON
    del n["node_id"]

    # Labeling convention: "-"-separated leaf names
    n["name"] = name = "-".join(sorted(map(str, leafNames)))

    return leafNames

#Makes labeled json tree for visulaization in d3, makes and returns cc object within label_tree
def make_tree_json(row_clusters, df_by_gene, path_filename):
    T= hierarchy.to_tree(row_clusters)

    # Create dictionary for labeling nodes by their IDs
    labels = list(df_by_gene.index)
    id2name = dict(zip(range(len(labels)), labels))

    # Initialize nested dictionary for d3, then recursively iterate through tree
    d3Dendro = dict(children=[], name="Root1")
    add_node( T, d3Dendro )
    label_tree( d3Dendro["children"][0], id2name )
    # Output to JSON
    json.dump(d3Dendro, open(os.path.join(path_filename,"d3-dendrogram.json"), "w"), sort_keys=True, indent=4)

    return cc


#finds significant genes between subclusters
def find_twobytwo(cc, df_by_cell, full_by_cell_df, path_filename, cluster_size=20):
    gene_list = full_by_cell_df.index.tolist()
    by_gene_df = full_by_cell_df.transpose()
    pair_dict = {}
    parent = cc[0][1]
    p_num = cc[0][0]
    l_nums = [x[0] for x in cc]
    c_lists = [c[1] for c in cc[1:]]
    unique_count = 1
    pair_list = []
    for i, c in enumerate(c_lists):
        for i2, c2 in enumerate(c_lists):
            overlap = [i for i in c if i in c2]
            if not overlap and len(c)>=cluster_size and len(c2)>=cluster_size:
                if (c,c2) not in pair_list:
                    pair_list.append((c,c2))
                    pair_list.append((c2,c))
                    pair_dict[str(len(c))+'cells_vs_'+str(len(c2))+'cells'+str(unique_count)]= [c, c2]
                    unique_count+=1

    for v, k in pair_dict.items():
        g_pvalue_dict = {}
        index_list = []
        sig_gene_list = []
        cell_list1 = [x.strip('\n') for x in k[0]]
        cell_list2 = [xx.strip('\n') for xx in k[1]]
        group1 = str(len(cell_list1))
        group2 = str(len(cell_list2))
        df_by_cell_1 = full_by_cell_df[cell_list1]
        df_by_cell_2 = full_by_cell_df[cell_list2]
        df_by_gene_1 = df_by_cell_1.transpose()
        df_by_gene_2 = df_by_cell_2.transpose()
        for g in gene_list:
            g_pvalue = scipy.stats.f_oneway(df_by_gene_1[g], df_by_gene_2[g])
            if g_pvalue[0] > 0 and g_pvalue[1] <= 1:
                g_pvalue_dict[g] = g_pvalue
                if g not in [s[0] for s in sig_gene_list]:
                    sig_gene_list.append([g, g_pvalue[1]])

        sig_gene_list.sort(key=lambda tup: tup[1])
        pvalues = [p[1] for p in sig_gene_list]
        gene_index = [ge[0] for ge in sig_gene_list]
        mean_log2_exp_list = []
        sig_1_2_list = []
        mean1_list = []
        mean2_list = []
        for sig_gene in gene_index:
            sig_gene_df = by_gene_df[sig_gene]
            mean_log2_exp_list.append(sig_gene_df.mean())
            sig_cell_df = sig_gene_df.transpose()
            mean_cell1 = sig_cell_df[cell_list1].mean()
            mean1_list.append(mean_cell1)
            mean_cell2 = sig_cell_df[cell_list2].mean()
            mean2_list.append(mean_cell2)
            ratio_1_2 = (mean_cell1+1)/(mean_cell2+1)
            sig_1_2_list.append(ratio_1_2)
        sig_df = pd.DataFrame({'pvalues':pvalues,'mean_all':mean_log2_exp_list,'mean_group1':mean1_list, 'mean_group2':mean2_list, 'ratio_1_2':sig_1_2_list}, index=gene_index)
        cell_names_df = pd.DataFrame({'cells1':pd.Series(cell_list1, index=range(len(cell_list1))), 'cells2':pd.Series(cell_list2, index=range(len(cell_list2)))})
        sig_df.to_csv(os.path.join(path_filename,'sig_'+v+'_pvalues.txt'), sep = '\t')
        cell_names_df.to_csv(os.path.join(path_filename,'sig_'+v+'_cells.txt'), sep = '\t')

def ellip_enclose(points, color, inc=1, lw=2, nst=2):
    """
    Plot the minimum ellipse around a set of points.

    Based on:
    https://github.com/joferkington/oost_paper_code/blob/master/error_ellipse.py
    """

    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    x = points[:,0]
    y = points[:,1]
    cov = np.cov(x, y)
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    w, h = 2 * nst * np.sqrt(vals)
    center = np.mean(points, 0)
    ell = patches.Ellipse(center, width=inc*w, height=inc*h, angle=theta,
                          facecolor=color, alpha=0.2, lw=0)
    edge = patches.Ellipse(center, width=inc*w, height=inc*h, angle=theta,
                          facecolor='none', edgecolor=color, lw=lw)
    return ell, edge

def return_top_pca_gene(df_by_gene, num_genes=100):
    gene_pca = skPCA(n_components=3)
    np_by_gene = np.asarray(df_by_gene)

    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(gene_pca.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=df_by_gene.columns.tolist())
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len(df_by_gene.columns.tolist()))
    top_pca_list = Pc_sort_df.index.tolist()
    new_gene_matrix = df_by_gene[top_pca_list[0:num_genes]]
    return new_gene_matrix, top_pca_list[0:num_genes]

def plot_PCA(args, df_by_gene, path_filename, num_genes=100, gene_list_filter=False, title='', plot=False, label_map=False, gene_map = False):
    gene_list = df_by_gene.columns.tolist()
    sns.set_palette("RdBu_r", 10, 1)
    if gene_list_filter:
        sig_by_gene = df_by_gene[gene_list_filter]
        sig_by_cell = sig_by_gene.transpose()
    else:
        sig_by_gene = df_by_gene
        sig_by_cell = sig_by_gene.transpose()
    gene_pca = skPCA(n_components=3)
    np_by_gene = np.asarray(sig_by_gene)

    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(gene_pca.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=sig_by_gene.columns.tolist())
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len(sig_by_gene.columns.tolist()))
    top_pca_list = Pc_sort_df.index.tolist()

    top_by_gene = df_by_gene[top_pca_list[0:num_genes]]
    gene_top = skPCA(n_components=2)
    cell_pca = skPCA(n_components=2)
    top_by_cell = top_by_gene.transpose()
    np_top_gene = np.asarray(top_by_cell)
    np_top_cell = np.asarray(top_by_gene)
    top_cell_trans = cell_pca.fit_transform(np_top_cell)
    top_gene_trans = gene_top.fit_transform(np_top_gene)
    if not np.isnan(top_cell_trans).any():
        fig, (ax_cell, ax_gene) = plt.subplots(2, 1, figsize=(15, 30), sharex=False)
        rect_cell = ax_cell.patch
        rect_gene = ax_gene.patch
        rect_cell.set_facecolor('white')
        rect_gene.set_facecolor('white')
        ax_cell.grid(b=True, which='major', color='grey', linestyle='--', linewidth=0.3)
        ax_gene.grid(b=True, which='major', color='grey', linestyle='--', linewidth=0.3)
        if label_map:
            annotate = args.annotate_cell_pca
            X = [x for x in top_cell_trans[:, 0]]
            Y = [y for y in top_cell_trans[:, 1]]
            labels = [label_map[cell][2] for cell in top_by_cell.columns.tolist()]
            markers = [label_map[cell][1] for cell in top_by_cell.columns.tolist()]
            colors = [label_map[cell][0] for cell in top_by_cell.columns.tolist()]
            label_done = []
            xy_by_color_dict = {}
            for c in set(colors):
                xy_by_color_dict[c] = []
            for X_pos, Y_pos, m, color, l in zip(X, Y, markers, colors, labels):
                if l in label_done:
                    lab = ''
                else:
                    lab= l
                    label_done.append(l)
                xy_by_color_dict[color].append([X_pos,Y_pos])
                ax_cell.scatter(X_pos, Y_pos, marker=m, c=color, label=lab, s=30)
            if args.add_ellipse:
                for c in set(colors):
                    ell, edge = ellip_enclose(np.asarray(xy_by_color_dict[c]), c)
                    ax_cell.add_artist(ell)
                    ax_cell.add_artist(edge)
        else:
            ax_cell.scatter(top_cell_trans[:, 0], top_cell_trans[:, 1], alpha=0.75)
            annotate = args.annotate_cell_pca
        ax_cell.set_xlim([min(top_cell_trans[:, 0])-1, max(top_cell_trans[:, 0]+1)])
        ax_cell.set_ylim([min(top_cell_trans[:, 1])-1, max(top_cell_trans[:, 1]+2)])
        ax_cell.set_title(title+'_cell')
        if label_map:
            handles, labs = ax_cell.get_legend_handles_labels()
            # sort both labels and handles by labels
            labs, handles = zip(*sorted(zip(labs, handles), key=lambda t: t[0]))
            ax_cell.legend(handles, labs, loc='best', ncol=1, prop={'size':12}, markerscale=1.5, frameon=True)
        ax_cell.set_xlabel('PC1')
        ax_cell.set_ylabel('PC2')
        if annotate:
            for label, x, y in zip(top_by_cell.columns, top_cell_trans[:, 0], top_cell_trans[:, 1]):
                ax_cell.annotate(label, (x+0.1, y+0.1))

        if gene_map:
            X = [x for x in top_gene_trans[:, 0]]
            Y = [y for y in top_gene_trans[:, 1]]
            labels = top_by_gene.columns.tolist()
            markers = [gene_map[gene][1] for gene in top_by_gene.columns.tolist()]
            colors = [gene_map[gene][0] for gene in top_by_gene.columns.tolist()]
            xy_by_color_dict = {}
            for c in set(colors):
                xy_by_color_dict[c] = []
            for X_pos, Y_pos, m, color, l in zip(X, Y, markers, colors, labels):
                xy_by_color_dict[color].append([X_pos,Y_pos])
                ax_gene.scatter(X_pos, Y_pos, marker=m, c=color, label = l, s=30)
            if args.add_ellipse:
                for c in set(colors):
                    ell, edge = ellip_enclose(np.asarray(xy_by_color_dict[c]), c)
                    ax_gene.add_artist(ell)
                    ax_gene.add_artist(edge)

        else:
            ax_gene.scatter(top_gene_trans[:, 0], top_gene_trans[:, 1], alpha=0.75)
        ax_gene.set_xlim([min(top_gene_trans[:, 0])-1, max(top_gene_trans[:, 0])+1])
        ax_gene.set_ylim([min(top_gene_trans[:, 1])-1, max(top_gene_trans[:, 1])+2])
        ax_gene.set_title(title+'_gene')
        ax_gene.set_xlabel('PC1')
        ax_gene.set_ylabel('PC2')
        if args.annotate_gene_subset:
            plot_subset_path = os.path.join(os.path.dirname(args.filepath),args.annotate_gene_subset)
            genes_plot = pd.read_csv(plot_subset_path, sep='\t', index_col=False)
            for label, x, y in zip(top_by_gene.columns, top_gene_trans[:, 0], top_gene_trans[:, 1]):
                if label in genes_plot['GeneID'].tolist():
                    if '_' in label:
                        label = label.split('_')[0]
                    ax_gene.annotate(label, (x+0.1, y+0.1))
        else:
            for label, x, y in zip(top_by_gene.columns, top_gene_trans[:, 0], top_gene_trans[:, 1]):
                if '_' in label:
                    label = label.split('_')[0]
                ax_gene.annotate(label, (x+0.1, y+0.1))
        if plot:
            plt.show()
        if title != '':
            save_name = '_'.join(title.split(' ')[0:2])
            plt.savefig(os.path.join(path_filename,save_name+'_skpca.pdf'), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(path_filename,'Group0_skpca.pdf'), bbox_inches='tight')
        plt.close('all')
        return top_pca_list
    else:
        return []

def plot_SVD(args,df_by_gene, path_filename, num_genes=100, gene_list_filter=False, title='', plot=False, label_map=False, gene_map = False):
    gene_list = df_by_gene.columns.tolist()
    sns.set_palette("RdBu_r", 10, 1)
    if gene_list_filter:
        sig_by_gene = df_by_gene[gene_list_filter]
        sig_by_cell = sig_by_gene.transpose()
    else:
        sig_by_gene = df_by_gene
        sig_by_cell = sig_by_gene.transpose()
    gene_pca = TruncatedSVD(n_components=3)
    np_by_gene = np.asarray(sig_by_gene)

    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(gene_pca.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=sig_by_gene.columns.tolist())
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len(sig_by_gene.columns.tolist()))
    top_pca_list = Pc_sort_df.index.tolist()

    top_by_gene = df_by_gene[top_pca_list[0:num_genes]]
    gene_top = TruncatedSVD(n_components=2)
    cell_pca = TruncatedSVD(n_components=2)
    top_by_cell = top_by_gene.transpose()
    np_top_gene = np.asarray(top_by_cell)
    np_top_cell = np.asarray(top_by_gene)
    top_cell_trans = cell_pca.fit_transform(np_top_cell)
    top_gene_trans = gene_top.fit_transform(np_top_gene)
    if not np.isnan(top_cell_trans).any():
        fig, (ax_cell, ax_gene) = plt.subplots(2, 1, figsize=(15, 30), sharex=False)
        rect_cell = ax_cell.patch
        rect_gene = ax_gene.patch
        rect_cell.set_facecolor('white')
        rect_gene.set_facecolor('white')
        ax_cell.grid(b=True, which='major', color='grey', linestyle='--', linewidth=0.3)
        ax_gene.grid(b=True, which='major', color='grey', linestyle='--', linewidth=0.3)
        if label_map:
            annotate = args.annotate_cell_pca
            X = [x for x in top_cell_trans[:, 0]]
            Y = [y for y in top_cell_trans[:, 1]]
            labels = [label_map[cell][2] for cell in top_by_cell.columns.tolist()]
            markers = [label_map[cell][1] for cell in top_by_cell.columns.tolist()]
            colors = [label_map[cell][0] for cell in top_by_cell.columns.tolist()]
            label_done = []
            xy_by_color_dict = {}
            for c in set(colors):
                xy_by_color_dict[c] = []
            for X_pos, Y_pos, m, color, l in zip(X, Y, markers, colors, labels):
                if l in label_done:
                    lab = ''
                else:
                    lab= l
                    label_done.append(l)
                xy_by_color_dict[color].append([X_pos,Y_pos])
                ax_cell.scatter(X_pos, Y_pos, marker=m, c=color, label=lab, s=30)
            if args.add_ellipse:
                for c in set(colors):
                    ell, edge = ellip_enclose(np.asarray(xy_by_color_dict[c]), c)
                    ax_cell.add_artist(ell)
                    ax_cell.add_artist(edge)
        else:
            ax_cell.scatter(top_cell_trans[:, 0], top_cell_trans[:, 1], alpha=0.75)
            annotate = args.annotate_cell_pca
        ax_cell.set_xlim([min(top_cell_trans[:, 0])-1, max(top_cell_trans[:, 0]+1)])
        ax_cell.set_ylim([min(top_cell_trans[:, 1])-1, max(top_cell_trans[:, 1]+2)])
        ax_cell.set_title(title+'_cell')
        if label_map:
            handles, labs = ax_cell.get_legend_handles_labels()
            # sort both labels and handles by labels
            labs, handles = zip(*sorted(zip(labs, handles), key=lambda t: t[0]))
            ax_cell.legend(handles, labs, loc='best', ncol=1, prop={'size':12}, markerscale=1.5, frameon=True)
        ax_cell.set_xlabel('PC1')
        ax_cell.set_ylabel('PC2')
        if annotate:
            for label, x, y in zip(top_by_cell.columns, top_cell_trans[:, 0], top_cell_trans[:, 1]):
                ax_cell.annotate(label, (x+0.1, y+0.1))

        if gene_map:
            X = [x for x in top_gene_trans[:, 0]]
            Y = [y for y in top_gene_trans[:, 1]]
            labels = top_by_gene.columns.tolist()
            markers = [gene_map[gene][1] for gene in top_by_gene.columns.tolist()]
            colors = [gene_map[gene][0] for gene in top_by_gene.columns.tolist()]
            xy_by_color_dict = {}
            for c in set(colors):
                xy_by_color_dict[c] = []
            for X_pos, Y_pos, m, color, l in zip(X, Y, markers, colors, labels):
                xy_by_color_dict[color].append([X_pos,Y_pos])
                ax_gene.scatter(X_pos, Y_pos, marker=m, c=color, label = l, s=30)
            if args.add_ellipse:
                for c in set(colors):
                    ell, edge = ellip_enclose(np.asarray(xy_by_color_dict[c]), c)
                    ax_gene.add_artist(ell)
                    ax_gene.add_artist(edge)

        else:
            ax_gene.scatter(top_gene_trans[:, 0], top_gene_trans[:, 1], alpha=0.75)
        ax_gene.set_xlim([min(top_gene_trans[:, 0])-1, max(top_gene_trans[:, 0])+1])
        ax_gene.set_ylim([min(top_gene_trans[:, 1])-1, max(top_gene_trans[:, 1])+2])
        ax_gene.set_title(title+'_gene')
        ax_gene.set_xlabel('PC1')
        ax_gene.set_ylabel('PC2')
        if args.annotate_gene_subset:
            plot_subset_path = os.path.join(os.path.dirname(args.filepath),args.annotate_gene_subset)
            genes_plot = pd.read_csv(plot_subset_path, sep='\t', index_col=False)
            for label, x, y in zip(top_by_gene.columns, top_gene_trans[:, 0], top_gene_trans[:, 1]):
                if label in genes_plot['GeneID'].tolist():
                    if '_' in label:
                        label = label.split('_')[0]
                    ax_gene.annotate(label, (x+0.1, y+0.1))
        else:
            for label, x, y in zip(top_by_gene.columns, top_gene_trans[:, 0], top_gene_trans[:, 1]):
                if '_' in label:
                    label = label.split('_')[0]
                ax_gene.annotate(label, (x+0.1, y+0.1))
        if plot:
            plt.show()
        if title != '':
            save_name = '_'.join(title.split(' ')[0:2])
            plt.savefig(os.path.join(path_filename,save_name+'_TruncatedSVD.pdf'), bbox_inches='tight')
            #plot_url = py.plot_mpl(fig)
        else:
            #plot_url = py.plot_mpl(fig)
            plt.savefig(os.path.join(path_filename,'Group0_TruncatedSVD.pdf'), bbox_inches='tight')
        plt.close('all')
        return top_pca_list
    else:
        return []

#create cell and gene TSNE scatter plots (one pdf)
def plot_TSNE(args,df_by_gene, path_filename, num_genes=100, gene_list_filter=False, title='', plot=False, label_map=False, gene_map = False):
    gene_list = df_by_gene.columns.tolist()
    sns.set_palette("RdBu_r", 10, 1)
    if gene_list_filter:
        sig_by_gene = df_by_gene[gene_list_filter]
        sig_by_cell = sig_by_gene.transpose()
    else:
        sig_by_gene = df_by_gene
        sig_by_cell = sig_by_gene.transpose()
    gene_pca = TruncatedSVD(n_components=3)
    np_by_gene = np.asarray(sig_by_gene)

    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(gene_pca.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=sig_by_gene.columns.tolist())
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len(sig_by_gene.columns.tolist()))
    top_pca_list = Pc_sort_df.index.tolist()

    top_by_gene = df_by_gene[top_pca_list[0:num_genes]]
    gene_top = TSNE(n_components=2, init='pca', random_state=0)
    cell_pca = TSNE(n_components=2, init='pca', random_state=0)
    top_by_cell = top_by_gene.transpose()
    np_top_gene = np.asarray(top_by_cell)
    np_top_cell = np.asarray(top_by_gene)
    top_cell_trans = cell_pca.fit_transform(np_top_cell)
    top_gene_trans = gene_top.fit_transform(np_top_gene)
    if not np.isnan(top_cell_trans).any():
        fig, (ax_cell, ax_gene) = plt.subplots(2, 1, figsize=(15, 30), sharex=False)
        rect_cell = ax_cell.patch
        rect_gene = ax_gene.patch
        rect_cell.set_facecolor('white')
        rect_gene.set_facecolor('white')
        ax_cell.grid(b=True, which='major', color='grey', linestyle='--', linewidth=0.3)
        ax_gene.grid(b=True, which='major', color='grey', linestyle='--', linewidth=0.3)
        if label_map:
            annotate = args.annotate_cell_pca
            X = [x for x in top_cell_trans[:, 0]]
            Y = [y for y in top_cell_trans[:, 1]]
            labels = [label_map[cell][2] for cell in top_by_cell.columns.tolist()]
            markers = [label_map[cell][1] for cell in top_by_cell.columns.tolist()]
            colors = [label_map[cell][0] for cell in top_by_cell.columns.tolist()]
            label_done = []
            xy_by_color_dict = {}
            for c in set(colors):
                xy_by_color_dict[c] = []
            for X_pos, Y_pos, m, color, l in zip(X, Y, markers, colors, labels):
                if l in label_done:
                    lab = ''
                else:
                    lab= l
                    label_done.append(l)
                xy_by_color_dict[color].append([X_pos,Y_pos])
                ax_cell.scatter(X_pos, Y_pos, marker=m, c=color, label=lab, s=30)
            if args.add_ellipse:
                for c in set(colors):
                    ell, edge = ellip_enclose(np.asarray(xy_by_color_dict[c]), c)
                    ax_cell.add_artist(ell)
                    ax_cell.add_artist(edge)
        else:
            ax_cell.scatter(top_cell_trans[:, 0], top_cell_trans[:, 1], alpha=0.75)
            annotate = args.annotate_cell_pca
        ax_cell.set_xlim([min(top_cell_trans[:, 0])-1, max(top_cell_trans[:, 0]+1)])
        ax_cell.set_ylim([min(top_cell_trans[:, 1])-1, max(top_cell_trans[:, 1]+2)])
        ax_cell.set_title(title+'_cell')
        if label_map:
            handles, labs = ax_cell.get_legend_handles_labels()
            # sort both labels and handles by labels
            labs, handles = zip(*sorted(zip(labs, handles), key=lambda t: t[0]))
            ax_cell.legend(handles, labs, loc='best', ncol=1, prop={'size':12}, markerscale=1.5, frameon=True)
        ax_cell.set_xlabel('PC1')
        ax_cell.set_ylabel('PC2')
        if annotate:
            for label, x, y in zip(top_by_cell.columns, top_cell_trans[:, 0], top_cell_trans[:, 1]):
                ax_cell.annotate(label, (x+0.1, y+0.1))

        if gene_map:
            X = [x for x in top_gene_trans[:, 0]]
            Y = [y for y in top_gene_trans[:, 1]]
            labels = top_by_gene.columns.tolist()
            markers = [gene_map[gene][1] for gene in top_by_gene.columns.tolist()]
            colors = [gene_map[gene][0] for gene in top_by_gene.columns.tolist()]
            xy_by_color_dict = {}
            for c in set(colors):
                xy_by_color_dict[c] = []
            for X_pos, Y_pos, m, color, l in zip(X, Y, markers, colors, labels):
                xy_by_color_dict[color].append([X_pos,Y_pos])
                ax_gene.scatter(X_pos, Y_pos, marker=m, c=color, label = l, s=30)
            if args.add_ellipse:
                for c in set(colors):
                    ell, edge = ellip_enclose(np.asarray(xy_by_color_dict[c]), c)
                    ax_gene.add_artist(ell)
                    ax_gene.add_artist(edge)

        else:
            ax_gene.scatter(top_gene_trans[:, 0], top_gene_trans[:, 1], alpha=0.75)
        ax_gene.set_xlim([min(top_gene_trans[:, 0])-1, max(top_gene_trans[:, 0])+1])
        ax_gene.set_ylim([min(top_gene_trans[:, 1])-1, max(top_gene_trans[:, 1])+2])
        ax_gene.set_title(title+'_gene')
        ax_gene.set_xlabel('PC1')
        ax_gene.set_ylabel('PC2')
        if args.annotate_gene_subset:
            plot_subset_path = os.path.join(os.path.dirname(args.filepath),args.annotate_gene_subset)
            genes_plot = pd.read_csv(plot_subset_path, sep='\t', index_col=False)
            for label, x, y in zip(top_by_gene.columns, top_gene_trans[:, 0], top_gene_trans[:, 1]):
                if label in genes_plot['GeneID'].tolist():
                    if '_' in label:
                        label = label.split('_')[0]
                    ax_gene.annotate(label, (x+0.1, y+0.1))
        else:
            for label, x, y in zip(top_by_gene.columns, top_gene_trans[:, 0], top_gene_trans[:, 1]):
                if '_' in label:
                    label = label.split('_')[0]
                ax_gene.annotate(label, (x+0.1, y+0.1))
        if plot:
            plt.show()
        if title != '':
            save_name = '_'.join(title.split(' ')[0:2])
            plt.savefig(os.path.join(path_filename,save_name+'_TSNE.pdf'), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(path_filename,'Group0_TSNE.pdf'), bbox_inches='tight')
        plt.close('all')
        return top_pca_list
    else:
        return []

#create cell and gene TSNE scatter plots (one pdf)
def plot_kmeans(args, df_by_gene, path_filename, kmeans_range, num_genes=100, gene_list_filter=False, title='', plot=False, label_map=False, gene_map = False, run_sig_test=False):
    from sklearn.metrics import silhouette_samples, silhouette_score
    import matplotlib.cm as cm
    gene_list = df_by_gene.columns.tolist()
    sns.set_palette("RdBu_r", 10, 1)
    if gene_list_filter:
        sig_by_gene = df_by_gene[gene_list_filter]
        sig_by_cell = sig_by_gene.transpose()
    else:
        sig_by_gene = df_by_gene
        sig_by_cell = sig_by_gene.transpose()
    gene_pca = TruncatedSVD(n_components=3)
    np_by_gene = np.asarray(sig_by_gene)

    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(gene_pca.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=sig_by_gene.columns.tolist())
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len(sig_by_gene.columns.tolist()))
    top_pca_list = Pc_sort_df.index.tolist()

    top_by_gene = df_by_gene[top_pca_list[0:num_genes]]

    gene_top = TruncatedSVD(n_components=2)
    cell_pca = TruncatedSVD(n_components=2)
    top_by_cell = top_by_gene.transpose()
    np_top_gene = np.asarray(top_by_cell)
    np_top_cell = np.asarray(top_by_gene)
    top_cell_trans = cell_pca.fit_transform(np_top_cell)
    top_gene_trans = gene_top.fit_transform(np_top_gene)
    range_n_clusters = range(kmeans_range[0],kmeans_range[1])

    for n_clusters in range_n_clusters:
        # Create a subplot with 1 row and 2 columns
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(18, 7)

        ax1.set_xlim([-0.1, 1])
        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax1.set_ylim([0, len(np_top_cell) + (n_clusters + 1) * 10])
        #cluster cell PCA
        cell_clusterer = KMeans(n_clusters=n_clusters)
        top_cell_pred = cell_clusterer.fit_predict(top_cell_trans)
        #cluster gene PCA
        gene_clusterer = KMeans(n_clusters=n_clusters)
        top_gene_pred = gene_clusterer.fit_predict(top_gene_trans)

        pred_dict = {'SampleID':top_by_cell.columns, 'GroupID':['kmeans_'+str(p) for p in top_cell_pred]}
        df_pred = pd.DataFrame(pred_dict)
        cell_group_path = os.path.join(path_filename,'kmeans_cell_groups_'+str(n_clusters)+'.txt')
        df_pred.to_csv(cell_group_path, sep = '\t')
        #compute silouette averages and values
        silhouette_avg_cell = silhouette_score(top_cell_trans, top_cell_pred)
        print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg_cell)
        silhouette_avg_gene = silhouette_score(top_gene_trans, top_gene_pred)


        sample_silhouette_values_cell = silhouette_samples(top_cell_trans, top_cell_pred)
        sample_silhouette_values_gene = silhouette_samples(top_gene_trans, top_gene_pred)

        y_lower = 10
        for i in range(n_clusters):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            ith_cluster_silhouette_values = \
                sample_silhouette_values_cell[top_cell_pred == i]

            ith_cluster_silhouette_values.sort()

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.spectral(float(i) / n_clusters)
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                              0, ith_cluster_silhouette_values,
                              facecolor=color, edgecolor=color, alpha=0.7)

            # Label the silhouette plots with their cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples

        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")

        # The vertical line for average silhoutte score of all the values
        ax1.axvline(x=silhouette_avg_cell, color="red", linestyle="--")

        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

        # 2nd Plot showing the actual clusters formed
        colors = cm.spectral(top_cell_pred.astype(float) / n_clusters)
        ax2.scatter(top_cell_trans[:, 0], top_cell_trans[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                    c=colors)

        # Labeling the clusters
        centers = cell_clusterer.cluster_centers_
        # Draw white circles at cluster centers
        ax2.scatter(centers[:, 0], centers[:, 1],
                    marker='o', c="white", alpha=1, s=200)

        for i, c in enumerate(centers):
            ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1, s=50)

        ax2.set_title("The visualization of the clustered data.")
        ax2.set_xlabel("Feature space for the 1st feature")
        ax2.set_ylabel("Feature space for the 2nd feature")

        plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                      "with n_clusters = %d" % n_clusters),
                     fontsize=14, fontweight='bold')

        plt.savefig(os.path.join(path_filename,'Group0_kmeans_'+str(n_clusters)+'_clusters.pdf'), bbox_inches='tight')
        plt.close('all')

        #use colors to make label map compatable with heatmap
        color_dict ={}
        markers = ['o', 'v','D','*','x','h', 's','p','8','^','>','<', 'd','o', 'v','D','*','x','h', 's','p','8','^','>','<', 'd']
        color_dict =dict(zip(top_by_cell.columns, zip(colors,[markers[pred] for pred in top_cell_pred],['kmeans_'+str(p) for p in top_cell_pred])))
        group_color_dict = dict(zip(['kmeans_'+str(p) for p in top_cell_pred],zip(colors,[markers[pred] for pred in top_cell_pred])))
        #run heatmap with kmeans clustering and colors
        top_pca_by_gene, top_pca = return_top_pca_gene(df_by_gene, num_genes=args.gene_number)
        top_pca_by_cell = top_pca_by_gene.transpose()
        cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, top_pca, top_pca_by_gene, path_filename, num_to_plot=args.gene_number, title= 'kmeans_label_with_'+str(n_clusters)+'_clusters',label_map=color_dict)
        if run_sig_test:
            multi_group_sig(args, df_by_gene.transpose(), cell_group_path, path_filename, group_color_dict, from_kmeans=str(n_clusters))



def clust_heatmap(args, gene_list, df_by_gene, path_filename, num_to_plot, title='', plot=False, label_map=False, gene_map=False, fontsize=18):
    cell_list = df_by_gene.index.tolist()
    cell_num =len(cell_list)
    longest_side = max(num_to_plot,cell_num*2)
    if longest_side == num_to_plot:
        sns.set(context= 'poster', font_scale = .4*(num_to_plot/100))
        width_heatmap = min(28+round(cell_num/50),42+round(cell_num/40))
        len_heatmap = min(43+round(num_to_plot/10),58+round(num_to_plot/30))
        title_set = 1.15
    else:
        sns.set(context= 'poster', font_scale = .6*(cell_num/120))
        width_heatmap = min(42+round(cell_num/9),68+round(cell_num/40))
        len_heatmap = min(47+round(num_to_plot/8),50+round(num_to_plot/30))
        title_set = 1.12
    font = {'size'   : fontsize}

    plt.rc('font', **font)
    if len(str(args.z_direction)) > 1:
        z_choice = str(args.z_direction)
        if z_choice != 'None':
            sys.exit('Please enter a valid option (0, 1, or None) for z_direction')
    else:
        z_choice = int(args.z_direction)
        if z_choice != 0 and z_choice != 1:
            sys.exit('Please enter a valid option (0, 1, or None) for z_direction')
    cmap = sns.diverging_palette(255, 10, s=99, sep=1, as_cmap=True)

    cluster_df = df_by_gene[gene_list[0:num_to_plot]].transpose()
    cluster_df[abs(cluster_df)<3e-12] = 0.0
    try:
        cg = sns.clustermap(cluster_df, method=args.method, metric=args.metric, z_score=z_choice, figsize=(width_heatmap, len_heatmap), cmap =cmap)
        col_order = cg.dendrogram_col.reordered_ind
        row_order = cg.dendrogram_row.reordered_ind
        if label_map and gene_map:
            Xlabs = [cell_list[i] for i in col_order]
            Xcolors = [label_map[cell][0] for cell in Xlabs]
            col_colors = pd.DataFrame({'Cell Groups': Xcolors},index=Xlabs)
            Xgroup_labels = [label_map[cell][2] for cell in Xlabs]
            Ylabs = [gene_list[i] for i in row_order]
            Ycolors = [gene_map[gene][0] for gene in Ylabs]
            Ygroup_labels= [gene_map[gene][2] for gene in Ylabs]
            row_colors = pd.DataFrame({'Gene Groups': Ycolors},index=Ylabs)
            cg = sns.clustermap(cluster_df, method=args.method, metric=args.metric, z_score=z_choice,row_colors=row_colors, col_colors=col_colors, figsize=(width_heatmap, len_heatmap), cmap =cmap)
        elif label_map:
            Xlabs = [cell_list[i] for i in col_order]
            Xcolors = [label_map[cell][0] for cell in Xlabs]
            Xgroup_labels = [label_map[cell][2] for cell in Xlabs]
            col_colors = pd.DataFrame({'Cell Groups': Xcolors},index=Xlabs)
            cg = sns.clustermap(cluster_df, method=args.method, metric=args.metric, z_score=z_choice, col_colors=col_colors, figsize=(width_heatmap, len_heatmap), cmap =cmap)
        elif gene_map:
            Ylabs = [gene_list[i] for i in row_order]
            Ycolors = [gene_map[gene][0] for gene in Ylabs]
            Ygroup_labels= [gene_map[gene][2] for gene in Ylabs]
            row_colors = pd.DataFrame({'Gene Groups': Ycolors},index=Ylabs)
            cg = sns.clustermap(cluster_df, method=args.method, metric=args.metric, z_score=z_choice,row_colors=row_colors, figsize=(width_heatmap, len_heatmap), cmap =cmap)

        cg.ax_heatmap.set_title(title, y=title_set)
        cg.cax.set_title('Z-score')
        if label_map:
            leg_handles_cell =[]
            group_seen_cell = []
            for xtick, xcolor, xgroup_name in zip(cg.ax_heatmap.get_xticklabels(), Xcolors, Xgroup_labels):
                xtick.set_color(xcolor)
                xtick.set_rotation(270)
                xtick.set_fontsize(fontsize)
                if xgroup_name not in group_seen_cell:
                    leg_handles_cell.append(patches.Patch(color=xcolor, label=xgroup_name))
                    group_seen_cell.append(xgroup_name)
        else:
            for xtick in cg.ax_heatmap.get_xticklabels():
                xtick.set_rotation(270)
                xtick.set_fontsize(fontsize)
        if gene_map:
            leg_handles_gene =[]
            group_seen_gene = []
            for ytick, ycolor, ygroup_name in zip(cg.ax_heatmap.get_yticklabels(), list(reversed(Ycolors)), list(reversed(Ygroup_labels))):
                ytick.set_color(ycolor)
                ytick.set_rotation(0)
                ytick.set_fontsize(fontsize)
                if ygroup_name not in group_seen_gene:
                    leg_handles_gene.append(patches.Patch(color=ycolor, label=ygroup_name))
                    group_seen_gene.append(ygroup_name)
        else:
            for ytick in cg.ax_heatmap.get_yticklabels():
                ytick.set_rotation(0)
                ytick.set_fontsize(fontsize)
        if gene_map and label_map:
            gene_legend = cg.ax_heatmap.legend(handles=leg_handles_gene, loc=2, bbox_to_anchor=(1.04, 0.8), title='Gene groups', prop={'size':fontsize})
            plt.setp(gene_legend.get_title(),fontsize=fontsize)
            cg.ax_heatmap.add_artist(gene_legend)
            cell_legend = cg.ax_heatmap.legend(handles=leg_handles_cell, loc=2, bbox_to_anchor=(1.04, 1), title='Cell groups', prop={'size':fontsize})
            plt.setp(cell_legend.get_title(),fontsize=fontsize)
            #cg.ax_heatmap.add_artist(cell_legend)
        elif label_map:
            cell_legend = cg.ax_heatmap.legend(handles=leg_handles_cell, loc=2, bbox_to_anchor=(1.04, 1), title='Cell groups', prop={'size':fontsize})
            plt.setp(cell_legend.get_title(),fontsize=fontsize)
        elif gene_map:
            gene_legend = cg.ax_heatmap.legend(handles=leg_handles_gene, loc=2, bbox_to_anchor=(1.04, 0.8), title='Gene groups', prop={'size':fontsize})
            plt.setp(gene_legend.get_title(),fontsize=fontsize)
        if plot:
            plt.show()
        cell_linkage = cg.dendrogram_col.linkage

        link_mat = pd.DataFrame(cell_linkage,
                    columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
                    index=['cluster %d' %(i+1) for i in range(cell_linkage.shape[0])])
        if title != '':
            save_name = '_'.join(title.split(' ')[0:2])
            #plot_url = py.plot_mpl(cg)
            cg.savefig(os.path.join(path_filename, save_name+'_heatmap.pdf'), bbox_inches='tight')
        else:
            #plot_url = py.plot_mpl(cg)
            cg.savefig(os.path.join(path_filename,'Group0_Heatmap_all_cells.pdf'), bbox_inches='tight')
        plt.close('all')
        return cell_linkage, df_by_gene[gene_list[0:num_to_plot]], col_order
    except FloatingPointError:
        print('Linkage distance has too many zeros. Filter to remove non-expressed genes in order to produce heatmap. Heatmap with '+ str(len(cell_list))+' will not be created.')
        return False, False, False

def make_subclusters(args, cc, log2_expdf_cell, log2_expdf_cell_full, path_filename, base_name, gene_corr_list, label_map=False, gene_map=False, cluster_size=20, group_colors=False):
    '''
    Walks a histogram branch map 'cc' and does PCA (SVD), heatmap and correlation search for each non-overlapping
    tree branch. Stops at defined cluster_size (default is 20).
    '''
    #initial cell group is parent
    parent = cc[0][1]
    #p_num is the number of cells in the parent group
    p_num = cc[0][0]
    #l_nums is number of members of each leaf of tree
    l_nums = [x[0] for x in cc]
    #cell list is the list of list of cells in each leaf of tree
    c_lists = [c[1] for c in cc]
    #Group ID will increment with each group so that each subcluter has a unique ID
    group_ID = 0

    for num_members, cell_list in zip(l_nums, c_lists):
        #run all cell groups that are subgroups of the parent and greater than or equal to the cutoff cluster_size
        if num_members < p_num and num_members >= cluster_size:
            group_ID+=1
            #save name for all files generated within this cluster i.e. 'Group_2_with_105_cells_heatmap.pdf'
            current_title = 'Group_'+str(group_ID)+'_with_'+str(num_members)+'_cells'
            cell_subset = log2_expdf_cell[list(set(cell_list))]
            gene_subset = cell_subset.transpose()
            gene_subset = gene_subset.loc[:,(gene_subset!=0).any()]
            full_cell_subset = log2_expdf_cell_full[list(set(cell_list))]
            full_gene_subset = full_cell_subset.transpose()
            full_gene_subset = full_gene_subset.loc[:,(full_gene_subset!=0).any()]
            norm_df_cell1 = np.exp2(full_cell_subset)
            norm_df_cell = norm_df_cell1 -1
            norm_df_cell.to_csv(os.path.join(path_filename, base_name+'_'+current_title+'_matrix.txt'), sep = '\t', index_col=0)
            if gene_map:
                top_pca_by_gene, top_pca = return_top_pca_gene(gene_subset, num_genes=args.gene_number)
                plot_SVD(args,gene_subset, path_filename, num_genes=len(gene_subset.columns.tolist()), title=current_title, plot=False, label_map=label_map, gene_map=gene_map)
            else:
                top_pca_by_gene, top_pca = return_top_pca_gene(full_gene_subset, num_genes=args.gene_number)
                plot_SVD(args,full_gene_subset, path_filename, num_genes=int(args.gene_number), title=current_title, plot=False, label_map=label_map)
            if len(top_pca)<args.gene_number:
                plot_num = len(top_pca)
            else:
                plot_num = args.gene_number
            if top_pca != []:
                top_pca_by_cell = top_pca_by_gene.transpose()
                #if no_corr flag is provided (False) no correlation plots will be made
                if args.no_corr:

                    if gene_corr_list != []:
                        top_genes_search = top_pca[0:50]
                        corr_plot(top_genes_search, full_gene_subset, path_filename, num_to_plot=3, gene_corr_list= gene_corr_list, title = current_title, label_map=label_map)
                    else:
                        top_genes_search = top_pca[0:50]
                        corr_plot(top_genes_search, full_gene_subset, path_filename, num_to_plot=3, title = current_title, label_map=label_map)
                if gene_map:
                    cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, top_pca, top_pca_by_gene, path_filename, num_to_plot=plot_num, title=current_title, plot=False, label_map=label_map, gene_map = gene_map)
                else:
                    cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, top_pca, top_pca_by_gene, path_filename,num_to_plot=plot_num, title=current_title, plot=False, label_map=label_map)
                plt.close('all')
            else:
                print('Search for top genes by PCA failed in '+current_title+'. No plots will be generated for this subcluster. ')
                pass

def clust_stability(args, log2_expdf_gene, path_filename, iterations, label_map=False):
    sns.set(context='poster', font_scale = 1)
    sns.set_palette("RdBu_r")
    stability_ratio = []
    total_genes = len(log2_expdf_gene.columns.tolist())
    end_num = 1000
    iter_list = range(100,int(round(end_num)),int(round(end_num/iterations)))
    for gene_number in iter_list:
        title= str(gene_number)+' genes plot.'
        top_pca_by_gene, top_pca = return_top_pca_gene(df_by_gene, num_genes=gene_number)
        top_pca_by_cell = top_pca_by_gene.transpose()
        cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, top_pca, top_pca_by_gene, num_to_plot=gene_number, title=title, label_map=label_map)
        if gene_number == 100:
            s1 = col_order
            s0 = col_order
        else:
            s2= col_order
            sm_running = difflib.SequenceMatcher(None,s1,s2)
            sm_first = difflib.SequenceMatcher(None,s0,s2)
            stability_ratio.append((sm_running.ratio(), sm_first.ratio()))
            s1=col_order
        plt.close()
    x= iter_list[1:]
    f, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    y1= [m[0] for m in stability_ratio]
    y2= [m[1] for m in stability_ratio]
    sns.barplot(x, y1, palette="RdBu_r", ax=ax1)
    ax1.set_ylabel('Running ratio (new/last)')
    sns.barplot(x, y2, palette="RdBu_r", ax=ax2)
    ax2.set_ylabel('Ratio to 100')
    plt.savefig(os.path.join(path_filename,'clustering_stability.pdf'), bbox_inches='tight')
    plt.show()
    plt.close('all')
    return stability_ratio

#run correlation matrix and save only those above threshold
def run_corr(df_by_gene, title, path_filename, method_name='pearson', sig_threshold= 0.5, run_new=True, min_period=3, save_corrs=False):
    if run_new:
        if len(df_by_gene.columns.tolist())>5000:
            df_by_gene, top_pca_list = return_top_pca_gene(df_by_gene, num_genes=5000)
        if method_name != 'kendall':
            corr_by_gene = df_by_gene.corr(method=method_name, min_periods=min_period)
        else:
            corr_by_gene = df_by_gene.corr(method=method_name)
        df_by_cell = df_by_gene.transpose()
        corr_by_cell = df_by_cell.corr()

        cor = corr_by_gene
        cor.loc[:,:] =  np.tril(cor.values, k=-1)
        cor = cor.stack()
        corr_by_gene_pos = cor[cor >=sig_threshold]
        corr_by_gene_neg = cor[cor <=(sig_threshold*-1)]
    else:
        corr_by_g_pos =  open(os.path.join(path_filename,'gene_correlations_sig_pos_'+method_name+'.p'), 'rb')
        corr_by_g_neg =  open(os.path.join(path_filename,'gene_correlations_sig_neg_'+method_name+'.p'), 'rb')
        corr_by_gene_pos = pickle.load(corr_by_g_pos)
        corr_by_gene_neg = pickle.load(corr_by_g_neg)
    if save_corrs:
        with open(os.path.join(path_filename,'gene_correlations_sig_neg_'+method_name+'.p'), 'wb') as fp:
            pickle.dump(corr_by_gene_neg, fp)
        with open(os.path.join(path_filename,'gene_correlations_sig_pos_'+method_name+'.p'), 'wb') as fp0:
            pickle.dump(corr_by_gene_pos, fp0)
        with open(os.path.join(path_filename,'by_gene_corr.p'), 'wb') as fp1:
            pickle.dump(corr_by_gene, fp1)
        with open(os.path.join(path_filename,'by_cell_corr.p'), 'wb') as fp2:
            pickle.dump(corr_by_cell, fp2)


    cor_pos_df = pd.DataFrame(corr_by_gene_pos)
    cor_neg_df = pd.DataFrame(corr_by_gene_neg)
    sig_corr = cor_pos_df.append(cor_neg_df)
    sig_corrs = pd.DataFrame(sig_corr[0], columns=["corr"])

    if run_new:
        sig_corrs.to_csv(os.path.join(path_filename, title+'_counts_corr_sig_'+method_name+'.txt'), sep = '\t')
    return sig_corrs

#finds most correlated gene groups that are not overlapping
def find_top_corrs(terms_to_search, sig_corrs, num_to_return, gene_corr_list = []):
    all_corrs_list = []
    best_corrs_list = []
    for term_to_search in terms_to_search:
        corr_tup = [(term_to_search, 1)]
        for index, row in sig_corrs.iterrows():
            if term_to_search in index:
                if index[0]==term_to_search:
                    corr_tup.append((index[1],row['corr']))
                else:
                    corr_tup.append((index[0],row['corr']))
        all_corrs_list.append(corr_tup)
    all_corrs_list.sort(key=len, reverse=True)
    good_count = 0
    corr_genes_seen = []
    while good_count <= num_to_return:
        for i, corrs in enumerate(all_corrs_list):
            if corrs[0][0] not in corr_genes_seen:
                best_corrs_list.append(corrs)
                good_count+=1
            for g, c in corrs:
                if g not in corr_genes_seen and '-' not in str(c):
                    corr_genes_seen.append(g)
    if gene_corr_list != []:
        search_corrs = []
        for term in gene_corr_list:
            corr_tup = [(term, 1)]
            for index, row in sig_corrs.iterrows():
                if term in index:
                    if index[0]==term:
                        corr_tup.append((index[1],row['corr']))
                    else:
                        corr_tup.append((index[0],row['corr']))
            search_corrs.append(corr_tup)
        best_corrs_list = search_corrs+best_corrs_list
        return best_corrs_list[0:num_to_return+len(gene_corr_list)+1]
    else:
        return best_corrs_list[0:num_to_return]


#corr_plot finds and plots all correlated genes, log turns on log scale, sort plots the genes in the rank order of the gene searched
def corr_plot(terms_to_search, df_by_gene_corr, path_filename, title, num_to_plot, gene_corr_list = [], label_map=False, log=False, sort=True, sig_threshold=0.5):
    size_cells = len(df_by_gene_corr.index.tolist())
    figlen=int(size_cells/12)
    if figlen < 15:
        figlen = 15
    ncol = int(figlen/3.2)

    if size_cells <100:
        sig_threshold = -0.137*math.log(size_cells)+1.1322
    sig_corrs = run_corr(df_by_gene_corr, title, path_filename, sig_threshold=sig_threshold)
    corr_list = find_top_corrs(terms_to_search, sig_corrs, num_to_plot, gene_corr_list=gene_corr_list)
    for corr_tup in corr_list:
        term_to_search = corr_tup[0][0]
        corr_tup.sort(key=itemgetter(1), reverse=True)
        corr_df = pd.DataFrame(corr_tup, columns=['GeneID', 'Correlation'])
        corr_df.to_csv(os.path.join(path_filename, title+'_Corr_w_'+term_to_search+'_list.txt'), sep = '\t', index=False)

        to_plot = [x[0] for x in corr_tup]
        sns.set_palette(sns.cubehelix_palette(len(to_plot), start=1, rot=-.9, reverse=True))
        sns.set_context("notebook", font_scale=.8, rc={"lines.linewidth": 1})
        try:
            sorted_df = df_by_gene_corr.sort_values(by=[term_to_search])
            log2_df = np.log2(df_by_gene_corr[to_plot])
            sorted_log2_df=np.log2(sorted_df[to_plot])
            ylabel='Counts (log2)'
            if sort and log:
                ax = sorted_log2_df.plot(figsize = (figlen,10))
                xlabels = sorted_log2_df[to_plot].index.values
            elif sort:
                ax = sorted_df[to_plot].plot(figsize = (figlen,10))
                xlabels = sorted_df[to_plot].index.values
            elif log:
                ax = log2_df.plot(figsize = (figlen,10))
                ylabel= 'log2 FPKM'
                xlabels = log2_df.index.values
            else:
                ax = df_by_gene_corr[to_plot].plot(figsize = (figlen,10))
                xlabels = df_by_gene_corr[to_plot].index.values
            ax.set_xlabel('Cell #')
            ax.set_ylabel(ylabel)
            ax.set_title('Correlates with '+term_to_search, loc='right')
            ax.xaxis.set_minor_locator(LinearLocator(numticks=len(xlabels)))
            if label_map:
                ax.set_xticklabels(xlabels, minor=True, rotation='vertical', fontsize=3)
                Xcolors = [label_map[cell][0] for cell in xlabels]
                group_labels = [label_map[cell][2] for cell in xlabels]
                group_seen = []
                leg_handles = []
                for xtick, xcolor, group_name in zip(ax.get_xticklabels(which='minor'), Xcolors, group_labels):
                    xtick.set_color(xcolor)
                    xtick.set_rotation(90)
                    if group_name not in group_seen:
                        leg_handles.append(patches.Patch(color=xcolor, label=group_name))
                        group_seen.append(group_name)

            else:
                ax.set_xticklabels(xlabels, minor=True, rotation='vertical', fontsize=3)
            ax.set_ylim([0, df_by_gene_corr[to_plot].values.max()])
            ax.xaxis.set_major_formatter(ticker.NullFormatter())
            ax.tick_params(axis='x', which ='minor', labelsize=9)
            #scale bbox anchoring to account for number of correlated genes and plot size
            if len(corr_tup)>1:
                bbox_height = float(1E-13)*pow(len(corr_tup),6) - float(7E-11)*pow(len(corr_tup),5) + float(1E-8)*pow(len(corr_tup),4) - float(8E-7)*pow(len(corr_tup),3) - float(3E-5)*pow(len(corr_tup),2) + 0.0086*len(corr_tup) + 1.0042
            else:
                bbox_height = 1.05
            l_labels = [str(x[0])+' '+"%.2f" % x[1] for x in corr_tup]
            if label_map:
                first_legend = ax.legend(l_labels, loc='upper left', bbox_to_anchor=(0.01, bbox_height), ncol=ncol, prop={'size':10})
                ax = plt.gca().add_artist(first_legend)
                plt.legend(handles=leg_handles, loc='upper right', bbox_to_anchor=(0.9, bbox_height+.1))
            else:
                ax.legend(l_labels, loc='upper left', bbox_to_anchor=(0.01, bbox_height), ncol=ncol, prop={'size':10})
            fig = plt.gcf()
            fig.subplots_adjust(bottom=0.08, top=0.95, right=0.98, left=0.03)
            plt.savefig(os.path.join(path_filename, title+'_corr_with_'+term_to_search+'.pdf'), bbox_inches='tight')
            plt.close('all')
        except KeyError:
            print(term_to_search+' not in this matrix.')
            pass

'''Compares each defined cell group to each other group and returns all genes with p-value and adjusted p-value.
Also creates a best_gene_list with the top genes between each group, by adjusted p-value.
Also creates a barplot of top significant genes between groups (can be unique or not).
'''
def multi_group_sig(args, full_by_cell_df, cell_group_filename, path_filename, color_dict_cell, sig_to_plot = 20, from_kmeans=''):
    #create seperate file for all of the significance files
    if from_kmeans == '':
        multi_sig_filename = os.path.join(path_filename,'user_defined_group_pairwise_significance_files')
    else:
        multi_sig_filename = os.path.join(path_filename,from_kmeans+'_group_pairwise_significance_files')
    try:
        os.mkdir(multi_sig_filename)
    except OSError:
        print(multi_sig_filename+' already exists. Files will be overwritten.')
    plot_pvalue = False
    from rpy2.robjects.packages import importr
    from rpy2.robjects.vectors import FloatVector
    stats = importr('stats')
    cell_groups_df = pd.read_csv(open(cell_group_filename,'rU'), sep=None, engine='python')
    group_name_list = list(set(cell_groups_df['GroupID']))
    group_pairs = list(set(itertools.permutations(group_name_list,2)))
    gene_list = full_by_cell_df.index.tolist()
    cell_group_ident_0 = zip(cell_groups_df['SampleID'],cell_groups_df['GroupID'])
    cell_group_ident= [c for c in cell_group_ident_0]
    barplot_dict = {}
    by_gene_df = full_by_cell_df.transpose()
    best_gene_list = []
    best_gene_groups = []
    best_vs_list =[]
    best_pvalue_list = []
    for name in group_name_list:
        barplot_dict[name] = {'genes':[], 'pvalues':[], 'fold_change':[], 'Vs':[]}
    for gp in group_pairs:
        index_list = []
        sig_gene_list = []
        gp_vs_all_seen = []
        sig_vs_all_gene_list = []
        cell_list1 = [c[0] for c in cell_group_ident if c[1] == gp[0]]
        cell_list2 = [c[0] for c in cell_group_ident if c[1] == gp[1]]
        cell1_present = [c for c in set(cell_list1) if c in set(full_by_cell_df.columns.tolist())]
        cell2_present = [c for c in set(cell_list2) if c in set(full_by_cell_df.columns.tolist())]
        all_other_cells =  [c for c in set(full_by_cell_df.columns.tolist()) if c not in set(cell_list1)]
        if cell1_present != [] and cell2_present != []:
            df_by_cell_1 = full_by_cell_df.ix[:,cell1_present]
            df_by_cell_2 = full_by_cell_df.ix[:,cell2_present]
            df_by_cell_other = full_by_cell_df.ix[:,all_other_cells]
            df_by_gene_1 = df_by_cell_1.transpose()
            df_by_gene_2 = df_by_cell_2.transpose()
            df_by_gene_other = df_by_cell_other.transpose()
            for g in gene_list:
                g_pvalue = scipy.stats.f_oneway(df_by_gene_1[g], df_by_gene_2[g])
                if gp[0] not in gp_vs_all_seen:
                    g_pvalue_all = scipy.stats.f_oneway(df_by_gene_1[g], df_by_gene_other[g])
                if g_pvalue[0] > 0 and g_pvalue[1] <= 1:
                    if g not in [s[0] for s in sig_gene_list]:
                        sig_gene_list.append([g, g_pvalue[1]])
                        if gp[0] not in gp_vs_all_seen:
                            sig_vs_all_gene_list.append([g, g_pvalue_all[1]])

            sig_gene_list.sort(key=lambda tup: tup[1])
            if gp[0] not in gp_vs_all_seen:
                sig_vs_all_gene_list.sort(key=lambda tup: tup[1])
                pvalues_all = [p[1] for p in sig_vs_all_gene_list]

                p_adjust_all = stats.p_adjust(FloatVector(pvalues_all), method = 'BH')
                gene_index_all = [ge[0] for ge in sig_vs_all_gene_list]
                mean_log2_exp_list_all = []
                sig_1_2_list_all = []
                mean1_list_all = []
                mean2_list_all = []
                for sig_gene in gene_index_all:
                    sig_gene_df_all = by_gene_df.ix[:,sig_gene]
                    mean_log2_exp_list_all.append(sig_gene_df_all.mean())
                    sig_cell_df_all = sig_gene_df_all.transpose()
                    mean_cell1_all = sig_cell_df_all[cell1_present].mean()
                    mean1_list_all.append(mean_cell1_all)
                    mean_cell_other = sig_cell_df_all[all_other_cells].mean()
                    mean2_list_all.append(mean_cell_other)
                    ratio_1_other = (mean_cell1_all+1)/(mean_cell_other+1)
                    sig_1_2_list_all.append(ratio_1_other)
                sig_df_vs_other = pd.DataFrame({'pvalues':pvalues_all,'adjusted_p_values':p_adjust_all,'mean_all':mean_log2_exp_list_all, 'mean_'+gp[0]:mean1_list_all, 'mean_all_other':mean2_list_all, 'ratio '+gp[0]+' to everything':sig_1_2_list_all}, index=gene_index_all)
                sig_df_vs_other.to_csv(os.path.join(multi_sig_filename,'sig_'+gp[0]+'_VS_all_other_pvalues.txt'), sep = '\t')
                gp_vs_all_seen.append(gp[0])
            pvalues = [p[1] for p in sig_gene_list]
            p_adjust = stats.p_adjust(FloatVector(pvalues), method = 'BH')
            gene_index = [ge[0] for ge in sig_gene_list]

            mean_log2_exp_list = []
            sig_1_2_list = []
            mean1_list = []
            mean2_list = []
            for sig_gene in gene_index:
                sig_gene_df = by_gene_df.ix[:,sig_gene]
                mean_log2_exp_list.append(sig_gene_df.mean())
                sig_cell_df = sig_gene_df.transpose()
                mean_cell1 = sig_cell_df[cell1_present].mean()
                mean1_list.append(mean_cell1)
                mean_cell2 = sig_cell_df[cell2_present].mean()
                mean2_list.append(mean_cell2)
                ratio_1_2 = (mean_cell1+1)/(mean_cell2+1)
                sig_1_2_list.append(ratio_1_2)
            sig_df = pd.DataFrame({'pvalues':pvalues,'adjusted_p_values':p_adjust,'mean_all':mean_log2_exp_list,'mean_'+gp[0]:mean1_list, 'mean_'+gp[1]:mean2_list, 'ratio '+gp[0]+' to '+gp[1]:sig_1_2_list}, index=gene_index)
            cell_names_df = pd.DataFrame({gp[0]+'_cells':pd.Series(cell1_present, index=range(len(cell1_present))), gp[1]+'_cells2':pd.Series(cell2_present, index=range(len(cell2_present)))})
            sig_df.to_csv(os.path.join(multi_sig_filename,'sig_'+gp[0]+'_VS_'+gp[1]+'_pvalues.txt'), sep = '\t')
            cell_names_df.to_csv(os.path.join(multi_sig_filename,'sig_'+gp[0]+'_VS_'+gp[1]+'_cells.txt'), sep = '\t')


            top_fc_df1 = sig_df.loc[(sig_df['ratio '+gp[0]+' to '+gp[1]]>1.3)]
            new_fc_list = [-1/x if x <0.3 else x for x in top_fc_df1['ratio '+gp[0]+' to '+gp[1]]]
            top_fc_df1.loc[:,'ratio '+gp[0]+' to '+gp[1]] = new_fc_list
            top_fc_df = top_fc_df1.sort_values(by='adjusted_p_values',axis=0, ascending=True)
            genes = top_fc_df.index.tolist()
            pvalues = top_fc_df['adjusted_p_values'].tolist()
            fc = top_fc_df['ratio '+gp[0]+' to '+gp[1]].tolist()
            z = zip(genes,pvalues,fc)
            z_all = [s for s in z]
            if args.sig_unique:
                top_t = [g for g in z_all if g[0] not in barplot_dict[gp[0]]['genes']]
            else:
                top_t = [g for g in z_all if g[0]]
            if args.exclude_genes:
                hu_cc_gene_df = pd.DataFrame.from_csv(os.path.join(os.path.dirname(args.filepath),args.exclude_genes), sep='\t', header=0, index_col=False)
                exclude_list = hu_cc_gene_df['GeneID'].tolist()
                top_t2 = [g for g in top_t if g[0] not in barplot_dict[gp[0]]['genes'] and g[0] not in exclude_list]
            else:
                top_t2 = top_t
            top = [list(t) for t in zip(*top_t2)]
            sig_to_plot = min(len(top[0]),sig_to_plot)
            if sig_to_plot != 0:
                barplot_dict[gp[0]]['genes']= barplot_dict[gp[0]]['genes']+[str(gene.strip(' ')) for gene in top[0][0:sig_to_plot]]
                barplot_dict[gp[0]]['pvalues']= barplot_dict[gp[0]]['pvalues']+top[1][0:sig_to_plot]
                barplot_dict[gp[0]]['fold_change']= barplot_dict[gp[0]]['fold_change']+top[2][0:sig_to_plot]
                barplot_dict[gp[0]]['Vs']= barplot_dict[gp[0]]['Vs']+ ['significance vs '+gp[1] for x in range(0,len(top[0][0:sig_to_plot]))]
                best_gene_list = best_gene_list+top[0][0:sig_to_plot]
                best_gene_groups = best_gene_groups+[gp[0] for x in range(0,len(top[0][0:sig_to_plot]))]
                best_vs_list = best_vs_list + [gp[1] for x in range(0,len(top[0][0:sig_to_plot]))]
                best_pvalue_list = best_pvalue_list + top[1][0:sig_to_plot]
        else:
            if cell1_present == []:
                print(gp[1], 'not present in cell matrix')
            else:
                print(gp[0], 'not present in cell matrix')
    fig, axs = plt.subplots(1, len(group_name_list), figsize=(23+len(group_name_list),13), sharex=False, sharey=False)
    axs = axs.ravel()
    color_map = {}

    #plot top significant genes for each group compared to all other groups
    for i, name in enumerate(group_name_list):
        to_plot= barplot_dict[name]
        for v in set(to_plot['Vs']):
            color_map[v] = color_dict_cell[v.split("significance vs ")[-1]][0]
        if color_map != {}:
            g = sns.barplot(x='pvalues', y='genes', hue='Vs', data=to_plot, ax = axs[i], palette=color_map)
            axs[i].set_xscale("log", nonposx='clip')
            bar_list = []
            if plot_pvalue:
                for p in axs[i].patches:
                    height = p.get_height()
                    width = p.get_width()
                    bar_list.append(width)
                max_bar = max(bar_list)
                for p in axs[i].patches:
                    height = p.get_height()
                    width = p.get_width()
                    axs[i].text(max_bar*50,p.get_y()+(height), "{:.2e}".format(width))
            rect = axs[i].patch
            rect.set_facecolor('white')
            #sns.despine(left=True, bottom=True, top=True)
            axs[i].invert_xaxis()
            axs[i].xaxis.set_ticks_position('none')
            axs[i].yaxis.tick_right()
            axs[i].set_title(name)
            axs[i].legend(loc='upper left', bbox_to_anchor=(0.01, 1.11+(0.01*len(group_name_list))), ncol=1, prop={'size':15})
            axs[i].set_xlabel('adjusted p-value')
            for xmaj in axs[i].xaxis.get_majorticklocs():
                axs[i].axvline(x=xmaj,ls='--', lw = 0.5, color='grey', alpha=0.3)
            axs[i].xaxis.grid(True, which="major", linestyle='-')
            plt.subplots_adjust(left=.08, wspace=.3)
    if from_kmeans == '':
        plt.savefig(os.path.join(path_filename,'differential_genes_foldchanges.pdf'), bbox_inches='tight')
        best_gene_df = pd.DataFrame({'GeneID':best_gene_list, 'GroupID':best_gene_groups, 'Vs':best_vs_list, 'adjusted_pvalue':best_pvalue_list})
        best_gene_df.to_csv(os.path.join(path_filename,'Best_Gene_list.txt'), sep = '\t')
    else:
        plt.savefig(os.path.join(path_filename,'kmeans_'+from_kmeans+'_differential_genes_foldchanges.pdf'), bbox_inches='tight')
        best_gene_df = pd.DataFrame({'GeneID':best_gene_list, 'GroupID':best_gene_groups, 'Vs':best_vs_list, 'adjusted_pvalue':best_pvalue_list})
        best_gene_df.to_csv(os.path.join(path_filename,'kmeans_'+from_kmeans+'_Best_Gene_list.txt'), sep = '\t')

'''takes cell groups and creates dictionay 'label_map' that has attached color and marker
if the same cell is assigned to multiple groups it assigns it to the first groupID
'''
def cell_color_map(args, cell_group_filename, cell_list, color_dict):
    cell_groups_df = pd.read_csv(open(os.path.join(os.path.dirname(args.filepath), cell_group_filename),'rU'), sep=None, engine='python')
    cell_list_1 = list(set(cell_groups_df['SampleID'].tolist()))
    group_set = list(set(cell_groups_df['GroupID'].tolist()))
    if len(cell_groups_df['SampleID']) == len(cell_groups_df['GroupID']):
        group_seen = []
        label_map = {}
        cells_seen = []
        for cell, group in list(set(zip(cell_groups_df['SampleID'].tolist(), cell_groups_df['GroupID'].tolist()))):
            if cell not in cells_seen:
                group_count = group_set.index(group)
                label_map[cell] = (color_dict[group][0],color_dict[group][1],group)
                cells_seen.append(cell)
        non_group_cells = [c for c in cell_list if c not in cells_seen]
        if non_group_cells != []:
            from matplotlib import colors
            all_color_list = list(colors.cnames.keys())
            markers = ['o', 'v','D','*','x','h', 's','p','8','^','>','<', 'd','o', 'v','D','*','x','h', 's','p','8','^','>','<', 'd']
            color_list = ['b', 'm', 'r', 'c', 'g', 'orange', 'darkslateblue']+all_color_list
            for cell in non_group_cells:
                label_map[cell] = (color_list[group_count+1],markers[group_count+1],'No_Group')
    else:
        label_map = False
    return cell_list_1, label_map

#takes cell groups and creates dictionay 'label_map' that has attached color and marker
def gene_list_map(args, gene_list_file, gene_list, color_dict, exclude_list = []):
    gene_df1 = pd.read_csv(open(os.path.join(os.path.dirname(args.filepath), gene_list_file),'rU'), sep=None, engine='python')
    if exclude_list != []:
        gene_df1 = gene_df1.ix[~gene_df1['GeneID'].isin(exclude_list)]
    gene_df = gene_df1.copy()
    gene_list_1 = [g for g in list(set(gene_df['GeneID'].tolist())) if g in gene_list]
    if len(gene_df['GeneID']) == len(gene_df['GroupID']):
        gene_label_map = {}
        group_pos = 0
        group_seen = ['xyz' for i in range(len(set(gene_df['GroupID'].tolist())))]
        genes_seen = []
        for gene, group in zip(gene_df['GeneID'].tolist(), gene_df['GroupID'].tolist()):
            #if no GroupIDs are provided replace withe empty strings
            try:
                if math.isnan(float(group)):
                    group = ' '
            except ValueError:
                pass
            if gene not in genes_seen:
                if str(group) in group_seen:
                    pos = group_seen.index(str(group))
                else:
                    group_seen[group_pos] = str(group)
                    pos = group_pos
                    group_pos += 1
                gene_label_map[gene] = (color_dict[str(group)][0],color_dict[str(group)][1],str(group))
                genes_seen.append(gene)
        non_group_genes = [g for g in gene_list_1 if g not in genes_seen]
        if non_group_genes != []:
            all_color_list = list(colors.cnames.keys())
            markers = ['o', 'v','D','*','x','h', 's','p','8','^','>','<', 'd','o', 'v','D','*','x','h', 's','p','8','^','>','<', 'd']
            color_list = ['b', 'm', 'r', 'c', 'g', 'orange', 'darkslateblue']+all_color_list
            for cell in non_group_genes:
                gene_label_map[gene] = (color_list[group_pos+1],markers[group_pos+1],'No_ID')
    else:
        gene_label_map = False
    return gene_list_1, gene_label_map

#this script calls qgraph R package using rpy2, for gene or cell qgraph gene or cell groups must be provided (either or both)
def run_qgraph(data, cell_group_filename, gene_filename, label_map, gene_map, path_filename, gene_or_cell, minimum = 0.25, cut = 0.4, vsize = 1.5, legend = True, borders = False):
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    robjects.r.setwd(os.path.dirname(path_filename))
    qgraph = importr('qgraph')
    psych = importr('psych')
    if gene_or_cell=='cell':
        r_dataframe = pandas2ri.py2ri(data.transpose())
        cell_groups_df = pd.read_csv(open(cell_group_filename,'rU'), sep=None, engine='python')
        cell_list_1 = list(set(cell_groups_df['SampleID'].tolist()))
        group_set = list(set(cell_groups_df['GroupID'].tolist()))
        d = defaultdict(list)
        cell_list_all = data.columns.tolist()
        for i, cell in enumerate(cell_list_all):
            group = label_map[cell][2]
            d[group].append(i+1)

        label_list = robjects.vectors.StrVector(cell_list_all)
    elif gene_or_cell=='gene':
        r_dataframe = pandas2ri.py2ri(data.transpose())
        gene_groups_df = pd.read_csv(open(gene_filename,'rU'), sep=None, engine='python')
        gene_list_1 = list(set(gene_groups_df['GeneID'].tolist()))
        group_set = list(set(gene_groups_df['GroupID'].tolist()))
        d = defaultdict(list)
        d_color = []
        #color_dict2 = {'r':'red', 'm':'magenta', 'b':'blue', 'g':'green', 'c':'cyan'}
        gene_list_all = data.transpose().columns.tolist()
        for i, gene in enumerate(gene_list_all):
            group = gene_map[gene][2]
            color = gene_map[gene][0]
            d[group].append(i+1)
            #d_color.append(color_dict2[color])

        label_list = robjects.vectors.StrVector(gene_list_all)
        #d_color_r = robjects.vectors.StrVector(d_color)
    group_num = len(d)
    from rpy2.robjects.vectors import FloatVector
    for colname in d:
        d[colname] = FloatVector(d[colname])

    # data frame
    from rpy2.robjects.vectors import ListVector
    group_data = ListVector(d)
    pca = psych.principal(robjects.r.cor(r_dataframe),group_num, rotate = "promax")
    robjects.r.setwd(path_filename)
    qpca = qgraph.qgraph(pca, groups = group_data, layout = "circle", rotation = "promax",
    minimum = 0.2, cut = 0.4, vsize = FloatVector([1.5, 15]), labels= label_list, borders = False,
    vTrans = 200,filename='graph_pca_'+gene_or_cell, filetype = "pdf", height = 15, width = 15)
    if gene_or_cell == 'gene':
        Q = qgraph.qgraph(robjects.r.cor(r_dataframe), minimum = 0.25, cut = 0.4, vsize = 1.5, groups = group_data,
        legend = True, borders = False, labels = label_list, filename='graph_'+gene_or_cell, filetype = "pdf", height = 15, width = 15)
    elif gene_or_cell =='cell':
        Q = qgraph.qgraph(robjects.r.cor(r_dataframe), minimum = 0.25, cut = 0.4, vsize = 1.5, groups = group_data,
        legend = True, borders = False, labels = label_list, filename='graph_'+gene_or_cell, filetype = "pdf", height = 15, width = 15)
    Q = qgraph.qgraph(Q, layout = "spring", overlay=True)
    robjects.r.setwd(os.path.dirname(path_filename))


def main(args):
    try:
        new_file = os.path.join(os.path.dirname(args.filepath),args.base_name+'_scicast')
    except AttributeError:
        sys.exit('Please provide a valid path to a file.')
    if args.verbose:
        print('Making new folder for results of SCICAST clustering: '+new_file)
    try:
        os.mkdir(new_file)
    except OSError:
        print(new_file+' already exists. Files will be overwritten.')

    if args.gene_list_filename:
        if os.path.isfile(args.gene_list_filename):
            gene_list_file = args.gene_list_filename
        elif os.path.isfile(os.path.join(os.path.dirname(args.filepath),args.gene_list_filename)):
            gene_list_file = os.path.join(os.path.dirname(args.filepath),args.gene_list_filename)
        else:
            sys.exit('Error: Cannot find gene list file. Please place the gene list file in the same directory or provide a full path.')


    if args.cell_list_filename:
        if os.path.isfile(args.cell_list_filename):
            cell_file = args.cell_list_filename
        elif os.path.isfile(os.path.join(os.path.dirname(args.filepath),args.cell_list_filename)):
            cell_file = os.path.join(os.path.dirname(args.filepath),args.cell_list_filename)
        else:
            sys.exit('Error: Cannot find cell list file. Please place the gene list file in the same directory or provide a full path.')

    try:
        by_cell = pd.DataFrame.from_csv(args.filepath, sep="\t")
    except AttributeError:
        sys.exit('Please provide a valid path to a file.')

    dup_gene_list = [item for item, count in collections.Counter(by_cell.index).items() if count > 1]
    if len(dup_gene_list) >0:
        by_cell.drop(dup_gene_list, inplace=True)
    by_gene = by_cell.transpose()
    #create list of genes
    gene_list_inital = by_cell.index.tolist()
    #create cell list
    cell_list = [x for x in list(by_cell.columns.values)]


    df_by_gene1 = pd.DataFrame(by_gene, columns=gene_list_inital, index=cell_list)
    if args.exclude_genes:
        hu_cc_gene_df = pd.DataFrame.from_csv(os.path.join(os.path.dirname(args.filepath),args.exclude_genes), sep='\t', header=0, index_col=False)
        exclude_list = [ex_gene for ex_gene in hu_cc_gene_df['GeneID'].tolist() if ex_gene in df_by_gene1.index.tolist()]
        try:
            df_by_gene1.drop(exclude_list, axis=1, inplace=True)
        except ValueError:
            print(str(exclude_list)+' already filtered from matrix.')
    df_by_cell1 = df_by_gene1.transpose()
    if args.limit_cells and args.cell_list_filename:
        df_by_cell2, df_by_gene2 = make_new_matrix_cell(df_by_cell1, cell_file)
        log2_expdf_cell, log2_expdf_gene = log2_oulierfilter(df_by_cell2, plot=False, already_log2=args.already_log2)
    else:
        log2_expdf_cell, log2_expdf_gene = log2_oulierfilter(df_by_cell1, plot=False, already_log2=args.already_log2)

    #get kmeans range from input
    kmeans_range1 = [int(x) for x in args.kmeans_cluster_range.split(',')]
    kmeans_range = [kmeans_range1[0],kmeans_range1[1]+1]

    markers = ['o', 'v','D','*','x','h', 's','p','8','^','>','<', 'd','o', 'v','D','*','x','h', 's','p','8','^','>','<', 'd']
    from matplotlib import colors
    all_color_list = list(colors.cnames.keys())
    cell_color_list = ['b', 'm', 'r', 'c', 'g','y','k']+all_color_list
    color_dict_cells= {}
    color_dict_genes= {}
    if args.color_cells:
        color_list1 = args.color_cells.split(' ')
        for i, c in enumerate(color_list1):
            c_pair = c.split(',')
            if len(c_pair) == 2:
                color_dict_cells[c_pair[0]] = [c_pair[1],markers[i]]
            elif len(cc_pair == 3):
                color_dict_cells[c_pair[0]] = [c_pair[1],c_pair[2]]
    elif args.cell_list_filename:
        cell_groups_df = pd.read_csv(open(cell_file,'rU'), sep=None, engine='python')
        group_set = list(set(cell_groups_df['GroupID'].tolist()))
        for g,c,m in zip(group_set, cell_color_list[0:len(group_set)],markers[0:len(group_set)]):
            color_dict_cells[g] =[c,m]
    if args.color_genes =='same':
        color_dict_genes = color_dict_cells
    elif args.color_genes:
        color_list1 = args.color_genes.split(' ')
        for i, c in enumerate(color_list1):
            c_pair = c.split(',')
            if len(c_pair) == 2:
                color_dict_genes[c_pair[0]] = [c_pair[1],markers[i]]
            elif len(c_pair == 3):
                color_dict_genes[c_pair[0]] = [c_pair[1],c_pair[2]]
    elif args.gene_list_filename:
        gene_groups_df = pd.read_csv(open(gene_list_file,'rU'), sep=None, engine='python')
        group_set = list(set(gene_groups_df['GroupID'].tolist()))
        for i, group in enumerate(group_set):
            try:
                if math.isnan(float(group)):
                    group_set[i] = ' '
            except ValueError:
                pass
        for g,c,m in zip(group_set, cell_color_list[0:len(group_set)],markers[0:len(group_set)]):
            color_dict_genes[g] =[c,m]



    if args.gene_list_filename and args.cell_list_filename:
        if args.exclude_genes:
            hu_cc_gene_df = pd.DataFrame.from_csv(os.path.join(os.path.dirname(args.filepath),args.exclude_genes), sep='\t', header=0, index_col=False)
            exclude_list = hu_cc_gene_df['GeneID'].tolist()
            df_by_cell3, df_by_gene3 = make_new_matrix_gene(log2_expdf_gene, gene_list_file, exclude_list=exclude_list)
            df_by_cell, df_by_gene = make_new_matrix_cell(df_by_cell3, cell_file)
            gene_list, gene_color_map = gene_list_map(args, gene_list_file, df_by_gene.columns.values, color_dict_genes, exclude_list=exclude_list)
            cell_list, label_map = cell_color_map(args, cell_file, df_by_cell.columns.values, color_dict_cells)
        else:
            df_by_cell3, df_by_gene3 = make_new_matrix_gene(log2_expdf_gene, gene_list_file)
            df_by_cell, df_by_gene = make_new_matrix_cell(df_by_cell3, cell_file)
            gene_list, gene_color_map = gene_list_map(args, gene_list_file, df_by_gene.columns.values, color_dict_genes)
            cell_list, label_map = cell_color_map(args, cell_file, df_by_cell.columns.values, color_dict_cells)
    elif args.gene_list_filename:
        if args.exclude_genes:
            hu_cc_gene_df = pd.DataFrame.from_csv(args.exclude_genes, sep='\t', header=0, index_col=False)
            exclude_list = hu_cc_gene_df['GeneID'].tolist()
            df_by_cell, df_by_gene = make_new_matrix_gene(log2_expdf_gene, gene_list_file,exclude_list=exclude_list)
            gene_list, gene_color_map = gene_list_map(args, gene_list_file, df_by_gene.columns.values, color_dict_genes,exclude_list=exclude_list)
        else:
            df_by_cell, df_by_gene = make_new_matrix_gene(log2_expdf_gene, gene_list_file)
            gene_list, gene_color_map = gene_list_map(args, gene_list_file, df_by_gene.columns.values, color_dict_genes)
        label_map = False
    elif args.cell_list_filename:
        if args.exclude_genes:
            hu_cc_gene_df = pd.DataFrame.from_csv(os.path.join(os.path.dirname(args.filepath),args.exclude_genes), sep='\t', header=0, index_col=False)
            exclude_list = hu_cc_gene_df['GeneID'].tolist()
            df_by_cell3, df_by_gene3 = make_new_matrix_gene(log2_expdf_gene, log2_expdf_gene.columns.tolist(), exclude_list=exclude_list)
            df_by_cell, df_by_gene = make_new_matrix_cell(df_by_cell3, cell_file)
            cell_list, label_map = cell_color_map(args, cell_file, log2_expdf_cell.columns.values, color_dict_cells)
            gene_color_map = False
        else:
            df_by_cell, df_by_gene = make_new_matrix_cell(log2_expdf_cell, cell_file)
            cell_list, label_map = cell_color_map(args, cell_file, log2_expdf_cell.columns.values, color_dict_cells)
            gene_color_map = False
    else:
        if args.exclude_genes:
            hu_cc_gene_df = pd.DataFrame.from_csv(os.path.join(os.path.dirname(args.filepath),args.exclude_genes), sep='\t', header=0, index_col=False)
            exclude_list = hu_cc_gene_df['GeneID'].tolist()
            df_by_cell, df_by_gene = make_new_matrix_gene(log2_expdf_gene, log2_expdf_gene.columns.tolist(), exclude_list=exclude_list)
            label_map = False
            gene_color_map = False
        else:
            df_by_cell, df_by_gene = log2_expdf_cell, log2_expdf_gene
            label_map = False
            gene_color_map = False




    #run heatmap clustering stability function
    stability_iteration_num = int(args.test_clust_stability)
    if  stability_iteration_num != 0:
        if args.exclude_genes:
            hu_cc_gene_df = pd.DataFrame.from_csv(args.exclude_genes, sep='\t', header=0, index_col=False)
            exclude_list = hu_cc_gene_df['GeneID'].tolist()
            df_by_cell, df_by_gene = make_new_matrix_gene(log2_expdf_gene, log2_expdf_gene.columns.tolist(),exclude_list=exclude_list)
            stability_ratio = clust_stability(args, df_by_gene, new_file, stability_iteration_num ,label_map=label_map)
        else:
            stability_ratio = clust_stability(args, log2_expdf_gene, new_file, stability_iteration_num , label_map=label_map)

    #if there are genes supplied with genes_corr flag process them to a list for correlation search
    if args.genes_corr != '':
        corr_gene_list = args.genes_corr.split(',')
    #otherwise pass an empty list
    else:
        corr_gene_list = []
    #run heatmaps and PCA only if no_heatmaps flag is not provided
    if not args.no_heatmaps:
        if label_map != False and gene_color_map != False: #if both cell and gene lists are provided
            plot_SVD(args,df_by_gene, new_file, num_genes=len(gene_list), title='all_cells_pca', plot=False, label_map=label_map, gene_map=gene_color_map)
            plot_TSNE(args,df_by_gene, new_file, num_genes=len(gene_list), title='all_cells_pca', plot=False, label_map=label_map, gene_map=gene_color_map)
            plot_kmeans(args, df_by_gene, new_file, num_genes=args.gene_number, kmeans_range=kmeans_range, title='all_cells_pca', plot=False, label_map=label_map, run_sig_test=args.kmeans_sig_test)
            top_pca_gene_df, top_pca = return_top_pca_gene(df_by_gene, num_genes=args.gene_number)
            if not args.no_corr:
                if corr_gene_list != []:
                    top_genes_search = top_pca[0:50]
                    corr_plot(top_genes_search, df_by_gene, new_file, num_to_plot=3, gene_corr_list= corr_gene_list, title = 'All_cells', label_map=label_map)
                else:
                    top_genes_search = top_pca[0:50]
                    corr_plot(top_genes_search, df_by_gene, new_file, num_to_plot=3, title = 'All_cells', label_map=label_map)
            cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, gene_list, df_by_gene, new_file, num_to_plot=len(gene_list), label_map=label_map, gene_map=gene_color_map)
            cc = make_tree_json(cell_linkage, plotted_df_by_gene, new_file)
            make_subclusters(args, cc, df_by_cell, df_by_cell, gene_corr_list=corr_gene_list ,path_filename=new_file, base_name=args.base_name, group_colors=True, label_map=label_map, gene_map=gene_color_map, cluster_size=args.depth_of_clustering)
            if args.qgraph_plot == 'both':
                run_qgraph(df_by_cell, cell_file, gene_list_file, label_map, gene_color_map, new_file, gene_or_cell='gene')
                run_qgraph(df_by_cell, cell_file, gene_list_file, label_map, gene_color_map, new_file, gene_or_cell='cell')
            elif args.qgraph_plot == 'gene':
                run_qgraph(df_by_cell, cell_file, gene_list_file, label_map, gene_color_map, new_file, gene_or_cell='gene')
            elif args.qgraph_plot == 'cell':
                run_qgraph(df_by_cell, cell_file, gene_list_file, label_map, gene_color_map, new_file, gene_or_cell='cell')
            if args.all_sig:
                find_twobytwo(cc, df_by_cell, df_by_cell, new_file, cluster_size=args.depth_of_clustering)
        elif label_map != False: #if only cell list is provided
            plot_SVD(args,df_by_gene, new_file, num_genes=int(args.gene_number), title='all_cells_pca', plot=False, label_map=label_map)
            plot_TSNE(args,df_by_gene, new_file, num_genes=int(args.gene_number), title='all_cells_pca', plot=False, label_map=label_map)
            plot_kmeans(args, df_by_gene, new_file, num_genes=args.gene_number, kmeans_range=kmeans_range, title='all_cells_pca', plot=False, label_map=label_map, run_sig_test=args.kmeans_sig_test)
            top_pca_by_gene, top_pca = return_top_pca_gene(df_by_gene, num_genes=args.gene_number)
            top_pca_by_cell = top_pca_by_gene.transpose()
            if not args.no_corr:
                if corr_gene_list != []:
                    top_genes_search = top_pca[0:50]
                    corr_plot(top_genes_search, df_by_gene, new_file, num_to_plot=3, gene_corr_list= corr_gene_list, title = 'All_cells', label_map=label_map)
                else:
                    top_genes_search = top_pca[0:50]
                    corr_plot(top_genes_search, df_by_gene, new_file, num_to_plot=3, title = 'All_cells', label_map=label_map)
            cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, top_pca, top_pca_by_gene, new_file, num_to_plot=args.gene_number, label_map=label_map)
            cc = make_tree_json(cell_linkage, plotted_df_by_gene, new_file)
            make_subclusters(args, cc, top_pca_by_cell, df_by_cell, new_file, base_name=args.base_name, gene_corr_list=corr_gene_list, group_colors=True, label_map=label_map, cluster_size=args.depth_of_clustering)
            if args.all_sig:
                find_twobytwo(cc, df_by_cell, log2_expdf_cell, new_file, cluster_size=args.gene_number)
        elif gene_color_map != False: #if only gene list is provided
            plot_SVD(args,df_by_gene, new_file, num_genes=int(args.gene_number), title='all_cells_pca', plot=False, label_map=label_map, gene_map=gene_color_map)
            plot_TSNE(args,df_by_gene, new_file, num_genes=int(args.gene_number), title='all_cells_pca', plot=False, label_map=label_map, gene_map=gene_color_map)
            plot_kmeans(args,df_by_gene, new_file, num_genes=args.gene_number, kmeans_range=kmeans_range, title='all_cells_pca', plot=False, label_map=label_map, run_sig_test=args.kmeans_sig_test)
            top_pca_by_gene, top_pca = return_top_pca_gene(df_by_gene, num_genes=args.gene_number)
            top_pca_by_cell = top_pca_by_gene.transpose()
            if not args.no_corr:
                if corr_gene_list != []:
                    top_genes_search = top_pca[0:50]
                    corr_plot(top_genes_search, df_by_gene, new_file, num_to_plot=3, gene_corr_list= corr_gene_list, title = 'All_cells', label_map=label_map, gene_map=gene_color_map)
                else:
                    top_genes_search = top_pca[0:50]
                    corr_plot(top_genes_search, df_by_gene, new_file, num_to_plot=3, title = 'All_cells', label_map=label_map)
            cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, top_pca, top_pca_by_gene, new_file, num_to_plot=args.gene_number, label_map=label_map, gene_map=gene_color_map)
            cc = make_tree_json(cell_linkage, plotted_df_by_gene, new_file)
            make_subclusters(args, cc, top_pca_by_cell, df_by_cell, new_file, base_name=args.base_name, gene_corr_list=corr_gene_list, group_colors=True, label_map=label_map, gene_map=gene_color_map, cluster_size=args.depth_of_clustering)
            if args.all_sig:
                find_twobytwo(cc, df_by_cell, log2_expdf_cell, new_file, cluster_size=args.depth_of_clustering)
        else: #if neither cell or gene lists are provided
            label_map=False
            group_colors = False
            plot_SVD(args, log2_expdf_gene, new_file, num_genes=args.gene_number, title='all_cells_pca', plot=False, label_map=label_map)
            plot_TSNE(args, log2_expdf_gene, new_file, num_genes=args.gene_number, title='all_cells_pca', plot=False, label_map=label_map)
            plot_kmeans(args,log2_expdf_gene, new_file, num_genes=args.gene_number, kmeans_range=kmeans_range, title='all_cells_pca', plot=False, label_map=label_map, run_sig_test=args.kmeans_sig_test)
            top_pca_by_gene, top_pca = return_top_pca_gene(df_by_gene, num_genes=args.gene_number)
            top_pca_by_cell = top_pca_by_gene.transpose()
            if not args.no_corr:
                if corr_gene_list != []:
                    top_genes_search = top_pca[0:50]
                    corr_plot(top_genes_search, df_by_gene, new_file, num_to_plot=3, gene_corr_list= corr_gene_list, title = 'All_cells', label_map=label_map)
                else:
                    top_genes_search = top_pca[0:50]
                    corr_plot(top_genes_search, df_by_gene, new_file, num_to_plot=3, title = 'All_cells', label_map=label_map)
            cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, top_pca, top_pca_by_gene, new_file, num_to_plot=args.gene_number, label_map=label_map)
            cc = make_tree_json(cell_linkage, plotted_df_by_gene, new_file)
            make_subclusters(args,cc, top_pca_by_cell, df_by_cell, new_file, base_name=args.base_name, gene_corr_list=corr_gene_list, cluster_size=args.depth_of_clustering)
            if args.all_sig:
                find_twobytwo(cc, top_pca_by_cell, log2_expdf_cell, new_file, cluster_size=args.depth_of_clustering)
    if args.group_sig_test and args.cell_list_filename:
        multi_group_sig(args, log2_expdf_cell, cell_file, new_file, color_dict_cells)


if __name__ == "__main__":
    from scicast_argparse import check_gui_parser
    try:
        gui = check_gui_parser()
        run_gui = gui.gui_true
    except:
        run_gui = False
    if run_gui:
        from tkinter_scicast import Window
        scicast_window = Window()
        scicast_window.mainloop()
        from sci_load import Sci_load
        scil = Sci_load()
        try:
            opts_all = scil.load_options(all_options_dict = scicast_window.all_dict)
        except AttributeError:
            sys.exit('Please provide (at a minimum) a valid path to a file and click Run scicast.')
        main(opts_all)

    else:
        from scicast_argparse import get_parser
        opts_all = get_parser()
        main(opts_all)
