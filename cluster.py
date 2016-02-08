import pickle as pickle
import numpy as np
import pandas as pd
import os
import sys
from subprocess import call
import matplotlib
matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
import scipy
import json
from sklearn.decomposition import PCA as skPCA
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram, set_link_color_palette, to_tree, inconsistent
import seaborn as sns
from matplotlib.colors import rgb2hex, colorConverter
from pprint import pprint
import difflib
from operator import itemgetter
import itertools
from functools import reduce
import matplotlib.ticker as ticker



def make_new_matrix_gene(org_matrix_by_gene, gene_list_file):
    split_on='_'
    gene_df = pd.read_csv(gene_list_file, delimiter= '\t')
    gene_list = list(set(gene_df['GeneID'].tolist()))
    group_list = gene_df['GroupID'].tolist()
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

def make_new_matrix_cell(org_matrix_by_cell, cell_list_file):
    cell_df = pd.read_csv(cell_list_file, delimiter= '\t')
    cell_list_new = [cell.strip('\n') for cell in cell_df['SampleID'].tolist()]
    new_cmatrix_df = org_matrix_by_cell[cell_list_new]
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



def run_cluster(by_gene_matrix, path_filename):
    cell_list = [x for x in list(by_gene_matrix.index.values)]
    cell_dist = pdist(np.array(by_gene_matrix), metric='euclidean')
    row_dist = pd.DataFrame(squareform(cell_dist), columns=cell_list, index=cell_list)
    row_clusters = linkage(cell_dist, metric=metric, method='average')
    link_mat = pd.DataFrame(row_clusters,
                 columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
                 index=['cluster %d' %(i+1) for i in range(row_clusters.shape[0])])
    row_dendr = dendrogram(row_clusters, labels=cell_list, leaf_rotation=90, leaf_font_size=8)

    plt.savefig(os.path.join(path_filename,'dendrogram_gene.png'))
    plt.clf()
    return cell_dist, row_dist, row_clusters, link_mat, row_dendr

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

#Makes labeled json tree for visulaization in d3
def make_tree_json(row_clusters, df_by_gene, path_filename):
    T= to_tree(row_clusters)

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


def plot_PCA(df_by_gene, path_filename, num_genes=100, gene_list_filter=False, title='', plot=False, label_map=False, gene_map = False):
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
    if args.verbose:
        print(top_pca_list[0:num_genes], 'top_pca_list')
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
            annotate = False
            X = [x for x in top_cell_trans[:, 0]]
            Y = [y for y in top_cell_trans[:, 1]]
            labels = [label_map[cell][2] for cell in top_by_cell.columns.tolist()]
            markers = [label_map[cell][1] for cell in top_by_cell.columns.tolist()]
            colors = [label_map[cell][0] for cell in top_by_cell.columns.tolist()]
            label_done = []
            for X_pos, Y_pos, m, color, l in zip(X, Y, markers, colors, labels):
                if l in label_done:
                    lab = ''
                else:
                    lab= l
                    label_done.append(l)
                ax_cell.scatter(X_pos, Y_pos, marker=m, c=color, label=lab, s=30)

        else:
            ax_cell.scatter(top_cell_trans[:, 0], top_cell_trans[:, 1], alpha=0.75)
            annotate = True
        ax_cell.set_xlim([min(top_cell_trans[:, 0])-1, max(top_cell_trans[:, 0]+1)])
        ax_cell.set_ylim([min(top_cell_trans[:, 1])-1, max(top_cell_trans[:, 1]+2)])
        ax_cell.set_title(title+'_cell')
        ax_cell.legend(loc='best', ncol=1, prop={'size':12}, markerscale=1.5, frameon=True)
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
            for X_pos, Y_pos, m, color, l in zip(X, Y, markers, colors, labels):
                ax_gene.scatter(X_pos, Y_pos, marker=m, c=color, label = l, s=30)
        else:
            ax_gene.scatter(top_gene_trans[:, 0], top_gene_trans[:, 1], alpha=0.75)
        ax_gene.set_xlim([min(top_gene_trans[:, 0])-1, max(top_gene_trans[:, 0])+1])
        ax_gene.set_ylim([min(top_gene_trans[:, 1])-1, max(top_gene_trans[:, 1])+2])
        ax_gene.set_title(title+'_gene')
        ax_gene.set_xlabel('PC1')
        ax_gene.set_ylabel('PC2')
        for label, x, y in zip(top_by_gene.columns, top_gene_trans[:, 0], top_gene_trans[:, 1]):
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

def clust_heatmap(gene_list, df_by_gene, path_filename, num_to_plot, title='', plot=False, label_map=False, gene_map=False):
    if num_to_plot >175:
        sns.set(context= 'poster', font_scale = 0.65/(num_to_plot/100))
    else:
        sns.set(context= 'poster', font_scale = .80, font ='Verdana')

    if len(str(args.z_direction)) > 1:
        z_choice = str(args.z_direction)
        if z_choice != 'None':
            sys.exit('Please enter a valid option (0, 1, or None) for z_direction')
    else:
        z_choice = int(args.z_direction)
        if z_choice != 0 and z_choice != 1:
            sys.exit('Please enter a valid option (0, 1, or None) for z_direction')
    cmap = sns.diverging_palette(255, 10, s=99, sep=1, as_cmap=True)
    cell_list = df_by_gene.index.tolist()
    cg = sns.clustermap(df_by_gene[gene_list[0:num_to_plot]].transpose(), metric=args.metric, method=args.method, z_score=z_choice, figsize=(30, 30), cmap =cmap)
    col_order = cg.dendrogram_col.reordered_ind
    row_order = cg.dendrogram_row.reordered_ind
    cg.ax_heatmap.set_title(title)
    if label_map:
        Xlabs = [cell_list[i] for i in col_order]
        Xcolors = [label_map[cell][0] for cell in Xlabs]
        for xtick, xcolor in zip(cg.ax_heatmap.get_xticklabels(), Xcolors):
            xtick.set_color(xcolor)
            xtick.set_rotation(270)
    else:
        for xtick in cg.ax_heatmap.get_xticklabels():
            xtick.set_rotation(270)
    if gene_map:
        Ylabs = [gene_list[i] for i in row_order]
        Ycolors = [gene_map[gene][0] for gene in Ylabs]
        for ytick, ycolor in zip(cg.ax_heatmap.get_yticklabels(), list(reversed(Ycolors))):
            ytick.set_color(ycolor)
            ytick.set_rotation(0)
    else:
        for ytick in cg.ax_heatmap.get_yticklabels():
            ytick.set_rotation(0)
    if plot:
        plt.show()
    cell_linkage = cg.dendrogram_col.linkage

    link_mat = pd.DataFrame(cell_linkage,
                columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
                index=['cluster %d' %(i+1) for i in range(cell_linkage.shape[0])])
    if title != '':
        save_name = '_'.join(title.split(' ')[0:2])
        cg.savefig(os.path.join(path_filename, save_name+'_heatmap.pdf'), bbox_inches='tight')
    else:
        cg.savefig(os.path.join(path_filename,'Group0_Heatmap_all_cells.pdf'), bbox_inches='tight')
    plt.close('all')
    return cell_linkage, df_by_gene[gene_list[0:num_to_plot]], col_order

def make_subclusters(cc, log2_expdf_cell, log2_expdf_cell_full, path_filename, base_name, gene_corr_list, label_map=False, gene_map=False, cluster_size=20, group_colors=False):
    parent = cc[0][1]
    p_num = cc[0][0]
    l_nums = [x[0] for x in cc]
    c_lists = [c[1] for c in cc]
    group_ID = 0

    for num_members, cell_list in zip(l_nums, c_lists):
        if num_members < p_num and num_members >= cluster_size:
            group_ID+=1
            current_title = 'Group_'+str(group_ID)+'_with_'+str(num_members)+'_cells'
            cell_subset = log2_expdf_cell[list(set(cell_list))]
            full_cell_subset = log2_expdf_cell_full[list(set(cell_list))]
            gene_subset = cell_subset.transpose()
            norm_df_cell1 = np.exp2(full_cell_subset)
            norm_df_cell = norm_df_cell1 -1
            norm_df_cell.to_csv(os.path.join(path_filename, base_name+'_'+current_title+'_matrix.txt'), sep = '\t', index_col=0)
            if gene_map:
                top_pca = plot_PCA(gene_subset, path_filename, num_genes=len(gene_subset.columns.tolist()), title=current_title, plot=False, label_map=label_map, gene_map=gene_map)
            else:
                top_pca = plot_PCA(gene_subset, path_filename, num_genes=int(args.gene_number), title=current_title, plot=False, label_map=label_map)
            if len(top_pca)<args.gene_number:
                plot_num = len(top_pca)
            else:
                plot_num = args.gene_number
            if top_pca != []:
                top_pca_by_gene = gene_subset[top_pca]
                top_pca_by_cell = top_pca_by_gene.transpose()
                if gene_corr_list != []:
                    top_genes_search = gene_corr_list+top_pca[0:30]
                    corr_plot(top_genes_search, gene_subset, path_filename, num_to_plot=3, title = current_title, label_map=label_map)
                else:
                    top_genes_search = top_pca[0:30]
                    corr_plot(top_genes_search, gene_subset, path_filename, num_to_plot=3, title = current_title, label_map=label_map)
                if gene_map:
                    cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(top_pca, top_pca_by_gene, path_filename, num_to_plot=plot_num, title=current_title, plot=False, label_map=label_map, gene_map = gene_map)
                else:
                    cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(top_pca, top_pca_by_gene, path_filename,num_to_plot=plot_num, title=current_title, plot=False, label_map=label_map)
                plt.close('all')
            else:
                pass

def clust_stability(log2_expdf_gene, iterations, path_filename):
    sns.set(context='poster', font_scale = 1)
    sns.set_palette("RdBu_r")
    stability_ratio = []
    total_genes = len(log2_expdf_gene.columns.tolist())
    end_num = 1000
    iter_list = range(100,int(round(end_num)),int(round(end_num/iterations)))
    for gene_number in iter_list:
        title= str(gene_number)+' genes plot.'
        top_pca = plot_PCA(log2_expdf_gene, num_genes=gene_number, title=title)
        top_pca_by_gene = log2_expdf_gene[top_pca]
        top_pca_by_cell = top_pca_by_gene.transpose()
        cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(top_pca, top_pca_by_gene, num_to_plot=gene_number, title=title)
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
def run_corr(df_by_gene, title, path_filename, method_name='pearson', sig_threshold= 0.5, run_new=True, min_period=3):
    if run_new:
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

        with open(os.path.join(path_filename,'gene_correlations_sig_neg_'+method_name+'.p'), 'wb') as fp:
            pickle.dump(corr_by_gene_neg, fp)
        with open(os.path.join(path_filename,'gene_correlations_sig_pos_'+method_name+'.p'), 'wb') as fp0:
            pickle.dump(corr_by_gene_pos, fp0)
        with open(os.path.join(path_filename,'by_gene_corr.p'), 'wb') as fp1:
            pickle.dump(corr_by_gene, fp1)
        with open(os.path.join(path_filename,'by_cell_corr.p'), 'wb') as fp2:
            pickle.dump(corr_by_cell, fp2)
    else:
        corr_by_g_pos =  open(os.path.join(path_filename,'gene_correlations_sig_pos_'+method_name+'.p'), 'rb')
        corr_by_g_neg =  open(os.path.join(path_filename,'gene_correlations_sig_neg_'+method_name+'.p'), 'rb')
        corr_by_gene_pos = pickle.load(corr_by_g_pos)
        corr_by_gene_neg = pickle.load(corr_by_g_neg)

    cor_pos_df = pd.DataFrame(corr_by_gene_pos)
    cor_neg_df = pd.DataFrame(corr_by_gene_neg)
    sig_corr = cor_pos_df.append(cor_neg_df)
    sig_corrs = pd.DataFrame(sig_corr[0], columns=["corr"])

    if run_new:
        sig_corrs.to_csv(os.path.join(path_filename, title+'_counts_corr_sig_'+method_name+'.txt'), sep = '\t')
    return sig_corrs
def find_top_corrs(terms_to_search, sig_corrs, num_to_return):
    all_corrs_list = []
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

    return all_corrs_list[0:num_to_return]


#corr_plot finds and plots all correlated genes, log turns on log scale, sort plots the genes in the rank order of the gene searched
def corr_plot(terms_to_search, df_by_gene, path_filename, title, num_to_plot, label_map=False, log=False, sort=True, sig_threshold=0.5):
    sig_corrs = run_corr(df_by_gene, title, path_filename, sig_threshold=sig_threshold)
    corr_list = find_top_corrs(terms_to_search, sig_corrs, num_to_plot)
    for corr_tup in corr_list:
        term_to_search = corr_tup[0][0]
        corr_tup.sort(key=itemgetter(1), reverse=True)
        corr_df = pd.DataFrame(corr_tup, columns=['GeneID', 'Correlation'])
        corr_df.to_csv(os.path.join(path_filename, title+'_Corr_w_'+term_to_search+'_list.txt'), sep = '\t', index=False)
        if args.verbose:
            print('Correlations:')
            for c in corr_tup:
                print(c)
        to_plot = [x[0] for x in corr_tup]
        sns.set_palette(sns.cubehelix_palette(len(to_plot), start=1, rot=-.9, reverse=True))
        try:
            sorted_df = df_by_gene.sort_values(by=[term_to_search])
            log2_df = np.log2(df_by_gene[to_plot])
            sorted_log2_df=np.log2(sorted_df[to_plot])
            ylabel='FPKM (log2)'
            if sort and log:
                ax = sorted_log2_df.plot()
                xlabels = sorted_log2_df[to_plot].index.values
            elif sort:
                ax =sorted_df[to_plot].plot()
                xlabels = sorted_df[to_plot].index.values
            elif log:
                ax = log2_df.plot()
                ylabel= 'log2 FPKM'
                xlabels = log2_df.index.values
            else:
                ax = df_by_gene[to_plot].plot()
                xlabels = df_by_gene[to_plot].index.values
            ax.set_xlabel('Cell #')
            ax.set_ylabel(ylabel)
            ax.set_title('Correlates with '+term_to_search, loc='right')
            ax.xaxis.set_minor_locator(LinearLocator(numticks=len(xlabels)))
            if label_map:
                ax.set_xticklabels(xlabels, minor=True, rotation='vertical', fontsize=6)
                Xcolors = [label_map[cell][0] for cell in xlabels]
                for xtick, xcolor in zip(ax.get_xticklabels(which='minor'), Xcolors):
                    xtick.set_color(xcolor)
                    xtick.set_rotation(90)
            else:
                ax.set_xticklabels(xlabels, minor=True, rotation='vertical', fontsize=6)
            ax.set_ylim([0, df_by_gene[to_plot].values.max()])
            ax.xaxis.set_major_formatter(ticker.NullFormatter())
            ax.tick_params(axis='x', which ='minor', labelsize=10)
            if len(corr_tup) > 15:
                l_labels = [str(x[0])+' '+"%.2f" % x[1] for x in corr_tup]
                ax.legend(l_labels, loc='upper left', bbox_to_anchor=(0.01, 1.05), ncol=6, prop={'size':6})
            else:
                l_labels = [str(x[0])+' '+"%.2f" % x[1] for x in corr_tup]
                ax.legend(l_labels, loc='upper left', bbox_to_anchor=(0.01, 1.05), ncol=4, prop={'size':8})
            fig = plt.gcf()
            fig.subplots_adjust(bottom=0.08, top=0.95, right=0.98, left=0.03)
            plt.savefig(os.path.join(path_filename, title+'_corr_with_'+term_to_search+'.pdf'), bbox_inches='tight')
            plt.close('all')
        except KeyError:
            if args.verbose:
                print(term_to_search+' not in this matrix')
            pass


def multi_group_sig(full_by_cell_df, cell_group_filename, path_filename):
    cell_groups_df = pd.read_csv(cell_group_filename, delimiter= '\t')
    group_name_list = list(set(cell_groups_df['GroupID']))
    group_pairs = list(set(itertools.permutations(group_name_list,2)))
    gene_list = full_by_cell_df.index.tolist()
    cell_group_ident_0 = zip(cell_groups_df['SampleID'],cell_groups_df['GroupID'])
    cell_group_ident= [c for c in cell_group_ident_0]
    for gp in group_pairs:
        g_pvalue_dict = {}
        index_list = []
        sig_gene_list = []
        cell_list1 = [c[0] for c in cell_group_ident if c[1] == gp[0]]
        cell_list2 = [c[0] for c in cell_group_ident if c[1] == gp[1]]
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
        by_gene_df = full_by_cell_df.transpose()
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
        sig_df.to_csv(os.path.join(path_filename,'sig_'+gp[0]+'_'+gp[1]+'_pvalues.txt'), sep = '\t')
        cell_names_df.to_csv(os.path.join(path_filename,'sig_'+gp[0]+'_'+gp[1]+'_cells.txt'), sep = '\t')

def cell_color_map(cell_group_filename, cell_list, color_list, markers):
    cell_groups_df = pd.read_csv(cell_group_filename, delimiter= '\t')
    cell_list_1 = list(set(cell_groups_df['SampleID'].tolist()))
    if len(cell_groups_df['SampleID']) == len(cell_groups_df['GroupID']):
        group_count = -1
        group_seen = []
        label_map = {}
        cells_seen = []
        for cell, group in zip(cell_groups_df['SampleID'].tolist(), cell_groups_df['GroupID'].tolist()):
            if cell not in cells_seen:
                if group not in group_seen:
                    group_count += 1
                    group_seen.append(group)
                label_map[cell] = (color_list[group_count],markers[group_count],group)
                cells_seen.append(cell)
        non_group_cells = [c for c in cell_list if c not in cells_seen]
        if non_group_cells != []:
            for cell in non_group_cells:
                label_map[cell] = (color_list[group_count+1],markers[group_count+1],'No_Group')
    else:
        label_map = False
    return cell_list_1, label_map

def gene_list_map(gene_list_file, gene_list, color_list, markers):
    gene_df = pd.read_csv(os.path.join(os.path.dirname(args.filepath), gene_list_file), delimiter= '\t')
    gene_list_1 = list(set(gene_df['GeneID'].tolist()))
    if len(gene_df['GeneID']) == len(gene_df['GroupID']):
        gene_label_map = {}
        group_pos = 0
        group_seen = ['xyz' for i in range(len(set(gene_df['GroupID'].tolist())))]
        genes_seen = []
        for gene, group in zip(gene_df['GeneID'].tolist(), gene_df['GroupID'].tolist()):
            if gene not in genes_seen:
                if str(group) in group_seen:
                    pos = group_seen.index(str(group))
                else:
                    group_seen[group_pos] = str(group)
                    pos = group_pos
                    group_pos += 1
                gene_label_map[gene] = (color_list[pos],markers[pos],group)
                genes_seen.append(gene)
        non_group_genes = [g for g in gene_list if g not in genes_seen]
        if non_group_genes != []:
            for cell in non_group_genes:
                gene_label_map[gene] = (color_list[group_pos+1],markers[group_pos+1],'No_ID')
    else:
        gene_label_map = False
    return gene_list_1, gene_label_map

def main(args):
    new_file = os.path.join(os.path.dirname(args.filepath),args.base_name+'_subgroups')
    if args.verbose:
        print('Making new folder for results of SCICAST clustering: '+new_file)
    call('mkdir -p '+new_file, shell=True)

    by_cell = pd.DataFrame.from_csv(args.filepath, sep='\t')

    by_gene = by_cell.transpose()
    #create list of genes
    gene_list_inital = by_cell.index.tolist()
    #create cell list
    cell_list = [x for x in list(by_cell.columns.values)]


    df_by_gene1 = pd.DataFrame(by_gene, columns=gene_list_inital, index=cell_list)
    df_by_cell1 = pd.DataFrame(by_cell, columns=cell_list, index=gene_list_inital)
    log2_expdf_cell, log2_expdf_gene = log2_oulierfilter(df_by_cell1, plot=False)
    color_list = ['r', 'm', 'c', 'b', 'g', 'orange', 'darkslateblue']
    markers = ['o', 'v','D','*','x','h', 's']


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

    if args.gene_list_filename and args.cell_list_filename:
        df_by_cell3, df_by_gene3 = make_new_matrix_gene(log2_expdf_gene, gene_list_file)
        df_by_cell, df_by_gene = make_new_matrix_cell(df_by_cell3, cell_file)
        gene_list, gene_color_map = gene_list_map(gene_list_file, df_by_gene.columns.values, color_list, markers)
        cell_list, label_map = cell_color_map(cell_file, df_by_cell.columns.values, color_list, markers)
    elif args.gene_list_filename:
        df_by_cell, df_by_gene = make_new_matrix_gene(log2_expdf_gene, gene_list_file)
        gene_list, gene_color_map = gene_list_map(gene_list_file, df_by_gene.columns.values, color_list, markers)
    elif args.cell_list_filename:
        df_by_cell, df_by_gene = make_new_matrix_cell(log2_expdf_cell, cell_file)
        cell_list, label_map = cell_color_map(cell_file, df_by_cell.columns.values, color_list, markers)
        gene_color_map = False
    else:
        df_by_cell, df_by_gene = df_by_cell1, df_by_gene1
        label_map = False
        gene_color_map = False



    if args.test_clust_stability != 0:
        stability_ratio = clust_stability(log2_expdf_gene, iterations = test_clust_stability)
    if args.genes_corr != '':
        corr_gene_list = args.genes_corr.split(',')
    else:
        corr_gene_list = []

    if label_map != False:
        if gene_color_map != False:
            plot_PCA(df_by_gene, new_file, num_genes=len(gene_list), title='all_cells_pca', plot=False, label_map=label_map, gene_map=gene_color_map)
            cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(gene_list, df_by_gene, new_file, num_to_plot=len(gene_list), label_map=label_map, gene_map=gene_color_map)
            cc = make_tree_json(cell_linkage, plotted_df_by_gene, new_file)
            make_subclusters(cc, df_by_cell, log2_expdf_cell, gene_corr_list=corr_gene_list ,path_filename=new_file, base_name=args.base_name, group_colors=True, label_map=label_map, gene_map=gene_color_map, cluster_size=args.depth_of_clustering)
        else:
            top_pca = plot_PCA(df_by_gene, new_file, num_genes=int(args.gene_number), title='all_cells_pca', plot=False, label_map=label_map)
            top_pca_by_gene = log2_expdf_gene[top_pca]
            top_pca_by_cell = top_pca_by_gene.transpose()
            cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(top_pca, top_pca_by_gene, new_file, num_to_plot=args.gene_number, label_map=label_map)
            cc = make_tree_json(cell_linkage, plotted_df_by_gene, new_file)
            make_subclusters(cc, top_pca_by_cell, log2_expdf_cell, new_file, base_name=args.base_name, gene_corr_list=corr_gene_list, group_colors=True, label_map=label_map, cluster_size=args.depth_of_clustering)
        if args.all_sig:
            find_twobytwo(cc, df_by_cell, log2_expdf_cell, new_file, cluster_size=args.depth_of_clustering)
    else:
        label_map=False
        group_colors = False
        top_pca = plot_PCA(log2_expdf_gene, new_file, num_genes=args.gene_number, title='all_cells_pca', plot=False, label_map=label_map)
        top_pca_by_gene = log2_expdf_gene[top_pca]
        top_pca_by_cell = top_pca_by_gene.transpose()
        cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(top_pca, top_pca_by_gene, new_file, num_to_plot=args.gene_number, label_map=label_map)
        cc = make_tree_json(cell_linkage, plotted_df_by_gene, new_file)
        make_subclusters(cc, top_pca_by_cell, log2_expdf_cell, new_file, base_name=args.base_name, gene_corr_list=corr_gene_list, cluster_size=args.depth_of_clustering)
        if args.all_sig:
            find_twobytwo(cc, top_pca_by_cell, log2_expdf_cell, new_file, cluster_size=args.depth_of_clustering)
    if args.group_sig_test and args.cell_list_filename:
        multi_group_sig(log2_expdf_cell, cell_file, new_file)

    #cell_dist, row_dist, row_clusters, link_mat, row_dendr = run_cluster(top_pca_by_gene)

    #augmented_dendrogram(row_clusters, labels=top_pca_by_cell.columns.tolist(), leaf_rotation=90, leaf_font_size=8)

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
