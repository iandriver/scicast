import numpy as np
import pandas as pd
import os
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy
import json

from scipy.cluster import hierarchy
import seaborn as sns
from functools import reduce
import matplotlib.patches as patches
import warnings





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
def find_twobytwo(cc, args, matrix_data):
    if not args.verbose:
        warnings.filterwarnings("ignore", category=RuntimeWarning)
    gene_list = matrix_data.gene_list
    by_gene_df = matrix_data.log2_df_cell.transpose()
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
            if not overlap and len(c)>=args.depth_of_clustering and len(c2)>=args.depth_of_clustering:
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
        df_by_cell_1 = matrix_data.log2_df_cell[cell_list1]
        df_by_cell_2 = matrix_data.log2_df_cell[cell_list2]
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
        sig_df.to_csv(os.path.join(matrix_data.new_filepath,'sig_'+v+'_pvalues.txt'), sep = '\t')
        cell_names_df.to_csv(os.path.join(matrix_data.new_filepath,'sig_'+v+'_cells.txt'), sep = '\t')







def clust_heatmap(args, matrix_data, title= '', matrix_subset = None, fontsize=18, stablity_plot_num = 0, kmeans_color_map= {}, plot=True):
    try:
        from .dim_reduction import return_top_pca_gene
    except (SystemError, ValueError, ImportError):
        from dim_reduction import return_top_pca_gene
    if isinstance(matrix_subset,pd.DataFrame):
        cell_list = matrix_subset.index.tolist()
        cell_num =len(cell_list)
        df_by_gene = matrix_subset

        num_to_plot = min(args.gene_number, len(matrix_subset.columns.tolist()))
    else:
        if isinstance(matrix_data.log2_df_cell_gene_restricted,pd.DataFrame):
            cell_list = matrix_data.cell_list
            cell_num =len(cell_list)
            df_by_gene = matrix_data.log2_df_cell_gene_restricted.transpose()
            num_to_plot = min(args.gene_number, len(matrix_data.short_gene_list))
        else:
            cell_list = matrix_data.cell_list
            cell_num =len(cell_list)
            df_by_gene = matrix_data.log2_df_cell.transpose()
            num_to_plot = min(args.gene_number, len(matrix_data.gene_list))

    if kmeans_color_map:
        cell_color_map = kmeans_color_map
        first_group_name = 'Kmeans_groups'
        groups_list = ['Kmeans_groups']
    elif matrix_data.cell_label_map:
        cell_color_map = matrix_data.cell_label_map
        first_group_name = matrix_data.cell_group_names[0]
        groups_list = matrix_data.cell_group_names
    else:
        cell_color_map = {}

    if stablity_plot_num != 0:
        num_to_plot = stablity_plot_num
    top_pca_by_gene, top_pca_gene_list = return_top_pca_gene(args, df_by_gene.transpose(), user_num_genes=num_to_plot)
    longest_side = max(num_to_plot,cell_num*2)
    if longest_side > 500:
        fontsize=10
    elif longest_side >250:
        fontsize=14
    if longest_side == num_to_plot:
        sns.set(context= 'poster', font_scale = min(.4*(num_to_plot/100),3))
        width_heatmap = min(28+round(cell_num/50),42+round(cell_num/40))
        len_heatmap = min(43+round(num_to_plot/8),68+round(num_to_plot/20))
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
    gene_list = top_pca_gene_list[0:num_to_plot]
    cluster_df = df_by_gene[gene_list].transpose()
    cluster_df[abs(cluster_df)<3e-12] = 0.0
    try:
        cg = sns.clustermap(cluster_df, method=args.method, metric=args.metric, z_score=z_choice, figsize=(width_heatmap, len_heatmap), cmap =cmap)
        col_order = cg.dendrogram_col.reordered_ind
        row_order = cg.dendrogram_row.reordered_ind

        if cell_color_map and matrix_data.gene_label_map and stablity_plot_num == 0:
            Xlabs = [cell_list[i] for i in col_order]
            Xcolors = {}
            for group_index , group_name in enumerate(groups_list):
                Xcolors[group_name]=[cell_color_map[cell][group_index][0] for cell in Xlabs]
            col_colors = pd.DataFrame(Xcolors,index=Xlabs)
            Xgroup_labels = [cell_color_map[cell][0][2] for cell in Xlabs]
            Ylabs = [gene_list[i] for i in row_order]
            Ycolors = [matrix_data.gene_label_map[gene][0] for gene in Ylabs]
            Ygroup_labels= [matrix_data.gene_label_map[gene][2] for gene in Ylabs]
            row_colors = pd.DataFrame({'Gene Groups': Ycolors},index=Ylabs)
            cg = sns.clustermap(cluster_df, method=args.method, metric=args.metric, z_score=z_choice, row_colors=row_colors, col_colors=col_colors, figsize=(width_heatmap, len_heatmap), cmap =cmap)
        elif cell_color_map:
            Xlabs = [cell_list[i] for i in col_order]
            Xcolors = {}
            for group_index , group_name in enumerate(groups_list):
                Xcolors[group_name]=[cell_color_map[cell][group_index][0] for cell in Xlabs]
            col_colors = pd.DataFrame(Xcolors,index=Xlabs)
            Xgroup_labels = [cell_color_map[cell][0][2] for cell in Xlabs]
            cg = sns.clustermap(cluster_df, method=args.method, metric=args.metric, z_score=z_choice, col_colors=col_colors, figsize=(width_heatmap, len_heatmap), cmap =cmap)
        elif matrix_data.gene_label_map and stablity_plot_num == 0:
            Ylabs = [gene_list[i] for i in row_order]
            Ycolors = [matrix_data.gene_label_map[gene][0] for gene in Ylabs]
            Ygroup_labels= [matrix_data.gene_label_map[gene][2] for gene in Ylabs]
            row_colors = pd.DataFrame({'Gene Groups': Ycolors},index=Ylabs)
            cg = sns.clustermap(cluster_df, method=args.method, metric=args.metric, z_score=z_choice,row_colors=row_colors, figsize=(width_heatmap, len_heatmap), cmap =cmap)

        cg.ax_heatmap.set_title(title, y=title_set)
        cg.cax.set_title('Z-score')
        if cell_color_map:

            for xtick, xcolor, xgroup_name in zip(cg.ax_heatmap.get_xticklabels(), [cell_color_map[cell][0][0] for cell in Xlabs], Xgroup_labels):
                xtick.set_color(xcolor)
                xtick.set_rotation(270)
                xtick.set_fontsize(fontsize)

            leg_handles_cell = {}

            for group_index , group_name in enumerate(groups_list):
                group_seen_cell = []
                leg_handles_cell[group_name] = []
                for xcolor, xgroup_name in zip([cell_color_map[cell][group_index][0] for cell in Xlabs], [cell_color_map[cell][group_index][2] for cell in Xlabs]):
                    if xgroup_name not in group_seen_cell:
                        leg_handles_cell[group_name].append(patches.Patch(color=xcolor, label=xgroup_name))
                        group_seen_cell.append(xgroup_name)
        else:
            for xtick in cg.ax_heatmap.get_xticklabels():
                xtick.set_rotation(270)
                xtick.set_fontsize(fontsize)
        if matrix_data.gene_label_map and stablity_plot_num == 0:
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
        if matrix_data.gene_label_map and cell_color_map and stablity_plot_num == 0:
            gene_legend = cg.ax_heatmap.legend(handles=leg_handles_gene, loc=2, bbox_to_anchor=(1.04, 1-(0.2*len(groups_list))), title='Gene groups', prop={'size':fontsize})
            plt.setp(gene_legend.get_title(),fontsize=fontsize)
            cg.ax_heatmap.add_artist(gene_legend)
            for group_count ,group_name in enumerate(groups_list):
                cell_legend = cg.ax_heatmap.legend(handles=leg_handles_cell[group_name], loc=2, bbox_to_anchor=(1.04, 1-(0.2*group_count)), title='Cell groups '+group_name, prop={'size':fontsize})
                if group_count == 0:
                    plt.setp(cell_legend.get_title(),fontsize=fontsize)
                    #cg.ax_heatmap.add_artist(cell_legend)
                else:
                    #plt.setp(cell_legend.get_title(),fontsize=fontsize)
                    cg.ax_heatmap.add_artist(cell_legend)
        elif cell_color_map:
            for group_count ,group_name in enumerate(groups_list):
                cell_legend = cg.ax_heatmap.legend(handles=leg_handles_cell[group_name], loc=2, bbox_to_anchor=(1.04, 1-(0.2*group_count)), title='Cell groups '+group_name, prop={'size':fontsize})
                if group_count == 0:
                    plt.setp(cell_legend.get_title(),fontsize=fontsize)

                else:
                    cg.ax_heatmap.add_artist(cell_legend)

        elif matrix_data.gene_label_map and stablity_plot_num == 0:
            gene_legend = cg.ax_heatmap.legend(handles=leg_handles_gene, loc=2, bbox_to_anchor=(1.04, 0.8), title='Gene groups', prop={'size':fontsize})
            plt.setp(gene_legend.get_title(),fontsize=fontsize)

        cell_linkage = cg.dendrogram_col.linkage

        link_mat = pd.DataFrame(cell_linkage,
                    columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
                    index=['cluster %d' %(i+1) for i in range(cell_linkage.shape[0])])
        if plot:
            if title != '':
                save_name = '_'.join(title.split(' ')[0:2])
                #plot_url = py.plot_mpl(cg)
                cg.savefig(os.path.join(matrix_data.new_filepath, save_name+'_heatmap.'+args.image_format), bbox_inches='tight')
            else:
                #plot_url = py.plot_mpl(cg)
                cg.savefig(os.path.join(matrix_data.new_filepath,'Group0_Heatmap_all_cells.'+args.image_format), bbox_inches='tight')
            plt.close('all')
        return cell_linkage, df_by_gene[gene_list], col_order
    except FloatingPointError:
        print('Linkage distance has too many zeros. Filter to remove non-expressed genes in order to produce heatmap. Heatmap with '+ str(len(cell_list))+' will not be created.')
        return False, False, False

def make_subclusters(args, cc, matrix_data):
    '''
    Walks a histogram branch map 'cc' and does PCA (SVD), heatmap and correlation search for each non-overlapping
    tree branch. Stops at defined cluster_size (default is 20).
    '''
    try:
        from .dim_reduction import return_top_pca_gene, plot_SVD
    except (SystemError, ValueError, ImportError):
        from dim_reduction import return_top_pca_gene, plot_SVD

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
        #run all cell groups that are subgroups of the parent and greater than or equal to the cutoff args.depth_of_clustering
        if num_members < p_num and num_members >= args.depth_of_clustering:
            group_ID+=1
            #save name for all files generated within this cluster i.e. 'Group_2_with_105_cells_heatmap.pdf'
            current_title = 'Group_'+str(group_ID)+'_with_'+str(num_members)+'_cells'

            if isinstance(matrix_data.log2_df_cell_gene_restricted,pd.DataFrame):
                cell_subset = matrix_data.log2_df_cell_gene_restricted.loc[:,list(set(cell_list))]
                gene_subset = cell_subset.transpose()
                gene_subset = gene_subset.loc[:,(gene_subset!=0).any()]

                full_cell_subset = matrix_data.log2_df_cell.loc[:,list(set(cell_list))]
                full_gene_subset = full_cell_subset.transpose()
                full_gene_subset = full_gene_subset.loc[:,(full_gene_subset!=0).any()]
                num_genes = min(args.gene_number, len(matrix_data.short_gene_list))
                top_pca_by_gene, top_pca = return_top_pca_gene(args, cell_subset)
                plot_SVD(args, matrix_data, gene_subcluster_matrix=gene_subset, title=current_title)
            else:
                full_cell_subset = matrix_data.log2_df_cell.loc[:,list(set(cell_list))]
                full_gene_subset = full_cell_subset.transpose()
                full_gene_subset = full_gene_subset.loc[:,(full_gene_subset!=0).any()]
                top_pca_by_gene, top_pca = return_top_pca_gene(args, full_cell_subset)
                plot_SVD(args,matrix_data, gene_subcluster_matrix=full_gene_subset, title=current_title)

            norm_df_cell1 = np.exp2(full_cell_subset)
            norm_df_cell = norm_df_cell1 -1
            norm_df_cell.to_csv(os.path.join(matrix_data.new_filepath, args.base_name+'_'+current_title+'_matrix.txt'), sep = '\t', index_label=0)

            if len(top_pca)<args.gene_number:
                plot_num = len(top_pca)
            else:
                plot_num = args.gene_number

            if top_pca != []:
                top_pca_by_cell = top_pca_by_gene.transpose()
                #if no_corr flag is provided (False) no correlation plots will be made
                if not args.no_corr:
                    try:
                        from .correlation import corr_plot
                    except (SystemError, ValueError, ImportError):
                        from correlation import corr_plot
                    top_genes_search = top_pca[0:50]
                    corr_plot(top_genes_search, full_gene_subset, args, matrix_data, title = current_title)

                    cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, matrix_data, matrix_subset = top_pca_by_gene, title=current_title)
                plt.close('all')
            else:
                print('Search for top genes by PCA failed in '+current_title+'. No plots will be generated for this subcluster. ')
                pass
