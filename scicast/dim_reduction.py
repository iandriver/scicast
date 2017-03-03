import numpy as np
import pandas as pd
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA as skPCA
import seaborn as sns
import matplotlib.patches as patches
from sklearn.manifold import TSNE
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans



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

def return_top_pca_gene(args, by_cell_matrix, user_num_genes = None):
    gene_pca = skPCA(n_components=3)
    np_by_gene = np.asarray(by_cell_matrix.transpose())
    gene_index = by_cell_matrix.index.tolist()

    if user_num_genes is not None:
        num_genes = user_num_genes
    else:
        num_genes = min(args.gene_number, len(gene_index))
    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(gene_pca.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=gene_index)
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len(gene_index))
    top_pca_list = Pc_sort_df.index.tolist()
    new_cell_matrix = by_cell_matrix.ix[top_pca_list[0:num_genes],:]
    return new_cell_matrix.transpose(), top_pca_list[0:num_genes]

def plot_PCA(args, matrix_data, title = '', gene_subcluster_matrix=False):
    sns.set_palette("RdBu_r", 10, 1)
    if gene_subcluster_matrix:
        gene_list = gene_subcluster_matrix.columns.tolist()
        len_gene_list = len(gene_list)
        num_genes = min(args.gene_number, len_gene_list)

        df_by_gene = gene_subcluster_matrix.transpose()

    elif isinstance(matrix_data.log2_df_cell_gene_restricted,pd.DataFrame):
        gene_list = matrix_data.short_gene_list
        len_gene_list = len(gene_list)
        num_genes = min(args.gene_number, len_gene_list)

        df_by_gene = matrix_data.log2_df_cell_gene_restricted.transpose()

    else:
        gene_list = matrix_data.gene_list
        len_gene_list = len(gene_list)
        num_genes = min(args.gene_number, len_gene_list)

        df_by_gene = matrix_data.log2_df_cell.transpose()

    gene_pca = skPCA(n_components=3)
    np_by_gene = np.asarray(df_by_gene)

    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(gene_pca.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=gene_list)
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len_gene_list)
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
        if matrix_data.cell_label_map:
            annotate = args.annotate_cell_pca
            X = [x for x in top_cell_trans[:, 0]]
            Y = [y for y in top_cell_trans[:, 1]]
            for group_index , group_name in enumerate(matrix_data.cell_group_names):
                if group_index == 0:
                    labels = [matrix_data.cell_label_map[cell][group_index][2] for cell in top_by_cell.columns.tolist()]
                    colors = [matrix_data.cell_label_map[cell][group_index][0] for cell in top_by_cell.columns.tolist()]
                elif group_index > 0:
                    labels = [a+'_'+b for a, b in zip(labels, [matrix_data.cell_label_map[cell][group_index][2] for cell in top_by_cell.columns.tolist()])]

            markers = [matrix_data.cell_label_map[cell][0][1] for cell in top_by_cell.columns.tolist()]

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
        if matrix_data.cell_label_map:
            handles, labs = ax_cell.get_legend_handles_labels()
            # sort both labels and handles by labels
            labs, handles = zip(*sorted(zip(labs, handles), key=lambda t: t[0]))
            ax_cell.legend(handles, labs, loc='best', ncol=1, prop={'size':12}, markerscale=1.5, frameon=True)
        ax_cell.set_xlabel('PC1')
        ax_cell.set_ylabel('PC2')
        if annotate:
            for label, x, y in zip(top_by_cell.columns, top_cell_trans[:, 0], top_cell_trans[:, 1]):
                ax_cell.annotate(label, (x+0.1, y+0.1))

        if matrix_data.gene_label_map:
            X = [x for x in top_gene_trans[:, 0]]
            Y = [y for y in top_gene_trans[:, 1]]
            labels = top_by_gene.columns.tolist()
            markers = [matrix_data.gene_label_map[gene][1] for gene in top_by_gene.columns.tolist()]
            colors = [matrix_data.gene_label_map[gene][0] for gene in top_by_gene.columns.tolist()]
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

        if title != '':
            save_name = '_'.join(title.split(' ')[0:2])
            plt.savefig(os.path.join(matrix_data.new_filepath,save_name+'_skpca.'+args.image_format), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(matrix_data.new_filepath,'Group0_skpca.'+args.image_format), bbox_inches='tight')
        plt.close('all')
        return top_pca_list
    else:
        return []

def plot_SVD(args, matrix_data, gene_subcluster_matrix=False, title = ''):
    sns.set(context= 'poster', font_scale = 1.2)
    sns.set_palette("RdBu_r", 10, 1)
    if isinstance(gene_subcluster_matrix, pd.DataFrame):
        gene_list = gene_subcluster_matrix.columns.tolist()
        len_gene_list = len(gene_list)
        num_genes = min(args.gene_number, len_gene_list)

        df_by_gene = gene_subcluster_matrix

    elif isinstance(matrix_data.log2_df_cell_gene_restricted, pd.DataFrame):
        gene_list = matrix_data.short_gene_list
        len_gene_list = len(gene_list)
        num_genes = min(args.gene_number, len_gene_list)

        df_by_gene = matrix_data.log2_df_cell_gene_restricted.transpose()

    else:
        gene_list = [x for x in matrix_data.gene_list]
        len_gene_list = len(gene_list)
        num_genes = min(args.gene_number, len_gene_list)

        df_by_gene = matrix_data.log2_df_cell.transpose()
    gene_pca = TruncatedSVD(n_components=3)
    np_by_gene = np.asarray(df_by_gene)
    gene_list = df_by_gene.columns.tolist()
    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    pca_t = gene_pca.components_.T
    Pc_df = pd.DataFrame(pca_t, columns=['PC-1', 'PC-2', 'PC-3'], index=df_by_gene.columns.tolist())
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len_gene_list)
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
        if matrix_data.cell_label_map:
            annotate = args.annotate_cell_pca
            X = [x for x in top_cell_trans[:, 0]]
            Y = [y for y in top_cell_trans[:, 1]]
            for group_index , group_name in enumerate(matrix_data.cell_group_names):
                if group_index == 0:
                    labels = [matrix_data.cell_label_map[cell][group_index][2] for cell in top_by_cell.columns.tolist()]
                    colors = [matrix_data.cell_label_map[cell][group_index][0] for cell in top_by_cell.columns.tolist()]
                elif group_index > 0:
                    labels = [a+'_'+b for a, b in zip(labels, [matrix_data.cell_label_map[cell][group_index][2] for cell in top_by_cell.columns.tolist()])]
            label_set = list(set(labels))
            markers_dict = {}
            for i, lab_set in enumerate(label_set):
                markers_dict[lab_set] = matrix_data.markers[i]
            markers = [markers_dict[x] for x in labels]

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
        if matrix_data.cell_label_map:
            handles, labs = ax_cell.get_legend_handles_labels()
            # sort both labels and handles by labels
            labs, handles = zip(*sorted(zip(labs, handles), key=lambda t: t[0]))
            ax_cell.legend(handles, labs, loc='best', ncol=1, prop={'size':12}, markerscale=1.5, frameon=True)
        ax_cell.set_xlabel('PC1')
        ax_cell.set_ylabel('PC2')
        if annotate:
            for label, x, y in zip(top_by_cell.columns, top_cell_trans[:, 0], top_cell_trans[:, 1]):
                ax_cell.annotate(label, (x+0.1, y+0.1))

        if matrix_data.gene_label_map:
            X = [x for x in top_gene_trans[:, 0]]
            Y = [y for y in top_gene_trans[:, 1]]
            labels = top_by_gene.columns.tolist()
            markers = [matrix_data.gene_label_map[gene][1] for gene in top_by_gene.columns.tolist()]
            colors = [matrix_data.gene_label_map[gene][0] for gene in top_by_gene.columns.tolist()]
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
        if title != '':
            save_name = '_'.join(title.split(' ')[0:2])
            plt.savefig(os.path.join(matrix_data.new_filepath,save_name+'_TruncatedSVD.'+args.image_format), bbox_inches='tight')
            #plot_url = py.plot_mpl(fig)
        else:
            #plot_url = py.plot_mpl(fig)
            plt.savefig(os.path.join(matrix_data.new_filepath,'Group0_TruncatedSVD.'+args.image_format), bbox_inches='tight')
        plt.close('all')
        return top_pca_list
    else:
        return []

#create cell and gene TSNE scatter plots (one file for both)
def plot_TSNE(args, matrix_data, title= ''):
    sns.set(context= 'poster', font_scale = 1.2)
    sns.set_palette("RdBu_r", 10, 1)
    if isinstance(matrix_data.log2_df_cell_gene_restricted,pd.DataFrame):
        gene_list = matrix_data.short_gene_list
        len_gene_list = len(gene_list)
        num_genes = min(args.gene_number, len_gene_list)

        df_by_gene = matrix_data.log2_df_cell_gene_restricted.transpose()
    else:
        gene_list = matrix_data.gene_list
        len_gene_list = len(gene_list)
        num_genes = min(args.gene_number, len_gene_list)

        df_by_gene = matrix_data.log2_df_cell.transpose()

    gene_pca = TruncatedSVD(n_components=3)
    np_by_gene = np.asarray(df_by_gene)

    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(gene_pca.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=df_by_gene.columns.tolist())
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len_gene_list)
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
        if matrix_data.cell_label_map:
            annotate = args.annotate_cell_pca
            X = [x for x in top_cell_trans[:, 0]]
            Y = [y for y in top_cell_trans[:, 1]]
            for group_index , group_name in enumerate(matrix_data.cell_group_names):
                if group_index == 0:
                    labels = [matrix_data.cell_label_map[cell][group_index][2] for cell in top_by_cell.columns.tolist()]
                    colors = [matrix_data.cell_label_map[cell][group_index][0] for cell in top_by_cell.columns.tolist()]
                elif group_index > 0:
                    labels = [a+'_'+b for a, b in zip(labels, [matrix_data.cell_label_map[cell][group_index][2] for cell in top_by_cell.columns.tolist()])]

            markers = [matrix_data.cell_label_map[cell][0][1] for cell in top_by_cell.columns.tolist()]

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
        if matrix_data.cell_label_map:
            handles, labs = ax_cell.get_legend_handles_labels()
            # sort both labels and handles by labels
            labs, handles = zip(*sorted(zip(labs, handles), key=lambda t: t[0]))
            ax_cell.legend(handles, labs, loc='best', ncol=1, prop={'size':12}, markerscale=1.5, frameon=True)
        ax_cell.set_xlabel('PC1')
        ax_cell.set_ylabel('PC2')
        if annotate:
            for label, x, y in zip(top_by_cell.columns, top_cell_trans[:, 0], top_cell_trans[:, 1]):
                ax_cell.annotate(label, (x+0.1, y+0.1))

        if matrix_data.gene_label_map:
            X = [x for x in top_gene_trans[:, 0]]
            Y = [y for y in top_gene_trans[:, 1]]
            labels = top_by_gene.columns.tolist()
            markers = [matrix_data.gene_label_map[gene][1] for gene in top_by_gene.columns.tolist()]
            colors = [matrix_data.gene_label_map[gene][0] for gene in top_by_gene.columns.tolist()]
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

        if title != '':
            save_name = '_'.join(title.split(' ')[0:2])
            plt.savefig(os.path.join(matrix_data.new_filepath,save_name+'_TSNE.'+args.image_format), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(matrix_data.new_filepath,'Group0_TSNE.'+args.image_format), bbox_inches='tight')
        plt.close('all')
        return top_pca_list
    else:
        return []

#create kmeans cluster and silhouette_score plot
def plot_kmeans(args, matrix_data, kmeans_range, title=''):
    sns.set(context= 'poster', font_scale = 1.2)
    from sklearn.metrics import silhouette_samples, silhouette_score
    try:
        from .heatmaps import clust_heatmap
        from .significance_testing import multi_group_sig
    except (SystemError, ValueError, ImportError):
        from heatmaps import clust_heatmap
        from significance_testing import multi_group_sig
    import matplotlib.cm as cm
    sns.set_palette("RdBu_r", 10, 1)
    path_filename = matrix_data.new_filepath
    if isinstance(matrix_data.log2_df_cell_gene_restricted,pd.DataFrame):
        gene_list = matrix_data.short_gene_list
        len_gene_list = len(gene_list)
        num_genes = min(args.gene_number, len_gene_list)
        df_by_gene = matrix_data.log2_df_cell_gene_restricted.transpose()

    else:
        gene_list = matrix_data.gene_list
        len_gene_list = len(gene_list)
        num_genes = min(args.gene_number, len_gene_list)
        df_by_gene = matrix_data.log2_df_cell.transpose()


    gene_pca = TruncatedSVD(n_components=3)
    np_by_gene = np.asarray(df_by_gene)
    #find top genes by PCA
    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(gene_pca.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=gene_list)
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len_gene_list)
    top_pca_list = Pc_sort_df.index.tolist()

    top_by_gene = df_by_gene[top_pca_list[0:num_genes]]

    #use 2D t-SNE for kmeans clustering if option is selected
    if args.use_TSNE:
        gene_top = TSNE(n_components=2, init='pca', random_state=0)
        cell_pca = TSNE(n_components=2, init='pca', random_state=0)
    else:
        gene_top = TruncatedSVD(n_components=2)
        cell_pca = TruncatedSVD(n_components=2)
    top_by_cell = top_by_gene.transpose()
    np_top_gene = np.asarray(top_by_cell)
    np_top_cell = np.asarray(top_by_gene)
    top_cell_trans = cell_pca.fit_transform(np_top_cell)
    top_gene_trans = gene_top.fit_transform(np_top_gene)

    range_n_clusters = range(kmeans_range[0],kmeans_range[1]+1)

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
        if args.verbose:
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

        plt.savefig(os.path.join(path_filename,'Group0_kmeans_'+str(n_clusters)+'_clusters.'+args.image_format), bbox_inches='tight')
        plt.close('all')

        #use colors to make label map compatable with heatmap
        color_dict ={}
        for cell, m_color_tup in zip(top_by_cell.columns, zip(colors,[matrix_data.markers[pred] for pred in top_cell_pred],['kmeans_'+str(p) for p in top_cell_pred])):
            color_dict[cell] = []
            color_dict[cell].append(m_color_tup)

        group_color_dict = dict(zip(['kmeans_'+str(p) for p in top_cell_pred],zip(colors,[matrix_data.markers[pred] for pred in top_cell_pred])))
        #run heatmap with kmeans clustering and colors

        top_pca_by_gene, top_pca = return_top_pca_gene(args, df_by_gene.transpose())
        cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, matrix_data, matrix_subset=top_pca_by_gene, title= 'kmeans_label_with_'+str(n_clusters)+'_clusters', kmeans_color_map=color_dict)
        if args.kmeans_sig_test:
            sns.set(context= 'poster', font_scale = 1.2-(n_clusters*.08))
            multi_group_sig(args, matrix_data, from_kmeans='kmeans_'+str(n_clusters), alt_color_dict=group_color_dict, kmeans_groups_file=cell_group_path)
