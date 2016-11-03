import numpy as np
import pandas as pd
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import itertools
import matplotlib.patches as patches
import warnings


def group_gene_expression(args, matrix_data, threshold_of_expression = 1, from_kmeans='', alt_color_dict = {}, kmeans_groups_file = ''):
    path_filename = matrix_data.new_filepath
    full_by_cell_df = matrix_data.data_by_cell

    if alt_color_dict:
        color_dict_cell = alt_color_dict
        cell_groups_df = pd.read_csv(open(kmeans_groups_file,'rU'), sep=None, engine='python')
        group_name_list = list(set(cell_groups_df['GroupID']))

    else:
        color_dict_cell = matrix_data.color_dict_cells
        cell_groups_df = pd.read_csv(open(matrix_data.cell_list_filepath,'rU'), sep=None, engine='python')
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
        barplot_dict[name] = {'# genes expressed':[], 'pvalues':[], 'Vs':[]}
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
            df_by_cell_1 = full_by_cell_df.loc[:,cell1_present]
            df_by_cell_2 = full_by_cell_df.loc[:,cell2_present]
            df_by_cell_other = full_by_cell_df.loc[:,all_other_cells]
            df_by_gene_1 = df_by_cell_1.transpose()
            df_by_gene_2 = df_by_cell_2.transpose()
            df_by_gene_other = df_by_cell_other.transpose()
            for g in gene_list:
                g_pvalue = scipy.stats.f_oneway(df_by_gene_1[g], df_by_gene_2[g])




'''Compares each defined cell group to each other group and returns all genes with p-value and adjusted p-value.
Also creates a best_gene_list with the top genes between each group, by adjusted p-value.
Also creates a barplot of top significant genes between groups (can be unique or not).
'''
def multi_group_sig(args, matrix_data, sig_to_plot = 20, from_kmeans='', alt_color_dict = {}, kmeans_groups_file = ''):
    if not args.verbose:
        warnings.filterwarnings("ignore", category=RuntimeWarning)
    path_filename = matrix_data.new_filepath
    #create seperate file for all of the significance files
    if from_kmeans == '':
        multi_sig_filename = os.path.join(path_filename,'user_defined_group_pairwise_significance_files')
    else:
        multi_sig_filename = os.path.join(path_filename,from_kmeans+'_group_pairwise_significance_files')
    try:
        os.mkdir(multi_sig_filename)
    except OSError:
        if args.verbose:
            print(multi_sig_filename+' already exists. Files will be overwritten.')

    full_by_cell_df = matrix_data.log2_df_cell

    if alt_color_dict:
        color_dict_cell = alt_color_dict
        cell_groups_df = pd.read_csv(open(kmeans_groups_file,'rU'), sep=None, engine='python')
        group_name_list = list(set(cell_groups_df['GroupID']))

    else:
        color_dict_cell = matrix_data.color_dict_cells
        cell_groups_df = pd.read_csv(open(matrix_data.cell_list_filepath,'rU'), sep=None, engine='python')
        group_name_list = list(set(cell_groups_df['GroupID']))
    plot_pvalue = False
    from rpy2.robjects.packages import importr
    from rpy2.robjects.vectors import FloatVector
    stats = importr('stats')

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
            df_by_cell_1 = full_by_cell_df.loc[:,cell1_present]
            df_by_cell_2 = full_by_cell_df.loc[:,cell2_present]
            df_by_cell_other = full_by_cell_df.loc[:,all_other_cells]
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
                    sig_gene_df_all = by_gene_df.loc[:,sig_gene]
                    mean_log2_exp_list_all.append(sig_gene_df_all.mean())
                    cell1_all = by_gene_df.loc[cell1_present,sig_gene]
                    mean_cell1_all = cell1_all.mean()
                    mean1_list_all.append(mean_cell1_all)
                    cell_other = by_gene_df.loc[all_other_cells,sig_gene]
                    mean_cell_other = cell_other.mean()
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
                sig_gene_df = by_gene_df.loc[:,sig_gene]
                mean_log2_exp_list.append(sig_gene_df.mean())

                cell1 = by_gene_df.loc[cell1_present,sig_gene]
                mean_cell1 = cell1.mean()
                mean1_list.append(mean_cell1)
                cell2 = by_gene_df.loc[cell2_present,sig_gene]
                mean_cell2 = cell2.mean()
                mean2_list.append(mean_cell2)
                ratio_1_2 = (mean_cell1+1)/(mean_cell2+1)
                sig_1_2_list.append(ratio_1_2)
            sig_df = pd.DataFrame({'pvalues':pvalues,'adjusted_p_values':p_adjust,'mean_all':mean_log2_exp_list,'mean_'+gp[0]:mean1_list, 'mean_'+gp[1]:mean2_list, 'ratio '+gp[0]+' to '+gp[1]:sig_1_2_list}, index=gene_index)
            cell_names_df = pd.DataFrame({gp[0]+'_cells':pd.Series(cell1_present, index=range(len(cell1_present))), gp[1]+'_cells2':pd.Series(cell2_present, index=range(len(cell2_present)))})
            sig_df.to_csv(os.path.join(multi_sig_filename,'sig_'+gp[0]+'_VS_'+gp[1]+'_pvalues.txt'), sep = '\t')
            cell_names_df.to_csv(os.path.join(multi_sig_filename,'sig_'+gp[0]+'_VS_'+gp[1]+'_cells.txt'), sep = '\t')


            top_fc_df1 = sig_df.loc[(sig_df['ratio '+gp[0]+' to '+gp[1]]>1.3)]

            top_fc_df = top_fc_df1.sort_values(by='adjusted_p_values',axis=0, ascending=True)
            genes = top_fc_df.index.tolist()
            pvalues = top_fc_df['adjusted_p_values'].tolist()
            fc = top_fc_df.loc[:,'ratio '+gp[0]+' to '+gp[1]].tolist()
            z = zip(genes,pvalues,fc)
            z_all = [s for s in z]
            if args.sig_unique:
                top_t = [g for g in z_all if g[0] not in barplot_dict[gp[0]]['genes']]
            else:
                top_t = [g for g in z_all if g[0]]
            if matrix_data.exclude_list != None:
                top_t2 = [g for g in top_t if g[0] not in barplot_dict[gp[0]]['genes'] and g[0] not in matrix_data.exclude_list]
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
                if args.verbose:
                    print(gp[1], 'not present in cell matrix')
            else:
                if args.verbose:
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
        plt.savefig(os.path.join(path_filename,from_kmeans+'_differential_genes_foldchanges.pdf'), bbox_inches='tight')
        best_gene_df = pd.DataFrame({'GeneID':best_gene_list, 'GroupID':best_gene_groups, 'Vs':best_vs_list, 'adjusted_pvalue':best_pvalue_list})
        best_gene_df.to_csv(os.path.join(path_filename,from_kmeans+'_Best_Gene_list.txt'), sep = '\t')
