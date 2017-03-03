
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
import numpy as np
from collections import Counter



def group_gene_expression(args, matrix_data, threshold_of_expression = 1, from_kmeans='', alt_color_dict = {}, kmeans_groups_file = ''):
    path_filename = matrix_data.new_filepath
    full_by_cell_df = matrix_data.data_by_cell

    if alt_color_dict:
        color_dict_cell = alt_color_dict
        cell_groups_df = pd.read_table(open(kmeans_groups_file,'rU'), sep='\s+', engine='python')
        group_name_list = [str(gp) for gp in list(set(cell_groups_df['GroupID']))]

    else:
        color_dict_cell = matrix_data.color_dict_cells
        cell_groups_df = pd.read_table(open(matrix_data.cell_list_filepath,'rU'), sep='\s+', engine='python')
        group_name_list = [str(gp) for gp in list(set(cell_groups_df['GroupID']))]


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
        barplot_dict[str(name)] = {'# genes expressed':[], 'pvalues':[], 'Vs':[]}
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
        cell_groups_df = pd.read_table(open(kmeans_groups_file,'rU'), sep='\s+', engine='python')
        primary_group_name = 'GroupID'
        group_name_list = [str(gp) for gp in list(set(cell_groups_df[primary_group_name]))]

    else:
        color_dict_cell = matrix_data.color_dict_cells
        cell_groups_df = matrix_data.cell_groups_df
        group_name_list = matrix_data.cell_group_names
        primary_group_name = group_name_list[0]
        group_name_list = [str(gp) for gp in list(set(cell_groups_df[primary_group_name]))]
    if len(list(set(group_name_list))) ==1:
        print("MultiGroupSig Error: only one group exists")
        return(None)
    plot_pvalue = False
    from rpy2.robjects.packages import importr
    from rpy2.robjects.vectors import FloatVector
    stats = importr('stats')

    group_pairs = list(set(itertools.combinations(group_name_list,2)))
    gene_list = full_by_cell_df.index.tolist()
    cell_group_ident_0 = zip(cell_groups_df['SampleID'],cell_groups_df[primary_group_name])
    cell_group_ident= [c for c in cell_group_ident_0]
    barplot_dict = {}
    vs_dict = {}
    by_gene_df = full_by_cell_df.transpose()
    best_gene_list = []
    best_gene_groups = []
    best_vs_list =[]
    best_pvalue_list = []
    all_gene_plot_dict = {}
    for name in group_name_list:
        barplot_dict[str(name)] = {'genes':[], 'pvalues':[], 'fold_change':[], 'Vs':[]}
        all_gene_plot_dict[str(name)] = []
        vs_dict[str(name)] = []

    gp_present_list = []
    for gp in group_pairs:
        gp = [str(x) for x in gp]
        index_list = []
        sig_gene_list = []
        gp_vs_all_seen = []
        sig_vs_all_gene_list = []
        cell_list1 = [c[0] for c in cell_group_ident if str(c[1]) == gp[0]]
        cell_list2 = [c[0] for c in cell_group_ident if str(c[1]) == gp[1]]
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
                pvalues_all = np.asarray([p[1] for p in sig_vs_all_gene_list], dtype=float)

                p_adjust_all = np.asarray(stats.p_adjust(FloatVector(pvalues_all), method = 'BH'), dtype=float)
                gene_index_all = [ge[0] for ge in sig_vs_all_gene_list]
                mean_log2_exp_list_all = []
                sig_1_2_list_all = []
                mean1_list_all = []
                mean2_list_all = []
                for sig_gene in gene_index_all:
                    sig_gene_df_all = by_gene_df.loc[:,sig_gene]
                    mean_log2_exp_list_all.append(np.asarray(sig_gene_df_all, dtype=float).mean())
                    cell1_all = by_gene_df.loc[cell1_present,sig_gene]
                    mean_cell1_all = np.asarray(cell1_all, dtype=float).mean()
                    mean1_list_all.append(mean_cell1_all)
                    cell_other = by_gene_df.loc[all_other_cells,sig_gene]
                    mean_cell_other = np.asarray(cell_other, dtype=float).mean()
                    mean2_list_all.append(mean_cell_other)
                    ratio_1_other = (mean_cell1_all+1)/(mean_cell_other+1)
                    sig_1_2_list_all.append(ratio_1_other)
                sig_vs_other_dict = {'pvalues':pvalues_all,'adjusted_p_values':p_adjust_all,'mean_all':np.asarray(mean_log2_exp_list_all,dtype=float), 'mean_'+str(gp[0]):np.asarray(mean1_list_all, dtype=float), 'mean_all_other':np.asarray(mean2_list_all,dtype=float), 'ratio '+str(gp[0])+' to everything':np.asarray(sig_1_2_list_all, dtype=float)}

                sig_df_vs_other = pd.DataFrame(sig_vs_other_dict, index=gene_index_all)

                sig_df_vs_other.to_csv(os.path.join(multi_sig_filename,'sig_'+str(gp[0])+'_VS_all_other_pvalues.txt'), sep = '\t')
                gp_vs_all_seen.append(gp[0])
            pvalues = [p[1] for p in sig_gene_list]
            p_adjust = stats.p_adjust(FloatVector(pvalues), method = 'BH')
            gene_index = [ge[0] for ge in sig_gene_list]

            mean_log2_exp_list = []
            sig_1_2_list = []
            sig_2_1_list = []
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
                ratio_2_1 = (mean_cell2+1)/(mean_cell1+1)
                sig_1_2_list.append(ratio_1_2)
                sig_2_1_list.append(ratio_2_1)
            sig_df_1 = pd.DataFrame({'pvalues':pvalues,'adjusted_p_values':p_adjust,'mean_all':mean_log2_exp_list,'mean_'+gp[0]:mean1_list, 'mean_'+gp[1]:mean2_list, 'ratio '+gp[0]+' to '+gp[1]:sig_1_2_list}, index=gene_index)
            sig_df_2 = pd.DataFrame({'pvalues':pvalues,'adjusted_p_values':p_adjust,'mean_all':mean_log2_exp_list,'mean_'+gp[1]:mean2_list, 'mean_'+gp[1]:mean1_list, 'ratio '+gp[1]+' to '+gp[0]:sig_2_1_list}, index=gene_index)
            cell_names_df_1 = pd.DataFrame({gp[0]+'_cells':pd.Series(cell1_present, index=range(len(cell1_present))), gp[1]+'_cells2':pd.Series(cell2_present, index=range(len(cell2_present)))})
            cell_names_df_2 = pd.DataFrame({gp[1]+'_cells':pd.Series(cell2_present, index=range(len(cell2_present))), gp[0]+'_cells2':pd.Series(cell1_present, index=range(len(cell1_present)))})
            sig_df_1.to_csv(os.path.join(multi_sig_filename,'sig_'+gp[0]+'_VS_'+gp[1]+'_pvalues.txt'), sep = '\t')
            sig_df_2.to_csv(os.path.join(multi_sig_filename,'sig_'+gp[1]+'_VS_'+gp[0]+'_pvalues.txt'), sep = '\t')
            cell_names_df_1.to_csv(os.path.join(multi_sig_filename,'sig_'+gp[0]+'_VS_'+gp[1]+'_cells.txt'), sep = '\t')
            cell_names_df_2.to_csv(os.path.join(multi_sig_filename,'sig_'+gp[1]+'_VS_'+gp[0]+'_cells.txt'), sep = '\t')

            #select only genes that have a fold change greater than 1.3 and are significant by adusted-p-value
            top_fc_df1 = sig_df_1.loc[(sig_df_1['ratio '+gp[0]+' to '+gp[1]]>1.3) & (sig_df_1['adjusted_p_values']<=0.05)]
            top_fc_df2 = sig_df_2.loc[(sig_df_2['ratio '+gp[1]+' to '+gp[0]]>1.3) & (sig_df_2['adjusted_p_values']<=0.05)]
            two_df_list = [top_fc_df1 , top_fc_df2]

            #make dictionary of values to plot for significant genes, for both 1vs2 and 2vs1, eg.
            for i, fc_df in enumerate(two_df_list):
                top_fc_df = fc_df.sort_values(by='adjusted_p_values',axis=0, ascending=True)
                genes = top_fc_df.index.tolist()
                pvalues = top_fc_df['adjusted_p_values'].tolist()
                if i == 0:
                    fc = top_fc_df.loc[:,'ratio '+gp[0]+' to '+gp[1]].tolist()
                    #create list of tuples with (gene, pvalue, foldchange) for each significant gene
                    z_all = zip(genes,pvalues,fc)

                    if matrix_data.exclude_list != None:
                        top_t2 = [g for g in z_all if g[0] not in matrix_data.exclude_list]
                    else:
                        top_t2 = z_all
                    top = [list(t) for t in zip(*top_t2)]
                    try:
                        if len(top[0]) > 0:
                            all_gene_plot_dict[gp[0]].append(top)
                            vs_dict[gp[0]].append(gp[1])
                            gp_present_list.append(gp[0])
                    except IndexError:
                        print(gp[0]+' vs. '+gp[1]+' sig test failed.')
                #same operation on the inverse comparison
                elif i == 1:
                    fc = top_fc_df.loc[:,'ratio '+gp[1]+' to '+gp[0]].tolist()
                    #create list of tuples with (gene, pvalue, foldchange) for each significant gene
                    z_all = zip(genes,pvalues,fc)

                    if matrix_data.exclude_list != None:
                        top_t2 = [g for g in z_all if g[0] not in matrix_data.exclude_list]
                    else:
                        top_t2 = z_all
                    top = [list(t) for t in zip(*top_t2)]
                    try:
                        if len(top[0]) > 0:
                            all_gene_plot_dict[gp[1]].append(top)
                            vs_dict[gp[1]].append(gp[0])
                            gp_present_list.append(gp[1])
                    except IndexError:
                        print(gp[1]+' vs. '+gp[0]+' sig test failed.')
        else:
            if cell1_present == []:
                if args.verbose:
                    print(gp[1], 'not present in cell matrix')
            else:
                if args.verbose:
                    print(gp[0], 'not present in cell matrix')
    plt.rc('font', weight='bold')
    fig, axs = plt.subplots(1, len(group_name_list), figsize=(42+len(group_name_list),16), sharex=False, sharey=False)
    axs = axs.ravel()
    color_map = {}
    best_gene_df_list = []
    list_of_gene_lists = []
    for g_index, gp in enumerate(group_name_list):
        for l in all_gene_plot_dict[str(gp)]:
            if l != []:
                list_of_gene_lists.append(l[0])
        if list_of_gene_lists != []:
            flat_list = [item for sublist in list_of_gene_lists for item in sublist]
            gene_counter = Counter(itertools.chain.from_iterable(flat_list))
            common_genes = []
            for k,v in gene_counter.most_common():
                if v > 1:
                    common_genes.append(k)
            vs_list = vs_dict[str(gp)]
            plot_vs_df_list = []
            for i, vs_name in enumerate(vs_list):
                top_df1 = pd.DataFrame(all_gene_plot_dict[gp][i])
                top_df = top_df1.T
                top_df.columns = ['GeneID', 'adjusted_pvalue', 'fold_change']
                top_df.sort_values(by='adjusted_pvalue', inplace=True)
                if not args.sig_unique:
                    common_df = top_df[top_df['GeneID'].isin(common_genes)]
                    if not common_df.empty:
                        final_gp_df = common_df.sort_values(by='adjusted_pvalue')
                        if len(common_df['GeneID'] <= sig_to_plot):
                            final_gp_df = pd.concat([common_df,top_df])
                    else:
                        final_gp_df = top_df
                else:
                    final_gp_df = top_df.drop_duplicates('GeneID')


                final_gp_df['GroupID'] = pd.Series([gp for c in top_df['GeneID']])
                final_gp_df['Vs'] = pd.Series(["Significance vs. "+vs_name for c in top_df['GeneID']])


                if len(final_gp_df['GeneID']) > sig_to_plot:
                    plot_gp_df_vs = final_gp_df.loc[list(range(0,sig_to_plot))]
                else:
                    plot_gp_df_vs = final_gp_df

                plot_vs_df_list.append(plot_gp_df_vs)
                best_gene_df_list.append(final_gp_df)

                color_map["Significance vs. "+vs_name] = color_dict_cell[vs_name][0]
            try:
                final_plot_df = pd.concat(plot_vs_df_list)
            except ValueError:
                break
            if args.sig_unique:
                final_plot_df = pd.concat(plot_vs_df_list).sort_values(by='adjusted_pvalue', inplace=True)
                final_plot_df.drop_duplicates('GeneID', inplace=True)
            if color_map != {}:
                g = sns.barplot(x='adjusted_pvalue', y='GeneID', hue='Vs', data=final_plot_df, ax = axs[g_index], palette=color_map)
                axs[g_index].set_xscale("log", nonposx='clip')
                bar_list = []
                if plot_pvalue:
                    for p in axs[g_index].patches:
                        height = p.get_height()
                        width = p.get_width()
                        bar_list.append(width)
                    max_bar = max(bar_list)
                    for p in axs[g_index].patches:
                        height = p.get_height()
                        width = p.get_width()
                        axs[g_index].text(max_bar*50,p.get_y()+(height), "{:.2e}".format(width))
                rect = axs[g_index].patch
                rect.set_facecolor('white')
                #sns.despine(left=True, bottom=True, top=True)
                axs[g_index].invert_xaxis()
                axs[g_index].xaxis.set_ticks_position('none')
                axs[g_index].yaxis.tick_right()
                axs[g_index].set_title(gp)
                axs[g_index].legend(loc='upper left', bbox_to_anchor=(0.01, 1.11+(0.01*len(group_name_list))), ncol=1, prop={'size':18})
                axs[g_index].set_xlabel('adjusted p-value')
                for xmaj in axs[g_index].xaxis.get_majorticklocs():
                    axs[g_index].axvline(x=xmaj,ls='--', lw = 0.5, color='grey', alpha=0.3)
                axs[g_index].xaxis.grid(True, which="major", linestyle='-')
                plt.subplots_adjust(left=.08, wspace=.3)
    best_gene_df = pd.concat(best_gene_df_list)
    if from_kmeans == '':
        plt.savefig(os.path.join(path_filename,'differential_genes_foldchanges.'+args.image_format), bbox_inches='tight')
        best_gene_df.to_csv(os.path.join(path_filename,'Best_Gene_list.txt'), sep = '\t')
    else:
        plt.savefig(os.path.join(path_filename,from_kmeans+'_differential_genes_foldchanges.'+args.image_format), bbox_inches='tight')
        best_gene_df.to_csv(os.path.join(path_filename,from_kmeans+'_Best_Gene_list.txt'), sep = '\t')
