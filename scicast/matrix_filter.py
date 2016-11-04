
import numpy as np
import pandas as pd
import os
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import math






#Takes raw matrix data and returns matrix object
class Matrix_filter(object):

    def __init__(self, cell_df, args, cell_list_filepath, gene_list_filepath):

        self.cell_list_filepath = cell_list_filepath
        self.gene_list_filepath = gene_list_filepath
        self.cell_names = cell_df.columns.tolist()
        self.gene_names = cell_df.index.tolist()
        self.data_by_cell = cell_df
        if args.base_name:
            self.new_filepath = os.path.join(os.path.dirname(args.filepath),args.base_name+'_scicast_analysis')
        else:
            self.new_filepath = os.path.join(os.path.dirname(args.filepath),'scicast_analysis')


        self.markers = ['o', 'v','D','*','x','h', 's','p','8','^','>','<', 'd','o', 'v','D','*','x','h', 's','p','8','^','>','<', 'd']
        from matplotlib import colors
        all_color_list = list(colors.cnames.keys())
        self.cell_color_list = ['b', 'm', 'r', 'c', 'g','y','k']+all_color_list

        if args.exclude_genes:
            excluded_genes_df = pd.read_csv(os.path.join(os.path.dirname(args.filepath),args.exclude_genes), sep='\t', header=0, index_col=False)
            try:
                self.exclude_list = [g1 for g1 in excluded_genes_df['GeneID'].tolist() if g1 in self.gene_names]
                not_in_mat = [g1 for g1 in excluded_genes_df['GeneID'].tolist() if g1 not in self.gene_names]
            except KeyError:
                sys.exit("Error: Please provide Gene Exclude list file with 'GeneID' as header.")

            self.data_by_cell.drop(self.exclude_list, axis=0, inplace=True)
            if args.verbose:
                print(str(not_in_mat)+' already filtered from matrix.')
        else:
            self.exclude_list = []

        #if cells are to be excluded filter them out and generate log2 matrix otherwise just generate log2_matrix
        if args.limit_cells and self.cell_list_filepath:
            self.make_new_matrix_cell()
        else:
            self.cell_list = self.cell_names
        self.gene_list = self.data_by_cell.index.tolist()

        #if a gene list is provided create the new matrix with only those genes included
        if self.gene_list_filepath:
            self.make_new_matrix_gene(self.gene_list_filepath)
        else:
            self.gene_group_list = []
            self.short_gene_list, self.short_gene_matrix_cell =[],[]


        self.log2_outlierfilter(args)

        #if user input color and or marker values are provided generate color_dict
        if args.color_cells:
            self.color_dict_cells= {}
            color_list1 = args.color_cells.split(' ')
            for i, c in enumerate(color_list1):
                c_pair = c.split(',')
                if len(c_pair) == 2:
                    self.color_dict_cells[c_pair[0]] = [c_pair[1],self.markers[i]]
                elif len(cc_pair == 3):
                    self.color_dict_cells[c_pair[0]] = [c_pair[1],c_pair[2]]
        #otherwise generate generic color dict from colors and markers
        elif self.cell_list_filepath:
            self.color_dict_cells= {}
            cell_groups_df = pd.read_csv(open(self.cell_list_filepath,'rU'), sep=None, engine='python')
            self.cell_group_set = list(set(cell_groups_df['GroupID'].tolist()))
            for g,c,m in zip(self.cell_group_set, self.cell_color_list[0:len(self.cell_group_set)],self.markers[0:len(self.cell_group_set)]):
                self.color_dict_cells[g] =[c,m]

        #if gene user color inputs are provided
        if args.color_genes:
            self.color_dict_genes= {}
        #if the gene group names are the same as cell groups
        if args.color_genes =='same':
            self.color_dict_genes = self.color_dict_cells
        elif args.color_genes:
            color_list1 = args.color_genes.split(' ')
            for i, c in enumerate(color_list1):
                c_pair = c.split(',')
                if len(c_pair) == 2:
                    self.color_dict_genes[c_pair[0]] = [c_pair[1],self.markers[i]]
                elif len(c_pair == 3):
                    self.color_dict_genes[c_pair[0]] = [c_pair[1],c_pair[2]]
        elif self.gene_list_filepath:
            no_groups = False
            self.color_dict_genes= {}
            gene_groups_df = pd.read_csv(open(self.gene_list_filepath,'rU'), sep=None, engine='python')

            gene_group_set = list(set(gene_groups_df['GroupID'].tolist()))
            for i, group in enumerate(gene_group_set):
                try:
                    if math.isnan(float(group)):
                        no_groups = True
                        gene_group_set[i] = ' '
                except ValueError:
                    pass
            if not no_groups:
                for g,c,m in zip(gene_group_set, self.cell_color_list[0:len(gene_group_set)],self.markers[0:len(gene_group_set)]):
                    self.color_dict_genes[g] =[c,m]
            else:
                self.color_dict_genes[' '] =['b','o']

        self.cell_color_map(args)
        self.gene_color_map(args)

        self.log2_df_cell.to_csv(os.path.join(self.new_filepath, 'log2_matrix_after_filtering.txt'), sep= '\t', index_label=0)
        self.data_by_cell.to_csv(os.path.join(self.new_filepath, 'count_matrix_after_filtering.txt'), sep= '\t', index_label=0)

    def threshold_genes(self, by_cell, row_sum=1, cell_sum=1):
        return by_cell.loc[by_cell.sum(1) > row_sum, by_cell.sum(0) > cell_sum]

    def index_to_label(self,index):
        """Convert a pandas index or multiindex to an axis label."""
        if isinstance(index, pd.MultiIndex):
            return "-".join(map(str, index.names))
        else:
            return index.name

    def make_new_matrix_gene(self, gene_list_or_file):
        from itertools import compress
        if isinstance(gene_list_or_file,str):
            gene_df = pd.read_csv(open(gene_list_or_file,'rU'), sep=None, engine='python')
            try:
                self.short_gene_list = list(set(gene_df['GeneID'].tolist()))
                if self.exclude_list:

                    self.short_gene_list = [g for g in self.short_gene_list if g not in self.exclude_list and g in self.gene_list]
                else:
                    self.short_gene_list = gene_list
            except KeyError:
                sys.exit("Error: Please provide Gene list file with 'GeneID' as header.")
            try:
                if self.exclude_list:
                    mask = [True if g not in self.exclude_list and g in self.gene_list else False for g in self.short_gene_list]
                    self.gene_group_list = list(compress(gene_df['GroupID'].tolist(), mask))
                else:
                    self.gene_group_list = gene_df['GroupID'].tolist()
            except KeyError:
                sys.exit("Error: Please provide Gene list file with 'GroupID' as header.")
            try:
                self.short_gene_matrix_cell = self.threshold_genes(self.data_by_cell.loc[self.short_gene_list,:])

            except KeyError as error_gene:
                cause1 = error_gene.args[0].strip(' not in index')
                cause = [v.strip('\n\' ') for v in cause1.strip('[]').split(' ')]
                absent_gene = cause
                if args.verbose:
                    print(' '.join(absent_gene)+' not in matrix file.')
                new_list = [x for x in gene_list if x not in absent_gene]
                self.short_gene_matrix_cell = self.data_by_cell.loc[new_list,:]

        elif isinstance(gene_list_or_file,list):
            if self.exclude_list == []:
                self.short_gene_list = gene_list_or_file
            else:
                self.short_gene_list = [g for g in gene_list_or_file if g not in self.exclude_list and g in self.gene_list]
            try:
                self.short_gene_matrix_cell = self.threshold_genes(self.data_by_cell.loc[self.short_gene_list,:])
            except KeyError as error_gene:
                cause = error_gene.args[0]
                absent_gene = cause.split('\'')[1]
                if args.verbose:
                    print(absent_gene+' not in matrix file.')
                new_list = [x for x in gene_list if x not in [absent_gene]]
                gmatrix_df = org_matrix_by_gene[new_list]
        else:
            sys.exit("Error: gene list must be filepath or a list.")

    def make_new_matrix_cell(self):
        cell_df = pd.read_csv(open(self.cell_list_filepath,'rU'), sep=None, engine='python')
        cell_list_new = list(set([cell.strip('\n') for cell in cell_df['SampleID'].tolist()]))
        cell_list_old = self.data_by_cell.columns.tolist()
        overlap = [c for c in cell_list_new if c in cell_list_old]
        not_in_matrix = [c for c in cell_list_new if c not in cell_list_old]
        if not_in_matrix != []:
            print('These cells were in the cell list provided, but not found in the matrix provided:')
            print(not_in_matrix)
        self.data_by_cell = self.threshold_genes(self.data_by_cell.ix[:,overlap])
        self.cell_list = self.data_by_cell.columns.tolist()






    def find_top_common_genes(self):
        top_common_list = []
        count = 0
        done = False
        log2_df_by_gene = self.log2_df_cell.transpose()
        log2_df2_gene = log2_df_by_gene.apply(pd.to_numeric,errors='coerce')
        log_mean = log2_df2_gene.mean(axis=0).sort_values(ascending=False)
        try:
            log2_sorted_gene = log2_df_by_gene.reindex_axis(log2_df_by_gene.mean(axis=0).sort_values(ascending=False).index, axis=1)
        except ValueError:
            overlap_list = [item for item, count in collections.Counter(self.log2_df_cell.index).items() if count > 1]
            print(overlap_list, len(overlap_list))
            sys.exit('Error: Duplicate GeneIDs are present.')
        for gene in log2_sorted_gene.columns.tolist():
            if sum(genes < 1 for genes in log2_df_by_gene[gene])<6:
                if count < 20:
                    count+=1
                    top_common_list.append(gene)
            if count == 20:
                done = True
                break
        if done:
            return log2_df_by_gene[top_common_list].transpose()
        else:
            return [0]

    def log2_outlierfilter(self, args):
        if not args.already_log2:
            self.log2_df_cell = np.log2(self.data_by_cell+1)
            if isinstance(self.short_gene_matrix_cell, pd.DataFrame):
                self.log2_df_cell_gene_restricted = np.log2(self.short_gene_matrix_cell+1)
            else:
                self.log2_df_cell_gene_restricted = []
        else:
            self.log2_df_cell = self.data_by_cell
            if isinstance(self.short_gene_matrix_cell, pd.DataFrame):
                self.log2_df_cell_gene_restricted = self.short_gene_matrix_cell
            else:
                self.log2_df_cell_gene_restricted = []

        top_log2 = self.find_top_common_genes()
        if all(top_log2) != 0:
            log2_df= self.log2_df_cell.apply(pd.to_numeric,errors='coerce')
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

            filtered_df_by_cell = self.data_by_cell.ix[:,keep_col]
            if isinstance(self.short_gene_matrix_cell, pd.DataFrame):
                filtered_df_by_cell_gene_restricted = self.short_gene_matrix_cell.ix[:,keep_col]

            if not args.already_log2:
                filtered_df_by_cell +=1
                filtered_log2 = np.log2(filtered_df_by_cell[filtered_df_by_cell>0])
                if isinstance(self.short_gene_matrix_cell, pd.DataFrame):
                    filtered_df_by_cell_gene_restricted += 1
                    filtered_df_by_cell_gene_restricted_log2 = np.log2(filtered_df_by_cell_gene_restricted[filtered_df_by_cell_gene_restricted>0])

            else:
                filtered_log2 = filtered_df_by_cell[filtered_df_by_cell>0]
                if isinstance(self.short_gene_matrix_cell, pd.DataFrame):
                    filtered_df_by_cell_gene_restricted_log2 = filtered_df_by_cell_gene_restricted[filtered_df_by_cell_gene_restricted>0]
            fig1, ax1 = plt.subplots(figsize=(35,15))
            sns.boxplot(data=filtered_log2, whis= .75, notch=True, ax=ax1)
            sns.stripplot(x=filtered_log2.columns.values, y=filtered_log2.mean(axis=0), size=4, jitter=True, edgecolor="gray", ax=ax1)
            xtickNames = plt.setp(ax1, xticklabels=xticks)
            plt.setp(xtickNames, rotation=90, fontsize=9)
            plt.savefig(os.path.join(self.new_filepath, 'top_log2_genes_plot.pdf'), bbox_inches='tight')
            plt.clf()
            fig2, ax2 = plt.subplots(figsize=(12,10))
            sns.distplot(filtered_log2.mean(), color='r', label='After filtering', ax=ax2)
            sns.distplot(self.log2_df_cell.mean(), color='b', hist=False, label='Before filtering', ax=ax2)
            plt.savefig(os.path.join(self.new_filepath, 'distribution_log2_genes_plot.pdf'), bbox_inches='tight')

            self.log2_df_cell = filtered_log2
            if isinstance(self.short_gene_matrix_cell, pd.DataFrame):
                self.log2_df_cell_gene_restricted = self.threshold_genes(filtered_df_by_cell_gene_restricted_log2)
            else:
                self.log2_df_cell_gene_restricted = []
        else:
            if args.verbose:
                print("no common genes found")

    '''takes cell groups and creates dictionay 'label_map' that has attached color and marker
    if the same cell is assigned to multiple groups it assigns it to the first groupID
    '''
    def cell_color_map(self, args):
        if self.cell_list_filepath:
            cell_groups_df = pd.read_csv(open(os.path.join(os.path.dirname(args.filepath), self.cell_list_filepath),'rU'), sep=None, engine='python')
            cell_list_1 = list(set(cell_groups_df['SampleID'].tolist()))
            if len(cell_groups_df['SampleID']) == len(cell_groups_df['GroupID']):
                group_seen = []
                self.cell_label_map = {}
                cells_seen = []
                for cell, group in list(set(zip(cell_groups_df['SampleID'].tolist(), cell_groups_df['GroupID'].tolist()))):
                    if cell not in cells_seen:
                        group_count = self.cell_group_set.index(group)
                        self.cell_label_map[cell] = (self.color_dict_cells[group][0],self.color_dict_cells[group][1],group)
                        cells_seen.append(cell)
                non_group_cells = [c for c in cell_list_1 if c not in cells_seen]
                if non_group_cells != []:
                    for cell in non_group_cells:
                        self.cell_label_map[cell] = (self.color_list[group_count+1],self.markers[group_count+1],'No_Group')
            else:
                self.cell_label_map = False
        else:
            self.cell_label_map = False

    #takes cell groups and creates dictionay 'label_map' that has attached color and marker
    def gene_color_map(self, args):

        if self.gene_list_filepath:
            gene_df1 = pd.read_csv(open(os.path.join(os.path.dirname(args.filepath), self.gene_list_filepath),'rU'), sep=None, engine='python')
            if self.exclude_list != []:
                gene_df = gene_df1.ix[~gene_df1['GeneID'].isin(self.exclude_list)]
            gene_list_1 = [g for g in list(set(gene_df['GeneID'].tolist())) if g in self.gene_list]
            if len(gene_df['GeneID']) == len(gene_df['GroupID']):
                self.gene_label_map ={}
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
                        self.gene_label_map[gene] = (self.color_dict_genes[str(group)][0],self.color_dict_genes[str(group)][1],str(group))
                        genes_seen.append(gene)
                non_group_genes = [g for g in gene_list_1 if g not in genes_seen]
                if non_group_genes != []:
                    for cell in non_group_genes:
                        self.gene_label_map[gene] = (self.color_list[group_pos+1],self.markers[group_pos+1],'No_ID')
            else:
                self.gene_label_map = False
        else:
            self.gene_label_map = False
