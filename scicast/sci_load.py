
class Sci_load(object):
    def __init__(self):
        pass


    def load_options(self, all_options_dict):

        self.filepath = all_options_dict['filepath']
        self.base_name = all_options_dict['output_name']
        self.method = all_options_dict['method']
        self.metric = all_options_dict['metric']
        self.gene_number = all_options_dict['gene_number']
        self.depth_of_clustering = all_options_dict['depth']
        self.cell_list_filename = all_options_dict['cell_file']
        self.gene_list_filename = all_options_dict['gene_file']
        self.z_direction = all_options_dict['zdir']
        self.qgraph_plot = all_options_dict['qgraph']
        self.exclude_genes = all_options_dict['exclude_genes']
        self.kmeans_cluster_range = all_options_dict['kmeans_cluster_range']
        self.color_cells = all_options_dict['color_cells']
        self.color_genes = all_options_dict['color_genes']
        self.genes_corr = all_options_dict['genes_corr']
        self.test_clust_stability = all_options_dict['test_clust_stability']
        self.annotate_gene_subset = all_options_dict['annotate_gene_subset']
        self.no_heatmaps = all_options_dict["Don't Run Heatmaps"]
        self.no_corr = all_options_dict["Don't Run Correlation"]
        self.verbose = all_options_dict["Verbose"]
        self.group_sig_test = all_options_dict["Test Significance by Groups (User Defined)"]
        self.all_sig = all_options_dict["Test Significance by Unbiased Clusters"]
        self.limit_cells= all_options_dict["Exclude Cells Not in User Cell Groups"]
        self.add_ellipse= all_options_dict["Add Ellipse"]
        self.annotate_cell_pca= all_options_dict["Add Cell Names to PCA"]
        self.sig_unique= all_options_dict["Display Only Unique Signifcant Genes"]
        self.kmeans_sig_test = all_options_dict["Run Significance Test for kmeans clusters"]
        self.already_log2 = all_options_dict["Input Matrix is already log2"]
        self.use_TSNE = all_options_dict["use t-SNE (for kmeans clustering)"]
        self.image_format = all_options_dict["image_format"]
        return(self)
