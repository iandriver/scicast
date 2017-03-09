

try:
    import tkinter as tk
    from tkinter.filedialog import askopenfilename

    class Window(tk.Tk):
        def __init__(self):
            tk.Tk.__init__(self)
            self.title("scicast")
            self.path = tk.StringVar()
            self.cell_path = tk.StringVar()
            self.gene_path = tk.StringVar()
            self.gene_label_path = tk.StringVar()
            self.exclude_gene_path = tk.StringVar()
            self.asset = tk.StringVar()
            self.gene_number = tk.IntVar(value=200)
            self.depth_number = tk.IntVar(value=20)
            self.kmeans_cluster_range = tk.StringVar(value='2,4')
            self.color_cells = tk.StringVar()
            self.color_genes = tk.StringVar()
            self.test_clust_stability = tk.IntVar(value=0)
            self.genes_corr = tk.StringVar()
            self.annotate_gene_subset = tk.StringVar()

            #type or choose gene matrix file
            dir_label = tk.Label(self, text="Browse or type path to gene cell matrix file:")
            path_entry = tk.Entry(self, textvariable=self.path, width=40)
            browse_button = tk.Button(self, text="Browse for gene cell matrix file", command=self.browse)

            #type or choose cell group file
            cell_label = tk.Label(self, text="Browse or type path to cell group file:")
            cell_path_entry = tk.Entry(self, textvariable=self.cell_path, width=40)
            cell_browse_button = tk.Button(self, text="Browse for cell group file", command=self.browse_cellp)

            #type or choose gene group file
            gene_label = tk.Label(self, text="Browse or type path to gene group file:")
            gene_path_entry = tk.Entry(self, textvariable=self.gene_path, width=40)
            gene_browse_button = tk.Button(self, text="Browse for gene group file", command=self.browse_genep)

            #type or choose file of genes to exclude from all analysis
            exclude_gene_label = tk.Label(self, text="Browse or type path to file with genes to exclude from analysis (i.e cell cycle):")
            exclude_gene_path_entry = tk.Entry(self, textvariable=self.exclude_gene_path, width=40)
            exclude_gene_browse_button = tk.Button(self, text="Browse for exclude genes file", command=self.browse_excludeg)

            #type or choose file of genes which will be labeled in gene PCA plot
            annotate_gene_subset_label = tk.Label(self, text="Annotate only these genes in gene PCA:")
            annotate_gene_subset_path_entry = tk.Entry(self, textvariable=self.annotate_gene_subset, width=40)
            annotate_gene_subset_browse_button = tk.Button(self, text="Browse gene annotation file", command=self.browse_annotateg)

            #define file extensions
            self.file_opt = options = {}
            options['defaultextension'] = '.txt'
            options['filetypes'] = [('all files', '.*'), ('text files', '.txt'),('csv files', '.csv'), ('cufflinks counts', '.count_table'), ('cufflinks fpkm', 'fpkm_table'), ('gene matrix', '.matrix')]

            #setup metric menu options
            self.metric_menu_var = tk.StringVar()
            self.metric_menu_var.set("seuclidean")
            metric_menu_label = tk.Label(self, text="Choose Metric:")
            metric_option_menu = tk.OptionMenu(self, self.metric_menu_var, 'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule')

            #setup method option menu
            self.method_menu_var = tk.StringVar()
            self.method_menu_var.set("ward")
            method_menu_label = tk.Label(self, text="Choose Method:")
            method_option_menu = tk.OptionMenu(self, self.method_menu_var, 'single', 'complete', 'average', 'weighted', 'centroid', 'median')

            #setup qgraph option menu
            self.qgraph_menu_var = tk.StringVar()
            self.qgraph_menu_var.set("none")
            qgraph_menu_label = tk.Label(self, text="Choose which qgraph networks to generate:")
            qgraph_option_menu = tk.OptionMenu(self, self.qgraph_menu_var, 'gene','cell','both')

            #setup image format selection menu
            self.image_format_menu_var = tk.StringVar()
            self.image_format_menu_var.set("pdf")
            image_format_menu_label = tk.Label(self, text="Select image format for output files:")
            image_format_option_menu = tk.OptionMenu(self, self.image_format_menu_var, 'tif', 'png', 'jpeg')

            #setup z-direction option menu
            self.zdir_menu_var = tk.IntVar()
            self.zdir_menu_var.set(0)
            zdir_menu_label = tk.Label(self, text="Choose z:")
            zdir_option_menu = tk.OptionMenu(self, self.zdir_menu_var, 1,'None')



            self.flags = ["Don't Run Heatmaps","Don't Run Correlation", "Verbose", "Test Significance by Groups (User Defined)", "Test Significance by Unbiased Clusters", "Exclude Cells Not in User Cell Groups", "Add Ellipse", "Add Cell Names to PCA", "Display Only Unique Signifcant Genes", "Run Significance Test for kmeans clusters", "Input Matrix is already log2", "use t-SNE (for kmeans clustering)"]
            self.variables = []



            asset_label = tk.Label(self, text="Output File Name:")
            asset_entry = tk.Entry(self, textvariable=self.asset, width=40)

            gene_number_label = tk.Label(self, text="Number of genes to include")
            gene_number_entry = tk.Entry(self, textvariable=self.gene_number, width=10)

            kmeans_range_label = tk.Label(self, text="Range of cluster for kmeans (inclusive):")
            kmeans_range_entry = tk.Entry(self, textvariable=self.kmeans_cluster_range, width=10)

            depth_number_label = tk.Label(self, text="Depth at which subclustering will stop")
            depth_number_entry = tk.Entry(self, textvariable=self.depth_number, width=10)

            color_cells_label = tk.Label(self, text="Provide specific colors and markers for each cell group.")
            color_cells_entry = tk.Entry(self, textvariable=self.color_cells, width=20)

            color_genes_label = tk.Label(self, text="Provide specific colors and markers for each gene group.")
            color_genes_entry = tk.Entry(self, textvariable=self.color_cells, width=20)

            test_clust_stability_label = tk.Label(self, text="Number of iterations to test cluster stability with varying gene numbers.")
            test_clust_stability_entry = tk.Entry(self, textvariable=self.test_clust_stability, width=10)

            genes_corr_label = tk.Label(self, text="Comma seperated list of genes to add to correlation search.")
            genes_corr_entry = tk.Entry(self, textvariable=self.genes_corr, width=20)


            create_button = tk.Button(self, text="Run scicast", command=self.genAsset)

            dir_label.grid(row=1, column=1, columnspan=2, sticky='w')
            path_entry.grid(row=2, column=1, columnspan=2, sticky='w')
            browse_button.grid(row=3, column=1, columnspan=2, sticky='w')

            cell_label.grid(row=4, column=1, columnspan=2, sticky='w')
            cell_path_entry.grid(row=5, column=1, columnspan=2, sticky='w')
            cell_browse_button.grid(row=6, column=1, columnspan=2, sticky='w')

            gene_label.grid(row=7, column=1, columnspan=2, sticky='w')
            gene_path_entry.grid(row=8, column=1, columnspan=2, sticky='w')
            gene_browse_button.grid(row=9, column=1, columnspan=2, sticky='w')



            exclude_gene_label.grid(row=10, column=1, columnspan=2, sticky='w')
            exclude_gene_path_entry.grid(row=11, column=1, columnspan=2, sticky='w')
            exclude_gene_browse_button.grid(row=12, column=1, columnspan=2, sticky='w')

            annotate_gene_subset_label.grid(row=9, column=4, columnspan=2, sticky='w')
            annotate_gene_subset_path_entry.grid(row=10, column=4, columnspan=2, sticky='w')
            annotate_gene_subset_browse_button.grid(row=11, column=4, columnspan=2, sticky='w')

            gene_number_label.grid(row=13, column=1, columnspan=2, sticky='w')
            gene_number_entry.grid(row=14, column=1, columnspan=2, sticky='w')

            depth_number_label.grid(row=15, column=1, columnspan=2, sticky='w')
            depth_number_entry.grid(row=16, column=1, columnspan=2, sticky='w')

            for i, flag in enumerate(self.flags):
                var = tk.BooleanVar()
                tk.Checkbutton(self, text=flag, variable=var).grid(row=1+i, column=3, columnspan=1, sticky='w')
                self.variables.append(var)

            metric_menu_label.grid(row=2+len(self.flags), column=3, columnspan=1, sticky='w')
            metric_option_menu.grid(row=3+len(self.flags), column=3, columnspan=1, sticky='w')

            method_menu_label.grid(row=4+len(self.flags), column=3, columnspan=1, sticky='w')
            method_option_menu.grid(row=5+len(self.flags), column=3, columnspan=1, sticky='w')

            qgraph_menu_label.grid(row=15, column=4, columnspan=1, sticky='w')
            qgraph_option_menu.grid(row=16, column=4, columnspan=1, sticky='w')

            image_format_menu_label.grid(row=17, column=4, columnspan=1, sticky='e')
            image_format_option_menu.grid(row=18, column=4, columnspan=1, sticky='e')

            zdir_menu_label.grid(row=17, column=4, columnspan=1, sticky='w')
            zdir_option_menu.grid(row=18, column=4, columnspan=1, sticky='w')

            kmeans_range_label.grid(row=12, column=4, columnspan=1, sticky='w')
            kmeans_range_entry.grid(row=13, column=4, columnspan=1, sticky='w')

            color_cells_label.grid(row=1, column=4, columnspan=1, sticky='w')
            color_cells_entry.grid(row=2, column=4, columnspan=1, sticky='w')

            color_genes_label.grid(row=3, column=4, columnspan=1, sticky='w')
            color_genes_entry.grid(row=4, column=4, columnspan=1, sticky='w')

            test_clust_stability_label.grid(row=5, column=4, columnspan=1, sticky='w')
            test_clust_stability_entry.grid(row=6, column=4, columnspan=1, sticky='w')

            genes_corr_label.grid(row=7, column=4, columnspan=1, sticky='w')
            genes_corr_entry.grid(row=8, column=4, columnspan=1, sticky='w')

            asset_label.grid(row=17, column=1, columnspan=1, sticky='w')
            asset_entry.grid(row=18, column=1, columnspan=1, sticky='w')

            create_button.grid(row=24, column=2, columnspan=2)


        def browse(self):
            file_path= askopenfilename(**self.file_opt)
            if file_path:
                self.path.set(file_path)
        def browse_cellp(self):
            file_path= askopenfilename(**self.file_opt)
            if file_path:
                self.cell_path.set(file_path)
        def browse_genep(self):
            file_path= askopenfilename(**self.file_opt)
            if file_path:
                self.gene_path.set(file_path)
        def browse_excludeg(self):
            file_path= askopenfilename(**self.file_opt)
            if file_path:
                self.exclude_gene_path.set(file_path)
        def browse_annotateg(self):
            file_path= askopenfilename(**self.file_opt)
            if file_path:
                self.annotate_gene_subset.set(file_path)

        def genAsset(self):
            all_options_dict = {}
            asset_path = self.path.get()
            asset_name = self.asset.get()
            asset_metric_menu_option = self.metric_menu_var.get()
            asset_method_menu_option = self.method_menu_var.get()
            asset_gene_number = self.gene_number.get()
            asset_depth = self.depth_number.get()
            asset_cell_path = self.cell_path.get()
            asset_gene_path = self.gene_path.get()
            asset_zdir = self.zdir_menu_var.get()
            asset_qgraph = self.qgraph_menu_var.get()
            asset_image_format = self.image_format_menu_var.get()
            asset_kmeans_cluster_range = self.kmeans_cluster_range.get()
            asset_exclude_gene_path = self.exclude_gene_path.get()
            asset_color_cells = self.color_cells.get()
            asset_color_genes = self.color_genes.get()
            asset_test_clust_stability = self.test_clust_stability.get()
            asset_genes_corr = self.genes_corr.get()
            asset_annotate_gene_subset = self.annotate_gene_subset.get()
            for var, flag in zip(self.variables, self.flags):
                all_options_dict[flag] = var.get()
            all_options_dict['filepath'] = asset_path
            all_options_dict['output_name'] = asset_name
            all_options_dict['method'] = asset_method_menu_option
            all_options_dict['metric'] =asset_metric_menu_option
            all_options_dict['gene_number'] =asset_gene_number
            all_options_dict['depth'] = asset_depth
            all_options_dict['cell_file'] = asset_cell_path
            all_options_dict['gene_file'] = asset_gene_path
            all_options_dict['zdir'] = asset_zdir
            all_options_dict['qgraph'] = asset_qgraph
            all_options_dict['image_format'] = asset_image_format
            all_options_dict['kmeans_cluster_range'] = asset_kmeans_cluster_range
            all_options_dict['exclude_genes'] = asset_exclude_gene_path
            all_options_dict['color_cells'] = asset_color_cells
            all_options_dict['color_genes'] = asset_color_genes
            all_options_dict['test_clust_stability'] = asset_test_clust_stability
            all_options_dict['genes_corr'] = asset_genes_corr
            all_options_dict['annotate_gene_subset'] = asset_annotate_gene_subset
            self.all_dict = all_options_dict
            self.destroy()
except ImportError:
    import Tkinter as tk
    import tkFileDialog

    class Window(tk.Frame):
        def __init__(self):
            tk.Frame.__init__(self)
            #self.title("scicast")
            self.path = tk.StringVar()
            self.cell_path = tk.StringVar()
            self.gene_path = tk.StringVar()
            self.gene_label_path = tk.StringVar()
            self.exclude_gene_path = tk.StringVar()
            self.asset = tk.StringVar()
            self.gene_number = tk.IntVar(value=200)
            self.depth_number = tk.IntVar(value=20)
            self.kmeans_cluster_range = tk.StringVar(value='2,4')
            self.color_cells = tk.StringVar()
            self.color_genes = tk.StringVar()
            self.test_clust_stability = tk.IntVar(value=0)
            self.genes_corr = tk.StringVar()
            self.annotate_gene_subset = tk.StringVar()

            #type or choose gene matrix file
            dir_label = tk.Label(self, text="Browse or type path to gene cell matrix file:")
            path_entry = tk.Entry(self, textvariable=self.path, width=40)
            browse_button = tk.Button(self, text="Browse for gene cell matrix file", command=self.browse)

            #type or choose cell group file
            cell_label = tk.Label(self, text="Browse or type path to cell group file:")
            cell_path_entry = tk.Entry(self, textvariable=self.cell_path, width=40)
            cell_browse_button = tk.Button(self, text="Browse for cell group file", command=self.browse_cellp)

            #type or choose gene group file
            gene_label = tk.Label(self, text="Browse or type path to gene group file:")
            gene_path_entry = tk.Entry(self, textvariable=self.gene_path, width=40)
            gene_browse_button = tk.Button(self, text="Browse for gene group file", command=self.browse_genep)

            #type or choose file of genes to exclude from all analysis
            exclude_gene_label = tk.Label(self, text="Browse or type path to file with genes to exclude from analysis (i.e cell cycle):")
            exclude_gene_path_entry = tk.Entry(self, textvariable=self.exclude_gene_path, width=40)
            exclude_gene_browse_button = tk.Button(self, text="Browse for exclude genes file", command=self.browse_excludeg)

            #type or choose file of genes which will be labeled in gene PCA plot
            annotate_gene_subset_label = tk.Label(self, text="Annotate only these genes in gene PCA:")
            annotate_gene_subset_path_entry = tk.Entry(self, textvariable=self.annotate_gene_subset, width=40)
            annotate_gene_subset_browse_button = tk.Button(self, text="Browse gene annotation file", command=self.browse_annotateg)

            #define file extensions
            self.file_opt = options = {}
            options['defaultextension'] = '.txt'
            options['filetypes'] = [('all files', '.*'), ('text files', '.txt'),('csv files', '.csv')]

            #setup metric menu options
            self.metric_menu_var = tk.StringVar()
            self.metric_menu_var.set("seuclidean")
            metric_menu_label = tk.Label(self, text="Choose Metric:")
            metric_option_menu = tk.OptionMenu(self, self.metric_menu_var, 'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule')

            #setup method option menu
            self.method_menu_var = tk.StringVar()
            self.method_menu_var.set("ward")
            method_menu_label = tk.Label(self, text="Choose Method:")
            method_option_menu = tk.OptionMenu(self, self.method_menu_var, 'single', 'complete', 'average', 'weighted', 'centroid', 'median')

            #setup qgraph option menu
            self.qgraph_menu_var = tk.StringVar()
            self.qgraph_menu_var.set("none")
            qgraph_menu_label = tk.Label(self, text="Choose which qgraph networks to generate:")
            qgraph_option_menu = tk.OptionMenu(self, self.qgraph_menu_var, 'gene','cell','both')

            #setup image format selection menu
            self.image_format_menu_var = tk.StringVar()
            self.image_format_menu_var.set("pdf")
            image_format_menu_label = tk.Label(self, text="Select image format for output files:")
            image_format_option_menu = tk.OptionMenu(self, self.image_format_menu_var, 'tif', 'png', 'jpeg')

            #setup z-direction option menu
            self.zdir_menu_var = tk.IntVar()
            self.zdir_menu_var.set(0)
            zdir_menu_label = tk.Label(self, text="Choose z:")
            zdir_option_menu = tk.OptionMenu(self, self.zdir_menu_var, 1,'None')

            self.flags = ["Don't Run Heatmaps","Don't Run Correlation", "Verbose", "Test Significance by Groups (User Defined)", "Test Significance by Unbiased Clusters", "Exclude Cells Not in User Cell Groups", "Add Ellipse", "Add Cell Names to PCA", "Display Only Unique Signifcant Genes", "Run Significance Test for kmeans clusters", "Input Matrix is already log2", "use t-SNE (for kmeans clustering)"]
            self.variables = []



            asset_label = tk.Label(self, text="Output File Name:")
            asset_entry = tk.Entry(self, textvariable=self.asset, width=40)

            gene_number_label = tk.Label(self, text="Number of genes to include")
            gene_number_entry = tk.Entry(self, textvariable=self.gene_number, width=10)

            kmeans_range_label = tk.Label(self, text="Range of cluster for kmeans (inclusive):")
            kmeans_range_entry = tk.Entry(self, textvariable=self.kmeans_cluster_range, width=10)

            depth_number_label = tk.Label(self, text="Depth at which subclustering will stop")
            depth_number_entry = tk.Entry(self, textvariable=self.depth_number, width=10)

            color_cells_label = tk.Label(self, text="Provide specific colors and markers for each cell group.")
            color_cells_entry = tk.Entry(self, textvariable=self.color_cells, width=20)

            color_genes_label = tk.Label(self, text="Provide specific colors and markers for each gene group.")
            color_genes_entry = tk.Entry(self, textvariable=self.color_cells, width=20)

            test_clust_stability_label = tk.Label(self, text="Number of iterations to test cluster stability with varying gene numbers.")
            test_clust_stability_entry = tk.Entry(self, textvariable=self.test_clust_stability, width=10)

            genes_corr_label = tk.Label(self, text="Comma seperated list of genes to add to correlation search.")
            genes_corr_entry = tk.Entry(self, textvariable=self.genes_corr, width=20)


            create_button = tk.Button(self, text="Run scicast", command=self.genAsset)

            dir_label.grid(row=1, column=1, columnspan=2, sticky='w')
            path_entry.grid(row=2, column=1, columnspan=2, sticky='w')
            browse_button.grid(row=3, column=1, columnspan=2, sticky='w')

            cell_label.grid(row=4, column=1, columnspan=2, sticky='w')
            cell_path_entry.grid(row=5, column=1, columnspan=2, sticky='w')
            cell_browse_button.grid(row=6, column=1, columnspan=2, sticky='w')

            gene_label.grid(row=7, column=1, columnspan=2, sticky='w')
            gene_path_entry.grid(row=8, column=1, columnspan=2, sticky='w')
            gene_browse_button.grid(row=9, column=1, columnspan=2, sticky='w')



            exclude_gene_label.grid(row=10, column=1, columnspan=2, sticky='w')
            exclude_gene_path_entry.grid(row=11, column=1, columnspan=2, sticky='w')
            exclude_gene_browse_button.grid(row=12, column=1, columnspan=2, sticky='w')

            annotate_gene_subset_label.grid(row=9, column=4, columnspan=2, sticky='w')
            annotate_gene_subset_path_entry.grid(row=10, column=4, columnspan=2, sticky='w')
            annotate_gene_subset_browse_button.grid(row=11, column=4, columnspan=2, sticky='w')

            gene_number_label.grid(row=13, column=1, columnspan=2, sticky='w')
            gene_number_entry.grid(row=14, column=1, columnspan=2, sticky='w')

            depth_number_label.grid(row=15, column=1, columnspan=2, sticky='w')
            depth_number_entry.grid(row=16, column=1, columnspan=2, sticky='w')

            for i, flag in enumerate(self.flags):
                var = tk.BooleanVar()
                tk.Checkbutton(self, text=flag, variable=var).grid(row=1+i, column=3, columnspan=1, sticky='w')
                self.variables.append(var)

            metric_menu_label.grid(row=2+len(self.flags), column=3, columnspan=1, sticky='w')
            metric_option_menu.grid(row=3+len(self.flags), column=3, columnspan=1, sticky='w')

            method_menu_label.grid(row=4+len(self.flags), column=3, columnspan=1, sticky='w')
            method_option_menu.grid(row=5+len(self.flags), column=3, columnspan=1, sticky='w')

            qgraph_menu_label.grid(row=15, column=4, columnspan=1, sticky='w')
            qgraph_option_menu.grid(row=16, column=4, columnspan=1, sticky='w')

            image_format_menu_label.grid(row=17, column=5, columnspan=1, sticky='w')
            image_format_option_menu.grid(row=18, column=5, columnspan=1, sticky='w')

            zdir_menu_label.grid(row=17, column=4, columnspan=1, sticky='w')
            zdir_option_menu.grid(row=18, column=4, columnspan=1, sticky='w')

            kmeans_range_label.grid(row=12, column=4, columnspan=1, sticky='w')
            kmeans_range_entry.grid(row=13, column=4, columnspan=1, sticky='w')

            color_cells_label.grid(row=1, column=4, columnspan=1, sticky='w')
            color_cells_entry.grid(row=2, column=4, columnspan=1, sticky='w')

            color_genes_label.grid(row=3, column=4, columnspan=1, sticky='w')
            color_genes_entry.grid(row=4, column=4, columnspan=1, sticky='w')

            test_clust_stability_label.grid(row=5, column=4, columnspan=1, sticky='w')
            test_clust_stability_entry.grid(row=6, column=4, columnspan=1, sticky='w')

            genes_corr_label.grid(row=7, column=4, columnspan=1, sticky='w')
            genes_corr_entry.grid(row=8, column=4, columnspan=1, sticky='w')

            asset_label.grid(row=17, column=1, columnspan=1, sticky='w')
            asset_entry.grid(row=18, column=1, columnspan=1, sticky='w')

            create_button.grid(row=24, column=2, columnspan=2)


        def browse(self):
            file_path= tkFileDialog.askopenfilename(**self.file_opt)
            if file_path:
                self.path.set(file_path)
        def browse_cellp(self):
            file_path= tkFileDialog.askopenfilename(**self.file_opt)
            if file_path:
                self.cell_path.set(file_path)
        def browse_genep(self):
            file_path= tkFileDialog.askopenfilename(**self.file_opt)
            if file_path:
                self.gene_path.set(file_path)
        def browse_excludeg(self):
            file_path= tkFileDialog.askopenfilename(**self.file_opt)
            if file_path:
                self.exclude_gene_path.set(file_path)
        def browse_annotateg(self):
            file_path= tkFileDialog.askopenfilename(**self.file_opt)
            if file_path:
                self.annotate_gene_subset.set(file_path)

        def genAsset(self):
            all_options_dict = {}
            asset_path = self.path.get()
            asset_name = self.asset.get()
            asset_metric_menu_option = self.metric_menu_var.get()
            asset_method_menu_option = self.method_menu_var.get()
            asset_gene_number = self.gene_number.get()
            asset_depth = self.depth_number.get()
            asset_cell_path = self.cell_path.get()
            asset_gene_path = self.gene_path.get()
            asset_zdir = self.zdir_menu_var.get()
            asset_qgraph = self.qgraph_menu_var.get()
            asset_image_format = self.image_format_menu_var.get()
            asset_kmeans_cluster_range = self.kmeans_cluster_range.get()
            asset_exclude_gene_path = self.exclude_gene_path.get()
            asset_color_cells = self.color_cells.get()
            asset_color_genes = self.color_genes.get()
            asset_test_clust_stability = self.test_clust_stability.get()
            asset_genes_corr = self.genes_corr.get()
            asset_annotate_gene_subset = self.annotate_gene_subset.get()
            for var, flag in zip(self.variables, self.flags):
                all_options_dict[flag] = var.get()
            all_options_dict['filepath'] = asset_path
            all_options_dict['output_name'] = asset_name
            all_options_dict['method'] = asset_method_menu_option
            all_options_dict['metric'] =asset_metric_menu_option
            all_options_dict['gene_number'] =asset_gene_number
            all_options_dict['depth'] = asset_depth
            all_options_dict['cell_file'] = asset_cell_path
            all_options_dict['gene_file'] = asset_gene_path
            all_options_dict['zdir'] = asset_zdir
            all_options_dict['qgraph'] = asset_qgraph
            all_options_dict['image_format'] = asset_image_format
            all_options_dict['kmeans_cluster_range'] = asset_kmeans_cluster_range
            all_options_dict['exclude_genes'] = asset_exclude_gene_path
            all_options_dict['color_cells'] = asset_color_cells
            all_options_dict['color_genes'] = asset_color_genes
            all_options_dict['test_clust_stability'] = asset_test_clust_stability
            all_options_dict['genes_corr'] = asset_genes_corr
            all_options_dict['annotate_gene_subset'] = asset_annotate_gene_subset
            self.all_dict = all_options_dict
            self.destroy()
