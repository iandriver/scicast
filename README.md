# scicast
Single Cell Iterative Clustering and Significance Testing
------
**This is still in development. If it is useful please let me know.  I will be updating with a proper manual and vignette with sample data set as I get closer to publishing. -Ian (ian.driver@ucsf.edu)**


Requires
--------
- Python 2.7 or 3.4+ (preferred)

- [numpy](http://www.numpy.org/)

- [scipy](http://www.scipy.org/)

- [matplotlib](http://matplotlib.org/)

- [pandas >= 0.19.0](http://pandas.pydata.org/)

- [seaborn >= 0.7.1](http://seaborn.pydata.org/)

- [R >3.3](https://www.r-project.org/)
        -must have [qgraph](http://sachaepskamp.com/qgraph) installed in R for qgraph_plot to work

            install.packages("qgraph")

- [rpy2](http://rpy2.bitbucket.org/)

- [fastcluster](https://pypi.python.org/pypi/fastcluster)

pip or brew or conda install any that are missing, prior to scicast installation.

Note on rpy2
--------------
Install or update R, and then pip install rpy2. Otherwise it won't link properly.

Installation
------------

To install the released version, just type:

    pip install scicast

  Example command
--------------
Simple Command:
There is now a gui (tkinter based):
```python
          scicast -gui
```
will launch the gui interface:

scicast GUI window
-------------
![scicast GUI Window Image](scicast_gui.png)

otherwise:

```python
          scicast -f path/to/DESeq_krasnow_185_At2_rsem_matrix_norm.txt -n krasnow_example_data
```
Command with cell and gene files:
```python
          scicast -f /Path/to/gene_cell_matrix -n my_gene_cell_matrix_analyis -v -cell_group sc_cell_groups.csv -gene_list sc_gene_groups.csv -depth 100 -z 0 -genes_corr Myc,Nanog,Vegf -method ward -metric seuclidean -annotate_gene_subset gene_annotation_subset.csv -add_ellipse
```
Example Dataset
-------------
For the soon to be vignette and for trying out scicast an example single cell dataset is provided: DESeq_krasnow_185_At2_rsem_matrix_norm.txt

From the paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4145853/
Raw sequencing was mapped and aligned to GRCm38 using rsem and normalized using Deseq2.
Only E18.5 and Adult At2 cells are included.

=======

Usage
-------
usage: scicast

      -gui        To run gui instead of command line. (default: False)

else:

                  [-h] -f FILEPATH [-n BASE_NAME] [-g GENE_NUMBER] [-v]
                  [-cell_group CELL_LIST_FILENAME]
                  [-gene_list GENE_LIST_FILENAME] [-group_sig_test]
                  [-stability TEST_CLUST_STABILITY] [-metric METRIC]
                  [-method METHOD] [-depth DEPTH_OF_CLUSTERING]
                  [-genes_corr GENES_CORR] [-all_sig] [-already_log2]
                  [-z Z_DIRECTION] [-no_corr] [-no_heatmaps]
                  [-exclude_genes EXCLUDE_GENES] [-limit_cells] [-add_ellipse]
                  [-annotate_gene_subset ANNOTATE_GENE_SUBSET]
                  [-annotate_cell_pca] [-color_cells COLOR_CELLS]
                  [-color_genes COLOR_GENES] [-qgraph_plot QGRAPH_PLOT]
                  [-sig_unique] [-kmeans_cluster_range KMEANS_CLUSTER_RANGE]
                  [-kmeans_sig_test]

optional arguments:

  -h, --help            show this help message and exit

  -f FILEPATH, -filepath FILEPATH

                        Filepath to cell-gene matrix file (default: None)

  -n BASE_NAME, -name BASE_NAME

                        The base name for files that will be generated.
                        (default: scicast_analysis)

  -g GENE_NUMBER, -gene_number GENE_NUMBER

                        The number of genes that will be selected from top pca
                        genes on each iteration. (default: 200)

  -v, -verbose          

                        Print more (default: False)

  -cell_group CELL_LIST_FILENAME

                        Provide a filename with two columns: 'SampleID' and
                        'GroupID'. If GroupID is empty the SampleID list will
                        be used to restrict cells prior to analysis. (default:
                        False)

  -gene_list GENE_LIST_FILENAME

                        Path to a file with a two columns with headers
                        'GeneID' and 'GroupID' (Singular format). GeneID list
                        will be used to create a new matrix file with only
                        those genes included. (default: False)

  -group_sig_test       

                        If cell groups are provided, perform significance
                        testing between all groups (independent of any cluster
                        groups). (default: False)

  -stability TEST_CLUST_STABILITY

                        Provide a number of iterations to test how stable
                        clustering is as the number of top PCA genes changes
                        from 100-1000. Output will be clustering heatmaps for
                        each iteration a summary of changes as gene number
                        varies. (default: 0)

  -metric METRIC        

                        The distance metric to use. The distance function can
                        be: braycurtis, canberra, chebyshev, cityblock,
                        correlation, cosine, dice, euclidean, hamming,
                        jaccard, kulsinski, mahalanobis, matching, minkowski,
                        rogerstanimoto, russellrao, seuclidean, sokalmichener,
                        sokalsneath, sqeuclidean, yule. (default: euclidean)

  -method METHOD        

                        The linkage method to use. The linkage algorithm can
                        be: single, complete, average, weighted, centroid,
                        median or ward. (default: weighted)

  -depth DEPTH_OF_CLUSTERING

                        The size in cell number at which sub-clustering
                        analysis will stop clustering, pca and correlation
                        analysis. (default: 20)

  -genes_corr GENES_CORR

                        If you want to look for correlation on a specific gene
                        or set of genes enter them as a comma seperated list
                        i.e. 'Gapdh,Actb'. (default: )

  -all_sig              

                        Do significance testing on all hierarchical clustering
                        groups. The minimum number of cells in a group is set
                        by --depth. (default: False)

  -already_log2         

                        If the gene matrix is already log2 transformed i.e.
                        rld or vsd DESeq2, this will disable log2 transform.
                        (default: False)

  -z Z_DIRECTION        

                        Either 0 (rows) or 1 (columns) or None. Whether or not
                        to calculate z-scores for the rows or the columns. Z
                        scores are: z = (x - mean)/std, so values in each row
                        (column) will get the mean of the row (column)
                        subtracted, then divided by the standard deviation of
                        the row (column). This ensures that each row (column)
                        has mean of 0 and variance of 1. (default: 0)

  -no_corr              

                        Don't run correlation search. Default is on. (default:
                        False)

  -no_heatmaps          

                        Don't run heatmaps and pca. Default is will generate
                        heatmaps and pca. (default: False)

  -exclude_genes EXCLUDE_GENES

                        Filepath to list of genes to be excluded from analysis
                        (at all levels of analysis). Header must be 'GeneID'.
                        (default: False)

  -limit_cells          

                        With cell group file, will exclude all other cells
                        from analysis (at all levels of analysis). Header must
                        be 'SampleID'. (default: False)

  -add_ellipse          

                        When present colored ellipses will be added to cell
                        and gene PCA plots. Must provide gene and/or cell
                        groups. (default: False)

  -annotate_gene_subset ANNOTATE_GENE_SUBSET

                        Provide path or filename (if in same file) to file
                        with genes to be annotated on gene PCA. Must have
                        'GeneID' header. (default: False)

  -annotate_cell_pca    

                        Option will annotate cell PCA with cell names. Default
                        is off (False). (default: False)

  -color_cells COLOR_CELLS

                        Follow with cell group names and pairs (matching cell
                        group names). Forces color and marker style on group.
                        Just color or color and marker can be provided.
                        Example: Group1,color,marker Group2,red,v (default:
                        False)

  -color_genes COLOR_GENES

                        Can be 'same' if names are the same as cell groups.
                        Follow with gene group names and pairs (matching gene
                        group names). Forces color and marker style on group.
                        Just color or color and marker can be provided.
                        Example: Group1,color,marker Group2,red,v (default:
                        False)

  -qgraph_plot QGRAPH_PLOT

                        Can be 'cell' 'gene' or 'both' will plot qgraph gene
                        and/or cell network and pca correlation network plot
                        to pdf. Requires both gene and cell groups to be
                        provided. (default: False)

  -sig_unique           

                        group_sig_test will by default find unique sets of
                        significant genes, so that the whole list has no
                        duplicates. This switches the top significant genes
                        amoung groups to allow repeats. (default: True)

  -kmeans_cluster_range KMEANS_CLUSTER_RANGE

                        range of cluster sizes to run kmeans clustering
                        (inclusive) (default: 2,4)

  -kmeans_sig_test      

                        Run significance test between groups defined by kmeans
                        clusters. (default: False)
