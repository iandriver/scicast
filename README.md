# SCICAST
Single Cell Iterative Clustering and Significance Testing

usage: cluster.py [-h] -f FILEPATH [-n BASE_NAME] [-g GENE_NUMBER] [-v]
                  [-cell_group CELL_GROUP_FILE] [-group_sig_test]
                  [-gene_list MAKE_GENE_MATRIX] [-cell_list MAKE_CELL_MATRIX]
                  [-stability TEST_CLUST_STABILITY] [--metric METRIC]
                  [--method METHOD]

optional arguments:
  -h, --help            show this help message and exit
  -f FILEPATH, --filepath FILEPATH
                        Filepath to cell-gene matrix file (default: None)
  -n BASE_NAME, --name BASE_NAME
                        The base name for files that will be generated
                        (default: )
  -g GENE_NUMBER, --gene_number GENE_NUMBER
                        The number of genes that will be selected from top pca
                        genes on each iteration. (default: 200)
  -v, --verbose         Print more (default: False)
  -cell_group CELL_GROUP_FILE
                        Optional: Provide path to file with group names as
                        headers and columns with cells in that group.
                        (default: False)
  -group_sig_test       True or Flase. When True, if cell groups are provided,
                        perform significance testing between all groups.
                        (default: False)
  -gene_list MAKE_GENE_MATRIX
                        Path to a file with a two columns with headers
                        'GeneID' and 'GroupID' (Singular format). GeneID list
                        will be used to create a new matrix file with only
                        those genes included. (default: False)
  -cell_list MAKE_CELL_MATRIX
                        Path to a file with a column 'SampleID' (Singular
                        format). SampleID list will be used to create a new
                        matrix file with only those cells included. (default:
                        False)
  -stability TEST_CLUST_STABILITY
                        Provide a number of iterations to test how stable
                        clustering is as the number of top PCA genes changes
                        from 100-1000. Output will be clustering heatmaps for
                        each iteration a summary of changes as gene number
                        varies. (default: 0)
  --metric METRIC       The distance metric to use. The distance function can
                        be: braycurtis, canberra, chebyshev, cityblock,
                        correlation, cosine, dice, euclidean, hamming,
                        jaccard, kulsinski, mahalanobis, matching, minkowski,
                        rogerstanimoto, russellrao, seuclidean, sokalmichener,
                        sokalsneath, sqeuclidean, yule. (default: euclidean)
  --method METHOD       The linkage method to use. The linkage algorithm can
                        be: single, complete, average, weighted, centroid
                        median or ward. (default: average)
