scicast (Single Cell Iterative Clustering and Significance Testing) what does it do?
------------
When approaching large RNA sequencing datasets (single cell or bulk with replicates) the first question is often what are the subsets/groups/populations and what are the defines those groups? This software is designed to do that rapidly and provide the results in multiple visual and data formats that can be used for further analysis.

How do I use it?
------------
After following the installation instructions, it is usable through a gui and through the command line. Both the gui and the command line do the same thing and the command is logged in "scicast_command_log.txt" so that it can be re-run by scicast "copy from log".

**scicast vignette**
==================
Gui and command line for each step will be shown.

Using the provided "krasnow_AT2_185_rsem_deseq_counts_norm.txt" file.

![GUI: Load the dataset](scicast_with_parameters1.png)

Command Line:

```bash
        scicast -f path_to_file/krasnow_AT2_185_rsem_deseq_counts_norm.txt -n
        vignette -method ward -metric seuclidean -g 200 -depth 100 -z 0
        -qgraph_plot none -kmeans_cluster_range 2,4 -kmeans_sig_test
```

![Files output:](scicast_filelist_ouput1.png)

1.  "scicast_command_log.txt" - log of command(s)  
  Records exact command used to generate the output in the folder. Can be copy pasted to re-run and serves as a record of how results were generated.

  Example:

  ```bash
          2017-01-07 14:11
          -f /Users/iandriver/scicast/krasnow_AT2_185_rsem_deseq_counts_norm.txt
          -n vignette -method ward -metric seuclidean -g 200 -depth 100 -z 0
          -qgraph_plot none -kmeans_cluster_range 2,4 -kmeans_sig_test
          -image_format png
  ```  
2.  "top_log2_genes_plot.png" - plot of expression in log2 of top genes for each sample.  
  An overview of expression in all samples and a way to see what samples are included after filtering.

  Example:

  ![top_log2](vignette_scicast_analysis/top_log2_genes_plot.png)  
3.  "distribution_log2_genes_plot.png" - distribution of mean expression before and after filtering.  
  The distribution of mean log2 expression of all genes before and after filtering (often the same, especially if filtered prior to scicast).

  Example:

  ![log2_dist](vignette_scicast_analysis/distribution_log2_genes_plot.png)  
