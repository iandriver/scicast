import os
import seaborn as sns
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import difflib

def clust_stability(args, matrix_data, iterations):
    try:
        from .heatmaps import clust_heatmap
        from .dim_reduction import return_top_pca_gene
    except (SystemError, ValueError, ImportError):
        from heatmaps import clust_heatmap
        from dim_reduction import return_top_pca_gene


    stability_ratio = []
    total_genes = len(matrix_data.log2_df_cell.index.tolist())
    end_num = min(1000, total_genes)
    iter_list = list(range(100,int(round(end_num)),int(round(end_num/iterations))))
    for gene_number in iter_list:
        title= str(gene_number)+' genes plot.'
        top_pca_by_gene, top_pca = return_top_pca_gene(args, matrix_data.log2_df_cell, user_num_genes=gene_number)
        top_pca_by_cell = top_pca_by_gene.transpose()
        cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(args, matrix_data, matrix_subset= top_pca_by_gene, title=title, stablity_plot_num=gene_number)
        if gene_number == 100:
            s1 = col_order
            s0 = col_order
        else:
            s2= col_order
            sm_running = difflib.SequenceMatcher(None,s1,s2)
            sm_first = difflib.SequenceMatcher(None,s0,s2)
            stability_ratio.append((sm_running.ratio(), sm_first.ratio()))
            s1=col_order
        plt.close()
    sns.set(context='poster', font_scale = 1)
    sns.set_palette("RdBu_r")
    x= iter_list[1:]
    f, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    y1= [m[0] for m in stability_ratio]
    y2= [m[1] for m in stability_ratio]
    sns.barplot(x, y1, palette="RdBu_r", ax=ax1)
    ax1.set_ylabel('Running ratio (new/last)')
    sns.barplot(x, y2, palette="RdBu_r", ax=ax2)
    ax2.set_ylabel('Ratio to 100')
    plt.savefig(os.path.join(matrix_data.new_filepath,'clustering_stability.'+args.image_format), bbox_inches='tight')
    plt.close('all')
    return stability_ratio
