from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP
import rpy2.robjects as robjects
import os


with open('/Users/idriver/RockLab-files/SCICAST/normalize_functions.R', 'r') as f:
    r_file = f.read()
norm_functions = STAP(r_file, "normalize_functions")

#path to fpkm file (usually cuffnorm output)
path_to_file = '/Volumes/Seq_data/count-picard_Pdgfra_ctrl'
#default file name will use genes.fpkm_table from cuffnorm
file_name = 'Pdgfra_ctrl_all_count_table.txt'
#provide base name for output files
base_name ='Pdgfra_ctrl_counts'
file_path = path_to_file+'/'+file_name
group1_terms = 'ctrl1'
group1_name = 'all'
group2_terms = 'low'
group2_name = 'low'
norm_functions.make_cpm(robjects.r(file_path) , robjects.r(base_name))
norm_functions.run_deseq2(robjects.r(file_path), robjects.r(base_name), robjects.r(group1_terms), robjects.r(group1_name), robjects.r(group2_terms), robjects.r(group2_name))
