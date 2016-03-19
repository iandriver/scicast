from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP
import rpy2.robjects as robjects
from rpy2 import rinterface as ri
import os


norm_func_path = '/Users/idriver/RockLab-files/SCICAST/SCICAST/normalize_functions.R'
r_source = robjects.r['source']
r_source(norm_func_path)
r_make_cpm = robjects.globalenv['make_cpm']
r_run_deseq2 = robjects.globalenv['run_deseq2']

#path to fpkm file (usually cuffnorm output)
path_to_file = '/Volumes/Seq_data/count-picard_Pdgfra_ctrl'
#default file name will use genes.fpkm_table from cuffnorm
file_name = 'Pdgfra_ctrl_all_count_table.txt'
#provide base name for output files
base_name = robjects.r('"Pdgfra_ctrl_counts"')
file_path = robjects.r('"/Volumes/Seq_data/count-picard_Pdgfra_ctrl/Pdgfra_ctrl_all_count_table.txt"')
robjects.globalenv['file_path'] = file_path
robjects.globalenv['base_name'] = base_name
print(file_path)
group1_terms = ['ctrl1']
group1_name = ['all']
group2_terms = ['low']
group2_name = ['low']
r_make_cpm(file_path , base_name)
dds = r_run_deseq2(file_path, base_name, robjects.vectors.StrVector(group1_terms), robjects.vectors.StrVector(group1_name), robjects.vectors.StrVector(group2_terms), robjects.vectors.StrVector(group2_name))
