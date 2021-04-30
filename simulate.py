from SERGIO.sergio import sergio
import numpy as np
import os
import pandas as pd

our_input_targets = 'C:/Users/Alex/Box/MSCBIO2041/Ecoli_1200_net4_input_GRN.txt'
our_input_MRs = 'C:/Users/Alex/Box/MSCBIO2041/Ecoli_1200_net4_input_MRs.txt'

their_input_targets = 'data_sets/De-noised_100G_9T_300cPerT_4_DS1/Interaction_cID_4_testing.txt'
their_input_MRs = 'data_sets/De-noised_100G_9T_300cPerT_4_DS1/Regs_cID_4.txt'

AR_GRN_1 = 'Input_with_Autoregulation/Ecoli_100_net1_input_GRN.txt'
AR_MR_1 = 'Input_with_Autoregulation/Ecoli_100_net1_input_MRs.txt'
AR_GRN_2 = 'Input_with_Autoregulation/Ecoli_1200_net4_input_GRN.txt'
AR_MR_2 = 'Input_with_Autoregulation/Ecoli_1200_net4_input_MRs.txt'

# Simulate data for 100 genes
sim = sergio(number_genes=100, number_bins=3, number_sc=300, noise_params=1, decays=0.8, sampling_state=15, noise_type='dpd')
sim.build_graph(input_file_taregts=AR_GRN_1, input_file_regs=AR_MR_1, shared_coop_state=2)
sim.simulate()
expr = sim.getExpressions()
expr_clean = np.concatenate(expr, axis=1)


# Add technical noise
"""
Add outlier genes
"""
expr_O = sim.outlier_effect(expr, outlier_prob = 0.01, mean = 0.8, scale = 1)

"""
Add Library Size Effect
"""
libFactor, expr_O_L = sim.lib_size_effect(expr_O, mean = 4.6, scale = 0.4)

"""
Add Dropouts
"""
binary_ind = sim.dropout_indicator(expr_O_L, shape = 6.5, percentile = 82)
expr_O_L_D = np.multiply(binary_ind, expr_O_L)

"""
Convert to UMI count
"""
count_matrix = sim.convert_to_UMIcounts(expr_O_L_D)

"""
Make a 2d gene expression matrix
"""
count_matrix = np.concatenate(count_matrix, axis = 1)

# Write out final count matrix
path = os.path.expanduser('Output_with_Autoregulation/ds1_ar_counts.txt')
np.savetxt(path, count_matrix)

# Simulate data for 1200 genes
sim = sergio(number_genes=1200, number_bins=9, number_sc=300, noise_params=1, decays=0.8, sampling_state=15, noise_type='dpd')
sim.build_graph(input_file_taregts=AR_GRN_2, input_file_regs=AR_MR_2, shared_coop_state=2)
sim.simulate()
expr = sim.getExpressions()
expr_clean = np.concatenate(expr, axis=1)

# Add technical noise
"""
Add outlier genes
"""
expr_O = sim.outlier_effect(expr, outlier_prob = 0.01, mean = 0.8, scale = 1)

"""
Add Library Size Effect
"""
libFactor, expr_O_L = sim.lib_size_effect(expr_O, mean = 4.6, scale = 0.4)

"""
Add Dropouts
"""
binary_ind = sim.dropout_indicator(expr_O_L, shape = 6.5, percentile = 82)
expr_O_L_D = np.multiply(binary_ind, expr_O_L)

"""
Convert to UMI count
"""
count_matrix = sim.convert_to_UMIcounts(expr_O_L_D)

"""
Make a 2d gene expression matrix
"""
count_matrix = np.concatenate(count_matrix, axis = 1)

# Write out final count matrix
path = os.path.expanduser('Output_with_Autoregulation/ds4_ar_counts.txt')
np.savetxt(path, count_matrix)