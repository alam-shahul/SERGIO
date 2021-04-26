from SERGIO.sergio import sergio
import numpy as np
import pandas as pd

our_input_targets = 'C:/Users/Alex/Box/MSCBIO2041/Ecoli_1200_net4_input_GRN.txt'
our_input_MRs = 'C:/Users/Alex/Box/MSCBIO2041/Ecoli_1200_net4_input_MRs.txt'

their_input_targets = 'data_sets/De-noised_100G_9T_300cPerT_4_DS1/Interaction_cID_4_testing.txt'
their_input_MRs = 'data_sets/De-noised_100G_9T_300cPerT_4_DS1/Regs_cID_4.txt'

sim = sergio(number_genes=100, number_bins = 9, number_sc = 300, noise_params = 1, decays=0.8, sampling_state=15, noise_type='dpd')
sim.build_graph(input_file_taregts =our_input_targets, input_file_regs=our_input_MRs, shared_coop_state=2)
sim.simulate()
expr = sim.getExpressions()
expr_clean = np.concatenate(expr, axis = 1)