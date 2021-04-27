#!/usr/bin/env python
# coding: utf-8


import numpy as np
import networkx as nx #Only needed for breaking cycles
import os


class DotFileProcessor:
    def __init__(self, filename, input_dir=''):
        self.dotfilepath = indir + os.sep + filename
        self._read_file()

    def _read_file(self):
        with open(self.dotfilepath, 'r') as f:
            self.file_content = np.array([line.split('\t\n') for line in f.readlines()]).squeeze()
        # file_content is a numpy array of shape (number of lines, )
        # file_content[0] and file_content[1] contains metadata not important for us 
        self.number_of_lines = len(self.file_content)

    def get_gene_names(self):
        # Create gene_names list which will work as a dictionary for index to gene name mapping
        gene_names = []
        for i in range(2, self.number_of_lines):
            if 'value' in self.file_content[i]:
                break
            gene_names.append(self.file_content[i].split('"', 2)[1])
        return gene_names

    def get_gene_regulation_graph_matrices(self, gene_names):
        number_of_genes = len(gene_names)
        # Create gene regulation dictionary 
        gene_adj_mat = np.zeros((number_of_genes, number_of_genes))
        gene_sign_mat = np.zeros((number_of_genes, number_of_genes))

        for i in range(self.number_of_lines):
            if 'value' in self.file_content[i]:
                segments = self.file_content[i].split('->')
                regulator_gene = segments[0].split('"', 2)[1]
                regulated_gene = segments[1].split('=', 2)[0].split('"', 2)[1]
                sign = segments[1].split('=', 2)[1][1]
                regulator_gene_id = gene_names.index(regulator_gene)
                regulated_gene_id = gene_names.index(regulated_gene)
                if sign == '-':
                    gene_adj_mat[regulator_gene_id][regulated_gene_id] = 1
                    gene_sign_mat[regulator_gene_id][regulated_gene_id] = -1
                elif sign == '+':
                    gene_adj_mat[regulator_gene_id][regulated_gene_id] = 1
                    gene_sign_mat[regulator_gene_id][regulated_gene_id] = 1

        return gene_adj_mat, gene_sign_mat


# In[55]:


class Dot2Sergio:
    """
    Input: The filename, input and output directory as strings
           no_master_autoregular_in_GRN and no_cycle as boolean denoting whether to include master regulators in the GRN file and whether to remove cycle
    Output: Can generate GRN and MR input files
    """
    def __init__(self, filename, input_dir='', output_dir='', no_master_autoregular_in_GRN=True, no_cycle=True):
        self.filename = filename
        self.input_dir = input_dir
        self.output_dir = output_dir
        DFP = DotFileProcessor(filename, input_dir)
        self.gene_names = DFP.get_gene_names()
        self.number_of_genes = len(self.gene_names)
        self.gene_adj_mat, self.gene_sign_mat = DFP.get_gene_regulation_graph_matrices(self.gene_names)
        if no_cycle is True:
            self._break_cycles()
        self.no_master_autoregular_in_GRN = no_master_autoregular_in_GRN
        self.set_number_of_regulators_per_gene()
        self.set_master_and_ar_genes()

    def set_number_of_regulators_per_gene(self):
        self.number_of_regulators_per_gene = self.gene_adj_mat.sum(
            axis=0)  # Column wise sum, expresses how many regulatory genes (in edges) each genes have, those having 0 are master regulators

    def set_master_and_ar_genes(self):
        self.no_ar_master_gene_ids = []
        for i in range(self.number_of_genes):
            if self.number_of_regulators_per_gene[i] == 0:
                self.no_ar_master_gene_ids.append(i)

        self.autoregulator_master_gene_ids = []
        self.autoregulator_target_gene_ids = []
        for i in range(self.number_of_genes):
            if self.gene_adj_mat[i, i] == 1:
                if self.number_of_regulators_per_gene[i] == 1:
                    self.autoregulator_master_gene_ids.append(i)
                if self.number_of_regulators_per_gene[i] > 1:
                    self.autoregulator_target_gene_ids.append(i)

        self.master_gene_ids = self.no_ar_master_gene_ids + self.autoregulator_master_gene_ids

    def _break_cycles(self):  # From each cycle present, remove one edge from that cycle so that the cycle breaks
        # Fist remove the autoregulations and keep track of the autoregulatory nodes
        autoreg_genes = []
        for i in range(self.gene_adj_mat.shape[0]):
            if self.gene_adj_mat[i][i] == 1:
                autoreg_genes.append(i)
                self.gene_adj_mat[i][i] = 0

        # Now break the cycles
        while (True):
            nx_graph = nx.from_numpy_matrix(self.gene_adj_mat, create_using=nx.DiGraph)
            try:
                cycle_edge = nx.find_cycle(nx_graph)[0]
            except:
                break
            self.gene_adj_mat[cycle_edge[0]][cycle_edge[1]] = 0

        # Again store the autoregulatory edges
        for gene in autoreg_genes:
            self.gene_adj_mat[gene][gene] = 1

    def generate_GRN_File(self, hillnumber=2):  # Create input GRN with autoregulation
        # Remove autoregulation for master genes 
        if self.no_master_autoregular_in_GRN is True:
            for ids in self.autoregulator_master_gene_ids:
                self.gene_sign_mat[ids][ids] = 0
            master_gene_ids = self.master_gene_ids
        else:
            master_gene_ids = self.no_ar_master_gene_ids

        assert self.gene_sign_mat.shape[0] == self.number_of_genes
        assert self.gene_sign_mat.shape[1] == self.number_of_genes

        np.random.seed(0)
        self.hillnumber = hillnumber
        line_ls = []

        for target_gene_id in range(self.number_of_genes):
            if target_gene_id not in master_gene_ids:
                regulators = self.number_of_regulators_per_gene[target_gene_id]
                K_queue = []
                line = '{}'.format(float(target_gene_id))
                line += ',' + '{}'.format(float(self.number_of_regulators_per_gene[target_gene_id]))
                for i in range(self.number_of_genes):
                    if self.gene_adj_mat[i][target_gene_id] == 1:
                        line += ',' + '{}'.format(float(i))
                        K_i = self.gene_sign_mat[i][target_gene_id] * np.random.uniform(1, 5)
                        K_queue.append(K_i)
                # write all the K values
                for i in range(len(K_queue)):
                    line += ',' + '{}'.format(float(K_queue[i]))

                # write the hill numbers
                for i in range(len(K_queue)):
                    line += ',' + '{}'.format(float(hillnumber))
                line_ls.append(line)
        # shuffle_the_list
        np.random.shuffle(line_ls)

        output_file_path = self.output_dir + os.sep + self.filename.split('.')[0] + '_input_GRN.txt'
        with open(output_file_path, 'w') as f:
            for line in line_ls:
                f.write(line + '\n')

    def generate_MR_file(self, number_of_cell_types=9, high_min=0.6, high_max=1.0, low_min=0.1, low_max=0.4):
        self.number_of_mrs = len(self.master_gene_ids)
        self.number_of_cell_types = number_of_cell_types
        number_of_expression_states = 2 ** self.number_of_cell_types - 1
        taken_states = []
        assert self.number_of_cell_types <= number_of_expression_states  # Otherwise not possible with only 2 states for each MR
        # Generate States
        while (True):
            if (len(taken_states)) == self.number_of_mrs:
                break
            x = np.random.randint(number_of_expression_states)
            if x not in taken_states:
                taken_states.append(x)

        MR_mat = np.zeros((self.number_of_mrs, self.number_of_cell_types))
        for i, state in enumerate(taken_states):
            state_val = format(state, '0' + str(self.number_of_cell_types) + 'b')
            MR_mat[i, :] = list(state_val)

        expression_range = dict()
        expression_range[1] = (high_min, high_max)
        expression_range[0] = (low_min, low_max)
        for i in range(MR_mat.shape[0]):
            for j in range(MR_mat.shape[1]):
                min_val, max_val = expression_range[MR_mat[i][j]]
                MR_mat[i][j] = np.random.uniform(min_val, max_val)

        line_ls = []
        for i, master_gene in enumerate(self.master_gene_ids):
            line = '{}'.format(master_gene) + ','
            line += ','.join('{}'.format(x) for x in MR_mat[i])
            line_ls.append(line)

        np.random.shuffle(line_ls)

        output_file_path = self.output_dir + os.sep + self.filename.split('.')[0] + '_input_MRs.txt'
        with open(output_file_path, 'w') as f:
            for line in line_ls:
                f.write(line + '\n')


if __name__ == '__main__':
    curdir = os.getcwd()
    indir = curdir + os.sep + 'GNW_sampled_GRNs'
    outdir = curdir + os.sep + 'Input_with_Autoregulation'

    ds = Dot2Sergio(filename='Ecoli_100_net1.dot', input_dir=indir, output_dir=outdir)
    ds.generate_GRN_File(hillnumber=2)
    ds.generate_MR_file(number_of_cell_types=3, high_min=0.6, high_max=1.0, low_min=0.1, low_max=0.4)
