from __future__ import division, print_function
import os
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import islice
sns.set_style('ticks')


class RubinovAlgorithm:
    def __init__(self, hic_data, chr_number, n_iterations=1000):
        self.hic_file = hic_data
        self.chromosome_number = chr_number
        self.n_iterations = n_iterations
        # self.start = start
        # self.end = end
        self.get_hic_matrix()
        self.get_tad_orientation()
        self.get_first_level_tads()
        self.get_tads()


    def get_hic_matrix(self):
        print('Loading data...')
        hic_file = h5py.File(self.hic_file, 'r')
        chromosome = str(self.chromosome_number) + ' ' + str(self.chromosome_number)
        hic_matrix = hic_file[chromosome].value
        df = pd.DataFrame(hic_matrix)
        df = df.loc[:, (df != 0).any(axis=0)]
        df = df.loc[(df != 0).any(axis=1), :]
        hic_matrix = np.array(df)
        np.fill_diagonal(hic_matrix[1:, :], np.zeros(len(hic_matrix) - 1))
        np.fill_diagonal(hic_matrix[:, 1:], np.zeros(len(hic_matrix) - 1))
        self.hic_for_vizualization = hic_matrix
        self.hic_matrix = np.power(np.array(hic_matrix), 4)
        return self.hic_matrix

    def get_tad_orientation(self):
        left_summs = np.sum(np.triu(self.hic_matrix), axis=0)
        right_summs = np.sum(np.triu(self.hic_matrix), axis=1)
        self.pairs = list(zip(left_summs, right_summs))

        orientation = []
        for pair in self.pairs:
            if pair[0] < pair[1]:
                orientation.append(1)
            elif pair[0] > pair[1]:
                orientation.append(-1)
        self.orientation = orientation

    def get_first_level_tads(self):
        tads_all = []
        tads_left_pos = []
        tads_length = []
        borders = []
        new_pairs = []

        first_tads = []
        first_tads_left_pos = []
        first_tads_length = []
        brdr = []
        resupd = []
        numbers = iter(range(len(self.orientation)-1))
        for i in numbers:
            if self.orientation[i] == 1 and self.orientation[i + 1] == -1:
                brdr.append([i, i + 2])
                if self.pairs[i][0] + self.pairs[i + 1][0] > self.pairs[i][1] + self.pairs[i + 1][1]:
                    first_tads.append(-1)
                elif self.pairs[i][0] + self.pairs[i + 1][0] < self.pairs[i][1] + self.pairs[i + 1][1]:
                    first_tads.append(1)

                nres = [self.pairs[i][0] + self.pairs[i + 1][0], self.pairs[i][1] + self.pairs[i + 1][1]]
                resupd.append(nres)
                first_tads_left_pos.append(i)
                first_tads_length.append(2)
                next(islice(numbers, 1, 1), None)

            else:
                first_tads.append(self.orientation[i])
                first_tads_left_pos.append(i)
                first_tads_length.append(1)
                resupd.append(self.pairs[i])

        new_pairs.append(resupd)
        self.new_pairs = new_pairs

        tads_all.append(self.orientation)
        tads_all.append(first_tads)
        self.tads_all = tads_all

        borders.append(brdr)
        self.borders = borders

        tads_left_pos.append(first_tads_left_pos)
        self.tads_left_pos = tads_left_pos

        tads_length.append(first_tads_length)
        self.tads_length = tads_length

    def get_tads(self):
        print('Getting TADs...')
        for level in range(self.n_iterations):
            tads = []
            tads_left = []
            tads_len = []
            border = []
            updated_pairs = []
            numbers = iter(range(len(self.tads_all[level + 1])-1))

            if len(self.borders[level]) == 1 and self.borders[level][0][0] == 0 and\
                self.borders[level][0][1] == len(self.hic_matrix):
                max_level = level
                print('Highest level', max_level)
                break

            for i in numbers:
                if self.tads_all[level + 1][i] == 1 and self.tads_all[level + 1][i + 1] == -1:
                    border.append([self.tads_left_pos[level][i], self.tads_left_pos[level][i] +
                                   self.tads_length[level][i] + self.tads_length[level][i + 1]])

                    if self.new_pairs[level][i][0] + self.new_pairs[level][i + 1][0] > \
                            self.new_pairs[level][i][1] + self.new_pairs[level][i + 1][1]:
                        tads.append(-1)
                    elif self.new_pairs[level][i][0] + self.new_pairs[level][i + 1][0] < \
                            self.new_pairs[level][i][1] + self.new_pairs[level][i + 1][1]:
                        tads.append(1)

                    n_pairs = [self.new_pairs[level][i][0] + self.new_pairs[level][i + 1][0],
                            self.new_pairs[level][i][1] + self.new_pairs[level][i + 1][1]]

                    updated_pairs.append(n_pairs)
                    tads_left.append(self.tads_left_pos[level][i])
                    tads_len.append(self.tads_length[level][i] + self.tads_length[level][i + 1])
                    next(islice(numbers, 1, 1), None)

                else:
                    tads.append(self.tads_all[level + 1][i])
                    tads_left.append(self.tads_left_pos[level][i])
                    tads_len.append(self.tads_length[level][i])
                    updated_pairs.append(self.new_pairs[level][i])

            self.new_pairs.append(updated_pairs)
            self.tads_all.append(tads)
            self.borders.append(border)
            self.tads_left_pos.append(tads_left)
            self.tads_length.append(tads_len)
        print('Done')

    def mtxplot(self, letter, color):
        for j in range(len(letter)):
            bgn = letter[j][0] - self.startmtx
            end = letter[j][1] - self.startmtx
            plt.plot([bgn + 2, end + 2], [bgn, bgn], color=color)
            plt.plot([end + 2, end + 2], [bgn, end], color=color)

    def visualise_tree(self, start=None, end=None):
        print('Building hierarchical tree...')
        self.startmtx = start + 1
        self.endmtx = end
        plt.figure(figsize=(20, 20))
        sns.heatmap(self.hic_for_vizualization[start:end, start:end], cmap='Reds')
        for i in range(self.startmtx, self.endmtx):
            val = self.orientation[i]
            plt.text(i - self.startmtx, i + 0.8 - self.startmtx, val, {'color': 'red' if val > 0 else 'blue'})
        if len(self.borders) % 3 == 0:
            colors = ['red', 'blue', 'green'] * int(len(self.borders) / 3)
        else:
            colors = ['red', 'blue', 'green'] * int(len(self.borders) / 3)
            colors.append(['orange', 'blue', 'green'][: len(self.borders) % 3])
        for letter, color in zip(self.borders, colors):
            self.mtxplot(letter, color)
        plt.savefig(str(os.path.splitext(self.hic_file)[0]) + '_hierarhical_tree' + '.pdf', format='pdf')
        print('Done')

