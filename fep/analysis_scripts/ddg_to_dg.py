#!/usr/bin/env python
# coding=utf-8

import networkx as nx
import numpy as np
import math
from collections import defaultdict

class CycleCloseCorrect(object):
    """ Cycle Closure Correction for calculated ddG of all pairs
        The implementation details of the CCC algorithm can not be shown.
    """

    def __init__(self, ddg_dic):

    def get_ccc_ddg(self):

        return ccc_ddg_dic, ccc_ddg_err_dic


class DDG2DG(object):
    """ convert the ddG after cycle closure correction to dG according to the user-specified exp_dG using the MLE method
        The code was adapted from https://github.com/OpenFreeEnergy/cinnabar/blob/main/cinnabar/stats.py
        Reference: https://pubs.acs.org/doi/abs/10.1021/acs.jcim.9b00528
    """

    def __init__(self, ccc_ddg_dic, ccc_ddg_err_dic, exp_dg_dic):
        """ init the DDG2DG

        :param ccc_ddg_dic: ddG dict after cycle closure correction for all pairs
        :param ccc_ddg_err_dic: ddG cycle closure correction error dict for all pairs
        :param exp_dg_dic: user-specified exp_dG dict for certain ligands
        """
        self._exp_dg_dic = exp_dg_dic
        self._ccc_ddg_dic = ccc_ddg_dic
        self._ccc_ddg_err_dic = ccc_ddg_err_dic

    def _construct_data(self):
        lig_set = set()
        for pair in self._ccc_ddg_dic:
            lig1, lig2 = pair.split("~")
            lig_set.add(lig1)
            lig_set.add(lig2)
        self._N = len(lig_set)

        self._ddg_matrix = np.zeros((self._N, self._N))
        self._ddg_err_matrix = np.zeros((self._N, self._N))
        self._F = np.zeros((self._N, self._N))
        self._z = np.zeros([self._N])

        self._lig_list = list(lig_set)
        for pair in self._ccc_ddg_dic:
            lig1, lig2 = pair.split("~")
            i = self._lig_list.index(lig1)
            j = self._lig_list.index(lig2)
            self._ddg_matrix[i, j] = -self._ccc_ddg_dic[pair]
            self._ddg_matrix[j, i] = self._ccc_ddg_dic[pair]
            self._ddg_err_matrix[i, j] = self._ccc_ddg_err_dic[pair]
            self._ddg_err_matrix[j, i] = self._ccc_ddg_err_dic[pair]

        for pair in self._ccc_ddg_dic:
            lig1, lig2 = pair.split("~")
            i = self._lig_list.index(lig1)
            j = self._lig_list.index(lig2)
            if abs(self._ddg_err_matrix[i, j]) > 1e-6:
                self._F[i, j] = -self._ddg_err_matrix[i, j] ** (-2)
                self._F[j, i] = -self._ddg_err_matrix[i, j] ** (-2)
        for i, lig in enumerate(self._lig_list):
            self._F[i, i] = -np.sum(self._F[i, :])

        for pair in self._ccc_ddg_dic:
            lig1, lig2 = pair.split("~")
            i = self._lig_list.index(lig1)
            j = self._lig_list.index(lig2)
            if abs(self._ddg_err_matrix[i, j]) > 1e-6:
                self._z[i] += self._ddg_matrix[i, j] * self._ddg_err_matrix[i, j] ** (-2)
                self._z[j] += self._ddg_matrix[j, i] * self._ddg_err_matrix[j, i] ** (-2)

    def get_predict_dg(self):
        self._construct_data()

        Finv = np.linalg.pinv(self._F)
        f_i = np.matmul(Finv, self._z)

        # Compute uncertainty
        variance = np.diagonal(Finv)
        variance = variance ** 0.5
 
        predict_dg_dic = {}
        predict_dg_err_dic = {}

        for idx, lig in enumerate(self._lig_list):
            predict_dg_dic[lig] = f_i[idx]
            predict_dg_err_dic[lig] = variance[idx]

        # Shift to have lowest RMSE
        exp_data = [self._exp_dg_dic[lig] for lig in self._lig_list]
        shift = calc_shift(exp_data, f_i)
        for lig in predict_dg_dic:
            predict_dg_dic[lig] -= shift

        return predict_dg_dic, predict_dg_err_dic


def ddg_to_dg(ddg_dic, ddg_err_dic, exp_dg_dic):
    """ correct the calculated ddG for all pairs, and convert to dG for all ligands

    :param ddg_dic: calculated ddG dict of all pairs
    :param ddg_err_dic: calculated ddG error dict of all pairs
    :param exp_dg_dic: user-specified exp_dG dict for certain ligands
    :return: corrected ddG for all pairs and converted dG for all ligands and their associated errors
    """
    # Cycle closure correction
    cycle_close_correct = CycleCloseCorrect(ddg_dic)
    ccc_ddg_dic, ccc_ddg_err_dic = cycle_close_correct.get_ccc_ddg()
    # For those edges not involved in a cycle, using the original fep error
    for pair in ccc_ddg_err_dic:
        if abs(ccc_ddg_err_dic[pair]) < 1e-6:
            ccc_ddg_err_dic[pair] = ddg_err_dic[pair]

    # CCC ddG to dG
    ddg2dg = DDG2DG(ccc_ddg_dic, ccc_ddg_err_dic, exp_dg_dic)
    predict_dg_dic, dg_err_dic = ddg2dg.get_predict_dg()

    return ccc_ddg_dic, predict_dg_dic, ccc_ddg_err_dic, dg_err_dic

def calc_shift(exp, sim):
    """ shift two sets of data to have the lowest RMSE
    """
    sim_aver=np.average(sim)
    exp_aver=np.average(exp)
    bsim = np.subtract(sim, sim_aver)
    bexp = np.subtract(exp, exp_aver)
    c_temp = np.sum(bsim - bexp) / (float(len(bsim)))
    c = sim_aver - exp_aver + c_temp
    
    return c

if __name__ == "__main__":
    # Test data
    import pandas as pd

    data = pd.read_csv("pfkfb3_Merck_ddG_test.csv")
    dg_exp_dict = {}
    fep_ddg_dict = {}
    fep_err_dict = {}
    for idx, row in data.iterrows():
        lig1 = str(int(row["Ligand1"]))
        lig2 = str(int(row["Ligand2"]))
        lig_pair = lig1+'~'+lig2
        fep_ddg_dict[lig_pair] = row["FEP"]
        fep_err_dict[lig_pair] = row["FEP Error"]

    data2 = pd.read_csv("pfkfb3_Merck_dG_test.csv")
    for idx, row in data2.iterrows():
        lig = str(int(row["Ligand"]))
        dg_exp_dict[lig] = row["Pred. dG"]

    ccc_ddg_dic, predict_dg_dic, ccc_ddg_err_dic, dg_err_dic = ddg_to_dg(fep_ddg_dict, fep_err_dict, dg_exp_dict)
    print(ccc_ddg_dic, predict_dg_dic, ccc_ddg_err_dic, dg_err_dic)

