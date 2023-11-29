#!/usr/bin/env python
# coding=utf-8
from __future__ import print_function

from builtins import range
from rdkit import Chem
from itertools import combinations
import networkx as nx
from base import BaseMixin
from base_fragment import BaseFragment
from fun_group_table import bigring_bond_score

class Compound(BaseMixin):
    def __init__(self, rdmol):
        '''Class for molecule fragmentation
        '''
        self.rdmol = rdmol
        self.basefrag = None
        self.broken_bonds = []
        self.bigRing_broken_bonds = []

    @classmethod
    def comFromMolFile(cls, filename):
        '''
        Args:
            filename (str): path to the mol file

        Returns:
            :obj:`Compound`
        '''
        rdmol = Chem.MolFromMolFile(filename, removeHs=False)

        return cls(rdmol)

    def cutCompound(self):
        '''Cut the compound into elementary fragments
           Return all elementary fragments of the compound'''
        rings, big_rings = self.findRings()
        rdmol = self.rdmol
        if big_rings:
            rdmol, uncutted_rings = self.cutBigRing(big_rings)
            for ring in uncutted_rings:
                rings.append(list(ring))
        basefrag = BaseFragment.BaseFragmentFromRdmol(rdmol)
        basefrag.cutAllChain(rings)
        self.basefrag = basefrag

        return self.basefrag

    def cutBigRing(self, rings):
        '''Cut big rings'''
        broken_bonds_idx = []
        G = self.getGraph(removeHs=True)
        degree = G.degree()
        big_rings = {tuple(ring) for ring in rings if len(ring) > 10}
        cut_shared = False

        while big_rings:
            shared_edges = {}
            for a, b in combinations(big_rings, 2):
                share = set(a) & set(b)
                for i, j in combinations(share, 2):
                    if G.has_edge(i, j):
                        edge = (i, j) if i < j else (j, i)
                        if edge in shared_edges:
                            shared_edges[edge].add(a)
                            shared_edges[edge].add(b)
                        else:
                            shared_edges[edge] = set((a, b))
            if cut_shared:
                if not shared_edges:
                    break
                min_shared_edges = min(len(v) for v in shared_edges.values())
            broken_rings = set()
            for br in big_rings:
                bond_score_dict = {}
                for edge in G.edges():
                    if G.get_edge_data(*edge)['bond_type'] != 'SINGLE':
                        continue
                    elif (not cut_shared and edge in shared_edges) or \
                            (cut_shared and len(shared_edges.get(edge, [])) > min_shared_edges):
                        continue
                    elif edge[0] in br and edge[1] in br:
                        s_idx, e_idx = edge[0], edge[1]
                        s_atom = G.node[s_idx]['atom'] if G.node[s_idx]['atom'] == 'C' else 'Y'
                        e_atom = G.node[e_idx]['atom'] if G.node[e_idx]['atom'] == 'C' else 'Y'
                        bond_score = bigring_bond_score[(
                            s_atom, e_atom, degree[s_idx], degree[e_idx])]
                        if bond_score == 0:
                            broken_bonds_idx.append(edge)
                            broken_rings.add(br)
                            bond_score_dict = {}
                            if cut_shared:
                                for bri in shared_edges[edge]:
                                    broken_rings.add(bri)
                            break
                        else:
                            bond_score_dict[bond_score] = edge
                if bond_score_dict:
                    scores = sorted(bond_score_dict.keys())
                    broken_bonds_idx.append(bond_score_dict[scores[0]])
                    broken_rings.add(br)
                    if cut_shared:
                        cut_shared = False
                        for bri in shared_edges[broken_bonds_idx[-1]]:
                            broken_rings.add(bri)
                        break
                if cut_shared and br in broken_rings:
                    cut_shared = False
                    break
            if broken_rings:
                for br in broken_rings:
                    big_rings.remove(br)
            else:
                cut_shared = True

        if not broken_bonds_idx:
            return self.rdmol, big_rings
        self.bigRing_broken_bonds = [
            self.rdmol.GetBondBetweenAtoms(s, e) for s, e in broken_bonds_idx]
        rw = Chem.RWMol(self.rdmol)
        for s, e in broken_bonds_idx:
            rw.RemoveBond(s, e)
        new_rdmol = rw.GetMol()

        return new_rdmol, big_rings

    def findRings(self):
        '''Find all rings in the molecule and classify them into rings and big rings (>= 10 member)'''
        big_rings, rings = [], []
        smallest_rings = list(self.findSmallestRings())
        if len(smallest_rings) == 0:
            return [], []
        elif len(smallest_rings) == 1:
            if len(smallest_rings[0]) > 10:
                return [], [list(smallest_rings[0])]
            else:
                return [list(smallest_rings[0])], []
        set_big = {i for i, ringa in enumerate(smallest_rings)
                   if len(ringa) > 10}
        G = nx.Graph()
        G.add_nodes_from(set(range(len(smallest_rings))) - set_big)
        if not hasattr(G, 'node'):
            G.node = G.nodes
        for n, m in combinations(G.node, 2):
            set_n, set_m = set(smallest_rings[n]), set(smallest_rings[m])
            if set_n & set_m:
                G.add_edge(n, m)
        group = nx.connected_components(G)
        for grp in group:
            ring = set()
            for i in grp:
                for idx in smallest_rings[i]:
                    ring.add(idx)
            rings.append(list(ring))

        # Deal with two rings connected by a double bond
        atom_ring_idx = {}
        for idx, grp in enumerate(rings):
            for i in grp:
                atom_ring_idx[i] = idx
        G = nx.Graph()
        G.add_nodes_from(range(self.rdmol.GetNumAtoms()))
        for bond in self.rdmol.GetBonds():
            begin_idx, end_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE or \
                    (begin_idx in atom_ring_idx and atom_ring_idx.get(end_idx, -1) == atom_ring_idx[begin_idx]):
                G.add_edge(begin_idx, end_idx)
        if G.edges():
            tmp_rings = []
            group = nx.connected_components(G)
            for grp in group:
                ring = []
                has_ring = False
                for idx in grp:
                    ring.append(idx)
                    if idx in atom_ring_idx:
                        has_ring = True
                if has_ring:
                    tmp_rings.append(ring)
            rings = tmp_rings

        ring_atoms = {i for ring in rings for i in ring}
        for br in set_big:
            if any((i in ring_atoms for i in smallest_rings[br])):
                continue
            big_rings.append(list(smallest_rings[br]))
        return rings, big_rings


if __name__ == "__main__":
    '''Input: MOL/SDF file of the input molecule
    Output: BaseFragment object, which contains the following two attributes
    BaseFragment.frag_rings (list of list of int): a list of ring elementary fragments cut from the input molecule. 
                                                   Each inner list records the indices from the original molecule that constitue the fragment.
    BaseFragment.frag_chains (list of list of int): a list of chain elementary fragments cut from the input molecule.
                                                   Each inner list records the indices from the original molecule that constitue the fragment.
    '''
    cmpd = Compound.comFromMolFile("ejm-31.sdf")
    frags = cmpd.cutCompound()
    print(frags.frag_rings)
    print(frags.frag_chains)