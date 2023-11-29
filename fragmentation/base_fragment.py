#!/usr/bin/env python
# coding=utf-8
from __future__ import print_function, unicode_literals

from builtins import range
from builtins import object
import networkx as nx
from itertools import combinations
from rdkit import Chem

CHAIN_NUM = 0

class BaseAtom(object):
    '''base atom class

    Args:
        symbol (str): element symbol
        idx (int): atom index
        score (int, optional): score defined in fun_group_table, default 0
        chain_num (int, optional): chain number the atom belongs to, default -1
    '''
    def __init__(self, symbol, idx, score=0, chain_num=-1):
        self.symbol = symbol
        self.idx = idx
        self.score = score
        self.chain_num = chain_num

class BaseBond(object):
    '''base bond class

    Args:
        beginIdx, endIdx (int, int): indices of two end atoms
        bondtype (int or :obj:`rdkit.Chem.rdchem.BondType`, optional): bond type, default -1
    '''
    def __init__(self, beginIdx, endIdx, bondtype=1):
        self.beginIdx = beginIdx
        self.endIdx = endIdx
        self.bondtype = bondtype

class BaseFragment(object):
    '''Base fragment class, provide methods for molecule fragmentation

    Args:
        baseAtoms (list of :obj:`BaseAtom`): atom list
        baseBonds (list of :obj:`BaseBond`): bond list

    Attributes:
        baseAtoms (list of :obj:`BaseAtom`): atom list
        baseBonds (list of :obj:`BaseBond`): bond list
        frag_rings (list of list of int): ring fragments list
        frag_chains (list of list of int): chain fragments list
    '''

    def __init__(self, baseAtoms, baseBonds):
        self.baseAtoms = baseAtoms
        self.baseBonds = baseBonds
        self.frag_rings = []
        self.frag_chains = []

    def createGraph(self, removeHs=True):
        '''create graph for the molecule

        Args:
            removeHs (bool, optional): If remove H, default True

        Returns:
            :obj:`networkx.classes.graph.Graph`
        '''
        G = nx.Graph()
        for atom in self.baseAtoms:
            G.add_node(atom.idx, atom=atom)
        for bond in self.baseBonds:
            G.add_edge(bond.beginIdx, bond.endIdx, bondtype=bond.bondtype)
        if removeHs:
            for atom in self.baseAtoms:
                if atom.symbol == 'H':
                    G.remove_node(atom.idx)
        if not hasattr(G, 'node'):
            G.node = G.nodes
        return G

    def getAllAtomDegree(self):
        '''return a dict of heavy atoms and their degree'''
        G = self.createGraph(removeHs=False)
        degree = dict(G.degree())
        for node in G.nodes():
            if G.node[node]['atom'].symbol == 'H':
                degree.pop(node)
        return degree

    @classmethod
    def BaseFragmentFromRdmol(cls, rdmol):
        '''Creat a BaseFragment class from rdmol

        Args:
            rdmol (:obj:`rdkit.Chem.rdchem.Mol`): RDKit molecule

        Returns:
            :obj:`BaseFragment`
        '''
        baseAtoms = [BaseAtom(atom.GetSymbol(), atom.GetIdx()) for atom in rdmol.GetAtoms()]
        baseBonds = [BaseBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bondtype=bond.GetBondType()) for bond in rdmol.GetBonds()]

        return cls(baseAtoms, baseBonds)

    def removeGroup(self, groups):
        '''Remove groups of atoms from molecule graph. Also remove halogen atoms.

        Args:
            groups (list of list of int): list of atom groups to remove

        Returns:
            :obj:`networkx.classes.graph.Graph`: new networkx graph after removing atom groups
        '''
        G = self.createGraph()
        degree = G.degree()
        for grp in groups:
            for idx in grp:
                G.remove_node(idx)
        for node in list(G.nodes()):
            if G.node[node]['atom'].symbol in halogen and degree[node] == 1:
                G.remove_node(node)
        return G

    def getRingAndChainFragment(self, rings):
        '''Cut the molecule into ring fragments and chain fragments that wait for further cutting

        A ring fragment contains all atoms on the ring, all halogen atoms attached to the ring, and all atoms double bonded to the ring.

        Args:
            rings (list of list of int): list of rings fragments

        Returns:
            list of list of int: list of chain fragments waiting for further cutting
        '''
        group_chain = []
        self.mergeAtomToRing(rings)
        self.checkRings(rings)
        G = self.removeGroup(rings)
        group_list = [list(grp) for grp in nx.connected_components(G)]
        for grp in group_list:
            if len(grp) < 4:
                self.frag_chains.append(grp)
            else:
                group_chain.append(grp)
        self.frag_rings = rings

        return group_chain

    def checkRings(self, rings):
        '''Check whether multiple rings contain repeat atoms.'''
        total_count = sum([len(ring) for ring in rings])
        tmp_ring = []
        for ring in rings:
            tmp_ring += ring
        set_count = len(set(tmp_ring))
        assert total_count == set_count, 'duplicated ring atoms'

    def mergeAtomToRing(self, rings):
        '''Merge all halogen atoms that attached to a ring and atoms double bonded to a ring into the ring.

        Args:
            rings (list of list of int): list of ring fragments
        '''
        if not rings:
            return
        assert all((isinstance(ring, list) for ring in rings))
        G = self.createGraph()
        for ring in rings:
            tmp = []
            for r in ring:
                for nei in G.neighbors(r):
                    if nei not in ring:
                        if G.get_edge_data(r, nei)['bondtype'] != Chem.rdchem.BondType.SINGLE or (
                            G.degree(nei) == 1 and G.node[nei]['atom'].symbol in halogen):
                            tmp.append(nei)
            ring += tmp

    def getAllChainScore(self, group_chain):
        '''calculate the score for each functional group in a chain fragment

        Args:
            group_chain (list of list of int): list of chain fragments that need to cut

        Returns:
            list of :obj:`networkx.classes.graph.Graph`: the graph for each chain with calculated scores
        '''
        G = self.createGraph()
        G = self.nonC2Y(G)
        chain_graphs = []
        for grp in group_chain:
            subG = G.subgraph(grp)
            self.chainScore(subG)
            chain_graphs.append(subG)

        return chain_graphs

    def nonC2Y(self, G):
        '''convert non carbon element in the graph to Y'''
        for node in G.node.values():
            if node['atom'].symbol != 'C':
                atom = BaseAtom('Y', node['atom'].idx)
            else:
                atom = BaseAtom('C', node['atom'].idx)
            node['atom'] = atom
        return G

    def getSubFromGraph(self, G, cent, max_deep=2):
        '''traverse all vertexes and their pathes with depth less than max_deep, return the corresponding subgraph

        Args:
            G (:obj:`networkx.classes.graph.Graph`): NetworkX graph
            cent (int): indices of the vertexes
            max_deep (int): search depth, default 2

        Returns:
            :obj:`networkx.classes.graph.Graph`: subgraph of G
        '''
        assert max_deep > 0, "max_deep must > 0"
        nodes_times = dict(G.degree())
        per = set(G.neighbors(cent))
        visit = list(per) + [cent]
        for i in visit:
            nodes_times[i] = -1
        cur_deep = 1
        while cur_deep < max_deep:
            tmp = []
            while per:
                neig = list(G.neighbors(per.pop()))
                tmp += neig
            for t in tmp:
                if nodes_times[t] >= 0:
                    nodes_times[t] = -1
                    per.add(t)
                    visit.append(t)
            cur_deep += 1

        return G.subgraph(visit)

    def isMatchBondType(self, basefrag):
        '''return whether the fragment contains non-single bond

        Args:
            basefrag (:obj:`BaseFragment`): fragment to check

        Returns:
            bool: whether the fragment contains non-single bond
        '''
        for bond in basefrag.baseBonds:
            if bond.bondtype > 1:
                return True
        return False

    def chainScore(self, G, debug=False):
        '''calculate score for the functional groups in the chain and save them in the score attribute

        Args:
            G (:obj:`networkx.classes.graph.Graph`): the graph for the chain that need to calculate score
        '''
        import fun_group_table
        fun_group_dict = fun_group_table.fun_grop_dict
        if not hasattr(G, 'node'):
            G.node = G.nodes
        degree = dict(G.degree())
        for key, value in degree.items():
            if value == 1:
                for bf in fun_group_dict[value]:
                    if bf.baseAtoms[0].symbol == G.node[key]['atom'].symbol:
                        G.node[key]['atom'].score = bf.baseAtoms[0].score
            elif value > 1 and value < 5:
                for bf in fun_group_dict[value]:
                    sub_g = self.getSubFromGraph(G, key, max_deep=2)
                    is_match_bondtype = self.isMatchBondType(bf)
                    if is_match_bondtype:
                        g = bf.createGraph()
                        gm = nx.algorithms.isomorphism.GraphMatcher(sub_g, g, node_match=lambda n1, n2 : n1['atom'].symbol == n2['atom'].symbol,
                                                                    edge_match=lambda b1, b2 : b1['bondtype'] == b2['bondtype'])
                        if self.choiceMatch(gm, key):
                            G.node[key]['atom'].score = bf.baseAtoms[0].score
                            break
                    else:
                        g = bf.createGraph()
                        gm = nx.algorithms.isomorphism.GraphMatcher(sub_g, g, node_match=lambda n1, n2 : n1['atom'].symbol == n2['atom'].symbol)
                        if self.choiceMatch(gm, key):
                            G.node[key]['atom'].score = bf.baseAtoms[0].score
                            break
            else:
                G.node[key]['atom'].score = 36
        return G

    def choiceMatch(self, gm, idx):
        '''judge if the functional group graph match is success

        Args:
            gm (:obj:`networkx.algorithms.isomorphism.GraphMatcher`): the graph matcher object
            idx (int): the index of the atom to be scored

        Returns:
            match failure -- False;
            match success -- dict (int: int): return match result
        '''
        for match in gm.subgraph_isomorphisms_iter():
            if idx not in match.keys():
                continue
            if match[idx] == 0:
                return match
        return False

    def getSingleChain(self, chain_graph):
        '''Cut a given chain graph into several single chains

        After the chain is cut, the chain number of single chains are saved in the chain_num attribute of each atom; 
        Merge single atom chain into the main chain;
        Chains with atoms ranging from 2 to 4 are considered single chains. Longer chains will continue to be cut.

        Args:
            chain_graph (:obj:`networkx.classes.graph.Graph`): the graph for each chain
        '''
        global CHAIN_NUM
        if not hasattr(chain_graph, 'node'):
            chain_graph.node = chain_graph.nodes
        degree = dict(chain_graph.degree())
        if not self.checkDegree(chain_graph, degree):
            CHAIN_NUM += 1
            for idx in chain_graph.nodes():
                chain_graph.node[idx]['atom'].chain_num = CHAIN_NUM
            return
        terminal_atom = [key for key, value in degree.items() if value == 1]
        mainchain = self.findMainChain(chain_graph, terminal_atom)
        g_copy = chain_graph.copy()
        for idx in mainchain:
            g_copy.remove_node(idx)
        group_list = [list(grp) for grp in nx.connected_components(g_copy)]
        single_chains = []
        many_chains = []
        for grp in group_list:
            if len(grp) == 1:
                mainchain.append(grp[0])
            elif len(grp) <=4:
                single_chains.append(grp)
            else:
                many_chains.append(grp)
        CHAIN_NUM += 1
        for idx in mainchain:
            chain_graph.node[idx]['atom'].chain_num = CHAIN_NUM
        for sgl in single_chains:
            CHAIN_NUM += 1
            for idx in sgl:
                chain_graph.node[idx]['atom'].chain_num = CHAIN_NUM
        many_chains_graph = [chain_graph.subgraph(grp) for grp in many_chains]
        for onegrp in many_chains_graph:
            self.getSingleChain(onegrp)

    def getAllSingleChains(self, rings):
        '''Cut the entire molecule into several single chains and record the chain number into the chain_num attribute of each atom. 

        Args:
            rings (list of list of int): list of ring fragments

        Returns:
            list of :obj:`networkx.classes.graph.Graph`: list of graphes for single chains
        '''
        group_chain = self.getRingAndChainFragment(rings)
        chain_graphs = self.getAllChainScore(group_chain)
        all_single_chains = []
        for g in chain_graphs:
            self.getSingleChain(g)
        for g in chain_graphs:
            chain_type_set = set([g.node[idx]['atom'].chain_num for idx in g.nodes()])
            singchain_idx_dict = {num:[] for num in chain_type_set}
            for node in g.nodes():
                singchain_idx_dict[g.node[node]['atom'].chain_num].append(node)
            for value in singchain_idx_dict.values():
                chain = g.subgraph(value)
                all_single_chains.append(chain.copy())

        global CHAIN_NUM
        CHAIN_NUM = 0
        return all_single_chains

    def findMainChain(self, G, terminal_atom):
        '''Find the main chain in the fragment

        A main chain is the longest path in the molecule. If there are two pathes with equal
        length, pick the one with the highest sum of atom scores. Merge all atoms that are
        not single bonded to the main chain into the main chain

        Args:
            G (:obj:`networkx.classes.graph.Graph`): networkx graph
            terminal_atom (list of int): list of indices for the terminal atoms

        Returns:
            list of int: list of atoms for the main chain
        '''
        all_paths = []
        maxlen = 0
        if not hasattr(G, 'node'):
            G.node = G.nodes
        path = dict(nx.all_pairs_shortest_path(G))
        for s, e in combinations(terminal_atom, 2):
            if len(path[s][e]) > maxlen:
                maxlen = len(path[s][e])
            all_paths.append(path[s][e])
        same_long = []
        for cha in all_paths:
            if len(cha) == maxlen:
                same_long.append(cha)
        for sam in same_long:
            total_score = 0
            tmp_score = sum([G.node[idx]['atom'].score for idx in sam])
            if tmp_score > total_score:
                mainchain = sam
                total_score = tmp_score
        extra_main = []
        for main_idx in mainchain:
            for end_idx in G.neighbors(main_idx):
                if end_idx in mainchain:
                    continue
                if G.get_edge_data(main_idx, end_idx)['bondtype'] > 1:
                    extra_main.append(end_idx)
        mainchain += extra_main
        return mainchain

    def checkDegree(self, G, degree):
        '''Judge wheter the fragment contains branches

        Args:
            G (:obj:`networkx.classes.graph.Graph`): graph for the fragment
            degree (dict (int: int)): degree for each atom in the graph, G.degree()

        Returns:
            bool: whether the fragment contains branches
        '''
        for key, value in degree.items():
            if value == 3:
                if all([degree[n] != 1 for n in G.neighbors(key)]):
                    return True
            elif value == 4:
                if sum([degree[n] == 1 for n in G.neighbors(key)]) < 2:
                    return True
        return False

    def cutChain(self, chain):
        '''Cut the chain into fragments

        1. Do not further cut chains with less than four heavy atoms

        2. When there are four heavy atoms, we do not further cut if there is a tertiary
           functional group (a functional group that has an atom connected to three other
           atoms) or there is a non-carbon heavy atom that is not at the end of the chain,
           otherwise we cut from the center.

        3. When there are more than four heavy atoms:
            3a. For the bonds that connect to a 1-point score atom, cut the bond that connects 
                to a 0-point atom and retain all the rest. For example in the fragment
                C-C(0)..Y(1)-C=Y, “..” is the bond to cut while “-” are the bonds to retain.
                The numbers in the parenthesis are the scores of the corresponding atoms.
            3b. Check all atoms (with score higher than 1 point) from the highest score to the
                lowest score. Retain all non-single bonds connected to them. For all single
                bonds connected to them, cut the bond that connects to another high score
                atom (larger than 1 point).
            3c. Cut another single bond of a two-connection low score (score 0 or 1) atom
                that is connected to a high score atom. For example, in Y=C(7)-C(0)..C cut “..”.
            3d. When a low score atom is simultaneously connected to two high score atoms,
                cut the bond that connects to the high-score atom with the lowest score.
            3e. Merge the end atom to the nearest fragment. If an end atom connects to
                another 0-point atom, then these two atoms become a standalone elementary
                fragment.
            3f. Merge any left over single high score atom to the nearby highest score atom.
            3g. For the remaining uncut chains, cut every 3 heavy atoms so as to obtain an
                elementary fragment. If the number of atoms is not divisible by 3, first cut 2
                heavy atoms and continue with the procedure just mentioned.
            3h. Merge any non-end single-atom fragment to the nearest fragment. If there is
                a 0-point atom connected to it, then these two atoms become a standalone
                elementary fragment.

        Args:
            chain (:obj:`networkx.calsses.graph.Graph`): graph for the chain

        Returns:
            list of list of int: list of cut chain fragments
        '''
        chain_frags = []

        # 1.
        if len(chain.nodes()) < 4:
            chain_frags.append(list(chain.nodes()))
            return chain_frags

        # 2.
        degree = dict(chain.degree())
        if len(chain.nodes()) == 4:
            terminal_atom = []
            for k, v in degree.items():
                if v > 2:
                    chain_frags.append(list(chain.nodes()))
                    return chain_frags
                if v == 2 and chain.node[k]['atom'].score == 1:
                    chain_frags.append(list(chain.nodes()))
                    return chain_frags
                elif v == 1:
                    terminal_atom.append(k)
            one_path = nx.shortest_path(chain, terminal_atom[0], terminal_atom[1])
            chain_frags.append(one_path[:2])
            chain_frags.append(one_path[2:])
            return chain_frags

        # 3.
        # Modify and save bondtype property: -1, can be cut; 0, cannot be cut; >0, undetermined 
        index_score_dict = {node:chain.node[node]['atom'].score for node in chain.nodes() if chain.node[node]['atom'].score > 1}
        index_score_lt1_dict = {node:chain.node[node]['atom'].score for node in chain.nodes() if chain.node[node]['atom'].score < 2}
        index_score_list = sorted(index_score_dict.items(), key=lambda e:e[1], reverse=True)

        # 3a.
        for idx, val in index_score_lt1_dict.items():
            if val == 0:
                continue
            for nei in chain.neighbors(idx):
                if chain.node[nei]['atom'].score == 1:
                    chain.get_edge_data(idx, nei)['bondtype'] == -1
                elif chain.node[nei]['atom'].score == 0:
                    chain.get_edge_data(idx, nei)['bondtype'] == 0
        # 3b.
        for idx, _ in index_score_list:
            for end_idx in chain.neighbors(idx):
                if chain.get_edge_data(idx, end_idx)['bondtype'] == 1:
                    if chain.node[end_idx]['atom'].score < 2:
                        chain.get_edge_data(idx, end_idx)['bondtype'] = -1
                    else:
                        chain.get_edge_data(idx, end_idx)['bondtype'] = 0
                elif chain.get_edge_data(idx, end_idx)['bondtype'] > 1:
                    chain.get_edge_data(idx, end_idx)['bondtype'] = -1
        # 3c.
        set0_bonds = []
        for idx in index_score_lt1_dict:
            if degree[idx] == 2:
                neigs = list(chain.neighbors(idx))
                if chain.get_edge_data(idx, neigs[0])['bondtype'] == 1 and chain.get_edge_data(idx, neigs[1])['bondtype'] == -1:
                    set0_bonds.append((idx, neigs[0]))
                elif chain.get_edge_data(idx, neigs[0])['bondtype'] == -1 and chain.get_edge_data(idx, neigs[1])['bondtype'] == 1:
                    set0_bonds.append((idx, neigs[1]))
        for start, end in set0_bonds:
            chain.get_edge_data(start, end)['bondtype'] = 0
        # 3d.
        for idx in index_score_lt1_dict:
            neis = list(chain.neighbors(idx))
            if len(neis) == 2:
                per, bak = chain.node[neis[0]]['atom'].score, chain.node[neis[1]]['atom'].score
                if per > 1 and bak > 1:
                    if per > bak:
                        chain.get_edge_data(idx, neis[1])['bondtype'] = 0
                    else:
                        chain.get_edge_data(idx, neis[0])['bondtype'] = 0
                    
        # 3e.
        terminal_atom = [k for k, v in degree.items() if v == 1]
        for tm in terminal_atom:
            end_idx = list(chain.neighbors(tm))[0]
            if chain.get_edge_data(tm, end_idx)['bondtype'] == 0:
                chain.get_edge_data(tm, end_idx)['bondtype'] = -1
                # A terminal atom forms a fragment with a nearby 0-point C atom.
                if chain.node[end_idx]['atom'].score == 0:
                    for eend_idx in chain.neighbors(end_idx):
                        if eend_idx != tm:
                            chain.get_edge_data(end_idx, eend_idx)['bondtype'] = 0
        # 3f.
        for k in index_score_dict:
            neis = list(chain.neighbors(k))
            bond_flag = [chain.get_edge_data(k, nei)['bondtype'] == 0 for nei in neis]
            if all(bond_flag):
                topscore = -1
                topscore_idx = 0
                for nei in neis:
                    if chain.node[nei]['atom'].score > topscore:
                        topscore = chain.node[nei]['atom'].score
                        topscore_idx = nei
                chain.get_edge_data(k, topscore_idx)['bondtype'] = -1
        
        no_flag_atom = []
        for idx in chain.nodes():
            if chain.node[idx]['atom'].score < 2:
                check_bond = [chain.get_edge_data(idx, nei)['bondtype'] == 1 for nei in chain.neighbors(idx)]
                if any(check_bond):
                    no_flag_atom.append(idx)

        # 3g.
        sub_g = chain.subgraph(no_flag_atom)
        groups = list(nx.connected_components(sub_g))
        for grp in groups:
            chain_len = len(grp)
            if chain_len < 4:
                chain_frags.append(list(grp))
                continue
            tmp_subg = chain.subgraph(list(grp))
            tmp_degree = dict(tmp_subg.degree())
            terminal_atom = [k for k, v in tmp_degree.items() if v == 1]
            assert len(terminal_atom) == 2, "terminal_atom must be 2"
            one_path = nx.shortest_path(tmp_subg, terminal_atom[0], terminal_atom[1])
            if chain_len == 4:
                chain_frags += [one_path[:2], one_path[2:]]
                continue

            num, mod = divmod(chain_len, 3)
            if mod > 0:
                chain_frags.append(one_path[0:2])
            for i in range(num):
                if mod == 0:
                    tmp_frag = one_path[3*i:3*(i+1)]
                else:
                    tmp_frag = one_path[3*i+2:3*(i+1)+2]
                chain_frags.append(tmp_frag)

        # 3h.
        for node in chain.nodes():
            if chain.node[node]['atom'].score > 1:
                continue
            neigs = list(chain.neighbors(node))
            if len(neigs) != 2:
                continue
            if chain.get_edge_data(node, neigs[0])['bondtype'] == 0 and chain.get_edge_data(node, neigs[1])['bondtype'] == 0:
                for nei in neigs:
                    if chain.node[nei]['atom'].score == 0:
                        chain.get_edge_data(node, nei)['bondtype'] = -1
                        for nei2 in chain.neighbors(nei):
                            if nei2 != node:
                                chain.get_edge_data(nei, nei2)['bondtype'] = 0
                        break
                else:
                    chain.get_edge_data(node, neigs[0])['bondtype'] = -1

        # remove single chain nodes
        for idx in no_flag_atom:
            chain.remove_node(idx)
        # cut the bonds that are marked as "can be cut" 
        for start, end in list(chain.edges()):
            if chain.get_edge_data(start, end)['bondtype'] == 0:
                chain.remove_edge(start, end)
        for grp in nx.connected_components(chain):
            chain_frags.append(list(grp))

        return chain_frags

    def cutAllChain(self, rings):
        '''cut all chains in the molecule

        Ring and chain fragments are saved in the frag_rings and frag_chains attributes respectively.

        Args:
            rings (list of list of int): list of rings
        '''
        chains = self.getAllSingleChains(rings)
        for chain in chains:
            self.chainScore(chain)
            re_frags = self.cutChain(chain)
            self.frag_chains += re_frags

        # Merge halogen atoms into the chain fragments
        G = self.createGraph()
        for cha in self.frag_chains:
            for idx in list(cha):
                for nei in G.neighbors(idx):
                    if G.degree(nei) == 1 and G.node[nei]['atom'].symbol in halogen:
                        cha.append(nei)

halogen = {'F', 'Cl', 'Br', 'I', 'At'}
