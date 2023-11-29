#!/usr/bin/env python
# coding=utf-8

from builtins import next
from builtins import object

class BaseMixin(object):
    def findSmallestRings(self, get_ring_bonds=False):
        '''return the smallest rings set'''
        import networkx as nx
        G = self.getGraph()
        all_rings = set()
        if get_ring_bonds:
            ring_bonds = set()
        for edge in G.edges():
            G.remove_edge(*edge)
            try:
                for path in nx.all_shortest_paths(G, *edge):
                    idx0 = path.index(min(path))
                    idx1 = (idx0+1) % len(path)
                    if path[idx0-1] > path[idx1]:
                        all_rings.add(tuple(path[idx0:]+path[:idx0]))
                    else:
                        all_rings.add(tuple(reversed(path[idx1:]+path[:idx1])))
            except nx.NetworkXNoPath:
                pass
            else:
                if get_ring_bonds:
                    ring_bonds.add(edge)
            G.add_edge(*edge)
        if get_ring_bonds:
            return all_rings, ring_bonds
        return all_rings

    def getSubstructMatch(self, rdmol, matchBondType=True):
        if matchBondType:
            return self.rdmol.GetSubstructMatch(rdmol)
        else:
            import networkx as nx
            G1 = createGraph(nx, self.rdmol)
            G2 = createGraph(nx, rdmol)
            gm = nx.algorithms.isomorphism.GraphMatcher(
                G1, G2, node_match=lambda n1, n2: n1['atom'] == n2['atom'])
            match = []
            if gm.subgraph_is_isomorphic():
                match_dict = next(gm.subgraph_isomorphisms_iter())
                match = [k for k, _ in sorted(match_dict.items(),
                                              key=lambda x: x[1])]

            return tuple(match)

    def getGraph(self, removeHs=False):
        import networkx as nx
        return createGraph(nx, self.rdmol, removeHs=removeHs)

def createGraph(nx, rdmol, removeHs=False, match_charge=False):
    G = nx.Graph()
    for atm in rdmol.GetAtoms():
        if match_charge:
            G.add_node(atm.GetIdx(), atom="%s_%s" %
                       (atm.GetSymbol(), atm.GetFormalCharge()))
        else:
            G.add_node(atm.GetIdx(), atom=atm.GetSymbol())
    for bond in rdmol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(),
                   bond_type="%s" % bond.GetBondType())
    if removeHs:
        for atm in rdmol.GetAtoms():
            if atm.GetSymbol() == 'H':
                G.remove_node(atm.GetIdx())
    if not hasattr(G, 'node'):
        G.node = G.nodes

    return G
