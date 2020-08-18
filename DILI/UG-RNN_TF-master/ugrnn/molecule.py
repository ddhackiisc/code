from __future__ import absolute_import
from __future__ import division

import logging
from rdkit import Chem

import networkx as nx
import numpy as np
import six

logger = logging.getLogger(__name__)


class Molecule:
    max_number_of_parents = 4

    def __init__(self, smile, logp=None, contract_rings=False):
        self.smile = smile
        self.logp = logp
        # logger.info("Parsing Molecule {:},contract rings: {:}".format(smile, contract_rings))
        self.atoms = []
        m = Chem.MolFromSmiles(smile)  # makes a mol object
        # Chem.Kekulize(self.m)
        self.no_of_atoms = m.GetNumAtoms()
        self.graph = nx.Graph()  # creates an empty graph for the molecule

        for i in xrange(self.no_of_atoms):
            atom = m.GetAtomWithIdx(i)
            # adds a node for each atom and atom features returns [1x5] array
            self.graph.add_node(i, attr_dict={"atom_features": Molecule.atom_features(atom)})
            for neighbour in atom.GetNeighbors():  # iterates through each neighbour of the atom
                neighbour_idx = neighbour.GetIdx()
                # returns bond index of the coressponding bond
                bond = m.GetBondBetweenAtoms(i, neighbour_idx)
                # adds an edge to the graph and bond features returns [1x6] array with bond info
                self.graph.add_edge(i, neighbour_idx,
                                    attr_dict={"bond_features": Molecule.bond_features(bond)})

        if contract_rings:
            # makes the rings of the grpahs into a single node as explaines in the paper and relabels the nodes
            self.reduce_graph_rings()

        # creates an array of directed graphs from the undirected graph  taking each node as sink node
        self.create_directed_graphs()
        #  creates a local input vector with all teh graph atom and bond features
        self.create_feature_vectors()

    def create_directed_graphs(self):
        '''
        :return:
        '''
        self.directed_graphs = np.empty(
            (self.no_of_atoms, self.no_of_atoms - 1, 3), dtype=int)

        # parse all the atoms one by one and get directed graph to that atom
        # as the sink node
        for idx in xrange(self.no_of_atoms):
            # get shortest path from the root to all the other atoms and then reverse the edges.
            # return a dictionary of paths keyed by target name. Each path is am iterable with the nodes of the path
            # mentioned in order
            path = nx.single_source_dijkstra_path(self.graph, idx)
            G = nx.DiGraph()
            # goes to each of the nodes of the graph and add it t the current directed graph
            for i in xrange(self.no_of_atoms):
                temp = path[i]
                temp.reverse()  # reverses egde direction of that particular   path
                G.add_path(temp)

            # do a topological sort to get a order of atoms with all edges pointing to the root
            topological_order = nx.topological_sort(G)

            sorted_path = np.empty((self.no_of_atoms - 1, 3))

            no_of_incoming_edges = {}
            # slightly comlicated loop but in effect it returns a 3D array with frat dimenstion
            # itersting through the root node second dimension iterating through each node in the coressponding
            # directed graph and the thord dimesion is 1x3 array
            for i in xrange(self.no_of_atoms - 1):
                node = topological_order[i]
                edge = (nx.edges(G, node))[0]
                if edge[1] in no_of_incoming_edges:
                    index = no_of_incoming_edges[edge[1]]
                    no_of_incoming_edges[edge[1]] += 1
                else:
                    index = 0
                    no_of_incoming_edges[edge[1]] = 1
                sorted_path[i, :] = [node, edge[1], index]
            self.directed_graphs[idx, :, :] = sorted_path

    def create_feature_vectors(self):
        '''
        :return:
        '''
        # create a three dimesnional matrix I,
        # such that Iij is the local input vector for jth vertex in ith DAG

        length_of_bond_features = Molecule.num_bond_features()
        length_of_atom_features = Molecule.num_atom_features()

        self.local_input_vector = np.zeros(
            (self.no_of_atoms, self.no_of_atoms, Molecule.num_of_features()))

        for idx in xrange(self.no_of_atoms):
            sorted_path = self.directed_graphs[idx, :, :]  # selects each graph based on a root node

            self.local_input_vector[idx, idx, :length_of_atom_features] = \
                self.get_atom_features(
                    idx)  # adds the atom features of the root node to the locl input cevtor

            no_of_incoming_edges = {}
            for i in xrange(self.no_of_atoms - 1):
                node1 = sorted_path[i, 0]  # node1 is the particular node
                node2 = sorted_path[i, 1]  # node2 is the parent node of node1

                self.local_input_vector[idx, node1, :length_of_atom_features] = \
                    self.get_atom_features(node1)  # adds the atom features of the particular node

                if node2 in no_of_incoming_edges:
                    index = no_of_incoming_edges[node2]
                    no_of_incoming_edges[node2] += 1
                    if index >= Molecule.max_number_of_parents:
                        continue
                else:
                    index = 0
                    no_of_incoming_edges[node2] = 1

                start = length_of_atom_features + index * length_of_bond_features
                end = start + length_of_bond_features

                self.local_input_vector[idx, node2, start:end] = \
                    self.get_bond_features(node1, node2)

    def get_cycle(self):
        try:
            return nx.find_cycle(self.graph)
        except:
            return []

    def reduce_graph_rings(self):
        '''
        :return:
        '''
        cycle_name_format = "R_{:}"  # format to name the contracted rings
        index = 0
        # return array of egdes of cycle. each edge is represented by an order pair of node numbers
        cycle = self.get_cycle()

        while cycle:  # while the cycle array is not empty
            cycle_name = cycle_name_format.format(index)
            self.graph.add_node(cycle_name)

            # ebunch = zip(cycle, (cycle[1:] + cycle[:1]))
            self.graph.remove_edges_from(cycle)

            for node1, node2 in cycle:
                if isinstance(node1, six.string_types):  # checks if node1 is a ring
                    self.graph.add_edge(node1, cycle_name,  # adds an adge between rings and bond_feaures returns an array with first entry 1 and others 0
                                        attr_dict={"bond_features": Molecule.bond_features_between_contract_rings()})
                    continue

                neighbours = self.graph.neighbors(node1)
                if not neighbours:  # means neighbours are also part of cycle
                    continue
                for neighbour in neighbours:
                    edge_attrs = self.get_bond_features(neighbour, node1)
                    self.graph.add_edge(neighbour, cycle_name, attr_dict={
                        "bond_features": edge_attrs})  # adds an edge between the cycele and neighbours
                    # removes the edge between neighbour and cycle members
                    self.graph.remove_edge(node1, neighbour)

            nx.set_node_attributes(self.graph, "atom_features",  # adds
                                   values={cycle_name: Molecule.atom_features_of_contract_rings(0)})

            for node1, node2 in cycle:
                if not isinstance(node1, six.string_types):  # check ifnode1 is not a ring,
                    self.graph.remove_node(node1)  # then removes it from the graph
            index += 1
            cycle = self.get_cycle()  # reinitializes the cycle array

        self.graph = nx.convert_node_labels_to_integers(self.graph,  # relabels thegraph nodes consecutively
                                                        first_label=0)

        nx.draw(self.graph)
        self.no_of_atoms = len(self.graph)

    def get_atom_features(self, node_id):
        attrs = nx.get_node_attributes(self.graph, "atom_features")
        return attrs[node_id]

    def get_bond_features(self, node1, node2):
        attrs = self.graph.get_edge_data(node1, node2)
        return attrs["bond_features"]

    # returns an array [atom type, no.of neighbours, No.of H, implicit valence, aromaticity]
    @staticmethod
    def atom_features(atom):
        return np.array(Molecule.one_of_k_encoding_unk(atom.GetSymbol(),
                                                       ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl',
                                                        'Br', 'Mg', 'Na',
                                                        'Ca', 'Fe', 'As', 'Al', 'I', 'B', 'V', 'K',
                                                        'Tl', 'Yb',
                                                        'Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti',
                                                        'Zn', 'H',  # H?
                                                        'Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In',
                                                        'Mn',
                                                        'Zr', 'Cr', 'Pt', 'Hg', 'Pb', 'Unknown']) +
                        Molecule.one_of_k_encoding(atom.GetDegree(), [0, 1, 2, 3, 4, 5]) +
                        Molecule.one_of_k_encoding_unk(atom.GetTotalNumHs(), [0, 1, 2, 3, 4]) +
                        Molecule.one_of_k_encoding_unk(atom.GetImplicitValence(),
                                                       [0, 1, 2, 3, 4, 5]) +
                        [atom.GetIsAromatic()])

    @staticmethod
    def atom_features_of_contract_rings(degree):
        return np.array(Molecule.one_of_k_encoding_unk('Unknown',
                                                       ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl',
                                                        'Br', 'Mg', 'Na',
                                                        'Ca', 'Fe', 'As', 'Al', 'I', 'B', 'V', 'K',
                                                        'Tl', 'Yb',
                                                        'Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti',
                                                        'Zn', 'H',  # H?
                                                        'Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In',
                                                        'Mn', 'Zr',
                                                        'Cr', 'Pt', 'Hg', 'Pb', 'Unknown']) +
                        Molecule.one_of_k_encoding(degree, [0, 1, 2, 3, 4, 5]) +
                        Molecule.one_of_k_encoding_unk(0, [0, 1, 2, 3, 4]) +
                        Molecule.one_of_k_encoding_unk(0, [0, 1, 2, 3, 4, 5]) +
                        [0])

    @staticmethod
    def bond_features_between_contract_rings():
        return np.array([1, 0, 0, 0, 0, 0])

    # Returns array with bond info
    @staticmethod
    def bond_features(bond):
        bt = bond.GetBondType()
        return np.array([bt == Chem.rdchem.BondType.SINGLE,
                         bt == Chem.rdchem.BondType.DOUBLE,
                         bt == Chem.rdchem.BondType.TRIPLE,
                         bt == Chem.rdchem.BondType.AROMATIC,
                         bond.GetIsConjugated(),
                         bond.IsInRing()])

    @staticmethod
    def num_of_features():
        return Molecule.max_number_of_parents*Molecule.num_bond_features() + Molecule.num_atom_features()

    @staticmethod
    def one_of_k_encoding(x, allowable_set):
        if x not in allowable_set:
            raise Exception(
                "input {0} not in allowable set{1}:".format(x, allowable_set))
        return map(lambda s: x == s, allowable_set)

    @staticmethod
    def one_of_k_encoding_unk(x, allowable_set):
        """Maps inputs not in the allowable set to the last element."""
        if x not in allowable_set:
            x = allowable_set[-1]
        return map(lambda s: x == s, allowable_set)

    @staticmethod
    def num_atom_features():
        # Return length of feature vector using a very simple molecule.
        m = Chem.MolFromSmiles('CC')
        alist = m.GetAtoms()
        a = alist[0]
        return len(Molecule.atom_features(a))

    @staticmethod
    def num_bond_features():
        # Return length of feature vector using a very simple molecule.
        simple_mol = Chem.MolFromSmiles('CC')
        Chem.SanitizeMol(simple_mol)
        return len(Molecule.bond_features(simple_mol.GetBonds()[0]))


if __name__ == '__main__':
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_format)
    logger = logging.getLogger(__name__)
    m = Molecule("c2(Cl)c(Cl)c(Cl)c1nccnc1c2(Cl)", True)
