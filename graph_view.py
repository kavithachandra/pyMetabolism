"""
Multiple ways of creating a graph representation of a L{Metabolism} instance.
Generation of a stoichiometric matrix is also possible.

@todo: set up logging, complete docs
"""


import networkx
from numpy import hstack, vstack, zeros


class BipartiteMetabolicNetwork(networkx.DiGraph):
    """
A bipartite graph representation of a L{Metabolism} instance.

@todo: many things
    """
    
    def __init__(self, reactions=None, name='', weighted=True, node_info=False,
        *args, **kwargs):
        """
@param reactions: A list of L{Reaction}s from which the graph can be populated.
@type reactions: C{list}

@param name: A simple string to identify the graph if desired.
@type name: C{str}

@param weighted: By default the stoichiometric factors are used as edge weights,
                this can be turned off.
@type weighted: C{bool}

@param node_info: If true, information will be added to the graph's nodes.
                 For compound nodes the following information will be tried in
                 order:
                 
                    1. L{Compound.charge}
                    2. L{Compound.mass}
                 
                 For reaction nodes L{Reaction.rate} will be used.
@type node_info: C{bool}
        """
        super(BipartiteMetabolicNetwork, self).__init__(name=str(name), *args,
            **kwargs)
        if reactions:
            self._populate_graph_on_init(reactions, weighted, node_info)
    
    def _populate_graph_on_init(self, reactions, weighted, node_info):
        for rxn in reactions:
            self.add_reaction(rxn)
    
    def add_reaction(self, rxn, weighted=True, node_info=False):
        """
docstring for add_reaction
        """
        subEdges = [(elem, reaction, reaction.getStoichiometricFactor(elem)) for elem in reaction.substrates]
        prodEdges = [(reaction, elem, reaction.getStoichiometricFactor(elem)) for elem in reaction.products]
        if reaction.reversibleQ:
            return subEdges + prodEdges + [(j, i, d) for i, j, d in subEdges] + [(j, i, d) for i, j, d in prodEdges]
        else:
            return subEdges + prodEdges


class MetaboliteCentricNetwork(networkx.DiGraph):
    """"""
    def __init__(self, reactions=None, name='', weighted=True, *args, **kwargs):
        super(MetaboliteCentricNetwork, self).__init__(name=str(name), *args, **kwargs)
        if reactions:
            self.__populate_graph_on_init(reactions)
    
    def __populate_graph_on_init(self, reactions):
        for reac in reactions:
            self.add_reaction(reac)
    
    def add_reaction(self, reaction):
        """docstring for add_reaction"""
        edges = [(substr, prod, reaction) for prod in reaction.products for substr in reaction.substrates]
        if reaction.reversible:
            self.add_edges_from(edges + [(j, i, d) for i, j, d in edges])
        else:
            self.add_edges_from(edges)


class StoichiometricMatrix(object):
    """A class representing a stoichiometric matrix.
    
    Columns represent reactions
    Rows represent compounds
    Coefficients ...
    """
    def __init__(self, *args, **kwargs):
        super(StoichiometricMatrix, self).__init__(*args, **kwargs)
        self.matrix = None
        self.compound_map = dict()
        self.reaction_map = dict()
        
    def add_stoichiometry_from(self, metabolism):
        """"""
        j = len(self.reaction_map)
        i = len(self.compound_map)
        for reaction in metabolism.reactions:
            if reaction not in self.reaction_map:
                if self.matrix == None:
                    i = self._init_matrix(reaction)
                    j = 1
                    continue
                self.reaction_map[reaction] = j
                j += 1
                self.matrix = hstack((self.matrix, zeros((self.matrix.shape[0],
                    1))))
                for compound in reaction:
                    if compound not in self.compound_map:
                        self.compound_map[compound] = i
                        i += 1
                        self.matrix = vstack((self.matrix, zeros((1,
                            self.matrix.shape[1]))))
                    self.matrix[self.compound_map[compound]]\
                        [self.reaction_map[reaction]] =\
                        reaction.get_stoich_coeff(compound)
    
    def _init_matrix(self, reaction):
        num = len(reaction)
        self.matrix = zeros((num, 1))
        for i, compound in enumerate(reaction):
            self.compound_map[compound] = i
            self.matrix[i][0] = reaction.get_stoich_coeff(compound)
        self.reaction_map[reaction] = 0
        return num
