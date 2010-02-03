"""
Basic classes modelling compounds, reactions, and metabolism.

@author: Nikolaus Sonnenschein
@author: Moritz Beber
@contact: niko(at)foo.org
@contact: moritz(at)foo.org
@copyright: Jacobs University Bremen. All rights reserved.
@since: 2009-09-11
@todo: set up logging

"""

class Compound(object):
    """
    A class modeling a chemical compound. Primary identifier of L{Compound}s is a
    simple string C{identifier} although many synonymous identifiers may be set up.
    
    @cvar _compound_memory: A dictionary that stores instances of L{Compound} as
                            values to their C{identifier} key.
    @type _compound_memory: C{dict}
    
    @ivar identifier: The main name of the compound will be used for comparisons and
                      representation.
    @type identifier: C{str}
    
    @ivar compartment: Define a compartment where this compound occurs for a multi-
                       compartmentalised metabolism.
    @type compartment: C{str}
    
    @ivar synonyms: Alternative names for the compound.
    @type synonyms: C{list}
    
    @ivar formula: Elemental formula of the compound.
    @type formula: C{str}
    
    @ivar in_chi:
    @type in_chi: C{str}
    
    @ivar in_chey_key:
    @type in_chey_key: C{str}
    
    @ivar smiles:
    @type smiles: C{str}
    
    @ivar charge: Electric charge of the compound.
    @type charge: C{int}
    
    @ivar mass: Mass of the compound.
    @type mass: C{float}
    @note: Currently, no unit system for mass is implemented it is left up to the
           user.
    
    @todo: clarify compartment use, describe synonymous strings better
    """
    
    _compound_memory = dict()
    
    def __new__(cls, identifier, compartment=None, synonyms=None,
        formula=None, in_chi=None, in_chey_key=None, smiles=None, charge=None,
        mass=None, *args, **kwargs):
        """
        @return: Either returns an old L{Compound} instance if the name already exists
        or passes a new L{Compound} C{class instance} to be initialised.
        @rtype: L{Compound} C{class instance}

        @attention: This method is never called directly.
        """
        identifier = str(identifier)
        if identifier in cls._compound_memory:
            return cls._compound_memory[identifier]
        else:
            instance = super(Compound, cls).__new__(cls, *args, **kwargs)
            return instance
    
    def __init__(self, identifier, compartment=None, synonyms=None,
        formula=None, in_chi=None, in_chey_key=None, smiles=None, charge=None,
        mass=None, *args, **kwargs):
        """
        Either does nothing if the L{Compound} instance already exists or
        intialises a new L{Compound} instance.
        """
        if identifier in self.__class__._compound_memory:
            return None
        super(Compound, self).__init__(*args, **kwargs)
        self.identifier = identifier
        self.compartment = compartment
        self.synonyms = synonyms
        self.formula = formula
        self.in_chi = in_chi
        self.in_chey_key = in_chey_key
        self.smiles = smiles
        if charge:
            self.charge = int(charge)
        else:
            self.charge = None
        if mass:
            self.mass = float(mass)
        else:
            self.mass = None
        self.__class__._compound_memory[self.identifier] = self
        
    def __str__(self):
        """
        @rtype: C{str}
        """
        if self.compartment:
            return self.identifier + '<' + self.compartment + '>'
        else:
            return self.identifier
        
    def __contains__(self, element):
        """
        Checks for atomic element in compound.

        @raise NotImplementedError:
        """
        raise NotImplementedError
    
    def __hash__(self):
        """
        @rtype: C{int}
        """
        return hash(self.identifier)
    
    def __cmp__(self, other):
        """
        @rtype: C{int}
        """
        return cmp(id(self), id(other))
    
    def get_id(self):
        """
        @return: Returns identifier of compound.
        @rtype: C{str}        
        """
        return self.identifier


class Reaction(object):
    """
    A class modeling a chemical reaction.

    @cvar _reaction_memory: A dictionary that stores instances of L{Reaction} as
                            values to their C{identifier} key.
    @type _reaction_memory: C{dict}

    @ivar identifier: The main name of the reaction will be used for comparisons and
                      representation.
    @type identifier: C{str}

    @ivar substrates: L{Compound}s being consumed by this reaction.
    @type substrates: C{tuple}

    @ivar products: L{Compound}s being produced by this reaction.
    @type products: C{tuple}

    @ivar stoichiometry: Tuple containing the positive integer stoichiometric
                         factors of all substrates and products in the order
                         supplied to C{substrates} and C{products}.
    @type stoichiometry: C{tuple}
    @attention: The stoichiometric factors given in C{stoichiometry} should be
                positive integers or floats, their sign will be adjusted if necessary.

    @ivar synonyms: Synonymous names for this reaction.
    @type synonyms: C{list}

    @ivar rate: Define the value of the rate constant for this reaction.
    @type rate: C{float}

    @ivar stoichiometry_dict: Dictionary providing easy access to stoichiometric
                              factors (values) of all compounds (keys) in the
                              reaction.
    @type stoichiometry_dict: C{dict}

    @ivar reversible: Specify whether the reaction is considered reversible.
    @type reversible: C{bool}
    
    @todo: Fix stoichiometry issues ...
    """
    _reaction_memory = dict()
    
    def __new__(cls, identifier, substrates, products, stoichiometry,
        reversible=False, synonyms=None, rate=None, *args, **kwargs):
        """
        @return: Either returns an old L{Reaction} instance if the name already exists
        or passes a new L{Reaction} C{class instance} to be initialised.
        @rtype: L{Reaction} C{class instance}

        @attention: This method is never called directly.
        """
        identifier = str(identifier)
        if identifier in cls._reaction_memory:
            return cls._reaction_memory[identifier]
        else:
            instance = super(Reaction, cls).__new__(cls, *args, **kwargs)
            return instance
    
    def __init__(self, identifier, substrates, products, stoichiometry,
        reversible=False, synonyms=None, rate=None, *args, **kwargs):
        """
        Either does nothing if the L{Reaction} instance already exists or
        intialises a new L{Reaction} instance.
        """
        if identifier in self.__class__._reaction_memory:
            return None
        super(Reaction, self).__init__(*args, **kwargs)
        self.identifier = identifier
        self.substrates = tuple(substrates)
        self.products = tuple(products)
        self.stoichiometry = tuple([abs(coeff) for coeff in stoichiometry])
        self.stoichiometry_dict = dict(zip(list(self.substrates) + list(self.products), self.stoichiometry))
        self.reversible = bool(reversible)
        if synonyms:
            self.synonyms = synonyms
        else:
            self.synonyms = None
        if rate:
            self.rate = float(rate)
        else:
            self.rate = None
        self._consistency_check()
        self.__class__._reaction_memory[self.identifier] = self
    
    def _consistency_check(self):
        """
        Asserts some basic consistency of the L{Reaction} instance.

            1. The number of substrates plus products equals the number of
               stoichiometric factors.

            2. With enough meta data (L{Compound} formula, charge, or mass) stoichiometric
               balancing is checked.

        @raise AssertionError: If the reaction is not well balanced.

        @todo: Elemental balancing.
        """
        assert (len(self.products) + len(self.substrates)) ==\
            len(self.stoichiometry), "The number of stoichimetric factors does"\
            " not match the number of compounds."
        # elemental balancing
        check = True
        for compound in (self.substrates + self.products):
            if not compound.formula:
                check = False
                break
        if check:
            pass
        # mass balancing
        check = True
        sum_subs = 0.
        sum_prods = 0.
        for compound in self.substrates:
            if not compound.mass:
                check = False
                break
            else:
                sum_subs += float(stoichiometry_dict[compound]) * compound.mass
        if check:
            for compound in self.products:
                if not compound.mass:
                    check = False
                    break
                else:
                    sum_prods += float(stoichiometry_dict[compound]) *\
                        compound.mass
        if check:
            assert sum_subs == sum_prods, "There is a mass imbalance in"\
                " reaction '%s'" % self.identifier
        # charge balancing
        check = True
        sum_subs = 0
        sum_prods = 0
        for compound in self.substrates:
            if not compound.charge:
                check = False
                break
            else:
                sum_subs += stoichiometry_dict[compound] * compound.charge
        if check:
            for compound in self.products:
                if not compound.charge:
                    check = False
                    break
                else:
                    sum_prods += stoichiometry_dict[compound] * compound.charge
        if check:
            assert (sum_subs + sum_prods) == 0, "There is a charge imbalance in"\
                " reaction '%s'" % self.identifier
    
    def __iter__(self):
        return (compound for compound in self.substrates + self.products)
    
    def __str__(self):
        """
        @return: A C{str} representation of the reaction, e.g., 2 A + 4 B -> 1 C or
                 2 A + 4 B <=> 1 C for a reversible reaction.
        @rtype: C{str}
        """
        def util(compound_list):
            reaction_str = list()
            for compound in compound_list:
                reaction_str.append(str(abs(self.stoichiometry_dict[compound])))
                reaction_str.append(compound.__str__())
                if not (compound == compound_list[-1]):
                    reaction_str.append('+')
            return reaction_str
        reaction_str = list()
        reaction_str.extend(util(self.substrates))
        if self.reversible:
            reaction_str.append('<=>')
        else:
            reaction_str.append('->')
        reaction_str.extend(util(self.products))
        return ' '.join(reaction_str)
    
    def __contains__(self, compound):
        """
        Checks for the presence of compound in the reaction.

        @param compound: Presence tested for.
        @type compound: L{Compound} or C{str}
        @rtype: C{bool}
        """
        if isinstance(compound, str):
            return compound in [c.get_id() for c in self.get_compounds()]
        if isinstance(compound, Compound):
            return compound in self.get_compounds()
        else:
            return False
    
    def __len__(self):
        """
        @rtype: C{int}
        """
        return len(self.substrates) + len(self.products)
    
    def __cmp__(self, other):
        """
        @rtype: C{int}
        """
        return cmp(id(self), id(other))
        
    def get_id(self):
        """
        @return: Identifier or reaction
        @rtype: C{str}
        """
        return self.identifier
        
    def get_compounds(self):
        """
        @return: Return a list of all L{Compound}s participating in this reaction.
        @rtype: C{list}
        """
        return list(self.substrates + self.products)
    
    def get_stoich_coeff(self, compound):
        """
        @param compound: The compound whose stoichiometric factor is queried.
        @type compound: L{Compound} or C{str}
        @return: Return the stoichiometric coefficient of a compound.
        @rtype: C{int}
        @raise KeyError: If C{compound} is not contained in the reaction.
        """
        if compound in self:
            if isinstance(compound, str):
                compound = Compound(compound)
            if compound in self.substrates:
                return -self.stoichiometry_dict[compound]
            elif compound in self.products:
                return self.stoichiometry_dict[compound]
        else:
            msg = "'%s' is not participating in reaction '%s'" % (compound,
                reaction.identifier)
            raise KeyError
            
    def reversibleQ(self):
        """docstring for reversiblQ"""
        return self.reversible


class Metabolism(object):
    """
    Implements the representation of a metabolism.

    @cvar _metabolism_memory: A dictionary that stores instances of L{Metabolism} as
                            values to their C{name} key.
    @type _metabolism_memory: C{dict}

    @cvar _metabolism_counter: A counter of global existing L{Metabolism} instances.
    @type _metabolism_counter: C{int}

    @ivar reactions: A list of all reactions.
    @type reactions: C{list}

    @ivar name: Identifier of the metabolism, automatically set if not provided.
    @type name: C{str}

    @ivar compounds: A set of all compounds found in all the reactions.
    @type compounds: C{set}

    @ivar reactions_dict: Convenience dictionary allows access to reactions (values)
                          by their identifier (key).
    @type reactions_dict: C{dict}
    """
    _metabolism_memory = dict()
    _metabolism_counter = 0
    
    def __new__(cls, reactions=None, name='', *args, **kwargs):
        """
        @return: Either returns an old L{Metabolism} instance if the name already exists
        or passes a new L{Metabolism} C{class instance} to be initialised.
        @rtype: L{Metabolism} C{class instance}

        @attention: This method is never called directly.
        """
        if name:
            name = str(name)
        if name in cls._metabolism_memory:
            return cls._metabolism_memory[name]
        else:
            instance = super(Metabolism, cls).__new__(cls, *args, **kwargs)
            return instance
    
    def __init__(self, reactions=None, name=None, *args, **kwargs):
        """
        Either does nothing if the L{Metabolism} instance already exists or
        intialises a new L{Metabolism} instance.
        """
        if name in self.__class__._metabolism_memory:
            return None
        super(Metabolism, self).__init__(*args, **kwargs)
        if name:
            self.name = name
        else:
            self.name = "Metabolism-%d" % self.__class__._metabolism_counter
            self.__class__._metabolism_counter += 1
        if reactions:
            self.reactions = list(reactions)
        self.compounds = set()
        for rxn in self.reactions:
            self.compounds.update(rxn.get_compounds()) 
        self.reactions_dict = dict([(rxn.identifier, rxn) for rxn in
            self.reactions])
        self.__class__._metabolism_memory[self.name] = self
    
    def __str__(self):
        """
        @return: Provides some statistics about the system e.g. no. of reactions.
        @rtype: C{str}
        """
        info = "System name: '%s'\nNumber of reactions: %i\nNumber of compounds: %i" % (self.name, len(self), len(self.compounds))
        return info 
    
    def __len__(self):
        """
        @return: Returns the number of reactions.
        @rtype: C{int}
        """
        return len(self.reactions)
    
    def __contains__(self, reaction):
        """
        Checks for the presence of a reaction in the metabolism.

        @param reaction: Presence tested for.
        @type reaction: L{Reaction} or C{str}
        @rtype: C{bool}
        """
        if isinstance(reaction, str):
            return reaction in [r.get_id() for r in self.get_reactions()]
        if isinstance(reaction, Reaction):
            return reaction in self.get_reactions()
        else:
            return False
    
    def __getitem__(self, rxn):
        """
        @param rxn: The reaction to be returned.
        @type rxn: C{str} or C{int}

        @return: Returns a reaction either by string identifier or index.
        @rtype: L{Reaction}

        @raise IndexError: If C{rxn} is an integer out of bounds.
        @raise KeyError: If C{rxn} is a string and not present in the C{reactions_dict}.
        """
        if isinstance(rxn, int):
            return self.reactions[rxn]
        if isinstance(rxn, str):
            return self.reactions_dict[rxn]
    
    def __cmp__(self, other):
        """
        @rtype: C{int}
        """
        return cmp(id(self), id(other))
    
    def get_compounds(self):
        """
        @return: Return a list of all L{Compound}s.
        @rtype: C{list}
        """
        return list(self.compounds)
        
    def get_reactions(self):
        """
        @return: Return a list of all L{Reaction}s.
        @rtype: C{list}
        """
        return list(self.reactions)
        
