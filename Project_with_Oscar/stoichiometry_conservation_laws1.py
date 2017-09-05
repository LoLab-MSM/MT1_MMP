import numpy as np
import sympy
import re
from collections import OrderedDict
from pysb.bng import generate_equations

#to get the stoichiometry matrix
def stoichiometry_matrix(model):
    generate_equations(model)
    sm = np.zeros((len(model.species), len(model.reactions)))
    for i_s, sp in enumerate(model.species):
        for i_r, r in enumerate(model.reactions):
            sm[i_s][i_r] = r['products'].count(i_s) - r['reactants'].count(i_s)
    return sm

#?????
def stoichimetry_matrix_passengers(model, pruned_system):
    generate_equations(model)
    sm = np.zeros((len(pruned_system.keys()), len(model.reactions)))
    for i_s, pa in enumerate(pruned_system.keys()):
        for i_r, r in enumerate(model.reactions):
            if r['rate'] in pruned_system[pa].as_coefficients_dict().keys():
                sm[i_s][i_r] = r['products'].count(pa) - r['reactants'].count(pa)
    return sm

#to get the rhs of mcl
def conservation_laws_values(model, conser_laws):
    if not isinstance(conser_laws, list):
        conser_laws = [conser_laws]

    initial_conditions_expanded = {}
    y0 = np.zeros(len(model.species))
    for cp, value_obj in model.initial_conditions:
        value = value_obj.value
        si = model.get_species_index(cp)
        y0[si] = value

    for spp in range(len(model.species)):
        initial_conditions_expanded['__s%d' % spp] = y0[spp]

    value_constants = {}
    for conser in conser_laws:
        constant_to_solve = [atom for atom in conser.atoms(sympy.Symbol) if re.match(r'[a-d]', str(atom))]
        solution = sympy.solve(conser, constant_to_solve)
        solution_ready = solution[0]
        solution_ready = solution_ready.subs(initial_conditions_expanded)
        value_constants[constant_to_solve[0]] = solution_ready

    return value_constants

#MAIN CODE to get the lhs and rhs of MCL
def conservation_relations(model, pruned_system=None):
    if pruned_system is not None:
        stoichiometry = stoichimetry_matrix_passengers(model, pruned_system)
        model_species = [model.species[i] for i in pruned_system.keys()]
    else:
        stoichiometry = stoichiometry_matrix(model)
        model_species = model.species
    
    #print stoichiometry
            #[[-1.  1.]
            #[-1.  1.]
            #[ 1. -1.]]
    sto_rank = np.linalg.matrix_rank(stoichiometry)
    #print sto_rank 
                #1
    species_info = OrderedDict()
    for sp in model_species:
        #print str(sp)
        species_info[str(sp)] = sympy.Symbol('__s%d' % model.get_species_index(sp))
    '''Var Data'''
    #print species_info
                    #OrderedDict([('A1(a1=None, a2=None)', __s0), 
                    #('A2(a1=None, a3=None)', __s1), 
                    #('A3(a2=None)', __s2), 
                    #('A1(a1=1, a2=None) % A1(a1=1, a2=None)', __s3), 
                    #('A1(a1=None, a2=1) % A2(a1=1, a3=None)', __s4), 
                    #('A2(a1=None, a3=1) % A3(a2=1)', __s5), 
                    #('A1(a1=1, a2=2) % A1(a1=1, a2=None) % A2(a1=2, a3=None)', __s6), 
                    #('A1(a1=1, a2=2) % A1(a1=1, a2=3) % A2(a1=2, a3=None) % A2(a1=3, a3=None)', __s7), 
                    #('A1(a1=None, a2=1) % A2(a1=1, a3=2) % A3(a2=2)', __s8), 
                    #('A1(a1=1, a2=2) % A1(a1=1, a2=None) % A2(a1=2, a3=3) % A3(a2=3)', __s9), 
                    #('A1(a1=1, a2=2) % A1(a1=1, a2=3) % A2(a1=3, a3=4) % A2(a1=2, a3=None) % A3(a2=4)', __s10), 
                    #('A1(a1=1, a2=2) % A1(a1=1, a2=3) % A2(a1=2, a3=4) % A2(a1=3, a3=5) % A3(a2=4) % A3(a2=5)', __s11)])

    sto_augmented = np.concatenate((stoichiometry, np.identity(stoichiometry.shape[0])), axis=1) #adding I at next column
    #print sto_augmented
                        #[[-1.  1.  1.  0.  0.]
                        #[-1.  1.  0.  1.  0.]
                        #[ 1. -1.  0.  0.  1.]]
    sto_augmented = sympy.Matrix(sto_augmented)
    sto_reduced = sto_augmented.rref()[0] #getting the reduced echelon form matrix
    #print sto_reduced
                    #In Matrix[[]] format.
                    #[1, -1, 0, 0, 1], 
                    #[0, 0, 1, 0, 1], 
                    #[0, 0, 0, 1, 1]
    #print stoichiometry.shape[1]
    conservation_matrix = sto_reduced[sto_rank:, stoichiometry.shape[1]:] #sto_reduced['rank menentukan mulai row brp', 'shape[1] == number of reactions in matrix menentukan mulai column brp'] 
                                                                            #sto_reduced[1:,2:]
    conservation_matrix = conservation_matrix.applyfunc(sympy.Integer)
    '''Conservation Matrix'''
    #print conservation_matrix
                            #Matrix([[1, 0, 0, 2, 1, 0, 2, 2, 1, 2, 2, 2], 
                            #[0, 1, 0, 0, 1, 1, 1, 2, 1, 1, 2, 2], 
                            #[0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 2]])
    conservation_laws = conservation_matrix.dot(species_info.values())
    '''The MCL Eqs'''
    #print conservation_laws
                        #[__s0 + 2*__s10 + 2*__s11 + 2*__s3 + __s4 + 2*__s6 + 2*__s7 + __s8 + 2*__s9,
                        # __s1 + 2*__s10 + 2*__s11 + __s4 + __s5 + __s6 + 2*__s7 + __s8 + __s9, 
                        #__s10 + 2*__s11 + __s2 + __s5 + __s8 + __s9]
    if not isinstance(conservation_laws, list):
        conservation_laws = [conservation_laws]

    for ii, cl in enumerate(conservation_laws):
        conservation_laws[ii] = cl - sympy.Symbol('a%d' % ii)

    value_constants = conservation_laws_values(model, conservation_laws)

    return conservation_laws, value_constants