import mt1_mmp_model
import pysb.bng
import sympy
import networkx
import itertools
from stoichiometry_conservation_laws import conservation_relations, conservation_laws_values
from sympy.core.relational import Equality

'''???????'''
def do(self, e, i=None): 
    """do `e` to both sides of self using function given or
    model expression with a variable representing each side:
    >> eq.do(i + 2)  # add 2 to both sides of Equality
    >> eq.do(i + a, i)  # add `a` to both sides; 3rd parameter identifies `i`
    >> eq.do(lambda i: i + a)  # modification given as a function
    """
    if isinstance(e, (sympy.FunctionClass, sympy.Lambda, type(lambda: 1))):
        return self.applyfunc(e)
    e = sympy.S(e)
    i = (set([i]) if i else e.free_symbols) - self.free_symbols
    if len(i) != 1:
        raise ValueError('not sure what symbol is being used to represent a side')
    i = i.pop()
    f = lambda side: e.subs(i, side)
    return self.func(*[f(side) for side in self.args])


Equality.do = do

'''To import mt1_mmp system'''
model = mt1_mmp_model.return_model('original')
'''To generate ODEs'''
pysb.bng.generate_equations(model)
#print len(model.reactions)
#print len(model.species)
#print model.species
#print model.reactions

# this gets the conservation laws (cl) and their constant values (cl_values)
'''Go to stoichiometry_conservation_laws.py'''
cl, cl_values = conservation_relations(model) 
'''From this, we have MCL'''


# this loops over the model reactions and storage in rcts_rules the reverse and no reverse reaction
# for each rule
rcts_rules = {}
for rule in model.rules:
    '''To storage reactants of reaction for each rule'''
    #untuk setiap rule, reaksi yang yang memiliki rule yg sama akan disimpan di rule_reactants.
    rcts_no_reverse = 0
    rcts_reverse = 0
    rule_reactants = []
    for reaction in model.reactions:
        if reaction['rule'][0] == rule.name and reaction['reverse'][0] == False:
            rule_reactants.append(reaction['reactants'])
    '''Printed for each rule'''
    #print rule_reactants 
                        #Rule 'ab': [(0, 1), (0, 4), (0, 7), (0, 9), (0, 10)]
                        #Rule 'bc': [(1, 2), (1, 5), (2, 3), (3, 5), (1, 7), (1, 8), (3, 7), (3, 8)]
                        #Rule 'cc': [(2, 2), (2, 4), (4, 4), (2, 6), (4, 6), (6, 6)]
    
    '''?????'''
    G = networkx.Graph()
    G.add_edges_from(rule_reactants)
    sp_interactions = [G[node].keys() for node in networkx.nodes(G)]
    '''Printed for each rule'''
    #print sp_interactions 
                        #Rule 'ab': [[1, 10, 4, 9, 7], [0], [0], [0], [0], [0]]
                        #Rule 'bc': [[8, 2, 5, 7], [1, 3], [8, 2, 5, 7], [1, 3], [1, 3], [1, 3]]
                        #Rule 'cc': [[2, 4, 6], [2, 4, 6], [2, 4, 6]]
    vars_set = set(map(tuple, sp_interactions))
    final_vars = map(list, vars_set)
    '''Printed for each rule''' 
    #print final_vars
                    #Rule 'ab': [[1, 10, 4, 9, 7], [0]]
                    #Rule 'bc': [[1, 3], [8, 2, 5, 7]]
                    #Rule 'cc': [[2, 4, 6]]
    if len(final_vars) == 1:
        final_vars.append(final_vars[0])
    interactions = list(itertools.product(*final_vars))
    '''Printed for each rule'''
    #print interactions 
                        #Rule 'ab': [(1, 0), (10, 0), (4, 0), (9, 0), (7, 0)]
                        #Rule 'bc': [(1, 8), (1, 2), (1, 5), (1, 7), (3, 8), (3, 2), (3, 5), (3, 7)]
                        #Rule 'cc': [(2, 2), (2, 4), (2, 6), (4, 2), (4, 4), (4, 6), (6, 2), (6, 4), (6, 6)]
    '''From this, we have the list of reactions that involved in each rule'''
    

    '''Getting the reaction rate of interactions'''
    for interaction in interactions:
        for reaction in model.reactions:
            #for (k1*X2*X3+...)
            if reaction['rule'][0] == rule.name and reaction['reverse'][0] == False and sorted(
                    reaction['reactants']) == sorted(interaction):
                #for symmetry reactions
                if reaction['reactants'][0] == reaction['reactants'][1]:
                    rcts_no_reverse += 2 * reaction['rate']
                #for non symmetry reaction
                else:
                    rcts_no_reverse += reaction['rate']
            #for reverse term (l1*X1+l2*X2+...)
            if reaction['rule'][0] == rule.name and reaction['reverse'][0] == True and sorted(
                    reaction['products']) == sorted(interaction):
                #for symmetry reactions
                if reaction['products'][0] == reaction['products'][1]:
                    rcts_reverse += 2 * reaction['rate']
                #for non symmetry reaction
                else:
                    rcts_reverse += reaction['rate']
    '''ODE for each rule'''
    rcts_rules[rule.name] = {'no_reverse': sympy.factor(rcts_no_reverse), 'reverse': sympy.factor(rcts_reverse)}
    '''Printed for each rule'''
    #print rcts_rules 
                    #1st loop:
                    #{'ab': {'no_reverse': __s0*kab*(__s1 + __s10 + __s4 + __s7 + 2*__s9), 'reverse': 0}}
                    
                    #2nd loop:
                    #{'ab': {'no_reverse': __s0*kab*(__s1 + __s10 + __s4 + __s7 + 2*__s9), 'reverse': 0}, 
                    #'bc': {'no_reverse': kbc*(__s1 + __s3)*(__s2 + 2*__s5 + __s7 + __s8), 'reverse': lbc*(2*__s10 + 2*__s11 + __s4 + __s6 + __s7 + __s8 + 2*__s9)}}
                    
                    #3rd loop:
                    #{'cc': {'no_reverse': 1.0*kcc*(1.0*__s2 + 1.0*__s4 + 1.0*__s6)**2, 'reverse': 2*lcc*(__s10 + __s11 + __s5 + __s7 + __s8 + __s9)}, 
                    #'ab': {'no_reverse': __s0*kab*(__s1 + __s10 + __s4 + __s7 + 2*__s9), 'reverse': 0}, 
                    #'bc': {'no_reverse': kbc*(__s1 + __s3)*(__s2 + 2*__s5 + __s7 + __s8), 'reverse': lbc*(2*__s10 + 2*__s11 + __s4 + __s6 + __s7 + __s8 + 2*__s9)}}
'''From this, we have the group ODEs'''
    
    
# this line loops over the reactions and factor the monomials and that way  we get the new units
new_units = {}
for rc in rcts_rules:
    var = sympy.simplify(sympy.factor(rcts_rules[rc]['no_reverse'])).as_coeff_mul()[1]
    '''To separate Var on each term and storage them as element of tuple'''
    #print var
            #'cc': 
            #(1.00000000000000, kcc, (__s2 + __s4 + __s6)**2)
            #'ab': 
            #(__s0, kab, __s1 + __s10 + __s4 + __s7 + 2*__s9)
            #'bc': 
            #(kbc, __s1 + __s3, __s2 + 2*__s5 + __s7 + __s8)
    final_vars = [ex for ex in var if (str(ex).startswith('__') or type(ex) == sympy.Pow or type(ex) == sympy.Add)]
    '''Get and remember the Var only with s0,s1,...'''
    #print final_vars
            #'cc':
            #[(__s2 + __s4 + __s6)**2]
            #'ab': 
            #[__s0, __s1 + __s10 + __s4 + __s7 + 2*__s9]
            #'bc': 
            #[__s1 + __s3, __s2 + 2*__s5 + __s7 + __s8]
    new_units[rc] = final_vars
    '''Storage the var to new_units Dict for each rcts_rules'''
    #print new_units
                    #1st loop:
                    #{'cc': [(__s2 + __s4 + __s6)**2]}
                    
                    #2nd loop:
                    #{'cc': [(__s2 + __s4 + __s6)**2], 
                    #'ab': [__s0, __s1 + __s10 + __s4 + __s7 + 2*__s9]}
                    
                    #3rd loop:
                    #{'cc': [(__s2 + __s4 + __s6)**2], 
                    #'ab': [__s0, __s1 + __s10 + __s4 + __s7 + 2*__s9], 
                    #'bc': [__s1 + __s3, __s2 + 2*__s5 + __s7 + __s8]}
'''From this we have list of term that will be defined as new variable'''


'''To write Group ODEs'''
# this defines the odes for the new units
new_units_odes = {}
for r in new_units:
    for eq in new_units[r]:
        if type(eq) == sympy.Pow:
            '''get the inside of pow'''
            eq_name = eq.args[0]
            print eq.args[0] 
                            #__s2 + __s4 + __s6
            new_units_odes[eq_name] = sympy.simplify(rcts_rules[r]['reverse'] - rcts_rules[r]['no_reverse'])
            print new_units_odes
                                #{__s2 + __s4 + __s6: -1.0*kcc*(__s2 + __s4 + __s6)**2 + 2*lcc*(__s10 + __s11 + __s5 + __s7 + __s8 + __s9)}
        else:
            new_units_odes[eq] = rcts_rules[r]['reverse'] - rcts_rules[r]['no_reverse']
            print new_units_odes
                                #{__s2 + __s4 + __s6: -1.0*kcc*(__s2 + __s4 + __s6)**2 + 2*lcc*(__s10 + __s11 + __s5 + __s7 + __s8 + __s9), 
                                #__s0: -__s0*kab*(__s1 + __s10 + __s4 + __s7 + 2*__s9)}

                                #{__s2 + __s4 + __s6: -1.0*kcc*(__s2 + __s4 + __s6)**2 + 2*lcc*(__s10 + __s11 + __s5 + __s7 + __s8 + __s9), 
                                #__s1 + __s10 + __s4 + __s7 + 2*__s9: -__s0*kab*(__s1 + __s10 + __s4 + __s7 + 2*__s9), 
                                #__s0: -__s0*kab*(__s1 + __s10 + __s4 + __s7 + 2*__s9)}

                                #{__s2 + __s4 + __s6: -1.0*kcc*(__s2 + __s4 + __s6)**2 + 2*lcc*(__s10 + __s11 + __s5 + __s7 + __s8 + __s9), 
                                #__s1 + __s10 + __s4 + __s7 + 2*__s9: -__s0*kab*(__s1 + __s10 + __s4 + __s7 + 2*__s9), 
                                #__s0: -__s0*kab*(__s1 + __s10 + __s4 + __s7 + 2*__s9), 
                                #__s1 + __s3: -kbc*(__s1 + __s3)*(__s2 + 2*__s5 + __s7 + __s8) + lbc*(2*__s10 + 2*__s11 + __s4 + __s6 + __s7 + __s8 + 2*__s9)}

                                #{__s2 + __s4 + __s6: -1.0*kcc*(__s2 + __s4 + __s6)**2 + 2*lcc*(__s10 + __s11 + __s5 + __s7 + __s8 + __s9), 
                                #__s1 + __s10 + __s4 + __s7 + 2*__s9: -__s0*kab*(__s1 + __s10 + __s4 + __s7 + 2*__s9), 
                                #__s2 + 2*__s5 + __s7 + __s8: -kbc*(__s1 + __s3)*(__s2 + 2*__s5 + __s7 + __s8) + lbc*(2*__s10 + 2*__s11 + __s4 + __s6 + __s7 + __s8 + 2*__s9), 
                                #__s0: -__s0*kab*(__s1 + __s10 + __s4 + __s7 + 2*__s9), 
                                #__s1 + __s3: -kbc*(__s1 + __s3)*(__s2 + 2*__s5 + __s7 + __s8) + lbc*(2*__s10 + 2*__s11 + __s4 + __s6 + __s7 + __s8 + 2*__s9)}
'''From this, we have group of ODEs'''

'''TO BE CONTINUED'''
'''this adds more conservation laws from the reactions'''
for idx, nam in enumerate(new_units):
    if len(new_units[nam]) > 1:
        new_cons = new_units[nam][0] - new_units[nam][1] - sympy.symbols('b%d' % idx, real=True)
        new_cons_value = conservation_laws_values(model, new_cons)
        cl.append(new_cons)
        cl_values[new_cons_value.keys()[0]] = new_cons_value.values()[0]


# this uses the conservation laws to define the new units in terms of the other variables
equal_units = {}
for rul in new_units:
    for un in new_units[rul]:
        for cc in cl:
            if sympy.solve(cc, un):
                equal_units[un] = sympy.solve(cc, un)[0]

# this line defines the reverse reactions in term of the other variables (to have the same variable in the eq)
equal_const = {}
for r in rcts_rules:
    if not rcts_rules[r]['reverse'] == 0:
        for cc in cl:
            # print sympy.solve(sympy.collect_const(cc.rhs), rcts_rules[r]['reverse'].as_two_terms()[1])
            if sympy.solve(cc, rcts_rules[r]['reverse'].as_coeff_mul()[1][1]):
                equal_const[rcts_rules[r]['reverse'].as_coeff_mul()[1][1]] = sympy.collect_const(
                        sympy.solve(cc, rcts_rules[r]['reverse'].as_coeff_mul()[1][1])[0])
            if sympy.solve(sympy.collect_const(cc), rcts_rules[r]['reverse'].as_coeff_mul()[1][1]):
                equal_const[rcts_rules[r]['reverse'].as_coeff_mul()[1][1]] = sympy.collect_const(
                        sympy.solve(sympy.collect_const(cc), rcts_rules[r]['reverse'].as_coeff_mul()[1][1])[0].evalf())

equal_const_units = equal_units.copy()
equal_const_units.update(equal_const)

# this is a dictionary where the keys are the new units and the values are dictionaries that contain all the
# information to define that equation in terms of una variable.
variables_to_change = {}
for unit in new_units:
    if len(new_units[unit]) > 1:
        for idx in range(len(new_units[unit])):
            if rcts_rules[unit]['reverse'] != 0:
                tmp_dict = {new_units[unit][idx]: equal_const_units[new_units[unit][idx]],
                            rcts_rules[unit]['reverse'].as_coeff_mul()[1][1]: equal_const_units[
                                rcts_rules[unit]['reverse'].as_coeff_mul()[1][1]]}
            else:
                tmp_dict = {new_units[unit][idx]: equal_const_units[new_units[unit][idx]]}
            variables_to_change[new_units[unit][1 - idx]] = tmp_dict
    else:
        tmp_dict = {rcts_rules[unit]['reverse'].as_coeff_mul()[1][1]: equal_const_units[
            rcts_rules[unit]['reverse'].as_coeff_mul()[1][1]]}
        variables_to_change[new_units[unit][0].as_base_exp()[0]] = tmp_dict

# this is a simple change of variable to allow sympy to solve the differential equations.
new_ode_vars = {}
final_odes = {}
for num, ode in enumerate(new_units_odes):
    new_ode_vars[ode] = sympy.symbols('U%d' % num)
    final_odes[sympy.symbols('U%d' % num)] = new_units_odes[ode].subs(variables_to_change[ode]).subs(ode, sympy.symbols(
            'U%d' % num))

new_ode_vars_ic = {}
for idx, u in enumerate(new_ode_vars):
    new_ode_vars_ic[new_ode_vars[u](0)] = conservation_laws_values(model, u - sympy.Symbol('d%d' % idx)).values()[0]

# solving the differential equations
solutions = []
equations = []
solutions_eq = []
for nom in final_odes:
    s = sympy.var('s')
    t = sympy.symbols('t')
    equation = sympy.Eq(final_odes[nom].subs(nom, nom(t)), nom(t).diff(t))
    equations.append(equation)
    sol = sympy.dsolve(sympy.expand(equation), nom(t), simplify=False)
    sol = sympy.simplify(sol)
    sol_lhs = sol.lhs
    sol_factor = sol_lhs.as_coeff_mul()[1][0]
    sol_step1 = sol.do(s * (1 / sol_factor))
    sol_step2 = sol_step1.do(sympy.exp(s))
    ode_par = sol.subs(t, 0).subs(new_ode_vars_ic)
    explicit_sol = sympy.solve(sol_step2, nom(t), dict=True)
    for s in explicit_sol:
        
        s[s.keys()[0]] = s.values()[0].subs({ode_par.rhs: ode_par.lhs})
    solutions.append(explicit_sol)

for s in solutions:
    for q in s:
        solutions_eq.append(sympy.Eq(q.keys()[0], q.values()[0]))

eqs_to_evaluate = []
for idx, solu in enumerate(solutions_eq):
    sol_ic = solu.subs(cl_values)
    for p in model.parameters: sol_ic = sol_ic.subs(p.name, p.value)
    sol_ic_copy = sol_ic
    print sol_ic
    if sympy.simplify(sol_ic.subs(sympy.Symbol('t'), 0)).rhs > 0:
        eqs_to_evaluate.append(sol_ic_copy)


import matplotlib.pyplot as plt
import numpy
from pysb.integrate import odesolve


t = numpy.linspace(0,30000,10000)
t2=numpy.linspace(0,5,500)
x = odesolve(model,t2, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)

f = sympy.lambdify(sympy.Symbol('t'), eqs_to_evaluate[2].rhs, 'numpy')

plt.plot(t2,x['s1pluss3'], linewidth=3,label='numerical')
plt.plot(t2, f(t2), 'r--', linewidth=5, label='theoretical')
plt.legend(loc=0)
plt.show()