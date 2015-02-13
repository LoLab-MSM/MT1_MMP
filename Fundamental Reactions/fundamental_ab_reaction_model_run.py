from pysb.integrate import odesolve#, rdesolve
from rdesolver import rdesolve
#from rdesolver_test import rdesolve
import numpy as np
import pylab as pl
from fundamental_ab_reaction_model import model

##########
from pysb.bng import generate_equations
 
generate_equations(model, verbose=True)
# 
# for ic in model.initial_conditions:
#     print ic
# print
# for sp in model.species:
#     print sp
# print
# 
# x = np.zeros(len(model.species))
# print x
# index = [model.get_species_index(ic[0]) for ic in model.initial_conditions]
# x[index] = [ic[1].value for ic in model.initial_conditions]
# print x

# for rxn in model.reactions:
#     print rxn['rate']
#     for arg in rxn['rate'].args:
#         print '\t', arg, type(arg)

# print model.reactions[0]['rate']
# for arg in model.reactions[0]['rate'].args:
#     print arg
  
# print model.odes[0], type(model.odes[0])

# import re
# x = [10., 20., 30.]
# p = [1., 2.]
# p_names = ['k', 'l']
# 
# print model.odes[0]
# ode = str(model.odes[0])
# 
# # species
# for i in range(len(x)):
#     ode = re.sub(r'__s%d' % i, 'x[%d]' % i, ode)
# print ode
# 
# # rate constants
# for i in range(len(p)):
#     ode = re.sub(p_names[i], 'p[%d]' % i, ode)
# print ode
# 
# ode = 'xdot = ' + ode
# print ode
# 
# ode_py = compile(ode, '<test>', 'exec')
# print ode_py
# exec ode_py in locals()
# print xdot

# v[0] * v[1] * k[0]

for p in model.parameters_rules():
    print p.name, p.value
print

# k = []
# for rxn in model.reactions:
#     print rxn
#     rule = model.rules[rxn['rule']]
#     print rule
#     if rxn['reverse']:
#         val = rule.rate_reverse.value
#     else:
#         val = rule.rate_forward.value
#     print val
#     k.append(val)
# print k


# k = [model.rules[rxn['rule']].rate_reverse.value if rxn['reverse'] else model.rules[rxn['rule']].rate_forward.value
#      for rxn in model.reactions]
# 
# print k

# p = [k for k in model.reactions]

quit()
##############

t=np.linspace(0,40)

zout = odesolve(model,t)

for i in range(len(model.odes)):
     print model.species[i], ": __s", i, ":", model.odes[i]


##Solve PDE using fipy
#model, mesh, initc, Dirichlet = False, Neumann = False, Diffusivity
#mesh='1d' or '2d' or '3d' , Lx, nx, Ly, ny
#initial = [0.5, 0.5, 0]
#boundary condition (Dirichlet / Neumann = True or False) True means the Dirichlet boundary condition is not zero
#D=[1,1,0]

"""parameters"""
mesh = ['1d', 100, 1.]
initial = [0.5, 0.3, 0.1]
boundary = [1, 0]
pout = rdesolve(model, mesh, initial, boundary) #, initial

#print zout
#pl.ion()
#pl.figure()
#pl.plot(t, zout['dt'], label="product")
#pl.plot(t, zout['at'], label="mol a")
#pl.plot(t, zout['ct'], label="mol c")
#pl.legend()
#pl.xlabel("Time (s)")
#pl.ylabel("Molecules or Cells")
#pl.show()

if __name__ == '__main__':
         raw_input("Press <return> to proceed...")