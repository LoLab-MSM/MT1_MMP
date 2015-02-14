from pysb import *
import pysb.bng
from fipy.meshes.factoryMeshes import Grid1D
from scipy.integrate._ode import ode


#PYSB
Model()

Monomer('v0', ['a'])
Monomer('v1', ['a'])

Parameter('k', 1.)
Parameter('l', 1.)

Parameter('v0_init', 0.5)
Parameter('v1_init', 0.3)
Parameter('v2_init', 0.1)

Initial(v0(a=None), v0_init)
Initial(v1(a=None), v1_init)
Initial(v0(a=1)%v1(a=1), v2_init)

Rule('v2', v0(a=None) + v1(a=None) <> v0(a=1)%v1(a=1), k, l)

##############################
#species, reactions, reactions_bidirectional, observables
pysb.bng.generate_equations(model)
##### DIFFUSIVITIES
model.diffusivities = [0.1]*len(model.species)
#####
print
print '##################################~PYSB~##################################'
print
print 'Species:'
for spc in model.species:
    print spc
print
print 'Reactions:' 
for r in model.reactions:
    print r
print
print 'ODES:'
for ode in model.odes:
    print ode
print
print 'Initial Conditions'
for init in model.initial_conditions:
    print init
print
for d in model.diffusivities:
    print d
print
print '##################################~FIPY~##################################'
print
##############################
#FIPY 
import fipy
import numpy
import re

"""Create Mesh"""
#??? in 2d, 3d, and Sphere Coordinate?
m = Grid1D(nx=100, Lx=1.)

"""Call Initial Conditions"""
# a=[]
# for sp in model.initial_conditions:
#     a.append(sp[0])
# print a

##call index initial value
# print model.get_species_index(model.initial_conditions[1][0])
index_nonzero_init = [model.get_species_index(i[0]) for i in model.initial_conditions]
# print index_nonzero_init

##input initial value into list
initt = numpy.zeros(len(model.species))
initt[index_nonzero_init] = [i[1].value for i in model.initial_conditions]
# print initt

##reshape initt
# a=numpy.array([3,2,1])
# b = numpy.reshape(a,(3,1))
# print b
ic = numpy.reshape(initt, (len(initt),1))
# print ic
# print

"""Define CellVariables"""
v=fipy.CellVariable(mesh=m, hasOld=True, value=ic, elementshape=(len(model.species),))

"""Define Fixed-Boundary Conditions"""
#?????
v.constrain([[1], [1], [1]], where=m.facesLeft)
v.constrain([[0], [0], [0]], where=m.facesRight)

"""Define Source Matrix"""
#ODES will be a list of odes(type: fipy's variable)
ODES = []
for ode in model.odes:
    ode=str(ode)
    ##modify SPECIES
    for i in range(len(model.odes)):
        ode = re.sub('__s%d' % i, 'v[%d]' % i, ode)
    ##modify RATE CONSTANT
    ##call the name and values and store them to lists
    rc_name = [model.rules[rxn['rule']].rate_reverse.name if rxn['reverse'] else model.rules[rxn['rule']].rate_forward.name for rxn in model.reactions]
    r = [model.rules[rxn['rule']].rate_reverse.value if rxn['reverse'] else model.rules[rxn['rule']].rate_forward.value for rxn in model.reactions]
    ##modifying
    for i in range(len(rc_name)):
        ode = re.sub(rc_name[i],'r[%i]' % i,ode)
    ##create ydot=eqn
    ode = 'vdot = ' + ode 
    ##calculate the eqn
    ode_cal = compile(ode, '<$$$>', 'exec')
    exec ode_cal in locals()
    ODES.append(vdot)
print ODES
print type(ODES)

"""%%%How to call rate constant name and value and input them to lists%%%"""
# w=model.odes[0]
# w=str(w)
# print w
# print type(w)
# print model.rules['v2']
# print model.rules[0].rate_forward.name 
# print model.rules[0].rate_forward.value
# print list(enumerate(model.parameters_rules()))


##call name of rate constants
# q=[]
# for rxn in model.reactions:#3
#     if rxn['reverse']:
#         q.append(model.rules[rxn['rule']].rate_reverse.name)#1
#     else:
#         q.append(model.rules[rxn['rule']].rate_forward.name)#2
# print q

#model.rules[rxn['rule']].rate_reverse.name if rxn['reverse'] else model.rules[rxn['rule']].rate_forward.name for rxn in model.reactions

# q_name = [model.rules[rxn['rule']].rate_reverse.name if rxn['reverse'] else model.rules[rxn['rule']].rate_forward.name for rxn in model.reactions]
# print q
"""%%%How to call rate constant name and value and input them to lists%%%"""

"""Equations"""
M = numpy.identity(len(model.species))
#####
for i,d in enumerate(model.diffusivities):
    M[i,i] = d
#####

N = M
N = numpy.reshape(M,(len(model.species),len(model.species),1))
Q=int()
for i in range(len(ODES)):
    Q += ODES[i] * N[i]
print Q
print type(Q)
eqn = fipy.TransientTerm(M) == fipy.DiffusionTerm([M]) + (Q)

print eqn

"""Perform Integration"""
s=[]
for j in range(len(model.species)):
    s.append(v[j])
vi = fipy.Viewer(vars=s)
#time????
for t in range(10): 
     v.updateOld()
     eqn.solve(var=v, dt=1.e-3)
     if __name__ == '__main__':
         vi.plot()
if __name__ == '__main__':
     raw_input("Press <return> to proceed...")