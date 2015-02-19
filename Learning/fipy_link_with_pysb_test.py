from pysb import *
import pysb.bng
from fipy.meshes.factoryMeshes import Grid1D
from scipy.integrate._ode import ode

#PYSB
Model()

Monomer('v0', ['a'])
Monomer('v1', ['a'])

Parameter('k', 0.75)
Parameter('l', 0.5)

Parameter('v0_init', 100)
Parameter('v1_init', 200)
Parameter('v2_init', 50)

Initial(v0(a=None), v0_init)
Initial(v1(a=None), v1_init)
Initial(v0(a=1)%v1(a=1), v2_init)

Rule('v2', v0(a=None) + v1(a=None) <> v0(a=1)%v1(a=1), k, l)

##############################
#species, reactions, reactions_bidirectional, observables
pysb.bng.generate_equations(model)
##### DIFFUSIVITIES
model.diffusivities = [1.5]*len(model.species)
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
m = fipy.Grid1D(nx=100, Lx=1.)
mm = fipy.Grid2D(nx=100, ny=100, Lx=1., Ly=1.)

m3 = fipy.Gmsh2DIn3DSpace('''
     radius = 5.0;
     cellSize = 0.3;
     
     // create inner 1/8 shell
     Point(1) = {0, 0, 0, cellSize};
     Point(2) = {-radius, 0, 0, cellSize};
     Point(3) = {0, radius, 0, cellSize};
     Point(4) = {0, 0, radius, cellSize};
     Circle(1) = {2, 1, 3};
     Circle(2) = {4, 1, 2};
     Circle(3) = {4, 1, 3};
     Line Loop(1) = {1, -3, 2} ;
     Ruled Surface(1) = {1};
     
     // create remaining 7/8 inner shells
     t1[] = Rotate {{0,0,1},{0,0,0},Pi/2} {Duplicata{Surface{1};}};
     t2[] = Rotate {{0,0,1},{0,0,0},Pi} {Duplicata{Surface{1};}};
     t3[] = Rotate {{0,0,1},{0,0,0},Pi*3/2} {Duplicata{Surface{1};}};
     t4[] = Rotate {{0,1,0},{0,0,0},-Pi/2} {Duplicata{Surface{1};}};
     t5[] = Rotate {{0,0,1},{0,0,0},Pi/2} {Duplicata{Surface{t4[0]};}};
     t6[] = Rotate {{0,0,1},{0,0,0},Pi} {Duplicata{Surface{t4[0]};}};
     t7[] = Rotate {{0,0,1},{0,0,0},Pi*3/2} {Duplicata{Surface{t4[0]};}};
     
     // create entire inner and outer shell
     Surface Loop(100)={1,t1[0],t2[0],t3[0],t7[0],t4[0],t5[0],t6[0]};
''', order=2).extrude(extrudeFunc=lambda r: 1.1 * r)

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
noise = fipy.GaussianNoiseVariable(mesh=m3,mean=0.5,variance=0.01).value
v=fipy.CellVariable(mesh=m3, hasOld=True,elementshape=(len(model.species),)) # value=ic,
v[:]=noise
"""Define Fixed-Boundary Conditions"""
#?????
# v.constrain([[0], [0], [0]], where=m.facesLeft)
# v.constrain([[0], [0], [0]], where=m.facesRight)

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
time=numpy.linspace(0,10,100)
#########
for t in range(10):
     v.updateOld()
     eqn.solve(var=v, dt=1.e-3)
     vi.plot(filename = 'img%05d.png' % t)
video_file_name = 'video.wmv'
import os
if os.path.isfile(video_file_name):
    print "removing old video file. Check to see if this is correct"
    os.remove(video_file_name)
os.system('avcon -r 2 -i img\%05d.png %s' % video_file_name)


if __name__ == '__main__':
     raw_input("Press <return> to proceed...")