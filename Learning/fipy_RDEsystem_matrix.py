from fipy import *
import numpy

# v0 + v1 <-> v2  k0,k1

# Monomer('v0')
# Monomer('v1')
# Monomer('v2')

# Parameter('k0', 1)
# Parameter('k1', 1)
# Rule('v0v1_v2', v0() + v1() <> v2(), [k0,k1])

# Parameter('v0_init', 0.5)
# Parameter('v1_init', 0.3)
# Parameter('v2_init', 0.1)
# Initial(v0(), v0_init)
# Initial(v1(), v1_init)
# Initial(v2(), v2_init)

# -----------------------

# pysb.bng.generate_equations(model)

# DEFINE A MESH
# m = ???

# INITIAL CONDITIONS
# initial = numpy.zeros(len(model.species))
# index = [model.get_species_index(ic[0]) for ic in model.initial_conditions]
# initial[index] = [ic[1].value for ic in model.initial_conditions]
# numpy.reshape(initial, (len(initial),1))

# VARIABLES
# v = CellVariable(mesh=m, hasOld=True, value=initial, elementshape=(len(model.species),)

# DEFINE BOUNDARY CONDITIONS
# ???

# DEFINE SOURCE MATRIX
# U = []
# for ode in model.odes:
#    ode = str(ode)
    # SPECIES
#     for i in range(len(model.species)):
#         ode = re.sub(r'__s%d' % i, 'v[%d]' % i, ode)
    # RATE CONSTANTS
#     p = [model.rules[rxn['rule']].rate_reverse.value if rxn['reverse'] else model.rules[rxn['rule']].rate_forward.value
#          for rxn in model.reactions]
#     p_names = [model.rules[rxn['rule']].rate_reverse.name if rxn['reverse'] else model.rules[rxn['rule']].rate_forward.name
#          for rxn in model.reactions]
#     for i,p in enumerate(model.parameters_rules()):
#         ode = re.sub(p_names[i], 'p[%d]' % i, ode)
#     ode = 'vdot = ' + ode
#     ode_py = compile(ode, '<xxx>', 'exec')
#     exec ode_py in locals()
#     U.append[vdot]


m = Grid1D(nx=100, Lx=1.)
#m = Grid2D(nx=100, ny=100, Lx=1., Ly=1.0)
v = CellVariable(mesh=m, hasOld=True, value=[[0.5], [0.3], [0.1]], elementshape=(3,))

v.constrain([[1], [1], [1]], where=m.facesLeft)
v.constrain([[0], [0], [0]], where=m.facesRight)

k = numpy.array([1,1])

# M = [[1,0,0],[0,1,0],[0,0,1]]
nVars = 3
M = numpy.identity(nVars)
N = numpy.identity(nVars)
N = numpy.reshape(M, (len(M),len(M),1))

U = [-v[0] * v[1] * k[0] + v[2] * k[1], -v[0] * v[1] * k[0] + v[2] * k[1], v[0] * v[1] * k[0] - v[2] * k[1]]
# Q = U[0] * [[1], [0], [0]] + U[1] * [[0], [1], [0]] + U[2] * [[0], [0], [1]] #(-v[0]) * (1,1)#(-v[0] * v[1], v[0]**2) 
Q = U[0] * N[0] + U[1] * N[1] + U[2] * N[2]
# Q = U * M
# print Q
# quit()

# print v[0]
# print type(v[0])

eqn = TransientTerm(M) == DiffusionTerm([M]) + (Q)
# eqn2 = TransientTerm([[1.,0.,0.], [0.,1.,0.], [0.,0.,1.]]) == DiffusionTerm([[[1.,0.,0.], [0.,1.,0.], [0.,0.,1.]]])
# print eqn
# print eqn2

# t = U[0] * [[1], [0]]
# w = U[1] * [[0], [1]]
# print v
# print v[0]
# print U[0]
# print t
# # print w
# # print Q
# quit()
vi = Viewer((v[0],v[1], v[2]))

for t in range(10): 
     v.updateOld()
     eqn.solve(var=v, dt=1.e-3)
     if __name__ == '__main__':
         vi.plot()
if __name__ == '__main__':
     raw_input("Press <return> to proceed...")