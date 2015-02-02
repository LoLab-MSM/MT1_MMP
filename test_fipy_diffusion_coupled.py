from fipy import *
m = Grid1D(nx=100, Lx=1.)

#v0 = CellVariable(mesh=m, hasOld=True, value=0.5)
#v1 = CellVariable(mesh=m, hasOld=True, value=0.5)

#v0.constrain(0, m.facesLeft)
#v0.constrain(0, m.facesRight)
#v1.constrain(0, m.facesLeft)
#v1.constrain(0, m.facesRight)

#vi = Viewer((v0, v1))

#eqn0 = TransientTerm(var=v0) == DiffusionTerm(0.01, var=v0) - DiffusionTerm(1, var=v1)
#eqn1 = TransientTerm(var=v1) == DiffusionTerm(1, var=v0) + DiffusionTerm(0.01, var=v1)

#eqn0 = TransientTerm(var=v0) == DiffusionTerm(1, var=v1)
#eqn1 = TransientTerm(var=v1) == DiffusionTerm(1, var=v0) 

#eqn0 = TransientTerm(var=v0) == DiffusionTerm(1, var=v1) - v0*v1
#eqn1 = TransientTerm(var=v1) == DiffusionTerm(1, var=v0) - v0*v1

#eqn = eqn0 & eqn1

#for t in range(1): 
#     v0.updateOld()
#     v1.updateOld()
 #    eqn.solve(dt=1.e-3)
 #    vi.plot()
     
#solution as vector
#problem : how to put source term?
v = CellVariable(mesh=m, hasOld=True, value=[[0.5],[0.5]], elementshape=(2,))
v.constrain([[0], [0]], m.facesLeft)
v.constrain([[0], [0]], m.facesRight)
vi = Viewer(vars=(v[0],v[1]))
S0 = -v[0]*v[1]
S1 = -v[0]*v[1]
v0 = v[0]
v1 = v[1]
source = -v0*v1 *[1,1]
eqn = TransientTerm([[1,0],[0,1]]) == DiffusionTerm([[[1,0],[0,1]]]) #+ source

for t in range(1): 
     v.updateOld()
     eqn.solve(var=v, dt=1.e-3)
     vi.plot()