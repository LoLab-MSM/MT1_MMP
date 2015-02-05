from fipy import *
nx = ny = 100
#m = Grid2D(dx=0.25, dy=0.25, nx=nx, ny=ny)
m = Grid2D(nx=100, ny=100, Lx=1., Ly=1.) 
# v0 = CellVariable(mesh=m, hasOld=True, value=0.5)
# v1 = CellVariable(mesh=m, hasOld=True, value=0.5)
# 
# v0.constrain(0, m.facesLeft)
# v0.constrain(1, m.facesRight)
# v1.constrain(1, m.facesLeft)
# v1.constrain(0, m.facesRight)
# 
# vi = Viewer((v0, v1))
# 
# #eqn0 = TransientTerm(var=v0) == DiffusionTerm(0.01, var=v0) - DiffusionTerm(1, var=v1)
# #eqn1 = TransientTerm(var=v1) == DiffusionTerm(1, var=v0) + DiffusionTerm(0.01, var=v1)
# 
# eqn0 = TransientTerm(var=v0) == DiffusionTerm(1, var=v0) + 100
# eqn1 = TransientTerm(var=v1) == DiffusionTerm(1, var=v1) 
# 
# #eqn0 = TransientTerm(var=v0) == DiffusionTerm(1, var=v1) - v0*v1
# #eqn1 = TransientTerm(var=v1) == DiffusionTerm(1, var=v0) - v0*v1
# 
# eqn = eqn0 & eqn1
# 
# for t in range(10): 
#      v0.updateOld()
#      v1.updateOld()
#      eqn.solve(dt=1.e-3)
#      vi.plot()
     
#solution as vector
#problem : how to put source term?
v = CellVariable(mesh=m, hasOld=True, value=[[0.5], [0.5]], elementshape=(2,))
v1= v[0].value
print v1
v.constrain([[0], [1]], where=m.facesLeft)
v.constrain([[1], [0]], where=m.facesRight)
# print v.faceValue
S0 = -v[0] * v[1]
S1 = -v[0] * v[1]
source = (S0, S1)
#print source
# eqn = TransientTerm([[1,0],[0,1]]) == DiffusionTerm([[[0.01,-1],[1,0.01]]])
eqn2 = TransientTerm([[1,0],[0,1]]) == DiffusionTerm([[[1,0],[0,1]]])
M = [[0.01,-1],[1,0.01]]
M2 = [[1,0],[0,1]]
N = [[1,0],[0,1]]
P = [100,0]
Q = (v[0], v[1])  #+ [-v[0] * v[1], -v[0] * v[1]]
eqn3 = TransientTerm(N) == DiffusionTerm([M2]) +(Q)
# print eqn3
# quit()

vu = [v[0],v[1]]
if __name__ == '__main__':
#     vi = Viewer(vars=(v[0],v[1]))
     vr=[]
     for i in range(2):
         vr.append(vu[i])   
     vi = Viewer(vars=vr)
for t in range(100): 
     v.updateOld()
     eqn3.solve(var=v, dt=1.e-3)
     if __name__ == '__main__':
         vi.plot()
         

if __name__ == '__main__':
     raw_input("Explicit transient diffusion. Press <return> to proceed...")