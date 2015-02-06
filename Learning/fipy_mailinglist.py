from fipy import *
import numpy
m = Grid1D(nx=100, Lx=1.)
#m = Grid2D(nx=100, ny=100, Lx=1., Ly=1.0)
v = CellVariable(mesh=m, hasOld=True, value=[[0.5], [1], [2]], elementshape=(3,))

v.constrain([[0], [1], [0.5]], where=m.facesLeft)
v.constrain([[1], [0], [0]], where=m.facesRight)

M = [[1,0,0],[0,1,0],[0,0,1]]
# M = numpy.identity(3)
N = [[1,0,0],[0,1,0],[0,0,1]]
U = [-v[0] * v[1] + v[2], -v[0] * v[1] + v[2], v[0] * v[1] - v[2]]
Q = U[0] * [[1], [0], [0]] + U[1] * [[0], [1], [0]] + U[2] * [[0], [0], [1]] #(-v[0]) * (1,1)#(-v[0] * v[1], v[0]**2) 
print v[0]
print Q 
eqn = TransientTerm(M) == DiffusionTerm([M]) + (Q)

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

for t in range(100): 
     v.updateOld()
     eqn.solve(var=v, dt=1.e-3)
     if __name__ == '__main__':
         vi.plot()
if __name__ == '__main__':
     raw_input("Press <return> to proceed...")