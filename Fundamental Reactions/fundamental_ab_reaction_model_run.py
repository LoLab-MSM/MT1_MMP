from pysb.integrate import odesolve#, rdesolve
from rdesolver import rdesolve
#from notes_rdesolver import rdesolve
import numpy as np
import pylab as pl
from fundamental_ab_reaction_model import model


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
initial = [1, 0.5, 2]
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