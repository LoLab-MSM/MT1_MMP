from pysb.integrate import odesolve, rdesolve
import numpy as np
import pylab as pl
from fundamental_ab_reaction_model import model


t=np.linspace(0,40)

zout = odesolve(model,t)

print model.odes
print model.species

##Solve PDE using fipy
pout = rdesolve(model)

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
