from pysb.integrate import odesolve
import numpy as np
#import pylab as pl
from MT1-MMP import bcmodel

#declare monomers
monomer_bc_model()

#declare rate constant
rate


t=np.linspace(0,5,100)

zout = odesolve(model, t)
#zout = odesolve(model,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
print model.odes
#print zout