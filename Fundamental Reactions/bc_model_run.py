#import pylab as pl
import bc_model
from bc_model import model

import numpy as np
from pysb.integrate import odesolve
#declare monomers
bc_model.monomer_bc_model()

#declare rate constant
bc_model.rate_constant_bc_model()

#the rules
bc_model.rule_original_bc_model()

#initial conditions
bc_model.initial_condition_bc_model()

bc_model.observe_bc_model()
print model.monomers.keys()
quit()

t=np.linspace(0,5,100)

zout = odesolve(model, t)
#zout = odesolve(model,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
print model.odes