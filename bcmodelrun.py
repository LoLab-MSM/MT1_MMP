#import pylab as pl
from pysb import *
from MT1_MMP import bcmodel
from MT1_MMP.bcmodel import model

import numpy as np
from pysb.integrate import odesolve
#declare monomers
bcmodel.monomer_bc_model()

#declare rate constant
bcmodel.rate_constant_bc_model()

#the rules
bcmodel.rule_original_bc_model()

#initial conditions
bcmodel.initial_condition_bc_model()

bcmodel.observe_bc_model()

t=np.linspace(0,5,100)

zout = odesolve(model, t)
#zout = odesolve(model,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
print model.odes

