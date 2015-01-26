from pysb.integrate import odesolve
import numpy as np
from fundamental_aa_reaction_model import model

t=np.linspace(0,100)

zout = odesolve(model,t)
print model.odes