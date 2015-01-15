from pysb.integrate import odesolve
import numpy as np
import pylab as pl
from AAtoP import model

t=np.linspace(0,100)

zout = odesolve(model,t)
print model.odes