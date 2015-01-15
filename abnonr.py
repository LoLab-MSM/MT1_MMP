from pysb import *
from pysb.integrate import odesolve
import pylab as pl
import numpy as np

Model()
Monomer('A',['b'])
Monomer('B',['b'])

Parameter('k', 1e-2)

Initial(A(b=None), Parameter('Ao', 1000))
Initial(B(b=None), Parameter('Bo', 1500))

Rule('C', A(b=None) + B(b=None) >> A(b=1) % B(b=1), k)

Observable('tA', A(b=None))
Observable('tB', B(b=None))
Observable('tC', A(b=1)%B(b=1))

t= np.linspace(0,1)

sol = odesolve(model,t)

pl.ion()
pl.figure()
pl.plot(t, sol['tA'], label="mol A")
pl.plot(t, sol['tB'], label="mol B")
pl.plot(t, sol['tC'], label="mol C")
pl.legend()
pl.xlabel("Time")
pl.ylable("Mol")
pl.show()