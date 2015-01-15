from pysb import *
from pysb.integrate import odesolve
import pylab as pl

Model()

Monomer('L',['s'])
Monomer('R',['s'])

Parameter('kf',1e-3)
Parameter('kr',1e-3)

Initial(L(s=None),Parameter('L_0',100))
Initial(R(s=None),Parameter('R_0',200))

Rule('L_binds_R',L(s=None)+R(s=None) <> L(s=1)%R(s=1),kf,kr)

Observable('LR',L(s=1)%R(s=1))
Observable('tL',L(s=None))
Observable('tR',R(s=None))

time=pl.linspace(0,40)
x=odesolve(model,time)
pl.plot(time,x['LR'])

pl.ion()
pl.figure()
pl.plot(t, x['LR'], label="LR")
pl.plot(t, x['tL'], label="L")
pl.plot(t, x['tR'], label="R")
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("molecules/cell")
pl.show()