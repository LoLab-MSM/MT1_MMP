from pysb.integrate import odesolve
import numpy as np
import pylab as pl
from c8bid import model

t=np.linspace(0,20000)

yout = odesolve(model,t,verbose=True)
pl.figure()
pl.plot(t, yout['uBid'], label="Bid")
pl.plot(t, yout['obspBid'], label="pBid")
pl.plot(t, yout['obsC8'], label="C8")
pl.plot(t, yout['C8Bid'], label="c8bid")
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("molecules/cell")
pl.show()

