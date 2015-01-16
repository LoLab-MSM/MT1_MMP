from pysb.integrate import odesolve
import numpy as np
import pylab as pl
from MT1-MMP import mt1mmp

monomer_abc_model()
    
rate_constant_abc_model()

rule_original_abc_model()

initial_condition_abc_model()

observe_abc_model()

t=np.linspace(0,5,400)

zout = odesolve(model,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
#print model.odes
#print zout

#pl.ion()
#pl.figure()
#pl.title('Important Component in MT1-MMP Model')
#pl.plot(t, zout['abcc'], label="cancer substrate")
#pl.plot(t, zout['ta'], label="mmp2")
#pl.plot(t, zout['tb'], label="timp2")
#pl.plot(t, zout['tc'], label="mt1mmp")
#pl.legend()
#pl.xlabel("Time (s)")
#pl.ylabel("Molecules or Cells")
#pl.show()

#pl.ion()
#pl.figure()
#
#pl.title('Time Course Simulation')
#pl.plot(t, zout['ta'], label="a")
#pl.plot(t, zout['tb'], label="b")
#pl.plot(t, zout['tc'], label="c")
#pl.plot(t, zout['tab'], label="ab")
#pl.plot(t, zout['tbc'], label="bc")
#pl.plot(t, zout['tcc'], label="cc")
#pl.plot(t, zout['abc'], label="abc")
#pl.plot(t, zout['bcc'], label="bcc")
#pl.plot(t, zout['abcc'], label="abcc")
#pl.plot(t, zout['abccb'], label="abccb")
#pl.plot(t, zout['abccba'], label="abccba")
#pl.legend()
#pl.xlabel("Time (s)")
#pl.ylabel("Molecules")
#pl.show()
#
#pl.ion()
#pl.figure()
#pl.xlim(0,.4)
#pl.title('Type 1 (disappeared soon)')
#pl.plot(t, zout['tb'], label="b")
#pl.plot(t, zout['tbc'], label="bc")
#pl.plot(t, zout['bcc'], label="bcc")
#pl.plot(t, zout['bccb'], label="bccb")
#pl.plot(t, zout['abccb'], label="abccb")
#pl.legend()
#pl.xlabel("Time (s)")
#pl.ylabel("Molecules")
#pl.show()
#
#pl.ion()
#pl.figure()
#pl.title('Type 2 (disappeared after long time)')
#pl.plot(t, zout['tab'], label="ab")
#pl.plot(t, zout['tc'], label="c")
#pl.legend()
#pl.xlabel("Time (s)")
#pl.ylabel("Molecules")
#pl.show()
#
#pl.ion()
#pl.figure()
#pl.xlim(0,.4)
#pl.title('Type 3 (Remain at some values)')
#pl.plot(t, zout['ta'], label="a")
#pl.plot(t, zout['tcc'], label="cc")
#pl.plot(t, zout['abc'], label="abc")
#pl.plot(t, zout['abcc'], label="abcc")
#pl.plot(t, zout['abccba'], label="abccba")
#pl.legend()
#pl.xlabel("Time (s)")
#pl.ylabel("Molecules")
#pl.show()


initb = np.linspace(1e-9, 1e-5, 1000)
N=[]

for i in range (len(initb)):
    model.parameters['bo'].value = initb[i]
    zout = odesolve(model,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
    N.append(zout['abcc'][-1])

pl.ion()    
pl.figure()
#pl.plot(initb, N, 'o')
pl.semilogx(initb, N)
pl.show()
    