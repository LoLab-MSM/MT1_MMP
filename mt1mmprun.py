from pysb.integrate import odesolve
import numpy as np
import pylab as pl
from MT1_MMP import mt1mmp
from MT1_MMP import mt1mmp as mt1mmp2

model=mt1mmp.return_model('original')
model2=mt1mmp2.return_model('abremoved')


t=np.linspace(0,30000,10000)
t2=np.linspace(0,5,500)
zouty = odesolve(model,t2, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
#print model.odes
#print zout

initb = np.logspace(-9,-5, 100)
N=[]
N2=[]
for i in range (len(initb)):
    model.parameters['bo'].value = initb[i]
    model2.parameters['bo'].value = initb[i]
    zout = odesolve(model,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
    zout2 = odesolve(model2,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
    N.append(zout['tabcc'][-1])
    N2.append(zout2['tabcc'][-1])

#pl.ion()    
pl.figure()
pl.ylim(0,2.5e-7)
#pl.plot(initb, N, 'o')
#pl.plot(initb, N)
pl.semilogx(initb, N, label='original')
pl.semilogx(initb, N2, label='remove ab')
pl.legend(loc=0)
pl.xlabel('TIMP2 at initial')
pl.ylabel('abcc at eq. state')
pl.show()

#other simulations
#pl.ion()
pl.figure()
pl.title('Time Course Simulation')
pl.plot(t2, zouty['ta'], label="a")
pl.plot(t2, zouty['tb'], label="b")
pl.plot(t2, zouty['tc'], label="c")
pl.plot(t2, zouty['tab'], label="ab")
pl.plot(t2, zouty['tbc'], label="bc")
pl.plot(t2, zouty['tcc'], label="cc")
pl.plot(t2, zouty['tabc'], label="abc")
pl.plot(t2, zouty['tbcc'], label="bcc")
pl.plot(t2, zouty['tabcc'], label="abcc")
pl.plot(t2, zouty['tabccb'], label="abccb")
pl.plot(t2, zouty['tabccba'], label="abccba")
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules")
pl.show()

#pl.ion()
pl.figure()
pl.xlim(0,.4)
pl.title('Type 1 (disappeared soon)')
pl.plot(t2, zouty['tb'], label="b")
pl.plot(t2, zouty['tbc'], label="bc")
pl.plot(t2, zouty['tbcc'], label="bcc")
pl.plot(t2, zouty['tbccb'], label="bccb")
pl.plot(t2, zouty['tabccb'], label="abccb")
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules")
pl.show()

#pl.ion()
pl.figure()
pl.title('Type 2 (disappeared after long time)')
pl.plot(t2, zouty['tab'], label="ab")
pl.plot(t2, zouty['tc'], label="c")
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules")
pl.show()

#pl.ion()
pl.figure()
pl.xlim(0,.4)
pl.title('Type 3 (Remain at some values)')
pl.plot(t2, zouty['ta'], label="a")
pl.plot(t2, zouty['tcc'], label="cc")
pl.plot(t2, zouty['tabc'], label="abc")
pl.plot(t2, zouty['tabcc'], label="abcc")
pl.plot(t2, zouty['tabccba'], label="abccba")
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules")
pl.show()



    