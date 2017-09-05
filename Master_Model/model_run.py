import model
import numpy as np
import pylab as pl
from pysb.integrate import odesolve

model1 = model.return_model('pysb')
model2 = model.return_model('old')

t=np.linspace(0,30000,10000)
initb = np.logspace(-9,-5, 100)

N=[]
N2=[]
for i in range (len(initb)):
    model1.parameters['A2_0'].value = initb[i]
    model2.parameters['A2_0'].value = initb[i]
    zout = odesolve(model1,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
    zout2 = odesolve(model2,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
    N.append(zout['tabcc'][-1])
    N2.append(zout2['tabcc'][-1])

pl.figure()
pl.ylim(0,2.5e-7)
#pl.plot(initb, N, 'o')
#pl.plot(initb, N)
pl.semilogx(initb, N, label='pysb')
pl.semilogx(initb, N2, label='old')
pl.legend(loc=0)
pl.xlabel('TIMP2 at initial')
pl.ylabel('abcc at eq. state')

t2=np.linspace(0,10,2000)
z1 = odesolve(model1,t2, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
z2 = odesolve(model2,t2, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)

pl.figure()
pl.xlim(0,0.1)
pl.title('a')
pl.plot(t2, z1['ta'], label="pysb")
pl.plot(t2, z2['ta'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")

pl.figure()
pl.title('b')
pl.plot(t2, z1['tb'], label="pysb")
pl.plot(t2, z2['tb'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")

pl.figure()
pl.xlim(0,0.1)
pl.title('c')
pl.plot(t2, z1['tc'], label="pysb")
pl.plot(t2, z2['tc'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")

pl.figure()
pl.title('ab')
pl.plot(t2, z1['tab'], label="pysb")
pl.plot(t2, z2['tab'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")

pl.figure()
pl.title('bc')
pl.plot(t2, z1['tbc'], label="pysb")
pl.plot(t2, z2['tbc'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")

pl.figure()
pl.title('cc')
pl.plot(t2, z1['tcc'], label="pysb")
pl.plot(t2, z2['tcc'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")

pl.figure()
pl.title('abc')
pl.plot(t2, z1['tabc'], label="pysb")
pl.plot(t2, z2['tabc'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")

pl.figure()
pl.title('bcc')
pl.plot(t2, z1['tbcc'], label="pysb")
pl.plot(t2, z2['tbcc'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")

pl.figure()
pl.title('abcc')
pl.plot(t2, z1['tabcc'], label="pysb")
pl.plot(t2, z2['tabcc'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")

pl.figure()
pl.title('bccb')
pl.plot(t2, z1['tbccb'], label="pysb")
pl.plot(t2, z2['tbccb'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")

pl.figure()
pl.title('abccb')
pl.plot(t2, z1['tabccb'], label="pysb")
pl.plot(t2, z2['tabccb'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")

pl.figure()
pl.title('abccba')
pl.plot(t2, z1['tabccba'], label="pysb")
pl.plot(t2, z2['tabccba'], label="old")
pl.legend()
pl.xlabel("Time")
pl.ylabel("Concentrations")
pl.show()
