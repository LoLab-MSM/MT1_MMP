from pysb.integrate import odesolve
import numpy as np
import pylab as pl
import mt1_mmp_model
from pysb.bng import generate_equations
#from MT1_MMP import mt1_mmp_model as mt1mmp2

model=mt1_mmp_model.return_model('original')
model2=mt1_mmp_model.return_model('abremoved')
model3=mt1_mmp_model.return_model('abc1removed')
model4=mt1_mmp_model.return_model('abc2removed')

# generate_equations(model)
# print model.species
# 
# 
# import stoichiometry_analysis as sto
# conserv_laws = sto.conservation_relations(model)
# 
# 
# print conserv_laws


t=np.linspace(0,30000,10000)
# t2=np.linspace(0,5,500)
# zouty = odesolve(model,t2, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
#print model.odes
#print zout

#to simulate abcc(infinity) vs b(0)
initb = np.logspace(-9,-5, 100)
N=[]
N2=[]
N3=[]
N4=[]
for i in range (len(initb)):
    model.parameters['bo'].value = initb[i]
    model2.parameters['bo'].value = initb[i]
    model3.parameters['bo'].value = initb[i]
    model4.parameters['bo'].value = initb[i]
    zout = odesolve(model,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
    zout2 = odesolve(model2,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
    zout3 = odesolve(model3,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
    zout4 = odesolve(model4,t, integrator='vode', with_jacobian=True, rtol=1e-20, atol=1e-20)
    N.append(zout['tabcc'][-1])
    N2.append(zout2['tabcc'][-1])
    N3.append(zout3['tabcc'][-1])
    N4.append(zout4['tabcc'][-1])

#pl.ion()    
pl.figure()
pl.ylim(0,2.5e-7)
#pl.plot(initb, N, 'o')
#pl.plot(initb, N)
pl.semilogx(initb, N, label='original')
pl.semilogx(initb, N2, label='remove ab')
pl.semilogx(initb, N3, label='remove abc')
pl.semilogx(initb, N4, label='remove total of abc')
pl.legend(loc=0)
pl.xlabel('TIMP2 at initial')
pl.ylabel('abcc at eq. state')
pl.show()

##other simulations
##pl.ion()
#pl.figure()
#pl.title('Time Course Simulation')
#pl.plot(t2, zouty['ta'], label="a")
#pl.plot(t2, zouty['tb'], label="b")
#pl.plot(t2, zouty['tc'], label="c")
#pl.plot(t2, zouty['tab'], label="ab")
#pl.plot(t2, zouty['tbc'], label="bc")
#pl.plot(t2, zouty['tcc'], label="cc")
#pl.plot(t2, zouty['tabc'], label="abc")
#pl.plot(t2, zouty['tbcc'], label="bcc")
#pl.plot(t2, zouty['tabcc'], label="abcc")
#pl.plot(t2, zouty['tabccb'], label="abccb")
#pl.plot(t2, zouty['tabccba'], label="abccba")
#pl.legend()
#pl.xlabel("Time (s)")
#pl.ylabel("Molecules")
#pl.show()
#
##pl.ion()
#pl.figure()
#pl.xlim(0,.4)
#pl.title('Type 1 (disappeared soon)')
#pl.plot(t2, zouty['tb'], label="b")
#pl.plot(t2, zouty['tbc'], label="bc")
#pl.plot(t2, zouty['tbcc'], label="bcc")
#pl.plot(t2, zouty['tbccb'], label="bccb")
#pl.plot(t2, zouty['tabccb'], label="abccb")
#pl.legend()
#pl.xlabel("Time (s)")
#pl.ylabel("Molecules")
#pl.show()
#
##pl.ion()
#pl.figure()
#pl.title('Type 2 (disappeared after long time)')
#pl.plot(t2, zouty['tab'], label="ab")
#pl.plot(t2, zouty['tc'], label="c")
#pl.legend()
#pl.xlabel("Time (s)")
#pl.ylabel("Molecules")
#pl.show()
#
##pl.ion()
#pl.figure()
#pl.title('Type 3 (Remain at some values)')
#pl.plot(t2, zouty['ta'], label="a")
#pl.plot(t2, zouty['tcc'], label="cc")
#pl.plot(t2, zouty['tabc'], label="abc")
#pl.plot(t2, zouty['tabcc'], label="abcc")
#pl.plot(t2, zouty['tabccba'], label="abccba")
#pl.legend()
#pl.xlabel("Time (s)")
#pl.ylabel("Molecules")
#pl.show()
#
#print("Species Number")
#for i in range(len(model.species)):
#    print "S", i, model.species[i]
#
#print("ODE Model of Original System")
#for i in range(len(model.odes)):
#    print "d/dt S", i, model.odes[i]
#    
#print("Species Number")    
#for i in range(len(model2.species)):
#    print "S", i, model2.species[i]
#
#print("ODE Model of AB removed")
#for i in range(len(model2.odes)):
#    print "d/dt S", i, model2.odes[i]


    