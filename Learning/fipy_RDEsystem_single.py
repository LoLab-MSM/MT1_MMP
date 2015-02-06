import fipy
 
"""create mess on the domain solution"""
m=fipy.Grid1D(nx=100, Lx=1.)          
v0 = fipy.CellVariable(name="molecule a", mesh=m, hasOld=True, value = 0.5)
v1 = fipy.CellVariable(name="molecule b", mesh=m, hasOld=True, value = 0.3)
v2 = fipy.CellVariable(name="molecule p", mesh=m, hasOld=True, value = 0.1)
       
"""define boundary conditions"""
#Dirichlet Boundary
v0.constrain(1., m.facesLeft)
v0.constrain(0., m.facesRight)
      
v1.constrain(1, m.facesLeft)
v1.constrain(0, m.facesRight)
      
v2.constrain(1, m.facesLeft)
v2.constrain(0, m.facesRight)
      
si = fipy.Viewer((v0, v1, v2))
si.plot()
        
"""equations"""
eqn0 = fipy.TransientTerm(var=v0) == fipy.DiffusionTerm(1, var=v0) #- __s0*__s1 +__s2 #+ model.odes[0]
eqn1 = fipy.TransientTerm(var=v1) == fipy.DiffusionTerm(1, var=v1) #- __s0*__s1 +__s2 #+ model.odes[1]
eqn2 = fipy.TransientTerm(var=v2) == fipy.DiffusionTerm(1, var=v2) #+ __s0*__s1 -__s2 #+ model.odes[1]
    
eqn = eqn0 & eqn1 & eqn2
  
"""Solve the equation"""
for t in range(10): 
     v0.updateOld()
     v1.updateOld()
     v2.updateOld()
     eqn.solve(dt=1e-3)
     si.plot()
if __name__ == '__main__':
     raw_input("Press <return> to proceed...")
print v2
# from fipy import Grid1D, CellVariable, TransientTerm, DiffusionTerm, Viewer
# m = Grid1D(nx=100, Lx=1.)
# v0 = CellVariable(mesh=m, hasOld=True, value=0.5)
# v1 = CellVariable(mesh=m, hasOld=True, value=0.5)
# eqn0 = TransientTerm(var=v0) == DiffusionTerm(1, var=v0)
# eqn1 = TransientTerm(var=v1) == DiffusionTerm(1, var=v1)
# eqn = eqn0 & eqn1
# vi = Viewer((v0, v1))
# for t in range(100): 
#      v0.updateOld()
#      v1.updateOld()
#      eqn.solve(dt=1)
#      vi.plot()
# if __name__ == '__main__':
#      raw_input("Press <return> to proceed...")