#Solver for RDE using fipy (def or class?)
def rdesolve(model, mesh, initc, bound):
     """solve the system not using matrix system"""
     import fipy

     """create mess on the domain solution"""
     mesh = mesh
     print mesh
     if  mesh[0] == '1d':
         m=fipy.Grid1D(nx=mesh[1], Lx=mesh[2])
     elif  mesh[0] == '2d':
         m=fipy.Grid2D(nx=mesh[1], ny=mesh[2], Lx=mesh[3], Ly=mesh[4])
     else:
         print mesh[0], 'is not recognized as mesh'
#          
     """store solutions into variables"""
      #using different variable name?
     __s0 = fipy.CellVariable(name="molecule a", hasOld=True, mesh=m, value = 0.5)
     __s1 = fipy.CellVariable(name="molecule c", hasOld=True, mesh=m, value = 0.5)
     __s2 = fipy.CellVariable(name="molecule p", hasOld=True, mesh=m, value = 0.)
      
     """define boundary conditions"""
      #Dirichlet Boundary
     __s0.constrain(0, m.facesLeft)
     __s0.constrain(1, m.facesRight)
        
     __s1.constrain(0.5, m.facesLeft)
     __s1.constrain(0.8, m.facesRight)
      
     __s2.constrain(1, m.facesLeft)
     __s2.constrain(0, m.facesRight)
     
     si = fipy.Viewer((__s0, __s1, __s2))
     si.plot()
     print __s0  
     """equations"""
     eqn0 = fipy.TransientTerm(var=__s0) == fipy.DiffusionTerm(1, var=__s0) - 10*__s0*__s1 +__s2 #+ model.odes[0]
     eqn1 = fipy.TransientTerm(var=__s1) == fipy.DiffusionTerm(1, var=__s1) - 10*__s0*__s1 +__s2 #+ model.odes[1]
     eqn2 = fipy.TransientTerm(var=__s2) == fipy.DiffusionTerm(1, var=__s2) + 10*__s0*__s1 -__s2 #+ model.odes[1]
      
     eqn = eqn0 & eqn1 & eqn2
 
     """Solve the equation"""
     for t in range(20): 
         __s0.updateOld()
         __s1.updateOld()
         __s2.updateOld()
         eqn.solve(dt=1e-3)
         si.plot()