#Solver for RDE using fipy (def or class?)
def rdesolve(model, mesh, initc, bound): #Diffusifity, Dirichlet = None, Neumann = None):
     """Integrate a RDE model over a given time t using fipy package
     Parameters
     ----------
     model : pysb.Model
         ODE model as a part of RDE model.
     mesh : tuple
         in which dimension we want to create the mess. 1D or 2D or 3D or sphere
     initc : vector-valued
         initial condition for each variables
     Dirichlet : vector-valued
         Dirichlet condition for the boundary
     Neumann : vector-valued
         Neumann condition for the boundary
     Diffusivity : vector-valued
         Diffusion constant values for every equation.
     time : integer
         max time value.
     
     dirichlet_condition_left / dirichlet_condition_right : vector-valued
         boundary conditions
     
     """
     import fipy
     import numpy
     import sympy
     import re
     
     """TEST 1D 2 RDEs"""
#      nx = ny = 100
#      m = fipy.Grid2D(dx=0.25, dy=0.25, nx=nx, ny=ny)
#      m = fipy.Grid1D(nx=100, Lx=1.) 
#      v = fipy.CellVariable(mesh=m, hasOld=True, value=[[0.5], [0.5]], elementshape=(2,))
#      v.constrain([[0], [1]], m.facesLeft)
#      v.constrain([[1], [0]], m.facesRight)
#      S0 = -v[0] * v[1]
#      S1 = -v[0] * v[1]
#      #source = [S0, S1]
#      #print source
#      eqn = fipy.TransientTerm([[1,0],[0,1]]) == fipy.DiffusionTerm([[[0.01,-1],[1,0.01]]])
#      eqn2 = fipy.TransientTerm([[1,0],[0,1]]) == fipy.DiffusionTerm([[[1,0],[0,1]]]) #+source
#      if __name__ == '__main__':
#          vi = fipy.Viewer(vars=(v[0],v[1]))
#      for t in range(10): 
#          v.updateOld()
#          eqn.solve(var=v, dt=1.e-3)
#          if __name__ == '__main__':
#              vi.plot()
# #    for t in range(10): 
# #         v.updateOld()
# #         eqn2.solve(var=v, dt=1.e-3)
# #         if __name__ == '__main__':
# #         vi.plot()
# 
#      if __name__ == '__main__':
#          raw_input("Explicit transient diffusion. Press <return> to proceed...")


     """"create mess on the domain solution"""
     mesh = mesh
     if  mesh[0] == '1d':
         m=fipy.Grid1D(nx=mesh[1], Lx=mesh[2])
     elif  mesh[0] == '2d':
         m=fipy.Grid2D(nx=mesh[1], dy=mesh[2], Lx=mesh[3], Ly=mesh[4])
     else:
         print mesh[0], 'is not recognized as mesh, try 1d or 2d or 3d'
         quit()
          
     """store solutions into variables"""
     initc = initc    
     X = fipy.CellVariable(mesh=m, hasOld=True, value = initc, elementshape = (len(model.species),)) 
      
     """define boundary conditions"""
     if bound <> 'None':
         bound = bound
         if  mesh[0] == '1d':
             X.constrain(bound[0], where=m.facesLeft)
             X.constrain(bound[1], where=m.facesRight)
         elif  mesh[0] == '2d':
             X.constrain(bound[0], where=m.facesLeft)
             X.constrain(bound[1], where=m.facesBottom)
             X.constrain(bound[2], where=m.facesRight)
             X.constrain(bound[3], where=m.facesTop)                  

     """Create Viewer"""
     vi=[]
     for i in range(len(model.species)):
         vi.append(X[i])   
     view = fipy.Viewer(vars=vi)
               
     """equations"""
     stringodes = ';'.join(['%s' % sympy.ccode(model.odes[i]) for i in range(len(model.odes))])

     for j in range(len(model.odes)):
         c_odes = re.sub(r'_+s%d' % j,'X[%d]' % j, stringodes)
         stringodes = c_odes

     odes_m= c_odes.split(';')
     print odes_m
     print type(odes_m)
     
     
     M=numpy.identity(len(model.species))
#      M=[[1,0,0], [0,1,0], [0,0,1]]
     N=numpy.reshape(M, (len(model.species), len(model.species), 1))
#      odes = [-X[0]*X[1] + X[2], -X[0]*X[1] + X[2], X[0]*X[1] - X[2]]
#      source = odes_m[0] * N[0] + odes_m[1] * N[1] + odes_m[2] * N[2]
     equation = fipy.TransientTerm(M) == fipy.DiffusionTerm([M]) #+ (source)    
      
     """Solve the equation"""
     for t in range(10):
         X.updateOld()
         equation.solve(var=X, dt=1.e-3)
         view.plot()
      