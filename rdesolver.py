#Solver for RDE using fipy (def or class?)
def rdesolve(model, mesh, initc): #Diffusifity, Dirichlet = None, Neumann = None):
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
     
     """TEST 1D 2 RDEs"""
     m = fipy.Grid1D(nx=100, Lx=1.) 
     v = fipy.CellVariable(mesh=m, hasOld=True, value=[[0.5], [0.5]], elementshape=(2,))
     v.constrain([[0], [1]], m.facesLeft)
     v.constrain([[1], [0]], m.facesRight)
     S0 = -v[0] * v[1]
     S1 = -v[0] * v[1]
     #source = [S0, S1]
     #print source
     eqn = fipy.TransientTerm([[1,0],[0,1]]) == fipy.DiffusionTerm([[[0.01,-1],[1,0.01]]])
     eqn2 = fipy.TransientTerm([[1,0],[0,1]]) == fipy.DiffusionTerm([[[1,0],[0,1]]]) #+source
     if __name__ == '__main__':
         vi = fipy.Viewer(vars=(v[0],v[1]))
     for t in range(10): 
         v.updateOld()
         eqn.solve(var=v, dt=1.e-3)
         if __name__ == '__main__':
             vi.plot()
     for t in range(10): 
         v.updateOld()
         eqn2.solve(var=v, dt=1.e-3)
         if __name__ == '__main__':
             vi.plot()
#      """"create mess on the domain solution"""
#      mesh = mesh
#      print mesh
#      if  mesh[0] == '1d':
#          m=fipy.Grid1D(nx=mesh[1], Lx=mesh[2])
#      elif  mesh[0] == '2d':
#            m=fipy.Grid2D(nx=mesh[1], dy=mesh[2], Lx=mesh[3], Ly=mesh[4])
#      else:
#          print mesh[0], 'is not recognized as mesh'
#          
#      """store solutions into variables"""
#      initc = initc    
#      X = fipy.CellVariable(mesh=m, hasOld=True, value = initc, elementshape = (len(model.species),)) 
#      
#      """define boundary conditions"""
#      Xleft = [0,0,0] ## dirichlet_condition at left
#      Xright = [0,0,0] ## dirichlet_condition at right
#      li=[]
#      ri=[]
#      for i in range(len(model.species)):
#          li.append(Xleft[i])
#          ri.append(Xright[i])
#      print li
#      X.constrain(li, m.facesLeft)
#      X.constrain(ri, m.facesRight)
#      
#      """Create Viewer"""
#      vi=[]
#      for i in range(len(model.species)):
#          vi.append(X[i])   
#      view = fipy.Viewer(vi)
#      
#      print view
#           
#      """equations"""
#      M=numpy.identity(len(model.species))
#      print M
#      equation = fipy.TransientTerm(M) == fipy.DiffusionTerm([M]) #+ Source Term??     
#      print equation
#      """Solve the equation"""
#      for t in range(10):
#          X.updateOld()
#          #equation.solve(var=X, dt=0.1)
#          view.plot()

     if __name__ == '__main__':
         raw_input("Explicit transient diffusion. Press <return> to proceed...")