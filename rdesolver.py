#Solver for RDE using fipy (def or class?)
"""paste this code into pysb.integrate.py"""
#not finished
def rdesolve(model):
#Diffusion constant values for every equation. Zero value must be determined as zero
#(model, diffusivity, dirichlet_condition_left, dirichlet_condition_right, initial_condition, time, ...?)
     """Integrate a RDE model over a given time t using fipy package
     Parameters
     ----------
     model : pysb.Model
         ODE model as a part of RDE model.
     time : integer
         max time value.
     diffusivity : vector-valued
         Diffusion constant values for every equation.
     initial_condition : vector-valued
         initial condition for each variables
     dirichlet_condition_left / dirichlet_condition_right : vector-valued
         boundary conditions
     
     """
     from fipy import *
     """"create mess on the domain solution"""
     m=Grid1D(nx=100, Lx=1.)
     
     """store solutions into variables"""
     #using different variable name?
     __s0 = CellVariable(name="molecule a", hasOld=True, mesh=m, value = 0.5)
     __s1 = CellVariable(name="molecule c", hasOld=True, mesh=m, value = 0.5)
     __s2 = CellVariable(name="molecule p", hasOld=True, mesh=m, value = 0.)
     
     #X as vector solution #Better!
     initial = [0.5, 0.5, 0] #?? how to write this matrix  outside integrate.py and pass the information to pdesolve function
     xi=[]
     for i in range(len(model.species)):
         xi.insert(len(xi),initial[i-1])       
     X = CellVariable(mesh=m, hasOld=True, value = xi, elementshape = (len(model.species),)) 
                                        
     
     """define boundary conditions"""
     #Dirichlet Boundary
     __s0.constrain(0, m.facesLeft)
     __s0.constrain(0, m.facesRight)
         
     __s1.constrain(0, m.facesLeft)
     __s1.constrain(0, m.facesRight)
     
     __s2.constrain(0, m.facesLeft)
     __s2.constrain(0, m.facesRight)
     
     #*define Dirichlet Boundary Conditions as vector #Better!
     Xleft = [0,0] ## dirichlet_condition at left
     Xright = [0,0] ## dirichlet_condition at right
     li=[]
     ri=[]
     for i in range(len(model.species)):
         li.insert(len(li),Xleft[i-1])
         ri.insert(len(ri),Xright[i-1])
     X.constrain(li, m.facesLeft)
     X.constrain(ri, m.facesRight)
     
     """Create Viewer"""
     si = Viewer((__s0, __s1, __s2))
     
     #*another way to define Viewer ??how to add tuple
     vi=[]
     for i in range(len(model.species)):
         vi.insert(len(vi),X[i-1])   
     sol = Viewer(vi)
          
     """equations"""
     eqn0 = TransientTerm(var=__s0) == DiffusionTerm(1, var=__s0) - __s0*__s1 +__s2 #+ model.odes[0]
     eqn1 = TransientTerm(var=__s1) == DiffusionTerm(1, var=__s1) - __s0*__s1 +__s2#+ model.odes[1]
     eqn2 = TransientTerm(var=__s2) == DiffusionTerm(1, var=__s2) + __s0*__s1 -__s2#+ model.odes[1]
     
     eqn = eqn0 & eqn1 & eqn2
     
     #*define equations as vector ?? big system means big matrix. effective and efisient?
     #define identity matrix M with dimension n:len(model.species)
     M=numpy.identity(len(model.species))
     
     equation = TransientTerm(M) == DiffusionTerm([M]) #+ Source Term??     
     
     """Solve the equation"""
     for t in range(10): 
         __s0.updateOld()
         __s1.updateOld()
         __s2.updateOld()
         eqn.solve(dt=1)
         si.plot()
         
     #Solve the system using vector field
     for t in range(10):
         X.updateOld()
         equation.solve(var=X, dt=1)
         sol.plot()