from pysb import * #calling all pysb classes and files
from pysb.integrate import Solver #importing Solver from integrate file
Model() #calling Model class but all input is by default, so we don't need to put any arguments in it (in core.py)
Monomer('A') #calling Monomer class (in core.py)
Parameter('w',3.0) #calling Parameter class (in core.py)
Rule('sA', None >> A(),w) #calling Rule class (in core.py)
t=[0,10,20,30]
solver = Solver(model,t) #calling Solver class (integrate.py)
solver.run() #calling run function in the Solver class (integrate.py)
print solver.y[:,1]