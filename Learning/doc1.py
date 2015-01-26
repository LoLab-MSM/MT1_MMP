from pysb import *
from pysb.integrate import Solver
Model()
Monomer('A')
Parameter('w',3.0)
Rule('sA', None >> A(),w)
t=[0,10,20,30]
solver = Solver(model,t)
solver.run()
print solver.y[:,1]