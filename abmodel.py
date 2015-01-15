from pysb import *


Model()

Monomer('a',['b'])
Monomer('c',['b'])

Parameter('k',1e-3)
Parameter('l',1e-10)

Rule('d', a(b=None) + c(b=None) <> a(b=1)%c(b=1), k, l)

Parameter('ao',100)
Parameter('co',200)
Initial(a(b=None),ao)
Initial(c(b=None),co)

Observable('dt', a(b=1)%c(b=1))
Observable('at', a(b=None))
Observable('ct', c(b=None))

