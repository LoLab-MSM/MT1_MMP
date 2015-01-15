from pysb import *

Model()

Monomer('b', ['b1'])
Monomer('c', ['c1', 'c2'])

Parameter('kbc', 2.74e6)
Parameter('lbc', 2e-4)
Parameter('kcc', 2*2e6)
Parameter('lcc', 1e-2)

Rule('bc', b(b1=None) + c(c1=None) <> b(b1=1)%c(c1=1), kbc, lbc)
Rule('cc', c(c2=None) + c(c2=None) <> c(c2=1)%c(c2=1), kcc, lcc)



Initial(b(b1=None), Parameter('bo', 1.57e-7))
Initial(c(c1=None, c2=None), Parameter('co', 1e-6))

Observable('tb', b(b1=None))
Observable('tc', c(c1=None, c2=None))
Observable('tbc', b(b1=1)%c(c1=1,c2=None))
Observable('tcc', c(c1=None,c2=1)%c(c1=None,c2=1))
Observable('bcc', b(b1=1)%c(c1=1, c2=2)%c(c1=None,c2=2))
Observable('bccb', b(b1=1)%c(c1=1, c2=2)%c(c2=2, c1=3)%b(b1=3))