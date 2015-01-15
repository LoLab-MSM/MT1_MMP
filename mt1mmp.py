from pysb import *

Model()

Monomer('a',['a1'])
Monomer('b',['b1','b2'])
Monomer('c',['c1','c2'])

Parameter('kab', 2.1e7)
Parameter('kbc', 2.74e6)
Parameter('lbc', 2e-4)
Parameter('kcc', 2*2e6)
Parameter('lcc', 1e-2)

#binding criteria : (ab) b1 with a1, (bc) b2 with c1,(cc) c2 with itself
Rule('ab', a(a1=None) + b(b1=None) >> a(a1=1)%b(b1=1), kab)
Rule('bc', b(b2=None) + c(c1=None) <> b(b2=1)%c(c1=1), kbc, lbc)
Rule('cc', c(c2=None) + c(c2=None) <> c(c2=1)%c(c2=1), kcc, lcc)

#Rule('ab', a(a1=None) + b(b1=None, b2=None) >> a(a1=1)%b(b1=1), kab)

Initial(a(a1=None), Parameter('ao', 1e-6))
Initial(b(b1=None, b2=None), Parameter('bo', 1.57e-7))
Initial(c(c1=None, c2=None), Parameter('co', 1e-6))

Observable('ta', a(a1=None))
Observable('tb', b(b1=None, b2=None))
Observable('tc', c(c1=None, c2=None))
Observable('tab', a(a1=1)%b(b1=1, b2=None))
Observable('tbc', b(b1=None, b2=1)%c(c1=1, c2=None))
Observable('tcc', c(c1=None, c2=1)%c(c1=None, c2=1))
Observable('abc', a(a1=1)%b(b1=1, b2=2)%c(c1=2, c2=None))
Observable('bcc', b(b1=None, b2=1)%c(c1=1, c2=2)%c(c1=None, c2=2))
Observable('abcc', a(a1=1)%b(b1=1,b2=2)%c(c1=2,c2=3)%c(c1=None,c2=3))
Observable('bccb', b(b1=None, b2=1)%c(c1=1, c2=2)%c(c2=2, c1=3)%b(b2=3, b1=None))
Observable('abccb', a(a1=1)%b(b1=1, b2=2)%c(c1=2, c2=3)%c(c2=3, c1=4)%b(b2=4, b1=None))
Observable('abccba', a(a1=1)%b(b1=1, b2=2)%c(c1=2, c2=3)%c(c2=3, c1=4)%b(b2=4, b1=5)%a(a1=5))
