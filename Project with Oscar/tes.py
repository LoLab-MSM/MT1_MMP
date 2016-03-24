from sympy.abc import x, y, a, b, c, d
from sympy.solvers.ode import dsolve, checkodesol
from sympy import pprint, Function
f = Function('f')
y = f(x)
genform = y.diff(x) - (b*y**2 + c*y + d)
sol = dsolve(genform, y)
pprint(sol, wrap_line=False)