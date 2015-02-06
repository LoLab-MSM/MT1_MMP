from fipy import *
dx=1.
mesh=Grid1D(dx=dx, nx=50)
phi = CellVariable(name="solution", mesh=mesh)
phi.value = 0.5
print phi
D=1.
phi.constrain(0, where=mesh.facesRight)
phi.constrain(1, where=mesh.facesLeft)
# print phi.faceValue
equation = TransientTerm()==ExplicitDiffusionTerm(coeff=D) + phi * 0.5
timestep = 0.9 * dx * dx**2 / (2*D)
steps=100
# phiAnalytical = CellVariable(name="analytical value", mesh=mesh)
# if __name__ == '__main__':
#      viewer = Viewer(vars=(phi, phiAnalytical), datamin=0., datamax=1.)
#      viewer.plot()
# x = mesh.cellCenters[0]
# t = timestep * steps
# try:
#      from scipy.special import erf # doctest: +SCIPY
#      phiAnalytical.setValue(1 - erf(x / (2 * numerix.sqrt(D * t)))) # doctest: +SCIPY
# except ImportError:
#      print "The SciPy library is not available to test the solution to \the transient diffusion equation"
# for step in range(steps):
#      equation.solve(var=phi,
#                dt=timestep)
#      if __name__ == '__main__':
#          viewer.plot()
# print phi.allclose(phiAnalytical, atol = 7e-4) # doctest: +SCIPY
#  
# if __name__ == '__main__':
#      raw_input("Explicit transient diffusion. Press <return> to proceed...")

# if __name__ == '__main__':
viewer = Viewer(vars=phi)
viewer.plot()
for i in range(steps):
     equation.solve(var=phi, dt=timestep)
#     if __name__ == '__main__':
     viewer.plot()
if __name__ == '__main__':
     raw_input("Explicit transient diffusion. Press <return> to proceed...")