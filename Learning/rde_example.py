from fipy import *
nx = 20
ny = nx
dx = 1.
dy = dx
L = dx * nx
mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
phi = CellVariable(name = "solution variable",mesh = mesh,value = 0.)
D = 1.
eq = TransientTerm() == DiffusionTerm(coeff=D)
valueTopLeft = 0
valueBottomRight = 1
X, Y = mesh.faceCenters
facesTopLeft = ((mesh.facesLeft & (Y > L / 2))| (mesh.facesTop & (X < L / 2)))
facesBottomRight = ((mesh.facesRight & (Y < L / 2))| (mesh.facesBottom & (X > L / 2)))
phi.constrain(valueTopLeft, facesTopLeft)
phi.constrain(valueBottomRight, facesBottomRight)
if __name__ == '__main__':
     viewer = Viewer(vars=phi, datamin=0., datamax=1.)
     viewer.plot()
timeStepDuration = 10 * 0.9 * dx**2 / (2 * D)
steps = 200
for step in range(steps):
     eq.solve(var=phi,dt=timeStepDuration)
     if __name__ == '__main__':
         viewer.plot()
if __name__ == '__main__':
     raw_input("Coupled equations. Press <return> to proceed...")
# if __name__ == "__main__":
#      nx = ny = 200
# else:
#      nx = ny = 200
# mesh = Grid2D(nx=nx, ny=ny, dx=0.25, dy=0.25)
# phi = CellVariable(name=r"$\phi$", mesh=mesh)
# psi = CellVariable(name=r"$\psi$", mesh=mesh)
# noise = GaussianNoiseVariable(mesh=mesh,mean=0.5,variance=0.01).value
# phi[:] = noise
# if __name__ == "__main__":
#      viewer = Viewer(vars=(phi, psi)) # , datamin=0., datamax=1.)
# D = a = epsilon = 1.
# dfdphi = a**2 * 2 * phi * (1 - phi) * (1 - 2 * phi)
# dfdphi_ = a**2 * 2 * (1 - phi) * (1 - 2 * phi)
# d2fdphi2 = a**2 * 2 * (1 - 6 * phi * (1 - phi))
# eq1 = (TransientTerm(var=phi) == DiffusionTerm(coeff=D, var=psi))
# eq2 = (ImplicitSourceTerm(coeff=1., var=psi)== ImplicitSourceTerm(coeff=-d2fdphi2, var=phi) - d2fdphi2 * phi + dfdphi- DiffusionTerm(coeff=epsilon**2, var=phi))
# eq3 = (ImplicitSourceTerm(coeff=1., var=psi) == ImplicitSourceTerm(coeff=dfdphi_, var=phi)- DiffusionTerm(coeff=epsilon**2, var=phi))
# eq = eq1 & eq2
# dexp = -5
# elapsed = 0.
# if __name__ == "__main__":
#      duration = 1#.5e-1
# else:
#      duration = 1#.5e-1
# while elapsed < duration:
#      dt = min(100, numerix.exp(dexp))
#      elapsed += dt
#      dexp += 0.01
#      eq.solve(dt=dt)
#      if __name__ == "__main__":
#          viewer.plot()
# if __name__ == '__main__':
#      raw_input("Coupled equations. Press <return> to proceed...")