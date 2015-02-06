from fipy import *
if __name__ == "__main__":
     nx = ny = 20
else:
     nx = ny = 10
mesh = Grid2D(nx=nx, ny=ny, dx=0.25, dy=0.25)
phi = CellVariable(name=r"$\phi$", mesh=mesh)
psi = CellVariable(name=r"$\psi$", mesh=mesh)
noise = GaussianNoiseVariable(mesh=mesh, mean=0.5, variance=0.01).value
phi[:] = noise
if __name__ == "__main__":
     viewer = Viewer(vars=(phi, psi)) # , datamin=0., datamax=1.)
var = CellVariable(mesh=mesh, elementshape=(2,))
print var[1]
var[0] = noise

if __name__ == "__main__":
     viewer = Viewer(vars=(var[0], var[1]))

D = a = epsilon = 1.
v0 = var[0]
dfdphi = a**2 * 2 * v0 * (1 - v0) * (1 - 2 * v0)
dfdphi_ = a**2 * 2 * (1 - v0) * (1 - 2 * v0)
d2fdphi2 = a**2 * 2 * (1 - 6 * v0 * (1 - v0))
source = var[1] * (0, 1)
impCoeff = -d2fdphi2 * ((0, 0),(1., 0)) + ((0, 0),(0, -1.))
eq = TransientTerm(((1., 0.),(0., 0.))) == DiffusionTerm([((0.,D),(-epsilon**2, 0.))]) + ImplicitSourceTerm(impCoeff) + source
