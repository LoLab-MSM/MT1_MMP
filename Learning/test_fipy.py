from fipy import *
dx=1.
mesh=Grid1D(dx=dx, nx=50)
X_1 = CellVariable(name="solution", mesh=mesh, value=0.)
D=1.
X_1.constrain(0, mesh.facesRight)
X_1.constrain(1, mesh.facesLeft)
equation = TransientTerm()==ExplicitDiffusionTerm(coeff=D)
timestep = 0.9 * dx * dx**2 / (2*D)
steps=100
#X_1_analytical = CellVariable(name="analytical value", mesh=mesh)
if __name__ == '__main__':
    viewer = Viewer(vars=X_1, datamin=0., datamax=1.)
    viewer.plot()
#x = mesh.cellCenters[0]
#t = timestep * steps
#try:
#     from scipy.special import erf 
#     X_1_analytical.setValue(1 - erf(x / (2 * numerix.sqrt(D * t)))) 
#except ImportError:
#     print "The SciPy library is not available to test the solution to \
#     the transient diffusion equation"
for i in range(steps):
    equation.solve(var=X_1, dt=timestep)
    if __name__ == '__main__':
        viewer.plot()
#print X_1.allclose(X_1_analytical, atol = 7e-4)
#if __name__ == '__main__':
#     raw_input("Explicit transient diffusion. Press <return> to proceed...")
