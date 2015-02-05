import fipy as fp
 
nx = ny = 20
 
mesh = fp.Grid2D(nx=nx, ny=ny, dx=0.25, dy=0.25)
#mesh = fp.Grid1D(nx=nx, dx=0.25)
 
var = fp.CellVariable(mesh=mesh, value=[[0.5],[0.5]], elementshape=(2,))
D = a = epsilon = 1.
v0 = var[0]
dfdphi = a**2 * 2 * v0 * (1 - v0) * (1 - 2 * v0)
dfdphi_ = a**2 * 2 * (1 - v0) * (1 - 2 * v0)
d2fdphi2 = a**2 * 2 * (1 - 6 * v0 * (1 - v0))
 
source = (- d2fdphi2 * v0 + dfdphi) * (0, 1)
impCoeff = -d2fdphi2 * ((0, 0), 
                        (1., 0)) + ((0, 0), 
                                    (0, -1.))
 
eq = fp.TransientTerm(((1., 0.), 
                       (0., 1.))) == fp.DiffusionTerm([((1.,          0), 
                                                        (0, 1))]) + source
 
print 'class:',source.__class__
print 'shape:',source.shape
print 'rank:',source.rank
 
dexp = -5
elapsed = 0.
duration = .5e-1
 
while elapsed < duration:
    dt = min(100, fp.numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    eq.solve(var=var, dt=dt)
 