from ngsolve import *
from netgen.csg import *

box = OrthoBrick ( Pnt(-1,-1,-1), Pnt(1,1,1) ).bc("outer")
infcoil = Cylinder( Pnt(0,0,0.2), Pnt(0,0,0.4), 0.1 ) - Cylinder( Pnt(0,0,0.2), Pnt(0,0,0.4), 0.05 )
coil = infcoil * Plane ( Pnt(0,0,0.4), Vec(0,0,1) ) * Plane ( Pnt(0,0,0.2), Vec(0,0,-1) )

geo = CSGeometry()
geo.Add ( (box-coil).mat("air") )
geo.Add (coil.mat("coil"))

mesh = Mesh (geo.GenerateMesh(maxh=0.3) )

Draw (mesh)


mesh.Curve(3)

fes = HCurl(mesh, order=3)

u = fes.TrialFunction()
v = fes.TestFunction()

mu0 = 1.257e-6

# the 1*u*v is a regularisation as in example 11.3
a = BilinearForm(fes, symmetric=True)
a += SymbolicBFI ( 1/mu0 * curl(u) * curl(v) + 1 * u*v )

c = Preconditioner(a, "bddc")
a.Assemble()

current = 1 / sqrt(x*x+y*y) * CoefficientFunction ( (-y, x, 0) )

f = LinearForm(fes)
f += SymbolicLFI (current * v, definedon=mesh.Materials("coil"))
f.Assemble()


gfu = GridFunction (fes)

inva = CGSolver(a.mat, c.mat)
gfu.vec.data = inva * f.vec

Bfield = curl(gfu)

Draw (Bfield, mesh, "B-field")
