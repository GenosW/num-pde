from netgen.csg import *
from ngsolve import *
from ngsolve.solvers import CG

from ngsolve.internal import visoptions
from ngsolve.internal import viewoptions

def MakeGeometry():
    geometry = CSGeometry()
    box = OrthoBrick(Pnt(-1,-1,-1),Pnt(2,1,2)).bc("outer")

    core = OrthoBrick(Pnt(0,-0.05,0),Pnt(0.8,0.05,1))- \
           OrthoBrick(Pnt(0.1,-1,0.1),Pnt(0.7,1,0.9))- \
           OrthoBrick(Pnt(0.5,-1,0.4),Pnt(1,1,0.6)).maxh(0.2).mat("core")
    
    coil = (Cylinder(Pnt(0.05,0,0), Pnt(0.05,0,1), 0.3) - \
            Cylinder(Pnt(0.05,0,0), Pnt(0.05,0,1), 0.15)) * \
            OrthoBrick (Pnt(-1,-1,0.3),Pnt(1,1,0.7)).maxh(0.2).mat("coil")
    
    geometry.Add ((box-core-coil).mat("air"))
    geometry.Add (core)
    geometry.Add (coil)
    return geometry

ngmesh = MakeGeometry().GenerateMesh(maxh=0.5)
mesh = Mesh(ngmesh)

#ngsglobals.msg_level = 5

fes = HCurl(mesh, order=1, dirichlet="outer")

u = fes.TrialFunction()
v = fes.TestFunction()

mur = { "core" : 1000, "coil" : 1, "air" : 1 }
mu0 = 1.257e-6

nu_coef = [ 1/(mu0*mur[mat]) for mat in mesh.GetMaterials() ]
print ("nu_coef=", nu_coef)

eps = 1
nu = CoefficientFunction(nu_coef)

a = BilinearForm(fes, symmetric=True)
a += nu*curl(u)*curl(v)*dx + eps*nu*u*v*dx

# the jacobi preconditioner
c = Preconditioner(a, type="local")
# the direct inverse 
# c = Preconditioner(a, type="direct")

f = LinearForm(fes)
f += CoefficientFunction((y,0.05-x,0)) * v * dx("coil")

blocks = []
# create your blocks here

gfu = GridFunction(fes)

# stores the error and iteration number
data = []

# needed for this...
callback = lambda k,r: data.append((k,r))

with TaskManager():
    a.Assemble()    
    f.Assemble()
    
    solvers.CG(mat=a.mat, pre=c.mat, sol = gfu.vec, rhs = f.vec, tol = 1e-8, maxsteps = 10000, callback = callback, printrates = True)
    print(data [-1])
    # Use your AFW preconditioner instead of the Jacobi
    # solvers.CG(mat=a.mat, pre=AFW, sol = gfu.vec, rhs = f.vec, tol = 1e-8, maxsteps = 10000, callback = callback)

Draw (gfu.Deriv(), mesh, "B-field", draw_surf=False)

visoptions.clipsolution = 'vec'
viewoptions.clipping.ny= -1
viewoptions.clipping.enable = 1

