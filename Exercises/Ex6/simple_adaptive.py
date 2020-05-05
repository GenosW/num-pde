from ngsolve import *
from netgen.geom2d import SplineGeometry
from math import pi

geo = SplineGeometry()
geo.AddRectangle((-1,-1),(1,1),bc=1)
mesh = Mesh( geo.GenerateMesh(maxh=0.5))

V = H1(mesh, order=3, dirichlet=[1])

u = V.TrialFunction()
v = V.TestFunction()

r = sqrt(x*x+y*y)
r2 = x*x+y*y
R = 0.5

alpha = IfPos(r-R,2*log(1/R),1)

solution = IfPos(r-R,1.0/16.0*log(r2)/log(R*R),1.0/8.0-r2/4.0)
sol_flux = alpha * IfPos(r-R,1.0/16.0/log(R*R)*
                     CoefficientFunction((2*x/r2,2*y/r2)),
                     CoefficientFunction((-x/2.0,-y/2.0)))
rhs = IfPos(r-R,0,1)

a = BilinearForm(V, symmetric=False)
a += SymbolicBFI(alpha*grad(u)*grad(v))

f = LinearForm(V)
f += SymbolicLFI(rhs*v)

gfu = GridFunction(V)
flux = alpha * grad(gfu)
Draw (gfu,mesh,"u")
Draw (grad(gfu),mesh,"grad_u")

def SolveBVP():
    V.Update()
    gfu.Update()
    a.Assemble()
    f.Assemble()
    gfu.Set(solution)
    f.vec.data -= a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(V.FreeDofs(),"umfpack") * f.vec
    Redraw (blocking=True)

l = []

def CalcError():
    
    # compute error on every element:
    err = 1/alpha*(flux - sol_flux)*(flux - sol_flux)
    elerr = Integrate (err, mesh, VOL, element_wise=True)
    
    # sort elements (corresponding to error contribution)
    err_and_el_sorted = sorted([(entry,i) for i, entry in enumerate(elerr)], key= lambda x:-x[0])
    # reset marks
    marks = [False for el in mesh.Elements()]

    # mark element with largest error until 10% of the error is on marked elements:
    sumerr = sum(elerr)    
    accsum = 0
    for err,el in err_and_el_sorted:
        if accsum < 0.1 * sumerr:
            marks[el] = True
            accsum += err
        else:
            break
    
    print ("V.ndof = ", V.ndof)
    H1error = sqrt(Integrate (1/alpha*(flux - sol_flux)*(flux - sol_flux), mesh, VOL))
    print ("weighted H1 (semi norm) error = ", H1error)
    L2error = sqrt(Integrate ((gfu - solution)*(gfu - solution), mesh, VOL))
    print ("L2 error = ", L2error)
    
    # call the refinement according to the marks:
    for el in mesh.Elements():
        mesh.SetRefinementFlag(el, marks[el.nr])

with TaskManager():
    while V.ndof < 10000:  
        SolveBVP()
        CalcError()
        mesh.Refine()
    SolveBVP()
