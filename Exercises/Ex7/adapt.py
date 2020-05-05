from ngsolve import *
from netgen.geom2d import *

maxh = 0.1

# generate the geometry and initial mesh
geo = SplineGeometry()

cos = [(0.7,0.45), (0.75,0.45), (0.8,0.45), (0.8, 0.55), (0.75,0.55), (0.7,0.55), (0.75,0.5)]
for co in cos:
    geo.AddPoint(*co)
geo.Append(["line", 0,1], leftdomain=3, rightdomain=1, bc='none'  , maxh=maxh)
geo.Append(["line", 1,6], leftdomain=3, rightdomain=3, bc='goal_2', maxh=maxh)
geo.Append(["line", 6,4], leftdomain=3, rightdomain=3, bc='goal_2', maxh=maxh)
geo.Append(["line", 4,5], leftdomain=3, rightdomain=1, bc='none'  , maxh=maxh)
geo.Append(["line", 5,0], leftdomain=3, rightdomain=1, bc='none'  , maxh=maxh)
geo.Append(["line", 1,2], leftdomain=3, rightdomain=1, bc='none'  , maxh=maxh)
geo.Append(["line", 2,3], leftdomain=3, rightdomain=1, bc='none'  , maxh=maxh)
geo.Append(["line", 3,4], leftdomain=3, rightdomain=1, bc='none'  , maxh=maxh)

geo.AddRectangle( (0,0), (1,1), leftdomain=1, rightdomain=0, bc='outer')
geo.AddRectangle( (0.2,0.45), (0.3,0.55), leftdomain=2, rightdomain=1,
                  bcs=['inner_bot', 'inner_right', 'inner_top', 'inner_left'])
geo.SetMaterial(1, 'outer')
geo.SetMaterial(2, 'inner')
geo.SetMaterial(3, 'goal_1')

mesh = Mesh(geo.GenerateMesh(maxh=0.1))
Draw(mesh)

# The first goal functional is already defined
def goal_fun1_1(fes):
    u,v = fes.TnT()
    f = LinearForm(fes)
    f += SymbolicLFI(100*v, VOL, definedon=fes.mesh.Materials('goal_1'))
    f.Assemble()
    return f.vec

# For the second goal functional use
# SymbolicLFI( ..., BND, definedon=fes.mesh.Materials('goal_1'))

# For the third functional use the following loop:

# for k,p in enumerate(fes.mesh.ngmesh.Points()):
#        if tuple(p)==(0.75,0.5,0):
#            f.vec[k] = 1.0

# (why is this a proper implementation of the functional? Remember which fes we use...)

# error-estimator
def err_est(fes_flux, gfflux, gfsol):
    lam = 1
    flux = lam * grad(gfsol)
    gfflux.Set (flux)
    err = 1/lam*(flux-gfflux)*(flux-gfflux)
    elerr = Integrate (err, fes_flux.mesh, VOL, element_wise=True)
    return elerr

order = 2
fes = H1(mesh, order=order, dirichlet='outer')    

u,v = fes.TnT()
a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v))
a.Assemble()

rhs = LinearForm(fes)
factor = CoefficientFunction([0,100])
rhs += SymbolicLFI(factor * v)
rhs.Assemble()

#Solution of primal formulation
sol = GridFunction(fes)
#Solution of dual formulation (for the goal functional)
dual_sol = GridFunction(fes)

Draw(sol)
#For error estimator
fes_flux = HDiv(mesh, order=order-1)
gfflux = GridFunction(fes_flux)

ref_flag = False
nmarked = 1
toterr = 1

max_ndof=15000
tol=5e-15
maxerr = 1

energyestimator = True

while fes.ndof < max_ndof and toterr > tol and nmarked > 0:
    print('continue, ', fes.ndof, maxerr)
    # do not refine the first time...
    if ref_flag:
        print('refine!')
        mesh.Refine()
    else:
        ref_flag = True

    # With the update function we tell all the objects that the mesh has changed
    # By this precalculated stuff like ndof gets updated
    fes_flux.Update()
    fes.Update()
    print('ndof now: ', fes.ndof)
    rhs.Assemble()
    a.Assemble()
    sol.Update()

    gfflux.Update()
    
    ainv = a.mat.Inverse(fes.FreeDofs(), inverse='sparsecholesky')

    # Calculate the solution of the primal formulation
    sol.vec.data = ainv * rhs.vec
    err_sol = err_est(fes_flux, gfflux, sol)

    # Calculate solution of dual formulation here
    # Define the product of the error estimator
    dual_sol.Update()
        
    # Maybe you want to change this to the global error of the functional
    toterr = sum(err_sol)    
    maxerr = max(err_sol)
    
    nmarked = 0
    if energyestimator:
        max_en_err = max(err_sol)
        for el in mesh.Elements():
            mesh.SetRefinementFlag(el, err_sol[el.nr] > 0.25*max_en_err)
            if err_sol[el.nr] > 0.25*maxerr:
                nmarked = nmarked + 1
        print(' ---------- marked', nmarked, 'elemenets of ', mesh.ne)
    else:
        print("add your code here")
        
    Redraw()
    

# Make pretty plots here...
