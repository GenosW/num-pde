# pylint: disable=unused-wildcard-import
import netgen.gui
from ngsolve import *
from netgen.geom2d import SplineGeometry, MakeRectangle
from time import time_ns

print("Welcome, Peter!")
### Definition of mesh ###
h = 0.5
geo = SplineGeometry()

geo.AddRectangle((0,0), (1,1), bcs=['b','r','t','l'], leftdomain=1, rightdomain=0)
geo.AddRectangle((0.3,0.5), (0.5,0.7), bcs=['b1','r1','t1','l1'], leftdomain=2, rightdomain=1)
geo.SetMaterial(1, "d2")
geo.SetMaterial(2, "d1")
mesh = Mesh(geo.GenerateMesh(maxh=h)) # standard mesh available
print('#'*40)
print("Mesh generated:")
print('#'*40)
print("Number of vertices: ", mesh.nv)
print("Number of elements: ", mesh.ne)
print("Dimension of mesh: ", mesh.dim)
lams = {"d1" : 1, "d2" : 10}
Lambda = CoefficientFunction([lams[mat] for mat in mesh.GetMaterials()])
cs = {"d1" : 1, "d2" : 0}
C = CoefficientFunction([cs[mat] for mat in mesh.GetMaterials()])
Draw(mesh)

### Definition of FES ###


# define discrete conforming spaces of Sigma and V
# any idea why we need a different order here?
# try to use different combinations later as well

order_V = 1
order_Sigma = order_V + 2
Sigma = HDiv(mesh, order = order_Sigma)
V = L2(mesh, order = order_V, dirichlet='b|r|t|l') # define dirichlet on V...but since we have hom dirichlet everywhere, it actually isnt necessary.
# -) Using V(order=1) and Sigma(order=3) actually works quite well too! 
#       -> gives similiar (better?) Wtot and runs faster (1s vs 0.33s)


#This defines a "compound FESpace"
X = FESpace([Sigma, V])

#This gives you a trial function sigma in Sigma and u in V
(sigma, u) = X.TrialFunction()
# Similiar for test functions
(tau, v) = X.TestFunction()
# normal vector needed for the fluxes
n = specialcf.normal(mesh.dim) # mesh.dim == 2

# define a "big" bilinearform on X which includes all the integrals of (4)
# you can use div(sigma) to get the divergence
B = BilinearForm(X)
#   <           (A)                 > + <       B    >
B += (-1/Lambda*sigma*tau + div(tau)*u + div(sigma)*v)*dx

# and a linearform on X
F = LinearForm(X)
#    < B >  + <         A              > ... not sure if A part is necessary since it's zero --> probably not
F += C*v*dx #+  0*(tau.Trace()*n)*ds

# solution on X
sol = GridFunction(X, 'sol_mixed')
sigma_sol = sol.components[0] #this gives you the sigma solution
u_sol = sol.components[1] #this gives you the u solution

### Solving the system ###
# Assemble the system of equations
with TaskManager(): # with Multithreading
    F.Assemble()
    B.Assemble()

# Solve the system
# NOTE: the Bilinearform B is NOT coercive -> not SPD, hence we can not use a sparsecholesky solver
# use: ... B.mat.Inverse(X.FreeDofs(), inverse = "umfpack")
# Calculate the solution field (function)
start = time_ns() # super accurate *cough-cough*
sol.vec.data = B.mat.Inverse(X.FreeDofs(), inverse = "umfpack") * F.vec
end = time_ns()
runtime = end - start
print("Done with solving!")

print("Trying to draw what you cobbled together here...")
Draw(sigma_sol, mesh, 'sigma')
Draw(u_sol, mesh, 'u')

### Calculating stuff in the system ###
# evaluate the flux using the sigma solution!
# Similar as for the H1, the var framework from the cont setting transfers to the discrete setting
# For sigma you can now evaluate sigma * n (since the normal trace is a cont operator for the H(div)

# Ω1 = (0.3,0.5)×(0.5,0.7)
BND_D1 = 'b1|r1|t1|l1' # boundary of Ω1
Wl = -Integrate(sigma_sol*n, mesh, definedon=mesh.Boundaries('l1'))
Wr = -Integrate(sigma_sol*n, mesh, definedon=mesh.Boundaries('r1'))
Wt = -Integrate(sigma_sol*n, mesh, definedon=mesh.Boundaries('t1'))
Wb = -Integrate(sigma_sol*n, mesh, definedon=mesh.Boundaries('b1'))
Wtot = Wl + Wr + Wl + Wb
Wtot2 = -Integrate(sigma_sol*n, mesh, definedon=mesh.Boundaries(BND_D1))
print('-'*30)
print("Total flux= ", Wtot, "(should be 0.04)")
print("Total flux= ", Wtot2)

print()
print('-'*30)
print('order_V =',order_V)
print('order_Sigma =',order_Sigma)
print('Runtime = ', runtime/1e+9)

# Let's try to use Trace() somehow...
# print("Let's try to use Trace() somehow...")
# Flux = GridFunction(Sigma)
# Flux += sigma_sol*n*ds#, definedon=mesh.Boundaries(BND_D1))
# print(Integrate(Flux, mesh, definedon=mesh.Boundaries(BND_D1)))
