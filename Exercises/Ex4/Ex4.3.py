# pylint: disable=unused-wildcard-import
import netgen.gui
from ngsolve import *
from netgen.geom2d import SplineGeometry

print("Welcome, Peter!")
### Definition of mesh ###
h = 0.1
geo = SplineGeometry()

geo.AddRectangle((0,0), (1,1), bcs=['b','r','t','l'], leftdomain=1, rightdomain=0)
geo.SetMaterial(1, "D")
mesh = Mesh(geo.GenerateMesh(maxh=h)) # standard mesh available
print('#'*40)
print("Mesh generated:")
print('#'*40)
print("Number of vertices: ", mesh.nv)
print("Number of elements: ", mesh.ne)
print("Dimension of mesh: ", mesh.dim)
Draw(mesh)

### Definition of FES ###


# define discrete conforming spaces of Sigma and V
# any idea why we need a different order here?
# try to use different combinations later as well

# a)
order = 3
t = 0.1
t2 = t*t
fc = 1
BND = 'b|r|t|l'
Vw = H1(mesh, order = order, dirichlet=BND)
Vb = H1(mesh, order = order, dim=2, dirichlet=BND)

#This defines a "compound FESpace"
V = FESpace([Vw, Vb])

# Trial and test functions
(w, beta) = V.TrialFunction()
(v, delta) = V.TestFunction()

# BLF
B = BilinearForm(V)
B += (grad(beta) * grad(delta) + 1/t2 * (grad(w) - beta)*(grad(v)-delta) )*dx
# LF
F = LinearForm(V)
F += fc*v*dx

# solution on V
sol = GridFunction(V, 'sol_mixed')
w_sol = sol.components[0]
beta_sol = sol.components[1]

### Solving the system ###
# Assemble the system of equations
with TaskManager(): # with Multithreading
    F.Assemble()
    B.Assemble()

# Solve the system
sol.vec.data = B.mat.Inverse(V.FreeDofs(), inverse = "umfpack") * F.vec
print("Done with solving!")

print("Trying to draw what you cobbled together here...")
Draw(w_sol, mesh, 'w')
Draw(beta_sol, mesh, 'beta')

### Calculating stuff in the system ###
# evaluate the flux using the sigma solution!
# Similar as for the H1, the var framework from the cont setting transfers to the discrete setting
# For sigma you can now evaluate sigma * n (since the normal trace is a cont operator for the H(div)

# Ω1 = (0.3,0.5)×(0.5,0.7)
BND_D1 = 'b1|r1|t1|l1' # boundary of Ω1
Wtot = Wl + Wr + Wl + Wb
Wtot2 = -Integrate(grad(b)-beta, mesh)
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
