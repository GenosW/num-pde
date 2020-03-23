import netgen.gui
from ngsolve import *
import netgen.geom2d as geom2dim
import matplotlib.pyplot as plt
import numpy as np


# Setup the mesh
default = input("Use default settings (h=0.1,k=3)? (y/n): ")
if default == 'y' or default == 'Y':
    h = 1e-1
    k = 3
else:
    h = float(input("mesh_size -> h= "))
    k = int(input("order -> k= "))

geo = geom2dim.SplineGeometry()

geo.AddRectangle((0,0), (1,1), bcs=['b','r','t','l'], leftdomain=1)
geo.SetMaterial(1, "d1")
mesh = Mesh(geo.GenerateMesh(maxh=h)) # standard mesh available
print('#'*40)
print("Mesh generated:")
print('#'*40)
print("Number of vertices: ", mesh.nv)
print("Number of elements: ", mesh.ne)
print("Mesh size: ", h)
Draw(mesh)

####### Setup #######
print('#'*40)
print("Solving PDE...")
print('#'*40)
fes = H1(mesh, order=k, dirichlet='r|l')

## Set trial and test function
u, v = fes.TnT()

# Set right hand side
f = LinearForm(fes)
f += x*v*dx

# The bilinear form
A = BilinearForm(fes, symmetric=True)
A += grad(u)*grad(v)*dx

# Now assemble the system of equations
with TaskManager():
    f.Assemble()
    A.Assemble()
    
###### Calculations #######
# Calculate the solution field (function)
gf = GridFunction(fes)
gf.vec.data = A.mat.Inverse(fes.FreeDofs(), inverse='sparsecholesky') * f.vec

###### Draw final soltion ########
Draw(gf, mesh, 'grid_function')

temp = input("Calculate flux through borders? (y/n): ")
if temp == 'y' or temp == 'Y':
    print('#'*40)
    print("Calculating flux...")
    print('#'*40)
    
    # Option 1
    wr = GridFunction(fes)
    wl = GridFunction(fes)
    wr.Set(x)
    wl.Set(1-x)

    print("---- Option 1 ----")
    Wr = Integrate( grad(gf)*grad(wr) - x*wr, mesh)
    Wl = Integrate( grad(gf)*grad(wl) - x*wl, mesh)
    Wtot = Wl + Wr
    print("Wl= ", Wl)
    print("Wr= ", Wr)
    print("Total flux= ", Wtot)


    input("Go on?")
    # Option 2
    print()
    print("---- Option 1 ----")
    vector = gf.vec.CreateVector()
    vector.data = A.mat*vector - f.vec

    Wr2 = InnerProduct(vector, wr.vec)
    Wl2 = InnerProduct(vector, wl.vec)
    Wtot2 = Wl2 + Wr2
    print("Wl= ", Wl2)
    print("Wr= ", Wr2)
    print("Total flux= ", Wtot2)

    # Output to console
    print("2.) Total flux= ", Wtot2)
