import netgen.gui
from ngsolve import *
import netgen.geom2d as geom2dim
import matplotlib.pyplot as plt
import numpy as np

go = 'y'
while(go == 'y'):
    # Setup the mesh
    default = input("Use default settings (h=0.1,k=3)? (y/n): ")
    if default == 'y' or temp == 'Y':
        h = 0.1e-1
        k = 3
    else:
        h = float(input("mesh_size -> h= "))
        k = int(input("order -> k= "))

    geo = geom2dim.SplineGeometry()

    geo.AddRectangle((0,0), (1,1), bcs=['b','r','t','l'], leftdomain=1)
    geo.AddRectangle((0.3,0.5), (0.5,0.7), bcs=['b1','r1','t1','l1'], leftdomain=2, rightdomain=1)
    geo.SetMaterial(1, "d2")
    geo.SetMaterial(2, "d1")
    mesh = Mesh(geo.GenerateMesh(maxh=h)) # standard mesh available
    print('#'*40)
    print("Mesh generated:")
    print('#'*40)
    print("Number of vertices: ", mesh.nv)
    print("Number of elements: ", mesh.ne)
    Draw(mesh)
    
    ####### Setup #######
    print('#'*40)
    print("Setting up domain...")
    print('#'*40)
    lams = {"d1" : 1, "d2" : 10}
    Lambda = CoefficientFunction([lams[mat] for mat in mesh.GetMaterials()])
    Cs = {"d1" : 1, "d2" : 0}
    C = CoefficientFunction([Cs[mat] for mat in mesh.GetMaterials()])

    print("materials:", mesh.GetMaterials())
    print("Lambda: ", Lambda)
    print("f: ", C)
    print('#'*40)
    print("Solving PDE...")
    print('#'*40)
    fes = H1(mesh, order=k, dirichlet='b|r|t|l')

    ## Set trial and test function
    u, v = fes.TnT()

    # Set right hand side
    f = LinearForm(fes)
    f += C*v*dx(definedon=mesh.Materials("d1"))

    # The bilinear form
    A = BilinearForm(fes, symmetric=True)
    A += Lambda*grad(u)*grad(v)*dx

    # Now assemble the system of equations
    with TaskManager():
        f.Assemble()
        A.Assemble()
        
    ###### Calculations #######
    # Calculate the solution field (function)
    gf = GridFunction(fes)
    gf.vec.data = A.mat.Inverse(fes.FreeDofs(), inverse='sparsecholesky') * f.vec


    Draw(gf, mesh, 'grid_function')

    temp = input("Calculate flux through borders of d1? (y/n): ")
    if temp == 'y' or temp == 'Y':
        print('#'*40)
        print("Calculating flux...")
        print('#'*40)
        # Calculate flux through border of d1
        GradX = GridFunction(fes)
        GradY = GridFunction(fes)
        Flux = GridFunction(fes)

        GradX.Set(Lambda*grad(gf)[0])
        GradY.Set(Lambda*grad(gf)[1])

        # Ω1 = (0.3,0.5)×(0.5,0.7)
        Wl = Integrate(IfPos(x-0.4, -1,1)*GradX, mesh, definedon=mesh.Boundaries('l1'))
        Wr = Integrate(IfPos(x-0.4, -1,1)*GradX, mesh, definedon=mesh.Boundaries('r1'))
        Wt = Integrate(IfPos(y-0.6, -1,1)*GradY, mesh, definedon=mesh.Boundaries('t1'))
        Wb = Integrate(IfPos(y-0.6, -1,1)*GradY, mesh, definedon=mesh.Boundaries('b1'))
        Wtot = Wl + Wr + Wl + Wb
        print("Total flux= ", Wtot, "(should be 0.04)")

        # Draw the latest solution and mesh
        Draw(GradX, mesh, 'GradX')
        Draw(GradY, mesh, 'GradY')
    
    go = input("Go again? (y/n): ")



