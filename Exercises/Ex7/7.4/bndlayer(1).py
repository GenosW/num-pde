from ngsolve import *
import ngsolve.meshes as ngm

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp

ngsglobals.msg_level = 1

i = 5
mesh = ngm.MakeStructured2DMesh(quads=False, nx = 2**i ,ny = 2**i)
Draw(mesh)
V = H1(mesh, order = 1 ,dirichlet = [1,2,3,4])

IR = IntegrationRule(points=[(0,0), (1,0), (0,1)], weights=[1/6, 1/6, 1/6]) 

u = V.TrialFunction()
v = V.TestFunction()

lumping = 1

# setup solution vector
L = np.arange(0,1.005,0.005)
sol_line = np.zeros((7,len(L)))

for i in range(7):
    epsilon = 10**(-i)

    if lumping:
        f = LinearForm(V)
        f += v*dx(intrules = {TRIG: IR})

        A = BilinearForm(V)
        A += epsilon*grad(u)*grad(v)*dx(intrules = {TRIG: IR})
        A += u*v*dx(intrules = {TRIG: IR})

    else:
        
        f = LinearForm(V)
        f += v*dx

        A = BilinearForm(V)
        A += epsilon*grad(u)*grad(v)*dx
        A += u*v*dx

    f.Assemble()
    A.Assemble()

    gfu = GridFunction(V)

    gfu.vec.data = A.mat.Inverse(V.FreeDofs(), inverse="sparsecholesky") * f.vec

    Draw(gfu, mesh, "sol")

    sol_line[i] = gfu(mesh(L, 0.5))[:,0]



plt.plot(L, sol_line.T, '--')
plt.legend(range(7))
plt.show()


# check if M is diagonal
M = BilinearForm(V)
M += u*v*dx
M.Assemble()

rows,cols,vals = M.mat.COO()
M = sp.csr_matrix((vals,(rows,cols))).toarray()
non_diagonal_elements = np.count_nonzero(M - np.diag(np.diagonal(M)))

print("non diagonal elements: ", non_diagonal_elements)
