from ngsolve import *
mesh = Mesh("quads.vol.gz")
V = H1(mesh, order=1)
bfi = SymbolicBFI(6*grad(V.TrialFunction())*grad(V.TestFunction()))
lfi = SymbolicLFI(16*V.TrialFunction())

#storage to gather the global matrix:
a = Matrix(V.ndof,V.ndof)
for i in range(V.ndof):
    for j in range(V.ndof):
        a[i,j] = 0.0;

#storage to gather the global vector:
f = Vector(V.ndof)
for i in range(V.ndof):
    f[i] = 0.0

for el in V.Elements():
    print("---------------------------------------- ")
    print("currently on element: ", el.nr)
    print("the vertices of this elements are:\n", el.vertices)
    print("the degrees of freedom of this element are:\n", el.dofs)
    
    dofs = el.dofs
    elmat = bfi.CalcElementMatrix(el.GetFE(),el.GetTrafo())
    print("the element matrix of this element is:\n", elmat)
    elvec = lfi.CalcElementVector(el.GetFE(),el.GetTrafo())
    print("the element vector of this element is:\n", elvec)
    
print("the final assembled matrix is:\n",a)
print("the final assembled vector is:\n",f)

u = GridFunction(V)
u.vec[0:8]=1
u.vec[8]=2.5

print(u.vec)

Draw(u,mesh,"u")
