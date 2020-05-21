from ngsolve import *
from netgen.geom2d import SplineGeometry

geo = SplineGeometry()
p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [(0,0), (1,0), (1,1), (0,1)] ]
p5,p6 =  [ geo.AppendPoint(x,y) for x,y in [(2,0), (2,1)] ]

geo.Append (["line", p1, p2], leftdomain=1, rightdomain=0, bc='outer')
geo.Append (["line", p2, p3], leftdomain=1, rightdomain=2, bc='interface')
geo.Append (["line", p3, p4], leftdomain=1, rightdomain=0, bc='outer')
geo.Append (["line", p4, p1], leftdomain=1, rightdomain=0, bc='outer')
geo.Append (["line", p2, p5], leftdomain=2, rightdomain=0, bc='outer')
geo.Append (["line", p5, p6], leftdomain=2, rightdomain=0, bc='outer')
geo.Append (["line", p6, p3], leftdomain=2, rightdomain=0, bc='outer')

geo.SetMaterial(1, 'left')
geo.SetMaterial(2, 'right')

mesh = Mesh(geo.GenerateMesh(maxh=0.1))
Draw(mesh)

order = 3
Vl = H1(mesh, order = order, definedon = "left")
Vr = H1(mesh, order = order, definedon = "right")
Q = SurfaceL2(mesh, order = order-1)

X = FESpace([Vl,Vr,Q])

a = BilinearForm(X)

f = LinearForm(X)

gfu = GridFunction(X)

Draw (gfu.components[0], mesh, "ul")
Draw (gfu.components[1], mesh, "ur")


#to draw the global solution
V = L2(mesh, order = order)
gfu_total = GridFunction(V)
gfu_total.Set(gfu.components[0] + gfu.components[1])
Draw(gfu_total)
