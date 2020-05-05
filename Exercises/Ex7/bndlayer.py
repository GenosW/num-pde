from ngsolve import *
import ngsolve.meshes as ngm

ngsglobals.msg_level = 1

i = 5
mesh = ngm.MakeStructured2DMesh(quads=False, nx = 2**i ,ny = 2**i)
Draw(mesh)
#V = H1(mesh, order = 1 ,dirichlet = [1,2,3,4])

#For point 3
#IR = IntegrationRule(points = [], weights = []) 

#a = BilinearForm(...)
#a += (... )*dx(intrules = {TRIG: IR})

