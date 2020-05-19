from time import time_ns
from netgen.geom2d import unit_square
from ngsolve import *
from ngsolve import internal
import netgen.gui
#from IPython import display
import matplotlib.pyplot as plt
#%matplotlib inline
#internal.visoptions.autoscale = False

# mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
# #print(mesh.GetBoundaries())
# Draw(mesh)

# Some constants and functions
f = 1
ud = CoefficientFunction(x) 
g = CoefficientFunction(y)
alpha = 5

# Boundaries
Dirichlet = "bottom"
Neumann = "left"
Robin = "top|right"

# Runtimes
runtimes = {"primal": -1, "mixed": -1}

def primal(mesh, order=2, draw=False, rt=False):
    FESp = H1(mesh, order=order, dirichlet=Dirichlet)
    up, vp = FESp.TnT()

    # Forms
    ap = BilinearForm(FESp)
    #                           + Bilinear Robin-Term
    ap += grad(up)*grad(vp)*dx  + alpha*up*vp*ds(definedon=mesh.Boundaries(Robin))

    fp = LinearForm(FESp)
    fp += f*vp*dx + g*vp*ds(definedon=mesh.Boundaries(Neumann+"|"+Robin))

    with TaskManager():
        ap.Assemble()
        fp.Assemble()

    # Gridfunction
    gfup = GridFunction(FESp, "u_primal")
    gfup.Set(ud, BND) # BND =^ definedon=mesh.Boundaries("bottom|right|top|left")

    # Need to take of the Dirichlet BC
    with TaskManager():
        start_p = time_ns()
        # Homogenyze the system... "do-it-yoursolve" version
        # rhs = fp.vec.CreateVector()
        # rhs.data = fp.vec - ap.mat * gfup.vec
        # gfup.vec.data += ap.mat.Inverse(freedofs=FESp.FreeDofs()) * rhs
        solvers.BVP(bf=ap, lf=fp, gf=gfup)
        if rt:
            delta_p = (time_ns() - start_p) * 1e-9
            print(f"Runtime primal: {delta_p}s")
            runtimes["primal"] = round(delta_p, 3)

    flux = CoefficientFunction(grad(gfup))

    if draw:
        Draw (gfup)
        Draw (flux, mesh, "flux_primal")
        internal.visoptions.autoscale = False
    return gfup, flux

def mixed(mesh, order_flux=2, draw=False, rt=False):
    # Setup the mixed method
    V = HDiv(mesh, order=order_flux, dirichlet=Neumann) # "Neumann is the new Dirichlet" 
    Q = L2(mesh, order=order_flux-1)

    FESm = FESpace([V,Q]) # Compound space
    TnTm = FESm.TnT()
    (sigma, u), (tau, v) = TnTm
     # Setup the special function "normal" 
    normal = specialcf.normal(mesh.dim)

    # Set up forms again
    am = BilinearForm(FESm)
    #      a(sig,tau)+ b(sig,v)     + b(tau,u)       + bilinear Robin-Term
    am += (sigma*tau + div(sigma)*v + div(tau)*u)*dx + (1/alpha*(tau.Trace()*normal)*(sigma.Trace()*normal))*ds(definedon=mesh.Boundaries(Robin))
    am.Assemble()

    fm = LinearForm(FESm)
    #     source  + linear Dirichlet-term                                            + linear Robin-term
    fm += -f*v*dx + ud*(tau.Trace()*normal)*ds(definedon=mesh.Boundaries(Dirichlet)) + g/alpha*(tau.Trace()*normal)*ds(definedon=mesh.Boundaries(Robin)) 
    fm.Assemble()

    # Now on to the actual problem
    gfm = GridFunction(FESm, name="gf_mixed")
    gfsigma, gfum = gfm.components
    # Remember g from before? :D
    # It was defined as scalar valued --> "need to bring it into the vector world"
    gfsigma.Set(normal*g, BND) 
    with TaskManager():
        start_m = time_ns()
        # fm.vec.data -= am.mat * gfm.vec
        # gfm.vec.data += am.mat.Inverse(freedofs=FESm.FreeDofs(), inverse="umfpack") * fm.vec
        solvers.BVP(bf=am, lf=fm, gf=gfm)
        delta_m = (time_ns() - start_m) * 1e-9
        if rt:
            print(f"Runtime mixed: {delta_m}s")
            runtimes["mixed"] = round(delta_m, 3)

    if draw:
        Draw (gfsigma, mesh, "flux_mixed")
        Draw (gfum, mesh, "u_mixed")
        internal.visoptions.autoscale = False
    return gfm


hs = [0.5, 0.2, 0.1, 0.08, 0.05, 0.025, 0.01]
def compare_PvM(h_list=hs, order=2, ask=True):
    errors_flux = []
    errors_u = []
    for h, hnext in zip(h_list, h_list[1:]+["None"]):
        print("#"*40)
        mesh = Mesh(unit_square.GenerateMesh(maxh=h))
        print("h:", h)

        # Solve the problems
        #order = 1
        gfup, flux = primal(mesh, order=order, draw=False)
        gfm = mixed(mesh, order_flux=order+1, draw=False)
        gfsigma, gfum = gfm.components

        # Errors
        u_diff = gfup - gfum
        err_u = sqrt(Integrate( (u_diff)*(u_diff), mesh))
        errors_u.append(err_u)
        print ("err_u:", err_u)
        flux_diff = grad(gfup) - gfsigma
        err_flux = sqrt(Integrate(flux_diff*flux_diff, mesh))
        errors_flux.append(err_flux)
        print ("err_flux:", err_flux)
        
        Draw (gfup, mesh, "u_primal")
        Draw (flux, mesh, "flux_primal")
        Draw (gfum, mesh, "u_mixed")
        Draw (gfsigma, mesh, "flux_mixed")
        Draw (u_diff, mesh, "u_diff")
        Draw (flux_diff, mesh, "flux_diff")
        #internal.visoptions.linear = True
        #internal.visoptions.autoscale = False
        if ask:
            if not (hnext == "None"):
                print(f"Next h={hnext}")
            stop = input("Stop? (y/n): ")
            if stop == 'y':
                break
    return errors_u, errors_flux


hs = [0.5, 0.2, 0.1, 0.08, 0.05, 0.025, 0.01]

fig = plt.figure(figsize=(10,8))
orders = [1,2,3]
styles = ["x", "o", "s"]
for order, style in zip(orders, styles):
    print(f"Order= {order}")
    errors_u, errors_flux = compare_PvM(h_list=hs, order=order, ask=True)
    hs_used = hs[0:int(len(errors_u))]
    plt.semilogy(hs_used, errors_u, "-"+style, label=f"err(u)_o={order}")
    plt.semilogy(hs_used, errors_flux, "--"+style, label=f"err(flux)_o={order}")
plt.grid()
plt.legend()
pic = f"D:/Studium/Code/num-pde/Exercises/Ex9/Peter/Ex9.3_O=[{orders[0]}-{orders[-1]}].png"
plt.savefig(pic)
print(f"Saved plot at {pic}")
