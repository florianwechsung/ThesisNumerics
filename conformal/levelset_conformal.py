import firedrake as fd
import fireshape as fs
import fireshape.zoo as fsz
from cr_inner_product import CauchyRiemannAugmentation, distance_function
import ROL
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--use-cr", default=False, action="store_true")
parser.add_argument("--dirichlet", default=False, action="store_true")
parser.add_argument("--base-inner", type=str, default="elasticity",
                    choices=["elasticity", "laplace"])
parser.add_argument("--alpha", type=float, default=1e-2)
parser.add_argument("--clscale", type=float, default=0.1)
parser.add_argument("--maxiter", type=int, default=50)
args = parser.parse_args()


mesh = fs.DiskMesh(args.clscale, radius=3, smooth=10)

Q = fs.FeControlSpace(mesh)
d = distance_function(Q.get_space_for_inner()[0].mesh())
mu_base = 0.01 / (0.01 + d)


mu_base = fd.Constant(1.0)
if args.base_inner == "elasticity":
    inner = fs.ElasticityInnerProduct(Q, mu=mu_base, direct_solve=True)
elif args.base_inner == "laplace":
    inner = fs.LaplaceInnerProduct(Q, mu=mu_base, direct_solve=True)
else:
    raise NotImplementedError

if args.use_cr:
    mu_cr = mu_base/args.alpha
    inner = CauchyRiemannAugmentation(mu_cr, inner)


mesh_m = Q.mesh_m
(x, y) = fd.SpatialCoordinate(mesh_m)

a = 0.8
b = 2.0

expr = (fd.sqrt((x - a)**2 + b * y**2) - 1) \
    * (fd.sqrt((x + a)**2 + b * y**2) - 1) \
    * (fd.sqrt(b * x**2 + (y - a)**2) - 1) \
    * (fd.sqrt(b * x**2 + (y + a)**2) - 1) - 0.001

J = fsz.LevelsetFunctional(expr, Q, scale=0.1, quadrature_degree=5)
q = fs.ControlVector(Q, inner)

params_dict = {
    "General": {
        "Secant": {
            "Type": "Limited-Memory BFGS", "Maximum Storage": 1
        }
    },
    "Step": {
        "Type": "Line Search",
        "Line Search": {
            "Descent Method": {
                "Type": "Quasi-Newton Step"
            }
        }
    },
    "Status Test": {
        "Gradient Tolerance": 1e-7,
        "Step Tolerance": 1e-7,
        "Iteration Limit": args.maxiter,
    },
}

outdir = "output-%s/" % args.base_inner
out = fd.File(outdir + "domain.pvd")
def cb(*args):
    out.write(mesh_m.coordinates)
    # print(fd.assemble(abs(expr)*fd.ds))
J.cb = cb
params = ROL.ParameterList(params_dict, "Parameters")
problem = ROL.OptimizationProblem(J, q)
solver = ROL.OptimizationSolver(problem, params)
solver.solve(False)


def solve_somthing(mesh):
    V = fd.FunctionSpace(mesh, "CG", 2)
    u = fd.Function(V)
    v = fd.TestFunction(V)

    x, y = fd.SpatialCoordinate(mesh)
    # f = fd.sin(x) * fd.sin(y) + x**2 + y**2
    # uex = x**4 * y**4
    uex = x**3
    uex = uex - fd.assemble(uex * fd.dx)/fd.assemble(1 * fd.dx(domain=mesh))
    # f = fd.conditional(fd.ge(abs(x)-abs(y), 0), 1, 0)
    from firedrake import inner, grad, dx, ds, div
    f = uex - div(grad(uex))
    n = fd.FacetNormal(mesh)
    g = inner(grad(uex), n)
    # F = 0.1 * inner(u, v) * dx + inner(grad(u), grad(v)) * dx - f * v * dx - g * v * ds
    f = -div(grad(uex))
    F = inner(grad(u), grad(v)) * dx - f * v * dx

    bc = fd.DirichletBC(V, uex, "on_boundary")
    fd.solve(F == 0, u, bcs=bc, solver_parameters={"ksp_type": "cg", "pc_type": "jacobi", "snes_type": "ksponly", "ksp_converged_reason": None})
    print("||grad(u-uex)|| =", fd.norm(grad(u-uex)))
    err = fd.Function(fd.VectorFunctionSpace(mesh, "DG", V.ufl_element().degree())).interpolate(grad(u-uex))
    fd.File(outdir + "sln.pvd").write(u)
    fd.File(outdir + "err.pvd").write(err)

solve_somthing(mesh_m)
