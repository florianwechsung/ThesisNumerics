import firedrake as fd
import fireshape as fs
import fireshape.zoo as fsz
from cr_inner_product import CauchyRiemannAugmentation, distance_function
import ROL
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--base-inner", type=str, default="elasticity",
                    choices=["elasticity", "laplace"])
parser.add_argument("--alpha", type=float, default=None)
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

if args.alpha is not None:
    mu_cr = mu_base/args.alpha
    mu_div = mu_base/args.alpha
    mu_div = fd.Constant(0.)
    inner = CauchyRiemannAugmentation(mu_cr, mu_div, inner)


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

outdir = "./output/levelset/levelset-base-%s-cr-%s/" % (args.base_inner, args.alpha)
out = fd.File(outdir + "domain.pvd")

params = ROL.ParameterList(params_dict, "Parameters")
problem = ROL.OptimizationProblem(J, q)
solver = ROL.OptimizationSolver(problem, params)
boundary_derivatives = []
gradient_norms = []
objective_values = []

refinner = fs.ElasticityInnerProduct(Q, mu=fd.Constant(1.), direct_solve=True)
deriv = fs.ControlVector(Q, refinner)


def cb(*args):
    out.write(mesh_m.coordinates)
    J.derivative(deriv)
    deriv.apply_riesz_map()
    gradient_norms.append(deriv.norm())
    # gradient_norms.append(solver.getAlgorithmState().gnorm)
    objective_values.append(solver.getAlgorithmState().value)
    boundary_derivatives.append(fd.assemble(fd.inner(expr, expr) * fd.ds)**0.5)


J.cb = cb

solver.solve(False)

np.savetxt(outdir + "gradient_norms.txt", gradient_norms)
np.savetxt(outdir + "objective_values.txt", objective_values)
np.savetxt(outdir + "boundary_derivatives.txt", boundary_derivatives)
# import sys; sys.exit()

def solve_something(mesh):
    V = fd.FunctionSpace(mesh, "CG", 1)
    u = fd.Function(V)
    v = fd.TestFunction(V)

    x, y = fd.SpatialCoordinate(mesh)
    # f = fd.sin(x) * fd.sin(y) + x**2 + y**2
    # uex = x**4 * y**4
    uex = fd.sin(x)*fd.sin(y)#*(x*y)**3
    # def source(xs, ys):
    #     return 1/((x-xs)**2+(y-ys)**2 + 0.1)
    # uex = source(0, 0)
    uex = uex - fd.assemble(uex * fd.dx)/fd.assemble(1 * fd.dx(domain=mesh))
    # f = fd.conditional(fd.ge(abs(x)-abs(y), 0), 1, 0)
    from firedrake import inner, grad, dx, ds, div, sym
    eps = fd.Constant(0.0)
    f = uex - div(grad(uex)) + eps * div(grad(div(grad(uex))))
    n = fd.FacetNormal(mesh)
    g = inner(grad(uex), n)
    g1 = inner(grad(div(grad(uex))), n)
    g2 = div(grad(uex))
    # F = 0.1 * inner(u, v) * dx + inner(grad(u), grad(v)) * dx + inner(grad(grad(u)), grad(grad(v))) * dx - f * v * dx - g * v * ds
    F = inner(u, v) * dx + inner(grad(u), grad(v)) * dx - f * v * dx - g * v * ds
    F += eps * inner(div(grad(u)), div(grad(v))) * dx
    F += eps * g1 * v * ds
    F -= eps * g2 * inner(grad(v), n) * ds
    # f = -div(grad(uex))
    # F = inner(grad(u), grad(v)) * dx - f * v * dx

    # bc = fd.DirichletBC(V, uex, "on_boundary")
    bc = None
    fd.solve(F == 0, u, bcs=bc, solver_parameters={
        "ksp_type": "cg",
        "ksp_atol": 1e-13,
        "ksp_rtol": 1e-13,
        "ksp_dtol": 1e-13,
        "ksp_stol": 1e-13,
        "pc_type": "jacobi",
        "pc_factor_mat_solver_type": "mumps",
        "snes_type": "ksponly",
        "ksp_converged_reason": None
    })
    print("||u-uex||             =", fd.norm(u-uex))
    print("||grad(u-uex)||       =", fd.norm(grad(u-uex)))
    print("||grad(grad(u-uex))|| =", fd.norm(grad(grad(u-uex))))
    err = fd.Function(fd.TensorFunctionSpace(mesh, "DG", V.ufl_element().degree()-2)).interpolate(grad(grad(u-uex)))
    # err = fd.Function(fd.FunctionSpace(mesh, "DG", V.ufl_element().degree())).interpolate(u-uex)
    fd.File(outdir + "sln.pvd").write(u)
    fd.File(outdir + "err.pvd").write(err)

# solve_something(mesh_m)
print(outdir)
