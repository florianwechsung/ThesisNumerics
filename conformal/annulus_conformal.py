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
parser.add_argument("--weighted", default=False, action="store_true")
parser.add_argument("--rstar", type=float, default=0.79)
args = parser.parse_args()


mesh = fd.Mesh("annulus.msh")
R = 1.0
r = 0.5
print("Harmonic map exists for r^*/R^* = %.2f" % ((0.5*(R/r + r/R))**-1))
Rs = 1.0
rs = args.rstar

Q = fs.FeControlSpace(mesh)
d = distance_function(Q.get_space_for_inner()[0].mesh(), boundary_ids=[1, 2])
if args.weighted:
    mu_base = 0.01 / (0.01 + d)
else:
    mu_base = fd.Constant(1.)



# mu_base = fd.Constant(1.0)
if args.base_inner == "elasticity":
    inner = fs.ElasticityInnerProduct(Q, mu=mu_base, direct_solve=True)
elif args.base_inner == "laplace":
    inner = fs.LaplaceInnerProduct(Q, mu=mu_base, direct_solve=True)
else:
    raise NotImplementedError

if args.alpha is not None:
    mu_cr = mu_base/args.alpha
    inner = CauchyRiemannAugmentation(mu_cr, inner)


mesh_m = Q.mesh_m
(x, y) = fd.SpatialCoordinate(mesh_m)

r = fd.sqrt(x**2 + y**2)
expr = (r-fd.Constant(rs))*(r-fd.Constant(Rs))
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
        "Gradient Tolerance": 1e-8,
        "Step Tolerance": 1e-8,
        "Iteration Limit": args.maxiter,
    },
}

outdir = "./output/annulus/base-%s-cr-%s-weighted-%s-rstar-%.2f/" % (args.base_inner, args.alpha, args.weighted, args.rstar)
out = fd.File(outdir + "domain.pvd")

params = ROL.ParameterList(params_dict, "Parameters")
problem = ROL.OptimizationProblem(J, q)
solver = ROL.OptimizationSolver(problem, params)
boundary_derivatives = []
gradient_norms = []
objective_values = []


def cb(*args):
    out.write(mesh_m.coordinates)
    gradient_norms.append(solver.getAlgorithmState().gnorm)
    objective_values.append(solver.getAlgorithmState().value)
    boundary_derivatives.append(fd.assemble(fd.inner(expr, expr) * fd.ds)**0.5)


J.cb = cb

solver.solve(True)

np.savetxt(outdir + "gradient_norms.txt", gradient_norms)
np.savetxt(outdir + "objective_values.txt", objective_values)
np.savetxt(outdir + "boundary_derivatives.txt", boundary_derivatives)
# import sys; sys.exit()

print(outdir)
