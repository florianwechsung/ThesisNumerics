import firedrake as fd
import fireshape as fs
import fireshape.zoo as fsz
from cr_inner_product import CauchyRiemannAugmentation, distance_function
import ROL
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--use_cr", type=int, default=0)
parser.add_argument("--base_inner", type=str, default="elasticity")
args = parser.parse_args()

use_cr = bool(args.use_cr)
base_inner = args.base_inner

mesh = fd.Mesh("Sphere2D.msh")
Q = fs.FeControlSpace(mesh)
d = distance_function(Q.get_space_for_inner()[0].mesh(), eps=fd.Constant(0.1))
mu_base = 0.01 / (0.01 + d)

if base_inner == "elasticity":
    inner = fs.ElasticityInnerProduct(Q, fixed_bids=[1, 2, 3], mu=mu_base,
                                      direct_solve=True)
elif base_inner == "laplace":
    inner = fs.LaplaceInnerProduct(Q, fixed_bids=[1, 2, 3], mu=mu_base,
                                   direct_solve=True)
else:
    raise NotImplementedError

if use_cr:
    mu_cr = 100.0 * mu_base
    inner = CauchyRiemannAugmentation(mu_cr, inner)

mesh_m = Q.mesh_m
(x, y) = fd.SpatialCoordinate(mesh_m)
inflow_expr = fd.Constant((1.0, 0.0))
e = fsz.StokesSolver(
    mesh_m, inflow_bids=[1, 2, 3], inflow_expr=inflow_expr, noslip_bids=[4],
    direct=True, pin_pressure=True
)
# e.solve()
# fd.File("out.pvd").write(e.solution.split()[0])
# import sys; sys.exit()


directory = f"./output/stokes/base_{base_inner}_cr_{use_cr}/"
if not os.path.exists(directory):
    os.makedirs(directory, exist_ok=True)
out = fd.File(directory + "u.pvd")


vol = fsz.LevelsetFunctional(fd.Constant(1.0), Q)
baryx = fsz.LevelsetFunctional(x, Q)
baryy = fsz.LevelsetFunctional(y, Q)
econ = fs.EqualityConstraint([vol, baryx, baryy])
emul = ROL.StdVector(3)
econ_val = ROL.StdVector(3)

Je = fsz.EnergyObjective(e, Q, cb=None)
Jr = fs.ReducedObjective(Je, e)
# J = 2e-4 * Jr
J = 1e-2 * Jr
q = fs.ControlVector(Q, inner)


params_dict = {
    "General": {"Secant": {"Type": "Limited-Memory BFGS", "Maximum Storage": 5}},
    "Step": {
        "Type": "Augmented Lagrangian",
        "Line Search": {"Descent Method": {"Type": "Quasi-Newton Step"}},
        "Augmented Lagrangian": {
            "Subproblem Step Type": "Line Search",
            "Penalty Parameter Growth Factor": 10.0,
            # "Print Intermediate Optimization History": False,
            "Print Intermediate Optimization History": True,
            "Subproblem Iteration Limit": 20,
            "Use Default Initial Penalty Parameter": False,
            "Initial Penalty Parameter": 1.0,
            "Use Default Problem Scaling": False,
            # "Subproblem Iteration Limit": 100,
            # "Initial Optimality Tolerance": 1e-2,
            # "Initial Feasibility Tolerance": 1e-2,
        },
    },
    "Status Test": {
        "Gradient Tolerance": 1e-6,
        "Step Tolerance": 1e-6,
        "Contraint Tolerance": 1e-8,
        "Iteration Limit": 4,
    },
}
params = ROL.ParameterList(params_dict, "Parameters")
problem = ROL.OptimizationProblem(J, q, econ=econ, emul=emul)
solver = ROL.OptimizationSolver(problem, params)

g = q.clone()
gecon = q.clone()
J.update(q, None, 1)
# J.gradient(g, q, None)
# res = J.checkGradient(q, g, 8, 1)
# import sys; sys.exit()


Jvals = []
cnorms = []
gnorms = []
pde_solves = []

refinner = fs.LaplaceInnerProduct(Q, fixed_bids=[1, 2, 3], mu=mu_base,
                                  direct_solve=True)
refg = fs.ControlVector(Q, refinner)
refgecon = fs.ControlVector(Q, refinner)

fd.set_log_level(fd.INFO)


def cb(*args):
    econ.value(econ_val, None, None)
    econ.applyAdjointJacobian(refgecon, emul, None, None)
    J.gradient(refg, None, None)
    refg.plus(refgecon)
    Jvals.append(Je.value(None, None))
    cnorms.append(econ_val.norm())
    gnorms.append(refg.norm())
    pde_solves.append(e.num_solves)
    fd.info_blue(f"Lagrangian norm : {refg.norm()}")
    fd.info_blue(f"Constraint viol : {econ_val.norm()}")
    fd.info_blue(f"Lagrange multip : {[emul[i] for i in range(3)]}")
    out.write(e.solution.split()[0])
    # out_adj.write(e.solution_adj.split()[0])


# cb()
Jr.cb = cb
vol_before = vol.value(q, None)
solver.solve()
print(emul)
# print(Jvals)
np.save(directory + "Jvals", np.asarray(Jvals))
np.save(directory + "cnorms", np.asarray(cnorms))
np.save(directory + "gnorms", np.asarray(gnorms))
np.save(directory + "pde_solves", np.asarray(pde_solves))
# params_dict['Step']['Type'] = 'Composite Step'
# params_dict['Status Test']['Iteration Limit'] = 100
# print(vol.value(q, None)-vol_before)
# params = ROL.ParameterList(params_dict, "Parameters")
# problem = ROL.OptimizationProblem(J, q, econ=econ, emul=emul)
# solver = ROL.OptimizationSolver(problem, params)
# solver.solve()
print(vol.value(q, None) - vol_before)


e = fsz.StokesSolver(
    mesh_m, inflow_bids=[1, 2], inflow_expr=inflow_expr, noslip_bids=[4],
    direct=False
)
e.solve()
