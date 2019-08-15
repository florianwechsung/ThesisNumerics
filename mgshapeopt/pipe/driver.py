from pipeproblem import PipeProblem
from obstacleproblem import ObstacleProblem
import fireshape as fs
import fireshape.zoo as fsz
import firedrake as fd
import numpy as np
from distance_function import distance_function
from alfi import get_default_parser, get_solver, run_solver
from utils import C1Regulariser
import ROL
import sys

parser = get_default_parser()
parser.add_argument("--label", type=str, default=None)
parser.add_argument("--problem", type=str, default="obstacle")
parser.add_argument("--dim", type=int, default=2)
parser.add_argument("--order", type=int, default=1)
parser.add_argument("--opt-re", type=int, default=100)
parser.add_argument("--element-size", type=float, default=None)
parser.add_argument("--tikhonov", type=float, default=0)
parser.add_argument("--cr", type=float, default=0)
parser.add_argument("--htwo", type=float, default=0)
parser.add_argument("--surf", dest="surf", default=False,
                    action="store_true")
parser.add_argument("--spectral", dest="spectral", default=False,
                    action="store_true")
parser.add_argument("--smooth", dest="smooth", default=False,
                    action="store_true")
parser.add_argument("--lowstorage", dest="lowstorage", default=False,
                    action="store_true")

args, _ = parser.parse_known_args()
label = f"{args.problem}-{args.dim}d-{args.discretisation}-nref-{args.nref}-{args.solver_type}-{args.mh}-stab-{args.stabilisation_type}-stabw-{args.stabilisation_weight}-gamma-{args.gamma}-optre-{args.opt_re}-order-{args.order}-tikhonov-{args.tikhonov}-cr-{args.cr}-htwo-{args.htwo}-h-{args.element_size}"  # noqa
if args.label is not None:
    label = label + "-" + args.label
optre = args.opt_re
h = args.element_size
if args.problem == "pipe":
    if args.dim == 3:
        if h is None:
            h = 1./4.
    elif args.dim == 2:
        if h is None:
            h = 1./2.
    else:
        raise NotImplementedError
    problem = PipeProblem(order=args.order, dim=args.dim, element_size=h)
elif args.problem == "obstacle":
    if h is None:
        h = 1.
    problem = ObstacleProblem(order=args.order, dim=args.dim, element_size=h)
else:
    raise NotImplementedError
fixed_bids = problem.inflow_bids \
    + problem.noslip_fixed_bids \
    + problem.outflow_bids
Q = [None]
extension = [None]


def mh_cb(mh_r):
    fd.File("output/fine_init.pvd").write(mh_r[-1].coordinates)
    fd.File("output/coarse_init.pvd").write(mh_r[0].coordinates)
    # sys.exit()
    d = distance_function(mh_r[0])
    mu = 0.01/(d+0.01)
    extension[0] = fs.ElasticityForm(mu=mu)
    if args.cr > 0:
        mu_cr = args.cr * mu
        extension[0] = fs.CauchyRiemannAugmentation(extension[0], mu=mu_cr)
    if args.htwo > 0:
        mu_c1 = fd.Constant(args.htwo)
        extension[0] = C1Regulariser(extension[0], mu=mu_c1)
    if args.surf:
        Q[0] = fs.ScalarFeMultiGridControlSpace(
            mh_r, extension[0], order=args.order,
            fixed_bids=fixed_bids)
    else:
        Q[0] = fs.FeMultiGridControlSpace(mh_r, order=args.order)
    return Q[0].mh_m


solver = get_solver(args, problem, hierarchy_callback=mh_cb)
Q = Q[0]
extension = extension[0]
if args.surf:
    innerp = fs.ScalarSurfaceInnerProduct(
        Q, fixed_bids=fixed_bids)
else:
    innerp = fs.UflInnerProductFromForm(
        extension, Q, fixed_bids=fixed_bids, direct_solve=True)

# import IPython; IPython.embed()

res = [0, 1, 10, 50, 100, 150, 200, 250, 300, 400, 499, 625, 750, 875, 999]
res = [r for r in res if r <= optre-1]
if res[-1] != optre-1:
    res.append(optre-1)
results = run_solver(solver, res, args)
# import sys; sys.exit()

u, _ = fd.split(solver.z)
nu = solver.nu
objective_form = nu * fd.inner(fd.sym(fd.grad(u)), fd.sym(fd.grad(u))) * fd.dx
solver.setup_adjoint(objective_form)
solver.solver_adjoint.solve()

# import sys; sys.exit()


class Constraint(fs.PdeConstraint):

    def solve(self):
        super().solve()
        solver.solve(optre)


class Objective(fs.ShapeObjective):

    def __init__(self, *args_, **kwargs_):
        super().__init__(*args_, **kwargs_)
        self.fake_val = None

    def value(self, *args_):
        if self.fake_val is None:
            return super().value(*args_)
        return self.fake_val

    def value_form(self):
        return objective_form

    def derivative_form(self, w):
        if args.discretisation == "pkp0":
            fd.warning(fd.RED % "Using residual without grad-div")
            F = solver.F_nograddiv
        else:
            F = solver.F
        F = fd.replace(F, {F.arguments()[0]: solver.z_adj})
        L = F + self.value_form()
        X = fd.SpatialCoordinate(solver.z.ufl_domain())
        return fd.derivative(L, X, w)

    def derivative(self, out):
        super().derivative(out)
        if args.discretisation != "pkp0":
            return
        w = fd.TestFunction(self.V_m)
        u = solver.z.split()[0]
        v = solver.z_adj.split()[0]
        from firedrake import div, cell_avg, dx, tr, grad
        gamma = solver.gamma
        deriv = gamma * div(w) * cell_avg(div(u)) * div(v) * dx \
            + gamma * (cell_avg(div(u) * div(w) - tr(grad(u)*grad(w)))
                       - cell_avg(div(u)) * cell_avg(div(w))) * div(v) * dx \
            - gamma * cell_avg(div(u)) * tr(grad(v)*grad(w)) * dx
        fd.assemble(deriv, tensor=self.deriv_m,
                    form_compiler_parameters=self.params)
        outcopy = out.clone()
        outcopy.from_first_derivative(self.deriv_r)
        fd.warning(fd.RED % ("norm of extra term %e" % outcopy.norm()))
        out.plus(outcopy)

    def update(self, *args_):
        sys.stdout.flush()
        sys.stderr.flush()
        if super().update(*args_):
            for m_ in Q.mh_m:
                m_.coordinates.dat.data_ro
                m_._shared_data_cache["hierarchy_physical_node_locations"] = {}
            if hasattr(solver, 'vtransfer'):
                solver.vtransfer.force_rebuild()
            # for i, mesh in enumerate(Q.mh_m):
            #     if i == 0 or not args.lowstorage:
            #        fd.File("output/%s-mesh-%i.pvd" % (label, i)).write(mesh.coordinates)
            try:
                fd.warning(fd.BLUE % "Solve state")
                solver.solve(optre)
                fd.warning(fd.BLUE % "Solve adjoint")
                solver.solver_adjoint.solve()
                fd.warning(fd.BLUE % ("J(Omega)=%f" % fd.assemble(objective_form)))
            except:
                fd.warning(fd.RED % "Solver failed, let's try from scratch.")
                solver.z.assign(0)
                solver.z_adj.assign(0)
                res = list(range(0, optre+1, 25))
                run_solver(solver, res, args)


constraint = Constraint()
obj = Objective(Q)
J = obj
out = fd.File("output/%s.pvd" % label)


if args.spectral:
    Js = fsz.MoYoSpectralConstraint(1e3, fd.Constant(0.5), Q)
    J = J + Js
if args.tikhonov > 0:
    Jt = args.tikhonov * fsz.CoarseDeformationRegularization(extension, Q)
    J = J + Jt

if args.smooth:
    control_constraint = fs.InteriorControlConstraint(
        Q.V_r_coarse, form=extension)
else:
    dirichlet_extension = None
    control_constraint = None
q = fs.ControlVector(Q, innerp, control_constraint=control_constraint)
vol = fsz.LevelsetFunctional(fd.Constant(10.0), Q)
if args.problem == "pipe":
    econ_unscaled = fs.EqualityConstraint([vol])
    def wrap(f): return fs.DeformationCheckObjective(f, delta_threshold=0.25 if args.dim == 2 else 0.25,  # noqa
                                                     strict=False)
    scale = 1e1
    J = wrap(scale*J)
    volweight = 0.1 if args.dim == 2 else 1.
    vol = wrap(volweight * scale**0.5 * vol)
    econ = fs.EqualityConstraint([vol])
    emul = ROL.StdVector(1)
    econ_val = ROL.StdVector(1)
elif args.problem == "obstacle":
    if args.dim == 2:
        x, y = fd.SpatialCoordinate(Q.mesh_m)
    else:
        x, y, z = fd.SpatialCoordinate(Q.mesh_m)
    baryx = fsz.LevelsetFunctional(x, Q)
    baryy = fsz.LevelsetFunctional(y, Q)
    econ_unscaled = fs.EqualityConstraint([vol, baryx, baryy])
    if args.dim == 3:
        baryz = fsz.LevelsetFunctional(z, Q)
        econ_unscaled = fs.EqualityConstraint([vol, baryx, baryy, baryz])
    if args.surf:
        scale = 1e-3
    else:
        scale = 1e0

    def wrap(f): return fs.DeformationCheckObjective(f, delta_threshold=0.10,  # noqa
                                                  strict=False)
    J = wrap(scale*J)
    vol = wrap(scale**0.5 * vol)
    baryx = wrap(scale**0.5 * 1e1*baryx)
    baryy = wrap(scale**0.5 * 1e1*baryy)
    econ = fs.EqualityConstraint([vol, baryx, baryy])
    if args.dim == 3:
        baryz = wrap(scale**0.5 * 1e1*baryz)
        econ = fs.EqualityConstraint([vol, baryx, baryy, baryz])
    emul = ROL.StdVector(args.dim + 1)
    econ_val = ROL.StdVector(args.dim + 1)
else:
    raise NotImplementedError
params_dict = {
    'General': {
        'Secant': {
            'Type': 'Limited-Memory BFGS', 'Maximum Storage': 20
        }
    },
    'Step': {
        'Type': 'Augmented Lagrangian',
        'Line Search': {
            'Descent Method': {
                'Type': 'Quasi-Newton Step'
            },
        },
        'Augmented Lagrangian': {
            'Subproblem Step Type': 'Line Search',
            'Penalty Parameter Growth Factor': 4,
            'Print Intermediate Optimization History': True,
            'Subproblem Iteration Limit': 20,
            "Use Default Initial Penalty Parameter": False,
            "Initial Penalty Parameter": 1.0,
            "Use Default Problem Scaling": False,
            "Initial Optimality Tolerance": 1e-6,
            # "Initial Feasibility Tolerance": 1e-6,
        }
    },
    'Status Test': {
        'Gradient Tolerance': 1e-20,
        'Step Tolerance': 1e-20,
        'Iteration Limit': 12 if args.dim == 2 else 8,
    }
}


g = q.clone()
gecon = q.clone()
# J.update(q, None, 0)

# J.gradient(g, q, None)
# g.scale(1e-3)
# res = J.checkGradient(q, g, 5, 1)
# q.scale(0)
# J.update(q, None, 0)

params = ROL.ParameterList(params_dict, "Parameters")

problem = ROL.OptimizationProblem(J, q, econ=econ, emul=emul)
rolsolver = ROL.OptimizationSolver(problem, params)

data = {
    "iter": [],
    "state_snes_iters": [],
    "state_ksp_iters": [],
    "adjoint_snes_iters": [],
    "adjoint_ksp_iters": [],
    "drag": [],
    "Jval": [],
    "Jgrad": [],
    "cval": [],
}


def cb(*cbargs):
    sys.stdout.flush()
    sys.stderr.flush()
    if not args.lowstorage:
        out.write(solver.z.split()[0])

    state_snes_iters = solver.solver.snes.getIterationNumber()
    state_ksp_iters = solver.solver.snes.getLinearSolveIterations()
    state_ksp_iters *= 1/max(state_snes_iters, 1)
    adjoint_snes_iters = solver.solver_adjoint.snes.getIterationNumber()
    adjoint_ksp_iters = solver.solver_adjoint.snes.getLinearSolveIterations()
    adjoint_ksp_iters *= 1./max(adjoint_snes_iters, 1)

    if state_snes_iters == 0:
        return
    data["iter"].append(len(data["iter"]))
    data["state_snes_iters"].append(state_snes_iters)
    data["state_ksp_iters"].append(state_ksp_iters)
    data["adjoint_snes_iters"].append(adjoint_snes_iters)
    data["adjoint_ksp_iters"].append(adjoint_ksp_iters)

    econ_unscaled.value(econ_val, None, None)
    econ.applyAdjointJacobian(gecon, emul, None, None)
    J.gradient(g, None, None)
    g.plus(gecon)
    fd.warning("i: %d" % (len(data["iter"])-1))
    fd.warning("J: %e" % J.value(None, None))
    fd.warning("drag: %e" % obj.value(None, None))
    fd.warning("dL: %e" % g.norm())
    fd.warning("c : %e" % econ_val.norm())

    data["drag"].append(obj.value(None, None))
    data["Jval"].append(J.value(None, None))
    data["Jgrad"].append(g.norm())
    data["cval"].append(econ_val.norm())


obj.cb = cb

fd.warning(fd.BLUE % ("Initial volume: %f" % fd.assemble(fd.Constant(1) * fd.dx(domain=Q.mesh_m))))
rolsolver.solve(Q.mesh_m.mpi_comm().rank == 0)
fd.warning(fd.BLUE % ("Final volume: %f" % fd.assemble(fd.Constant(1) * fd.dx(domain=Q.mesh_m))))
if args.lowstorage:
    out.write(solver.z.split()[0])
fd.File("output/q-%s.pvd" % label).write(q.fun)
for k in data.keys():
    data[k] = np.asarray(data[k])
data["Jvalnormalised"] = data["Jval"]/data["Jval"][0]
data["Jgradnormalised"] = data["Jgrad"]/data["Jgrad"][0]
npdata = np.vstack([data[k] for k in data.keys()]).T
np.savetxt("output/" + label + ".txt", npdata, delimiter=";",
           header=";".join(data.keys()), comments='')
obj.cb = None
g = q.clone()
J.update(q, None, 0)
J.gradient(g, q, None)
g.scale(1e-2)
res = J.checkGradient(q, g, 7, 1)
