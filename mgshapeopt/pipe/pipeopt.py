from pipeproblem import PipeProblem
from obstacleproblem import ObstacleProblem
import fireshape as fs
import fireshape.zoo as fsz
import firedrake as fd
from distance_function import distance_function
from alfi import get_default_parser, get_solver, run_solver
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
parser.add_argument("--surf", dest="surf", default=False,
                    action="store_true")
parser.add_argument("--spectral", dest="spectral", default=False,
                    action="store_true")
parser.add_argument("--smooth", dest="smooth", default=False,
                    action="store_true")
args, _ = parser.parse_known_args()
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
    if args.dim == 3:
        if h is None:
            h = 1./4.
    elif args.dim == 2:
        if h is None:
            h = 1.
    else:
        raise NotImplementedError
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
    mu_cr = 1 * mu
    extension[0] = fs.CauchyRiemannAugmentation(extension[0], mu=mu_cr)
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
    # innerp = fs.ElasticityInnerProduct(
    #     Q, fixed_bids=fixed_bids, direct_solve=True)

# import IPython; IPython.embed()

res = [1, 10, 50, 100, 150, 200, 250, 300, 400, 500, 750, 1000]
res = [r for r in res if r <= optre]
if res[-1] != optre:
    res.append(optre)
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

    def update(self, *args):
        sys.stdout.flush()
        sys.stderr.flush()
        if super().update(*args):
            for m_ in Q.mh_m:
                m_.coordinates.dat.data_ro
                m_._shared_data_cache["hierarchy_physical_node_locations"] = {}
            if hasattr(solver, 'vtransfer'):
                solver.vtransfer.force_rebuild()
            for i, mesh in enumerate(Q.mh_m):
                fd.File("output/mesh-%i.pvd" % i).write(mesh.coordinates)
            fd.warning(fd.BLUE % "Solve state")
            solver.solve(optre)
            fd.warning(fd.BLUE % "Solve adjoint")
            solver.solver_adjoint.solve()
            fd.warning(fd.RED % ("J(Omega)=%f" % fd.assemble(objective_form)))


constraint = Constraint()
J = Objective(Q)
if args.label is not None:
    out = fd.File("output/u-%s.pvd" % args.label)
else:
    out = fd.File("output/u.pvd")

qactionout = fd.File("output/actionq.pvd")
temp = fs.InteriorControlConstraint(Q.V_r_coarse, form=extension)


def cb(*args):
    sys.stdout.flush()
    sys.stderr.flush()
    out.write(solver.z.split()[0])
    temp.temp2.assign(Q.intermediate_Ts[0])
    temp.temp2 -= Q.intermediate_Ids[0]
    temp.apply_action(temp.temp2, temp.temp1)
    temp.bc.apply(temp.temp1)
    qactionout.write(temp.temp1)


J.cb = cb
if args.spectral:
    Js = fsz.MoYoSpectralConstraint(1e3, fd.Constant(0.5), Q)
    J = J + Js
if args.tikhonov > 0:
    Jt = args.tikhonov * fsz.DeformationRegularization(
        Q, l2_reg=0., sym_grad_reg=1., skew_grad_reg=1.)
    J = J + Jt

if args.smooth:
    # dirichlet_extension = fs.DirichletExtension(Q.V_r_coarse, form=extension)
    control_constraint = fs.InteriorControlConstraint(Q.V_r_coarse, form=extension)
else:
    dirichlet_extension = None
    control_constraint = None
q = fs.ControlVector(Q, innerp, control_constraint=control_constraint)
vol = fsz.LevelsetFunctional(fd.Constant(10.0), Q)
x, y = fd.SpatialCoordinate(Q.mesh_m)
baryx = fsz.LevelsetFunctional(x, Q)
baryy = fsz.LevelsetFunctional(y, Q)
if args.problem == "pipe":
    wrap = lambda f: fs.DeformationCheckObjective(f, delta_threshold=0.10,  # noqa
                                                  strict=False)
    scale = 1e-1
    J = wrap(scale*J)
    vol = wrap(0.1 * scale**0.5 * vol)
    econ = fs.EqualityConstraint([vol])
    emul = ROL.StdVector(1)
elif args.problem == "obstacle":
    if args.surf:
        scale = 1e-3
    else:
        scale = 1e0
    # all points on lower part of obstacle need to have y<=0
    Jbox1 = fsz.MoYoBoxConstraint(
        fd.Constant(1e5), [5], Q, lower_bound=fd.Constant((-100, -100)),
        upper_bound=fd.Constant((100, 0)))
    # all points on upper part of obstacle need to have y>=0
    Jbox2 = fsz.MoYoBoxConstraint(
        fd.Constant(1e5), [6], Q, lower_bound=fd.Constant((-100, 0)),
        upper_bound=fd.Constant((100, 100)))

    wrap = lambda f: fs.DeformationCheckObjective(f, delta_threshold=0.50,  # noqa
                                                  strict=False)
    # wrap = lambda f: f

    # J = wrap(scale*J + Jbox1 + Jbox2)
    J = wrap(scale*J)
    vol = wrap(scale**0.5 * vol)
    baryx = wrap(scale**0.5 * 1e1*baryx)
    baryy = wrap(scale**0.5 * 1e1*baryy)
    econ = fs.EqualityConstraint([vol, baryx, baryy])
    emul = ROL.StdVector(3)
else:
    raise NotImplementedError
params_dict = {
    'General': {
        'Secant': {
            'Type': 'Limited-Memory BFGS', 'Maximum Storage': 5
        }
    },
    'Step': {
        'Type': 'Augmented Lagrangian',
        'Line Search': {
            'Descent Method': {
                'Type': 'Quasi-Newton Step'
            }
        },
        'Augmented Lagrangian': {
            'Subproblem Step Type': 'Line Search',
            'Penalty Parameter Growth Factor': 1.5,
            'Print Intermediate Optimization History': True,
            'Subproblem Iteration Limit': 100,
            "Use Default Initial Penalty Parameter": False,
            "Initial Penalty Parameter": 1.0,
            "Use Default Problem Scaling": False,
            "Initial Optimality Tolerance": 1e-6,
            "Initial Feasibility Tolerance": 1e-6,
        }
    },
    'Status Test': {
        'Gradient Tolerance': 1e-8,
        'Step Tolerance': 1e-6,
        'Iteration Limit': 1
    }
}


g = q.clone()
J.update(q, None, 0)
J.gradient(g, q, None)
g.scale(1e-3)
res = J.checkGradient(q, g, 5, 1)
q.scale(0)
J.update(q, None, 0)

params = ROL.ParameterList(params_dict, "Parameters")
if args.problem == "obstacle":
    lb = q.clone()
    lb.fun.assign(-100)
    ub = q.clone()
    ub.fun.assign(100)
    x, y = fd.SpatialCoordinate(lb.fun.ufl_domain())
    fd.DirichletBC(
        lb.fun.function_space(), fd.as_vector([-100, -y]), 6).apply(lb.fun)
    fd.DirichletBC(
        ub.fun.function_space(), fd.as_vector([100, -y]), 5).apply(ub.fun)
    bnd = ROL.Bounds(lb, ub, 1.0)
    bnd = None
else:
    bnd = None
problem = ROL.OptimizationProblem(J, q, econ=econ, emul=emul, bnd=bnd)
rolsolver = ROL.OptimizationSolver(problem, params)
vol_before = vol.value(q, None)
rolsolver.solve(Q.mesh_m.mpi_comm().rank == 0)
fd.File("output/q-%s.pvd" % args.label).write(q.fun)
g = q.clone()
J.update(q, None, 0)
J.gradient(g, q, None)
g.scale(1e-2)
res = J.checkGradient(q, g, 7, 1)
