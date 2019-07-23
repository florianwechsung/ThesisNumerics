from firedrake import UnitSquareMesh, SpatialCoordinate, VectorFunctionSpace, \
    Function, parameters, assemble, inner, grad, interpolate, dx, ds, sin, \
    cos, derivative, FunctionSpace, TestFunction, solve, adjoint, replace, \
    rhs, lhs, exp, ln, FacetNormal, DirichletBC
from math import log
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(linewidth=240)

"""
Example code to run taylor tests for shape derivatives for a range of
functionals with and without PDE constraint.
"""


class TaylorTest:

    """ Base class for Taylor tests """

    def __init__(self):
        N = 5
        self.mesh = UnitSquareMesh(N, N)
        self.T0 = self.mesh.coordinates.copy(deepcopy=True)
        # random perturmabion
        h = 1 / N
        self.X = SpatialCoordinate(self.mesh)
        W = VectorFunctionSpace(self.mesh, "CG", 1)
        self.w = Function(W)
        vec = self.w.vector()
        np.random.seed(1)
        vec.set_local(np.random.uniform(-h / 3., h /
                                        3., size=vec.get_local().shape))

        self.name = "missing name"
        self.bc = None

    def set_quadrature(self, degree):
        parameters["form_compiler"]["quadrature_degree"] = degree

    def eval_J(self):
        raise NotImplementedError

    def eval_dJdw(self):
        raise NotImplementedError

    def eval_ddJdw(self):
        raise NotImplementedError

    def perform_test(self):
        """
        #(name, mesh, T0, J0, w, dJdw, ddJdw=None, pre_assemble_callback=None):
        Evaluate J on the line T0+s*w, where T0 is the initial domain and w
        is a perturbation direction. Then, compare J(T0+s*w) with the Taylor
        expansion J0 + s*dJdw + 0.5*s*s*ddJdw and compute convergence rates.
        """

        print("".join(["#"] * 80))
        print("### " + self.name + " " +
              "".join(["#"] * (75 - len(self.name))))
        print("    err1 rate1    err2 rate2")
        self.mesh.coordinates.assign(self.T0)
        J0 = self.eval_J()
        dJdw = self.eval_dJdw()
        ddJdw = self.eval_ddJdw()
        steps = np.asarray([2 ** (-i) for i in range(1, 11)])
        err1 = np.zeros(steps.shape)
        err2 = np.zeros(steps.shape)
        for (i, s) in enumerate(steps):
            self.mesh.coordinates.assign(self.T0 + float(s) * self.w)

            if self.bc is not None:  # update DirichletBC
                self.bc.function_arg = self.bc._original_val

            # evaluate J on perturbed domain and its Taylor approximation
            Js = self.eval_J()
            pred1 = J0 + s * dJdw
            err1[i] = abs(Js - pred1)
            pred2 = pred1 + 0.5 * s * s * ddJdw
            err2[i] = abs(Js - pred2)
            if i is 0:
                rate1 = np.nan
                rate2 = np.nan
            else:
                rate1 = (log(err1[i]) - log(err1[i - 1])) / (
                    log(steps[i]) - log(steps[i - 1])
                )
                rate2 = (log(err2[i]) - log(err2[i - 1])) / (
                    log(steps[i]) - log(steps[i - 1])
                )
            print(
                "{0:1.2e}, {1:1.1f}, {2:1.1e}, {3:1.1f}".format(
                    err1[i], rate1, err2[i], rate2
                )
            )
        return steps, err1, err2


class Unconstrained(TaylorTest):
    """ Base class for functionals without a PDE constraint. """

    def eval_J(self):
        return assemble(self.J)

    def eval_dJdw(self):
        w = self.w
        J = self.J
        X = self.X
        return assemble(derivative(J, X, w))

    def eval_ddJdw(self):
        w = self.w
        J = self.J
        X = self.X
        return assemble(derivative(derivative(J, X, w), X, w))


class Constrained(TaylorTest):
    """ Base class for functionals with PDE constraint. """

    def __init__(self):
        super().__init__()

        self.params = {"ksp_type": "preonly", "pc_type": "lu"}
        V = FunctionSpace(self.mesh, "CG", 1)
        self.V = V
        u, v = Function(V), TestFunction(V)
        self.u = u
        self.v = v

    def eval_J(self):
        solve(self.F == 0, self.u, bcs=self.bc, solver_parameters=self.params)
        return assemble(self.J)

    def eval_dJdw(self):
        u = self.u
        v = self.v
        J = self.J
        F = self.F
        X = self.X
        w = self.w
        V = self.V
        params = self.params

        solve(self.F == 0, u, bcs=self.bc, solver_parameters=params)
        bil_form = adjoint(derivative(F, u))
        rhs = -derivative(J, u)
        u_adj = Function(V)
        solve(assemble(bil_form), u_adj, assemble(rhs), bcs=self.bc,
              solver_parameters=params)
        L = J + replace(self.F, {v: u_adj})
        self.L = L
        self.bil_form = bil_form
        return assemble(derivative(L, X, w))

    def eval_ddJdw(self):
        u = self.u
        v = self.v
        J = self.J
        F = self.F
        X = self.X
        w = self.w
        V = self.V
        L = self.L
        bil_form = self.bil_form
        params = self.params

        s = w
        y_s = Function(V)
        # follow p 65 of Hinze, Pinnau, Ulbrich, Ulbrich
        # Step 1:
        solve(
            assemble(derivative(F, u)),
            y_s,
            assemble(derivative(-F, X, s)),
            solver_parameters=params,
            bcs=self.bc,
        )
        # Step 2:
        Lyy_y_s = assemble(derivative(derivative(L, u), u, y_s))
        Lyu_s = assemble(derivative(derivative(L, u), X, s))

        h1 = Lyy_y_s
        h1 += Lyu_s

        Luy_y_s = assemble(derivative(derivative(L, X), u, y_s))
        Luu_s = assemble(derivative(derivative(L, X), X, s))
        h2 = Luy_y_s
        h2 += Luu_s
        h3_temp = Function(V)
        # Step 3:
        solve(assemble(bil_form), h3_temp, h1,
              bcs=self.bc, solver_parameters=params)
        F_h3_temp = replace(F, {v: h3_temp})
        h3 = assemble(derivative(-F_h3_temp, X))
        res = h2
        res += h3
        return res.vector().inner(w.vector())


class LevelsetExample1(Unconstrained):
    def __init__(self):
        super().__init__()
        X = self.X
        # simple case for plot
        self.name = "LevelSet Example 1"
        self.J = sin(X[0]) * cos(X[1]) * dx


class LevelsetExample2(Unconstrained):
    def __init__(self):
        super().__init__()
        X = self.X
        # volume and boundary integrals
        self.name = "LevelSet Example 2"
        self.J = (
            sin(X[0]) * cos(X[1]) * dx
            + pow(1.3 + X[0], 4.2) * pow(1.4 + X[1], 3.3) * dx
            + exp(sin(X[0]) + cos(X[1])) * dx
            + ln(5 + sin(X[0]) + cos(X[1])) * ds
        )
        self.set_quadrature(8)


class LevelsetExample3(Unconstrained):
    def __init__(self):
        super().__init__()
        X = self.X
        # volume and boundary integrals
        self.name = "LevelSet Example 3"
        V = FunctionSpace(self.mesh, "CG", 1)
        u = interpolate(sin(X[0]) * cos(X[1])**2, V)
        n = FacetNormal(self.mesh)
        self.J = (
            sin(X[0]) * cos(X[1]) * dx
            + pow(1.3 + X[0], 4.2) * pow(1.4 + X[1], 3.3) * dx
            + exp(sin(X[0]) + cos(X[1])) * dx
            + ln(5 + sin(X[0]) + cos(X[1])) * ds
            + inner(grad(u), n) * ds
        )
        self.set_quadrature(8)


class PdeConstraintExample1(Constrained):
    def __init__(self):
        super().__init__()
        u, v, X = self.u, self.v, self.X
        # simple case for plot
        self.name = "PDE constrained Example 1"
        f = X[1] * X[0]
        self.F = (u * v + inner(grad(u), grad(v)) - f * v) * dx
        self.J = u*dx


class PdeConstraintExample2(Constrained):
    def __init__(self):
        super().__init__()
        u, v, X = self.u, self.v, self.X
        # nonhomogeneous Neumann bc and nonlinear functional on bdry
        self.name = "PDE constrained Example 2"
        f = sin(X[1]) * cos(X[0])
        g = exp(f)
        self.F = (u * v + inner(grad(u), grad(v)) - f * v) * dx + g * v * ds
        self.J = u * dx + pow(1 + u * u, 2.5) * ds
        self.set_quadrature(10)


class NonlinearPdeConstraint(Constrained):
    def __init__(self):
        super().__init__()
        u, v, X = self.u, self.v, self.X
        self.name = "nonlinear PDE constraint"
        f = sin(X[1]) * cos(X[0])
        g = exp(f)
        self.F = (u * v + (1+u**2) * inner(grad(u),
                                           grad(v)) - f * v) * dx + g * v * ds
        self.J = u * dx + pow(1 + u * u, 2.5) * ds
        self.set_quadrature(10)


class DirichletBcConstraint(Constrained):
    def __init__(self):
        super().__init__()
        self.name = "PDE-constraint with DirBC"
        u, v, X = self.u, self.v, self.X
        f = sin(X[1]) * cos(X[0])
        g = f/2.
        self.g = g
        self.with_DirBC = 1  # special case for DirichletBC in TaylorTest
        self.bc = DirichletBC(self.V, 0., "on_boundary")
        self.F = (inner(grad(u+g), grad(v)) + (u+g)*v - f * v) * dx
        self.J = (u+g)*(u+g) * dx + pow(1 + u * u, 2.5) * ds
        self.set_quadrature(10)


if __name__ == "__main__":

    # Pick an Unconstrained example
    Unconstrained_t = LevelsetExample1
    # Unconstrained_t = LevelsetExample2
    # Unconstrained_t = LevelsetExample3

    level = Unconstrained_t()
    steps, err11, err12 = level.perform_test()

    # Pick a Constrained example
    Constrained_t = PdeConstraintExample1
    # Constrained_t = PdeConstraintExample2
    # Constrained_t = NonlinearPdeConstraint
    # Constrained_t = DirichletBcConstraint

    pde = Constrained_t()
    steps, err21, err22 = pde.perform_test()

    # create picture taylortest.pdf
    print("s,DeltaoneJone,DeltatwoJone,DeltaoneJtwo,DeltatwoJtwo\\\\")
    for i in range(len(steps)):
        print("%e,%e,%e,%e,%e\\\\" % (steps[i], err11[i], err12[i], err21[i], err22[i]))

    ys = [err11, err21, err12, err22]
    labels = [r"$|\delta_1(\rm J_1,s)|$", r"$|\delta_1(\rm J_2,s)|$",
              r"$|\delta_2(\rm J_2,s)|$", r"$|\delta_2(\rm J_2,s)|$"]
    markers = ["x", "o", "s", "d"]
    params = {'text.usetex': True,
              'font.size': 8.5,
              'font.family': 'serif',#Latin Modern Roman',
              'text.latex.unicode': True,
              'mathtext.fontset': 'cm',
              'legend.numpoints': 1,
              'legend.scatterpoints': 1,
              'legend.handlelength'  : 1.35,
              }
    plt.rcParams.update(params)
    #plt.figure(figsize=(14/2.54, 7/2.54))
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=.15, bottom=.2, right=.96, top=.88)
    #ax = plt.gca()
    y1 = steps**2
    y1 = y1 * 2*10**-3 / y1[0]
    plt.loglog(steps, y1, "-", label=r"$s^2$", linewidth=1.5)
    for i, y in enumerate(ys):
        if i==2:
            y2 = steps**3
            y2 = y2 * 3*10**-5 / y2[0]
            plt.loglog(steps, y2, ":", label=r"$s^3$", linewidth=1.5)
        plt.loglog(steps, y, label=labels[i], marker=markers[i],
                   markersize=2.5, markeredgewidth=1.5, linewidth=1.5)
    plt.legend(ncol=2, fontsize=5.5)#, bbox_to_anchor=(1, -0.3))
    plt.title(r"Taylor test", fontsize=8.5)
    plt.xlabel(r"$s$")
    #plt.ylabel(r"$\Delta_i(J, s)$")
    #plt.tight_layout()
    fig.set_size_inches(3.3, 1.85)
    plt.savefig("taylortest.pdf")
