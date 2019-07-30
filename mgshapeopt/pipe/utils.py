import fireshape as fs
import firedrake as fd


class C1Regulariser(fs.BilinearForm):

    def __init__(self, base_form, *args, mu=fd.Constant(1.0), **kwargs):
        super().__init__(*args, **kwargs)
        self.base_form = base_form
        self.mu = mu

    def get_form(self, V):
        u = fd.TrialFunction(V)
        v = fd.TestFunction(V)
        base = self.base_form.get_form(V)
        n = fd.FacetNormal(V.ufl_domain())

        def h(domain):
            if domain.topological_dimension() == 3:
                return (fd.FacetArea(domain)*6)**(1./3)
            else:
                return (fd.FacetArea(domain)*2)**(1./2)
        h = h(V.mesh())
        alpha = fd.Constant(50.)
        from firedrake import div, grad, dS, dx, inner, jump, avg

        def form(u, v):
            return inner(div(grad(u)), div(grad(v)))*dx \
                - inner(avg(div(grad(u))), jump(grad(v), n))*dS \
                - inner(jump(grad(u), n), avg(div(grad(v))))*dS \
                + alpha/h*inner(jump(grad(u), n), jump(grad(v), n))*dS

        return base + self.mu * (form(u[0], v[0]) + form(u[1], v[1]))

    def get_nullspace(self, V):
        return self.base_form.get_nullspace(V)

    def is_symmetric(self):
        return True
