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

        dim = V.mesh().topological_dimension()
        def h(domain):
            if dim == 3:
                return (fd.FacetArea(domain)*2)**(1./2)
            else:
                return fd.FacetArea(domain)
        h = h(V.mesh())
        alpha = fd.Constant(10)
        from firedrake import div, grad, dS, dx, inner, jump, avg, ds

        def form(u, v):
            return inner(div(grad(u)), div(grad(v)))*dx \
                - inner(avg(div(grad(u))), jump(grad(v), n))*dS \
                - inner(jump(grad(u), n), avg(div(grad(v))))*dS \
                + alpha/h*inner(jump(grad(u), n), jump(grad(v), n))*dS \
                - inner(div(grad(u)), inner(grad(v), n))*ds \
                - inner(inner(grad(u), n), div(grad(v)))*ds \
                + alpha/h*inner(grad(u), grad(v))*ds

        form = base + self.mu * (sum(form(u[i], v[i]) for i in range(dim)))
        return form

    def get_nullspace(self, V):
        return self.base_form.get_nullspace(V)

    def is_symmetric(self):
        return True
