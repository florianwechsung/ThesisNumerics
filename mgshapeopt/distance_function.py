import firedrake as fd


def distance_function(mesh, boundary_ids="on_boundary", order=1,
                      eps=fd.Constant(0.1)):
    V = fd.FunctionSpace(mesh, "CG", order)
    u = fd.Function(V)
    v = fd.TestFunction(V)
    a = eps * fd.inner(fd.grad(u), fd.grad(v)) * fd.dx \
        + fd.inner(fd.grad(u), fd.grad(u)) * v * fd.dx \
        - v * fd.dx
    fd.solve(a == 0, u, bcs=fd.DirichletBC(V, 0, boundary_ids))
    return u

    def get_nullspace(self, V):
        return self.inner.get_nullspace(V)
