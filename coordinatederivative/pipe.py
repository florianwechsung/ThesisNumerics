from firedrake import *

mesh = Mesh("pipe.msh")
coords = mesh.coordinates.vector()
X = SpatialCoordinate(mesh)

W = mesh.coordinates.function_space()
gradJ = Function(W)
phi, psi = TrialFunction(W), TestFunction(W)
bc = DirichletBC(W, 0, [1, 2, 3])
A_riesz = assemble(inner(grad(phi), grad(psi)) * dx, bcs=bc)

Z = VectorFunctionSpace(mesh, "CG", 2) \
    * FunctionSpace(mesh, "CG", 1)
z, z_adjoint = Function(Z), Function(Z)
u, p = split(z)
test = TestFunction(Z)
v, q = split(test)

nu = 1./400.
e = 2*nu*inner(sym(grad(u)), sym(grad(v)))*dx - p*div(v)*dx \
    + inner(dot(grad(u), u), v)*dx + div(u)*q*dx
uin = 6 * as_vector([(1-X[1])*X[1], 0])
bcs = [DirichletBC(Z.sub(0), 0., [3, 4]),
       DirichletBC(Z.sub(0), uin, 1)]
sp = {"pc_type": "lu", "mat_type": "aij",
      "pc_factor_mat_solver_type": "mumps"}

J = nu * inner(sym(grad(u)), sym(grad(u))) * dx
volume = Constant(1.) * dx(domain=mesh)
target_volume = assemble(volume)
dvol = derivative(volume, X)
c = 1/20
L = replace(e, {test: z_adjoint}) + J
dL = derivative(L, X)

out = File("output/u.pvd")
def solve_state_and_adjoint():
    solve(e==0, z, bcs=bcs, solver_parameters=sp)
    solve(derivative(L, z)==0, z_adjoint,
          bcs=homogenize(bcs), solver_parameters=sp)
    out.write(z.split()[0])

solve_state_and_adjoint()
print("i;J;dJ\\\\")
for i in range(200):
    dJ = assemble(dL).vector() \
        + assemble(dvol).vector() * c * 2\
        * (assemble(volume)-target_volume)
    solve(A_riesz, gradJ, dJ)
    print("%i;%.6f;%.6f\\\\"
          % (i, assemble(J), norm(grad(gradJ))))
    coords -= gradJ.vector()
    solve_state_and_adjoint()

