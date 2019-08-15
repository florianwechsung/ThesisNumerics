import fireshape as fs
import fireshape.zoo as fsz
import firedrake as fd
import ROL

h = 0.75
nref = 1
order = 1
mh = fd.OpenCascadeMeshHierarchy(
    "meshes/disk.step", element_size=h,
    levels=nref, order=order, cache=True, verbose=True,
    project_refinements_to_cad=False, reorder=True,
    gmsh="/home/wechsung/bin/gmsh-4.4.0-Linux64/bin/gmsh -algo del2d -smooth 5"
)

fd.File("output/initial-coarse.pvd").write(mh[0].coordinates)
fd.File("output/initial-fine.pvd").write(mh[-1].coordinates)

Q = fs.ScalarFeMultiGridControlSpace(
    mh, fs.CauchyRiemannAugmentation(fs.ElasticityForm()), order=order)
inner = fs.ScalarSurfaceInnerProduct(Q)
# Q = fs.FeMultiGridControlSpace(mh, order=order)
# inner = fs.ElasticityInnerProduct(Q, direct_solve=True)

a = 0.8
b = 2.0

x, y = fd.SpatialCoordinate(Q.mesh_m)
expr = (fd.sqrt((x - a)**2 + b * y**2) - 1) \
    * (fd.sqrt((x + a)**2 + b * y**2) - 1) \
    * (fd.sqrt(b * x**2 + (y - a)**2) - 1) \
    * (fd.sqrt(b * x**2 + (y + a)**2) - 1) - 0.001

J = 0.001 * fsz.LevelsetFunctional(expr, Q, quadrature_degree=30)

q = fs.ControlVector(Q, inner)
g = q.clone()
J.update(q, None, 0)
J.gradient(g, q, None)
res = J.checkGradient(q, g, 10, 1)
import sys
sys.exit()

out = fd.File("output/domain.pvd")


def cb(*args):
    out.write(Q.mesh_m.coordinates)


J.cb = cb
J = fs.DeformationCheckObjective(J, delta_threshold=0.10, strict=False)

params_dict = {
    "General": {
        "Secant": {
            "Type": "Limited-Memory BFGS", "Maximum Storage": 10
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
        "Gradient Tolerance": 1e-9,
        "Step Tolerance": 1e-9,
        "Iteration Limit": 100,
    },
}

params = ROL.ParameterList(params_dict, "Parameters")
problem = ROL.OptimizationProblem(J, q)
solver = ROL.OptimizationSolver(problem, params)
solver.solve()
