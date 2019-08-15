from firedrake import *
from firedrake.petsc import PETSc
from alfi import *
from alfi.bary import BaryMeshHierarchy
import numpy as np

# //Field[1] = Box;
# //Field[1].VIn = %f;
# //Field[1].VOut = %f;
# //Field[1].XMax = 1.1;
# //Field[1].XMin = -0.1;
# //Field[1].YMax = 0.05;
# //Field[1].YMin = -0.05;
# //Field[1].Thickness = 1;
# //Field[2] = Box;
# //Field[2].XMin = -0.01;
# //Field[2].XMax = 0.04;
# //Field[2].YMin = -0.01;
# //Field[2].YMax = 0.01;
# //
# //Field[2].Thickness = 0.5;
# //Field[2].VIn = %f;
# //Field[2].VOut = %f;
# //Field[3] = Min;
# //Field[3].FieldsList = {1, 2};

class ObstacleProblem(NavierStokesProblem):
    def __init__(self, order=2, dim=2, element_size=1./4):
        super().__init__()
        if dim == 2:
            self.inflow_bids = [1, 3, 4]
            self.noslip_fixed_bids = []
            self.noslip_free_bids = [5, 6]
            self.outflow_bids = [2]
        elif dim == 3:
            self.inflow_bids = [1, 3]
            self.noslip_fixed_bids = []
            self.noslip_free_bids = [4]
            self.outflow_bids = [2]

        self.order = order
        self.dim = dim
        h = element_size
        if dim == 2:
            self.element_size = """
    Characteristic Length {5, 6} = %f;
    Characteristic Length {1, 2, 3, 4} = %f;

    Field[1] = MathEval;
    Field[1].F = "%f + 0.6*sqrt(x*x+y*y)";
    Field[2] = MathEval;
    Field[2].F = "%f + 0.6*sqrt((x-1)*(x-1)+y*y)";
    Field[3] = Min;
    Field[3].FieldsList = {1, 2};
    Background Field = 3;
    """ % (h/20, h, h/200, h/200)
        else:
            self.element_size = """
Field[1] = Box;
Field[1].Thickness = 1.0;
Field[1].VIn = %f;
Field[1].VOut = %f;
Field[1].XMax = 0.6;
Field[1].XMin = -0.6;
Field[1].YMax = 0.15;
Field[1].YMin = -0.15;
Field[1].ZMax = 0.15;
Field[1].ZMin = -0.15;
Field[2] = Box;
Field[2].Thickness = 0.5;
Field[2].VIn = %f;
Field[2].VOut = %f;
Field[2].XMin = -0.55;
Field[2].XMax = -0.45;
Field[2].YMin = -0.05;
Field[2].YMax = 0.05;
Field[2].ZMin = -0.05;
Field[2].ZMax = 0.05;
Field[3] = Box;
Field[3].Thickness = 0.5;
Field[3].VIn = %f;
Field[3].VOut = %f;
Field[3].XMin = 0.45;
Field[3].XMax = 0.55;
Field[3].YMin = -0.05;
Field[3].YMax = 0.05;
Field[3].ZMin = -0.05;
Field[3].ZMax = 0.05;
Field[4] = Min;
Field[4].FieldsList = {1, 2, 3};
Background Field = 4;
""" % (h/20, h, h/40, h, h/40, h)


    def mesh(self, distribution_parameters):
        base = Mesh("meshes/obstacle%id.msh" % self.dim, distribution_parameters=distribution_parameters)
        return base

    def mesh_hierarchy(self, hierarchy, nref, callbacks, distribution_parameters):

        if hierarchy == "uniform":
            mh_constructor = MeshHierarchy
        elif hierarchy == "bary":
            mh_constructor = BaryMeshHierarchy
        else:
            raise NotImplemtedError("Pipe problem currently just works for uniform refinements")

        if self.dim == 2:
            fi = "meshes/airfoil2d.step"
        else:
            fi = "meshes/obstacle3d.step"
        warning(RED % "Using cached mesh for reproducability!")
        mh = OpenCascadeMeshHierarchy(
            fi, element_size=self.element_size,
            levels=nref, order=self.order, cache=False, verbose=True,
            distribution_parameters=distribution_parameters,
            callbacks=callbacks, project_refinements_to_cad=False,
            reorder=True, mh_constructor=mh_constructor, cache=True,
            gmsh="gmsh -algo front%id -optimize_netgen 5 -smooth 5 -clscale 0.5" % self.dim
        )
        return mh

    def bcs(self, Z):
        X = SpatialCoordinate(Z.mesh())
        if self.dim == 2:
            uin = Constant((1, 0))
        else:
            uin = Constant((1, 0, 0))

        bcs = [DirichletBC(Z.sub(0), uin, self.inflow_bids),
               DirichletBC(Z.sub(0), 0., self.noslip_free_bids + self.noslip_fixed_bids)]
        return bcs

    def has_nullspace(self): return False

    def char_length(self): return 1.0

    def relaxation_direction(self): return "0+:1+"

    def mesh_size(self, u, domain_type):
        # the factor below are chosen so that the calculation is exact for a regular simplex
        if domain_type == "facet":
            if self.dim == 2:
                return FacetArea(u.ufl_domain())
            elif self.dim == 3:
                return (2 * FacetArea(u.ufl_domain())**(1./2))  # area to length
        elif domain_type == "cell":
            if self.dim == 2:
                return (2*CellVolume(u.ufl_domain())**(1./2))  # area to length
            elif self.dim == 3:
                return (6*CellVolume(u.ufl_domain())**(1./3))  # volume to length

if __name__ == "__main__":


    parser = get_default_parser()
    args, _ = parser.parse_known_args()
    problem = ObstacleProblem()
    solver = get_solver(args, problem)

    start = 250
    end = 10000
    step = 250
    res = [0, 1, 10, 100] + list(range(start, end+step, step))
    res = [1, 10, 100]# + list(range(start, end+step, step))
    results = run_solver(solver, res, args)
