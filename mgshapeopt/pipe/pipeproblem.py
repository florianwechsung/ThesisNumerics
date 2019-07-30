from firedrake import *
from firedrake.petsc import PETSc
from alfi import *
import numpy as np


class PipeProblem(NavierStokesProblem):
    def __init__(self, order=2, dim=2, element_size=1./4):
        super().__init__()
        if dim == 2:
            self.inflow_bids = [2]
            self.noslip_fixed_bids = [1, 3, 8, 9]
            self.noslip_free_bids = [5, 6]
            self.outflow_bids = [10]
            # self.inflow_bids = [1]
            # self.noslip_fixed_bids = [3]
            # self.noslip_free_bids = [4]
            # self.outflow_bids = [2]
        elif dim == 3:
            self.inflow_bids = [1]
            self.noslip_fixed_bids = [2, 4]
            self.noslip_free_bids = [3]
            self.outflow_bids = [5]
            # self.inflow_bids = [1]
            # self.noslip_fixed_bids = [3]
            # self.noslip_free_bids = [4]
            # self.outflow_bids = [2]

        self.order = order
        self.dim = dim
        # self.element_size = element_size
        self.element_size = """
Field[1] = MathEval;
Field[1].F = "%f";
Background Field = 1;
""" % element_size

    def mesh(self, distribution_parameters):
        base = Mesh("meshes/pipe%id.msh" % self.dim, distribution_parameters=distribution_parameters)
        return base

    def mesh_hierarchy(self, hierarchy, nref, callbacks, distribution_parameters):
        if hierarchy != "uniform":
            raise NotImplemtedError("Pipe problem currently just works for uniform refinements")

        mh = OpenCascadeMeshHierarchy(
            "meshes/pipe%id.step" % self.dim, element_size=self.element_size,
            levels=nref, order=self.order, cache=False, verbose=True,
            distribution_parameters=distribution_parameters,
            callbacks=callbacks, project_refinements_to_cad=False,
            reorder=True,
            gmsh="/home/wechsung/bin/gmsh-4.4.0-Linux64/bin/gmsh -algo del%id -optimize_netgen 5 -smooth 5 " % self.dim
        )
        return mh

    def bcs(self, Z):
        X = SpatialCoordinate(Z.mesh())
        if self.dim == 2:
            uin = 4 * as_vector([(1/2-X[1])*(1/2+X[1]), 0])
        else:
            rsq = X[1]**2 + X[2]**2
            uin = as_vector([1-4*rsq, 0, 0])

        bcs = [DirichletBC(Z.sub(0), uin, self.inflow_bids),
               DirichletBC(Z.sub(0), 0., self.noslip_free_bids + self.noslip_fixed_bids)]
        return bcs

    def has_nullspace(self): return False

    def char_length(self): return 1.0

    def relaxation_direction(self): return "0+:1+"

    def mesh_size(self, u, domain):
        # the factor below are chosen so that the calculation is exact for a regular simplex
        if domain == "facet":
            if self.dim == 2:
                return FacetArea(u.ufl_domain())
            elif self.dim == 3:
                return (2 * FacetArea(u.ufl_domain())**(1./2))  # area to length
        elif domain == "cell":
            if self.dim == 2:
                return (2*CellVolume(u.ufl_domain())**(1./2))  # area to length
            elif self.dim == 3:
                return (6*CellVolume(u.ufl_domain())**(1./3))  # volume to length


if __name__ == "__main__":


    parser = get_default_parser()
    args, _ = parser.parse_known_args()
    problem = PipeProblem()
    solver = get_solver(args, problem)

    start = 250
    end = 10000
    step = 250
    res = [0, 1, 10, 100] + list(range(start, end+step, step))
    # res = [1, 10, 100]# + list(range(start, end+step, step))
    results = run_solver(solver, res, args)
