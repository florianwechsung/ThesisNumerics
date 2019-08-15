// Gmsh project created on Thu Apr 11 14:10:44 2019
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 2.25, 0, 2*Pi};
Line Loop(1) = {1};
Plane Surface(1) = {1};
