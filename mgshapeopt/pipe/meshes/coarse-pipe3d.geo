SetFactory("OpenCASCADE");
a() = ShapeFromFile("/fs4/e590/e590/fwe590/ThesisNumerics/mgshapeopt/pipe/meshes/pipe3d.step");

Field[1] = MathEval;
Field[1].F = "0.100000";
Background Field = 1;
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Volume("Combined volume", 6) = {a()};
