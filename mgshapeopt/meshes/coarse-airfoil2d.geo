SetFactory("OpenCASCADE");
a() = ShapeFromFile("/home/wechsung/ThesisNumerics/mgshapeopt/meshes/airfoil2d.step");

Characteristic Length {5, 6} = 0.070000;
Characteristic Length {1, 2, 3, 4} = 1.400000;

Field[1] = MathEval;
Field[1].F = "0.007000 + 0.6*sqrt(x*x+y*y)";
Field[2] = MathEval;
Field[2].F = "0.007000 + 0.6*sqrt((x-1)*(x-1)+y*y)";
Field[3] = Min;
Field[3].FieldsList = {1, 2};
Background Field = 3;
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Surface(7) = {1};
