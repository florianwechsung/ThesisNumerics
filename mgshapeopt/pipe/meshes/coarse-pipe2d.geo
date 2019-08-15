SetFactory("OpenCASCADE");
a() = ShapeFromFile("/home/wechsung/ThesisNumerics/mgshapeopt/pipe/meshes/pipe2d.step");

Field[1] = Box;
Field[1].Thickness = 0.25;
Field[1].VIn = 0.031250;
Field[1].VOut = 0.125000;
Field[1].XMax = 1.75;
Field[1].XMin = 1.25;
Field[1].YMax = 10;
Field[1].YMin = -10;
Field[1].ZMax = 10;
Field[1].ZMin = -10;
Field[2] = Box;
Field[2].Thickness = 0.25;
Field[2].VIn = 0.031250;
Field[2].VOut = 0.125000;
Field[2].XMax = 10.25;
Field[2].XMin = 9.75;
Field[2].YMax = 10;
Field[2].YMin = -10;
Field[2].ZMax = 10;
Field[2].ZMin = -10;
Field[3] = Min;
Field[3].FieldsList = {1, 2};
Background Field = 3;
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Line(7) = {7};
Physical Line(8) = {8};
Physical Line(9) = {9};
Physical Line(10) = {10};
Physical Line(11) = {11};
Physical Line(12) = {12};
Physical Line(13) = {13};
Physical Line(14) = {14};
Physical Line(15) = {15};
Physical Line(16) = {16};
Physical Line(17) = {17};
Physical Line(18) = {18};
Physical Line(19) = {19};
Physical Line(20) = {20};
Physical Line(21) = {21};
Physical Line(22) = {22};
Physical Line(23) = {23};
Physical Line(24) = {24};
Physical Line(25) = {25};
Physical Line(26) = {26};
Physical Line(27) = {27};
Physical Line(28) = {28};
Physical Line(29) = {29};
Physical Line(30) = {30};
Physical Line(31) = {31};
Physical Line(32) = {32};
Physical Line(33) = {33};
Physical Line(34) = {34};
Physical Surface(35) = {1};
Physical Surface(36) = {2};
Physical Surface(37) = {3};
BooleanUnion{ Surface{1}; Delete;}{Surface{2}; Surface{3}; Delete;}