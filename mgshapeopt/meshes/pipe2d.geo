// Gmsh project created on Tue Jan 22 11:40:52 2019
SetFactory("OpenCASCADE");

Point(5) =  {0,   0, 0.0, 1.0};
Point(6) =  {1.5, 0, 0.0, 1.0};
Point(7) =  {3.0, 0, 0.0, 1.0};
Point(8) =  {4.0, 0, 0.0, 1.0};
Point(9) =  {4.75, 0.5,0.0, 1.0};
Point(10) = {5.5,  5.0,0.0, 1.0};
Point(11) = {6.5,  5.0,0.0, 1.0};
Point(12) = {10.0, 5.0,0.0, 1.0};
Point(13) = {16.0, 5.0,0.0, 1.0};

Bezier(10) = {6, 7, 8, 9, 10, 11, 12};
Line(15) = {5, 6};
Line(16) = {12, 13};
Point(14) = {0, 0.5, 0.0, 1.0};
Point(15) = {0, -0.5, 0.0, 1.0};
Line(17) = {14, 15};
Wire(1) = {15, 10, 16};
Extrude { Curve{17}; } Using Wire {1}
Delete{ Curve{17, 16, 15, 10}; }
Field[1] = MathEval;
Field[1].F = "1.0";
Background Field = 1;

Physical Curve("Inflow") = {19};
Physical Curve("Outflow") = {27};
Physical Curve("WallFixed") = {18, 20, 25, 26};
Physical Curve("WallFree") = {22, 23};
Physical Surface("Pipe") = {1, 2, 3};
