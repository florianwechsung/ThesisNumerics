SetFactory("OpenCASCADE");
//Rectangle(1) = {-4, -4, 0, 10, 8, 0};
//Rectangle(1) = {-1, -1, 0, 4, 2, 0};

//Point(1) = {-4, 0, 0, 1};
//Point(2) = {-4, 4, 0, 1};
//Point(3) = {6, 4, 0, 1};
//Point(4) = {6, 0, 0, 1};
//
//Line(1) = {1, 2};
//Line(2) = {2, 3};
//Line(3) = {3, 4};
//
//
//Point(5) = {0, 0, 0,  0.1};
//Point(6) = {0.5, 0.10, 0, 0.1};
//Point(7) = {1.0, 0, 0, 0.1};
//Point(8) = {0.5, -0.10, 0, 0.1};
//Spline(5) = {5, 6, 7};
//
//Extrude {{1, 0, 0}, {0, 0, 0}, 2*Pi} {
//  Line{1}; Line{2}; Line{3}; Line{5}; 
//}
//Surface Loop(1) = {2, 1, 3};
//Surface Loop(2) = {4};
//Volume(1) = {1};
//Volume(2) = {2};
//BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
Cylinder(1) = {-4, 0, 0, 8, 0, 0, 2, 2*Pi};
Sphere(2) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
Dilate {{0, 0, 0}, {1, 0.25, 0.25}} {
  Volume{2}; 
}
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
Field[1] = Box;
Field[1].Thickness = 1.0;
Field[1].VIn = 0.05;
Field[1].VOut = 1;
Field[1].XMax = 0.6;
Field[1].XMin = -0.6;
Field[1].YMax = 0.15;
Field[1].YMin = -0.15;
Field[1].ZMax = 0.15;
Field[1].ZMin = -0.15;
Field[2] = Box;
Field[2].Thickness = 0.5;
Field[2].VIn = 0.025;
Field[2].VOut = 1;
Field[2].XMin = -0.55;
Field[2].XMax = -0.45;
Field[2].YMin = -0.05;
Field[2].YMax = 0.05;
Field[2].ZMin = -0.05;
Field[2].ZMax = 0.05;
Field[3] = Box;
Field[3].Thickness = 0.5;
Field[3].VIn = 0.025;
Field[3].VOut = 1;
Field[3].XMin = 0.45;
Field[3].XMax = 0.55;
Field[3].YMin = -0.05;
Field[3].YMax = 0.05;
Field[3].ZMin = -0.05;
Field[3].ZMax = 0.05;
Field[4] = Min;
Field[4].FieldsList = {1, 2, 3};
Background Field = 4;
