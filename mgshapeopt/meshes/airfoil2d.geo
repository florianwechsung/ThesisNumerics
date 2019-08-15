SetFactory("OpenCASCADE");
//Rectangle(1) = {-4, -4, 0, 10, 8, 0};
Rectangle(1) = {-1, -1, 0, 4, 2, 0};

Point(5) = {0, 0, 0,  0.1};
Point(6) = {0.5, 0.10, 0, 0.1};
Point(7) = {1.0, 0, 0, 0.1};
Point(8) = {0.5, -0.10, 0, 0.1};
Spline(5) = {5, 6, 7};
Spline(6) = {5, 8, 7};

Line Loop(2) = {5, -6};
Plane Surface(2) = {2};
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
