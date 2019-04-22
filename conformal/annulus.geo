Point(1) = {0, 0, 0, 0.1};

Point(2) = {0.5, 0, 0, 0.1};
Point(3) = {-0.5, 0, 0, 0.1};
Point(4) = {0, 0.5, 0, 0.1};
Point(5) = {0, -0.5, 0, 0.1};

Circle(1) = {4, 1, 2};
Circle(2) = {2, 1, 5};
Circle(3) = {5, 1, 3};
Circle(4) = {3, 1, 4};

Point(6) = {1, 0, 0, 0.1};
Point(7) = {-1, 0, 0, 0.1};
Point(8) = {0, 1, 0, 0.1};
Point(9) = {0, -1, 0, 0.1};

Circle(5) = {8, 1, 6};
Circle(6) = {6, 1, 9};
Circle(7) = {9, 1, 7};
Circle(8) = {7, 1, 8};
Line Loop(10) = {8, 5, 6, 7};
Line Loop(9) = {4, 1, 2, 3};
Plane Surface(11) = {10, 9};
Physical Line("Outer") = {8, 5, 6, 7};
Physical Line("Inner") = {4, 1, 2, 3};

Physical Surface("Annulus") = {11};
