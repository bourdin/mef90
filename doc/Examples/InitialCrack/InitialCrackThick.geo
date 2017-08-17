L = 1.;
l = 0.25;
h = 0.025;

Point(1) = {-L/2., -L/2., 0, h};
Point(2) = { L/2., -L/2., 0, h};
Point(3) = { L/2.,  L/2., 0, h};
Point(4) = {-L/2.,  L/2., 0, h};
Point(5) = {-l/2., -h/4., 0, h};
Point(6) = { l/2., -h/4., 0, h};
Point(7) = { l/2.,  h/4., 0, h};
Point(8) = {-l/2.,  h/4., 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9) = {1, 5};
Line(10) = {7, 3};

Line Loop(10) = {1, 2, -10, -6, -5, -9};
Line Loop(11) = {9, -8, -7, 10, 3, 4};
Line Loop(12) = {-5, -6, -7, -8};

Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};

Physical Surface(1) = {10, 11};
Physical Surface(2) = {12};
Physical Line(30) = {1};
Physical Line(40) = {-3};
