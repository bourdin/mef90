L = 1.;
l = 0.25;
h = 0.025;

Point(1) = {-L/2., -L/2., 0, h};
Point(2) = { L/2., -L/2., 0, h};
Point(3) = { L/2.,  L/2., 0, h};
Point(4) = {-L/2.,  L/2., 0, h};
Point(5) = {-l/2, 0., 0., h};
Point(6) = { l/2, 0., 0., h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {1, 5};
Line(6) = {5, 6};
Line(60) = {6, 5};
Line(7) = {6, 3};

Line Loop(10) = {1, 2, -7, 60, -5};
Line Loop(11) = {5, 6, 7, 3, 4};
Plane Surface(10) = {10};
Plane Surface(11) = {11};

Physical Surface(1) = {10,11};
Physical Line(30) = {1};
Physical Line(40) = {-3};
Physical Line(50) = {6,60};
