Mesh.Algorithm = 6;

L = 10.;
H = 1;
h = .2;

Point(1) = {-L/2, -H/2, 0, h};
Point(2) = { L/2, -H/2, 0, h};
Point(3) = { L/2,  H/2, 0, h};
Point(4) = {-L/2,  H/2, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Surface(1) = {1};
Physical Line(20) = {4};
Physical Line(30) = {2};
Physical Line(40) = {1};
Physical Line(50) = {3};