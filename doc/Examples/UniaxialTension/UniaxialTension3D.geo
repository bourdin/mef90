L = 1;
H = 0.1;
h = 0.005;

Point(1) = {-L/2., -H/2., -H/2., h};
Point(2) = { L/2., -H/2., -H/2., h};
Point(3) = { L/2.,  H/2., -H/2., h};
Point(4) = {-L/2.,  H/2., -H/2., h};
Point(5) = {-L/2., -H/2.,  H/2., h};
Point(6) = { L/2., -H/2.,  H/2., h};
Point(7) = { L/2.,  H/2.,  H/2., h};
Point(8) = {-L/2.,  H/2.,  H/2., h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {1, 5};
Line(6) = {5, 6};
Line(7) = {6, 2};
Line(8) = {6, 7};
Line(9) = {7, 8};
Line(10) = {8, 5};
Line(11) = {8, 4};
Line(12) = {7, 3};

Line Loop(16) = {1, 2, 3, 4};
Plane Surface(16) = {16};
Line Loop(18) = {5, 6, 7, -1};
Plane Surface(18) = {18};
Line Loop(20) = {10, 6, 8, 9};
Plane Surface(20) = {20};
Line Loop(22) = {3, -11, -9, 12};
Plane Surface(22) = {22};
Line Loop(24) = {12, -2, -7, 8};
Plane Surface(24) = {24};
Line Loop(26) = {10, -5, -4, -11};
Plane Surface(26) = {26};
Surface Loop(28) = {20, 26, 18, 24, 22, 16};
Volume(28) = {28};

Physical Surface(20) = {26};
Physical Surface(30) = {24};
Physical Volume(1) = {28};
Physical Point(400) = {5};
Physical Point(401) = {6};
Physical Point(402) = {2};
