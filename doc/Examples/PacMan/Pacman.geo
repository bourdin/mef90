R = 1.;
w = .5;
hf = 0.005;
hc = .1;


Point(1) = {0, 0, 0, hf};
Point(2) = {-R*Cos(w), -R*Sin(w), 0, hc};
//Point(3) = {0, -R, 0, hc};
Point(4) = {R, 0, 0, hc};
//Point(5) = {0, R, 0, hc};
Point(6) = {-R*Cos(w), R*Sin(w), 0, hc};

Line(10) = {1, 2};
Circle(11) = {2, 1, 4};
Circle(12) = {4, 1, 6};
Line(13) = {6, 1};

Line Loop(100) = {10, 11, 12, 13};
Plane Surface(100) = {100};

Physical Surface(1) = {100};
Physical Line(30) = {10, 13};
Physical Line(40) = {11, 12};
