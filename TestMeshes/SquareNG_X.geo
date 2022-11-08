SetFactory("OpenCASCADE");
L = 1;
h = 0.1;

Point(1) = {-L/2,-L/2,0,h};
Point(2) = { L/2,-L/2,0,h};
Point(3) = { L/2, L/2,0,h};
Point(4) = {-L/2, L/2,0,h};

For i In {1:4}
        Line(i) = {i,i%4+1};
EndFor

Line Loop(1)        = {1:4};
Plane Surface(1)    = {1};
Physical Surface(1) = {1};
Physical Line(10)   = {4};
Physical Line(20)   = {2};
Physical Line(30)   = {1,3};
Physical Point(100) = {1,4};
Physical Point(200) = {2,3};
