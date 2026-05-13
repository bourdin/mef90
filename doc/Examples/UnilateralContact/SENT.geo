SetFactory("OpenCASCADE");
Mesh.Algorithm = 5;
L = 1;
epsc = 0.01;

hf = 0.0025;  // Fine mesh size
hc = 0.1;    // Coarse mesh size

Rectangle(1) = {-L/2, -L/2, 0, L, L, 0};

Point(10) = {-L/2, -epsc/2, 0, hc};
Point(11) = {0, 0, 0, hf};
Point(12) = {-L/2, epsc/2, 0, hc};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 10};
Line Loop(10) = {10, 11, 12};
Plane Surface(10) = {10};

BooleanDifference(2) = {Surface{1}; Delete;} {Surface{10}; Delete;};

Field[1] = Box;
Field[1].VIn = hf;
Field[1].VOut = hc;
Field[1].XMin = -L/16;
Field[1].XMax = L;
Field[1].YMin = -L/2;
Field[1].YMax = L/2;
Field[1].Thickness = L/4;

Background Field = 1;

Physical Surface(1) = {2};
Physical Curve(20)  = {1};
Physical Curve(30)  = {4};
Physical Curve(40)  = {7};
Physical Curve(50)  = {5};
Physical Curve(60)  = {6};
Physical Curve(100) = {2, 3};