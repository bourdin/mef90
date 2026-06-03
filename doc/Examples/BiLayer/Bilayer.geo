SetFactory("OpenCASCADE");
Mesh.Algorithm = 5; //Delaunay
Lx  = 10;
Ly1 = 1;
Ly2 = .5;
lc  = .5;
wc  = 0.1;
h   =.05;

Rectangle(1) = {-Lx/2,0, 0, Lx, Ly1, 0};
Rectangle(2) = {-Lx/2,-Ly2, 0, Lx, Ly2, 0};

// Crack
Point(100) = {-wc/2, Ly1, 0, h};
Point(101) = {0, Ly1 - lc, 0, h};
Point(102) = {wc/2, Ly1, 0, h};
Line(100) = {100, 101};
Line(101) = {101, 102};
Line(102) = {102, 100};
Line Loop(100) = {100, 101, 102};
Plane Surface(100) = {100};

BooleanDifference {Surface{1}; Delete;} {Surface{100}; Delete;}

Coherence;

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Line(10) = {2, 5};
Physical Line(20) = {8};
Physical Line(30) = {3, 4};

Field[1] = Constant;
Field[1].IncludeBoundary = 1;
Field[1].VIn = h;
Field[1].VOut = h;
Field[1].SurfacesList = {1, 2};
Background Field = 1;