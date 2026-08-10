SetFactory("OpenCASCADE");
Mesh.Algorithm = 5;
W = 4;
D = 4;
R = 1;
h = 0.2;

Rectangle(1) = {-W/2, -D/2, 0, W, D/2, 0};
Rectangle(2) = {-W/2, 0, 0, W, D/2, 0};

Disk(10) = {0,0,0,R};
BooleanDifference(11) = {Surface{1}; Delete;} {Surface{10};};
BooleanDifference(12) = {Surface{2}; Delete;} {Surface{10}; Delete;};

Physical Surface (11) = {11};
Physical Surface (12) = {12};
Coherence;

Physical Line(20) = {8, 1};
Physical Line(30) = {10, 5};

Physical Line(40) = {3, 7};


Field[1]      =  Constant;
Field[1].VIn  =  h;
Field[1].VOut =  h;
Field[1].SurfacesList = {11, 12};

Background Field = 1;
