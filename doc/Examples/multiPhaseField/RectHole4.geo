SetFactory("OpenCASCADE");
Mesh.Algorithm = 5;
W = 4;
D = 4;
R = 1;
h = 0.2;

Rectangle(1) = {-W/2, -D/2, 0, W/2, D/2, 0};
Rectangle(2) = {0, -D/2, 0, W/2, D/2, 0};
Rectangle(3) = {-W/2, 0, 0, W/2, D/2, 0};
Rectangle(4) = {0, 0, 0, W/2, D/2, 0};

Disk(10) = {0,0,0,R};
BooleanDifference(11) = {Surface{1}; Delete;} {Surface{10};};
BooleanDifference(12) = {Surface{2}; Delete;} {Surface{10};};
BooleanDifference(13) = {Surface{3}; Delete;} {Surface{10};};
BooleanDifference(14) = {Surface{4}; Delete;} {Surface{10}; Delete;};

Physical Surface (11) = {11};
Physical Surface (12) = {12};
Physical Surface (13) = {13};
Physical Surface (14) = {14};
Coherence;

Physical Line(20) = {1, 10};
Physical Line(30) = {8, 16};

Physical Line(40) = {3, 6, 13, 14};


Field[1]      =  Constant;
Field[1].VIn  =  h;
Field[1].VOut =  h;
Field[1].SurfacesList = {11, 12, 13, 14};

Background Field = 1;
