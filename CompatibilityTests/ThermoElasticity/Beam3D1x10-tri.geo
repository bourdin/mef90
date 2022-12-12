SetFactory("OpenCASCADE");
Mesh.Algorithm = 6;

L = 10.;
H = 1;
h = .2;

Box(1) = {-L/2,-H/2,-H/2,L,H,H};
Physical Volume(1) = {1};

Physical Surface(20) = 1;
Physical Surface(30) = 2;

Physical Surface(40) = 3;
Physical Surface(50) = 4;
Physical Surface(60) = 5;
Physical Surface(70) = 6;
