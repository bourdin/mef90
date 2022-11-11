SetFactory("OpenCASCADE");
Mesh.Format=1;
Mesh.MshFileVersion=2;

L = 1;
h = 1;

Box(1) = {-L/2,-L/2,-L/2,L,L,L};

MeshSize{1:8} = h;
Physical Volume(1) = {1};
Physical Surface(20) = {1};
Physical Line(300)   = {5};
Physical Point(4000) = {2};
