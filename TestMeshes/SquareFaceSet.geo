SetFactory("OpenCASCADE");
Mesh.Algorithm = 5;
//Mesh.Format=1;
//Mesh.MshFileVersion=2;

L = 1;
h = 0.5;

Point(1) = {0,0,0,h};
Point(2) = {1,0,0,h};
Point(3) = {1,1,0,h};
Point(4) = {0,1,0,h};

For i In {1:4}
        Line(i) = {i,i%4+1};
EndFor

Line Loop(1) = {1:4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};
Physical Line(20)   = {1};
Physical Line(21)   = {2};
Physical Point(300) = {2};
Physical Point(301) = {4};
