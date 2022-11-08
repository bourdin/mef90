SetFactory("OpenCASCADE");
L = 1;
h = 1;

Box(1) = {0,0,0,L,L,L};

Physical Volume(1) = {1};
Physical Surface(10) = {1};
Physical Surface(20) = {2};

Field[1]             = Constant;
Field[1].VIn         = h;
Field[1].VolumesList = {1};
Background Field = 1;