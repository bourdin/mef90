SetFactory("OpenCASCADE");
Mesh.Algorithm = 5;

m = 1;
t1 = 1*m;
t2 = 5*m;
w = 40*m;

hf = 0.02;
hc = 20*hf;


Rectangle(1) = {-w/2, 0, 0, w, t1, 0};
Rectangle(2) = {-w/2, t1, 0, w, t2, 0};

Coherence;

Field[1] = Box;
Field[1].VIn = hf;
Field[1].VOut = hc;
Field[1].XMin = -w/2;
Field[1].XMax = w/2;
Field[1].YMin = 0;
Field[1].YMax = t1;
Field[1].Thickness = t1/2;

Background Field = 1;

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Curve(30)  = {1};
Physical Curve(40)  = {6};
Physical Curve(50)  = {4, 7};
Physical Curve(60)  = {2, 5};
Physical Point(100) = {6};