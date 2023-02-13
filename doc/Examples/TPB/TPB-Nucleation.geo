mm = 1./90;
l = 90*mm;
w = 20*mm;
h = 8*mm;
L = 80*mm;
a = 4.5*mm;
ratio=3/4;
wf = a+(w-a)*ratio;
epsc = .001*mm;
theta = Atan(1);

//ell = 0.3*mm;
//hf = ell/6;
//hc = .5*mm;

hf = 1e-3;
hc = 2e-2;


Point(1) = {0, a, 0, hf};
Point(2) = {a*Sin(theta), 0, 0, hc};
Point(3) = {l/2, 0, 0, hc};
Point(4) = {l/2, w, 0, hc};
Point(5) = {-l/2, w, 0, hc};
Point(6) = {-l/2, 0, 0, hc};
Point(7) = {-a*Sin(theta), 0, 0, hc};

Point(8) = {a*Sin(theta)/2, a/2, 0, hf};
Point(9) = {a*Sin(theta)/2, wf, 0, hf};
Point(10) = {-a*Sin(theta)/2, wf, 0, hf};
Point(11) = {-a*Sin(theta)/2, a/2, 0, hf};


Point(12) = {-L/2, 0, 0, hc};
Point(13) = {0, w, 0, hc};
Point(14) = {L/2, 0, 0, hc};

Line(1) = {1, 8};
Line(2) = {8, 9};
Line(3) = {9, 10};
Line(4) = {10, 11};
Line(5) = {11, 1};
Line Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};

Line(6) = {8, 2};
Line(7) = {2, 14};
Line(8) = {14, 3};
Line(9) = {3, 4};
Line(10) = {4, 13};
Line(11) = {13, 5};
Line(12) = {5, 6};
Line(13) = {6, 12};
Line(14) = {12, 7};
Line(15) = {7, 11};
Line Loop(2) = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, -4, -3, -2};
Plane Surface(2) = {2};

Physical Surface(1) = {1};
Physical Surface(2) = {2};

Physical Line(30) = {5, 1};
Physical Point(400) = {13};
Physical Point(401) = {12};
Physical Point(402) = {14};
Physical Point(500) = {1};