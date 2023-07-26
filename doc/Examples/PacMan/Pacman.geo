SetFactory("Built-in");
R=1;
Rf = .25;
w=Pi/8;
hf = 0.01;
hc = 0.1;


Point(1) = {0, 0, 0, hf};
Point(2) = {Rf,0,0,hf};
Point(3) = {R,0,0,hc};

Point(4) = {-Rf*Cos(w), -Rf*Sin(w), 0, hf};
Point(5) = {-Rf*Cos(w), Rf*Sin(w), 0, hf};
Point(6) = {-R*Cos(w), -R*Sin(w), 0, hc};
Point(7) = {-R*Cos(w), R*Sin(w), 0, hc};

Circle(1) = {4,1,2};
Circle(2) = {2,1,5};
Circle(3) = {6,1,3};
Circle(4) = {3,1,7};

Line(10) = {1,5};
Line(11) = {5,7};
Line(20) = {1,4};
Line(21) = {4,6};

Line Loop(1000) = {2,1,-10,20};
Plane Surface(1000) = {1000};
Line Loop(2000) = {3,4,-11,-2,-1,21};
Plane Surface(2000) = {2000};

Physical Surface(1) = {1000};
Physical Surface(2) = {2000};
Physical Line(30)   = {3,4};
Physical Line(40)   = {10,20};