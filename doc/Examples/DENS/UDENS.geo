SetFactory("Built-in");

H = 5;
W = 2;
a = .5;
Hf = 2.*a;
R = .1; // must be < a
hf = .025;
hc = .2;

Point(1) = {W/2,-H/2,0,hc};
Point(2) = {W/2,-Hf/2,0,hf};
Point(3) = {W/2,-R,0,hf};
Point(4) = {W/2-a+R,-R,0,hf};
Point(5) = {W/2-a,0,0,hf};
Point(6) = {W/2-a+R,R,0,hf};
Point(7) = {W/2,R,0,hf};
Point(8) = {W/2,Hf/2,0,hf};
Point(9) = {W/2,H/2,0,hc};

Point(10) = {-W/2,H/2,0,hc};
Point(11) = {-W/2,Hf/2,0,hf};
Point(12) = {-W/2,R,0,hf};
Point(13) = {-W/2+a-R,R,0,hf};
Point(14) = {-W/2+a,0,0,hf};
Point(15) = {-W/2+a-R,-R,0,hf};
Point(16) = {-W/2,-R,0,hf};
Point(17) = {-W/2,-Hf/2,0,hf};
Point(18) = {-W/2,-H/2,0,hc};

Point(200) = {-W/2+a-R,0,0,hf};
Point(201) = {W/2-a+R,0,0,hf};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Circle(4) = {4,201,5};
Circle(5) = {5, 201,6};
//Circle(4) = {W/2-a+R,0,0,R,Pi,3*Pi/2};
//Circle(5) = {W/2-a+R,0,0,R,Pi/2,Pi};

Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};

//Circle(13) = {-W/2+a-R,0,0,R,0,Pi/2};
//Circle(14) = {-W/2+a-R,0,0,R,-Pi/2,0};
Circle(13) = {14,200,13};
Circle(14) = {15,200,14};

Line(15) = {15,16};
Line(16) = {16,17};
Line(17) = {17,18};
Line(18) = {18,1};

Line(100) = {17,2};
Line(101) = {8,11};

//Coherence;

Line Loop(1000) = {1,-100,17,18};
Plane Surface(1000) = {1000};

Line Loop(2000) = {-101,8,9,10};
Plane Surface(2000) = {2000};

Line Loop(3000) = {100, 2, 3, 4, 5, 6, 7, 101, 11, 12, -13, -14, 15, 16};
Plane Surface(3000) = {3000};

Physical Surface(1) = {3000};
Physical Surface(2) = {1000,2000};

Physical Line(20)   = {9};
Physical Line(21)   = {18};
//
////Physical Line(20) = {101};
////Physical Line(21) = {100};
//
//
Physical Line(40)   = {3,4,5,6,12,13,14,15};
