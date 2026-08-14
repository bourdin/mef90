SetFactory("OpenCASCADE");
Mesh.Algorithm = 5;

// Parameters, all overridable from the command line, i.e.
//     gmsh -2 surfing-Bimat.geo -setnumber theta 20 -o surfing-Bimat-20.msh
DefineConstant[
    Lx    = {10,  Name "Parameters/1 Lx"},
    Ly    = {5,   Name "Parameters/2 Ly"},
    lc    = {0.5, Name "Parameters/3 crack length"},  // lc < Lx/2 keeps the crack tip in the left subdomain
    epsc  = {.01, Name "Parameters/4 crack height"},  // height of the initial crack at its mouth
    h     = {.1,  Name "Parameters/5 mesh size"},
    theta = {0,   Name "Parameters/6 theta"}          // The angle between the interface and the y-axis, in degrees
];

th  = theta * Pi / 180.;
R   = 2 * (Lx + Ly); // any length large enough to cover the whole domain
tol = epsc / 10.;    // tolerance of the bounding box queries

// The domain
Rectangle(1) = {0, -Ly/2, 0, Lx, Ly, 0};

// The initial crack: a wedge of length lc and height epsc, with its mouth on
// the left edge of the domain and its tip at (lc, 0)
Point(11) = {0,  epsc/2, 0};
Point(12) = {0, -epsc/2, 0};
Point(13) = {lc, 0, 0};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 11};
Curve Loop(11) = {11, 12, 13};
Plane Surface(10) = {11};
BooleanDifference(11) = { Surface{1}; Delete; } { Surface{10}; Delete; };

// A half-plane whose left edge is the interface, i.e. the line through
// (Lx/2, 0) with direction (Sin(th), Cos(th)): a large rectangle lying in
// x > Lx/2, rotated by -th about (Lx/2, 0, 0)
Rectangle(2) = {Lx/2, -R, 0, R, 2*R, 0};
Rotate {{0, 0, 1}, {Lx/2, 0, 0}, -th} { Surface{2}; }

// Split the notched domain in two subdomains, on the left and on the right of
// the interface. Each of them consists of more than one piece if the crack cuts
// all the way through it, hence the lists
Sl() = BooleanDifference  { Surface{11};         } { Surface{2};         };
Sr() = BooleanIntersection{ Surface{11}; Delete; } { Surface{2}; Delete; };
nl = #Sl();

// Fragment the subdomains so that they share the interface curve. The fragments
// are returned in the order they were given, so the first nl ones make up the
// left subdomain
S() = BooleanFragments{ Surface{Sl(), Sr()}; Delete; }{};
Left()  = {};
Right() = {};
For i In {0 : #S()-1}
    If (i < nl)
        Left()  += S(i);
    Else
        Right() += S(i);
    EndIf
EndFor

// The interface is what the two subdomains have in common, the lips of the
// crack are the curves within the bounding box of the wedge
outer()  = Abs(CombinedBoundary{ Surface{S()}; });
iface()  = Abs(Boundary{ Surface{Left()}; });
iface() -= outer();
crack()  = Curve In BoundingBox{-tol, -epsc/2-tol, -tol, lc+tol, epsc/2+tol, tol};
crack() -= iface();
outer() -= crack();

Physical Surface(1) = {Left()};  // left of the interface
Physical Surface(2) = {Right()}; // right of the interface
Physical Curve(20)  = {outer()}; // exterior boundary, crack lips excluded
Physical Curve(30)  = {iface()}; // material interface
Physical Curve(40)  = {crack()}; // lips of the initial crack

MeshSize{ PointsOf{ Surface{S()}; } } = h;
