#include <array>
#include <cmath>
#include "twovec.h"
#include "coordinate.h"
#include "shortestdistance.h"

std::array<double, 5> dist3D_Segment_to_Segment(Coordinate S1a, Coordinate S1b, Coordinate S2a, Coordinate S2b){
    double EPSILON = 0.00000001;
    TwoVec u = S1b - S1a; // Vector of parametric line
    TwoVec v = S2b - S2a;
    TwoVec w = S1a - S2a;
    double a = u*u; // For derivation, see http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment()
    double b = u*v;
    double c = v*v;
    double d = u*w;
    double e = v*w;
    double D = a*c - b*b;
    double sc, sN, sD = D;
    double tc, tN, tD = D;

    //Compute line parameters of two closest points
    if(D > EPSILON){ //Lines not parallel, most common case so in the if
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if(sN < 0.0){
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if(sN > sD){
            sN = sD;
            tN = e+b;
            tD = c;
        }
    }
    else{
        sN = 0.0; //Parallel, but you still need an answer, so force the first starting point as a closest point
        sD = 1.0;
        tN = e;
        tD = c;
    }
    if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d +  b);
            sD = a;
        }
    }
    sc = (abs(sN) < EPSILON ? 0.0 : sN / sD);
    tc = (abs(tN) < EPSILON ? 0.0 : tN / tD);
    TwoVec dP = w + (u*sc) - (v*tc);
    //Coordinate closestPoint1 = S1a + u*sc; keep these, for they are pretty useful parameters to know, but a downright disaster to return
    //Coordinate closestPoint2 = S2a + v*tc;
    double distance = sqrt(dP*dP);
    std::array<double, 5> ans = {distance, dP.x, dP.y, sc, tc};
    return ans;
}
