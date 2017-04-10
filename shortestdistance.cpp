#include "twovec.h" //Coordinate needs to know what a TwoVec is
#include "coordinate.h" //Particle needs to know what a Coordinate is
#include "global.h"
#include <iostream>
#include <list>
#include <cmath>
#include <vector>
#include <array>
#include <deque>
#include "particle.h"
#include "shortestdistance.h"


double determinant(TwoVec a, TwoVec b, TwoVec c){
    double x = b.y * c.z - c.y * b.z;
    double y = a.y * c.z - c.y * a.z;
    double z = a.y * b.z - b.y * a.z;
    return a.x*x - b.x*y + c.x*z;
}

//http://stackoverflow.com/questions/2824478/shortest-distance-between-two-line-segments#18994296
std::array<double, 5> shortest_distance(Coordinate a0, Coordinate a1, Coordinate b0, Coordinate b1){
    double EPSILON = 0.0001;
    std::array<Coordinate, 5> ans;

    TwoVec A = a1 - a0;
    TwoVec B = b1 - b0;
    double magA = sqrt(A*A);
    double magB = sqrt(B*B);
    A = A*(1/magA);
    B = B*(1/magB);
    TwoVec C = cross(A, B);
    double D = (sqrt(C.z*C.z))*(sqrt(C.z*C.z));
    // Case that lines are parallel
    if(fabs(D) < EPSILON){
        double d0 = A * (b0 - a0);
        double d1 = A * (b1 - a0);
        if(d0 <= 0 && 0 >= d1){
            if(fabs(d0) < fabs(d1)){
                return std::array<double, 5>{sqrt((a0 - b0)*(a0 - b0)), a0.x - b0.x, a0.y - b0.y, 0, 0};
            }
            return std::array<double, 5>{sqrt((a0 - b1)*(a0 - b1)), a0.x - b1.x, a0.y - b1.y, 0, 1};
        }
        else if(d0 >= magA && magA <= d1){
            if(fabs(d0) < fabs(d1)){
                return std::array<double, 5>{sqrt((a1 - b0)*(a1 - b0)), a1.x - b0.x, a1.y - b0.y, 1, 0};
            }
            return  std::array<double, 5>{sqrt((a1 - b1)*(a1 - b1)), a1.x - b1.x, a1.y - b1.y, 1, 1};
        }
        TwoVec whenOverlap = TwoVec((a0.x + A.x*d0) - b0.x, (a0.y + A.y*d0) - b0.y);
        return std::array<double, 5>{sqrt(whenOverlap*whenOverlap), 0, 0, -1, -1};
    }
    // Extended lines cross somewhere
    TwoVec t = b0 - a0;
    double detA = determinant(t, B, C);
    double detB = determinant(t, A, C);

    // These t's vary between 0 and magA or magB, not from 0 to 1;
    double t0 = detA/D;
    double t1 = detB/D;

    // Closest point on infinite line
    Coordinate pA = a0 + (A*t0);
    Coordinate pB = b0 + (B*t1);

    // Clamping, parameters not in [0, 1] but in [0, magA/B]
    if(t0 < 0) pA = a0;
    else if(t0 > magA) pA = a1;
    if(t1 <0) pB = b0;
    else if(t1 > magB) pB = b1;

    // Project A
    if(t0 < 0 || t0 > magB){
        double dot = B*(pA - b0);
        if(dot < 0) dot = 0;
        else if(dot > magB) dot = magB;
        pB = b0 + (B * dot);
    }

    // Project B
    if(t1 < 0 || t1 > magA){
        double dotdot = A*(pB - a0);
        if(dotdot < 0) dotdot = 0;
        else if(dotdot > magA) dotdot = magA;
        pA = a0 + (A*dotdot);
    }

    // Clamps coordinates, but not parameters
    if(t0 < 0) t0 = 0;
    if(t1 < 0) t1 = 0;
    if(t0 > magA) t0 = 1;
    if(t1 > magB) t1 = 1;

    return std::array<double, 5>{dist(pA, pB), pA.x - pB.x, pA.y - pB.y, t0/magA, t1/magB};

}

///Don't forget the next-nearest neighbor approach, it reduces to O(n^5/2)
void repulsiveForce(std::vector<Particle> &plist){
    double d;
    double Ftot;
    TwoVec Fp;
    TwoVec Fq;
    std::array<double, 5> st; //Contains scalar distance, coordinates of vectorial distance and values of parameter for parametrization line
    for(int p1 = 0; p1 < plist.size(); ++p1){
        for(int i = 0; i < npivot + 1; ++i){ //Loop over segments of p1
            for(int p2 = p1 + 1; p2 < plist.size(); ++p2){ //Loop over all particles
                for (int j = 0; j < npivot + 1; ++j){ //Loop over segments of other particles
                    st = shortest_distance(plist[p1].positions[i], plist[p1].positions[i+1], plist[p2].positions[j], plist[p2].positions[j+1]);
                    d = st[0];
                    if (d < plist[p1].D){ //If norm(d) smaller than diameter
                        Ftot = ko*(d-plist[p1].D);
                        Fp = TwoVec(st[1], st[2])*Ftot; //Force of a segment of particle p1 on p2
                        Fq = Fp*-1;
                        plist[p1].forces[i] += Fq*(1 - st[3]);
                        plist[p1].forces[i + 1] += Fq*st[3];
                        plist[p2].forces[j] += Fp*(1 - st[4]);
                        plist[p2].forces[j + 1] += Fp*st[4];
                    }
                }
            }
        }
    }
}
