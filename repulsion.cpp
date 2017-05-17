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
#include "repulsion.h"

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
        Coordinate ahalfway = a0 + 0.5*(a1 - a0);
        Coordinate bhalfway = b0 + 0.5*(b1 - b0);
        return std::array<double, 5>{sqrt(whenOverlap*whenOverlap), ahalfway.x - bhalfway.x, ahalfway.y - bhalfway.y, 0.5, 0.5};
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

///Gets the average position of a particle
Coordinate getCenter(Particle &p){
    Coordinate center(0, 0);
    for(Coordinate i : p.positions){
        center.x += i.x;
        center.y += i.y;
    }
    center.x /= (npivot+2);
    center.y /= (npivot+2);
    return center;
}


//std::vector<int> chooseForRepulsion(std::vector<Particle> &plist, Particle &pCompare, int indexSelf){
//    Coordinate cComp = getCenter(pCompare);
//    std::vector<int> chosenIndices;
//    Coordinate ci;
//    for(int pin = 0; pin < plist.size(); ++pin){
//        ci = getCenter(plist[pin]);
//        if(ci.x < cComp.x + dRepulse && ci.x > cComp.x - dRepulse && ci.y < cComp.y + dRepulse && ci.y > cComp.y - dRepulse){
//            if(pin != indexSelf){
//                chosenIndices.push_back(pin);
//            }
//        }
//    }
//    return chosenIndices;
//}
//
/////Returns a vector of pointers to all elements that need to be considered for repulsion with particle pCompare
//std::vector<int> selectForRepulsion(std::vector<Particle> &plist, Particle &pCompare){
//    std::vector<int> chosenIndices;
//    for(int i = 0; i < plist.size(); ++i){
//        if(dist(getCenter(plist[i]), getCenter(pCompare)) < repulsionSelection && dist(getCenter(plist[i]), getCenter(pCompare)) > 0.000001){
//            chosenIndices.push_back(i);
//        }
//    }
//    return chosenIndices;
//}

std::vector<std::pair<int, int>> findIndicesRepulsion(std::vector<Particle> &plist){
    std::vector<Coordinate> centers;
    for(Particle &p : plist){
        centers.push_back(getCenter(p));
    }
    std::vector<std::pair<int, int>> ans;
    int csize = centers.size();
    for(int i = 0; i < csize; ++i){
        for(int j = i+1; j < csize; ++j){
            if(centers[i].x < centers[j].x + dRepulse && centers[i].x > centers[j].x - dRepulse && centers[i].y < centers[j].y + dRepulse && centers[i].y > centers[j].y - dRepulse){
                ans.push_back(std::pair<int, int>{i, j});
            }
        }
    }
    return ans;

}

void applyRepulsiveForceParticleToParticle(Particle &p1, Particle &p2){
    double d;
    double Ftot;
    TwoVec Fp;
    TwoVec Fq;
    std::array<double, 5> st; //Contains scalar distance, coordinates of vectorial distance and values of parameter for parametrization line
    for(int i = 0; i < npivot + 1; ++i){
        for(int j = 0; j < npivot + 1; ++j){
            st = shortest_distance(p1.positions[i], p1.positions[i+1], p2.positions[j], p2.positions[j + 1]);
            d = st[0];
            if(d < p1.D){
                Ftot = ko*(d - p1.D);
                Fp = TwoVec(st[1], st[2])*Ftot;
                Fq = -1*Fp;
                p1.forces[i] += Fq*(1 - st[3]);
                p1.pressures[i] += norm(Fq*(1 - st[3]))/p1.len;
                p1.forces[i + 1] += Fq*st[3];
                p1.pressures[i + 1] += norm(Fq*st[3])/p1.len;
                p2.forces[j] += Fp*(1 - st[4]);
                p2.pressures[j] += norm(Fp*(1 - st[4]))/p2.len;
                p2.forces[j + 1] += Fp*st[4];
                p2.pressures[j + 1] += norm(Fp*st[4])/p2.len;
            }
        }
    }
}

void repulsiveForce(std::vector<Particle> &plist){
    std::vector<std::pair<int, int>> indices = findIndicesRepulsion(plist);
    for(std::pair<int, int> indexPair : indices){
        applyRepulsiveForceParticleToParticle(plist[indexPair.first], plist[indexPair.second]);
    }
}

/////Calculates the repulsive forces for all points that conform to selectForRepulsion()
//void newerRepulsiveForce(std::vector<Particle> &plist){
//    double d;
//    double Ftot;
//    TwoVec Fp;
//    TwoVec Fq;
//    std::array<double, 5> st; //Contains scalar distance, coordinates of vectorial distance and values of parameter for parametrization line
//    std::vector<int> chosen;
//    int chosenIndex;
//    for(int p1 = 0; p1 < plist.size(); ++p1){
//        for(int i = 0; i < npivot + 1; ++i){
//            chosen = chooseForRepulsion(plist, plist[p1], p1);
//            for(int p2 = 0; p2 < chosen.size(); ++p2){
//                chosenIndex = chosen[p2];
//                for(int j = 0; j < npivot + 1; ++j){
//                    st = shortest_distance(plist[p1].positions[i], plist[p1].positions[i+1], plist[chosenIndex].positions[j], plist[chosenIndex].positions[j+1]);
//                    d = st[0];
//                    if(d < plist[p1].D){
//                        Ftot = ko*(d - plist[p1].D);
//                        Fp = TwoVec(st[1], st[2])*Ftot;
//                        Fq = Fp*-1;
//                        plist[p1].forces[i] += Fq*(1 - st[3]);
//                        plist[p1].pressures[i] += norm(Fq*(1 - st[3]))/plist[p1].len;
//                        plist[p1].forces[i + 1] += Fq*st[3];
//                        plist[p1].pressures[i + 1] += norm(Fq*st[3])/plist[p1].len;
//                        plist[chosenIndex].forces[j] += Fp*(1 - st[4]);
//                        plist[chosenIndex].pressures[j] += norm(Fp*(1 - st[4]))/plist[chosenIndex].len;
//                        plist[chosenIndex].forces[j + 1] += Fp*st[4];
//                        plist[chosenIndex].pressures[j + 1] += norm(Fp*st[4])/plist[chosenIndex].len;
//                    }
//                }
//            }
//        }
//    }
//}

//void oldRepulsiveForce(std::vector<Particle> &plist){
//    double d;
//    double Ftot;
//    TwoVec Fp;
//    TwoVec Fq;
//    std::array<double, 5> st; //Contains scalar distance, coordinates of vectorial distance and values of parameter for parametrization line
//    for(int p1 = 0; p1 < plist.size(); ++p1){
//        for(int i = 0; i < npivot + 1; ++i){ //Loop over segments of p1
//            for(int p2 = p1 + 1; p2 < plist.size(); ++p2){ //Loop over all particles
//                for (int j = 0; j < npivot + 1; ++j){ //Loop over segments of other particles
//                    st = shortest_distance(plist[p1].positions[i], plist[p1].positions[i+1], plist[p2].positions[j], plist[p2].positions[j+1]);
//                    d = st[0];
//                    if (d < plist[p1].D){ //If norm(d) smaller than diameter
//                        Ftot = ko*(d-plist[p1].D);
//                        Fp = TwoVec(st[1], st[2])*Ftot; //Force of a segment of particle p1 on p2
//                        Fq = Fp*-1;
//                        plist[p1].forces[i] += Fq*(1 - st[3]);
//                        plist[p1].forces[i + 1] += Fq*st[3];
//                        plist[p2].forces[j] += Fp*(1 - st[4]);
//                        plist[p2].forces[j + 1] += Fp*st[4];
//                    }
//                }
//            }
//        }
//    }
//}
