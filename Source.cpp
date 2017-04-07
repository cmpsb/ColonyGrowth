#include <iostream>
#include <cmath>
#include <cstdlib>
#include <boost/random.hpp>
#include <ctime>
#include <vector>
#include <string>
#include "print.h"
#include <ctime> // time_t
#include <cstdio>
#include <array>
#include <fstream>
#include <tuple>
#include "twovec.h" //Coordinate needs to know what a TwoVec is
#include "coordinate.h" //Particle needs to know what a Coordinate is
#include "global.h" //Particle needs to access the global variables
#include "shortestdistance.h"
#include "particle.h" // All clean and tidy in its own file
#include "division.h"


	//Random number generator
boost::mt19937 generator(time(0)); //number generator from random list
boost::normal_distribution<> //setup distributions
normalDistGrowth(0.0, growthRateDev), //growth distribution added to daughter cells, sigma = 0.277
normalDistAngle(0.0, 0), //Noise in orientation for daughter cells
normalDistLength(maxLength, maxLengthDev); //distribution of maximum length that particles can reach, 4.54, sigma = 0.46
boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
randomMu(generator, normalDistGrowth), //generates deviation in growth rate
randomTheta(generator, normalDistAngle), //generates noise in orientation
randomLmax(generator, normalDistLength); //generates division lengths

///Don't forget the next-nearest neighbor approach, it reduces to O(n^5/2)
void repulsiveForce(vector<Particle> &plist){
    double d;
    double Ftot;
    TwoVec Fp;
    TwoVec Fq;
    std::array<double, 5> st; //Contains scalar distance, coordinates of vectorial distance and values of parameter for parametrization line
    for(int p1 = 0; p1 < plist.size(); ++p1){
        for(int i = 0; i < npivot + 1; ++i){ //Loop over segments of p1
            for(int p2 = p1 + 1; p2 < plist.size(); ++p2){ //Loop over all particles
                for (int j = 0; j < npivot + 1; ++j){ //Loop over segments of other particles
                    st = dist3D_Segment_to_Segment(plist[p1].positions[i], plist[p1].positions[i+1], plist[p2].positions[j], plist[p2].positions[j+1]);
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



///Master function growAll
void growAll(vector<Particle> &p){
    int lengthBefore = p.size(); //Since a particle cannot divide twice during the same growAll() call, only loop over the initial amount of particles
    for(int i = 0; i < lengthBefore; i++){ //Unable to be range based since the ranging object changes length, leading to an infinite loop
        p[i].grow(); //Why no straight lines without noise: double uncertainty for pi
        if(p[i].L > p[i].Lmax){
            Particle pnew = Particle(0, 0, 0, 0, diameter, 0); //Necessary copy of all old parameters
            p.push_back(pnew); //Add new particle to p
            divide(p[i], p[p.size() - 1]); //Set new properties of daughter particles
            p[p.size()-1].ID = p.size() - 1; //Set new particle ID
            std::cout << "DIVISION" << std::endl;
			p[i].str();
			p[p.size()-1].str();
        }
    }
}

///Master function moveAll
void moveAll(vector<Particle> &p){
    for(Particle &part : p){
        part.clear(); //Put on top to be able to monitor the forces after a time step
        part.forceInternal();
    }
    repulsiveForce(p);
    for(Particle &part : p){
        part.removeSelfOverlap();
        part.torsionForce();
        part.move();
    }
}

void writeAll(vector<Particle> &p, ofstream &outStream, int ts){
    for(Particle &part : p){
        outStream << ts << " ";
        outStream << part.ID << " ";
        outStream << part.D << " ";
        for(Coordinate &coord : part.positions){
            outStream << coord.x << "," << coord.y << " ";
        }
        for(TwoVec &tv : part.forces){
            outStream << tv.x << "," << tv.y << " ";
        }
		outStream << ";";
    }
    outStream << std::endl;
}

int main(){
    ofstream outStream;
    outStream.open("/home/romano/Documents/Workspace/Results.txt");
    std::vector<Particle> p;
    Particle Test = Particle(0, 0, 0, 1, 1, growthRate);
    Test.positions[1] = Coordinate(1,1);
    Test.positions[2] = Coordinate(2,1);
    p.push_back(Test);
    int ts = 0;
    while(p.size() < 16){
		if(ts % 50 == 0) growAll(p);
        if(ts % 1000 == 0) writeAll(p, outStream, ts);
        moveAll(p);
        ts++;
    }
    return 0;
}
    //        Output stream to .txt, which is then read in by Matplotlib

//    int ts = 0;
//    while(p.size() < 4){
//        if(ts % 20 == 0){
//            growAll(p);
//        }
//        moveAll(p);
//        for(Particle pk : p){
//            outStream << ts << std::endl;
//            outStream << pk.ID << std::endl;
//            outStream << pk.D << std::endl;
//            for(Coordinate c : pk.positions){
//                outStream << c.x << "," << c.y << std::endl;
//            }
//            for(TwoVec f : pk.forces){
//                outStream << f.x << "," << f.y << std::endl;
//            }
//        }
//        ++ts;
//    }

