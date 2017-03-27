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

void print(double in){
    std::cout << in << std::endl;
}

void divide(Particle &pOld, Particle &pNew){

    // Find middle point of particle
    int split = npivot/2 + 1;
    /*
    Got to add angular noise to the model, but not sure if on particle or spring scale
    */
    pNew.mu = pOld.mu + randomMu(); //Growth noise instead of length noise
    pNew.Lmax = pOld.Lmax;
    pNew.D = pOld.D;
    // Make new particle from the left half of the particle positions
    pNew.positions[0] = pOld.positions[0];
    // Correct for distance D/2
    double d = dist(pOld.positions[split], pOld.positions[split-1]);
    double phi = ang(pOld.positions[split], pOld.positions[split-1]);
    pNew.positions[npivot+1].x = pOld.positions[split-1].x + (d - pOld.D/2)*cos(phi);
    pNew.positions[npivot+1].y = pOld.positions[split-1].y + (d - pOld.D/2)*sin(phi);
    //make space for an interpolated particle
    for(int i = 0; i < (npivot-1)/2; i++){
        pNew.positions[2*i+2] = pOld.positions[i+1];
    }
    //interpolation
    for(int i = 1; i <= npivot; i+=2){
        pNew.positions[i].x = (pNew.positions[i-1].x + pNew.positions[i+1].x)/2;
        pNew.positions[i].y = (pNew.positions[i-1].y + pNew.positions[i+1].y)/2;
    }
    //! Set rest length equal to total particle length divided by number of springs
    double totalLength;
    for(int i = 0; i < npivot + 1; i++){
        totalLength += dist(pNew.positions[i], pNew.positions[i+1]);
    }
    pNew.len = totalLength;
    pNew.L = totalLength/npivot;
    // ---Update old particle from the right half of the particle positions---
    d = dist(pOld.positions[split+1], pOld.positions[split]);
    phi = ang(pOld.positions[split], pOld.positions[split+1]);
    pOld.positions[0].x = pOld.positions[split+1].x + (d - pOld.D/2)*cos(phi);
    pOld.positions[0].y = pOld.positions[split+1].y + (d - pOld.D/2)*sin(phi);
    //make space for an interpolated particle
    for(int i = 0; i < (npivot-1)/2; i++){
        pOld.positions[2*i+2] = pOld.positions[i+1+split];
    }
    //interpolation
    for(int i = 1; i <= npivot; i+=2){
        pOld.positions[i].x = (pOld.positions[i-1].x + pOld.positions[i+1].x)/2;
        pOld.positions[i].y = (pOld.positions[i-1].y + pOld.positions[i+1].y)/2;
    }
    //! Set rest length equal to total particle length divided by number of springs
    totalLength = 0;
    for(int i = 0; i < npivot + 1; i++){
        totalLength += dist(pOld.positions[i], pOld.positions[i+1]);
    }
    pOld.len = totalLength;
    pOld.L = totalLength/npivot;

}

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
        p[i].grow(); //Why no straight lines without noise
        if(p[i].L > p[i].Lmax){
            Particle pnew = Particle(0, 0, 0, 0, 0, 0); //Necessary copy of all old parameters
            p.push_back(pnew); //Add new particle to p
            divide(p[i], p[p.size() - 1]); //Set new properties of daughter particles
            p[p.size()-1].ID = p.size() - 1; //Set new particle ID
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
        part.move();
    }
}

int main(){
    Particle Test = Particle(0, -0, 0, 0.25, 1, 0.01);
    std::vector<Particle> p;
    p.push_back(Test);
    int ts = 0;
    ofstream outStream;
    outStream.open("/home/romano/Documents/Workspace/Results.txt");
    int foo;
    while(p.size() < 33){
        foo = p.size();
        if(ts % 100 == 0){
            growAll(p);
            if(p.size() != foo) print( (double)ts);
        }
        moveAll(p);

        //Output stream to .txt, which is then read in by Matplotlib
        for(Particle pk : p){
            outStream << ts << std::endl;
            outStream << pk.ID << std::endl;
            outStream << pk.D << std::endl;
            for(Coordinate c : pk.positions){
                outStream << c.x << "," << c.y << std::endl;
            }
            for(TwoVec f : pk.forces){
                outStream << f.x << "," << f.y << std::endl;
            }
        }
        ++ts;
    }

    return 0;
}
