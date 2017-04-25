#include <iostream>
#include <cmath>
#include <cstdlib>
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
#include "repulsion.h"
#include "particle.h" // All clean and tidy in its own file
#include "division.h"

///Master function growAll
void growAll(vector<Particle> &p, int ts){
    int lengthBefore = p.size(); //Since a particle cannot divide twice during the same growAll() call, only loop over the initial amount of particles
    for(int i = 0; i < lengthBefore; i++){ //Unable to be range based since the ranging object changes length, leading to an infinite loop
        p[i].grow();
        if(p[i].L > p[i].Lmax){
			std::cout << "before" << std::endl;
			p[i].str();
            Particle pnew = Particle(0, 0, 0, 0, diameter, 0); //Necessary copy of all old parameters
            p.push_back(pnew); //Add new particle to p
            divide(p[i], p[p.size() - 1]); //Set new properties of daughter particles
            p[p.size()-1].ID = p.size() - 1; //Set new particle ID
            std::cout << "Division number " << p.size() - 1 << " has occurred at time step " << ts << std::endl;
			p[i].str();
			p[1].str();
        }
    }
}

///Master function moveAll
void moveAll(vector<Particle> &p){
    for(Particle &part : p){
        part.clear(); //Put on top to be able to monitor the forces after a time step
    }
    repulsiveForce(p);
    for(Particle &part : p){
		part.forceInternal();
        part.removeSelfOverlap();
        part.move();
        part.clear();
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
    std::cout << ts << std::endl;
}

void run(std::string path){
	std::cout << path << std::endl;
    ofstream outStream;
    if(onCluster) outStream.open(path);
    else outStream.open("/home/romano/Documents/Workspace/1E-2.txt");
    Particle Test(0, 0, 0, startLength, diameter, growthRate);
    std::vector<Particle> p;
    p.push_back(Test);
    int ts = 0;
    while(p.size() < Nmax){
        if(ts % relaxTime == 0) growAll(p, ts);
        if(ts % writeTime == 0) writeAll(p, outStream, ts);
        moveAll(p);
        ts++;
    }
}

int main(int argc, char* argv[]){
run("/home/romano/Documents/STUFFTODAY/today.txt");
return 0;
}


