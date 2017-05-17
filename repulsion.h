#ifndef SHORTESTDISTANCE_H_INCLUDED
#define SHORTESTDISTANCE_H_INCLUDED

#include "particle.h"

/* THE INCLUSION ABOVE CAN CAUSE UGLINESS */

Coordinate getCenter(Particle &p);

std::array<double, 5> dist3D_Segment_to_Segment(Coordinate S1a, Coordinate S1b, Coordinate S2a, Coordinate S2b);

std::array<double, 5> shortest_distance(Coordinate a0, Coordinate a1, Coordinate b0, Coordinate b1);

std::vector<std::pair<int, int>> findIndicesRepulsion(std::vector<Particle> &plist);

void applyRepulsiveForceParticleToParticle(Particle &p1, Particle &p2);

void repulsiveForce(std::vector<Particle> &plist);

#endif // SHORTESTDISTANCE_H_INCLUDED
