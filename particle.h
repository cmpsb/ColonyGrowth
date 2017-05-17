#ifndef PARTICLE_H_INCLUDED
#define PARTICLE_H_INCLUDED

struct Particle{
    /// Constructor of the Particle structure
    Particle(double xstart, double ystart, double angle, double restlength, double diameter, double growth);

    std::vector<int> neighbours; //List of ID's of neighboring particles
    int ID; //ID number
    double D, mu, Lmax; //diameter, growth rate, maximum length of a spring
    double theta; //angular noise in initial orientation
    double L; //rest length of springs
    double len; //head to head length of the particle

    std::array<Coordinate, npivot + 2> positions;
    std::array<TwoVec, npivot + 2> forces;
    std::array<double, npivot + 2> pressures;

    int colour = 240;

	void grow();
    void str();
    void forceInternal();
    void removeSelfOverlap();
    void torsionForce();
    void move();
    void rotate();
    void clear();
	void straighten();
};

#endif // PARTICLE_H_INCLUDED
