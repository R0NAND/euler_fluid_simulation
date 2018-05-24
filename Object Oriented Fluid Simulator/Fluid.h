#pragma once
#include "Field.h"
#include <SFML\Graphics.hpp>
class Fluid
{
public:
	//constructor and destructor
	Fluid();
	Fluid(int f_xSize, int f_ySize, float f_diff, float f_visc, float f_dt, int f_solverIterations);
	~Fluid();

	//the public functions
	void addSource(int x0, int y0, int x, int y);
	void update();
	void reset();
	void render(sf::Image &image);


private:

	//Field property declarations
	Field u;
	Field v;
	Field density;
	Field u0;
	Field v0;
	Field density0;
	Field uSource;
	Field vSource;
	Field densitySource;
	Field div;
	Field p;

	//Physical Property declarations
	float diff;
	float visc;

	//Simulation parameters
	float dt;
	int solverIterations;
	int xSize;
	int ySize;

	//fluid physics functions
	void diffuseVelocity();
	void diffuseDensity();
	void advectVelocity();
	void advectDensity();
	void projectVelocity();
	void setBoundary(int d, Field x);
	void swapPointers(float *&x0, float *&x1);



};

