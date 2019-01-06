#pragma once
class Field
{
public:
	//constructor and destructor
	Field();
	void set(int xSize, int ySize);
	~Field();
	friend class Fluid;
private:
	//array containing the field parameters
	float *values;

	//field dimensions
	int x;
	int y;
	int xy;
};

