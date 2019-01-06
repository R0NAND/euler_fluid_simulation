#include "Field.h"


Field::Field()
{
}

void Field::set(int xSize, int ySize) {
	x = xSize;
	y = ySize;
	xy = (x + 2) * (y + 2);
	values = new float[(xSize+2)*(ySize+2)];
	for (int i = 0; i < xy; i++) {
		values[i] = 0;
	}
}


Field::~Field()
{
}
