#include "Fluid.h"
#include <iostream>

Fluid::Fluid()
{
}


Fluid::Fluid(int f_xSize, int f_ySize, float f_diff, float f_visc, float f_dt, int f_solverIterations) {
	xSize = f_xSize;
	ySize = f_ySize;
	diff = f_diff;
	visc = f_visc;
	dt = f_dt;
	solverIterations = f_solverIterations;
	u.set(xSize, ySize);
	v.set(xSize, ySize);
	density.set(xSize, ySize);
	u0.set(xSize, ySize);
	v0.set(xSize, ySize);
	density0.set(xSize, ySize);
	uSource.set(xSize, ySize);
	vSource.set(xSize, ySize);
	densitySource.set(xSize, ySize);
	div.set(xSize, ySize);
	p.set(xSize, ySize);
}

void Fluid::update() {

	for (int x = 1; x <= xSize; x++) {
		for (int y = 1; y <= ySize; y++) {
			u0.values[x + y * (xSize + 2)] = uSource.values[x + y * (xSize + 2)];
			v0.values[x + y * (xSize + 2)] = vSource.values[x + y * (xSize + 2)];
			uSource.values[x + y * (xSize + 2)] = 0;
			vSource.values[x + y * (xSize + 2)] = 0;

		}
	}

	diffuseVelocity();
	projectVelocity();
	swapPointers(u0.values, u.values);
	swapPointers(v0.values, v.values);
	advectVelocity();
	projectVelocity();
	
	for (int x = 1; x <= xSize; x++) {
		for (int y = 1; y <= ySize; y++) {
			density0.values[x + y * (xSize + 2)] += densitySource.values[x + y * (xSize + 2)];
			densitySource.values[x + y * (xSize + 2)] = 0;
		}
	}

	diffuseDensity();
	swapPointers(density0.values, density.values);
	advectDensity();

}

void Fluid::reset() {
	for (int x = 1; x <= xSize; x++) {
		for (int y = 1; y <= ySize; y++) {
			u0.values[x + y * (xSize + 2)] = 0;
			v0.values[x + y * (xSize + 2)] = 0;
			u.values[x + y * (xSize + 2)] = 0;
			v.values[x + y * (xSize + 2)] = 0;
			density0.values[x + y * (xSize + 2)] = 0;
			densitySource.values[x + y * (xSize + 2)] = 0;
		}
	}
}

void Fluid::render(int size, sf::Image &image) {

	//draw pixel by pixel
	for (int x = 0; x <= xSize + 1; x++) {
		for (int y = 0; y <= ySize + 1; y++) {
			if (density.values[x + y * (xSize + 2)] > 1) {
				density.values[x + y * (xSize + 2)] = 1;
			}
			for (int i = 0; i <= size; i++){
				for (int j = 0; j <= size; j++){
					image.setPixel(x*size + i, y*size + j, sf::Color(255 * density.values[x + y * (xSize + 2)], 255 * density.values[x + y * (xSize + 2)], 255 * density.values[x + y * (xSize + 2)]));
				}
			}
			density0.values[x + y * (xSize + 2)] = density.values[x + y * (xSize + 2)];
			u0.values[x + y * (xSize + 2)] = u.values[x + y * (xSize + 2)];
			v0.values[x + y * (xSize + 2)] = v.values[x + y * (xSize + 2)];
		}
	}
}
void Fluid::addSource(int x0, int y0, int x1, int y1) {
	int dx = x1 - x0;
	int dy = y1 - y0;
	if (abs(dx) >= abs(dy)) {
		int yi = 1;
		if (dy < 0) {
			yi = -1;
			dy = -dy;
		}
		int D = 2 * dy - dx;
		int y = y0;
		int elementID = 0;
		if (dx > 0) {
			for (int x = x0; x <= x1; x++) {
				elementID = x + (xSize + 2) * y;
				densitySource.values[elementID] = 1;
				uSource.values[elementID] = (x1 - x0) / dt;
				vSource.values[elementID] = (y1 - y0) / dt;
				if (D > 0) {
					y += yi;
					D -= 2 * dx;
				}
				D += 2 * dy;
			}
		}
		else if (dx < 0) {
			D = 2 * dy + dx;
			for (int x = x0; x >= x1; x--) {
				elementID = x + (xSize + 2) * y;
				densitySource.values[elementID] = 1;
				uSource.values[elementID] = (x1 - x0) / dt;
				vSource.values[elementID] = (y1 - y0) / dt;
				if (D > 0) {
					y += yi;
					D += 2 * dx;
				}
				D += 2 * dy;
			}
		}
	}
	else if (abs(dy) > abs(dx)) {
		int xi = 1;
		if (dx < 0) {
			xi = -1;
			dx = -dx;
		}
		int D = 2 * dx - dy;
		int x = x0;
		int elementID = 0;
		if (dy > 0) {
			for (int y = y0; y <= y1; y++) {
				elementID = x + (xSize + 2) * y;
				densitySource.values[elementID] = 1;
				uSource.values[elementID] = (x1 - x0) / dt;
				vSource.values[elementID] = (y1 - y0) / dt;
				if (D > 0) {
					x += xi;
					D -= 2 * dy;
				}
				D += 2 * dx;
			}
		}
		else if (dy < 0) {
			D = 2 * dx + dy;
			for (int y = y0; y >= y1; y--) {
				elementID = x + (xSize + 2) * y;
				densitySource.values[elementID] = 1;
				uSource.values[elementID] = (x1 - x0) / dt;
				vSource.values[elementID] = (y1 - y0) / dt;
				if (D > 0) {
					x += xi;
					D += 2 * dy;
				}
				D += 2 * dx;
			}
		}
	}
}

void Fluid::diffuseVelocity() {
	for (int x = 1; x <= xSize; x++) {
		for (int y = 1; y <= ySize; y++) {
			u.values[x + y * (xSize + 2)] = u0.values[x + y * (xSize + 2)] + visc*(0.25*u.values[x - 1 + y * (xSize + 2)] + 0.25*u.values[x + 1 + y * (xSize + 2)] +
				0.25*u.values[x + (y - 1)*(xSize + 2)] + 0.25*u.values[x + (y + 1)*(xSize + 2)] - u0.values[x + y * (xSize + 2)]);
			v.values[x + y * (xSize + 2)] = v0.values[x + y * (xSize + 2)] + visc*(0.25*v.values[x - 1 + y * (xSize + 2)] + 0.25*v.values[x + 1 + y * (xSize + 2)] +
				0.25*v.values[x + (y - 1)*(xSize + 2)] + 0.25*v.values[x + (y + 1)*(xSize + 2)] - v0.values[x + y * (xSize + 2)]);
		}
	}
	setBoundary(1, u);
	setBoundary(2, v);
}


void Fluid::advectVelocity() {
	for (int x = 1; x <= xSize; x++) {
		for (int y = 1; y <= ySize; y++) {
			double ax = 0;
			double ay = 0;
			int ax0 = 0;
			int ax1 = 0;
			int ay0 = 0;
			int ay1 = 0;
			double s0 = 0;
			double s1 = 0;
			double t0 = 0;
			double t1 = 0;
			float dt0 = dt * (xSize);

			ax = x - dt0 * u.values[x + xSize * y];
			ay = y - dt0 * v.values[x + xSize * y];
			if (ax < 0.5) {
				ax = 0.5;
			}
			if (ax > xSize + 0.5) {
				ax = xSize + 0.5;
			}
			ax0 = (int)ax;
			ax1 = ax0 + 1;
			if (ay < 0.5) {
				ay = 0.5;
			}
			if (ay > xSize + 0.5) {
				ay = xSize + 0.5;
			}
			ay0 = (int)ay;
			ay1 = ay0 + 1;
			s1 = ax - ax0;
			s0 = 1 - s1;
			t1 = ay - ay0;
			t0 = 1 - t1;
			u.values[x + (xSize + 2) * y] = s0 * (t0*u0.values[ax0 + ay0 * (xSize + 2)] + t1 * u0.values[ax0 + ay1 * (xSize + 2)]) +
				s1 * (t0*u0.values[ax1 + ay0 * (xSize + 2)] + t1 * u0.values[ax1 + ay1 * (xSize + 2)]);
			v.values[x + y * (xSize + 2)] = s0 * (t0*v0.values[ax0 + ay0 * (xSize + 2)] + t1 * v0.values[ax0 + ay1 * (xSize + 2)]) +
				s1 * (t0*v0.values[ax1 + ay0 * (xSize + 2)] + t1 * v0.values[ax1 + ay1 * (xSize + 2)]);
		}
	}
	setBoundary(1, u);
	setBoundary(2, v);
}


void Fluid::projectVelocity() {
	float h = 1.0 / xSize;

	//project
	for (int i = 1; i <= xSize; i++) {
		for (int j = 1; j <= ySize; j++) {
			div.values[i + j * (xSize + 2)] = -0.5*h*(u.values[i + 1 + j * (xSize + 2)] - u.values[i - 1 + j * (xSize + 2)] + v.values[i + (j + 1)*(xSize + 2)] - v.values[i + (j - 1)*(xSize + 2)]);
			p.values[i + j * (xSize + 2)] = 0;
		}
	}
	setBoundary(1, u);
	setBoundary(2, v);

	for (int k = 0; k<solverIterations; k++) {
		for (int i = 1; i <= xSize; i++) {
			for (int j = 1; j <= ySize; j++) {
				p.values[i + j * (xSize + 2)] = (div.values[i + j * (xSize + 2)] + p.values[i - 1 + j * (xSize + 2)] + p.values[i + 1 + j * (xSize + 2)] + p.values[i + (j - 1)*(xSize + 2)] + p.values[i + (j + 1)*(xSize + 2)]) / 4;
			}
		}
	}
	setBoundary(0, density);

	for (int i = 1; i <= xSize; i++) {
		for (int j = 1; j <= ySize; j++) {
			u.values[i + j * (xSize + 2)] -= 0.5*(p.values[i + 1 + j * (xSize + 2)] - p.values[i - 1 + j * (xSize + 2)]) / h;
			v.values[i + j * (xSize + 2)] -= 0.5*(p.values[i + (j + 1)*(xSize + 2)] - p.values[i + (j - 1)*(xSize + 2)]) / h;
		}
	}
	setBoundary(1, u);
	setBoundary(2, v);
}


void Fluid::diffuseDensity() {
	for (int x = 1; x <= xSize; x++) {
		for (int y = 1; y <= ySize; y++) {
			density.values[x + y * (xSize + 2)] = density0.values[x + y * (xSize + 2)] + diff*(0.25*density.values[x - 1 + y * (xSize + 2)] + 0.25*density.values[x + 1 + y * (xSize + 2)] +
				0.25*density.values[x + (y - 1)*(xSize + 2)] + 0.25*density.values[x + (y + 1)*(xSize + 2)] - density0.values[x + y * (xSize + 2)]);
		}
	}
	setBoundary(0, density);
}


void Fluid::advectDensity() {

	double ax = 0;
	double ay = 0;
	int ax0 = 0;
	int ax1 = 0;
	int ay0 = 0;
	int ay1 = 0;
	double s0 = 0;
	double s1 = 0;
	double t0 = 0;
	double t1 = 0;
	float dt0 = dt * (xSize);

	for (int x = 1; x <= xSize; x++) {
		for (int y = 1; y <= xSize; y++) {
			ax = x - dt0 * u.values[x + (xSize + 2) * y];
			ay = y - dt0 * v.values[x + (xSize + 2) * y];
			if (ax < 0.5) {
				ax = 0.5;
			}
			if (ax >(xSize + 2) - 1.5) {
				ax = (xSize + 2) - 1.5;
			}
			ax0 = (int)ax;
			ax1 = ax0 + 1;
			if (ay < 0.5) {
				ay = 0.5;
			}
			if (ay >(xSize + 2) - 1.5) {
				ay = (xSize + 2) - 1.5;
			}
			ay0 = (int)ay;
			ay1 = ay0 + 1;
			s1 = ax - ax0;
			s0 = 1 - s1;
			t1 = ay - ay0;
			t0 = 1 - t1;
			density.values[x + (xSize + 2) * y] = s0 * (t0*density0.values[ax0 + (xSize + 2) * ay0] + t1 * density0.values[ax0 + (xSize + 2) * ay1]) +
				s1 * (t0*density0.values[ax1 + (xSize + 2) * ay0] + t1 * density0.values[ax1 + (xSize + 2) * ay1]);
		}
	}
	setBoundary(0, density);
}

void Fluid::setBoundary(int d, Field f)
{
	for (int i = 1; i <= xSize; i++) {
		f.values[0 + (xSize + 2)*i] = (d == 1 ? -f.values[1+ (xSize + 2)*i] : f.values[1+ (xSize + 2)*i]);
		f.values[(xSize+1) + (xSize + 2)*i] = (d == 1 ? -f.values[xSize + (xSize + 2)*i] : f.values[xSize + (xSize + 2)* i]);
		f.values[i + (xSize + 2) * 0] = (d == 2 ? -f.values[i + (xSize + 2) * 1] : f.values[i + (xSize + 2) * 1]);
		f.values[i + (xSize + 2) * (xSize + 1)] = (d == 2 ? -f.values[i + (xSize + 2) * (xSize)] : f.values[i + (xSize + 2) * (xSize)]);
		/*v.values[0 + (xSize + 2)*i] = 0;
		v.values[(xSize+1) + (xSize + 2)*i] = 0;
		v.values[i + (xSize + 2) * 0] = 0;
		v.values[i + (xSize + 2) * (xSize + 1)] = 0;
		u.values[0 + (xSize + 2)*i] = 0;
		u.values[(xSize+1) + (xSize + 2)*i] = 0;
		u.values[i + (xSize + 2) * 0] = 0;
		u.values[i + (xSize + 2) * (xSize + 1)] = 0;*/
	}
	f.values[0 + (xSize + 2) * 0] = 0.5*(f.values[1 + (xSize + 2) * 0] + f.values[0 + (xSize + 2) * 1]);
	f.values[0 + (xSize + 2) * (xSize + 1)] = 0.5*(f.values[1 + (xSize + 2) * (xSize + 1)] + f.values[0 + (xSize + 2) * (xSize)]);
	f.values[xSize+1 + 0] = 0.5*(f.values[xSize + 0] + f.values[xSize + 1 + 1]);
	f.values[xSize+1 + (xSize + 2) * (xSize + 1)] = 0.5*(f.values[xSize + (xSize + 2) * (xSize + 1)] + f.values[xSize + 1 + (xSize + 2) * (xSize)]);

}

void Fluid::swapPointers(float *&x0, float *&x1) {
	float *memory = x0;
	x0 = x1;
	x1 = memory;
}

Fluid::~Fluid()
{
}