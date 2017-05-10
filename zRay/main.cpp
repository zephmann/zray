/*
 * zRay ~ A Monte Carlo Path Tracer
 * Zephyr Mann ~ 2013
 * SCAD ~ Studio 2
 * zephmann@gmail.com
 */

#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string.h>
#include <cstdlib>

using namespace std;

class Pixel {
public:
	Pixel() {
		r = g = b = 0.;
	}

	Pixel(double _r, double _g, double _b)
		: r(_r), g(_g), b(_b) {}

	double r,g,b;
};

inline Pixel operator+(Pixel &p1, Pixel &p2) {
	Pixel tp;
	tp.r = p1.r + p2.r;
	tp.g = p1.g + p2.g;
	tp.b = p1.b + p2.b;
	return tp;
}

inline Pixel operator*(Pixel &p, double f) {
	Pixel tp;
	tp.r = p.r * f;
	tp.g = p.g * f;
	tp.b = p.b * f;
	return tp;
}

inline Pixel operator*(Pixel &p1, Pixel &p2) {
	return Pixel(p1.r*p2.r, p1.g*p2.g, p1.b*p2.b);
}

inline double operator-(Pixel &p1, Pixel &p2) {
	double r = p1.r - p2.r;
	double g = p1.g - p2.g;
	double b = p1.b - p2.b;
	return (r*r + g*g + b*b);
}

inline double rand01() {
	return rand() / double(RAND_MAX);
}


const double PI = 3.1415926535;
const double PI_2 = PI * 0.5;

#include "Geometry.h"
#include "Shape.h"
#include "Light.h"
#include "Camera.h"
#include "Image.h"
#include "Material.h"
#include "Texture.h"
#include "Scene.h"

int main(int argc, const char* argv[])
{
	int width;
	int height;
	Camera myCam;

	bool fast;
	int samps;
	int aadivs;
	double noiseThreshold;

	if(argc != 2)
	{
		cout << "Usage: ./zray.out scenefile outputfile\n";
		return 1;
	}

	string inputfile = argv[1];
	string outputfile;

	Scene myScene;
	if(!myScene.loadScene(inputfile, outputfile, myCam, width, height, fast, samps, aadivs, noiseThreshold))
	{
		cout << "ERROR: cannot load scene!\n";
		return 0;
	}

	srand(time(0));
	time_t curr = time(0);
	cout << ctime(&curr) << endl;

	init_noise();

	Image test(width, height);

	noiseThreshold *= noiseThreshold;

	double invsamps = 1. / samps;

	int aadiv2 = aadivs*aadivs;
	double aainc = 1./(aadivs-1);
	double invaa = 1./(aadiv2);

	int size = width*height;
	double mside = 1. / double(max(width, height));
	double hwid = double(width) * 0.5f;
	double hhei = double(height) * 0.5f;
	double invwid = 1. / double(width);

	if(0) {
		Ray r;
		myCam.GetRay(0.,-0.2,r);
		Pixel p = myScene.Render(r);
		return 0;
	}

	// Use Open MP for parallelization, loop through all pixels
	#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < size; i++)
	{
		// calc the x and y of the pixel
		int x = i % width;
		int y = i * invwid;

		Ray r;
		Pixel tp;
		Pixel p;
		double rx,px,py;

		// used for rendering a test image
		if(fast)
		{
			// render a number of samples for the current pixel
			for (int s = 0; s < samps; s++)
			{
				// calc the direction and randomly offset
				rx = rand01()*2 - 1;
				px = (double(x-hwid)+rx) * mside;
				rx = rand01()*2 - 1;
				py = (double(y-hhei)+rx) * mside;

				// get a new ray
				myCam.GetRay(px, py, r);

				// get the color result from the current ray
				tp = myScene.Render(r);
				p = p + tp;
			}
			p = p * invsamps;
		}

		// else render a full image
		else
		{
			int s = 0;
			double noise = 1.1;
			Pixel *aaPix = new Pixel[aadiv2];
			int aaInd;
			Pixel subPix, oldAvg, curAvg;

			// loop while we're under the number of samples and are still above the noise threshold
			while (s < samps && noise > noiseThreshold)
			{
				subPix = Pixel(0,0,0);
				aaInd = 0;
				noise = 0;

				// loop through the choosen sub-pixels
				for(double j = -0.5; j <= 0.5; j+=aainc)
				{
					for(double k = -0.5; k <= 0.5; k+=aainc)
					{
						// calc the direction and randomly offset
						rx = (rand01()*2 - 1) * aainc;
						px = (double(x-hwid)+j+rx) * mside;
						rx = (rand01()*2 - 1) * aainc;
						py = (double(y-hhei)+k+rx) * mside;

						// get a new ray
						myCam.GetRay(px, py, r);

						// get the color result from the current ray
						Pixel tp = myScene.Render(r);

						subPix = subPix + tp;

						aaPix[aaInd] = tp;
						aaInd++;
					}
				}

				// add to the total for the sub-pixel and scale by the inverse of the number of anti-alias samples
				subPix = subPix * invaa;
				p = p + subPix;

				// inc number of samples
				s++;

				// calc the average noise compared to prev pixel vs prev sample
				if(s == 1)
				{
					for(int j = 0; j < aadiv2; j++)
					{
						noise = max(noise, subPix-aaPix[j]);
					}
					oldAvg = p;
				}
				else if(s < samps)
				{
					curAvg = p * (1./s);
					noise = curAvg - oldAvg;

					oldAvg = curAvg;
				}
			}

			delete aaPix;

			if(s < samps)
				p = p * (1./s);
			else
				p = p * invsamps;

			// color by num samples
			if(0)
				p.r = p.g = p.b = (s * invsamps);
		}

		// clamp color values to 1
		p.r = min(1., p.r);
		p.g = min(1., p.g);
		p.b = min(1., p.b);

		test.data[i] = p;

		if(x == 0)
			cout << y << endl;
	}

	curr = time(0);
	cout << ctime(&curr) << endl;

	test.write(outputfile);

	if (getenv("windir"))
		system("PAUSE");

	return 0;
}
