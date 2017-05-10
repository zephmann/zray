/*
 * zRay ~ A Monte Carlo Path Tracer
 * Zephyr Mann ~ 2013
 * SCAD ~ Studio 2
 * zephmann@gmail.com
 *
 * Texture.h
 */

#ifndef TEXTURE_H
#define TEXTURE_H

class Texture {
public:
	string filename;
	int width;
	int height;
	Pixel *data;

	Texture() {
		width = height = -1;
	}
	Texture(string f)
		: filename(f) {}

	bool ReadFile();


	double getGray(const double u, const double v);
	Pixel getPixel(const double u, const double v);

	~Texture() {
		delete data;
	}
};

bool Texture::ReadFile() {
	ifstream texFile;
	texFile.open(filename.c_str(), ios::in);

	string word1;

	if(texFile.is_open())
	{
		texFile >> word1;
		if (word1 != "P6")
			return false;

		texFile >> word1;

		if(word1 == "#")
		{
			getline(texFile, word1);
		}

		texFile >> word1;
		width = atoi(word1.c_str());

		if(word1 == "#")
		{
			getline(texFile, word1);
		}

		texFile >> word1;
		height = atoi(word1.c_str());

		if(word1 == "#")
		{
			getline(texFile, word1);
		}

		texFile >> word1;
		double max = 1. / atoi(word1.c_str());
		getline(texFile, word1);

		unsigned char r,g,b;
		int size = width*height;
		data = new Pixel[size];

		for(int i = 0; i < size; i++)
		{
			r = char(texFile.get());
			g = char(texFile.get());
			b = char(texFile.get());

			data[i] = Pixel(double(r)*max, double(g)*max, double(b)*max);
		}

		texFile.close();

		return true;
	}

	return false;
}

double Texture::getGray(const double u, const double v) {
	double du = fabs(fmod(u,1.));
	double dv = fabs(fmod(v,1.));

	int u1 = int(du*width);
	int u2 = u1+1;
	int v1 = int(dv*height);
	int v2 = v1+1;

	Pixel tp = data[v1*width + u1];


	return (tp.r*0.3 + tp.g*0.59 + tp.b*0.11);
}

Pixel Texture::getPixel(const double u, const double v) {
	double tu = u;
	while (tu < 0)
		tu++;
	double tv = v;
	while (tv < 0)
		tv++;

	double du = 1-fabs(fmod(tu,1.));
	double dv = 1-fabs(fmod(tv,1.));

	int u1 = int(du*width);
	int u2 = u1+1;
	int v1 = int(dv*height);
	int v2 = v1+1;

	return data[v1*width + u1];
}

#endif
