/*
 * zRay ~ A Monte Carlo Path Tracer
 * Zephyr Mann ~ 2013
 * SCAD ~ Studio 2
 * zephmann@gmail.com
 *
 * Image.h
 */

#ifndef IMAGE_H
#define IMAGE_H

class Image {
public:
	Pixel *data;
	int w;
	int h;

	Image(int nw, int nh)
		: w(nw), h(nh) { data = new Pixel[w*h]; }

	~Image() {
		delete [] data;
	}

	void write(string filename) const;
};

void Image::write(string filename) const {
	cout << filename << endl;

	int size = w*h;
	unsigned char *ch = new unsigned char[size*3];

	int index;
	int k = 0;
	for(int j = h-1; j > -1; j--)
	{
		for(int i = 0; i < w; i++)
		{
			index = j*w+i;
			ch[k++] = int(data[index].r*255) & 0xFF;
			ch[k++] = int(data[index].g*255) & 0xFF;
			ch[k++] = int(data[index].b*255) & 0xFF;
		}
	}

	FILE *f = fopen(filename.c_str(), "wb");         // Write image to PPM file.
	fprintf(f, "P6\n%d %d\n%d\n", w, h, 255);
	fwrite(ch, 1, size*3, f);
}

#endif
