/*
 * zRay ~ A Monte Carlo Path Tracer
 * Zephyr Mann ~ 2013
 * SCAD ~ Studio 2
 * zephmann@gmail.com
 *
 * Camera.h
 */

#ifndef CAMERA_H
#define CAMERA_H

#ifndef PI
#define PI 3.14159265359f
#endif

class Camera {
public:
	Camera() {}
	Camera(double nc, Point no, double nrx, double nry, double nrz, double nr, double ng, double nb)
		: c(nc), o(no), rx(nrx), ry(nry), rz(nrz), r(nr), g(ng), b(nb) {}

	void SetCam(double nc, Point no, double nrx, double nry, double nrz, double nr, double ng, double nb);

	void GetRay(double x, double y, Ray &r) const;

private:
	double c;          // cone angle in radians
	Point o;          // origin point
	double rx, ry, rz; // camera rotation in radians
	double r, g, b;    // background color
};

void Camera::SetCam(double nc, Point no, double nrx, double nry, double nrz, double nr, double ng, double nb) {
	c = nc;
	o = no;
	rx = nrx;
	ry = nry;
	rz = nrz;
	r = nr;
	g = ng;
	b = nb;
}

// x and y range from -1 to 1
void Camera::GetRay(double x, double y, Ray &r) const {
	double rad = tan(c) * 1.058201058f;
	double side = sqrtf(rad * rad * 0.5f);

	double nx = x*side;
	double ny = y*side;
	double nz = -1;

	// rotate around x-axis
	double tx = nx;
	double ty = ny*cos(rx) + sin(rx);
	double tz = ny*sin(rx) - cos(rx);

	// rotate around y-axis
	nx = tz*sin(ry) + tx*cos(ry);
	ny = ty;
	nz = tz*cos(ry) - tx*sin(ry);

	// rotate around z-axis
	tx = nx*cos(rz) - ny*sin(rz);
	ty = nx*sin(rz) + ny*cos(rz);
	tz = nz;

	Vector nd(tx, ty, tz);

	r = Ray(o, Normalize(nd), 0, INFINITY);
}

#endif
