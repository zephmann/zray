/*
 * zRay ~ A Monte Carlo Path Tracer
 * Zephyr Mann ~ 2013
 * SCAD ~ Studio 2
 * zephmann@gmail.com
 *
 * Light.h
 */

#ifndef LIGHT_H
#define LIGHT_H

class Light {
public:
	Light() {}
	Light(string nn, Point nc, Vector nd, double nr, double ng, double nb, double ni, int nt, bool ndif = true, bool nspec = true, bool nshad = false, double nsang = 0.f, int nshadr = 1, int nshadd = 1, int ndec = 0, double ncon = 0.f, double npen = 0.f, double ndrop = 0.f)
		: name(nn), c(nc), d(Normalize(nd)), r(nr), g(ng), b(nb), i(ni), type(nt), diffuse(ndif), specular(nspec), shadows(nshad), shadAngle(nsang), shadRays(nshadr), shadDepth(nshadd), decay(ndec), cone(ncon), full(npen), dropoff(ndrop) {
		angDif = full - cone;
	}

	Point c;
	Vector d;
	double r, g, b;
	double i;
	bool diffuse;
	bool specular;
	bool shadows;
	double shadAngle;
	int shadRays;
	int shadDepth;
	string name;
	int type; // 0 - dir light, 1 - point, 2 spot
	int decay;
	double cone; //cosine of cone angle
	double full; //cosing of cone + penum angle
	double angDif;
	double dropoff;
};

#endif
