/*
 * zRay ~ A Monte Carlo Path Tracer
 * Zephyr Mann ~ 2013
 * SCAD ~ Studio 2
 * zephmann@gmail.com
 *
 * Shape.h
 */

#ifndef SHAPE_H
#define SHAPE_H

#include <stdlib.h>

#ifndef EPSILON
#define EPSILON 0.000001
#endif

#ifndef length
#define length(a) ( sizeof (a) / sizeof (a[0]) )
#endif

// 7 plane bounding volume
class BBox
{
public:
	BBox() {
		double s = sqrtf(3) / 3.f;
		planes = new Vector[7];
		planes[0] = Vector(1,0,0);
		planes[1] = Vector(0,1,0);
		planes[2] = Vector(0,0,1);
		planes[3] = Vector(s,s,s);
		planes[4] = Vector(-s,s,s);
		planes[5] = Vector(-s,-s,s);
		planes[6] = Vector(s,-s,s);

		tVals = new double[14];
		for(int i = 0; i < 14; i+=2)
		{
			tVals[i] = INFINITY;
			tVals[i+1] = -INFINITY;
		}
	}

	~BBox() {
		delete tVals;
	}

	bool Intersect(const Ray &ray, double &tmin, double &tmax, int &planeIndex) const;

	double *tVals;
	Vector *planes;
};

bool BBox::Intersect(const Ray &ray, double &tmin, double &tmax, int &planeIndex) const {
	double tn, tf, temp, numer, denom;
	tmin = -INFINITY;
	tmax = INFINITY;

	/*for(int i=0; i < 7; i++)
	{
		int tValInd = i*2;
		numer = ray.numers[i];
		denom = ray.denoms[i];
		if(denom >= 0)
		{
			tn = (tVals[tValInd] - numer) * denom;
			tf = (tVals[tValInd+1] - numer) * denom;
		}
		else
		{
			tf = (tVals[tValInd] - numer) * denom;
			tn = (tVals[tValInd+1] - numer) * denom;
		}
		if(tn > tmin)
		{
			tmin = tn;
			planeIndex = i;
		}
		if(tf < tmax)
			tmax = tf;
		if(tmin > tmax)
			return false;
	}*/

	if(ray.denoms0 >= 0)
	{
		tn = (tVals[0] - ray.numers0) * ray.denoms0;
		tf = (tVals[1] - ray.numers0) * ray.denoms0;
	}
	else
	{
		tf = (tVals[0] - ray.numers0) * ray.denoms0;
		tn = (tVals[1] - ray.numers0) * ray.denoms0;
	}
	if(tn > tmin)
	{
		tmin = tn;
		planeIndex = 0;
	}
	if(tf < tmax)
		tmax = tf;
	if(tmin > tmax)
		return false;


	if(ray.denoms1 >= 0)
	{
		tn = (tVals[2] - ray.numers1) * ray.denoms1;
		tf = (tVals[3] - ray.numers1) * ray.denoms1;
	}
	else
	{
		tf = (tVals[2] - ray.numers1) * ray.denoms1;
		tn = (tVals[3] - ray.numers1) * ray.denoms1;
	}
	if(tn > tmin)
	{
		tmin = tn;
		planeIndex = 1;
	}
	if(tf < tmax)
		tmax = tf;
	if(tmin > tmax)
		return false;


	if(ray.denoms2 >= 0)
	{
		tn = (tVals[4] - ray.numers2) * ray.denoms2;
		tf = (tVals[5] - ray.numers2) * ray.denoms2;
	}
	else
	{
		tf = (tVals[4] - ray.numers2) * ray.denoms2;
		tn = (tVals[5] - ray.numers2) * ray.denoms2;
	}
	if(tn > tmin)
	{
		tmin = tn;
		planeIndex = 2;
	}
	if(tf < tmax)
		tmax = tf;
	if(tmin > tmax)
		return false;


	if(ray.denoms3 >= 0)
	{
		tn = (tVals[6] - ray.numers3) * ray.denoms3;
		tf = (tVals[7] - ray.numers3) * ray.denoms3;
	}
	else
	{
		tf = (tVals[6] - ray.numers3) * ray.denoms3;
		tn = (tVals[7] - ray.numers3) * ray.denoms3;
	}
	if(tn > tmin)
	{
		tmin = tn;
		planeIndex = 3;
	}
	if(tf < tmax)
		tmax = tf;
	if(tmin > tmax)
		return false;


	if(ray.denoms4 >= 0)
	{
		tn = (tVals[8] - ray.numers4) * ray.denoms4;
		tf = (tVals[9] - ray.numers4) * ray.denoms4;
	}
	else
	{
		tf = (tVals[8] - ray.numers4) * ray.denoms4;
		tn = (tVals[9] - ray.numers4) * ray.denoms4;
	}
	if(tn > tmin)
	{
		tmin = tn;
		planeIndex = 4;
	}
	if(tf < tmax)
		tmax = tf;
	if(tmin > tmax)
		return false;


	if(ray.denoms5 >= 0)
	{
		tn = (tVals[10] - ray.numers5) * ray.denoms5;
		tf = (tVals[11] - ray.numers5) * ray.denoms5;
	}
	else
	{
		tf = (tVals[10] - ray.numers5) * ray.denoms5;
		tn = (tVals[11] - ray.numers5) * ray.denoms5;
	}
	if(tn > tmin)
	{
		tmin = tn;
		planeIndex = 5;
	}
	if(tf < tmax)
		tmax = tf;
	if(tmin > tmax)
		return false;


	if(ray.denoms6 >= 0)
	{
		tn = (tVals[12] - ray.numers6) * ray.denoms6;
		tf = (tVals[13] - ray.numers6) * ray.denoms6;
	}
	else
	{
		tf = (tVals[12] - ray.numers6) * ray.denoms6;
		tn = (tVals[13] - ray.numers6) * ray.denoms6;
	}
	if(tn > tmin)
	{
		tmin = tn;
		planeIndex = 6;
	}
	if(tf < tmax)
		tmax = tf;
	if(tmin > tmax)
		return false;



	return true;
}

/* AABB class
class BBox
{
public:
	BBox() {
		x0 = y0 = z0 = 0.f;
		x1 = y1 = z1 = -1.f;
	}

	BBox(double nx0, double nx1, double ny0, double ny1, double nz0, double nz1)
		: x0(nx0), x1(nx1), y0(ny0), y1(ny1), z0(nz0), z1(nz1) {}

	~BBox() {}

	bool Intersect(const Ray &ray, double &tmin, double &tmax) const;

	double x0,x1,y0,y1,z0,z1;
};

bool BBox::Intersect(const Ray &ray, double &tmin, double &tmax) const {
	double temp;

	tmin = (x0 - ray.o.x) * ray.invd.x;
	tmax = (x1 - ray.o.x) * ray.invd.x;
	if(ray.invd.x < 0)
	{
		temp = tmin;
		tmin = tmax;
		tmax = temp;
	}

	double tymin = (y0 - ray.o.y) * ray.invd.y;
	double tymax = (y1 - ray.o.y) * ray.invd.y;
	if(ray.invd.y < 0)
	{
		temp = tymin;
		tymin = tymax;
		tymax = temp;
	}

	if((tmin > tymax) || (tymin > tmax))
		return false;
	if(tymin > tmin)
		tmin = tymin;
	if(tymax < tmax)
		tmax = tymax;

	double tzmin = (z0 - ray.o.z) * ray.invd.z;
	double tzmax = (z1 - ray.o.z) * ray.invd.z;
	if(ray.invd.z < 0)
	{
		temp = tzmin;
		tzmin = tzmax;
		tzmax = temp;
	}

	if((tmin > tzmax) || (tzmin > tmax))
		return false;
	if(tzmin > tmin)
		tmin = tzmin;
	if(tzmax < tmax)
		tmax = tzmax;

	return true;
}*/
//////////////// END BBOX CLASS //////////////////////////////


class Shape
{
public:
	Shape() {}

	virtual ~Shape();
	virtual Shape* clone() const;
	virtual bool Intersect(const Ray &ray, double &t, double &u, double &v, int &tri, bool pi, int pt) const;
	virtual Normal GetNormal(const Point &p, int i, double u, double v);
	virtual void GetUV(int i, double u, double v, double &texU, double &texV);
	virtual Point GetRandomPoint();
	virtual void Test();

	int matId;
	BBox box;
	bool primary;
	bool shadows;
};

Shape::~Shape() {}

Shape* Shape::clone() const {
	return new Shape(*this);
}

bool Shape::Intersect(const Ray &ray, double &t, double &u, double &v, int &tri, bool pi, int pt) const {
	cout << "SHAPE INTERSECT\n";
	return false;
}

Normal Shape::GetNormal(const Point &p, int i, double u, double v) {
	cout << "SHAPE GET NORMAL\n";
	return Normal(0,0,0);
}

void Shape::GetUV(int i, double u, double v, double &texU, double &texV) {
	cout << "SHAPE GET UV\n";
}

Point Shape::GetRandomPoint() {
	cout << "SHAPE RANDOM POINT\n";
	return Point(0,0,0);
}

void Shape::Test() {
	cout << "Shape\n";
}
//////////////// END SHAPE CLASS /////////////////////////////


class Sphere : public Shape
{
public:
	Sphere()
		: r(0) {}
	Sphere(Point &p, double nr)
		: c(p), r(nr) {

		/*box.x0 = p.x - nr;
		box.x1 = p.x + nr;
		box.y0 = p.y - nr;
		box.y1 = p.y + nr;
		box.z0 = p.z - nr;
		box.z1 = p.z + nr;*/
	}

	~Sphere() {}

	Sphere* clone() const {
		return new Sphere( *this );
	}

	bool Intersect(const Ray &ray, double &t, double &u, double &v, int &tri, bool pi, int pt) const;
	Normal GetNormal(const Point &p, int i, double u, double v);
	void GetUV(int i, double u, double v, double &texU, double &texV);
	Point GetRandomPoint();

	void Test() {
		cout << "Sphere\n";
	}

private:
	Point c;   // center of circle
	double r;   // radius
};

bool Sphere::Intersect(const Ray &ray, double &t, double &u, double &v, int &tri, bool pi, int pt) const {
	Point to = Point(ray.o - c);

	double A = Dot(ray.d, ray.d);
	double B = Dot(to, ray.d) * 2;
	double C = Dot(to, to) - (r*r);

	double disc = (B*B) - (4*A*C);
	if (disc < 0)
		return false;

	double discSqrt = sqrtf(disc);
	double Q;

	if (B < 0)
		Q = (-B + discSqrt) * 0.5f;
	else
		Q = (-B - discSqrt) * 0.5f;

	double t0 = Q / A;
	double t1 = C / Q;

	if (t0 > t1) {
		double temp = t0;
		t0 = t1;
		t1 = temp;
	}

	if (t1 < EPSILON)
		return false;

	if (t0 < EPSILON)
		t = t1;
	else
		t = t0;

	u = v = 0.f;
	tri = -1;

	return true;
}

// returns Normalized surface normal
Normal Sphere::GetNormal(const Point &p, int i, double u, double v) {
	return Normal(Normalize(p-c));
}

void Sphere::GetUV(int i, double u, double v, double &texU, double &texV) {
	texU = u;
	texV = v;
}

Point Sphere::GetRandomPoint() {
	return Point(0,0,0);
}
//////////////// END SPHERE CLASS ////////////////////////////


class Mesh : public Shape
{
public:
	Mesh() {}
	Mesh (int nt, int *newInds, Point *newVerts, Normal *newNorms, double *newUVs, int ni, int nv, bool s);
	Mesh (const Mesh &cp);
	bool Intersect(const Ray &ray, double &t, double &u, double &v, int &tri, bool pi, int pt) const;
	bool TriIntersect(const Ray &ray, int i, double &t, double &u, double &v) const;
	Normal GetNormal(const Point &p, int i, double u, double v);
	void GetUV(int i, double u, double v, double &texU, double &texV);
	Point GetRandomPoint();

	Mesh* clone() const {
		return (new Mesh(*this));
	}

	~Mesh() {
		delete[] index;
		delete[] verts;
		delete[] norms;
	}

	void Test() {
		cout << "Mesh " << numTris  << endl;
	}

private:
	int *index;
	Point *verts;
	Normal *norms;
	double *uvs;
	int numTris;
	int numInds;
	int numVerts;
	bool sided;
};

Mesh::Mesh (int nt, int *newInds, Point *newVerts, Normal *newNorms, double *newUVs, int ni, int nv, bool s) {
	numTris = nt;
	numInds = ni;
	numVerts = nv;
	sided = s;

	index = new int[numInds];
	for(int i = 0; i < numInds; i++)
	{
		index[i] = newInds[i];
	}

	verts = new Point[numVerts];
	Point nVert;
	double d;
	for(int i = 0; i < numVerts; i++)
	{
		nVert = newVerts[i];
		verts[i] = nVert;

		for(int j = 0; j < 7; j++)
		{
			double d = Dot(nVert, box.planes[j]);
			int tValInd = j*2;

			if(d < box.tVals[tValInd])
				box.tVals[tValInd] = d;
			if(d > box.tVals[tValInd+1])
				box.tVals[tValInd+1] = d;
		}
	}

	norms = new Normal[numInds];
	for(int i = 0; i < numInds; i++)
	{
		norms[i] = newNorms[i];
	}

	int numUvs = numInds*2;
	uvs = new double[numUvs];
	for(int i = 0; i < numUvs; i++)
	{
		uvs[i] = newUVs[i];
	}
}

Mesh::Mesh (const Mesh &cp) {

	numTris = cp.numTris;
	numInds = cp.numInds;
	numVerts  = cp.numVerts;
	sided = cp.sided;

	index = new int[numInds];
	for(int i = 0; i < numInds; i++)
	{
		index[i] = cp.index[i];
	}

	verts = new Point[numVerts];
	for(int i = 0; i < numVerts; i++)
	{
		verts[i] = cp.verts[i];
	}

	norms = new Normal[numInds];
	for(int i = 0; i < numInds; i++)
	{
		norms[i] = cp.norms[i];
	}

	int numUVs = numInds*2;
	uvs = new double[numUVs];
	for(int i = 0; i < numUVs; i++)
	{
		uvs[i] = cp.uvs[i];
	}

	for(int i = 0; i < 14; i++)
	{
		box.tVals[i] = cp.box.tVals[i];
	}

	matId = cp.matId;
	primary = cp.primary;
	shadows = cp.shadows;
}

bool Mesh::Intersect(const Ray &ray, double &t, double &u, double &v, int &tri, bool pi, int pt) const {
	bool hit = false;
	double tt,tu,tv;
	tt = t = INFINITY;

	for(int i = 0; i < numTris; i++)
	{
		if( !(pi && pt==i) && TriIntersect(ray,i,tt,tu,tv) && tt > 0.f && tt < t)
		{
			tri = i;
			t = tt;
			u = tu;
			v = tv;

			hit = true;
		}
	}

	return hit;
}

bool Mesh::TriIntersect(const Ray &ray, int i, double &t, double &u, double &v) const {

	Point vert0 = verts[index[i*3]];
	Point vert1 = verts[index[i*3 + 1]];
	Point vert2 = verts[index[i*3 + 2]];

	Vector ed0 = vert1 - vert0;
	Vector ed1 = vert2 - vert0;

	Vector pvec = Cross(ray.d, ed1);
	double det = Dot(ed0, pvec);

	//if((det < EPSILON && !sided) || (det < EPSILON && det > -EPSILON))
	if(det < EPSILON && (!sided || det > -EPSILON))
		return false;

	double invdet = 1.f / det;

	Vector tvec = ray.o - vert0;
	u = Dot(tvec, pvec) * invdet;
	if(u < 0.f || u > 1.f)
		return false;

	Vector qvec = Cross(tvec, ed0);
	v = Dot(ray.d, qvec) * invdet;
	if(v < 0.f || u+v > 1.f)
		return false;

	t = Dot(ed1, qvec) * invdet;

	/*double w = 1-(u+v);

	double u1,u2,u3,v1,v2,v3;
	int uvi = i*6;
	u1 = uvs[uvi];
	v1 = uvs[uvi+1];
	u2 = uvs[uvi+2];
	v2 = uvs[uvi+3];
	u3 = uvs[uvi+4];
	v3 = uvs[uvi+5];

	double tu = w*u1 + u*u2 + v*u3;
	double tv = w*v1 + u*v2 + v*v3;
	u = tu;
	v = tv;*/

	return true;
}

Normal Mesh::GetNormal(const Point &p, int i, double u, double v) {

	Normal n0 = norms[i*3];
	Normal n1 = norms[i*3 + 1];
	Normal n2 = norms[i*3 + 2];

	double w = 1 - (u+v);

	Normal n = n0*w + n1*u + n2*v;

	return Normalize(n);

	//return Normal(1,0,0);
}

void Mesh::GetUV(int i, double u, double v, double &texU, double &texV)
{
	double w = 1-(u+v);

	double u1,u2,u3,v1,v2,v3;
	int uvi = i*6;
	u1 = uvs[uvi];
	v1 = uvs[uvi+1];
	u2 = uvs[uvi+2];
	v2 = uvs[uvi+3];
	u3 = uvs[uvi+4];
	v3 = uvs[uvi+5];

	texU = w*u1 + u*u2 + v*u3;
	texV = w*v1 + u*v2 + v*v3;
}

Point Mesh::GetRandomPoint() {
	int randTri = (rand() % numTris)*3;

	Point vert0 = verts[index[randTri]];
	Point vert1 = verts[index[randTri + 1]];
	Point vert2 = verts[index[randTri + 2]];

	double u = rand01();
	double v = rand01() * (1-u);
	double w = 1 - (u + v);

	vert0 *= u;
	vert1 *= v;
	vert2 *= w;


	Point retP = (vert0 + Vector(vert1) + Vector(vert2));

	return (vert0 + Vector(vert1) + Vector(vert2));
}

#endif
