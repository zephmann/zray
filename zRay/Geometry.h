/*
 * zRay ~ A Monte Carlo Path Tracer
 * Zephyr Mann ~ 2013
 * SCAD ~ Studio 2
 * zephmann@gmail.com
 *
 * Geometry.h
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#ifndef INFINITY
#define INFINITY FLT_MAX
#endif

class Vector;
class Point;
class Ray;
class Normal;


class Vector
{
public:
	double x,y,z;

	Vector() {x = y = z = 0.;}
	Vector(double nx, double ny, double nz)
		: x(nx), y(ny), z(nz) {}
	Vector(const Pixel &p)
		: x(p.r), y(p.g), z(p.b) {}
	explicit Vector(const Normal &n);
	explicit Vector(const Point &p);

	Vector operator+(const Vector &v) const {
		return Vector(x + v.x, y + v.y, z + v.z);
	}

	Vector& operator+=(const Vector &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}

	Vector operator*(double f) const {
		return Vector(f*x, f*y, f*z);
	}

	Vector operator*=(double f) {
		x *= f; y *= f; z *= f;
		return *this;
	}

	Vector operator/(double f) const {
		if(f != 0)
		{
			double inv = 1. / f;
			return Vector(inv*x, inv*y, inv*z);
		}
		return Vector(0,0,0);
	}

	Vector& operator/=(double f) {
		if(f != 0)
		{
			double inv = 1. / f;
			x *= inv; y *= inv; z *= inv;
			return *this;
		}
		x = y = z = 0;
		return *this;
	}

	double LengthSquared() const { return x*x + y*y + z*z; }
	double Length() const { return sqrtf(x*x + y*y + z*z); }
};
////////////////////// END VECTOR ////////////////////////////////////////


class Point
{
public:
	double x,y,z;

	Point() {x = y = z = 0.;}
	Point(double nx, double ny, double nz)
		: x(nx), y(ny), z(nz) {}
	Point(const Vector &v)
		: x(v.x), y(v.y), z(v.z) {}

	Point operator+(const Vector &v) const {
		return Point(x + v.x, y + v.y, z + v.z);
	}

	Point& operator+=(const Vector &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}

	Vector operator-(const Point &p) const {
		return Vector(x - p.x, y - p.y, z - p.z);
	}

	Point operator*(double f) const {
		return Point(f*x, f*y, f*z);
	}

	Point& operator*=(double f) {
		x *= f; y *= f; z *= f;
		return *this;
	}

};
/////////////////////// END POINT /////////////////////////////////////////


inline double Dot(const Vector &v1, const Vector &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline double Dot(const Point &p1, const Vector &v2) {
	return p1.x * v2.x + p1.y * v2.y + p1.z * v2.z;
}

inline double Dot(const Vector &v1, const Point &p2) {
	return v1.x * p2.x + v1.y * p2.y + v1.z * p2.z;
}

inline double Dot(const Point &p1, const Point &p2) {
	return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}


class Ray {
public:
	Point o;
	Vector d;
	//double *numers;
	//double *denoms;

	double numers0;
	double numers1;
	double numers2;
	double numers3;
	double numers4;
	double numers5;
	double numers6;

	double denoms0;
	double denoms1;
	double denoms2;
	double denoms3;
	double denoms4;
	double denoms5;
	double denoms6;

	mutable double mint, maxt;
	int depth;
	int diffDepth;
	double opacDepth;
	int index;
	int tri;
	bool direct;

	Ray() : mint(0.f), maxt(INFINITY), depth(0), diffDepth(0), opacDepth(1.f) {}
	Ray(const Point &no, const Vector &nd, double start = 0, double end = INFINITY, int dd = 0, int dif = 0, double op = 1.f, int i = -1, int t = -1, bool _direct = false)
		: o(no), d(nd), mint(start), maxt(end), depth(dd), diffDepth(dif), opacDepth(op), index(i), tri(t), direct(_direct) {

		double s = sqrtf(3) / 3.f;
		Vector *planes = new Vector[7];
		planes[0] = Vector(1,0,0);
		planes[1] = Vector(0,1,0);
		planes[2] = Vector(0,0,1);
		planes[3] = Vector(s,s,s);
		planes[4] = Vector(-s,s,s);
		planes[5] = Vector(-s,-s,s);
		planes[6] = Vector(s,-s,s);
		
		numers0 = Dot(planes[0], no);
		numers1 = Dot(planes[1], no);
		numers2 = Dot(planes[2], no);
		numers3 = Dot(planes[3], no);
		numers4 = Dot(planes[4], no);
		numers5 = Dot(planes[5], no);
		numers6 = Dot(planes[6], no);

		denoms0 = 1 / Dot(planes[0], nd);
		denoms1 = 1 / Dot(planes[1], nd);
		denoms2 = 1 / Dot(planes[2], nd);
		denoms3 = 1 / Dot(planes[3], nd);
		denoms4 = 1 / Dot(planes[4], nd);
		denoms5 = 1 / Dot(planes[5], nd);
		denoms6 = 1 / Dot(planes[6], nd);

		delete planes;
	}

	Point operator()(double t) const { return o + d*t; }
};
////////////////////////////// END RAY /////////////////////////////////////////////


class Normal {
public:
	double x,y,z;

	Normal() { x=y=z=0.f; }
	Normal(double nx, double ny, double nz)
		: x(nx), y(ny), z(nz) {}

	Normal operator-() const {
		return Normal(-x, -y, -z);
	}

	Normal operator+(const Normal &n) const {
		return Normal(x + n.x, y + n.y, z + n.z);
	}

	Normal& operator+=(const Normal &n) {
		x += n.x; y += n.y; z += n.z;
		return *this;
	}

	Normal operator-(const Normal &n) const {
		return Normal(x - n.x, y - n.y, z - n.z);
	}

	Normal operator-=(const Normal &n) {
		x -= n.x; y -= n.y; z -= n.z;
		return *this;
	}

	Normal operator*(double f) const {
		return Normal(f*x, f*y, f*z);
	}

	Normal& operator*=(double f) {
		x *= f; y *= f; z *= f;
		return *this;
	}

	Normal operator/(double f) const {
		if(f != 0)
		{
			double inv = 1.f / f;
			return Normal(inv*x, inv*y, inv*z);
		}
		return Normal(0,0,0);
	}

	Normal& operator/=(double f) {
		if(f != 0)
		{
			double inv = 1.f / f;
			x *= inv; y *= inv; z *= inv;
			return *this;
		}
		x = y = z = 0;
		return *this;
	}

	double LengthSquared() const {
		return x*x + y*y + z*z;
	}

	double Length() const {
		return sqrtf(x*x + y*y + z*z);
	}

	double operator[](int i) const {
		return (&x)[i];
	}

	double& operator[](int i) {
		return (&x)[i];
	}

	explicit Normal(const Vector &v)
		: x(v.x), y(v.y), z(v.z) {}
};
//////////////////////////// END NORMAL ///////////////////////////////////////

inline Vector::Vector(const Point &p)
	: x(p.x), y(p.y), z(p.z) {}

inline Vector::Vector(const Normal &n)
	: x(n.x), y(n.y), z(n.z) {}

inline Vector operator*(double f, const Vector &v) {
	return v*f;
}

inline Vector operator-(const Normal &n, const Vector &v) {
	return Vector(n.x - v.x, n.y - v.y, n.z - v.z);
}

inline Vector operator-(const Vector &v, const Vector &n) {
	return Vector(v.x - n.x, v.y - n.y, v.z - n.z);
}

inline double Dot(const Normal &n1, const Vector &v2) {
	return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}

inline double Dot(const Vector &v1, const Normal &n2) {
	return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}

inline double Dot(const Normal &n1, const Normal &n2) {
	return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

inline double AbsDot(const Vector &v1, const Vector &v2) {
	return fabsf( v1.x * v2.x + v1.y * v2.y + v1.z * v2.z );
}

inline Vector Cross(const Vector &v1, const Vector &v2) {
	return Vector((v1.y * v2.z) - (v1.z * v2.y),
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x));
}

inline Vector Normalize(const Vector &v) {
	return v / v.Length();
}

inline Normal Normalize(const Normal &n) {
	return n / n.Length();
}

inline Point operator*(double f, const Point p) {
	return p*f;
}

inline Normal Faceforward(const Normal &n, const Vector &v) {
	return (Dot(n, v) < 0.f) ? -n : n;
}

inline Normal operator*(double f, const Normal n) {
	return n*f;
}

#endif
