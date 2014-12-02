/*
 * (C) 2014 Christian Parsons
 * www.chparsons.com.ar
 *
 * taken from:
 *
 * http://deusexmachina.googlecode.com/svn/trunk/DEM/Src/nebula2/inc/mathlib/plane.h
 * https://github.com/rgbdemo/nestk/
 *
 * @class plane
 * @ingroup NebulaMathDataTypes
 *
 * A plane in 3d space.
 *
 * (C) 2004 RadonLabsGmbH
 */

#pragma once
#include "ofxGeom.h"

class ofxPlane 
{

  public:

    ofxPlane();
    ofxPlane(float A, float B, float C, float D);
    ofxPlane(const ofxPlane& p);
    ofxPlane(const ofVec3f& v0, const ofVec3f& v1, const ofVec3f& v2);

    void set(float A, float B, float C, float D);
    void set(const ofVec3f& v0, const ofVec3f& v1, const ofVec3f& v2);

    float distance(const ofVec3f& v) const;
    ofVec3f normal() const;

    void intersect( const ofVec3f& p1, const ofVec3f& p2, ofVec3f& intersection ) const;
    bool intersect( const ofxLine3& line, ofVec3f& intersection ) const;
    bool intersect( const ofxPlane& plane, ofxLine3& line ) const;

    string toString();

    float a, b, c, d;
};

inline ofxPlane::ofxPlane() :
	a(0.0f), b(0.0f), c(0.0f), d(1.0f)
{}

inline ofxPlane::ofxPlane(float A, float B, float C, float D) :
	a(A), b(B), c(C), d(D)
{}

inline ofxPlane::ofxPlane(const ofxPlane& rhs) :
	a(rhs.a), b(rhs.b), c(rhs.c), d(rhs.d)
{}

inline string ofxPlane::toString()
{
  return "[ofxPlane a: "+ofToString(this->a)+" b: "+ofToString(this->b)+" c: "+ofToString(this->c)+" d: "+ofToString(this->d)+" ]";
}

inline void ofxPlane::set(float A, float B, float C, float D)
{
  this->a = A;
  this->b = B;
  this->c = C;
  this->d = D;
}

/**
  Constructs a Plane from 3 position vectors.
  */
inline void ofxPlane::set(const ofVec3f& v0, const ofVec3f& v1, const ofVec3f& v2)
{
  ofVec3f _normal = (v2 - v0).cross(v1 - v0);
  _normal.normalize();
  this->a = _normal.x;
  this->b = _normal.y;
  this->c = _normal.z;
  this->d = -(a * v0.x + b * v0.y + c * v0.z);
}

inline ofxPlane::ofxPlane(const ofVec3f& v0, const ofVec3f& v1, const ofVec3f& v2)
{
  this->set(v0, v1, v2);
}

/**
  Computes the distance of a point to the Plane. Return 0.0 if the 
  point is on the ofxPlane.
  */
inline float ofxPlane::distance(const ofVec3f& v) const
{
  return this->a * v.x + this->b * v.y + this->c * v.z + this->d;
}

/**
  Returns the Plane normal.
  */
inline ofVec3f ofxPlane::normal() const
{
  return ofVec3f(this->a, this->b, this->c);
}

/*
 * calc point of Plane/Line intersection
 * https://github.com/rgbdemo/nestk/blob/master/ntk%2Fgeometry%2Fplane.h
 */
inline void ofxPlane::intersect( const ofVec3f& p1, const ofVec3f& p2, ofVec3f& intersection ) const
{
  double u = a*p1.x + b*p1.y + c*p1.z + d;
  u /= a*(p1.x-p2.x) + b*(p1.y-p2.y) + c*(p1.z-p2.z);
  intersection = p1 + u * (p2-p1);
}

/*
 * calc point of Plane/Line intersection
 * one-sided Plane
 * returns false if the line is parallel to the Plane
 */
inline bool ofxPlane::intersect(const ofxLine3& l, ofVec3f& intersection) const
{
  float t;
  float f0 = this->a * l.b.x + this->b * l.b.y + this->c * l.b.z + this->d;
  float f1 = this->a * -l.m.x + this->b * -l.m.y + this->c * -l.m.z;

  if ((f1 < -0.0001f) || (f1 > 0.0001f))
  {
    t = f0 / f1;
    intersection = l.b + t*l.m;
    return true;
  }
  else
  {
    return false;
  }
}

/*
 * calc Plane/Plane intersection
 * returns false if Planes are parallel
 */

inline bool ofxPlane::intersect(const ofxPlane& p, ofxLine3& l) const
{
    ofVec3f n0 = this->normal();
    ofVec3f n1 = p.normal();
    float n00 = n0 % n0;
    float n01 = n0 % n1;
    float n11 = n1 % n1;
    float det = n00 * n11 - n01 * n01;
    const float tol = 1e-06f;
    if (fabs(det) < tol) 
    {
        return false;
    }
    else 
    {
        float inv_det = 1.0f/det;
        float c0 = (n11 * this->d - n01 * p.d) * inv_det;
        float c1 = (n00 * p.d - n01 * this->d) * inv_det;
        l.m = n0 * n1;
        l.b = n0 * c0 + n1 * c1;
        return true;
    }
}

