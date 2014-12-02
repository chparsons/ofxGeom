/*
 * (C) 2014 Christian Parsons
 * www.chparsons.com.ar
 *
 * taken from:
 *
 * http://deusexmachina.googlecode.com/svn/trunk/DEM/Src/nebula2/inc/mathlib/line.h
 *
 * @class line3
 * @ingroup NebulaMathDataTypes
 *
 * A 3-dimensional line class.
 *
 * (C) 2004 RadonLabs GmbH
 */

#pragma once
#include "ofxGeom.h"

class ofxLine3
{
  public:

    ofxLine3();
    ofxLine3(const ofxVec3f& v0, const ofxVec3f& v1);
    ofxLine3(const ofxLine3& l);

    /// set start and end point
    void set(const ofxVec3f& v0, const ofxVec3f& v1);
    /// get start point
    const ofxVec3f& start() const;
    /// get end point
    ofxVec3f end() const;

    /// get length
    float length() const;
    /// get squared length
    float length_squared() const;

    /// minimal distance of point to line
    float distance(const ofxVec3f& p) const;
    /// get point on line at t
    /// return p = b + m*t
    ofxVec3f interpolated(float t) const;

    /// intersect
    bool intersect(const ofxLine3& l, ofxVec3f& pa, ofxVec3f& pb) const;

    /// return t of the closest point on the line
    float closest_point(const ofxVec3f& p) const;

    /// point of origin
    ofxVec3f b;
    /// direction
    ofxVec3f m;
};

inline
ofxLine3::ofxLine3()
{}

inline
ofxLine3::ofxLine3(const ofxVec3f& v0, const ofxVec3f& v1) :
  b(v0),
  m(v1 - v0)
{}

inline
ofxLine3::ofxLine3(const ofxLine3& rhs) :
  b(rhs.b),
  m(rhs.m)
{}

inline
void
ofxLine3::set(const ofxVec3f& v0, const ofxVec3f& v1)
{
  this->b = v0;
  this->m = v1 - v0;
}

inline
const ofxVec3f&
ofxLine3::start() const
{
  return this->b;
}

inline
ofxVec3f
ofxLine3::end() const
{
  return this->b + this->m;
}

inline
float
ofxLine3::length() const
{
  return this->m.length();
}

inline
float
ofxLine3::length_squared() const
{
  return this->m.lengthSquared();
}

/**
  Returns a point on the line which is closest to a another point in space.
  This just returns the parameter t on where the point is located. If t is
  between 0 and 1, the point is on the line, otherwise not. To get the
  actual 3d point p:

  p = m + b*t
  */

inline
float
ofxLine3::closest_point(const ofxVec3f& p) const
{
  ofxVec3f diff(p - this->b);
  float l = (this->m % this->m);
  if (l > 0.0f)
  {
    float t = (this->m % diff) / l;
    return t;
  }
  else
  {
    return 0.0f;
  }
}

inline
float
ofxLine3::distance(const ofxVec3f& p) const
{
  ofxVec3f diff(p - this->b);
  float l = (this->m % this->m);
  if (l > 0.0f)
  {
    float t = (this->m % diff) / l;
    diff = diff - this->m * t;
    return diff.length();
  }
  else
  {
    // line is really a point...
    ofxVec3f v(p - this->b);
    return v.length();
  }
}

/**
  Returns p = b + m * t, given t. Note that the point is not on the line
  if 0.0 > t > 1.0
  */
inline
ofxVec3f
ofxLine3::interpolated(const float t) const
{
  return this->b + this->m * t;
  //return ofxVec3f( b + m*t );
}

/**
  Get line/line intersection. Returns the shortest line between two lines.
  */
inline
bool
ofxLine3::intersect(const ofxLine3& l, ofxVec3f& pa, ofxVec3f& pb) const
{
  const float EPS = 2.22e-16f;
  ofxVec3f p1 = this->b;
  ofxVec3f p2 = this->interpolated(10.0f);
  ofxVec3f p3 = l.b;
  ofxVec3f p4 = l.interpolated(10.0f);
  ofxVec3f p13, p43, p21;
  float d1343,d4321,d1321,d4343,d2121;
  float numer, denom;
  float mua, mub;

  p13.x = p1.x - p3.x;
  p13.y = p1.y - p3.y;
  p13.z = p1.z - p3.z;
  p43.x = p4.x - p3.x;
  p43.y = p4.y - p3.y;
  p43.z = p4.z - p3.z;
  if (n_abs(p43.x) < EPS && n_abs(p43.y) < EPS && n_abs(p43.z) < EPS) return false;
  p21.x = p2.x - p1.x;
  p21.y = p2.y - p1.y;
  p21.z = p2.z - p1.z;
  if (n_abs(p21.x)  < EPS && n_abs(p21.y) < EPS && n_abs(p21.z) < EPS) return false;

  d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
  d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
  d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
  d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
  d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;

  denom = d2121 * d4343 - d4321 * d4321;
  if (n_abs(denom) < EPS) return false;
  numer = d1343 * d4321 - d1321 * d4343;

  mua = numer / denom;
  mub = (d1343 + d4321 * (mua)) / d4343;

  pa.x = p1.x + mua * p21.x;
  pa.y = p1.y + mua * p21.y;
  pa.z = p1.z + mua * p21.z;
  pb.x = p3.x + mub * p43.x;
  pb.y = p3.y + mub * p43.y;
  pb.z = p3.z + mub * p43.z;

  return true;
}

