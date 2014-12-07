/*
 * 2014 Christian Parsons
 * www.chparsons.com.ar
 *
 * based on:
 *
 * http://deusexmachina.googlecode.com/svn/trunk/DEM/Src/nebula2/inc/mathlib/line.h
 *
 * A 2-dimensional line class.
 */

#pragma once
#include "ofxGeom.h"

class ofxLine2
{
  public:

    ofxLine2();
    ofxLine2(const ofVec2f& v0, const ofVec2f& v1);
    ofxLine2(const ofxLine2& rhs);

    /// return start point
    const ofVec2f& start() const;
    /// return end point
    ofVec2f end() const;

    /// return length
    float length() const;
    /// get point on line given t
    ofVec2f interpolated(const float t) const;

    bool intersect(const ofxLine2& l, ofVec2f& intersection) const;

    /// point of origin
    ofVec2f b;
    /// direction
    ofVec2f m;
};

inline
ofxLine2::ofxLine2()
{}

inline
ofxLine2::ofxLine2(const ofVec2f& v0, const ofVec2f& v1) :
  b(v0),
  m(v1 - v0)
{}

inline
ofxLine2::ofxLine2(const ofxLine2& rhs) :
  b(rhs.b),
  m(rhs.m)
{}

inline
const ofVec2f&
ofxLine2::start() const
{
  return this->b;
}

inline
ofVec2f
ofxLine2::end() const
{
  return this->b + this->m;
}

inline
float
ofxLine2::length() const
{
  return m.length();
}

inline
ofVec2f
ofxLine2::interpolated(const float t) const
{
  return this->b + this->m * t;
  //return ofVec2f(b + m * t);
}

/*
 * line/line intersection
 * Return FALSE if the lines don't intersect
 * from http://paulbourke.net/geometry/pointlineplane/
 */
inline
bool
ofxLine2::intersect(const ofxLine2& l, ofVec2f& intersection) const
{
  const float EPS = 2.22e-16f;

  ofVec2f _start = start();
  ofVec2f _end = end();

  ofVec2f _lstart = l.start();
  ofVec2f _lend = l.end();

  float x1 = _start.x;
  float y1 = _start.y;
  float x2 = _end.x;
  float y2 = _end.y;

  float x3 = _lstart.x;
  float y3 = _lstart.y;
  float x4 = _lend.x;
  float y4 = _lend.y;

  float mua,mub;
  float denom,numera,numerb;

  denom  = (y4-y3) * (x2-x1) - (x4-x3) * (y2-y1);
  numera = (x4-x3) * (y1-y3) - (y4-y3) * (x1-x3);
  numerb = (x2-x1) * (y1-y3) - (y2-y1) * (x1-x3);

  /* Are the line coincident? */
  if (fabs(numera) < EPS && fabs(numerb) < EPS && fabs(denom) < EPS) 
  {
    intersection.x = (x1 + x2) / 2;
    intersection.y = (y1 + y2) / 2;
    return true;
  }

  /* Are the line parallel */
  if (fabs(denom) < EPS) 
  {
    intersection.x = 0;
    intersection.y = 0;
    return false;
  }

  /* Is the intersection along the the segments */
  mua = numera / denom;
  mub = numerb / denom;

  if (mua < 0 || mua > 1 || mub < 0 || mub > 1) 
  {
    intersection.x = 0;
    intersection.y = 0;
    return false;
  }

  intersection.x = x1 + mua * (x2 - x1);
  intersection.y = y1 + mua * (y2 - y1);
  return true;
}

