/*
 * (C) 2014 Christian Parsons
 * www.chparsons.com.ar
 *
 * taken from:
 *
 * http://deusexmachina.googlecode.com/svn/trunk/DEM/Src/nebula2/inc/mathlib/line.h
 *
 * @class line2
 * @ingroup NebulaMathDataTypes
 *
 * A 2-dimensional line class.
 *
 * (C) 2004 RadonLabs GmbH
 */

#pragma once
#include "ofxGeom.h"

class ofxLine2
{
  public:

    ofxLine2();
    ofxLine2(const ofxVec2f& v0, const ofxVec2f& v1);
    ofxLine2(const ofxLine2& rhs);

    /// return start point
    const ofxVec2f& start() const;
    /// return end point
    ofxVec2f end() const;

    /// return length
    float length() const;
    /// get point on line given t
    ofxVec2f interpolated(const float t) const;

    /// point of origin
    ofxVec2f b;
    /// direction
    ofxVec2f m;
};

inline
ofxLine2::ofxLine2()
{}

inline
ofxLine2::ofxLine2(const ofxVec2f& v0, const ofxVec2f& v1) :
  b(v0),
  m(v1 - v0)
{}

inline
ofxLine2::ofxLine2(const ofxLine2& rhs) :
  b(rhs.b),
  m(rhs.m)
{}

inline
const ofxVec2f&
ofxLine2::start() const
{
  return this->b;
}

inline
ofxVec2f
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
ofxVec2f
ofxLine2::interpolated(const float t) const
{
  return this->b + this->m * t;
  //return ofxVec2f(b + m * t);
}
