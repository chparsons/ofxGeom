/*
 * 2014 Christian Parsons
 * www.chparsons.com.ar
 *
 * based on:
 *
 * http://deusexmachina.googlecode.com/svn/trunk/DEM/Src/nebula2/inc/mathlib/triangle.h
 *
 * A 3D triangle.
 */

#pragma once
#include "ofxGeom.h"

// Triangle points are tri(s,t) = b + s * e0 + t*e1 where
// 0 <= s <=1, 0 <= t <= 1 and 0 <= s+t <=1

class ofxTriangle 
{

  public:

    ofVec3f b, e0, e1;

    ofxTriangle() {};
    ofxTriangle(const ofVec3f& v0, const ofVec3f& v1, const ofVec3f& v2)
      : b(v0), e0(v1 - v0), e1(v2 - v0) {};
    ofxTriangle(const ofxTriangle& t)
      : b(t.b), e0(t.e0), e1(t.e1) {};

    string toString()
    {
      return "[ofxTriangle b: "+ofToString(this->b)+" e0: "+ofToString(this->e0)+" e1: "+ofToString(this->e1)+" ]";
    };

    void set(const ofVec3f& v0, const ofVec3f& v1, const ofVec3f& v2) 
    {
      b  = v0;
      e0 = v1 - v0;
      e1 = v2 - v0;
    };

    //--- get the face normal of the triangle ---------------------------------
    ofVec3f normal() const 
    {
      ofVec3f cross = e0.crossed(e1);
      cross.normalize();
      return cross;
    };

    //--- get the centroid (center of gravity) of the triangle ----------------
    ofVec3f centroid() const 
    {
      const float oneThird = 1.0f / 3.0f;
      return b + ((e0 + e1) * oneThird);
    };

    //--- get the plane of the triangle ---------------------------------------
    ofxPlane plane() const 
    {
      return ofxPlane(b, b + e0, b + e1);
    };

    //--- get one the edge points ---------------------------------------------
    ofVec3f vertex(int i) const
    {
      switch (i)
      {
        case 0: return b;
        case 1: return b + e0;
        case 2: return b + e1;
        default: return ofVec3f(0.0f, 0.0f, 0.0f);
      }
    };

    //--- check if and where line intersects triangle -------------------------
    //  Taken from Magic Software (http://www.cs.unc.edu/~eberly)
    //  Return false if line is parallel to triangle or hits its backside.

    bool intersect(const ofxLine3& line, float& ipos) 
    {

      // Compute plane of triangle, Dot(normal,X-tri.b) = 0 where 'normal' is
      // the plane normal.  If the angle between the line direction and normal
      // is small, then the line is effectively parallel to the triangle.

      const float fTolerance = 1e-04f;
      ofVec3f norm = e0.crossed(e1);
      float fDenominator = norm.dot(line.m);
      //float fLLenSqr     = line.m.dot(line.m);
      //float fNLenSqr     = norm.dot(norm);

      // check if intersecting backface or parallel...
      if (fDenominator >= -fTolerance) return false;

      //if ((fDenominator*fDenominator) <= (fTolerance*fLLenSqr*fNLenSqr)) {
      //    // line and triangle are parallel
      //    return false;
      //}

      // The line is X(t) = line.b + t*line.m.  Compute line parameter t for
      // intersection of line and plane of triangle.  Substitute in the plane
      // equation to get Dot(normal,line.b-tri.b) + t*Dot(normal,line.m)

      ofVec3f kDiff0(line.b - b);
      float fTime = -(norm.dot(kDiff0)) / fDenominator;
      if ((fTime < -fTolerance) || (fTime > (1.0f + fTolerance))) return false;

      // Find difference of intersection point of line with plane and vertex
      // of triangle.
      ofVec3f kDiff1(kDiff0 + line.m*fTime);

      // Compute if intersection point is inside triangle.  Write
      // kDiff1 = s0*E0 + s1*E1 and solve for s0 and s1.
      float fE00 = e0.dot(e0);
      float fE01 = e0.dot(e1);
      float fE11 = e1.dot(e1);
      float fDet = (float)fabs(fE00*fE11-fE01*fE01);     // = |normal|^2 > 0
      float fR0  = e0.dot(kDiff1);
      float fR1  = e1.dot(kDiff1);

      float fS0 = fE11*fR0 - fE01*fR1;
      float fS1 = fE00*fR1 - fE01*fR0;

      if ((fS0>=-fTolerance) && (fS1>=-fTolerance) && (fS0+fS1<=fDet+fTolerance)) {
        // intersection is inside triangle
        ipos = fTime;
        return true;
      }

      // intersection is outside triangle
      return false;
    };

    //--- check if and where line intersects triangle -------------------------
    //  Taken from Magic Software (http://www.cs.unc.edu/~eberly)
    //  Return false if line is parallel to triangle

    bool intersect_both_sides(const ofxLine3& line, float& ipos) {

      // Compute plane of triangle, Dot(normal,X-tri.b) = 0 where 'normal' is
      // the plane normal.  If the angle between the line direction and normal
      // is small, then the line is effectively parallel to the triangle.
      const float fTolerance = 1e-04f;
      ofVec3f norm = e0.crossed(e1);
      float fDenominator = norm.dot(line.m);
      float fLLenSqr     = line.m.dot(line.m);
      float fNLenSqr     = norm.dot(norm);

      // check if intersecting backface or parallel...
      if (fDenominator*fDenominator <= fTolerance*fLLenSqr*fNLenSqr) return false;

      //if ((fDenominator*fDenominator) <= (fTolerance*fLLenSqr*fNLenSqr)) {
      //    // line and triangle are parallel
      //    return false;
      //}

      // The line is X(t) = line.b + t*line.m.  Compute line parameter t for
      // intersection of line and plane of triangle.  Substitute in the plane
      // equation to get Dot(normal,line.b-tri.b) + t*Dot(normal,line.m)
      ofVec3f kDiff0(line.b - b);
      float fTime = -(norm.dot(kDiff0)) / fDenominator;
      if ((fTime<-fTolerance) || (fTime>(1.0f+fTolerance))) return false;

      // Find difference of intersection point of line with plane and vertex
      // of triangle.
      ofVec3f kDiff1(kDiff0 + line.m*fTime);

      // Compute if intersection point is inside triangle.  Write
      // kDiff1 = s0*E0 + s1*E1 and solve for s0 and s1.
      float fE00 = e0.dot(e0);
      float fE01 = e0.dot(e1);
      float fE11 = e1.dot(e1);
      float fDet = (float)fabs(fE00*fE11 - fE01*fE01);     // = |normal|^2 > 0
      float fR0  = e0.dot(kDiff1);
      float fR1  = e1.dot(kDiff1);

      float fS0 = fE11*fR0 - fE01*fR1;
      float fS1 = fE00*fR1 - fE01*fR0;

      if ((fS0 >= -fTolerance) && (fS1 >= -fTolerance) && (fS0 + fS1 <= fDet + fTolerance)) {
        // intersection is inside triangle
        ipos = fTime;
        return true;
      }

      // intersection is outside triangle
      return false;
    };
};

