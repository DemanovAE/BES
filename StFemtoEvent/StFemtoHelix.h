/**
 * \class StFemtoHelix
 * \brief Helix parametrization that uses TVector3
 *
 * Parametrization of a helix. Can also work with straight tracks, i.e. with
 * zero curvature. This represents only the mathematical model of a helix.
 *
 * \author Grigory Nigmatkulov; e-mail: nigmatkulov@gmail.com
 * \date May 07 2018
 */

#ifndef StFemtoHelix_h
#define StFemtoHelix_h

// C++ headers
#include <math.h>
#include <utility>
#include <algorithm>

// ROOT headers
#include "TVector3.h"

// FemtoDst headers
#ifdef _VANILLA_ROOT_
#include "SystemOfUnits.h"
#else
#include "StarClassLibrary/SystemOfUnits.h"
#endif

// Declare C++ namespaces
#if !defined(ST_NO_NAMESPACES)
using std::pair;
using std::swap;
using std::max;
#endif

//_________________
class StFemtoHelix {
  
 public:
  /// Default constructor
  StFemtoHelix();
  
  /// Constructor that takes next arguments:
  /// curvature, dip angle, phase, origin, h
  StFemtoHelix(Double_t c, Double_t dip, Double_t phase,
	      const TVector3& o, Int_t h=-1);

  /// Copy constructor
  StFemtoHelix(const StFemtoHelix&);

  // Assignment operator (will use the one, provided by compiler)
  //StFemtoHelix& operator=(const StFemtoHelix&);
  
  /// Destructor
  virtual ~StFemtoHelix();

  /// Return dip angle
  Double_t dipAngle()   const;
  /// Return curvature: 1/R in xy-plane
  Double_t curvature()  const;
  /// Return phase: aziumth in xy-plane measured from ring center
  Double_t phase()      const;
  /// Return x-center of circle in xy-plane
  Double_t xcenter()    const;
  /// Return y-center of circle in xy-plane
  Double_t ycenter()    const;
  /// Return -sign(q*B);
  Int_t    h()          const;

  /// Return origin of the helix = starting point
  const TVector3& origin() const;

  /// Set helix parameters
  void setParameters(Double_t c, Double_t dip, Double_t phase, const TVector3& o, Int_t h);

  /// coordinates of helix at point s
  Double_t x(Double_t s)  const;
  Double_t y(Double_t s)  const;
  Double_t z(Double_t s)  const;
  TVector3 at(Double_t s) const;

  /// pointing vector of helix at point s
  Double_t cx(Double_t s)  const;
  Double_t cy(Double_t s)  const;
  Double_t cz(Double_t s = 0)  const;
  TVector3 cat(Double_t s) const;

  /// returns period length of helix
  Double_t period()       const;
    
  /// path length at given r (cylindrical r)
  pair<Double_t, Double_t> pathLength(Double_t r)   const;
    
  /// path length at given r (cylindrical r, cylinder axis at x,y)
  pair<Double_t, Double_t> pathLength(Double_t r, Double_t x, Double_t y);
    
  /// path length at distance of closest approach to a given point
  Double_t pathLength(const TVector3& p, Bool_t scanPeriods = true) const;
    
  /// path length at intersection with plane
  Double_t pathLength(const TVector3& r,
		      const TVector3& n) const;

  /// path length at distance of closest approach in the xy-plane to a given point
  Double_t pathLength(Double_t x, Double_t y) const;

  /// path lengths at dca between two helices 
  pair<Double_t, Double_t> pathLengths(const StFemtoHelix&,
				       Double_t minStepSize = 10*micrometer,
				       Double_t minRange = 10*centimeter) const;
    
  /// minimal distance between point and helix
  Double_t distance(const TVector3& p, Bool_t scanPeriods = true) const;    
    
  /// checks for valid parametrization
  Bool_t valid(Double_t world = 1.e+5) const { return !bad(world); }
  Int_t  bad(Double_t world = 1.e+5) const;
    
  /// Move the origin along the helix to s which becomes then s=0
  virtual void moveOrigin(Double_t s);
  
  static const Double_t NoSolution;
    
 protected:
  
  /// Set curvature of the helix
  void setCurvature(Double_t);	/// performs also various checks
  /// Set phase of the helix
  void setPhase(Double_t);
  /// Set dip angle of the helix
  void setDipAngle(Double_t);
  /// Value of S where distance in x-y plane is minimal
  Double_t fudgePathLength(const TVector3&) const;

 protected:
  /// true for straight line case (B=0)
  Bool_t    mSingularity;  
  /// starting point of a helix
  TVector3  mOrigin;
  /// Dip angle
  Double_t  mDipAngle;
  /// Curvature = 1/R
  Double_t  mCurvature;
  /// Phase
  Double_t  mPhase;
  /// -sign(q*B);
  Int_t     mH;

  /// Cos of dip angle
  Double_t mCosDipAngle;
  /// Sin of dip angle
  Double_t mSinDipAngle;
  /// Cos of phase
  Double_t mCosPhase;
  //// Sin of phase
  Double_t mSinPhase;
    
  ClassDef(StFemtoHelix,1)
};

//
//     Non-member functions
//
Int_t operator== (const StFemtoHelix&, const StFemtoHelix&);
Int_t operator!= (const StFemtoHelix&, const StFemtoHelix&);
std::ostream& operator<<(std::ostream&, const StFemtoHelix&);

//
//     Inline functions
//
inline Int_t StFemtoHelix::h() const {return mH;}

inline Double_t StFemtoHelix::dipAngle() const {return mDipAngle;}

inline Double_t StFemtoHelix::curvature() const {return mCurvature;}

inline Double_t StFemtoHelix::phase() const {return mPhase;}

inline Double_t StFemtoHelix::x(Double_t s) const {
  if (mSingularity)
    return mOrigin.x() - s*mCosDipAngle*mSinPhase;
  else
    return mOrigin.x() + (cos(mPhase + s*mH*mCurvature*mCosDipAngle)-mCosPhase)/mCurvature;
}
 
inline Double_t StFemtoHelix::y(Double_t s) const {
  if (mSingularity)
    return mOrigin.y() + s*mCosDipAngle*mCosPhase;
  else
    return mOrigin.y() + (sin(mPhase + s*mH*mCurvature*mCosDipAngle)-mSinPhase)/mCurvature;
}

inline Double_t StFemtoHelix::z(Double_t s) const {
  return mOrigin.z() + s*mSinDipAngle;
}

inline Double_t StFemtoHelix::cx(Double_t s) const {
  if (mSingularity)
    return -mCosDipAngle*mSinPhase;
  else
    return -sin(mPhase + s*mH*mCurvature*mCosDipAngle)*mH*mCosDipAngle;
}

inline Double_t StFemtoHelix::cy(Double_t s) const {
  if (mSingularity)
    return mCosDipAngle*mCosPhase;
  else
    return cos(mPhase + s*mH*mCurvature*mCosDipAngle)*mH*mCosDipAngle;
}

inline Double_t StFemtoHelix::cz(Double_t /* s */)  const { return mSinDipAngle; }    
inline const TVector3& StFemtoHelix::origin() const { return mOrigin; }
inline TVector3 StFemtoHelix::at(Double_t s) const { return TVector3(x(s), y(s), z(s)); }
inline TVector3 StFemtoHelix::cat(Double_t s) const { return TVector3(cx(s), cy(s), cz(s)); }
inline Double_t StFemtoHelix::pathLength(Double_t X, Double_t Y) const { return fudgePathLength(TVector3(X, Y, 0)); }

inline Int_t StFemtoHelix::bad(Double_t WorldSize) const {

  Int_t ierr;
  if ( !::finite(mDipAngle) ) {
    return 11;
  }
  if ( !::finite(mCurvature) ) {
    return 12;
  }
  
  //ierr = mOrigin.bad(WorldSize);

  // The line above is commented and the StThreeVector::bad(double)
  // is rewritten here
  for(Int_t iIter=0; iIter<3; iIter++) {

    Double_t tmpVal;
    // Value StThreeVector.mX1[iter] ???
    switch(iIter) {
    case 0: tmpVal = mOrigin.X(); break;
    case 1: tmpVal = mOrigin.Y(); break;
    case 2: tmpVal = mOrigin.Z(); break;
    default: tmpVal = NAN;
    };
    
    if ( !::finite( tmpVal ) ) {
      ierr = 10 + iIter;
    }
    if ( ::fabs( tmpVal ) > WorldSize ) {
      ierr = 20 + iIter;
    }
  } //for(Int_t iIter=0; iIter<3; iIter+)
  
  if (ierr) {
    return (3+ierr*100);
  }
  if ( ::fabs(mDipAngle) > 1.58 ) {
    return 21;
  }
  
  Double_t qwe = ::fabs( ::fabs(mDipAngle) - M_PI/2 );
  if ( qwe < 1./WorldSize ) {
    return 31;
  }

  if ( ::fabs(mCurvature) > WorldSize )	{
    return 22;
  }
  if ( mCurvature < 0 )	{
    return 32;
  }
  if (abs(mH) != 1 ) {
    return 24;
  }

  return 0;
}

#endif
