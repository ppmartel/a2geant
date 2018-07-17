// AT Simulation of He Gas Scintillator Active Target
// J.R.M Annand, University of Glasgow
// Class LightColl
// Linear interpolate light collection efficiency at point x0,x1,x2
// from 3D voxel map
// 24/01/13 JRMA adapted from AcquMC original

// Interpolation method:
// Multidimensional Interpolation Methods
// InterpMethods.pdf version 04/09/2006, K.C. Johnson

#ifndef __LightColl_h__
#define __LightColl_h__

#include "G4ThreeVector.hh"
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

// For switching to higher-order interpolation schemes
enum { EDensLinear, EDensLinearC, EDensCubic, EDensCubicC };
// Switches for weighting of density fn.
// default no weight, bremsstrahlung (1/E), phase space sin(theta)
enum { EWgtUnity, EWgtBrem, EWgtSinTh, EWgtBremRes };

//-----------------------------------------------------------------------------
class LightColl{
 protected:
  G4double* fXarg;     // foam variables [0,1] interval
  G4double* fYN;       // evaluated points on grid
  G4double* fYpmt[4];
  G4double** fXN;      // grid point values for each dimension
  G4double* fUj;       // Linear "Cardinal" functions
  G4double* fXd;       // x[j] - x[fM0[j]]
  G4int* fSN;                           // # grid values for each dimension
  G4int* fM0;                           // nearest-&-less-than grid indices
  G4int* fM1;                           // fM1[j] = fM0[j] + s, s=0,1
  G4int* fM2;                           // for linear interp compensation
  G4int* fIwgt;                         // for weighting interp result
  G4int fNDim;                          // # dimensions
  G4int fNInterp;                       // # foam results
  G4int fDensOpt;                       // switch for density fn evaluation
  G4double fWgt0;                       // optional parms weighting fn.
  G4double fWgt1;
  G4double fWgt2;
  G4bool fIsInit;

  // N-dim interpolation...called by Density()
  virtual void LinIntpAccum( G4int, G4double&, Bool_t = kFALSE ); 
  virtual G4int IJK( G4int* );
  //  G4double GenFn(G4double, G4double, G4double, G4double, G4double);
  G4double YComp();
 public:
  LightColl(G4int, G4int = EDensLinear);
  virtual ~LightColl();
  virtual void ReadData( G4String, G4String );
  virtual G4double LinInterpN( G4double*, G4int, Bool_t = kFALSE );
  //  virtual G4double Density( G4int, G4double* );
  //  virtual G4double WgtDensity( G4double );
  //  virtual void Scale( G4double*, G4double* );
  virtual void SetWgt( G4int, G4int, G4double* = NULL );
  virtual G4double* GetXarg(){ return fXarg; }
  //  virtual G4double* GetScale(){ return fScale; }
  //  virtual G4double* GetXscaled(){ return fXscaled; }
  virtual G4double* GetYN(){ return fYN; }
  virtual G4double** GetXN(){ return fXN; }
  virtual G4int GetNDim(){ return fNDim; }
  virtual G4int GetNInterp(){ return fNInterp; }
  virtual G4int* GetIwgt(){ return fIwgt; }
  virtual G4int GetDensOpt(){ return fDensOpt; }
  virtual void GetGridLimits(G4double* );
};

//-----------------------------------------------------------------------------
inline void LightColl::LinIntpAccum( G4int j, G4double& yFit, Bool_t comp )
{
  // Evaluate cardinal functions recursively for each dimension j
  // and for sj = 0,1
  // Optionally apply compenstion for error in multi-linear interpolation
  for( G4int i=0; i<=1; i++ ){
    fUj[j] = (1-i) + (2*i-1)*fXd[j];
    fM1[j] = fM0[j] + i;
    if( j < (fNDim-1) ){
      LinIntpAccum( j+1, yFit, comp );
    }
    else{
      G4int jtot = IJK(fM1);
      G4double u = fYN[jtot];
      if( comp ) u -= YComp();
      for( G4int n=0; n<fNDim; n++ ) u *= fUj[n];
      yFit += u;
    }
  }
}

//-----------------------------------------------------------------------------
inline G4double LightColl::YComp()
{
  // 1st order compensation for error in linear interpolation
  // via "double differential" of grid points...check for off range

  G4double d2Y = 0.0;
  G4double Y[3];
  G4int j;
  G4int iy[3];
  for( j=0; j<fNDim; j++ )fM2[j] = fM1[j];
  for( j=0; j<fNDim; j++ ){
    if( fM1[j] == 0 ) iy[0] = 0;                       // off range low
    else if( fM1[j] == fSN[j] - 1) iy[0] = fSN[j] - 3; // off range high
    else iy[0] = fM1[j] - 1;                           // OK
    iy[1] = iy[0] + 1;
    iy[2] = iy[0] + 2;
    for( G4int i=0; i<3; i++ ){
      fM2[j] = iy[i];
      Y[i] = fYN[IJK(fM2)];
    }
    d2Y += Y[0] + Y[2] - 2*Y[1];
    fM2[j] = fM1[j];
  }
  return d2Y/12.0;
}

//-----------------------------------------------------------------------------
inline G4int LightColl::IJK(G4int* mj)
{
  // Convert N-dimensional i,j,k.... to an equivalent 1-D index
  G4int coeff = 1;
  G4int j1 = mj[fNDim-1];
  for( G4int j=fNDim-1; j>0; j-- ){
    coeff *= fSN[j];
    j1 += mj[j-1]*coeff;
  }
  return j1;
}

//-----------------------------------------------------------------------------
inline void LightColl::SetWgt( G4int iwgt, G4int i, G4double* opt  )
{
  // Initialise optional weighting for foam generator
  //
  if( i < fNDim ) fIwgt[i] = iwgt;
  if( opt ) {fWgt0 = opt[0]; fWgt1 = opt[1]; fWgt2 = opt[2];}
}

#endif
