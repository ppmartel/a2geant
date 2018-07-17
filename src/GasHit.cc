// AT Simulation of He Gas Scintillator Active Target
// J.R.M Annand, University of Glasgow
// Class GasHit
// Hit storage in He
// 21/01/13 JRMA adapted from SBS equivalent

#include "GasHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

G4Allocator<GasHit> GasHitAllocator;

//----------------------------------------------------------------------------
GasHit::GasHit()
{
  fMom   = G4ThreeVector(0,0,0);
  fPDef  = 0;
  fPID   = 0;
  fNPhot = 0;
  fEdep = 0;
  for(Int_t i=0; i<4; i++){ fEpmt[i] = 0.0; fLCpmt[i] = 0.0; }
  fPos.set(0,0,0);
  fID = 0;
  fTime = 0;
  fTimeE = 0;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
GasHit::~GasHit()
{;}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
GasHit::GasHit(const GasHit &right)
  : G4VHit()
{
  fMom    + right.fMom;
  fPDef   = right.fPDef;
  fPID    = right.fPID;
  fNPhot  = right.fNPhot;
  //
  fEdep = right.fEdep;
  for(Int_t i=0; i<4; i++){ 
    fEpmt[i] = right.fEpmt[i];
    fLCpmt[i] = right.fLCpmt[i];
  }
  fPos = right.fPos;
  fID = right.fID;
  fTime = right.fTime;
  fTimeE = right.fTimeE;
}

//-----------------------------------------------------------------------------
const GasHit& GasHit::operator=(const GasHit &right){
  fMom    + right.fMom;
  fPDef   = right.fPDef;
  fPID    = right.fPID;
  fNPhot  = right.fNPhot;
  fEdep = right.fEdep;
  for(Int_t i=0; i<4; i++){ 
    fEpmt[i] = right.fEpmt[i]; 
    fLCpmt[i] = right.fLCpmt[i]; 
  }
  fPos = right.fPos;
  fID = right.fID;
  fTime = right.fTime;
  fTimeE = right.fTimeE;
  return *this;
}

//----------------------------------------------------------------------------
void GasHit::Draw()
{;}

//-----------------------------------------------------------------------------
void GasHit::Print()
{;}









