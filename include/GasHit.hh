// AT Simulation of He Gas Scintillator Active Target
// J.R.M Annand, University of Glasgow
// Class GasHit
// Hit storage in He
// 21/01/13 JRMA adapted from SBS equivalent

#ifndef GasHit_h
#define GasHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "TMath.h"

class G4VTouchable;

class GasHit : public G4VHit
{
public:
  GasHit();
  ~GasHit();
  GasHit(const GasHit &right);
  const GasHit& operator=(const GasHit &right);
  G4int operator==(const GasHit &right) const;
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  void Draw();
  void Print();
protected: 
  G4ThreeVector            fMom;  
  G4int                    fPID; 
  G4ParticleDefinition*    fPDef; 
  G4double                 fNPhot;
  G4int                    fkID;
  G4double                 fLe;
  //
  G4double fEdep;       // Energy deposited in detector
  G4double fEpmt[4];    // weighted energy signal in PMT
  G4double fLCpmt[4];   // primary light coll eff in PMT
  G4ThreeVector fPos;   // Position of the hit (in what frame?)
  G4int fID;            // ID of detector hit
  G4double fTime;       // global time of hit
  G4double fTimeE;      // hit time over threshold
public:
  inline void AddEnergy(G4double de) {fEdep += de;};
  inline void AddEnergyPMT(G4double de, G4int i) {fEpmt[i] += de;};
  //
  inline void SetLC(G4double lc, G4int i) { fLCpmt[i] = lc; }
  inline void SetTrackID(G4int kTkID) { fkID  = kTkID;  } 
  inline void SetTrackLength(G4double kTkLe){ fLe   = kTkLe; }
  inline void SetEdep(G4double de) { fEdep   = de; } 
  inline void SetMomentum(G4ThreeVector p) { fMom   = p; }
  inline void SetPDef(G4ParticleDefinition*  pd) { fPDef  = pd; }
  inline void SetParentID (G4int j) { fPID = j; }
  inline void SetNPhot    (G4double np) { fNPhot = np; };
  inline void SetPos(G4ThreeVector pos) { fPos = pos; };
  inline void SetID(G4int i) { fID = i; };
  inline void SetTime(G4double t) { fTime = t;};
  inline void SetTimeE(G4double t) { if( !fTimeE )fTimeE = t; };
  //
  inline G4int    GetTrackID() { return fkID; };  
  inline G4double GetTrackLength() { return fLe; }; 
  inline G4double GetEdep() { return fEdep; }
  inline G4double GetEpmt(G4int i) { return fEpmt[i]; }
  inline G4double GetLCpmt(G4int i) { return fLCpmt[i]; }
  inline G4ThreeVector GetPos() { return fPos; }
  inline G4int GetID() { return fID; }
  inline G4double GetTime() { return fTime; }
  inline G4double GetTimeE() { return fTimeE; }
  inline G4int GetParentID() { return fPID; }
  inline G4ParticleDefinition* GetPDef() { return fPDef; }
  inline G4ThreeVector GetMom() { return fMom;  }
  inline G4double GetNPhot() { return fNPhot;}  
};

typedef G4THitsCollection<GasHit> GasHitsCollection;

extern G4Allocator<GasHit> GasHitAllocator;

//----------------------------------------------------------------------------
inline void* GasHit::operator new(size_t){
  void *aHit;
  aHit = (void *) GasHitAllocator.MallocSingle();
  return aHit;
}

//----------------------------------------------------------------------------
inline void GasHit::operator delete(void *aHit){
  GasHitAllocator.FreeSingle((GasHit*) aHit);
}


#endif


