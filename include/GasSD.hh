// AT Simulation of He Gas Scintillator Active Target
// J.R.M Annand, University of Glasgow
// Class GasSD
// Accumulate experimental signal while stepping through the sensitive volume
// 24/01/13 JRMA adapted from SBS equivalent

#ifndef GasSD_h
#define GasSD_h 1

#include "GasHit.hh"
#include "G4VSensitiveDetector.hh"
#include "LightColl.hh"

class G4Step;
class G4HCofThisEvent;

class GasSD : public G4VSensitiveDetector
{
public:
  GasSD(G4String, G4int, G4int);
  ~GasSD();
  void    Initialize(G4HCofThisEvent*);
  G4bool  ProcessHits(G4Step*, G4TouchableHistory*);
  void    EndOfEvent(G4HCofThisEvent*);
  void    clear();
  void    DrawAll();
  void    PrintAll();
  //
  G4double GetEffectiveEnergyDeposit(const G4Step*);
  void SetEthr(G4double ethr){ fEthr = ethr; }
  G4double GetEThr(){ return fEthr; }
  void SetTmax(G4double tm){ if(tm > 0.0) fTmax = tm; }
  G4double GetTmax(){ return fTmax; }
  void SetLightColl(G4String voxFile, G4double mapExtend=0){
    fVoxFile = voxFile;
    if( mapExtend )fVoxGridExtend = mapExtend;
    fLightColl->ReadData(fVoxFile,fVoxTree);
    fLightColl->GetGridLimits(fGridLim);
  }
  void SetVoxMapExt(G4double ext){ fVoxGridExtend = ext;}
private:
  LightColl* fLightColl;              // position dependent light collection
  GasHitsCollection* fCollection;
  G4int              fScintYield;
  G4int fHCID;
  G4int fNelements;                   // # target elements, usually 6
  G4int * fHitID;                     // element id already hit?
  G4int * fHits;                      // elements which have been hit
  G4int fNhits;                       // number of hits
  G4double fEthr;                     // threshold energy
  G4double fTmax;                     // max time to track
  G4String fVoxFile;                  // light collection map file
  G4String fVoxTree;                  // light collection map tree
  G4double fVoxGridExtend;            // extent beyond map grid
  G4double* fGridLim;                 // position limits map grid
};

#endif

