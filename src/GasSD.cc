// AT Simulation of He Gas Scintillator Active Target
// J.R.M Annand, University of Glasgow
// Class GasSD
// Accumulate experimental signal while stepping through the sensitive volume
// 24/01/13 JRMA adapted from SBS equivalent

#include "GasSD.hh"
#include "GasHit.hh"
#include "A2DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"

//----------------------------------------------------------------------------
GasSD::GasSD(G4String name, G4int ScintYield, G4int nElements)
  :G4VSensitiveDetector(name)
{
  collectionName.insert("GasSDHits");
  fCollection=NULL;
  fScintYield = ScintYield;
  fNelements=nElements+1;                     //numbering starts from 1 not 0
  fHitID=new G4int[fNelements];
  for(G4int i=0;i<fNelements;i++)fHitID[i]=-1;
  fHits=new G4int[fNelements];
  for(G4int i=0;i<fNelements;i++)fHits[i]=0;
  fNhits=0;
  fHCID=-1;
  fEthr = 1.0;
  fTmax = 500.0;                             // default max time to track
  fLightColl = new LightColl(3);             // 3D interpolation
  fVoxTree = "tree";
  fVoxFile = "root/VoxMap-2.root";
  fVoxGridExtend = 0.0;
  fGridLim = new G4double[6];
}

//-----------------------------------------------------------------------------
GasSD::~GasSD()
{}

//----------------------------------------------------------------------------
void GasSD::Initialize(G4HCofThisEvent*)
{
  fCollection = new GasHitsCollection
                      (SensitiveDetectorName,collectionName[0]);
  //A way to keep all the hits of this event in one place if needed
  ///  static G4int HCID = -1;
  //  if(HCID<0){ 
  //    HCID = GetCollectionID(0); 
  //  }
  //  HCE->AddHitsCollection( HCID, fCollection );
}

//----------------------------------------------------------------------------
G4bool GasSD::ProcessHits(G4Step* aStep,G4TouchableHistory* )
{
  //  
  // Check 1st that track is still alive
  G4Track* aTrack = aStep->GetTrack();
  if( aTrack->GetTrackStatus() != fAlive )
    return false;
  // Check that energy has been deposited
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;
  //G4int trID = aTrack->GetTrackID();
  // G4int pdgCode = aTrack->GetDynamicParticle()->GetPDGcode();
  //G4int parID = aTrack->GetParentID();
  //G4int stepID = aTrack->GetCurrentStepNumber();

  G4double edepW = 0;                                 // weighted energy
  G4double tHit = aStep->GetTrack()->GetGlobalTime(); // global time
  G4StepPoint*  thePrePoint         = aStep->GetPreStepPoint();
  G4StepPoint*  thePostPoint        = aStep->GetPostStepPoint();
  G4ThreeVector hitPos =
    thePrePoint->GetPosition() + thePostPoint->GetPosition();
  hitPos /= 2.;
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* volume = theTouchable->GetVolume(1);
  G4int id;
  // Get element copy number from the mother volume (body of AT cell)
  id = volume->GetCopyNo() - 43; // 44 = copy # He inside teflon
  // Calculate the light collection efficiency at the current hit position
  // Use 3-way linear interpolation between points on the voxel map
  // contained with fLightColl
  G4double xHit[3];            // coordinate vector for interpolator
  // Z coordinate offsets for the cells of the active target
  //static G4double zLim[] = {-200,-100,0,100,200};
  static G4double zLim[] = {-200,-133.33,-66.66,0,66.66,133.33,200};
  G4double z = hitPos.getZ();
  //G4int ii = (G4int)(z/100) + 2;
  for(G4int i=0; i<=6; i++){
    if( i == 6 )
      return false;
    if( (z >= zLim[i]) && (z < zLim[i+1]) ){
      id = i;
      break;
    }
  }  
  //G4double z = hitPos.getZ() + zOffset[id];
  // Check if z within interpolation limits
  if( fVoxGridExtend ){
    G4double diff;
    if( z < fGridLim[4] ){
      diff = fGridLim[4] - z;
      if( diff < fVoxGridExtend ) z = fGridLim[4];
    }
    else if( z > fGridLim[5] ){
      diff = z - fGridLim[5];
      if( diff < fVoxGridExtend ) z = fGridLim[5];
    }
  } 
  xHit[2] = z;
  xHit[1] = hitPos.getY();
  xHit[0] = hitPos.getX();
  //
  // Light collection efficiency into each PMT
  Double_t lEff[4] = {0.25,0.25,0.25,0.25};
  Int_t istop;                       // last PMT indices for this step
  //if( id < 4 ) istop = 4;            // main cells
  //else istop = 1;                    // end cells
  istop = 4;
  // Now get the light collection factors
  //  for(G4int i=0; i<istop; i++){
  //lEff[i] = fLightColl->LinInterpN(xHit,i);
  //}   
  GasHit* Hit;
  // Is this step connected with an existing track
  if (fHitID[id]==-1){
    Hit = new GasHit;                            // no its a new one
    Hit->SetPos(hitPos);                         // primary position
    Hit->SetID(id);                              // cell ID
    Hit->SetTime(tHit);                          // primary time
    fHitID[id] = fCollection->insert(Hit) -1;    // show track started
    fHits[fNhits]=id;                            // cell hits array
    fNhits++;                                    // no. hits
    for(Int_t i=0; i<istop; i++) Hit->SetLC(lEff[i],i); // Light Coll @1st step
  }
// Not a new track, accumulate pulse heights
  else Hit =  (*fCollection)[fHitID[id]];
  Hit->AddEnergy(edep);
  for(G4int i=0; i<istop; i++){
    edepW = edep*lEff[i];
    Hit->AddEnergyPMT(edepW,i);
  }   
  return true;
}

//----------------------------------------------------------------------------
void GasSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  //  GasHit* Hit;
  //  for(G4int i=0; i<fNhits; i++){
  //    Hit =  (*fCollection)[fHitID[fHits[i]]];
  //  }
  if(fHCID<0)
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  if(fNhits>0)  HCE->AddHitsCollection(fHCID,fCollection);
  for (G4int i=0;i<fNhits;i++) {
    fHitID[fHits[i]]=-1;
    fHits[i]=0;
  }
  fNhits=0;
}

//-----------------------------------------------------------------------------
void GasSD::clear(){
} 

//-----------------------------------------------------------------------------
void GasSD::DrawAll(){
} 

//-----------------------------------------------------------------------------
void GasSD::PrintAll(){
} 

