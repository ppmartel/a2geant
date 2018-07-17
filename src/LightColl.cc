// AT Simulation of He Gas Scintillator Active Target
// J.R.M Annand, University of Glasgow
// Class LightColl
// Linear interpolate light collection efficiency at point x0,x1,x2
// from 3D voxel map
// 24/01/13 JRMA adapted from AcquMC original

// Interpolation method:
// Multidimensional Interpolation Methods
// InterpMethods.pdf version 04/09/2006, K.C. Johnson

#include "LightColl.hh"

//-----------------------------------------------------------------------------
LightColl::LightColl( G4int n, G4int idens )
{
  // Initial setup of N-dim generator for TFoam
  fNDim = n;                          // # dimensions
  fNInterp = 1;                       // default 1 variable interpolated
  //  fScale = scale;                 // scaling factors for x[j]
  fM0 = new G4int[fNDim];             // nearest-&-less-than grid indices
  fM1 = new G4int[fNDim];             // fM1[j] = fM0[j] + s, s=0,1
  fM2 = new G4int[fNDim];             // for linear interpolation compensation
  fIwgt = new G4int[fNDim];           // weighting of interp output
  for( G4int i=0; i<fNDim; i++ ) fIwgt[i] = EWgtUnity; // default weighting off
  fUj = new G4double[fNDim];          // Linear "Cardinal" functions
  fXd = new G4double[fNDim];          // x[j] - x[fM0[j]]
  //  fXscaled = new G4double[fNDim]; // scaled variable vector
  fXarg = new G4double[fNDim];        // vector of variables [0,1] interval
  fDensOpt = idens;                   // evaluation option
  fXN = new G4double*[fNDim];
  fSN = new G4int[fNDim];
  fIsInit = false;
}


//---------------------------------------------------------------------------
LightColl::~LightColl()
{
}

//-----------------------------------------------------------------------------
G4double LightColl::LinInterpN( G4double* x, G4int n, Bool_t comp )
{
  // Multi-linear interpolation of N-dimensional function evaluated on
  // a regular grid of points.
  //  Multidimensional Interpolation Methods
  //  InterpMethods.pdf version 04/09/2006, K.C. Johnson
  // fYN    are the evaluations on the grid
  // fXN[j] is the array of grid point values for each dimension
  // fSN[j] is the number of grid points for each dimension
  // x[j]  is the array of values at which to evaluate the linear interpolation
  // N     is the number of dimensions

  fYN = fYpmt[n];
  G4double* xn;                    // grid values for particular dimension
  G4int j,m;
  // Determine the grid points which "enclose" the input value
  // for each dimension
  for( j=0; j<fNDim; j++ ){
    xn = fXN[j];
    for( m=0; m<fSN[j]; m++ ){
      if( x[j] < xn[m] ) break;
    }
    m--;
    if( (m < 0) || (m == fSN[j]) ){
      //printf("Variable value %g of dimension j = %d out of range\n",x[j], j );
      return 0.0;
    }
    fM0[j] = m;
    fXd[j] = (x[j] - xn[m])/(xn[m+1] - xn[m]);
  }
  G4double yFit = 0.0;
  LinIntpAccum(0,yFit,comp);
  return yFit;
}

//-----------------------------------------------------------------------------
void LightColl::ReadData( G4String sfname, G4String stname )
{
  // Read light collection voxel tree
  // Read n-dimensional grid of values of function Y = f(x0,x1,...xn)
  // Dimension should be 3, ie each point on grid is a position
  // Each dimension (x,y,z) should be 75
  // Total size 75*75*75
  //
  // delete any previous initialisation
  if( fIsInit ){
    for( Int_t i=0; i<fNDim; i++ )
      delete fXN[i];
    for( Int_t i=0; i<4; i++ )
      delete fYpmt[i];
  }
  Char_t* fname = (Char_t*)sfname.data();
  Char_t* tname = (Char_t*)stname.data();
  TFile* voxFile = new TFile(fname);
  TTree* voxTree = (TTree*)voxFile->Get(tname);
  Int_t nE = voxTree->GetEntries();
  printf(" Voxel tree entries: %d\n", nE);
  Double_t grid = TMath::Power((Double_t)nE,0.3333333333);
  Int_t nGrid = TMath::Nint(grid);
  printf(" Voxel grid dimension: %d\n", nGrid);
  for( Int_t i=0; i<fNDim; i++ ){
    fSN[i] = nGrid;
    fXN[i] = new G4double[nGrid];
  }
  for( Int_t i=0; i<4; i++ )
    fYpmt[i] = new G4double[nE];
  //
  Float_t x[3];
  Float_t xo[3];
  Int_t pmt[5];
  Double_t leff;
  voxTree->SetBranchAddress("Z",x+2);
  voxTree->SetBranchAddress("Y",x+1);
  voxTree->SetBranchAddress("X",x);
  voxTree->SetBranchAddress("Ntot",pmt);
  voxTree->SetBranchAddress("Ndet4",pmt+1);
  voxTree->SetBranchAddress("Ndet5",pmt+2);
  voxTree->SetBranchAddress("Ndet6",pmt+3);
  voxTree->SetBranchAddress("Ndet7",pmt+4);
  Int_t ii[3];
  xo[0] = xo[1] = xo[2] = 0;
  ii[0] = ii[1] = ii[2] = 0;
  Double_t* xxf;
  Double_t* ypmt;
  Int_t iii;
  //  for(Int_t i=0; i<3; i++)
  //    xf[i] = new G4double[nGrid];
  for(Int_t i=0; i<nE; i++){
    voxTree->GetEntry(i);
    for(Int_t j=0; j<3; j++){
      if( (ii[j] < nGrid) && (x[j] != xo[j]) ){
	xxf = fXN[j];
	iii = ii[j];
	xxf[iii] = x[j];
	ii[j]++;
	xo[j] = x[j];
      }
    }
    for(Int_t j=0; j<4; j++){
      if(pmt[0])
	leff = (Double_t)pmt[j+1]/pmt[0];
      else leff = 0;
      ypmt = fYpmt[j];
      ypmt[i] = leff;
    }
  }
  Int_t ng = nGrid-1;
  printf(" Min and max positions on grid:\n  x:%g %g  y:%g %g  z:%g %g\n",
	 fXN[0][0],fXN[0][ng],fXN[1][0],fXN[1][ng],fXN[2][0],fXN[2][ng]);
  //
  fIsInit = true;
  delete voxTree;
  voxFile->Close();
  delete voxFile;
  return;
}

//-----------------------------------------------------------------------------
void LightColl::GetGridLimits( G4double* glim)
{
  // supply limits of grid values for each dimension
  Int_t k = 0;
  for(Int_t i=0; i<fNDim; i++){
    Int_t ng = fSN[i] - 1;
    glim[k] = fXN[i][0];
    glim[k+1] = fXN[i][ng];
    k += 2;
  }
}
