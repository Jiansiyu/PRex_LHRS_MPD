//////////////////////////////////////////////////////////////////////////
//
// Bare-bones SBB BigBite spectrometer class 
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "PREXStand.h"
#include "THaTrack.h"
#include <map>
#include <string.h>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <cctype>
#include <cassert>
#include <vector>
// for test perpose
#include "iostream"

using namespace std;

ClassImp(PREXStand)

// Helper structure for parsing tensor data
typedef vector<PREXStand::THaMatrixElement> MEvec_t;
struct MEdef_t {
  MEdef_t() : npow(0), elems(0), isfp(false), fpidx(0) {}
  MEdef_t( Int_t npw, MEvec_t* elemp, Bool_t is_fp = false, Int_t fp_idx = 0 )
    : npow(npw), elems(elemp), isfp(is_fp), fpidx(fp_idx) {}
  MEvec_t::size_type npow; // Number of exponents for this element type
  MEvec_t* elems;          // Pointer to member variable holding data
  Bool_t isfp;             // This defines a focal plane matrix element
  Int_t fpidx;             // Index into fFPMatrixElems
};

//_____________________________________________________________________________
PREXStand::PREXStand( const char* name, const char* description ) :
  THaSpectrometer( name, description )
{
  // Constructor. Defines standard detectors

}

//_____________________________________________________________________________
PREXStand::~PREXStand()
{
  // Destructor
}

//_____________________________________________________________________________
Int_t PREXStand::FindVertices( TClonesArray& tracks/* tracks */ )
{
  // Reconstruct target coordinates for all tracks found.
  // TODO
  //Basically adapt the same class in the THaVDV.c for GEM detector tracking calculation
  // Siyu Added 2019-08-19
  Int_t n_exist = tracks.GetLast()+1;
  for(Int_t t=0; t < n_exist; t++){
	  THaTrack * theTrack = static_cast<THaTrack*>( tracks.At(t) );
	  CalcTargetCoords(theTrack);
  }
  if(n_exist!=0){
   fGoldenTrack= static_cast<THaTrack*>(tracks.At(0));
   fTrkIfo = *fGoldenTrack;
   fTrk    = fGoldenTrack;
  }else{
  fGoldenTrack=NULL;
  }
  return 0;
}

//_____________________________________________________________________________
void PREXStand::CalcTargetCoords(THaTrack *track){
	// Since the GEM detector already in Transport Coordination System, it does not need to transform from DCS to TCS
	// This is slightly different from the VDC
	// The CalcFocalPlaneCoord in the VDC is calculated in ConstructTracks,
	// It also would be fine to calculate here ?
	//  -- Siyu

	// Test the data in the buffer
// std::cout<<"In transport Plane ("<<(track->GetX())<<","<<(track->GetY())<<std::endl;
	
	// calculate the Focal plane coordination variables
 CalcFocalPlaneCoords(track);
	
	
  //Project the Focal Plane to the Target
  // calculates target coordinates from focal plane coordinates

  const Int_t kNUM_PRECOMP_POW = 10;

  Double_t x_fp, y_fp, th_fp, ph_fp;
  Double_t powers[kNUM_PRECOMP_POW][5];  // {(x), th, y, ph, abs(th) }
  Double_t x, y, theta, phi, dp, p, pathl;

  // first select the coords to use
  if( fCoordType == kTransport ) {
    x_fp = track->GetX();
    y_fp = track->GetY();
    th_fp = track->GetTheta();
    ph_fp = track->GetPhi();
  } else {  // kRotatingTransport
    x_fp = track->GetRX();
    y_fp = track->GetRY();
    th_fp = track->GetRTheta();
    ph_fp = track->GetRPhi();
  }
  
   // calculate the powers we need
  for(int i=0; i<kNUM_PRECOMP_POW; i++) {
    powers[i][0] = pow(x_fp, i);
    powers[i][1] = pow(th_fp, i);
    powers[i][2] = pow(y_fp, i);
    powers[i][3] = pow(ph_fp, i);
    powers[i][4] = pow(TMath::Abs(th_fp),i);
  }

  // calculate the matrices we need
  CalcMatrix(x_fp, fDMatrixElems);
  CalcMatrix(x_fp, fTMatrixElems);
  CalcMatrix(x_fp, fYMatrixElems);
  CalcMatrix(x_fp, fYTAMatrixElems);
  CalcMatrix(x_fp, fPMatrixElems);
  CalcMatrix(x_fp, fPTAMatrixElems);
  
  // calculate the coordinates at the target
  theta = CalcTargetVar(fTMatrixElems, powers);
  phi = CalcTargetVar(fPMatrixElems, powers)+CalcTargetVar(fPTAMatrixElems,powers);
  y = CalcTargetVar(fYMatrixElems, powers)+CalcTargetVar(fYTAMatrixElems,powers);

  //THaSpectrometer *app = static_cast<THaSpectrometer*>(GetApparatus());
  // calculate momentum
  dp = CalcTargetVar(fDMatrixElems, powers);
  p  = GetPcentral() * (1.0+dp);

  // pathlength matrix is for the Transport coord plane
  //pathl = CalcTarget2FPLen(fLMatrixElems, powers);

  //FIXME: estimate x ??
  x = 0.0;

  // Save the target quantities with the tracks
  track->SetTarget(x, y, theta, phi);
  track->SetDp(dp);
  track->SetMomentum(p);
  track->SetPathLen(pathl);
	
}
//_____________________________________________________________________________
Double_t PREXStand::CalcTargetVar(const vector<THaMatrixElement>& matrix,
                               const Double_t powers[][5])
{
  // calculates the value of a variable at the target
  // the x-dependence is already in the matrix, so only 1-3 (or np) used
  Double_t retval=0.0;
  Double_t v=0;
  for( vector<THaMatrixElement>::const_iterator it=matrix.begin();
       it!=matrix.end(); ++it )
    if(it->v != 0.0) {
      v = it->v;
      unsigned int np = it->pw.size(); // generalize for extra matrix elems.
      for (unsigned int i=0; i<np; ++i)
        v *= powers[it->pw[i]][i+1];
      retval += v;
  //      retval += it->v * powers[it->pw[0]][1]
  //                  * powers[it->pw[1]][2]
  //                  * powers[it->pw[2]][3];
    }

  return retval;
}

//_____________________________________________________________________________
void PREXStand::CalcMatrix(const Double_t x, vector<THaMatrixElement>& matrix){
  // calculates the values of the matrix elements for a given location
  // by evaluating a polynomial in x of order it->order with
  // coefficients given by it->poly

  for( vector<THaMatrixElement>::iterator it=matrix.begin();
       it!=matrix.end(); ++it ) {
    it->v = 0.0;

    if(it->order > 0) {
      for(int i=it->order-1; i>=1; --i)
        it->v = x * (it->v + it->poly[i]);
      it->v += it->poly[0];
    }
  }
}
//_____________________________________________________________________________
void PREXStand::CalcFocalPlaneCoords( THaTrack* track ){
	// for GEM, it is already in Transport Coordination System
	// Start from Transpot Coordination System, Project to the Focal plane
	// --Siyu 
	
	// read the Transport Coordination System Value
	double_t theta=track->GetTheta();
	double_t phi=track->GetPhi();
	double_t x=track->GetX();
	double_t y=track->GetY();
	
	// then calculate the rotating transport frame coordinates
	Double_t r_x = x;
	//
	CalcMatrix(r_x, fFPMatrixElems);
	Double_t r_y = y - fFPMatrixElems[Y000].v;  // Y000
	
	// TODO Questions ???
	// Calculate now the tan(rho) and cos(rho) of the local rotation angle.
    Double_t tan_rho_loc = fFPMatrixElems[T000].v;   // T000
    Double_t cos_rho_loc = 1.0/sqrt(1.0+tan_rho_loc*tan_rho_loc);
    
    // calculate the detector plan coordinations according to the Transport coordination
    double_t d_theta= (theta-tan_rho_loc)/(1.0 + theta*tan_rho_loc);
    double_t d_phi  = theta * cos_rho_loc *(1.00- d_theta*tan_rho_loc );

    Double_t r_phi   = (d_phi - fFPMatrixElems[P000].v /* P000 */)/
    	    (1.0-d_theta*tan_rho_loc) / cos_rho_loc;
    Double_t r_theta = (d_theta+tan_rho_loc) /
    	    (1.0-d_theta*tan_rho_loc);
    
    track->SetR(r_x, r_y, r_theta, r_phi);  // the roration transport
}

//_____________________________________________________________________________
Int_t PREXStand::TrackCalc()
{
  // Additioal track calculations

  // TODO

  return 0;
}

//____________________________________________________________________________
static Int_t ParseMatrixElements( const std::string& MEstring,
		                          std::map<std::string,MEdef_t>& matrix_map,
                                  const char* prefix ){
	// Parse the contents of MEstring (from the database) into the local
	// member variables holding the matrix elements
	const char* const here = "PREXStand::ParseMatrixElements";
	//TODO
	istringstream ist(MEstring.c_str());
	std::string word, w;
	bool findnext = true, findpowers = true;
	  Int_t powers_left = 0;
	  map<string,MEdef_t>::iterator cur = matrix_map.end();
	  PREXStand::THaMatrixElement ME;
	  while( ist >> word ) {
	    if( !findnext ) {
	      assert( cur != matrix_map.end() );
	      bool havenext = isalpha(word[0]);
	      if( findpowers ) {
	        assert( powers_left > 0 );
	        if( havenext || word.find_first_not_of("0123456789") != string::npos ||
	            atoi(word.c_str()) > 9 ) {
	          Error( Here(here,prefix), "Bad exponent = %s for matrix element \"%s\". "
	                 "Must be integer between 0 and 9. Fix database.",
	                 word.c_str(), w.c_str() );
	          return THaAnalysisObject::kInitError;
	        }
	        ME.pw.push_back( atoi(word.c_str()) );
	        if( --powers_left == 0 ) {
	          // Read all exponents
	          if( cur->second.isfp ) {
	            // Currently only the "000" focal plane matrix elements are supported
	            if( ME.pw[0] != 0 || ME.pw[1] != 0 || ME.pw[2] != 0 ) {
	              Error( Here(here,prefix), "Bad coefficients of focal plane matrix "
	                    "element %s = %d %d %d. Fix database.",
	                    w.c_str(), ME.pw[0], ME.pw[1], ME.pw[2] );
	              return THaAnalysisObject::kInitError;
	            } else {
	              findpowers = false;
	            }
	          }
	          else {
	            // Check if this element already exists, if so, skip
	            MEvec_t* mat = cur->second.elems;
	            assert(mat);
	            bool match = false;
	            for( MEvec_t::iterator it = mat->begin();
	                 it != mat->end() && !(match = it->match(ME)); ++it ) {}
	            if( match ) {
	              Warning( Here(here,prefix), "Duplicate definition of matrix element %s. "
	                      "Using first definition.", cur->first.c_str() );
	              findnext = true;
	            } else {
	              findpowers = false;
	            }
	          }
	        }
	      } else {
	        if( !havenext ) {
	          if( ME.poly.size() >= PREXStand::kPORDER )
	            havenext = true;
	          else {
	            ME.poly.push_back( atof(word.c_str()) );
	            if( ME.poly.back() != 0.0 ) {
	              ME.iszero = false;
	              ME.order = ME.poly.size();
	            }
	          }
	        }
	        if( havenext || ist.eof() ) {
	          if( ME.poly.empty() ) {
	            // No data read at all?
	            Error( Here(here,prefix), "Could not read in Matrix Element %s%d%d%d!",
	                  w.c_str(), ME.pw[0], ME.pw[1], ME.pw[2]);
	            return THaAnalysisObject::kInitError;
	          }
	          if( !ME.iszero ) {
	            MEvec_t* mat = cur->second.elems;
	            assert(mat);
	            // The focal plane matrix elements are stored in a single vector
	            if( cur->second.isfp ) {
	              PREXStand::THaMatrixElement& m = (*mat)[cur->second.fpidx];
	              if( m.order > 0 ) {
	                Warning( Here(here,prefix), "Duplicate definition of focal plane "
	                        "matrix element %s. Using first definition.",
	                        w.c_str() );
	              } else
	                m = ME;
	            } else
	              mat->push_back(ME);
	          }
	          findnext = true;
	          if( !havenext )
	            continue;
	        } // if (havenext)
	      } // if (findpowers) else
	    } // if (!findnext)
	    if( findnext ) {
	      cur = matrix_map.find(word);
	      if( cur == matrix_map.end() ) {
	        // Error( Here(here,prefix), "Unknown matrix element type %s. Fix database.",
	        //       word.c_str() );
	        // return THaAnalysisObject::kInitError;
	        continue;
	      }
	      ME.clear();
	      findnext = false; findpowers = true;
	      powers_left = cur->second.npow;
	      w = word;
	    }
	  } // while(word)

	  return THaAnalysisObject::kOK;

}


//_____________________________________________________________________________
Int_t PREXStand::ReadDatabase(const TDatime &date){

	// load the database for the GEM correction
	// TODO
	// add the rotation and correction matrix for the GEM detector
	const char * const here="ReadDatabase";

	//std::cout<<"\n\n\n This is a test for read the database \n\n\n"<<date.Print()<<std::endl;
	FILE* file = OpenFile( date );
	if( !file ) return kFileError;

// for GEM it does not needed
// Read fOrigin and fSize (currently unused)
//	Int_t err = ReadGeometry( file, date );
//	if( err ) {
//	   fclose(file);
//	   return err;
//	 }

	  // Read TRANSPORT matrices
	  //FIXME: move to HRS
	  fTMatrixElems.clear();
	  fDMatrixElems.clear();
	  fPMatrixElems.clear();
	  fPTAMatrixElems.clear();
	  fYMatrixElems.clear();
	  fYTAMatrixElems.clear();
	  fLMatrixElems.clear();

	  fFPMatrixElems.clear();
	  fFPMatrixElems.resize(3);

	  // For GEM not all of them are needed, will need to delete some unused ones
	  //TODO
	  std::map<std::string, MEdef_t> matrix_map;
	  //RANSPORT to focal-plane tensors
	  matrix_map["t"]   = MEdef_t( 3, &fFPMatrixElems, true, 0 );
	  matrix_map["y"]   = MEdef_t( 3, &fFPMatrixElems, true, 1 );
	  matrix_map["p"]   = MEdef_t( 3, &fFPMatrixElems, true, 2 );
	  // Standard focal-plane to target matrix elements (D=delta, T=theta, Y=y, P=phi)
	  matrix_map["D"]   = MEdef_t( 3, &fDMatrixElems );
	  matrix_map["T"]   = MEdef_t( 3, &fTMatrixElems );
	  matrix_map["Y"]   = MEdef_t( 3, &fYMatrixElems );
	  matrix_map["P"]   = MEdef_t( 3, &fPMatrixElems );
	  // Additional matrix elements describing the dependence of y-target and
	  // phi-target on the /absolute value/ of theta, found necessary in optimizing
	  // the septum magnet optics (R. Feuerbach, March 1, 2005)
	  matrix_map["YTA"] = MEdef_t( 4, &fYTAMatrixElems );
	  matrix_map["PTA"] = MEdef_t( 4, &fPTAMatrixElems );
	  // Matrix for calculating pathlength from z=0 (target) to focal plane (meters)
	  // (R. Feuerbach, October 16, 2003)
	  matrix_map["L"]   = MEdef_t( 4, &fLMatrixElems );

	  std::string MEstring;
	  DBRequest request1[] = {
	    { "matrixelem",  &MEstring, kString },
	    { 0 }
	  };
	 std::cout<<"fPrefix"<<fPrefix<<std::endl;
	 Int_t err = LoadDB( file, date, request1, fPrefix );
	  
	  if( err ) {
	    fclose(file);
	    return err;
	  }
	  if( MEstring.empty() ) {
	    Error( Here(here), "No matrix elements defined. Set \"maxtrixelem\" in database." );
	    fclose(file);
	    return kInitError;
	  }
	  // Parse the matrix elements
	  err = ParseMatrixElements( MEstring, matrix_map, fPrefix );
	  if( err ) {
	    fclose(file);
	    return err;
	  }
	  MEstring.clear();

	  // Ensure that we have all three focal plane matrix elements, else we cannot
	  // do anything sensible with the tracks
	  if( fFPMatrixElems[T000].order == 0 ) {
	    Error( Here(here), "Missing FP matrix element t000. Fix database." );
	    err = kInitError;
	  }
	  if( fFPMatrixElems[Y000].order == 0 ) {
	    Error( Here(here), "Missing FP matrix element y000. Fix database." );
	    err = kInitError;
	  }
	  if( fFPMatrixElems[P000].order == 0 ) {
	    Error( Here(here), "Missing FP matrix element p000. Fix database." );
	    err = kInitError;
	  }
	  if( err ) {
	    fclose(file);
	    return err;
	  }
	return kOK;
}

