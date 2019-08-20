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
	std::cout<<"In transport Plane ("<<(track->GetX())<<","<<(track->GetY())<<std::endl;

}

//_____________________________________________________________________________
void PREXStand::CalcFocalPlaneCoords( THaTrack* track ){

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
	  //std::cout<<"fPrefix"<<fPrefix<<std::endl;
	  //Int_t err = LoadDB( file, date, request1, fPrefix );
	  /*
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
*/


	return kOK;
}

