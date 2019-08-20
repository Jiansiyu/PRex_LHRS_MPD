#ifndef ROOT_TreeSearch_PREXStand
#define ROOT_TreeSearch_PREXStand

#include "THaSpectrometer.h"

class PREXStand : public THaSpectrometer {

public:
    PREXStand( const char *name, const char *description );
    virtual ~PREXStand();

    virtual Int_t FindVertices( TClonesArray& tracks );
    virtual Int_t TrackCalc();


    virtual Int_t  ReadDatabase( const TDatime& date );
    void CalcTargetCoords(THaTrack *the_track );
    void CalcFocalPlaneCoords( THaTrack* track );

    //_________________________________
    // used for the database
    enum { kPORDER = 7 };
    // Class for storing matrix element data
    class THaMatrixElement {
    public:
      THaMatrixElement() : iszero(true), order(0), v(0)
        { pw.reserve(5); poly.reserve(kPORDER); }
      bool match( const THaMatrixElement& rhs ) const
        { assert(pw.size() == rhs.pw.size()); return ( pw == rhs.pw ); }
      void clear()
        { iszero = true; pw.clear(); order = 0; v = 0.0; poly.clear(); }

      bool iszero;             // whether the element is zero
      std::vector<int> pw;     // exponents of matrix element
                               //   e.g. D100 = { 1, 0, 0 }
      int  order;
      double v;                // its computed value
      std::vector<double> poly;// the associated polynomial
    };

protected:
    enum ECoordType { kTransport, kRotatingTransport };
    enum EFPMatrixElemTag { T000 = 0, Y000, P000 };
    // Optics matrix elements (FIXME: move to HRS)
    std::vector<THaMatrixElement> fTMatrixElems;
    std::vector<THaMatrixElement> fDMatrixElems;
    std::vector<THaMatrixElement> fPMatrixElems;
    std::vector<THaMatrixElement> fPTAMatrixElems; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fYMatrixElems;
    std::vector<THaMatrixElement> fYTAMatrixElems; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fFPMatrixElems;  // matrix elements used in
                                              // focal plane transformations
                                              // { T, Y, P }

    std::vector<THaMatrixElement> fLMatrixElems;   // Path-length corrections (meters)


    protected:
    ClassDef(PREXStand,0) // BigBite spectrometer
};

#endif


