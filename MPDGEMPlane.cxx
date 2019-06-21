#include <iostream>
#include "MPDGEMPlane.h"
#include "TDatime.h"
#include "THaEvData.h"
#include <vector>
#include <numeric>

#include "GEMHit.h"


// for debug
#include <TCanvas.h>
#include <TH1F.h>
#include <map>
#include <TStyle.h>

#define ALL(c) (c).begin(), (c).end()

MPDGEMPlane::MPDGEMPlane( const char *name, const char *description,
    THaDetectorBase* parent ):
    GEMPlane(name,description,parent),
//    fNch(0),fStrip(NULL),fPedestal(NULL),fcommon_mode(NULL),
    fDefaultRMS(50.0), fNch(0), fStrip(0), 
    fADC0(0), fADC1(0), fADC2(0), fADC3(0), fADC4(0), fADC5(0), fADCSum(0),
    frADC0(0), frADC1(0), frADC2(0), frADC3(0), frADC4(0), frADC5(0),
    trigger_time(-1),ev_num(-1)

{
    fZeroSuppress    = kTRUE;
    fZeroSuppressRMS = 5.0;

    fMaxSamp = N_MPD_TIME_SAMP;

    for( UInt_t i = 0; i < fMaxSamp; i++ ){
      fADCForm[i] = NULL;
      frawADC[i] = NULL;
    }


    return;
}

MPDGEMPlane::~MPDGEMPlane() {
    if( fSigStrips.size() > 0 ){
        fADC0 = NULL;
        fADC1 = NULL;
        fADC2 = NULL;
        fADC3 = NULL;
        fADC4 = NULL;
        fADC5 = NULL;
	fADCSum = NULL;
        
	frADC0 = NULL;
        frADC1 = NULL;
        frADC2 = NULL;
        frADC3 = NULL;
        frADC4 = NULL;
        frADC5 = NULL;
        for( UInt_t i = 0; i < fMaxSamp; i++ ){
            delete fADCForm[i];
            fADCForm[i] = NULL;
            
	    delete frawADC[i];
            frawADC[i] = NULL;
        }
        delete fStrip;
        fStrip = NULL;
        /*
	delete fcommon_mode;
	fcommon_mode = NULL;
        delete fPedestal;
        fPedestal = NULL;
        */
	
    }

    return;
}

Int_t MPDGEMPlane::ReadDatabase( const TDatime& date ){
    std::cout << "[MPDGEMPlane::ReadDatabase]" << std::endl;

    // Read the database for the base class, quit if error

    Int_t status = ReadDatabaseCommon(date);
    if( status != kOK )
        return status;

    Int_t err;

    FILE* file = OpenFile( date );
    if( !file ) return kFileError;

    TString mapping;
    fMaxClusterSize = kMaxUInt;
    fMinAmpl   = 0.0;
    fSplitFrac = 0.0;
    fMapType   = kOneToOne;
    fMaxSamp   = 6;
    fChanMap.clear();
    fPed.clear();
    fAmplSigma = 0.36; // default, an educated guess


    std::vector<Double_t> rawped;
    std::vector<Double_t> rawrms;

    Int_t gbl = GetDBSearchLevel(fPrefix);

    Int_t do_noise = 1;

    const DBRequest request[] = {
        { "chanmap",        &fChanMapData,  kIntV,    0, 0},
        { "ped",            &rawped,        kDoubleV, 0, 1},
        { "rms",            &rawrms,        kDoubleV, 0, 1},

        { "strip.pos",      &fStart },
        { "strip.pitch",    &fPitch,          kDouble,  0, 0, gbl },
        { "maxclustsiz",    &fMaxClusterSize, kUInt,    0, 1, gbl },
        { "maxsamp",        &fMaxSamp,        kUInt,    0, 1, gbl },
        { "adc.min",        &fMinAmpl,        kDouble,  0, 1, gbl },
        { "split.frac",     &fSplitFrac,      kDouble,  0, 1, gbl },
        { "mapping",        &mapping,         kTString, 0, 1, gbl },
        { "do_noise",       &do_noise,        kInt,     0, 1, gbl },


        {}
    };
    err = LoadDB( file, date, request, fPrefix );

    if( err ) return err;

    fclose(file);

    SetBit( kDoNoise, do_noise );


    UInt_t nentry = fChanMapData.size()/MPDMAP_ROW_SIZE;
    for( UInt_t mapline = 0; mapline < nentry; mapline++ ){
        mpdmap_t thisdata;
        thisdata.crate  = fChanMapData[0+mapline*MPDMAP_ROW_SIZE];
        thisdata.slot   = fChanMapData[1+mapline*MPDMAP_ROW_SIZE];
        thisdata.mpd_id = fChanMapData[2+mapline*MPDMAP_ROW_SIZE];
        thisdata.gem_id = fChanMapData[3+mapline*MPDMAP_ROW_SIZE];
        thisdata.adc_id = fChanMapData[4+mapline*MPDMAP_ROW_SIZE];
        thisdata.i2c    = fChanMapData[5+mapline*MPDMAP_ROW_SIZE];
        thisdata.pos    = fChanMapData[6+mapline*MPDMAP_ROW_SIZE];
        thisdata.invert = fChanMapData[7+mapline*MPDMAP_ROW_SIZE];
        fMPDmap.push_back(thisdata);
    }

    fNelem = N_APV25_CHAN*fMPDmap.size();

    SafeDelete(fADCraw);
    SafeDelete(fADC);
    SafeDelete(fHitTime);
    SafeDelete(fADCcor);
    SafeDelete(fGoodHit);
    SafeDelete(fADCSum);

    std::cout << fName << " mapped to " << nentry << " APV25 chips" << std::endl;

    for( UInt_t i = 0; i < fMaxSamp; i++ ){
        fADCForm[i] = new Int_t [fNelem];
        frawADC[i] = new Int_t [fNelem];
        for( UInt_t j = 0; j < fMaxSamp; j++ ){
            fADCForm[i][j] = 0.0;
            frawADC[i][j] = 0.0;
        }
    }
    fADC0 = fADCForm[0];
    fADC1 = fADCForm[1];
    fADC2 = fADCForm[2];
    fADC3 = fADCForm[3];
    fADC4 = fADCForm[4];
    fADC5 = fADCForm[5];

    frADC0 = frawADC[0];
    frADC1 = frawADC[1];
    frADC2 = frawADC[2];
    frADC3 = frawADC[3];
    frADC4 = frawADC[4];
    frADC5 = frawADC[5];

    fADCSum = new Float_t[fNelem];
    fStrip  = new Int_t[fNelem];
//    fcommon_mode = new Int_t[N_APV25_CHAN*nentry];

    fPed.clear();
    fPed.resize(fNelem, 0.0);

    fRMS.clear();
    fRMS.resize(fNelem, fDefaultRMS);

    for( UInt_t i = 0; i < rawped.size(); i++ ){
        if( (i % 2) == 1 ) continue;
        UInt_t idx = (UInt_t) rawped[i];
	
    	if( idx < (UInt_t) fNelem ){
    		fPed[idx] = rawped[i+1];
    	} else {
		
    	    std::cout << "[MPDGEMPlane::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
    	}
    }

    for( UInt_t i = 0; i < rawrms.size(); i++ ){
        if( (i % 2) == 1 ) continue;
        UInt_t idx = (UInt_t) rawrms[i];

	if( idx < ((UInt_t) fNelem) ){
            if(  fRMS[idx] > 0 ){
		fRMS[idx] = rawrms[i+1];
            } else {
                std::cout << "[MPDGEMPlane::ReadDatabase]  WARNING: " << " strip " << idx  << " has invalid RMS: " << fRMS[idx] << std::endl;
            }
	} else {
	    std::cout << "[MPDGEMPlane::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
	}
    }

    fADCraw = new Float_t[fNelem];
    fADC = new Float_t[fNelem];
    fHitTime = new Float_t[fNelem];
    fADCcor = new Float_t[fNelem];
    fGoodHit = new Byte_t[fNelem];
    fSigStrips.reserve(fNelem);
    fStripsSeen.resize(fNelem);

    fIsInit = true;
    return 0;
}

Int_t MPDGEMPlane::DefineVariables( EMode mode ) {
    if( mode == kDefine and fIsSetup ) return kOK;
      fIsSetup = ( mode == kDefine );

      RVarDef vars[] = {
          { "nrawstrips",     "nstrips with decoder data",        "fNrawStrips" },
          { "nhitstrips",     "nstrips > 0",                      "fNhitStrips" },
          { "nstrips",        "Num strips with hits > adc.min",   "GetNsigStrips()" },
          { "hitocc",         "strips > 0 / n_all_strips",        "fHitOcc" },
          { "occupancy",      "nstrips / n_all_strips",           "fOccupancy" },
          { "strip.adcraw",   "Raw strip ADC sum",                "fADCraw" },
          { "strip.adc",      "Deconvoluted strip ADC sum",       "fADC" },
          { "strip.adc_c",    "Pedestal-sub strip ADC sum",       "fADCcor" },
          { "strip.time",     "Leading time of strip signal (ns)","fHitTime" },
          { "strip.good",     "Good pulse shape on strip",        "fGoodHit" },
          { "strip.number", "Strip number mapping", "fSigStrips" },
          { "nch", "Number of channels", "fNch" },
          { "strip_number", "Strip Number", "fStrip" },
          { "adc0", "ADC sample", "fADC0" },
          { "adc1", "ADC sample", "fADC1" },
          { "adc2", "ADC sample", "fADC2" },
          { "adc3", "ADC sample", "fADC3" },
          { "adc4", "ADC sample", "fADC4" },
          { "adc5", "ADC sample", "fADC5" },
          { "radc0", "raw ADC sample", "frADC0" },
          { "radc1", "raw ADC sample", "frADC1" },
          { "radc2", "raw ADC sample", "frADC2" },
          { "radc3", "raw ADC sample", "frADC3" },
          { "radc4", "raw ADC sample", "frADC4" },
          { "radc5", "raw ADC sample", "frADC5" },
          { "adc_sum", "ADC samples sum", "fADCSum" },
          { "nhits",          "Num hits (clusters of strips)",    "GetNhits()" },
          { "noise",          "Noise level (avg below adc.min)",  "fDnoise" },
          { "ncoords",        "Num fit coords",                   "GetNcoords()" },
          { "coord.pos",      "Position used in fit (m)",         "fFitCoords.TreeSearch::FitCoord.fPos" },
          { "coord.trkpos",   "Track pos from projection fit (m)","fFitCoords.TreeSearch::FitCoord.fTrackPos" },
          { "coord.trkslope", "Track slope from projection fit",  "fFitCoords.TreeSearch::FitCoord.fTrackSlope" },
          { "coord.resid",    "Residual of trkpos (m)",           "fFitCoords.TreeSearch::FitCoord.GetResidual()" },
          { "coord.3Dpos",    "Crossing position of fitted 3D track (m)", "fFitCoords.TreeSearch::FitCoord.f3DTrkPos" },
          { "coord.3Dresid",  "Residual of 3D trkpos (m)",        "fFitCoords.TreeSearch::FitCoord.Get3DTrkResid()" },
          { "coord.3Dslope",  "Slope of fitted 3D track wrt projection",  "fFitCoords.TreeSearch::FitCoord.f3DTrkSlope" },


          //{ "fNch",   "Number of channels",   "fNch" },
          //{ "common_mode", "Common Mode", "fcommon_mode" },
          { "trigger_time", "Trigger Time", "trigger_time" },
          { "ev_num","event counter","ev_num"},
          { 0 },
      };


      Int_t ret = DefineVarsFromList( vars, mode );

      if( ret != kOK )
          return ret;


      RVarDef nonmcvars[] = {
          { "hit.pos",  "Hit centroid (m)",      "fHits.TreeSearch::GEMHit.fPos" },
          { "hit.adc",  "Hit ADC sum",           "fHits.TreeSearch::GEMHit.fADCsum" },
          { "hit.size", "Num strips ",           "fHits.TreeSearch::GEMHit.fSize" },
          { "hit.type", "Hit analysis result",   "fHits.TreeSearch::GEMHit.fType" },
          { 0 }
      };
      ret = DefineVarsFromList( nonmcvars, mode );


      return kOK;

}

void MPDGEMPlane::Clear( Option_t* opt ){

    fNch = 0;
    TreeSearch::GEMPlane::Clear(opt);
    return;
}

Int_t MPDGEMPlane::Decode( const THaEvData& evdata ){

    //std::cout << "[MPDGEMPlane::Decode " << fName << "]" << std::endl;

	// used for debug perpose
	std::map<int,float> _rawHistoBuf;
	std::map<int,float> _cmSupHistoBuf;
	std::map<int,float> _pedSupHistoBuf;
	std::map<int,float> _zeroSupHistoBuf;
	std::map<int,float> _pedestalHistoBuf;
	std::map<int,float> _rmsHistoBuf;


	//   Time sample   RStrips ADC value
	std::map<int, std::map<int, float>> _TSrawHistoBuf;
	std::map<int, std::map<int, float>> _TScmSupHistoBuf;
	std::map<int, std::map<int, float>> _TSpedSupHistoBuf;
	std::map<int, std::map<int, float>> _TSzeroSupHistoBuf;

	Bool_t _above_th=kFALSE;

	fNch = 0;
	for (std::vector<mpdmap_t>::iterator it = fMPDmap.begin();
			it != fMPDmap.end(); ++it) {
		// Find channel for trigger time first
		Int_t effChan = it->mpd_id << 5;  // Channel reserved for trigger time
		ULong_t coarse_time1 = evdata.GetData(it->crate, it->slot, effChan, 0);
		UInt_t coarse_time2 = evdata.GetData(it->crate, it->slot, effChan, 1);
		UInt_t fine_time = evdata.GetData(it->crate, it->slot, effChan, 2);
		trigger_time = ((coarse_time1 << 20) | coarse_time2) + fine_time / 6.0;

		effChan = it->mpd_id << 4;
		ev_num = evdata.GetData(it->crate, it->slot, effChan, 0);

		// Start reading data sample
		effChan = it->mpd_id << 8 | it->adc_id;
		// Find channel for this crate/slot

		Int_t fNchan = evdata.GetNumChan(it->crate, it->slot);


		for (Int_t ichan = 0; ichan < fNchan; ++ichan) {
			Int_t chan = evdata.GetNextChan(it->crate, it->slot, ichan);
			if (chan != effChan)
				continue; // not part of this detector

			UInt_t nsamp = evdata.GetNumHits(it->crate, it->slot, chan);

			assert(nsamp == N_APV25_CHAN*fMaxSamp);

			Double_t arrADCSum[128]; // Copy of ADC sum for CMN Calculation

			// Get the raw data
			Int_t fNchStartOfAPV = fNch; // count the number of fired strips

			Int_t rawADCBuff[N_MPD_TIME_SAMP][N_APV25_CHAN]; // used for temporary buffer the raw result

			// temporary usage
			Int_t cmSupADCBuff[N_MPD_TIME_SAMP][N_APV25_CHAN]; // used for temporary buffer the raw result

			Int_t ADCBuff[N_MPD_TIME_SAMP][N_APV25_CHAN]; // used for temporary buffer the result
			Int_t ADCBuffSum[N_APV25_CHAN];
			for (int j = 0; j < N_APV25_CHAN; j++) {
				for (int i = 0; i < fMaxSamp; i++) {
					rawADCBuff[i][j] = 0;
					ADCBuff[i][i] = 0;
					ADCBuffSum[j] = 0;
				}

			}
			for (uint8_t adc_samp = 0; adc_samp < fMaxSamp; adc_samp++) {
				// calculate the common mode
				// data is packed like this
				// [ts1 of 128 chan] [ts2 of 128chan] ... [ts6 of 128chan]
				std::vector<Int_t> rawBuf;
				for (Int_t strip = 0; strip < N_APV25_CHAN; ++strip) {
					UInt_t isamp = adc_samp * N_APV25_CHAN + strip;
					assert(isamp < nsamp);

					Int_t rawadc = evdata.GetData(it->crate, it->slot, chan,
							isamp);
					rawBuf.push_back(rawadc);
				}
				// calculate the CM
				std::sort(rawBuf.begin(), rawBuf.end());

				// get a sub-vector of the vector
				const std::vector<Int_t> cmVect(rawBuf.begin() + 20,
						rawBuf.end() - 20);

				auto tsCommonMode = std::accumulate(cmVect.begin(),
						cmVect.end(), 0.0) / cmVect.size();

				// subtract the common mode from the raw data
				for (Int_t strip = 0; strip < N_APV25_CHAN; ++strip) {
					UInt_t isamp = adc_samp * N_APV25_CHAN + strip;
					assert(isamp < nsamp);

					// get the phys pos
					Int_t RstripPos = GetRStripNumber(strip, it->pos, it->invert);

					Int_t rawadc = evdata.GetData(it->crate, it->slot, chan,
							isamp);

					rawADCBuff[adc_samp][strip]=rawadc;
					cmSupADCBuff[adc_samp][strip]=rawadc - tsCommonMode;
					ADCBuff[adc_samp][strip]=rawadc - tsCommonMode - fPed[RstripPos];
					ADCBuffSum[strip]=ADCBuffSum[strip]+ADCBuff[adc_samp][strip];

				} // end of loop on the 128 channels

			} // end of loop on the Time samples

			//zero suppression

			for(auto strip = 0; strip < N_APV25_CHAN; ++strip){
				Int_t RstripPos = GetRStripNumber(strip, it->pos, it->invert);

				Bool_t isAboveThreshold = (float)(ADCBuffSum[strip]/((Int_t)fMaxSamp)) > fZeroSuppressRMS * fRMS[RstripPos];


				// loop on the time samples
				for(auto adc_samp=0;adc_samp<fMaxSamp;adc_samp++){
					_TSrawHistoBuf[adc_samp][RstripPos]=rawADCBuff[adc_samp][strip];
					_TScmSupHistoBuf[adc_samp][RstripPos]=cmSupADCBuff[adc_samp][strip];
					_TSpedSupHistoBuf[adc_samp][RstripPos]=ADCBuff[adc_samp][strip];
					if(adc_samp==0){
						std::cout << " pedstal sub"<<_TSpedSupHistoBuf[0][RstripPos]<<std::endl;
					}

					// apply the 5-sigma cut
					if(isAboveThreshold){
						_TSzeroSupHistoBuf[adc_samp][RstripPos]=ADCBuff[adc_samp][strip];
					}

				}
				std::cout << "Check Above TH:  plane:" << fName.Data()
										<< " Pos:" << RstripPos << "  ADC raw value:"
										<< _TSrawHistoBuf[0][RstripPos]<<" commonsub:"
										<<_TScmSupHistoBuf[0][RstripPos]<<" pedestsup:"
										<<_TSpedSupHistoBuf[0][RstripPos]<<" pedestal:"
										<< fPed[RstripPos]
										<< "  rms: " << fRMS[RstripPos]
										<< std::endl;



				_rawHistoBuf[RstripPos]=(rawADCBuff[0][strip]+rawADCBuff[1][strip]+rawADCBuff[2][strip])/3.0;
				_cmSupHistoBuf[RstripPos]=(cmSupADCBuff[0][strip]+cmSupADCBuff[1][strip]+cmSupADCBuff[2][strip])/3.0;
				_pedSupHistoBuf[RstripPos]=(ADCBuff[0][strip]+ADCBuff[1][strip]+ADCBuff[2][strip])/3.0;
				_rmsHistoBuf[RstripPos]=fRMS[RstripPos];
				_pedestalHistoBuf[RstripPos]=fPed[RstripPos];


				if(isAboveThreshold){
					_zeroSupHistoBuf[RstripPos]=(ADCBuff[0][strip]+ADCBuff[1][strip]+ADCBuff[2][strip])/3.0;
					_above_th=kTRUE;
				}
				//std::cout<<"\n Check Above TH:  plane:"<< fName.Data()<<" Pos:"<< RstripPos <<"  ADC value:"<<ADCBuff[1][strip]<<"  rms: "<<fRMS[RstripPos]<<std::endl;


				Vflt_t samples;
				samples.clear();

				for (auto adc_samp=0; adc_samp<fMaxSamp;adc_samp++){

					frawADC[adc_samp][fNch]=rawADCBuff[adc_samp][strip];
					fADCForm[adc_samp][fNch]=ADCBuff[adc_samp][strip];

					samples.push_back((float_t)ADCBuff[adc_samp][strip]);

				}

				//Comments: change the value to be the net charge, need to confirm with Seamus
				//Does it need to be the 'signal', noise information does not have any meaning full information
				MPDStripData_t stripdata = ChargeDep(samples);
				++fNrawStrips;
				++fNhitStrips;

				fADCraw[RstripPos] = stripdata.adcraw;
				fADC[RstripPos] = stripdata.adc;
				fHitTime[RstripPos] = stripdata.time;
				fGoodHit[RstripPos] = stripdata.pass;

				fADCcor[RstripPos] = stripdata.adc;

				if (isAboveThreshold) {

					fSigStrips.push_back(RstripPos);
				}

				if(!fZeroSuppress || isAboveThreshold){
					fNch++;
				}
			}

        }// End ichan loop: fNchan = total APVs 

    }

    fHitOcc    = static_cast<Double_t>(fNhitStrips) / fNelem;
    fOccupancy = static_cast<Double_t>(GetNsigStrips()) / fNelem;


    // generate the plot and save it into jpg file
//    if(evdata.GetEvNum()>60)
	{

		// create the histogram
		std::map<int,TH1F *> rawHisto;
		std::map<int,TH1F *> cmSupHisto;
		std::map<int,TH1F *> pedSupHisto;
		std::map<int,TH1F *> zeroSupHisto;

		// Initialize all the histograms
		for (auto adc_samp = 0; adc_samp < fMaxSamp; adc_samp++) {
			rawHisto[adc_samp] = new TH1F(
					Form("Raw_histo_evt%d_%s_Tsample%d", evdata.GetEvNum(),
							fName.Data(), adc_samp),
					Form("Raw_histo_evt%d_%s_Tsample%d", evdata.GetEvNum(),
							fName.Data(), adc_samp), 2000, 0, 2000);
			rawHisto[adc_samp]->GetYaxis()->SetLabelSize(0.1);
			rawHisto[adc_samp]->GetXaxis()->SetLabelSize(0.1);
			rawHisto[adc_samp]->SetTitleSize(2);

			cmSupHisto[adc_samp] = new TH1F(
					Form("cmSup_histo_evt%d_%s_Tsample%d", evdata.GetEvNum(),
							fName.Data(), adc_samp),
					Form("cmSup_histo_evt%d_%s_Tsample%d", evdata.GetEvNum(),
							fName.Data(), adc_samp), 2000, 0, 2000);
			cmSupHisto[adc_samp]->GetYaxis()->SetLabelSize(0.1);
			cmSupHisto[adc_samp]->GetXaxis()->SetLabelSize(0.1);

			pedSupHisto[adc_samp] = new TH1F(
					Form("pedSup_histo_evt%d_%s_Tsample%d", evdata.GetEvNum(),
							fName.Data(), adc_samp),
					Form("pedSup_histo_evt%d_%s_Tsample%d", evdata.GetEvNum(),
							fName.Data(), adc_samp), 2000, 0, 2000);
			pedSupHisto[adc_samp]->GetYaxis()->SetLabelSize(0.1);
			pedSupHisto[adc_samp]->GetXaxis()->SetLabelSize(0.1);

//			pedSupHisto[adc_samp]->GetYaxis()->SetRangeUser(-50,4000);
			zeroSupHisto[adc_samp] = new TH1F(
					Form("zeroSup_histo_evt%d_%s_Tsample%d", evdata.GetEvNum(),
							fName.Data(), adc_samp),
					Form("zeroSup_histo_evt%d_%s_Tsample%d", evdata.GetEvNum(),
							fName.Data(), adc_samp), 2000, 0, 2000);
			zeroSupHisto[adc_samp]->GetYaxis()->SetLabelSize(0.1);
			zeroSupHisto[adc_samp]->GetXaxis()->SetLabelSize(0.1);

			gStyle->SetTitleFontSize(0.1);

		}

		for(auto adc_samp = 0; adc_samp<fMaxSamp;adc_samp++){
			rawHisto[adc_samp]->GetXaxis()->SetRangeUser(-10,
					_TSrawHistoBuf[adc_samp].size());
			cmSupHisto[adc_samp]->GetXaxis()->SetRangeUser(-10,
					_TSrawHistoBuf[adc_samp].size());
			pedSupHisto[adc_samp]->GetXaxis()->SetRangeUser(-10,
					_TSrawHistoBuf[adc_samp].size());
			zeroSupHisto[adc_samp]->GetXaxis()->SetRangeUser(-10,
					_TSrawHistoBuf[adc_samp].size());
		}
		// dump data to the histogram
		for(auto adc_samp = 0; adc_samp<fMaxSamp;adc_samp++){
			// dump the raw
			for(auto iter=_TSrawHistoBuf[adc_samp].begin();iter!=_TSrawHistoBuf[adc_samp].end();iter++){
				rawHisto[adc_samp]->Fill(iter->first,iter->second);
//				rawHisto[adc_samp]->GetXaxis()->SetRangeUser(-10,_TSrawHistoBuf[adc_samp].size())
			}
			// dump the cm
			for (auto iter = _TScmSupHistoBuf[adc_samp].begin();iter != _TScmSupHistoBuf[adc_samp].end(); iter++) {
				cmSupHisto[adc_samp]->Fill(iter->first, iter->second);
			}
			// dump the ped
			for (auto iter = _TSpedSupHistoBuf[adc_samp].begin();iter != _TSpedSupHistoBuf[adc_samp].end(); iter++) {
				pedSupHisto[adc_samp]->Fill(iter->first, iter->second);
			}
			// dump the zero
			for (auto iter = _TSzeroSupHistoBuf[adc_samp].begin();iter != _TSzeroSupHistoBuf[adc_samp].end(); iter++) {
				zeroSupHisto[adc_samp]->Fill(iter->first, iter->second);
			}
		}

		TCanvas *c = new TCanvas(
						Form("Canvas_event%d_%s", evdata.GetEvNum(), fName.Data()),
						Form("Canvas_event%d_%s", evdata.GetEvNum(), fName.Data()),
						500, 400);

		c->Divide(1, 4 );
		// plot the graph to the canvas
		for (auto adc_samp = 0; adc_samp < fMaxSamp; adc_samp++) {
			c->cd(1);
			rawHisto[adc_samp]->Draw("HISTO");
			c->cd(2);
			cmSupHisto[adc_samp]->Draw("HISTO");
			c->cd(3);
			pedSupHisto[adc_samp]->Draw("HISTO");
			c->cd(4);
			zeroSupHisto[adc_samp]->Draw("HISTO");
			//c->Draw();
			c->Update();
			c->SaveAs(Form("result/Canvas_event%d_TS%d_%s.jpg", evdata.GetEvNum(),adc_samp, fName.Data()));
		}
		c->Close();
		c->Delete();


		// release the memory
		for (auto adc_samp = 0; adc_samp < fMaxSamp; adc_samp++) {
			rawHisto[adc_samp]->Delete();
			cmSupHisto[adc_samp]->Delete();
			pedSupHisto[adc_samp]->Delete();
			zeroSupHisto[adc_samp]->Delete();
		}

    }
if(evdata.GetEvNum()>100) exit(0);


/*    TCanvas *cped=new TCanvas(Form("Canvas_ped%d_%s",evdata.GetEvNum(),fName.Data()),Form("Canvas_ped%d_%s",evdata.GetEvNum(),fName.Data()),1000,1000);
    TH1F *pedestalHisto=new TH1F(Form("Test_Pedestal_histo"),Form("Test_Pedestal_histo"),600,0,600);
    TH1F *rmsHisto=new TH1F(Form("Test_Pedestal_rms_histo"),Form("Test_Pedestal_rms_histo"),600,0,600);

    for (auto iter=_pedestalHistoBuf.begin();iter!=_pedestalHistoBuf.end();iter++){
    		pedestalHisto->Fill(iter->first,iter->second);
    	}
    for (auto iter=_rmsHistoBuf.begin();iter!=_rmsHistoBuf.end();iter++){
    		rmsHisto->Fill(iter->first,iter->second);
    	}

    cped->Divide(1,2);
    cped->cd(1);
    pedestalHisto->Draw("HISTO");
    cped->cd(2);
    rmsHisto->Draw("HISTO");
    cped->Draw();
    cped->Update();
    getchar();*/

/*
if (_above_th)
  {
		TCanvas *c = new TCanvas(
				Form("Canvas_event%d_%s", evdata.GetEvNum(), fName.Data()),
				Form("Canvas_event%d_%s", evdata.GetEvNum(), fName.Data()),
				1000, 1000);

		TH1F *rawHisto = new TH1F(Form("%s_raw_histo", fName.Data()),
				Form("%s_raw_histo", fName.Data()), 600, 0, 600);
		TH1F *cmSupHisto = new TH1F(Form("Test_commonSup_histo"),
				Form("Test_commonSup_histo"), 600, 0, 600);
		TH1F *pedSupHisto = new TH1F(Form("Test_pedSup_histo"),
				Form("Test_pedSup_histo"), 600, 0, 600);
		TH1F *zeroSupHisto = new TH1F(Form("Test_zeroSup_histo"),
				Form("Test_zeroSup_histo"), 600, 0, 600);

		for (auto iter = _rawHistoBuf.begin(); iter != _rawHistoBuf.end();
				iter++) {
			//std::cout<<"Pos"<<iter->first<<"  adc:"<<iter->second<<std::endl;
			rawHisto->Fill(iter->first, iter->second);
		}

		for (auto iter = _cmSupHistoBuf.begin(); iter != _cmSupHistoBuf.end();
				iter++) {
			cmSupHisto->Fill(iter->first, iter->second);
		}
		for (auto iter = _pedSupHistoBuf.begin(); iter != _pedSupHistoBuf.end();
				iter++) {
			pedSupHisto->Fill(iter->first, iter->second);
		}
		for (auto iter = _zeroSupHistoBuf.begin();
				iter != _zeroSupHistoBuf.end(); iter++) {
			zeroSupHisto->Fill(iter->first, iter->second);
		}

		c->Divide(1, 4);
		c->cd(1);
		rawHisto->Draw("histo");
		c->cd(2);
		cmSupHisto->Draw("histo");
		c->cd(3);
		pedSupHisto->Draw("histo");
		c->cd(4);
		zeroSupHisto->Draw("histo");
		c->SaveAs(Form("result/Canvas_event%d_%s.jpg", evdata.GetEvNum(), fName.Data()));
		c->Draw();
		c->Update();
		c->Modified();

		getchar();

		rawHisto->Delete();
		cmSupHisto->Delete();
		pedSupHisto->Delete();
		zeroSupHisto->Delete();
//    c->Delete();
	}*/

    return FindGEMHits();
}

Int_t MPDGEMPlane::GetRStripNumber( UInt_t strip, UInt_t pos, UInt_t invert ){
    Int_t RstripNb = APVMAP[strip];
    RstripNb=RstripNb+(127-2*RstripNb)*invert;
    Int_t RstripPos = RstripNb + 128*pos;

    return RstripPos;
}


Int_t MPDGEMPlane::FindGEMHits(){
  // Find and analyze clusters. Clusters of active strips are considered
  // a "Hit".
  //
  // The cluster analysis is a critical part of the GEM analysis. Various
  // things can and probably need to be done right here already: splitting
  // oversized clusters, detecting noise hits/bogus clusters, detecting and
  // fitting overlapping clusters etc.
  //
  // This analysis may even need to be re-done after preliminary tracking to
  // see if the clustering can be improved using candidate tracks.
  // Additionally, correlated amplitude information from a second readout
  // direction in the same readout plane could be used here. These advanced
  // procedures would require significant redesign of the code:
  // all raw strip info will have to be saved and prcessed at a later point,
  // similar to the finding of hit pairs in like-oriented planes of the MWDC.
  //
  // For the moment, we implement a very simple algorithm: any cluster of
  // strips larger than what a single cluster should be is assumed to be two or
  // more overlapping hits, and the cluster will be split as follows: anything
  // that looks like a local peak followed by a valley will be considered an
  // actual cluster. The parameter frac = fSplitFrac (0.0 ... 1.0) can
  // be used for some crude tuning. frac > 0.0 means that a peak is
  // only a peak if the amplitude drops below (1-frac), so
  // frac = 0.1 means: trigger on a drop below 90% etc. Likewise for the
  // following valley: the bottom is found if the amplitude rises again
  // by (1+frac), so frac = 0.1 means: trigger on a rise above 110% etc.

  // The active strip numbers must be sorted for the clustering algorithm


    UInt_t nHits = 0;

    sort( ALL(fSigStrips) );

#ifndef NDEBUG
    TreeSearch::GEMHit* prevHit = 0;
#endif



    Double_t frac_down = 1.0 - fSplitFrac, frac_up = 1.0 + fSplitFrac;

    typedef Vint_t::iterator viter_t;
    Vint_t splits;                        // Strips with ampl split between 2 clusters
    viter_t next = fSigStrips.begin();
    while( next != fSigStrips.end() ) {
        viter_t start = next, cur = next;
        ++next;
        assert( next == fSigStrips.end() or *next > *cur );
        while( next != fSigStrips.end() and (*next - *cur) == 1  ) {
            ++cur;     // current
            ++next;
        }
        // Now the cluster candidate is between start and cur
        assert( *cur >= *start );

        // The "type" parameter indicates the result of the cluster analysis:
        // 0: clean (i.e. smaller than fMaxClusterSize, no further analysis)
        // 1: large, maximum at right edge, not split
        // 2: large, no clear minimum on the right side found, not split
        // 3: split, well-defined peak found (may still be larger than maxsize)
        Int_t  type = 0;
        UInt_t size = *cur - *start + 1;
        if( size > fMaxClusterSize ) {
            Double_t maxadc = 0.0, minadc = kBig;
            viter_t it = start, maxpos = start, minpos = start;
            enum EStep { kFindMax = 1, kFindMin, kDone };
            EStep step = kFindMax;
            while( step != kDone and it != next ) {
                Double_t adc = fADCcor[*it];

                switch( step ) {
                    case kFindMax:
                        // Looking for maximum
                        if( adc > maxadc ) {
                            maxpos = it;
                            maxadc = adc;
                        } else if( adc < maxadc * frac_down ) {
                            assert( maxadc > 0.0 );
                            step = kFindMin;
                            continue;
                        }
                        break;
                    case kFindMin:
                        // Looking for minimum
                        if(  adc < minadc ) {
                            minpos = it;
                            minadc = adc;
                        } else if( adc > minadc * frac_up ) {
                            assert( minadc < kBig );
                            step = kDone;
                        }
                        break;
                    case kDone:
                        assert( false );  // should never get here
                        break;
                }
                ++it;
            }
            if( step == kDone ) {
                // Found maximum followed by minimum
                assert( minpos != start );
                assert( minpos != cur );
                assert( *minpos > *maxpos );
                // Split the cluster at the position of the minimum, assuming that
                // the strip with the minimum amplitude is shared between both clusters
                cur  = minpos;
                next = minpos;
                // In order not to double-count amplitude, we split the signal height
                // of that strip evenly between the two clusters. This is a very
                // crude way of doing what we really should be doing: "fitting" a peak
                // shape and using the area and centroid of the curve
                fADCcor[*minpos] /= 2.0;
                splits.push_back(*minpos);
            }
            type = step;
            size = *cur - *start + 1;
            assert( *cur >= *start );
        }
        assert( size > 0 );

        // Compute weighted position average. Again, a crude (but fast) substitute
        // for fitting the centroid of the peak.
        Double_t xsum = 0.0, adcsum = 0.0;
        for( ; start != next; ++start ) {
            Int_t istrip = *start;
            Double_t pos = GetStart() + istrip * GetPitch();
            Double_t adc = fADCcor[istrip];
            xsum   += pos * adc;
            adcsum += adc;
        }

        assert( adcsum > 0.0 );
        Double_t pos = xsum/adcsum;

        // The resolution (sigma) of the position measurement depends on the
        // cluster size. In particular, if the cluster consists of only a single
        // hit, the resolution is much reduced
        Double_t resolution = fResolution;
        if( size == 1 ) {
            resolution = TMath::Max( 0.25*GetPitch(), fResolution );
            // The factor of 1/2*pitch is just a guess. Since with real GEMs
            // there _should_ always be more than one strip per cluster, we must
            // assume that the other strip(s) did not fire due to inefficiency.
            // As a result, the error is bigger than it would be if only ever one
            // strip fired per hit.
            //       resolution = TMath::Max( 0.5*GetPitch(), 2.0*fResolution );
            //     } else if( size == 2 ) {
            //       // Again, this is a guess, to be quantified with Monte Carlo
            //       resolution = 1.2*fResolution;
    }

    // Make a new hit
#ifndef NDEBUG
    TreeSearch::GEMHit* theHit = 0;
#endif
#ifndef NDEBUG
    theHit =
#endif
        new( (*fHits)[nHits++] ) TreeSearch::GEMHit( pos,
                adcsum,
                size,
                type,
                resolution,
                this
                );
#ifndef NDEBUG
    // Ensure hits are ordered by position (should be guaranteed by std::map)
    assert( (prevHit == 0) or (theHit->Compare(prevHit) > 0) );
    prevHit = theHit;
#endif
    }

    // Undo amplitude splitting, if any, so fADCcor contains correct ADC values
    for( viter_t it = splits.begin(); it != splits.end(); ++it ) {
        fADCcor[*it] *= 2.0;
    }

    // Negative return value indicates potential problem
    if( nHits > fMaxHits )
        nHits = -nHits;

    return nHits;
}



Int_t   MPDGEMPlane::Begin( THaRunBase* run ){
    TreeSearch::GEMPlane::Begin(run);
    return 0;
}

Int_t   MPDGEMPlane::End( THaRunBase* run ){
    TreeSearch::GEMPlane::End(run);
    return 0;
}

MPDStripData_t MPDGEMPlane::ChargeDep( const std::vector<Float_t>& amp ) {
    // Deconvolute signal given by samples in 'amp', return approximate integral.
    // Currently analyzes exactly 3 samples.
    // From Kalyan Allada
    // NIM A326, 112 (1993)

    //FIXME: from database, proper value for Tp
    const Float_t delta_t = 25.0; // time interval between samples (ns)
    const Float_t Tp      = 50.0; // RC filter time constant (ns)

    assert( amp.size() >= 3 );

    Float_t adcraw = delta_t*(amp[0]+amp[1]+amp[2]);

    // Weight factors calculated based on the response of the silicon microstrip
    // detector:
    // v(t) = (delta_t/Tp)*exp(-delta_t/Tp)
    // Need to update this for GEM detector response(?):
    // v(t) = A*(1-exp(-(t-t0)/tau1))*exp(-(t-t0)/tau2)
    // where A is the amplitude, t0 the begin of the rise, tau1 the time
    // parameter for the rising edge and tau2 the for the falling edge.

    Float_t x = delta_t/Tp;

    Float_t w1 = TMath::Exp(x-1)/x;
    Float_t w2 = -2*TMath::Exp(-1)/x;
    Float_t w3 = TMath::Exp(-x-1)/x;

    // Deconvoluted signal samples, assuming measurements of zero before the
    // leading edge
    Float_t sig[3] = { amp[0]*w1,
        amp[1]*w1+amp[0]*w2,
        amp[2]*w1+amp[1]*w2+amp[0]*w3 };

    Float_t adc    = delta_t*(sig[0]+sig[1]+sig[2]);
    Float_t time   = 0;     // TODO

    Bool_t pass;
    // Calculate ratios for 3 samples and check for bad signals
    if( amp[2] > 0 ) {
        Float_t r1 = amp[0]/amp[2];
        Float_t r2 = amp[1]/amp[2];
        pass = (r1 < 1.0 and r2 < 1.0 and r1 < r2);
    } else
        pass = false;

    return MPDStripData_t(adcraw,adc,time,pass);
}

ClassImp(MPDGEMPlane)

