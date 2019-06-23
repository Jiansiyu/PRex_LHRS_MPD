#include <TString.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMatrix.h>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TLine.h>
#pragma cling load("libTreeSearch-GEM.so");
#pragma cling load("../libprexCounting.so");

#define _ROOTREADER_MAX_GEM_CHAMBER 7
#define _ROOTREADER_MAX_TSAMPLE 3
struct position3d{
	double_t x;
	double_t y;
	double_t z;
public:
	position3d(double_t _x,double_t _y,double_t _z):
		x(_x), y(_y),z(_z){};
};

enum dimension{
	X,
	Y,
	Z
};

struct prex_tracking_event{
public:
	prex_tracking_event(){};
	~prex_tracking_event(){};

	// load the the data
	void addGEMHit(){};

	double_t GetZpos(Int_t chamberID){
		switch(chamberID){
		case 1:
			return 0.0;
		case 2:
			return 0.6858;
		case 3:
			return 0.889;
		case 4:
			return 1.104;

		case 5:
			return 1.180;

		case 6:
			return 1.256;
		}
		return 0;
	}
};


struct gemHit{

	gemHit();
public:
	void addHit(int detID,double_t *fadc,double_t *rstrips){};
};

void rootReader(TString fname="test_20532.root"){
	if(fname.IsNull()){
		std::cout<<"Please input the file name"<<std::endl;
	}

	TFile *fileio=TFile::Open(fname.Data());
	assert(fileio);
	TTree *PRex_GEM_tree;
	fileio->GetObject("T",PRex_GEM_tree);
	if(PRex_GEM_tree->IsZombie()){
		std::cout<<"[Error]: can not find tree in the file !!!"<<std::endl;
	}

	TCanvas *eventCanvas=new TCanvas("CanvasDisplay","CanvasDisplay",1000,1000);

	// check the detector are listed in the tree
	std::vector<int16_t> chamberList;
	std::map<int16_t,std::vector<int16_t>> TsampleLit;

	std::cout<<"List of Chambers:"<<std::endl;
	for(int16_t chambercount=0; chambercount<=_ROOTREADER_MAX_GEM_CHAMBER;chambercount++){
		if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("prex.gems.x%d.adc1",chambercount))){
			std::cout<<"	->"<< chambercount<<std::endl;
			chamberList.push_back(chambercount);
		}
	}

	// load the data


	// initialize the buffers
	for (auto chamberID : chamberList){

		std::cout<<"Reading chamber :"<<chamberID<<"\n	sample:";
		//std::vector<int16_t> TsampleLit;
		for (int adc_sample =0; adc_sample <_ROOTREADER_MAX_TSAMPLE; adc_sample++){
			if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("prex.gems.x%d.adc%d",chamberID,adc_sample))){
				std::cout<<adc_sample<<"  ";
				//TsampleLit.push_back(adc_sample);
				TsampleLit[chamberID].push_back(adc_sample);
			}
		}
		std::cout<<std::endl;

	}

	// allocate the memory to the tree

	// search for the vdcs

	Int_t fvdcXNum=0;
	double_t fvdcX[100];
	Int_t fvdcYNum=0;
	double_t fvdcY[100];

	std::string fvdcXNumForm(Form("Ndata.R.tr.x"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcXNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcXNumForm.c_str(),&fvdcXNum);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay resuly"<<std::endl;
	}

	std::string fvdcXForm(Form("R.tr.x"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcXForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcXForm.c_str(),fvdcX);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}


    //  chamber / value
	std::map<int16_t,double_t *>fstrip;
	std::map<int16_t,Int_t> fstripNum;
	// chamber / adcID / value
	std::map<Int_t,std::map<Int_t,Int_t>> fadcNum;
	std::map<Int_t,std::map<Int_t,double_t *>> fadc;
	for (auto chamberID : chamberList){
		for(auto adc_sample : TsampleLit[chamberID]){

			std::string fstripNumFormat(Form("Ndata.prex.gems.x%d.strip.number",chamberID));
			if(PRex_GEM_tree->GetListOfBranches()->Contains(fstripNumFormat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fstripNumFormat.c_str(),&fstripNum[chamberID]);
			}
			// load the strip number informations
			fstrip[chamberID]=new double_t [100];//[(fstripNum[chamberID])];
			std::string fstripFormat(Form("prex.gems.x%d.strip.number",chamberID));
			if (PRex_GEM_tree->GetListOfBranches()->Contains(fstripFormat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fstripFormat.c_str(),fstrip[chamberID]);
			}

			std::string fadcNumformat(Form("Ndata.prex.gems.x%d.adc%d",chamberID,adc_sample));
			if(PRex_GEM_tree->GetListOfBranches()->Contains(fadcNumformat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fadcNumformat.c_str(),&fadcNum[chamberID][adc_sample]);
			}

			fadc[chamberID][adc_sample]=new double_t [100];//[(fadcNum[chamberID][adc_sample])];
			std::string fadcformat(Form("prex.gems.x%d.adc%d",chamberID,adc_sample));
			if(PRex_GEM_tree->GetListOfBranches()->Contains(fadcformat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fadcformat.c_str(),fadc[chamberID][adc_sample]);
			}
		}
	}

	// read the vdc values

	for(auto entry=0;entry<PRex_GEM_tree->GetEntries();entry++){
	// load the data to the buff
	PRex_GEM_tree->GetEntry(entry);
	for(auto chamberID:chamberList){
		std::cout<<"Chamber: "<<chamberID<<std::endl;

		// print out the fired strips
		std::cout<<"	fired strips("<<fstripNum[chamberID]<<") :: ";
		for(int i =0 ; i < fstripNum[chamberID];i++){
			std::cout<<fstrip[chamberID][i]<<" ";
		}
		std::cout<<std::endl;

		//print out the ADC value
		// print out the fired strips
		std::cout<<"	ADC0 ("<<fadcNum[chamberID][0]<<") 	:: ";
		for(int i =0 ; i < fadcNum[chamberID][0];i++){
			std::cout<<fadc[chamberID][0][i]<<" ";
		}
		std::cout<<std::endl;
		std::cout<<"	ADC1 ("<<fadcNum[chamberID][1]<<") 	:: ";
		for(int i =0 ; i < fadcNum[chamberID][1];i++){
			std::cout<<fadc[chamberID][1][i]<<" ";
		}
		std::cout<<std::endl;
		std::cout<<"	ADC0 ("<<fadcNum[chamberID][2]<<") 	:: ";
		for(int i =0 ; i < fadcNum[chamberID][2];i++){
			std::cout<<fadc[chamberID][2][i]<<" ";
		}
		std::cout<<std::endl;
	}

	// data verification

	// load the data to plot

	const double_t positionshift[]={0,256.0,256.0,256.0,768.0,768.0,768.0};
	double_t positionZpos[]={0.0,  0.9, 1.5858, 1.789, 2.34135, 2.44305, 2.54};   // updated version of the z-position. take the vdc to be 0,

	// chamber / hit histo
	std::map<int16_t,TH1F *> hitHisto;
	for(auto chamberID : chamberList){
		if(fstripNum[chamberID]>0){
			if(!(hitHisto.find(chamberID)!=hitHisto.end())){
				hitHisto[chamberID]= new TH1F(Form("chamber%d",chamberID),Form("chamber%d",chamberID),6/0.00004,-3,3);
//				hitHisto[chamberID]->GetYaxis()->SetRangeUser(0,2.1);
				hitHisto[chamberID]->SetMarkerStyle(20);
				hitHisto[chamberID]->SetMarkerSize(1);
			}

			// loop on strips
			for(auto strips_iter=0;strips_iter<fstripNum[chamberID];strips_iter++){
				double_t x=(double_t)(fstrip[chamberID][strips_iter]-positionshift[chamberID])*0.0004;
				double_t y=positionZpos[chamberID];

				// caution UVa GEM need to match the VDC coordination
				if(chamberID>=4){
					x=-(double_t)(fstrip[chamberID][strips_iter]-positionshift[chamberID])*0.0004;
					y=positionZpos[chamberID];
				}

				std::cout<<"chamber:"<<chamberID<<"  ("<<x<<", "<<y<<")"<<std::endl;

				// apply the rotation for the GEM detectors
//				x=0.7071*(x-y);
//				y=0.7071*(x+y);

				double_t rotation_x=0.7071*(x-y);
				double_t rotation_y=0.7071*(x+y);

//				double_t rotation_x=x;
//				double_t rotation_y=y;

				std::cout<<"chamber:"<<chamberID<<"rotation:  ("<<rotation_x<<", "<<rotation_y<<")"<<std::endl;
				hitHisto[chamberID]->Fill(rotation_x,rotation_y);
			}
		}
	}
	eventCanvas->cd();
	Bool_t flag=kFALSE;


	TH1F *trackingHut=new TH1F("Tracking detectors","Tracking detectors",3.8/0.00004,-2.5,1.3);
	trackingHut->GetYaxis()->SetRangeUser(-1,2.8);
	trackingHut->Draw();
	flag=kTRUE;
	// draw the detector positions
	// SBU gem detectors

	std::map<int16_t,TLine *> detectorline;

	for(auto histo=hitHisto.begin();histo!=hitHisto.end();histo++){
		if(!flag){
			(histo->second)->Draw("P");
			flag=kTRUE;
		}else{
			(histo->second)->Draw("Psame");
		}
	}

	// draw the detector planes
	// beam center reference line


	TLine *beamcenter=new TLine(-2,2,0,0);
	beamcenter->SetLineWidth(1);
	beamcenter->SetLineColor(45);
	beamcenter->Draw("same");


	TLine *vdcplane=new TLine(-1.059,0,1.059,0);
	vdcplane->SetLineWidth(2);
	vdcplane->Draw("same");
	double_t detector_start_pos[]={0,-256*0.0004,-256*0.0004,-256*0.0004,-128*6*0.0004,-128*6*0.0004,-128*6*0.0004};
	for (int i =1; i <=6;i++){
//			TLine *line = new TLine(0,10,300,900);

		double_t start_x=detector_start_pos[i];
		double_t start_y=positionZpos[i];
		double_t end_x=-detector_start_pos[i];
		double_t end_y=positionZpos[i];
		//std::cout<<"before rotation start point:("<<start_x<<", "<<start_y<<")"<<"sum:"<<start_x+start_y<<" end point:("<<end_x<<",  "<<end_y<<")"<<std::endl;
		// apply the rotation matrix

		double_t rotation_start_x=(start_x-start_y)*0.7071;
		double_t rotation_start_y=(start_x+start_y)*0.7071;

		double_t rotation_end_x=(end_x-end_y)*0.7071;
		double_t rotation_end_y=(end_x+end_y)*0.7071;
		//std::cout<<"start point:("<<start_x<<", "<<start_y<<") end point:("<<end_x<<",  "<<end_y<<")"<<std::endl;
		//std::cout<<"start point:("<<rotation_start_x<<", "<<rotation_start_y<<") end point:("<<rotation_end_x<<",  "<<rotation_end_y<<")"<<std::endl;
		detectorline[i]=new TLine(rotation_start_x,rotation_start_y,rotation_end_x,rotation_end_y);
		detectorline[i]->SetLineWidth(2);
		detectorline[i]->Draw("same");
	}

	// plot the result from VDC if there is any
	// check the size of the x-dimension
	if(fvdcXNum>0){
		TH1F *vdchisto=new TH1F("vdc","vdc",0.6/0.00004,-2.1,0);
		for(auto i =0; i<fvdcXNum;i++){
			vdchisto->Fill(fvdcX[i],positionZpos[0]);
		}
		vdchisto->SetMarkerStyle(20);
		vdchisto->SetMarkerSize(1);
		vdchisto->Draw("HISTPsame");
	}

//	eventCanvas->Draw();
	eventCanvas->Update();
	getchar();
	}
}
