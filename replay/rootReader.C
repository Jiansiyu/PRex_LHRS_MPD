#include <TString.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMatrix.h>
#include <vector>
#include <string>
#include <TCanvas.h>

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

	TCanvas *eventCanvas=new TCanvas("CanvasDisplay","CanvasDisplay",600,600);

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
	double_t positionZpos[]={0.0,0.1,0.6858+0.1,0.889+0.6858+0.1,1.104+0.889+0.6858+0.1,1.180+1.104+0.889+0.6858+0.1,1.256+1.180+1.104+0.889+0.6858+0.1};

	// write the result to histo
	// chamber / hit histo
	std::map<int16_t,TH1F *> hitHisto;
	for(auto chamberID : chamberList){
		if(fstripNum[chamberID]!=0){
			if(!(hitHisto.find(chamberID)!=hitHisto.end())){
				hitHisto[chamberID]= new TH1F(Form("chamber%d",chamberID),Form("chamber%d",chamberID),0.6/0.0004,-0.3,0.3);
				hitHisto[chamberID]->GetYaxis()->SetRangeUser(-0.1,6);
				hitHisto[chamberID]->SetMarkerStyle(20);
				hitHisto[chamberID]->SetMarkerSize(1);
			}

			// loop on strips
			for(auto strips_iter=0;strips_iter<fstripNum[chamberID];strips_iter++){
				double_t x=(double_t)(fstrip[chamberID][strips_iter]-positionshift[chamberID])*0.0004;
				double_t y=positionZpos[chamberID];
				std::cout<<"chamber:"<<chamberID<<"  ("<<x<<", "<<y<<")"<<std::endl;
				hitHisto[chamberID]->Fill(x,y);
			}
		}
	}
	eventCanvas->cd();
	Bool_t flag=kFALSE;

	for(auto histo=hitHisto.begin();histo!=hitHisto.end();histo++){
		if(!flag){
			(histo->second)->Draw("HISTP");
			flag=kTRUE;
		}else{
			(histo->second)->Draw("HISTPsame");
		}
	}
	eventCanvas->Draw();
	eventCanvas->Update();
	getchar();
	}
}
