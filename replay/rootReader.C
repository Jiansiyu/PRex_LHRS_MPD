#include <TString.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TMatrix.h>
#include <vector>

#pragma cling load("libTreeSearch-GEM.so");
#pragma cling load("../libprexCounting.so");

#define _ROOTREADER_MAX_GEM_CHAMBER 7

struct position3d{
	double_t x;
	double_t y;
	double_t z;
public:
	position3d(double_t _x,double_t _y,double_t _z):
		x(_x), y(_y),z(_z){};
};
//
//struct prex_gem_infor{
//public:
//	prex_gem_infor(){
//		gemdetectorpos[1]=position3d(10,10,10);
//
//	};
//	~prex_gem_infor();
//	double_t getDetectorPos(int8_t){
//		return 0;
//	}
//
//private:
//	position3d gemdetectorpos;
//};


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
	// load the data
	double_t **fadc;
	double_t **fRstrips;

	// check the detector are listed in the tree
	std::vector<int16_t> chamberList;
	std::cout<<"List of Chambers:"<<std::endl;
	for(int16_t chambercount=0; chambercount<=_ROOTREADER_MAX_GEM_CHAMBER;chambercount++){
		if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("prex.gems.x%d.adc1",chambercount))){
			std::cout<<"	->"<< chambercount<<std::endl;
			chamberList.push_back(chambercount);
		}
	}

	// initialize the buffers
	for (auto chamberID : chamberList){


	}





//
//	// need to
//	double_t *fadc0;
//
//	PRex_GEM_tree->SetBranchAddress("prex.gems.y6.adc0",fadc0);
//	PRex_GEM_tree->GetEntry(25);
//
//
//
//	std::cout<<fadc0[0]<<fadc0[1]<<std::endl;

}
