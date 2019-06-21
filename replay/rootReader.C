#include <TString.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>

#pragma cling load("libTreeSearch-GEM.so");
#pragma cling load("../libprexCounting.so");

void rootReader(TString fname="test_20532.root"){
	if(fname.IsNull()){
		std::cout<<"Please input the file name"<<std::endl;
	}

	TFile *fileio=TFile::Open(fname.Data());
	TTree *PRex_GEM_tree;
	fileio->GetObject("T",PRex_GEM_tree);
	assert(fileio);
	if(PRex_GEM_tree->IsZombie()){
	}

	// need to
	double_t *fadc0;

	PRex_GEM_tree->SetBranchAddress("prex.gems.y6.adc0",fadc0);
	PRex_GEM_tree->GetEntry(25);

	std::cout<<fadc0[0]<<fadc0[1]<<std::endl;

}
