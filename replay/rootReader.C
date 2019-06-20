#include <TString.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>

#pragma cling load("libTreeSearch-GEM.so");
#pragma cling load("../libprexCounting.so");

int rootReader(TString fname){
	if(fname.IsNull()){
		std::cout<<"Please inpu the file name"<<std::endl;
		return 0;
	}

	TFile *fileio=TFile::Open(fname.Data());
	TTree *PRex_GEM_tree;
	fileio->GetObject("T",PRex_GEM_tree);

return 1;
}
