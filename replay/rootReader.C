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
#include <string>
#include <stdlib.h>
#include <cctype>
#pragma cling load("libTreeSearch-GEM.so");
#pragma cling load("../libprexCounting.so");

#define _ROOTREADER_MAX_GEM_CHAMBER 7
#define _ROOTREADER_MAX_TSAMPLE 3

void rootReader(TString fname="test_20532.root"){
	TCanvas *eventCanvas=new TCanvas("CanvasDisplay","CanvasDisplay",1000,1000);
	eventCanvas->Divide(1,2);

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


	// check the detector are listed in the tree
	std::vector<int16_t> chamberList;
	std::map<int16_t,std::vector<int16_t>> TsampleLit;

	std::cout<<"List of Chambers:"<<std::endl;
	for(int16_t chambercount=0; chambercount<=_ROOTREADER_MAX_GEM_CHAMBER;chambercount++){
		if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("RGEM.rgems.x%d.adc1",chambercount))){
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
			if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("RGEM.rgems.x%d.adc%d",chamberID,adc_sample))){
				std::cout<<adc_sample<<"  ";
				//TsampleLit.push_back(adc_sample);
				TsampleLit[chamberID].push_back(adc_sample);
			}
		}
		std::cout<<std::endl;

	}

	// allocate the memory to the tree

	// search for the vdcs
	Int_t fEvtNum=0;
	Int_t fRun=0;
	Int_t fvdcXNum=0;
	double_t fvdcX[100];
	Int_t fvdcYNum=0;
	double_t fvdcY[100];
	double_t fvdc_th[10];
	double_t fvdc_ph[10];

	std::string fEvtNumForm("fEvtHdr.fEvtNum");
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fEvtNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fEvtNumForm.c_str(),&fEvtNum);
	}else{
		std::cout<<"[Warning]:: fEvtNum data did not find in the replay resuly"<<std::endl;
	}
	std::string fRunForm("Event_Branch/fEvtHdr/fEvtHdr.fRun");
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fRunForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fRunForm.c_str(),&fRun);
	}else{
		std::cout<<"[Warning]:: fRun data did not find in the replay resuly"<<std::endl;
	}

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
	std::string fvdc_th_Form(Form("R.tr.th"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdc_th_Form.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdc_th_Form.c_str(),fvdc_th);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::string fvdc_ph_Form(Form("R.tr.ph"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdc_ph_Form.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdc_ph_Form.c_str(),fvdc_ph);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::string fvdcYNumForm(Form("Ndata.R.tr.y"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcYNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcYNumForm.c_str(),&fvdcYNum);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay resuly"<<std::endl;
	}
	std::string fvdcYForm(Form("R.tr.y"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcYForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcYForm.c_str(),fvdcY);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}


    //  chamber / value
	std::map<int16_t,double_t *>fstrip;
	std::map<int16_t,Int_t> fstripNum;
	// chamber / adcID / value
	std::map<Int_t,std::map<Int_t,Int_t>> fadcNum;
	std::map<Int_t,std::map<Int_t,double_t *>> fadc;
	std::string gem_root_header("RGEM.rgems");
	for (auto chamberID : chamberList){
		for(auto adc_sample : TsampleLit[chamberID]){

			std::string fstripNumFormat(Form("Ndata.%s.x%d.strip.number",gem_root_header.c_str(),chamberID));
			if(PRex_GEM_tree->GetListOfBranches()->Contains(fstripNumFormat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fstripNumFormat.c_str(),&fstripNum[chamberID]);
			}
			// load the strip number informations
			fstrip[chamberID]=new double_t [100];//[(fstripNum[chamberID])];
			std::string fstripFormat(Form("%s.x%d.strip.number",gem_root_header.c_str(),chamberID));
			if (PRex_GEM_tree->GetListOfBranches()->Contains(fstripFormat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fstripFormat.c_str(),fstrip[chamberID]);
			}

			std::string fadcNumformat(Form("Ndata.%s.x%d.adc%d",gem_root_header.c_str(),chamberID,adc_sample));
			if(PRex_GEM_tree->GetListOfBranches()->Contains(fadcNumformat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fadcNumformat.c_str(),&fadcNum[chamberID][adc_sample]);
			}

			fadc[chamberID][adc_sample]=new double_t [100];//[(fadcNum[chamberID][adc_sample])];
			std::string fadcformat(Form("%s.x%d.adc%d",gem_root_header.c_str(),chamberID,adc_sample));
			if(PRex_GEM_tree->GetListOfBranches()->Contains(fadcformat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fadcformat.c_str(),fadc[chamberID][adc_sample]);
			}
		}
	}

	// read the data for Y dimension
    //  chamber / value
	std::map<int16_t,double_t *>fstrip_y;
	std::map<int16_t,Int_t> fstripNum_y;
	// chamber / adcID / value
	std::map<Int_t,std::map<Int_t,Int_t>> fadcNum_y;
	std::map<Int_t,std::map<Int_t,double_t *>> fadc_y;
	for (auto chamberID : chamberList){
		for(auto adc_sample : TsampleLit[chamberID]){
			std::string fstripNumFormat(Form("Ndata.%s.y%d.strip.number",gem_root_header.c_str(),chamberID));
			if(PRex_GEM_tree->GetListOfBranches()->Contains(fstripNumFormat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fstripNumFormat.c_str(),&fstripNum_y[chamberID]);
			}

			fstrip_y[chamberID]=new double_t [100];//[(fstripNum[chamberID])];
			std::string fstripFormat(Form("%s.y%d.strip.number",gem_root_header.c_str(),chamberID));
			if (PRex_GEM_tree->GetListOfBranches()->Contains(fstripFormat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fstripFormat.c_str(),fstrip_y[chamberID]);
			}
		}
	}

	// read the vdc values
	for(auto entry=0;entry<PRex_GEM_tree->GetEntries() && entry<2000;entry++){
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

	// load the data to plot
	const double_t positionshift[]={0,256.0,256.0,256.0,768.0,768.0,768.0};   //  the center of the detector
	double_t positionZpos[]={0.0,  0.9, 1.5858, 1.789, 2.34135, 2.44305, 2.54};   // updated version of the z-position. take the vdc to be 0,

	const double_t positionshift_x[]={0,256.0,256.0,256.0,768.0,768.0,768.0};   //  the center of the detector
	const double_t positionshift_y[]={0,128.0,128.0,128.0,640.0,640.0,640.0};

	// write the data for x-z dimension
	// write the data for y-z dimension
	std::map<int16_t,TH1F *> GEMHisto_xz;
	std::map<int16_t,TH1F *> GEMHisto_yz;
	for(auto chamberID : chamberList){

		if(fstripNum[chamberID]>0){ // if there is hit on this dimension
			if(!(GEMHisto_xz.find(chamberID)!=GEMHisto_xz.end())){
			GEMHisto_xz[chamberID]=new TH1F(Form("chamber%d_xz",chamberID),Form("chamber%d_xz",chamberID),6/0.00004,-3,3);
			GEMHisto_xz[chamberID]->SetMarkerStyle(20);
			GEMHisto_xz[chamberID]->SetMarkerSize(1);
			}
			for(auto strips_iter=0;strips_iter<fstripNum[chamberID];strips_iter++){
				double_t x=(double_t)(fstrip[chamberID][strips_iter]-positionshift[chamberID])*0.0004;
				double_t z=positionZpos[chamberID];
				if(chamberID>=4){   // for UVa GEM
					x=-((double_t)(fstrip[chamberID][strips_iter]-positionshift[chamberID])*0.0004);
					z=positionZpos[chamberID];
				}
				GEMHisto_xz[chamberID]->Fill(x,z);
			}
		}

		// Y-Z plane
		if(fstripNum_y[chamberID]>0){
			if(!(GEMHisto_yz.find(chamberID)!=GEMHisto_yz.end())){
				GEMHisto_yz[chamberID]=new TH1F(Form("chamber%d_yz",chamberID),Form("chamber%d_yz",chamberID),6/0.00004,-3,3);
				GEMHisto_yz[chamberID]->SetMarkerStyle(20);
				GEMHisto_yz[chamberID]->SetMarkerSize(1);
			}
			for(auto strips_iter=0;strips_iter<fstripNum_y[chamberID];strips_iter++){
				double_t y=(double_t)(fstrip_y[chamberID][strips_iter]-positionshift_y[chamberID])*0.0004;
				double_t z=positionZpos[chamberID];
				if(chamberID>=4){   // for UVa GEM
					y=-((double_t)(fstrip_y[chamberID][strips_iter]-positionshift_y[chamberID])*0.0004);
					z=positionZpos[chamberID];
				}
				GEMHisto_yz[chamberID]->Fill(y,z);
			}
		}
	}


	TLine *beamcenter=new TLine(0,positionZpos[6],0,0);
	beamcenter->SetLineWidth(1);
	beamcenter->SetLineColor(45);
	eventCanvas->SetTitle(Form("Tracking_Detector_Run%d_Evt%d",fEvtNum,fRun));
	eventCanvas->cd(1);
	TH1F *trackingHut_xz=new TH1F("Tracking X-Z ","Tracking X-Z",2/0.00004,-1.1,1.1);
	trackingHut_xz->GetYaxis()->SetRangeUser(0,2.8);
	trackingHut_xz->Draw("histp");

	// draw beam center
	beamcenter->Draw("same");
	// draw VDC
	TLine *vdcplane_xz=new TLine(-1.059,0,1.059,0);
	vdcplane_xz->SetLineWidth(2);
	vdcplane_xz->Draw("same");

	// draw GEM planestd::map<int16_t,TLine *>
	std::map<int16_t,TLine *>GEMPlane_xz;
	for (int i =1; i <=6;i++){
		GEMPlane_xz[i]=new TLine(-positionshift_x[i]*0.0004,positionZpos[i],positionshift_x[i]*0.0004,positionZpos[i]);
		GEMPlane_xz[i]->SetLineWidth(2);
		GEMPlane_xz[i]->Draw("same");
	}

	for(auto histo=GEMHisto_xz.begin();histo!=GEMHisto_xz.end();histo++){
		histo->second->Draw("HISTPsame");
	}

	// draw VDC
	if(fvdcXNum>0){
		TH1F *vdchisto=new TH1F("vdc_xz","vdc_xz",3.8/0.00004,-2.5,1.3);
		vdchisto->Fill(fvdcX[0],0.00001);     // change to the first track
		vdchisto->SetMarkerStyle(20);
		vdchisto->SetMarkerColor(2);
		vdchisto->SetMarkerSize(1);
		vdchisto->Draw("HISTPsame");

		TLine *vdcTrack=new TLine(fvdcX[0],0.001,fvdcX[0]+positionZpos[6]*fvdc_th[0],positionZpos[6]);
		vdcTrack->SetLineWidth(1);
		vdcTrack->SetLineColor(6);
		vdcTrack->Draw("HISTPsame");
	}


	// draw the y-z plane
	eventCanvas->cd(2);
	TH1F *trackingHut_yz=new TH1F("Tracking Y-Z ","Tracking Y-Z",2/0.00004,-1.0,1.0);
	trackingHut_yz->GetYaxis()->SetRangeUser(0,2.8);
	trackingHut_yz->Draw("histp");

	// draw beam center
	beamcenter->Draw("same");

	// draw VDC
	TLine *vdcplane_yz=new TLine(-0.2,0,0.2,0);
	vdcplane_yz->SetLineWidth(2);
	vdcplane_yz->Draw("same");

	// draw GEM planestd::map<int16_t,TLine *>
	std::map<int16_t,TLine *>GEMPlane_yz;
	for (int i =1; i <=6;i++){
		GEMPlane_yz[i]=new TLine(-positionshift_y[i]*0.0004,positionZpos[i],positionshift_y[i]*0.0004,positionZpos[i]);
		GEMPlane_yz[i]->SetLineWidth(2);
		GEMPlane_yz[i]->Draw("same");
	}

	//draw GEM
	for(auto histo=GEMHisto_yz.begin();histo!=GEMHisto_yz.end();histo++){
		histo->second->Draw("HISTPsame");
	}

	// draw VDC
	if(fvdcYNum>0){
		TH1F *vdchisto=new TH1F("vdc_yz","vdc_yz",3.8/0.00004,-2.5,1.3);
		vdchisto->Fill(fvdcY[0],0.00001);     // change to the first track
		vdchisto->SetMarkerStyle(20);
		vdchisto->SetMarkerColor(2);
		vdchisto->SetMarkerSize(1);
		vdchisto->Draw("HISTPsame");

		TLine *vdcTrack=new TLine(fvdcY[0],0.001,fvdcY[0]+positionZpos[6]*fvdc_ph[0],positionZpos[6]);
		vdcTrack->SetLineWidth(1);
		vdcTrack->SetLineColor(6);
		vdcTrack->Draw("HISTPsame");
	}

	eventCanvas->Update();
	if(fvdcXNum>0 || fvdcYNum>0)
	{
//		getchar();
		eventCanvas->SaveAs(Form("result/PRex_Evt%d.jpg",entry));
	}

	}
}
