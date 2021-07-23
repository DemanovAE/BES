#include <string>
#include <iostream>
#include <vector>
#include <TGraph.h>
#include <Math/SpecFuncMathMore.h>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include "TFile.h"
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>

// Alex class to draw picture
#include "./Analysis_Function.h"

typedef vector<Double_t> VEC;

const Char_t *outPath = "/home/demanov/New_work/MyAnalysis";
const std::map<Int_t, const Char_t *> Energy_scan_text={{7, "7.7"},{11, "11.5"},{14, "14.5"},{19, "19.6"},{27, "27"},{39, "39"},{62, "62.4"}};
TString particle[]={"#pi^{+}","#pi^{-}","K^{+}","K^{-}","p","#bar{p}","#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}","Hadrons","","PID"};
TString bukva[] = {"(a)","(b)","(c)","(d)","(i)", "(f)", "(g)","(h)", "(i)", "(j)","(k)","(l)","(m)","(n)","(o)","(p)", "(q)"};
const Int_t cent[10]={80,70,60,50,40,30,20,10,5,0};

//Double_t pt_bin_ranges[npt_bins]  = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.2,3.6,4.0,4.4,4.8,5.2};
void makeplotstyle();
void DrawArticle(Int_t Energy, VEC pt_publ, VEC v2_publ, VEC v2_publ_err, VEC pt_publ_err, Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString PID, TString charge, Double_t max_y,Double_t max_x);
void DrawSravNiseem(Int_t Energy, VEC pt_publ, VEC v2_publ, VEC v2_publ_err, VEC pt_publ_err, Int_t CentBinMin, Int_t CentBinMax, Int_t harmonic, TString EtaGap, TString PID, TString charge, Double_t max_y,Double_t max_x);
void Plot_pict(Int_t Energy);
int GetNumberParticle(TString PID, TString charge);
void ParAntFlowVsPtForBESarticle( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, Double_t min_x_axis, Double_t max_x_axis,  vector<VEC> pt_publ_par, vector<VEC> v2_publ_par, vector<VEC> pt_publ_err_par, vector<VEC> v2_publ_err_par, vector<VEC> pt_publ_ant, vector<VEC> v2_publ_ant, vector<VEC> pt_publ_err_ant, vector<VEC> v2_publ_err_ant);


TString PhysRef(TString PID);
VEC RebinArticle( Int_t Energy, TString PID, TString charge, Int_t CentBinMin, Int_t CentBinMax);


int Article(Int_t GeV){
	

	makeplotstyle();
	TFile *fstyle = new TFile("style.root");
	TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
	tsty->cd();

	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
      
	VEC Niseem_GeV11_v3_0_40 = { 0.003015, 0.011383, 0.015285, 0.020407, 0.027482, 0.041795, 0.04471 , 0.067144, 0.079029};
	VEC Niseem_GeV11_v3_err_0_40 = { 0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
	VEC Niseem_GeV11_v3_pt_0_40 = { 0.2918, 0.4985, 0.6993, 0.8999, 1.1977, 1.592 , 2.0072, 2.4134, 2.7949};

	VEC Niseem_GeV14_v3_0_40 = { 0.002776, 0.00848 , 0.015829, 0.020706, 0.033391, 0.043138, 0.054956, 0.079994, 0.062053 };
	VEC Niseem_GeV14_v3_err_0_40 = { 0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
	VEC Niseem_GeV14_v3_pt_0_40 = { 0.2814, 0.4782, 0.6949, 0.8916, 1.1939, 1.6069, 2.0072, 2.4157, 2.8316};

	VEC Niseem_GeV19_v3_0_40 = { 0.003962, 0.010819, 0.015229, 0.021675, 0.032571, 0.049141, 0.061224, 0.072081, 0.074368 };
	VEC Niseem_GeV19_v3_err_0_40 = { 0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
	VEC Niseem_GeV19_v3_pt_0_40 = { 0.2922, 0.4935, 0.6883, 0.8961, 1.2013, 1.6039, 2.00  , 2.4026, 2.7987  };

	VEC Niseem_GeV27_v3_0_40 = { 0.004918, 0.011475, 0.018033, 0.023361, 0.035246, 0.05, 0.060246, 0.067623, 0.079918 };
	VEC Niseem_GeV27_v3_err_0_40 = { 0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
	VEC Niseem_GeV27_v3_pt_0_40 = { 0.3117, 0.5065, 0.6948, 0.9091, 1.2143, 1.6039, 2.013 , 2.4221, 2.8247  };

	VEC Niseem_GeV27_v2_0_40 = { 0.021162 , 0.040249 , 0.056432 , 0.070954 , 0.091286 , 0.114523 , 0.13278 , 0.145643 , 0.156846};
	VEC Niseem_GeV27_v2_err_0_40 = {0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
	VEC Niseem_GeV27_v2_pt_0_40 = { 0.2987, 0.5, 0.7013 , 0.9026 , 1.1948 , 1.6039 , 2.0065, 2.3961 , 2.8117 };

	VEC Niseem_GeV39_v3_0_40 = { 0.003584, 0.0105, 0.018235, 0.024737, 0.036138, 0.050377, 0.06379, 0.073919, 0.076215};
	VEC Niseem_GeV39_v3_err_0_40 = { 0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
	VEC Niseem_GeV39_v3_pt_0_40 = { 0.2913, 0.4843, 0.6901 , 0.8897 , 1.1854 , 1.5843 , 1.9899, 2.3769 , 2.7984   };

	VEC Niseem_GeV39_v2_0_40 = { 0.022092, 0.04094, 0.058966, 0.073282, 0.092502, 0.117025, 0.134553, 0.146731, 0.153151};
	VEC Niseem_GeV39_v2_err_0_40 = {0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002};
	VEC Niseem_GeV39_v2_pt_0_40 = {0.3166, 0.5152, 0.7136, 0.9246, 1.2212, 1.6297, 2.0372, 2.444, 2.8435 };

	VEC ptGeV_err ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	//DrawSravNiseem( 39, Niseem_GeV39_v2_pt_0_40, Niseem_GeV39_v2_0_40, Niseem_GeV39_v2_err_0_40, ptGeV_err, 4, 8, 2, "Eta15","Hadrons", "", 0.33,3.6);
	//DrawSravNiseem( 39, Niseem_GeV39_v3_pt_0_40, Niseem_GeV39_v3_0_40, Niseem_GeV39_v3_err_0_40, ptGeV_err, 4, 8, 3, "Eta15","Hadrons", "", 0.18,3.6);

	//DrawSravNiseem( 27, Niseem_GeV27_v2_pt_0_40, Niseem_GeV27_v2_0_40, Niseem_GeV27_v2_err_0_40, ptGeV_err, 4, 8, 2, "Eta15","Hadrons", "", 0.33,3.6);
	//DrawSravNiseem( 27, Niseem_GeV27_v3_pt_0_40, Niseem_GeV27_v3_0_40, Niseem_GeV27_v3_err_0_40, ptGeV_err, 4, 8, 3, "Eta15","Hadrons", "", 0.18,3.6);

	//DrawSravNiseem( 19, Niseem_GeV19_v3_pt_0_40, Niseem_GeV19_v3_0_40, Niseem_GeV19_v3_err_0_40, ptGeV_err, 4, 8, 3, "Eta15","Hadrons", "", 0.18,3.6);
	//DrawSravNiseem( 14, Niseem_GeV14_v3_pt_0_40, Niseem_GeV14_v3_0_40, Niseem_GeV14_v3_err_0_40, ptGeV_err, 4, 8, 3, "Eta15","Hadrons", "", 0.18,3.6);
	//DrawSravNiseem( 11, Niseem_GeV11_v3_pt_0_40, Niseem_GeV11_v3_0_40, Niseem_GeV11_v3_err_0_40, ptGeV_err, 4, 8, 3, "Eta15","Hadrons", "", 0.18,3.6);


	Plot_pict(GeV);
	
	return 0;

}


void DrawArticle(Int_t Energy, VEC pt_publ, VEC v2_publ, VEC v2_publ_err, VEC pt_publ_err ,Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString PID, TString charge, Double_t max_y,Double_t max_x){

	VEC rebin_vec = RebinArticle( Energy, PID, charge, CentBinMin, CentBinMax);
	TString prefix;
	
	if( strncmp(PID, "Hadrons",7)==0){
		prefix="";
    }else{
    	prefix="TPCandTOF";
    	//prefix="combPID";
    }
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<Draw_Object *> flow;
	std::vector<Draw_Object *> ratio;
		
	TFile *file_read;
	//TFile *file_read2 = new TFile(Form("/home/demanov/New_work/MyAnalysis/macro_draw/povarov/Flow_%iGeV_PID_Qw_Corr.root", Energy),"READ");
	//TFile *file_read2 = new TFile(Form("/home/demanov/New_work/MyAnalysis/macro_draw/povarov/newCorr/Flow_%iGeV_run11_PID_Qw.root", Energy),"READ");
	TFile *file_read2 = new TFile(Form("/home/demanov/New_work/MyAnalysis/macro_draw/povarov/auau_%igev_run11_flow_hadrons_effCorr.root", Energy),"READ");

	TFile *file_out = new TFile("OUT.root","RECREATE");	
	if( strncmp(PID, "Hadrons",7)==0){
		//file_read = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath,Energy),"READ");	
		file_read = new TFile(Form("%s/file_hadrons/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2.root",outPath,Energy),"READ");	
    }else{
		//file_read = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("/home/demanov/New_work/MyAnalysis/macro_draw/povarov/newCorr/Flow_%iGeV_run11_PID_Qw_newCorr.root", Energy),"READ");
		
		//file_read = new TFile(Form("%s/file_pid/Flow_%iGeVRun10_PID_10binYesTofdca11_sys.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file_pid/Flow_%iGeV_CombPIDall.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file_pid/Flow_%iGeV_comb_good.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file_pid/Flow_%iGeV_comb_p95.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/macroPID/fitGausRun18/Flow_%iGeV_Run18_comb.root",outPath, Energy),"READ");
		file_read = new TFile(Form("%s/macro_draw/povarov/auau_%igev_run11_flow_hadrons_correctNormalization.root",outPath, Energy),"READ");
		
		//file_read = new TFile(Form("%s/file_pid/Flow_%iGeVper1_PID_10binYesTofdca22_sys.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Energy),"READ");
    }
	if( strncmp(PID, "Proton",6)==0 ){
		for(Int_t p=0; p<pt_publ.size()-1;p++){
			pt_publ[p]=pt_publ[p+1];
			v2_publ[p]=v2_publ[p+1];
			pt_publ_err[p]=pt_publ_err[p+1];
			v2_publ_err[p]=v2_publ_err[p+1];
		}
		pt_publ.pop_back();
		pt_publ_err.pop_back();
		v2_publ.pop_back();
		v2_publ_err.pop_back();

	}

    //Published line
    flow.push_back(new Draw_Object);
    flow[0] -> SetParametrsGraph("graph","#font[42]{Published}",2, 24, 2.0);
    flow[0] -> SetVectorsGraph(pt_publ,v2_publ,pt_publ_err,v2_publ_err);

    //Data line
    flow.push_back(new Draw_Object);
    Analysis -> SetParametrs(Energy, file_read,file_out,2,PID,charge,EtaGap,prefix);
    flow[1] -> SetParametrsGraph("data","#font[42]{Old (multiplicity)}", 4, 22, 2.0);
    Analysis -> FlowVsPt_EtaSub(flow[1], CentBinMin, CentBinMax, rebin_vec);
    
    flow.push_back(new Draw_Object);
    Analysis -> SetParametrs(Energy, file_read2,file_out,2,PID,charge,EtaGap,"");
    flow[2] -> SetParametrsGraph("data","#font[42]{New (summ of pt)}", 6, 23, 2.0);
    Analysis -> FlowVsPt_EtaSub(flow[2], CentBinMin, CentBinMax, rebin_vec);
	  /* 
    flow.push_back(new Draw_Object);
    Analysis -> SetParametrs(Energy, file_read,file_out,2,PID,charge,EtaGap,"CombPID3");
    flow[3] -> SetParametrsGraph("data","#font[42]{Data (w=1)}", 1, 28, 2.0);
    Analysis -> FlowVsPt_EtaSub(flow[3], CentBinMin, CentBinMax, rebin_vec);
    */

	//flow.push_back(new Draw_Object);
    //Analysis -> SetParametrs(Energy, file_read2,file_out,2,PID,charge,EtaGap,prefix);
    //flow[2] -> SetParametrsGraph("data","#font[42]{Data, w=#frac{1}{f_{Corr}}}", 1, 28, 2.0);
    //Analysis -> FlowVsPt_EtaSub(flow[3], CentBinMin, CentBinMax, rebin_vec);

    //ratio line
    ratio.push_back(new Draw_Object);
    Analysis -> RatioGraphEVAL(ratio[0],flow[1],flow[0],"");

	ratio.push_back(new Draw_Object);
    Analysis -> RatioGraphEVAL(ratio[1],flow[2],flow[0],"");

	//ratio.push_back(new Draw_Object);
   	//Analysis -> RatioGraphEVAL(ratio[2],flow[3],flow[0],"");


    //canvas setening
    Draw_Picture_new *canvas = new Draw_Picture_new;

	Double_t min_x=0.;
    
    if( strncmp(PID, "Proton",6)==0){
    	min_x=0.29;
    }
    canvas->SetTLine(min_x, 0., max_x, 0., 1, 2, 2, "Up");
    canvas->SetTLine(min_x, 0., max_x, 0., 1, 2, 2, "Up");
    canvas->SetTLine(min_x, 0.9, max_x, 0.9, 1, 2, 2, "Down");
    canvas->SetTLine(min_x, 1.1, max_x, 1.1, 1, 2, 2, "Down");
    canvas->SetTLine(min_x, 1., max_x, 1., flow[0] -> GetColor(), 2, 1, "Down");
	canvas->SetAxisToCanvsWithBottomPanel(min_x, max_x, -0.01, max_y, 0.83, 1.17);
	canvas->SetDrawObject(flow,"Up");
	canvas->SetDrawObject(ratio,"Down");
	canvas->SetLegend(0.1,0.55,0.4,0.88,1);
	TCanvas *result = canvas->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.1]{ #splitline{Au+Au #sqrt{s_{NN}} = %s GeV Run11}{%s, centrality %i-%i %%}}}",Energy_scan_text.at(Energy), particle[GetNumberParticle(PID,charge)].Data(),cent[CentBinMax+1], cent[CentBinMin]),min_x + 0.36*max_x,0.8*max_y,
													"#font[42]{p_{T} [GeV/c]}","v_{2}","#font[42]{#frac{Data}{Published}}",PhysRef(PID));
	//result->SaveAs(Form("%s/picture/1_comparison_with_article/%iGeV/CombPID/NewPubl_pt_%s%s_%i_%i_%iGeVRun18_%s_new.png",outPath, Energy,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin], Energy, EtaGap.Data()));
	result->SaveAs(Form("%s/macro_draw/povarov/hadron/NewPubl_pt_%s%s_%i_%i_%iGeVRun18_%s_new.png",outPath, PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin], Energy, EtaGap.Data()));
	//result->SaveAs(Form("%s/picture/1_comparison_with_article/%iGeV/weight/3CorrPubl_pt_%s%s_%i_%i_%iGeVRun11_%s_new.png",outPath, Energy,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin], Energy, EtaGap.Data()));
	//result->SaveAs(Form("%s/picture/pdf/1_comparison_with_article/%iGeV/weight/newSysPubl_pt_%s%s_%i_%i_%iGeVRun10_%s.pdf",outPath, Energy,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin], Energy, EtaGap.Data()));
	delete result;
	file_read->Close();
	file_read2->Close();
	file_out->Close();
}

void DrawSravNiseem(Int_t Energy, VEC pt_publ, VEC v2_publ, VEC v2_publ_err, VEC pt_publ_err ,Int_t CentBinMin, Int_t CentBinMax, Int_t harmonic, TString EtaGap, TString PID, TString charge, Double_t max_y,Double_t max_x){

	VEC rebin_vec = RebinArticle( Energy, PID, charge, CentBinMin, CentBinMax);
	TString prefix;
	
	if( strncmp(PID, "Hadrons",7)==0){
		prefix="";
    }else{
    	prefix="TPCandTOF";
    }
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<Draw_Object *> flow;
	std::vector<Draw_Object *> ratio;
		
	TFile *file_read;
	TFile *file_out = new TFile("OUT.root","RECREATE");	
	if( strncmp(PID, "Hadrons",7)==0){
		file_read = new TFile(Form("%s/file_hadrons/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2.root",outPath,Energy),"READ");	
    }else{
		//file_read = new TFile(Form("%s/file/Flow_New%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Energy),"READ");
		file_read = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Energy),"READ");
    }

    //Published line
    flow.push_back(new Draw_Object);
    flow[0] -> SetParametrsGraph("graph","#font[42]{#scale[1.1]{STAR Preliminary}}",4, 24, 2.5);
    flow[0] -> SetVectorsGraph(pt_publ,v2_publ,pt_publ_err,v2_publ_err);

    //Data line
    flow.push_back(new Draw_Object);
    Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,charge,EtaGap,prefix);
    flow[1] -> SetParametrsGraph("data","#font[42]{#scale[1.1]{Data}}", 2, 33, 2.5);
    Analysis -> FlowVsPt_EtaSub(flow[1], CentBinMin, CentBinMax, rebin_vec);

    //ratio line
    ratio.push_back(new Draw_Object);
    Analysis -> RatioGraphEVAL(ratio[0],flow[1],flow[0],"");

    //canvas setening
    Draw_Picture_new *canvas = new Draw_Picture_new;
    
    canvas->SetTLine(0., 0., max_x, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0., max_x, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0.9, max_x, 0.9, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1.1, max_x, 1.1, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1., max_x, 1., flow[0] -> GetColor(), 2, 1, "Down");
	canvas->SetAxisToCanvsWithBottomPanel(0., max_x, -0.01, max_y, 0.83, 1.17);
	canvas->SetDrawObject(flow,"Up");
	canvas->SetDrawObject(ratio,"Down");
	canvas->SetLegend(0.1,0.65,0.4,0.85,1);
	TCanvas *result = canvas->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.1]{ #splitline{Au+Au #sqrt{s_{NN}} = %s GeV Run10}{%s, centrality %i-%i %%}}}",Energy_scan_text.at(Energy), particle[GetNumberParticle(PID,charge)].Data(),cent[CentBinMax+1], cent[CentBinMin]),0.36*max_x,0.8*max_y,
													"#font[42]{p_{T} [GeV/c]}",Form("v_{%i}",harmonic),"#font[42]{#frac{Data}{STAR}}","");
	result->SaveAs(Form("%s/picture/1_comparison_with_article/Niseem_v%i_pt_%s%s_%i_%i_%iGeV_Run10_%s.png",outPath, harmonic, PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin], Energy, EtaGap.Data()));
	delete result;
	file_read->Close();
	file_out->Close();
}


void ParAntFlowVsPtForBESarticle( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, Double_t min_x_axis, Double_t max_x_axis, vector<VEC> pt_publ_par, vector<VEC> v2_publ_par, vector<VEC> pt_publ_err_par, vector<VEC> v2_publ_err_par, vector<VEC> pt_publ_ant, vector<VEC> v2_publ_ant, vector<VEC> pt_publ_err_ant, vector<VEC> v2_publ_err_ant){

	VEC rebin_vec;
	Int_t point = 1;
	TString prefix = "TPCandTOF";

	Int_t Arr_energy[7]={27,11,14,19,27,39,62};
	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={22, 20, 23, 8, 25, 34, 8, 29};
	const Int_t style_data[]={26, 24, 32, 8, 25, 34, 8, 29};
	const Double_t size[]={1.8,1.6,1.6,1.4,1.3,1.3,1.3};
	TString arr_par[]={"Pion","Kaon","Proton"};

	TFile *file_1[7];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasDif = new Draw_Picture_new;
	Draw_Picture_new *canvasDifPubl = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;
	Draw_Picture_new *canvasRatioPubl = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow_ratio;
	std::vector<std::vector<Draw_Object *>> flow_ratio_publ;
	std::vector<std::vector<Draw_Object *>> flow_difference;
	std::vector<std::vector<Draw_Object *>> flow_difference_publ;
	std::vector<Draw_Object *> object;
	
	Double_t min_y_axis_ratio[3]={0.78,0.88,0.8};
	Double_t max_y_axis_ratio[3]={1.06,1.18,2.23};
	Double_t min_y_axis_dif[3]={-0.006, -0.011, -0.007};
	//Double_t min_y_axis_dif[3]={-0.022, -0.022, -0.012};
	Double_t max_y_axis_dif[3]={0.007, 0.006, 0.023};
	//Double_t max_y_axis_dif[3]={0.017, 0.022, 0.044};

	Double_t min_y_axis_ratio_publ[3]={0.88,0.88,0.84};
	Double_t max_y_axis_ratio_publ[3]={1.16,1.12,1.16};
	Double_t min_y_axis_dif_publ[3]={0.05, 0.05, 0.05};
	Double_t max_y_axis_dif_publ[3]={2.14, 2.14, 2.14};

	Double_t coord_y_draw_dif[3];
	Double_t coord_y_draw_ratio[3];

	Draw_Object *parOb = new Draw_Object;
	Draw_Object *antOb = new Draw_Object;

	Int_t b=0;
	Int_t b2=0;


	for(Int_t p=0; p<pt_publ_par[2].size()-1;p++){
		pt_publ_par[2][p]=pt_publ_par[2][p+1];
		v2_publ_par[2][p]=v2_publ_par[2][p+1];
		pt_publ_err_par[2][p]=pt_publ_err_par[2][p+1];
		v2_publ_err_par[2][p]=v2_publ_err_par[2][p+1];
	}
	pt_publ_par[2].pop_back();
	v2_publ_par[2].pop_back();
	for(Int_t p=0; p<pt_publ_ant[2].size()-1;p++){
		pt_publ_ant[2][p]=pt_publ_ant[2][p+1];
		v2_publ_ant[2][p]=v2_publ_ant[2][p+1];
		pt_publ_err_ant[2][p]=pt_publ_err_ant[2][p+1];
		v2_publ_err_ant[2][p]=v2_publ_err_ant[2][p+1];
	}
	pt_publ_ant[2].pop_back();
	v2_publ_ant[2].pop_back();
	

	for(Int_t par=0; par<3; par++){

		coord_y_draw_dif[par]=min_y_axis_dif[par] + (max_y_axis_dif[par] - min_y_axis_dif[par])*0.9;
		coord_y_draw_ratio[par]=min_y_axis_ratio[par] + (max_y_axis_ratio[par] - min_y_axis_ratio[par])*0.95;
		
		Int_t harm=0;

		parOb -> SetParametrsGraph("graph","#font[42]{#scale[1.1]{Par}}",color[par], style[0], size[0]);
    	parOb -> SetVectorsGraph(pt_publ_par[par],v2_publ_par[par],pt_publ_err_par[par],v2_publ_err_par[par]);
    	antOb -> SetParametrsGraph("graph","#font[42]{#scale[1.1]{Ant}}",color[par], style[0], size[0]);
    	antOb -> SetVectorsGraph(pt_publ_ant[par],v2_publ_ant[par],pt_publ_err_ant[par],v2_publ_err_ant[par]);
		
		for(Int_t j=0; j<1; j++){

			rebin_vec = RebinArticle( Arr_energy[j], arr_par[par], "Pos", CentBinMin, CentBinMax);
			//rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"","Eta01",prefix);
			object.push_back(new Draw_Object);
			Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.0]{Data}}"), color[par], style_data[j], size[j]);
		}
		//flow_difference.push_back(object);
		//object.clear();
		
		object.push_back(new Draw_Object);
    	Analysis -> DifferentGraphEVAL(object.back(),parOb,antOb,"STAR Published");
    	flow_difference.push_back(object);
    	object.clear();

    	object.push_back(new Draw_Object);
    	Analysis -> RatioGraphEVAL(object.back(),flow_difference[par][0],flow_difference[par][1],"");
    	flow_difference_publ.push_back(object);
    	object.clear();

		canvasDif->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.95*coord_y_draw_dif[par], 0.08, 0.0);
		canvasDifPubl->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.8*max_y_axis_dif_publ[par], 0.08, 0.0);
		b++;			

			
		for(Int_t j=0; j<1; j++){
			rebin_vec = RebinArticle( Arr_energy[j], arr_par[par], "Pos", CentBinMin, CentBinMax);
			//rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"","Eta01",prefix);
			object.push_back(new Draw_Object);
			Analysis -> RatioParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.0]{Data}}"), color[par], style_data[j], size[j]);
		}
		//flow_ratio.push_back(object);
		//object.clear();
		
		object.push_back(new Draw_Object);
    	Analysis -> RatioGraphEVAL(object.back(),parOb,antOb,"STAR Published");
    	flow_ratio.push_back(object);
    	object.clear();

    	object.push_back(new Draw_Object);
    	Analysis -> RatioGraphEVAL(object.back(),flow_ratio[par][1],flow_ratio[par][0],"");
    	flow_ratio_publ.push_back(object);
    	object.clear();


		canvasRatio->SetNumberPad(Form("#font[42]{%s}", bukva[b2].Data()), "Pad", 0.04, 0.98*coord_y_draw_ratio[par], 0.08, 0.0);
		canvasRatioPubl->SetNumberPad(Form("#font[42]{%s}", bukva[b2].Data()), "Pad", 0.04, 0.96*max_y_axis_ratio_publ[par], 0.08, 0.0);
		b2++;
	
	}

	canvasDif->SetTextInfo(Form("#font[42]{#scale[1.2]{ v_{2}(%s) - v_{2}(%s) }}",particle[GetNumberParticle("Pion","Pos")].Data(),particle[GetNumberParticle("Pion","Neg")].Data() ),"Y",0.4, 0.7,0.4,90);
	canvasDif->SetTextInfo(Form("#font[42]{#scale[1.2]{ v_{2}(%s) - v_{2}(%s) }}",particle[GetNumberParticle("Kaon","Pos")].Data(),particle[GetNumberParticle("Kaon","Neg")].Data() ),"Y",0.4, 0.4,0.4,90);
	canvasDif->SetTextInfo(Form("#font[42]{#scale[1.2]{ v_{2}(%s) - v_{2}(%s) }}",particle[GetNumberParticle("Proton","Pos")].Data(),particle[GetNumberParticle("Proton","Neg")].Data() ),"Y",0.4, 0.07,0.4,90);
	canvasDif->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{Au+Au #sqrt{s_{NN}}=27 GeV, %i-%i%%}}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.3,0.96*coord_y_draw_dif[0],0.08,0.0);
	canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{#splitline{empty: This analysis}{full: STAR Published}}}"),"Info",0.3,0.65*coord_y_draw_dif[0],0.08,0.0);
	//canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.6,0.95*coord_y_draw_dif[0],0.08,0.0);
	canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary }}"),"Prel",0.6,0.95*coord_y_draw_dif[2],0.08,0.0);
	canvasDif->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasDif->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[0], max_y_axis_dif[0], 0.06);
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[1], max_y_axis_dif[1], 0.06);
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[2], max_y_axis_dif[2], 0.06);
	canvasDif->SetDrawObject(flow_difference[0],"OnePad");
	canvasDif->SetLegend(0.07,0.08,0.98,0.2,3);
	TCanvas *result = canvasDif->CanvasNxM(640,1280, 3, 1, 0.29, 0.35, flow_difference,3,0,1);

	result->SaveAs(Form("%s/picture/1_comparison_with_article/good/only27BESupdateDif_vn_%s_cent_%i_%i_2.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/1_comparison_with_article/good/only27BESupdateDif_vn_%s_cent_%i_%i_2.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result;

	canvasRatio->SetTextInfo(Form("#font[42]{#scale[1.2]{#frac{ v_{2}(%s) }{ v_{2}(%s) } }}", particle[GetNumberParticle("Pion","Pos")].Data(), particle[GetNumberParticle("Pion","Neg")].Data()),"Y",0.07,0.8,0.3,0);
	canvasRatio->SetTextInfo(Form("#font[42]{#scale[1.2]{#frac{ v_{2}(%s) }{ v_{2}(%s) } }}", particle[GetNumberParticle("Kaon","Pos")].Data(), particle[GetNumberParticle("Kaon","Neg")].Data()),"Y",0.07,0.5,0.3,0);
	canvasRatio->SetTextInfo(Form("#font[42]{#scale[1.2]{#frac{ v_{2}(%s) }{ v_{2}(%s) } }}", particle[GetNumberParticle("Proton","Pos")].Data(), particle[GetNumberParticle("Proton","Neg")].Data()),"Y",0.07,0.15,0.3,0);
	canvasRatio->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{Au+Au #sqrt{s_{NN}}=27 GeV, %i-%i%%}}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.3,0.82*coord_y_draw_ratio[0],0.08,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{#splitline{empty: This analysis}{full: STAR Published}}}"),"Info",0.3,0.78*coord_y_draw_ratio[0],0.08,0.0);
	//canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.6,0.98*coord_y_draw_ratio[0],0.08,0.0);
	//canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary }}"),"Prel",0.6,0.96*coord_y_draw_ratio[2],0.08,0.0);
	canvasRatio->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasRatio->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[0], max_y_axis_ratio[0], 0.06);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[1], max_y_axis_ratio[1], 0.06);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[2], max_y_axis_ratio[2], 0.06);
	canvasRatio->SetDrawObject(flow_ratio[0],"OnePad");
	canvasRatio->SetLegend(0.04,0.7,0.98,0.85,3);
	TCanvas *result2 = canvasRatio->CanvasNxM(640,1280, 3, 1, 0.29, 0.35, flow_ratio,3,0,0);

	result2->SaveAs(Form("%s/picture/1_comparison_with_article/good/only27BESupdateRatio_vn_%s_cent_%i_%i_2.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result2->SaveAs(Form("%s/picture/1_comparison_with_article/good/only27BESupdateRatio_vn_%s_cent_%i_%i_2.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result2;

	
	canvasRatioPubl->SetTextInfo(Form("#font[42]{#scale[1.2]{#frac{ f(v_{2})[STAR] }{ f(v_{2})[Data] } (%s) }}", particle[GetNumberParticle("Pion","")].Data()),"Y",0.6,0.8,0.3,90);
	canvasRatioPubl->SetTextInfo(Form("#font[42]{#scale[1.2]{#frac{ f(v_{2})[STAR] }{ f(v_{2})[Data] } (%s) }}", particle[GetNumberParticle("Kaon","")].Data()),"Y",0.6,0.5,0.3,90);
	canvasRatioPubl->SetTextInfo(Form("#font[42]{#scale[1.2]{#frac{ f(v_{2})[STAR] }{ f(v_{2})[Data] } (%s) }}", particle[GetNumberParticle("Proton","")].Data()),"Y",0.6,0.15,0.3,90);
	canvasRatioPubl->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasRatioPubl->SetTextInfo(Form("#font[42]{ #scale[1.0]{#splitline{Au+Au #sqrt{s_{NN}}=27 GeV, %i-%i%%}{f(v_{2})=#frac{ v_{2}(X) }{ v_{2}(#bar{X}) }}}}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.3,0.95*max_y_axis_ratio_publ[0],0.08,0.0);
	//canvasRatioPubl->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.6,0.98*coord_y_draw_ratio[0],0.08,0.0);
	//canvasRatioPubl->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary }}"),"Prel",0.6,0.96*coord_y_draw_ratio[2],0.08,0.0);
	canvasRatioPubl->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasRatioPubl->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasRatioPubl->SetAxisToCanvsNM(min_x_axis, max_x_axis, min_y_axis_ratio_publ[0], max_y_axis_ratio_publ[0], 0.06);
	canvasRatioPubl->SetAxisToCanvsNM(min_x_axis, max_x_axis, min_y_axis_ratio_publ[1], max_y_axis_ratio_publ[1], 0.06);
	canvasRatioPubl->SetAxisToCanvsNM(min_x_axis, max_x_axis, min_y_axis_ratio_publ[2], max_y_axis_ratio_publ[2], 0.06);
	canvasRatioPubl->SetDrawObject(flow_ratio_publ[0],"OnePad");
	canvasRatioPubl->SetLegend(0.04,0.7,0.98,0.85,3);
	TCanvas *result3 = canvasRatioPubl->CanvasNxM(640,1280, 3, 1, 0.29, 0.35, flow_ratio_publ,3,0,0);

	result3->SaveAs(Form("%s/picture/1_comparison_with_article/good/PublOnlynly27BESupdateRatio_vn_%s_cent_%i_%i_2.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result3->SaveAs(Form("%s/picture/1_comparison_with_article/good/PublOnlynly27BESupdateRatio_vn_%s_cent_%i_%i_2.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result3;

	canvasDifPubl->SetTextInfo(Form("#font[42]{#scale[1.2]{#frac{ f(v_{2})[Data] }{ f(v_{2})[STAR] } (%s) }}", particle[GetNumberParticle("Pion","")].Data()),"Y",0.6,0.8,0.3,90);
	canvasDifPubl->SetTextInfo(Form("#font[42]{#scale[1.2]{#frac{ f(v_{2})[Data] }{ f(v_{2})[STAR] } (%s) }}", particle[GetNumberParticle("Kaon","")].Data()),"Y",0.6,0.5,0.3,90);
	canvasDifPubl->SetTextInfo(Form("#font[42]{#scale[1.2]{#frac{ f(v_{2})[Data] }{ f(v_{2})[STAR] } (%s) }}", particle[GetNumberParticle("Proton","")].Data()),"Y",0.6,0.15,0.3,90);
	canvasDifPubl->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasDifPubl->SetTextInfo(Form("#font[42]{ #scale[1.0]{#splitline{Au+Au #sqrt{s_{NN}}=27 GeV, %i-%i%%}{f(v_{2})=v_{2}(X) - v_{2}(#bar{X})}}}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.3,0.8*max_y_axis_dif_publ[0],0.08,0.0);
	//canvasDifPubl->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.6,0.98*coord_y_draw_ratio[0],0.08,0.0);
	//canvasDifPubl->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary }}"),"Prel",0.6,0.96*coord_y_draw_ratio[2],0.08,0.0);
	canvasDifPubl->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasDifPubl->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasDifPubl->SetAxisToCanvsNM(min_x_axis, max_x_axis, min_y_axis_dif_publ[0], max_y_axis_dif_publ[0], 0.06);
	canvasDifPubl->SetAxisToCanvsNM(min_x_axis, max_x_axis, min_y_axis_dif_publ[1], max_y_axis_dif_publ[1], 0.06);
	canvasDifPubl->SetAxisToCanvsNM(min_x_axis, max_x_axis, min_y_axis_dif_publ[2], max_y_axis_dif_publ[2], 0.06);
	canvasDifPubl->SetDrawObject(flow_difference_publ[0],"OnePad");
	canvasDifPubl->SetLegend(0.04,0.7,0.98,0.85,3);
	TCanvas *result4 = canvasDifPubl->CanvasNxM(640,1280, 3, 1, 0.29, 0.35, flow_difference_publ,3,0,0);

	result4->SaveAs(Form("%s/picture/1_comparison_with_article/good/PublOnlynly27BESupdateDif_vn_%s_cent_%i_%i_2.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result4->SaveAs(Form("%s/picture/1_comparison_with_article/good/PublOnlynly27BESupdateDif_vn_%s_cent_%i_%i_2.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result4;
	
	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}

/*
VEC RebinArticle(Int_t Energy, TString PID, TString charge, Int_t CentBinMin, Int_t CentBinMax){

	VEC rebin_vector;

	if(strncmp(PID, "Hadrons",7)==0 ){
		rebin_vector = {0.2, 0.3,0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 10.0};
	}

	if(Energy==62){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2 , 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 3.6, 10.0};
				}
			}	
		}
	}

	if(Energy==39){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2 , 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 3.6, 10.0};
				}
			}	
		}
	}


	if(Energy==27){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2, 4.0, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.2, 2.6, 3.2, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 3.6, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.6, 10.0};
				}
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2 , 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 3.2, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2 , 3.6, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.6, 10.0};
				}
			}	
		}
	}


	if(Energy==19){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.2, 2.6, 3.2, 4.0,10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.4, 2.8, 3.2, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.2, 10.0};
				}
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 3.8,10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 4.0, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 10.0};
				}
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.2, 2.6, 3.2, 3.8,10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.4, 2.8, 3.2, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 3.0, 10.0};
				}
			}
		}
	}

	if(Energy==14){
		if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
			rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.4, 2.8, 5.0};
		}
		if(strncmp(PID, "Proton",6)==0 ){
			rebin_vector = {0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.4, 2.8, 5.0};
		}

	}

	if(Energy==11){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.6, 3.0, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.8, 1.6, 2.4, 10.0};
				}
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 5.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 5.0};
				}
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.2, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 5.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.0, 5.0};
				}
			}
		}
	}

	if(Energy==7){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,2.2, 2.4, 2.8, 5.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.6, 1.0, 1.6, 2.4, 5.0};
				}
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.4, 2.8, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 1.0, 2.0, 10.0};
				}
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.4, 2.8, 5.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.6, 1.4, 1.8, 5.0};
				}
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 5.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 1.0, 1.5, 5.0};
				}
			}
		}
	}
	///////
	if(Energy==27){
		rebin_vector = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.2,2.6,3.0,5.0};
		//rebin_vector = {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.6,3.0,5.0};
	}
	if(Energy==39){
		//rebin_vector = {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.6,3.0,5.0};
		rebin_vector = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.2,2.6,3.0,5.0};
	}
	if(Energy<25){
		rebin_vector = {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.6,3.0,5.0};
	}
	/////

	return rebin_vector;
}
*/

VEC RebinArticle(Int_t Energy, TString PID, TString charge, Int_t CentBinMin, Int_t CentBinMax){

	VEC rebin_vector;

	if(strncmp(PID, "Hadrons",7)==0 ){
		//rebin_vector = {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 10.0};
		rebin_vector = {0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 10.0};
	}

	if(Energy==62){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2 , 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 3.6, 10.0};
				}
			}	
		}
	}

	if(Energy==39){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2 , 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 3.6, 10.0};
				}
			}	
		}
	}


	if(Energy==27){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = { 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2, 4.0, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.2, 2.6, 3.2, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = { 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 3.6, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = { 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.6, 10.0};
				}
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2 , 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = { 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 3.2, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2 , 3.6, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.6, 10.0};
				}
			}	
		}
	}


	if(Energy==19){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.2, 2.6, 3.2, 4.0,10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.4, 2.8, 3.2, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.2, 10.0};
				}
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 3.8,10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 4.0, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 10.0};
				}
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.2, 2.6, 3.2, 3.8,10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.4, 2.8, 3.2, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 3.0, 10.0};
				}
			}
		}
	}

	if(Energy==14){
		if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
			rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.4, 2.8, 5.0};
		}
		if(strncmp(PID, "Proton",6)==0 ){
			rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.4, 2.8, 5.0};
		}

	}

	if(Energy==11){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 4.0, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.6, 3.0, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.5, 0.8, 1.6, 2.4, 10.0};
				}
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 5.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 5.0};
				}
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.2, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 5.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = { 0.6, 0.8, 1.0, 1.4, 2.0, 5.0};
				}
			}
		}
	}

	if(Energy==7){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,2.2, 2.4, 2.8, 5.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.6, 1.0, 1.6, 2.4, 5.0};
				}
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.4, 2.8, 10.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 1.0, 2.0, 10.0};
				}
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 ,2.0, 2.2, 2.4, 2.8, 5.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 0.6, 1.4, 1.8, 5.0};
				}
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				if(strncmp(charge, "Pos",3)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 5.0};
				}
				if(strncmp(charge, "Neg",3)==0){
					rebin_vector = {0.2, 1.0, 1.5, 5.0};
				}
			}
		}
	}
	/*
	if(Energy==27){
		rebin_vector = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.2,2.6,3.0,5.0};
		//rebin_vector = {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.6,3.0,5.0};
	}
	if(Energy==39){
		//rebin_vector = {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.6,3.0,5.0};
		rebin_vector = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.2,2.6,3.0,5.0};
	}
	if(Energy<25){
		rebin_vector = {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.6,3.0,5.0};
	}
	*/

	return rebin_vector;
}



int GetNumberParticle(TString PID, TString charge){
	
	Int_t nPID;
    
    if( strncmp(PID, "_proton",7)==0 || strncmp(PID, "Proton",6)==0){
        if(strncmp(charge, "Pos",3)==0){
          nPID=4;
        }
        else{
          nPID=5;
        }
        if(strncmp(charge, "",1)==0){
          nPID=8;
        }
    }
    if( strncmp(PID, "_kaon",5)==0 || strncmp(PID, "Kaon",4)==0){
        if(strncmp(charge, "Pos",3)==0){
          nPID=2;
        }
        else{
          nPID=3;
        }
        if(strncmp(charge, "",1)==0){
          nPID=7;
        }
    }
    if( strncmp(PID, "_pion",5)==0 || strncmp(PID, "Pion",4)==0){
        if(strncmp(charge, "Pos",3)==0){
          nPID=0;
        }
        else{
          nPID=1;
        }
        if(strncmp(charge, "",1)==0){
          nPID=6;
        }
    }
    if( strncmp(PID, "Hadrons",7)==0 ){
        nPID=9;
    }
    if( strncmp(PID, "",1)==0){
        nPID=10;
    }
    if( strncmp(PID, "PID",3)==0){
        nPID=11;
    }
    return nPID;
}

TString PhysRef(TString PID){
	
	if(strncmp(PID, "Hadrons",7)==0){
		return "#font[42]{ *DOI: 10.1103/PhysRevC.86.054908 }";
	}
	if(strncmp(PID, "Pion",4)==0 || strncmp(PID, "Kaon",4)==0 || strncmp(PID, "Proton",6)==0){
		return "#font[42]{ *DOI: 10.1103/PhysRevC.93.014907 }";
	}
	return "";

}

void Plot_pict(Int_t Energy){

	VEC pt_bin_err={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	VEC ptGeV_err ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	VEC ptHadrons  = {0.3, 0.5, 0.7, 0.9, 1.1, 1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0};
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////39 GeV /////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

		// ---------------- Particle species: Proton ----------------
	VEC GeV39_pt_bin_center_Pr_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV39_v2_values_Pr_0_80 = {0.00683153,0.0158907,0.0315793,0.0516096,0.0727302,0.0925249,0.110851,0.126435,0.140123,0.150104,0.160199,0.165746,0.172611,0.179436,0.184446};
	VEC GeV39_v2_stat_error_Pr_0_80 = {0.000520759,0.000195399,0.000163827,0.000164381,0.00018273,0.000216225,0.000267061,0.000342627,0.000453176,0.000557974,0.000736763,0.00100473,0.00138938,0.00155757,0.00303019};  
	// ---------------- Particle species: antiProton --------------
	VEC GeV39_pt_bin_center_Ant_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV39_v2_values_Ant_0_80 = {0.00494108,0.0114479,0.0261197,0.0449316,0.0663671,0.0860504,0.103827,0.120229,0.132341,0.142772,0.151479,0.165265,0.16601,0.17546,0.178903};
	VEC GeV39_v2_stat_error_Ant_0_80 = {0.000897342,0.000371594,0.000306912,0.000307998,0.000345265,0.000415456,0.000522919,0.000689985,0.000821512,0.00110148,0.00148064,0.00206483,0.00291789,0.00336361,0.00684634};
	// ---------------- Particle species: PiM ----------------
	VEC GeV39_pt_bin_center_PiM_0_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV39_v2_values_PiM_0_80 = {0.0223896,0.042742,0.0624241,0.0798173,0.0944707,0.106976,0.11663,0.125184,0.131548,0.135949,0.138445,0.14286,0.148789,0.143203,0.151072};
	VEC GeV39_v2_stat_error_PiM_0_80 = {4.02884e-05,4.85174e-05,6.65846e-05,9.43872e-05,0.000134239,0.000189959,0.000271414,0.000385102,0.00053152,0.000717481,0.000968061,0.00129079,0.00172756,0.00185626,0.00332027};
	// ---------------- Particle species: PiP ----------------
	VEC GeV39_pt_bin_center_PiP_0_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV39_v2_values_PiP_0_80 = {0.0208012,0.041395,0.0611387,0.0783936,0.0932306,0.105414,0.115887,0.123129,0.129545,0.134079,0.136608,0.141293,0.145453,0.144127,0.147837};
	VEC GeV39_v2_stat_error_PiP_0_80 = {4.06899e-05,4.87531e-05,6.64991e-05,9.38308e-05,0.000132824,0.000187112,0.000265895,0.000375892,0.000515683,0.000692813,0.000918648,0.00122249,0.00162044,0.00172633,0.00303448};
	// ---------------- Particle species: KM ----------------
	VEC GeV39_pt_bin_center_KM_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.685,2.955,3.345};
	VEC GeV39_v2_values_KM_0_80 = {0.0129623,0.02816,0.0468495,0.0647472,0.0797767,0.0927023,0.102755,0.109907,0.117535,0.120576,0.121677,0.124428,0.126282,0.144784,0.139207};
	VEC GeV39_v2_stat_error_KM_0_80 = {0.000275854,0.000175113,0.00017688,0.000208595,0.000262253,0.000341829,0.000461195,0.000635014,0.000904939,0.00131148,0.00196627,0.00289826,0.00425974,0.0049056,0.0110844};
	// ---------------- Particle species: KP ----------------
	VEC GeV39_pt_bin_center_KP_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.685,2.955,3.345};
	VEC GeV39_v2_values_KP_0_80 = {0.0131131,0.0280332,0.0473805,0.0650401,0.0797504,0.0919515,0.101452,0.108697,0.114705,0.120665,0.121731,0.128768,0.123222,0.127912,0.134311};
	VEC GeV39_v2_stat_error_KP_0_80 = {0.000244437,0.000153924,0.000154119,0.000180181,0.000224963,0.000290565,0.000386755,0.000523833,0.000719747,0.0010333,0.00151358,0.002247,0.00334946,0.00396006,0.0063682};


	// ---------------- Particle species: PiM ----------------
	VEC GeV39_pt_bin_center_PiM_0_10 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.805,3.465};
	VEC GeV39_v2_values_PiM_0_10 = {0.0147436,0.0249519,0.0351931,0.0444924,0.0515255,0.0572386,0.0624405,0.0665682,0.0680159,0.0696394,0.0726135,0.0706268,0.0619339};
	VEC GeV39_v2_stat_error_PiM_0_10 = {7.87931e-05,9.36376e-05,0.000127759,0.000180173,0.000255988,0.000362511,0.000519879,0.000740854,0.00102086,0.00138664,0.00150644,0.00248854,0.00578472};
	// ---------------- Particle species: PiP ----------------
	VEC GeV39_pt_bin_center_PiP_0_10 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.805,3.465};
	VEC GeV39_v2_values_PiP_0_10 = {0.0118138,0.0229824,0.0337497,0.0432688,0.0511433,0.0577776,0.0629835,0.0666356,0.0688837,0.0707009,0.0758873,0.0793376,0.0710317};
	VEC GeV39_v2_stat_error_PiP_0_10 = {7.97649e-05,9.40754e-05,0.000127369,0.000178608,0.000252476,0.000355599,0.00050783,0.000718169,0.000987558,0.00132172,0.00141911,0.0022635,0.00528366};;
	// ---------------- Particle species: KM ----------------
	VEC GeV39_pt_bin_center_KM_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.435};
	VEC GeV39_v2_values_KM_0_10 = {0.00757911,0.0157954,0.026121,0.0363908,0.0452336,0.0522223,0.058466,0.0605987,0.0681485,0.0680222,0.0658601,0.0699299,0.10366};
	VEC GeV39_v2_stat_error_KM_0_10 = {0.000547625,0.000341637,0.000339806,0.000396405,0.000494154,0.000642514,0.000865204,0.00118888,0.00169966,0.00245065,0.00300416,0.00624675,0.0188689};
	// ---------------- Particle species: KP ----------------
	VEC GeV39_pt_bin_center_KP_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.435};
	VEC GeV39_v2_values_KP_0_10 = {0.00555827,0.0135945,0.0257326,0.0364142,0.045515,0.0526201,0.0574833,0.0609207,0.0644966,0.0688554,0.0702958,0.0619705,0.0345815};
	VEC GeV39_v2_stat_error_KP_0_10 = {0.000483909,0.000298923,0.00029486,0.000341234,0.000422957,0.000545761,0.00072527,0.000981056,0.00134703,0.00191111,0.00236173,0.00523082,0.0185284};
	// ---------------- Particle species: Proton ----------------
	VEC GeV39_pt_bin_center_Pr_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV39_v2_values_Pr_0_10 = {0.000344971,0.0042966,0.0113989,0.0231444,0.0362371,0.0487195,0.0597448,0.0688017,0.0768114,0.0834786,0.0895703,0.0933591,0.0928367,0.0985281,0.0996354};
	VEC GeV39_v2_stat_error_Pr_0_10 = {0.00107229,0.000396078,0.000324797,0.000318868,0.000347991,0.00040676,0.00049864,0.000637553,0.000842688,0.00103321,0.00136462,0.00186161,0.0025729,0.00289198,0.00563872};
	// ---------------- Particle species: antiProton ----------------
	VEC GeV39_pt_bin_center_Ant_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV39_v2_values_Ant_0_10 = {0.00104822,0.00044334,0.00784047,0.0174989,0.0296773,0.042122,0.0521295,0.0608011,0.0658994,0.0712058,0.0747417,0.0849594,0.0807271,0.0906798,0.0848712};
	VEC GeV39_v2_stat_error_Ant_0_10 = {0.00192614,0.000785758,0.000629482,0.000613868,0.000671617,0.000795191,0.000988842,0.00129467,0.00151703,0.00203264,0.00271891,0.00378032,0.00532838,0.00609857,0.0123108};


	// ---------------- Particle species: Proton ----------------
	VEC GeV39_pt_bin_center_Pr_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV39_v2_values_Pr_10_40 = {0.00625568,0.0163721,0.0351739,0.0584888,0.0837165,0.107311,0.129351,0.148297,0.164633,0.176239,0.18757,0.195966,0.205624,0.212297,0.216848};
	VEC GeV39_v2_stat_error_Pr_10_40 = {0.000544837,0.000204691,0.000172133,0.000173382,0.000193466,0.000229269,0.000283342,0.000363145,0.000479106,0.000593701,0.000783826,0.00106789,0.00147689,0.00165122,0.00321044};
	// ---------------- Particle species: antiProton ----------------
	VEC GeV39_pt_bin_center_Ant_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV39_v2_values_Ant_10_40 = {0.00367436,0.0106725,0.0266682,0.0486995,0.0735285,0.097254,0.118446,0.1377,0.154429,0.16567,0.177735,0.194259,0.197949,0.205333,0.221208};
	VEC GeV39_v2_stat_error_Ant_10_40 = {0.000925613,0.00038232,0.000316434,0.000318902,0.000358874,0.000432563,0.000545001,0.000718362,0.000865585,0.00115941,0.00155836,0.00217497,0.00308067,0.00356517,0.00727714};
	// ---------------- Particle species: PiM ----------------
	VEC GeV39_pt_bin_center_PiM_10_40 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.835,3.465};
	VEC GeV39_v2_values_PiM_10_40 = {0.024917,0.0485882,0.0711891,0.0911594,0.108139,0.122315,0.133424,0.142324,0.149065,0.15394,0.156085,0.159299,0.163019,0.16689};
	VEC GeV39_v2_stat_error_PiM_10_40 = {4.20321e-05,5.07775e-05,6.97073e-05,9.88055e-05,0.000140512,0.000198458,0.000282845,0.000400946,0.000554477,0.000747931,0.00100679,0.0013418,0.00132684,0.00296483};
	// ---------------- Particle species: PiP ----------------
	VEC GeV39_pt_bin_center_PiP_10_40 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.835,3.465};
	VEC GeV39_v2_values_PiP_10_40 = {0.0236042,0.0474268,0.0701121,0.0900108,0.106989,0.120903,0.132197,0.140385,0.146647,0.152065,0.151277,0.157466,0.158328,0.155246};
	VEC GeV39_v2_stat_error_PiP_10_40 = {4.24396e-05,5.10583e-05,6.97025e-05,9.83866e-05,0.000139269,0.000195935,0.000277684,0.000392379,0.000541797,0.00072417,0.000959308,0.00127277,0.00122297,0.00272259};
	// ---------------- Particle species: KM ---------------
	VEC GeV39_pt_bin_center_KM_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.835,3.465};
	VEC GeV39_v2_values_KM_10_40 = {0.0132178,0.0312904,0.0537999,0.0745947,0.0931274,0.108946,0.121914,0.131265,0.138464,0.142195,0.146256,0.150271,0.154251,0.144235};
	VEC GeV39_v2_stat_error_KM_10_40 = {0.000290698,0.000185182,0.000187752,0.000222014,0.000279484,0.000362775,0.000485722,0.000669697,0.000946231,0.00137214,0.00202761,0.00303479,0.00317906,0.00722992};
	// ---------------- Particle species: KP ----------------
	VEC GeV39_pt_bin_center_KP_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.805,3.435};
	VEC GeV39_v2_values_KP_10_40 = {0.0143249,0.0324692,0.054843,0.075677,0.0935934,0.108188,0.120058,0.129584,0.137597,0.142903,0.145207,0.147236,0.143582,0.124331};
	VEC GeV39_v2_stat_error_KP_10_40 = {0.000257904,0.00016312,0.000164087,0.000192338,0.000240327,0.000309352,0.00040976,0.000551598,0.00076111,0.00108504,0.00158232,0.00233087,0.00256974,0.00503966};


	// ---------------- Particle species: Proton ----------------
	VEC GeV39_pt_bin_center_Pr_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV39_v2_values_Pr_40_80 = {0.0173274,0.0329489,0.0579302,0.0899572,0.122273,0.153005,0.181775,0.205587,0.22685,0.241718,0.259529,0.254186,0.274926,0.281757,0.299217};
	VEC GeV39_v2_stat_error_Pr_40_80 = {0.00163998,0.000645123,0.000578342,0.000621392,0.000734358,0.000911515,0.00116596,0.00153552,0.00207608,0.00253886,0.00335845,0.00459636,0.00633159,0.00709883,0.013688};
	// ---------------- Particle species: antiProton ----------------
	VEC GeV39_pt_bin_center_Ant_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805};
	VEC GeV39_v2_values_Ant_40_80 = {0.0110292,0.0256727,0.0503999,0.0802053,0.115991,0.143638,0.174016,0.205972,0.22776,0.253882,0.27293,0.298666};
	VEC GeV39_v2_stat_error_Ant_40_80 = {0.00254538,0.00110941,0.000989204,0.00107262,0.00129246,0.00164622,0.00217035,0.00298469,0.00362504,0.00493413,0.00557032,0.0103886};
	// ---------------- Particle species: PiM ----------------
	VEC GeV39_pt_bin_center_PiM_40_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.835};
	VEC GeV39_v2_values_PiM_40_80 = {0.0289077,0.0593903,0.0895478,0.116196,0.140427,0.159371,0.172953,0.185262,0.197665,0.205783,0.21454,0.225593};
	VEC GeV39_v2_stat_error_PiM_40_80 = {0.000142567,0.000177444,0.0002491,0.000359443,0.000519995,0.000745305,0.00106672,0.00149978,0.00204763,0.0027489,0.00295397,0.00481477};
	// ---------------- Particle species: PiP ----------------
	VEC GeV39_pt_bin_center_PiP_40_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.835};
	VEC GeV39_v2_values_PiP_40_80 = {0.028888,0.0589122,0.0884838,0.114564,0.138449,0.156339,0.173457,0.181754,0.19405,0.19931,0.210363,0.227665};
	VEC GeV39_v2_stat_error_PiP_40_80 = {0.000143401,0.000178332,0.000249537,0.000358916,0.000517534,0.000736543,0.00105173,0.00148376,0.00201896,0.00269214,0.00284612,0.0045173};
	// ---------------- Particle species: KM ----------------
	VEC GeV39_pt_bin_center_KM_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.835};
	VEC GeV39_v2_values_KM_40_80 = {0.0225864,0.0437617,0.0699159,0.0982409,0.119146,0.139389,0.153179,0.167905,0.181927,0.188784,0.178602,0.210977};
	VEC GeV39_v2_stat_error_KM_40_80 = {0.000963893,0.000642692,0.000679603,0.000830741,0.00107209,0.00142409,0.00193787,0.0026898,0.00384267,0.00561045,0.00660844,0.0116723};
	// ---------------- Particle species: KP ----------------
	VEC GeV39_pt_bin_center_KP_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.835};
	VEC GeV39_v2_values_KP_40_80 = {0.0239543,0.044623,0.0718769,0.0978395,0.117474,0.136554,0.15272,0.165516,0.17048,0.181233,0.184436,0.167675};
	VEC GeV39_v2_stat_error_KP_40_80 = {0.0008657,0.000573558,0.000600633,0.000725502,0.000927351,0.00121505,0.00162533,0.00219898,0.00304011,0.00429892,0.00495433,0.00975087};

	VEC GeV39_v2_10_20  = {0.0210,0.0393,0.0558,0.0706,0.0846,0.0975,0.1099,0.1205,0.1302,0.1371,0.1444,0.1487,0.1524,0.1592};
	VEC GeV39_v2_10_20_err = {.0000,.0000,.0001,.0001,.0001,.0001,.0002,.0002,.0003,.0004,.0006,.0008,.0011,.0011};

	VEC GeV39_v2_20_30  = {0.0270,0.0509,0.0725,0.0919,0.1102,0.1270,0.1422,0.1556,0.1674,0.1770,0.1840,0.1913,0.1948,0.2023};
	VEC GeV39_v2_20_30_err = {.0001,.0001,.0001,.0002,.0002,.0003,.0004,.0006,.0008,.0010,.0014,.0020,.0027,.0029};

	VEC GeV39_v2_30_40  = {0.0304,0.0576,0.0824,0.1047,0.1251,0.1436,0.1604,0.1751,0.1868,0.1969,0.2032,0.2099,0.2182, 0.2168};
	VEC GeV39_v2_30_40_err = {.0001,.0001,.0001,.0001,.0002,.0002,.0003,.0004,.0005,.0007,.0009,.0013,.0017,.0018};

	VEC GeV39_v2_50_60  = {0.0298,0.0571,0.0839,0.1074,0.1286,0.1469,0.1639,0.1749,0.1885,0.1994,0.1999,0.2116,0.2119, 0.2163};
	VEC GeV39_v2_50_60_err = {.0001,.0002,.0002,.0003,.0005,.0006,.0008,.0011,.0015,.0020,.0027,.0037,.0050,.0052};


	if(Energy==39){
	
		DrawArticle( Energy, GeV39_pt_bin_center_Ant_0_80, GeV39_v2_values_Ant_0_80, GeV39_v2_stat_error_Ant_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Neg", 0.33,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_Pr_0_80, GeV39_v2_values_Pr_0_80, GeV39_v2_stat_error_Pr_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_KP_0_80, GeV39_v2_values_KP_0_80, GeV39_v2_stat_error_KP_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Pos", 0.28,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_KM_0_80, GeV39_v2_values_KM_0_80, GeV39_v2_stat_error_KM_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Neg", 0.28,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_PiP_0_80, GeV39_v2_values_PiP_0_80, GeV39_v2_stat_error_PiP_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Pos", 0.28,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_PiM_0_80, GeV39_v2_values_PiM_0_80, GeV39_v2_stat_error_PiM_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Neg", 0.28,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_Ant_0_10, GeV39_v2_values_Ant_0_10, GeV39_v2_stat_error_Ant_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Neg", 0.22,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_Pr_0_10, GeV39_v2_values_Pr_0_10, GeV39_v2_stat_error_Pr_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Pos", 0.22,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_KP_0_10, GeV39_v2_values_KP_0_10, GeV39_v2_stat_error_KP_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Pos", 0.16,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_KM_0_10, GeV39_v2_values_KM_0_10, GeV39_v2_stat_error_KM_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Neg", 0.16,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_PiP_0_10, GeV39_v2_values_PiP_0_10, GeV39_v2_stat_error_PiP_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Pos", 0.14,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_PiM_0_10, GeV39_v2_values_PiM_0_10, GeV39_v2_stat_error_PiM_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Neg", 0.14,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_Pr_10_40, GeV39_v2_values_Pr_10_40, GeV39_v2_stat_error_Pr_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_Ant_10_40, GeV39_v2_values_Ant_10_40, GeV39_v2_stat_error_Ant_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Neg", 0.33,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_KP_10_40, GeV39_v2_values_KP_10_40, GeV39_v2_stat_error_KP_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Pos", 0.28,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_KM_10_40, GeV39_v2_values_KM_10_40, GeV39_v2_stat_error_KM_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Neg", 0.28,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_PiP_10_40, GeV39_v2_values_PiP_10_40, GeV39_v2_stat_error_PiP_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Pos", 0.25,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_PiM_10_40, GeV39_v2_values_PiM_10_40, GeV39_v2_stat_error_PiM_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Neg", 0.25,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_Pr_40_80, GeV39_v2_values_Pr_40_80, GeV39_v2_stat_error_Pr_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Pos", 0.42,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_Ant_40_80, GeV39_v2_values_Ant_40_80, GeV39_v2_stat_error_Ant_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Neg", 0.42,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_KP_40_80, GeV39_v2_values_KP_40_80, GeV39_v2_stat_error_KP_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_KM_40_80, GeV39_v2_values_KM_40_80, GeV39_v2_stat_error_KM_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Neg", 0.33,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_PiP_40_80, GeV39_v2_values_PiP_40_80, GeV39_v2_stat_error_PiP_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV39_pt_bin_center_PiM_40_80, GeV39_v2_values_PiM_40_80, GeV39_v2_stat_error_PiM_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Neg", 0.33,3.6);
			
		//DrawArticle( Energy, ptHadrons, GeV39_v2_10_20, GeV39_v2_10_20_err, ptGeV_err, 6, 6, "Eta15","Hadrons", "", 0.33,3.6);
		//DrawArticle( Energy, ptHadrons, GeV39_v2_20_30, GeV39_v2_20_30_err, ptGeV_err, 5, 5, "Eta15","Hadrons", "", 0.33,3.6);
		//DrawArticle( Energy, ptHadrons, GeV39_v2_30_40, GeV39_v2_30_40_err, ptGeV_err, 4, 4, "Eta15","Hadrons", "", 0.33,3.6);
		//DrawArticle( Energy, ptHadrons, GeV39_v2_50_60, GeV39_v2_50_60_err, ptGeV_err, 2, 2, "Eta15","Hadrons", "", 0.33,3.6);


	}

	/*
	if(Energy==39){
		//Draw_particle_antipar( Energy, GeV39_pt_bin_center_Pr_0_80, GeV39_v2_values_Pr_0_80, GeV39_v2_stat_error_Pr_0_80, pt_bin_err,
		//							   GeV39_pt_bin_center_Ant_0_80, GeV39_v2_values_Ant_0_80, GeV39_v2_stat_error_Ant_0_80, pt_bin_err,0, 8, "Eta01","Proton", 0.33,3.6);

		Draw_particle_antipar( Energy, GeV39_pt_bin_center_KP_0_80, GeV39_v2_values_KP_0_80, GeV39_v2_stat_error_KP_0_80, pt_bin_err,
									   GeV39_pt_bin_center_KM_0_80, GeV39_v2_values_KM_0_80, GeV39_v2_stat_error_KM_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", 0.28,3.6);

		Draw_particle_antipar( Energy, GeV39_pt_bin_center_PiP_0_80, GeV39_v2_values_PiP_0_80, GeV39_v2_stat_error_PiP_0_80, pt_bin_err,
									   GeV39_pt_bin_center_PiM_0_80, GeV39_v2_values_PiM_0_80, GeV39_v2_stat_error_PiM_0_80, pt_bin_err, 0, 8, "Eta01","Pion", 0.28,3.6);
		
		//Draw_particle_antipar( Energy, GeV39_pt_bin_center_Pr_0_10, GeV39_v2_values_Pr_0_10, GeV39_v2_stat_error_Pr_0_10, pt_bin_err,
		//							   GeV39_pt_bin_center_Ant_0_10, GeV39_v2_values_Ant_0_10, GeV39_v2_stat_error_Ant_0_10, pt_bin_err, 7, 8, "Eta01","Proton", 0.22,3.6);
		
		Draw_particle_antipar( Energy, GeV39_pt_bin_center_KP_0_10, GeV39_v2_values_KP_0_10, GeV39_v2_stat_error_KP_0_10, pt_bin_err,
		 							   GeV39_pt_bin_center_KM_0_10, GeV39_v2_values_KM_0_10, GeV39_v2_stat_error_KM_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", 0.16,3.6);

		Draw_particle_antipar( Energy, GeV39_pt_bin_center_PiP_0_10, GeV39_v2_values_PiP_0_10, GeV39_v2_stat_error_PiP_0_10, pt_bin_err,
							      	   GeV39_pt_bin_center_PiM_0_10, GeV39_v2_values_PiM_0_10, GeV39_v2_stat_error_PiM_0_10, pt_bin_err, 7, 8, "Eta01","Pion", 0.14,3.6);

		//Draw_particle_antipar( Energy, GeV39_pt_bin_center_Pr_10_40, GeV39_v2_values_Pr_10_40, GeV39_v2_stat_error_Pr_10_40, pt_bin_err,
		//							   GeV39_pt_bin_center_Ant_10_40, GeV39_v2_values_Ant_10_40, GeV39_v2_stat_error_Ant_10_40, pt_bin_err, 4, 6, "Eta01","Proton", 0.33,3.6);

		Draw_particle_antipar( Energy, GeV39_pt_bin_center_KP_10_40, GeV39_v2_values_KP_10_40, GeV39_v2_stat_error_KP_10_40, pt_bin_err, 
									   GeV39_pt_bin_center_KM_10_40, GeV39_v2_values_KM_10_40, GeV39_v2_stat_error_KM_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", 0.28,3.6);

		Draw_particle_antipar( Energy, GeV39_pt_bin_center_PiP_10_40, GeV39_v2_values_PiP_10_40, GeV39_v2_stat_error_PiP_10_40, pt_bin_err,
							     	   GeV39_pt_bin_center_PiM_10_40, GeV39_v2_values_PiM_10_40, GeV39_v2_stat_error_PiM_10_40, pt_bin_err, 4, 6, "Eta01","Pion", 0.25,3.6);

		//Draw_particle_antipar( Energy, GeV39_pt_bin_center_Pr_40_80, GeV39_v2_values_Pr_40_80, GeV39_v2_stat_error_Pr_40_80, pt_bin_err,
		//							   GeV39_pt_bin_center_Ant_40_80, GeV39_v2_values_Ant_40_80, GeV39_v2_stat_error_Ant_40_80, pt_bin_err, 0, 3, "Eta01","Proton", 0.42,3.6);

		Draw_particle_antipar( Energy, GeV39_pt_bin_center_KP_40_80, GeV39_v2_values_KP_40_80, GeV39_v2_stat_error_KP_40_80, pt_bin_err,
									   GeV39_pt_bin_center_KM_40_80, GeV39_v2_values_KM_40_80, GeV39_v2_stat_error_KM_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", 0.33,3.6);

		Draw_particle_antipar( Energy, GeV39_pt_bin_center_PiP_40_80, GeV39_v2_values_PiP_40_80, GeV39_v2_stat_error_PiP_40_80, pt_bin_err,
									   GeV39_pt_bin_center_PiM_40_80, GeV39_v2_values_PiM_40_80, GeV39_v2_stat_error_PiM_40_80, pt_bin_err, 0, 3, "Eta01","Pion", 0.33,3.6);
	}
	*/


		// Event plane method: eta-sub
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////27 GeV /////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ---------------- Particle species: Proton ----------------
	VEC GeV27_pt_bin_center_Pr_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.435};
	VEC GeV27_v2_values_Pr_0_80 = {0.00794503,0.016472,0.0328142,0.0523477,0.0729255,0.0922382,0.109006,0.125281,0.137722,0.148279,0.153932,0.163041,0.168032,0.175865,0.180928};
	VEC GeV27_v2_stat_error_Pr_0_80 = {0.000682238,0.000265782,0.000226141,0.000230019,0.000259161,0.000310804,0.000388614,0.000502932,0.00066813,0.000856793,0.00115243,0.00159573,0.00225457,0.00260966,0.00503393};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV27_pt_bin_center_Ant_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.435};
	VEC GeV27_v2_values_Ant_0_80 = {0.000415428,0.00988674,0.0241081,0.0414571,0.0628851,0.0811606,0.0985156,0.11376,0.125846,0.138131,0.145702,0.153858,0.159991,0.162908,0.169893};
	VEC GeV27_v2_stat_error_Ant_0_80 = {0.0015074,0.000638496,0.000538509,0.000549987,0.000621507,0.000753979,0.000951449,0.00124761,0.00153575,0.00206931,0.00280249,0.00391116,0.00555564,0.00648414,0.0125137};

	// ---------------- Particle species: PiM ----------------
	VEC GeV27_pt_bin_center_PiM_0_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.805,3.465};
	VEC GeV27_v2_values_PiM_0_80 = {0.0218633,0.0414492,0.0604151,0.0773144,0.0910853,0.1024,0.11233,0.119006,0.124619,0.132218,0.131911,0.140354,0.138048,0.143187};
	VEC GeV27_v2_stat_error_PiM_0_80 = {6.25539e-05,7.62958e-05,0.00010536,0.000150461,0.000215751,0.000307269,0.000435508,0.000621287,0.000881245,0.0012254,0.00168957,0.00230782,0.00214512,0.00552546};

	// ---------------- Particle species: PiP ----------------
	VEC GeV27_pt_bin_center_PiP_0_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.805,3.465};
	VEC GeV27_v2_values_PiP_0_80 = {0.0204075,0.0404926,0.0592754,0.0759822,0.0897911,0.100898,0.110136,0.117802,0.121748,0.127693,0.133237,0.13354,0.138953,0.139932};
	VEC GeV27_v2_stat_error_PiP_0_80 = {6.33703e-05,7.70403e-05,0.000105998,0.000151183,0.000217215,0.000309484,0.000439386,0.000629185,0.000895574,0.00124855,0.00171713,0.00236994,0.00242954,0.00593679};

	// ---------------- Particle species: KM ----------------
	VEC GeV27_pt_bin_center_KM_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.805,3.435};
	VEC GeV27_v2_values_KM_0_80 = {0.012537,0.026693,0.044647,0.0608321,0.0762363,0.087412,0.0989232,0.105773,0.114272,0.119576,0.123284,0.127067,0.124691,0.125902};
	VEC GeV27_v2_stat_error_KM_0_80 = {0.000436839,0.000280824,0.000287405,0.000344485,0.000440569,0.000578052,0.000775283,0.00107023,0.00151357,0.0021815,0.00317164,0.00465061,0.00525826,0.0145206};

	// ---------------- Particle species: KP ----------------
	VEC GeV27_pt_bin_center_KP_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.805,3.465};
	VEC GeV27_v2_values_KP_0_80 = {0.0139658,0.0284135,0.0466836,0.0633981,0.0765928,0.0883424,0.0980818,0.104583,0.11152,0.117711,0.119196,0.118682,0.124314,0.106932};
	VEC GeV27_v2_stat_error_KP_0_80 = {0.000369429,0.000234884,0.000238474,0.000283337,0.000360682,0.000470827,0.00062698,0.000855642,0.00118904,0.00167376,0.00240471,0.00356718,0.00423466,0.0126745};
	// ---------------- Particle species: Proton ----------------
	VEC GeV27_pt_bin_center_Pr_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.435};
	VEC GeV27_v2_values_Pr_0_10 = {0.00563779,0.00550533,0.0138872,0.0260189,0.0368498,0.0493384,0.058616,0.0670935,0.0766293,0.0810413,0.0835657,0.0862228,0.0873898,0.0903064,0.088792};
	VEC GeV27_v2_stat_error_Pr_0_10 = {0.00139069,0.000533488,0.000444706,0.000443043,0.000490364,0.00058069,0.000720115,0.000928111,0.00123031,0.00157893,0.00212597,0.00295178,0.00418774,0.00485497,0.009385};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV27_pt_bin_center_Ant_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.535,3.075};
	VEC GeV27_v2_values_Ant_0_10 = {-0.00517344,-0.00238412,0.00694304,0.0132629,0.0269888,0.0419234,0.0476994,0.0579478,0.0637756,0.0724395,0.0718708,0.0828054,0.0791406};
	VEC GeV27_v2_stat_error_Ant_0_10 = {0.00323958,0.00134574,0.00110358,0.00109676,0.00120568,0.00143464,0.00178431,0.00231266,0.00281885,0.00378589,0.00510906,0.00590989,0.0117722};

	// ---------------- Particle species: PiM ----------------
	VEC GeV27_pt_bin_center_PiM_0_10 = {0.285,0.495,0.675,0.885,1.065,1.275,1.545,1.935,2.355,2.835};
	VEC GeV27_v2_values_PiM_0_10 = {0.0138583,0.0237411,0.034084,0.0428316,0.0497145,0.0554956,0.0608235,0.0645608,0.0704744,0.0609363};
	VEC GeV27_v2_stat_error_PiM_0_10 = {0.000121131,0.00014642,0.000201454,0.000287281,0.000410347,0.000585106,0.000679412,0.0013799,0.00265953,0.00457086};

	// ---------------- Particle species: PiP ----------------
	VEC GeV27_pt_bin_center_PiP_0_10 = {0.285,0.495,0.675,0.885,1.065,1.275,1.545,1.935,2.355,2.835};
	VEC GeV27_v2_values_PiP_0_10 = {0.0117471,0.0228517,0.0329,0.042002,0.0487612,0.0553591,0.0595984,0.0637105,0.0678853,0.066606};
	VEC GeV27_v2_stat_error_PiP_0_10 = {0.00012301,0.000147786,0.00020224,0.000287824,0.00041135,0.000586629,0.000681158,0.00139804,0.00267615,0.00466203};

	// ---------------- Particle species: KM ----------------
	VEC GeV27_pt_bin_center_KM_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.575,1.965,2.355,2.805};
	VEC GeV27_v2_values_KM_0_10 = {0.00542912,0.0142381,0.0251752,0.0335317,0.0441924,0.0490159,0.0559731,0.0661542,0.0619125,0.0664626};
	VEC GeV27_v2_stat_error_KM_0_10 = {0.000861928,0.000545604,0.000550471,0.000653157,0.000829444,0.00108271,0.00117228,0.00219677,0.00486547,0.00975503};

	// ---------------- Particle species: KP ----------------
	VEC GeV27_pt_bin_center_KP_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.575,1.965,2.355,2.805};
	VEC GeV27_v2_values_KP_0_10 = {0.00716901,0.0152748,0.0258303,0.0364002,0.0430943,0.0504399,0.0553866,0.0645685,0.0756379,0.0636638};
	VEC GeV27_v2_stat_error_KP_0_10 = {0.000726268,0.000454171,0.000454641,0.000535686,0.000677907,0.000882976,0.000946797,0.00180842,0.00375147,0.00832323};
	// ---------------- Particle species: Proton ----------------
	VEC GeV27_pt_bin_center_Pr_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.435};
	VEC GeV27_v2_values_Pr_10_40 = {0.00724341,0.0175996,0.0364746,0.0595057,0.0845299,0.107969,0.129523,0.148708,0.163259,0.176651,0.1848,0.195875,0.201553,0.209258,0.214492};
	VEC GeV27_v2_stat_error_Pr_10_40 = {0.000703917,0.000274745,0.000234747,0.000239978,0.000271603,0.000326546,0.000408775,0.000528588,0.000701306,0.000901556,0.00121072,0.00167138,0.00235442,0.00272466,0.00525693};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV27_pt_bin_center_Ant_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.435};
	VEC GeV27_v2_values_Ant_10_40 = {3.41673e-05,0.0096609,0.0249018,0.0458325,0.0693557,0.0897841,0.112537,0.131383,0.14558,0.159333,0.173659,0.179558,0.186008,0.196745,0.199973};
	VEC GeV27_v2_stat_error_Ant_10_40 = {0.00152348,0.000645128,0.000545297,0.000559015,0.000634942,0.000772649,0.000977222,0.00128147,0.00159098,0.00214354,0.00290067,0.00404856,0.0057128,0.00668967,0.0129286};

	// ---------------- Particle species: PiM ----------------
	VEC GeV27_pt_bin_center_PiM_10_40 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.805,3.465};
	VEC GeV27_v2_values_PiM_10_40 = {0.0250653,0.0475695,0.069042,0.0883132,0.104646,0.117747,0.12904,0.136938,0.143521,0.15096,0.151491,0.156554,0.138774};
	VEC GeV27_v2_stat_error_PiM_10_40 = {6.43381e-05,7.85067e-05,0.000108488,0.000155235,0.000222099,0.000316629,0.000447165,0.000636118,0.000900244,0.00125199,0.00140095,0.00236064,0.00564326};

	// ---------------- Particle species: PiP ----------------
	VEC GeV27_pt_bin_center_PiP_10_40 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.805,3.465};
	VEC GeV27_v2_values_PiP_10_40 = {0.0235111,0.046407,0.0682154,0.0871932,0.103428,0.116498,0.12745,0.134928,0.141224,0.146375,0.151377,0.154863,0.158313};
	VEC GeV27_v2_stat_error_PiP_10_40 = {6.51565e-05,7.9295e-05,0.00010925,0.000156206,0.000223992,0.000319636,0.000452404,0.000645997,0.000917705,0.0012787,0.00142484,0.0024869,0.00608305};

	// ---------------- Particle species: KM ----------------
	VEC GeV27_pt_bin_center_KM_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.465};
	VEC GeV27_v2_values_KM_10_40 = {0.0136614,0.0309157,0.0513795,0.0710171,0.0892261,0.10274,0.116668,0.126808,0.133538,0.140852,0.146276,0.146791,0.184954};
	VEC GeV27_v2_stat_error_KM_10_40 = {0.00045353,0.000292826,0.000300815,0.00036131,0.000462572,0.000607057,0.000813044,0.00112272,0.00158201,0.00226817,0.00267507,0.00549795,0.0139823};

	// ---------------- Particle species: KP ----------------
	VEC GeV27_pt_bin_center_KP_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.465};
	VEC GeV27_v2_values_KP_10_40 = {0.0147664,0.0324572,0.0541136,0.0742409,0.0902704,0.104575,0.115792,0.124146,0.130864,0.138143,0.136981,0.142757,0.0985712};
	VEC GeV27_v2_stat_error_KP_10_40 = {0.000384169,0.000245618,0.000250251,0.000297922,0.000379676,0.000495558,0.000659117,0.000898014,0.0012436,0.00174332,0.00205821,0.00428589,0.0118752};
	// ---------------- Particle species: Proton ----------------
	VEC GeV27_pt_bin_center_Pr_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.805,3.435};
	VEC GeV27_v2_values_Pr_40_80 = {0.0134984,0.0316417,0.058468,0.0875636,0.122787,0.151113,0.175259,0.20741,0.219541,0.239953,0.240922,0.261432,0.280528,0.319028};
	VEC GeV27_v2_stat_error_Pr_40_80 = {0.00218832,0.000894602,0.000812187,0.000881703,0.00105444,0.00132476,0.0017155,0.00228284,0.00309833,0.00395661,0.00536363,0.00744742,0.00813994,0.0234046};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV27_pt_bin_center_Ant_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.535,3.045};
	VEC GeV27_v2_values_Ant_40_80 = {0.00662859,0.0237488,0.0454148,0.0751553,0.112109,0.137055,0.169079,0.187734,0.218416,0.239201,0.242115,0.26207,0.287476};
	VEC GeV27_v2_stat_error_Ant_40_80 = {0.00422311,0.00187958,0.00170792,0.00188585,0.00229729,0.00296268,0.0039251,0.00542688,0.00680821,0.00937933,0.0117719,0.0152434,0.030083};

	// ---------------- Particle species: PiM ----------------
	VEC GeV27_pt_bin_center_PiM_40_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.965,2.355,2.835};
	VEC GeV27_v2_values_PiM_40_80 = {0.0280839,0.057612,0.0868324,0.113901,0.13458,0.151988,0.167185,0.183013,0.194763,0.208071,0.228593};
	VEC GeV27_v2_stat_error_PiM_40_80 = {0.00022182,0.000279632,0.000394499,0.000574355,0.000836336,0.00120471,0.00171964,0.00246,0.00283214,0.00536875,0.00908207};

	// ---------------- Particle species: PiP ----------------
	VEC GeV27_pt_bin_center_PiP_40_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.965,2.355,2.835};
	VEC GeV27_v2_values_PiP_40_80 = {0.0277103,0.0571959,0.0852322,0.111648,0.13395,0.148716,0.166455,0.181484,0.187652,0.210592,0.236303};
	VEC GeV27_v2_stat_error_PiP_40_80 = {0.000223711,0.000282544,0.00039842,0.000580809,0.000849794,0.00122565,0.00175536,0.00251897,0.00293369,0.00553557,0.00955107};

	// ---------------- Particle species: KM ----------------
	VEC GeV27_pt_bin_center_KM_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.965,2.355,2.805};
	VEC GeV27_v2_values_KM_40_80 = {0.0226684,0.0393333,0.0664529,0.092143,0.110992,0.132704,0.149874,0.155613,0.178867,0.219579,0.203398};
	VEC GeV27_v2_stat_error_KM_40_80 = {0.00153391,0.00103893,0.00111166,0.00138229,0.00181587,0.00241437,0.00327718,0.00454959,0.0052254,0.0112432,0.0197207};

	// ---------------- Particle species: KP ----------------
	VEC GeV27_pt_bin_center_KP_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.965,2.355,2.805};
	VEC GeV27_v2_values_KP_40_80 = {0.0250122,0.0441759,0.0711787,0.0930094,0.114819,0.129191,0.15154,0.159513,0.174123,0.17201,0.227863};
	VEC GeV27_v2_stat_error_KP_40_80 = {0.00132176,0.000885837,0.000938049,0.00115164,0.00149898,0.00197796,0.00265705,0.00364565,0.00409399,0.00831972,0.0170094};


	VEC pt27GeV  = {0.3, 0.5, 0.7, 0.9, 1.1, 1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0};
	VEC pt27GeV_err = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	VEC GeV27_v2_10_20  = { 0.0207,0.0383,0.0538,0.0683,0.0818,0.0949,0.1069,0.1178,0.1289,0.1359,0.1411,0.1454,0.1562,0.1555 };
	VEC GeV27_v2_10_20_err = {.0000,.0000,.0001,.0001,.0001,.0001,.0002,.0002,.0003,.0004,.0006,.0008,.0011,.0011};

	VEC GeV27_v2_20_30  = {0.0269,0.0496,0.0701,0.0887,0.1066,0.1232,0.1380,0.1521,0.1619,0.1736,0.1800,0.1883,0.1939,0.2015};
	VEC GeV27_v2_20_30_err = {0.0001,0.0001,0.0001,0.0002,0.0002,0.0003,0.0004,0.0006,0.0008,0.0010,0.0014,0.0020,0.0027,0.0029};

	VEC GeV27_v2_30_40  = {0.0299,0.0559,0.0794,0.1012,0.1214,0.1399,0.1565,0.1705,0.1822,0.1929,0.2032,0.2068,0.2214,0.2202};
	VEC GeV27_v2_30_40_err = {0.0001,0.0001,0.0002,0.0002,0.0003,0.0004,0.0006,0.0008,0.0011,0.0014,0.0019,0.0028,0.0038,0.0041};

	VEC GeV27_v2_50_60  = {0.0289,0.0558,0.0803,0.1042,0.1232,0.1399,0.1566,0.1697,0.1795,0.1928,0.2001,0.2100,0.2075,0.2243};
	VEC GeV27_v2_50_60_err = {0.0003,0.0004,0.0005,0.0007,0.0009,0.0013,0.0017,0.0024,0.0032,0.0044,0.0061,0.0087,0.0119,0.0131};


	if(Energy==27){
		/*
		DrawArticle( Energy, GeV27_pt_bin_center_Ant_0_80, GeV27_v2_values_Ant_0_80, GeV27_v2_stat_error_Ant_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Neg", 0.33,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_Pr_0_80, GeV27_v2_values_Pr_0_80, GeV27_v2_stat_error_Pr_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_KP_0_80, GeV27_v2_values_KP_0_80, GeV27_v2_stat_error_KP_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Pos", 0.28,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_KM_0_80, GeV27_v2_values_KM_0_80, GeV27_v2_stat_error_KM_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Neg", 0.28,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_PiP_0_80, GeV27_v2_values_PiP_0_80, GeV27_v2_stat_error_PiP_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Pos", 0.28,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_PiM_0_80, GeV27_v2_values_PiM_0_80, GeV27_v2_stat_error_PiM_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Neg", 0.28,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_Ant_0_10, GeV27_v2_values_Ant_0_10, GeV27_v2_stat_error_Ant_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Neg", 0.22,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_Pr_0_10, GeV27_v2_values_Pr_0_10, GeV27_v2_stat_error_Pr_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Pos", 0.22,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_KP_0_10, GeV27_v2_values_KP_0_10, GeV27_v2_stat_error_KP_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Pos", 0.16,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_KM_0_10, GeV27_v2_values_KM_0_10, GeV27_v2_stat_error_KM_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Neg", 0.16,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_PiP_0_10, GeV27_v2_values_PiP_0_10, GeV27_v2_stat_error_PiP_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Pos", 0.14,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_PiM_0_10, GeV27_v2_values_PiM_0_10, GeV27_v2_stat_error_PiM_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Neg", 0.14,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_Pr_10_40, GeV27_v2_values_Pr_10_40, GeV27_v2_stat_error_Pr_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_Ant_10_40, GeV27_v2_values_Ant_10_40, GeV27_v2_stat_error_Ant_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Neg", 0.33,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_KP_10_40, GeV27_v2_values_KP_10_40, GeV27_v2_stat_error_KP_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Pos", 0.28,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_KM_10_40, GeV27_v2_values_KM_10_40, GeV27_v2_stat_error_KM_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Neg", 0.28,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_PiP_10_40, GeV27_v2_values_PiP_10_40, GeV27_v2_stat_error_PiP_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Pos", 0.25,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_PiM_10_40, GeV27_v2_values_PiM_10_40, GeV27_v2_stat_error_PiM_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Neg", 0.25,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_Pr_40_80, GeV27_v2_values_Pr_40_80, GeV27_v2_stat_error_Pr_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Pos", 0.42,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_Ant_40_80, GeV27_v2_values_Ant_40_80, GeV27_v2_stat_error_Ant_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Neg", 0.42,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_KP_40_80, GeV27_v2_values_KP_40_80, GeV27_v2_stat_error_KP_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_KM_40_80, GeV27_v2_values_KM_40_80, GeV27_v2_stat_error_KM_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Neg", 0.33,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_PiP_40_80, GeV27_v2_values_PiP_40_80, GeV27_v2_stat_error_PiP_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV27_pt_bin_center_PiM_40_80, GeV27_v2_values_PiM_40_80, GeV27_v2_stat_error_PiM_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Neg", 0.33,3.6);
		*/
		DrawArticle( Energy, pt27GeV, GeV27_v2_10_20, GeV27_v2_10_20_err, pt27GeV_err, 6, 6, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, pt27GeV, GeV27_v2_20_30, GeV27_v2_20_30_err, pt27GeV_err, 5, 5, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, pt27GeV, GeV27_v2_30_40, GeV27_v2_30_40_err, pt27GeV_err, 4, 4, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, pt27GeV, GeV27_v2_50_60, GeV27_v2_50_60_err, pt27GeV_err, 2, 2, "Eta15","Hadrons", "", 0.33,3.6);
		/*

		std::vector<VEC> pt_par_0_80 = {GeV27_pt_bin_center_PiP_0_80, GeV27_pt_bin_center_KP_0_80, GeV27_pt_bin_center_Pr_0_80};
		std::vector<VEC> vn_par_0_80 = {GeV27_v2_values_PiP_0_80, GeV27_v2_values_KP_0_80, GeV27_v2_values_Pr_0_80};
		std::vector<VEC> pt_err_par_0_80 = {pt_bin_err,pt_bin_err,pt_bin_err};
		std::vector<VEC> vn_err_par_0_80 = {GeV27_v2_stat_error_PiP_0_80, GeV27_v2_stat_error_KP_0_80, GeV27_v2_stat_error_Pr_0_80};
		
		std::vector<VEC> pt_ant_0_80 = {GeV27_pt_bin_center_PiM_0_80, GeV27_pt_bin_center_KM_0_80, GeV27_pt_bin_center_Ant_0_80};
		std::vector<VEC> vn_ant_0_80 = {GeV27_v2_values_PiM_0_80, GeV27_v2_values_KM_0_80, GeV27_v2_values_Ant_0_80};
		std::vector<VEC> pt_err_ant_0_80 = {pt_bin_err,pt_bin_err,pt_bin_err};
		std::vector<VEC> vn_err_ant_0_80 = {GeV27_v2_stat_error_PiM_0_80, GeV27_v2_stat_error_KM_0_80, GeV27_v2_stat_error_Ant_0_80};

   		ParAntFlowVsPtForBESarticle(0, 8, "Eta01",-0.04, 3.16,pt_par_0_80,vn_par_0_80,pt_err_par_0_80,vn_err_par_0_80,pt_ant_0_80,vn_ant_0_80,pt_err_ant_0_80,vn_err_ant_0_80);
	
		std::vector<VEC> pt_par_0_10 = {GeV27_pt_bin_center_PiP_0_10, GeV27_pt_bin_center_KP_0_10, GeV27_pt_bin_center_Pr_0_10};
		std::vector<VEC> vn_par_0_10 = {GeV27_v2_values_PiP_0_10, GeV27_v2_values_KP_0_10, GeV27_v2_values_Pr_0_10};
		std::vector<VEC> pt_err_par_0_10 = {pt_bin_err,pt_bin_err,pt_bin_err};
		std::vector<VEC> vn_err_par_0_10 = {GeV27_v2_stat_error_PiP_0_10, GeV27_v2_stat_error_KP_0_10, GeV27_v2_stat_error_Pr_0_10};
		
		std::vector<VEC> pt_ant_0_10 = {GeV27_pt_bin_center_PiM_0_10, GeV27_pt_bin_center_KM_0_10, GeV27_pt_bin_center_Ant_0_10};
		std::vector<VEC> vn_ant_0_10 = {GeV27_v2_values_PiM_0_10, GeV27_v2_values_KM_0_10, GeV27_v2_values_Ant_0_10};
		std::vector<VEC> pt_err_ant_0_10 = {pt_bin_err,pt_bin_err,pt_bin_err};
		std::vector<VEC> vn_err_ant_0_10 = {GeV27_v2_stat_error_PiM_0_10, GeV27_v2_stat_error_KM_0_10, GeV27_v2_stat_error_Ant_0_10};

   		ParAntFlowVsPtForBESarticle(7, 8, "Eta01",-0.04, 3.16,pt_par_0_10,vn_par_0_10,pt_err_par_0_10,vn_err_par_0_10,pt_ant_0_10,vn_ant_0_10,pt_err_ant_0_10,vn_err_ant_0_10);
	
		std::vector<VEC> pt_par_10_40 = {GeV27_pt_bin_center_PiP_10_40, GeV27_pt_bin_center_KP_10_40, GeV27_pt_bin_center_Pr_10_40};
		std::vector<VEC> vn_par_10_40 = {GeV27_v2_values_PiP_10_40, GeV27_v2_values_KP_10_40, GeV27_v2_values_Pr_10_40};
		std::vector<VEC> pt_err_par_10_40 = {pt_bin_err,pt_bin_err,pt_bin_err};
		std::vector<VEC> vn_err_par_10_40 = {GeV27_v2_stat_error_PiP_10_40, GeV27_v2_stat_error_KP_10_40, GeV27_v2_stat_error_Pr_10_40};
		
		std::vector<VEC> pt_ant_10_40 = {GeV27_pt_bin_center_PiM_10_40, GeV27_pt_bin_center_KM_10_40, GeV27_pt_bin_center_Ant_10_40};
		std::vector<VEC> vn_ant_10_40 = {GeV27_v2_values_PiM_10_40, GeV27_v2_values_KM_10_40, GeV27_v2_values_Ant_10_40};
		std::vector<VEC> pt_err_ant_10_40 = {pt_bin_err,pt_bin_err,pt_bin_err};
		std::vector<VEC> vn_err_ant_10_40 = {GeV27_v2_stat_error_PiM_10_40, GeV27_v2_stat_error_KM_10_40, GeV27_v2_stat_error_Ant_10_40};

   		ParAntFlowVsPtForBESarticle(4, 6, "Eta01",-0.04, 3.16,pt_par_10_40,vn_par_10_40,pt_err_par_10_40,vn_err_par_10_40,pt_ant_10_40,vn_ant_10_40,pt_err_ant_10_40,vn_err_ant_10_40);
	
		std::vector<VEC> pt_par_40_80 = {GeV27_pt_bin_center_PiP_40_80, GeV27_pt_bin_center_KP_40_80, GeV27_pt_bin_center_Pr_40_80};
		std::vector<VEC> vn_par_40_80 = {GeV27_v2_values_PiP_40_80, GeV27_v2_values_KP_40_80, GeV27_v2_values_Pr_40_80};
		std::vector<VEC> pt_err_par_40_80 = {pt_bin_err,pt_bin_err,pt_bin_err};
		std::vector<VEC> vn_err_par_40_80 = {GeV27_v2_stat_error_PiP_40_80, GeV27_v2_stat_error_KP_40_80, GeV27_v2_stat_error_Pr_40_80};
		
		std::vector<VEC> pt_ant_40_80 = {GeV27_pt_bin_center_PiM_40_80, GeV27_pt_bin_center_KM_40_80, GeV27_pt_bin_center_Ant_40_80};
		std::vector<VEC> vn_ant_40_80 = {GeV27_v2_values_PiM_40_80, GeV27_v2_values_KM_40_80, GeV27_v2_values_Ant_40_80};
		std::vector<VEC> pt_err_ant_40_80 = {pt_bin_err,pt_bin_err,pt_bin_err};
		std::vector<VEC> vn_err_ant_40_80 = {GeV27_v2_stat_error_PiM_40_80, GeV27_v2_stat_error_KM_40_80, GeV27_v2_stat_error_Ant_40_80};

   		ParAntFlowVsPtForBESarticle(0, 3, "Eta01",-0.04, 3.16,pt_par_40_80,vn_par_40_80,pt_err_par_40_80,vn_err_par_40_80,pt_ant_40_80,vn_ant_40_80,pt_err_ant_40_80,vn_err_ant_40_80);
		*/

	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////19.6 GeV /////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////// 7.7 GeV //////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	VEC pt7GeV  = {0.3, 0.5, 0.7, 0.9, 1.1, 1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0};
	VEC pt7GeV_err = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,};

	Double_t v2_cent_7GeV[9] = {0.01459, 0.02467, 0.03611, 0.04675, 0.05162, 0.05127, 0.04743, 0.03990, 0.04460};
	Double_t v2_cent_7GeV_err[9] = {0.00056,0.00041,0.00025,0.00029,0.00039,0.00062,0.00113,0.00239,0.00401};

	VEC GeV7_v2_10_20  = {0.0201,0.0322,0.0425,0.0521,0.0659,0.0763,0.0911,0.0992,0.1115,0.1077,0.1346,0.1392,0.1265,0.1416};
	VEC GeV7_v2_10_20_err = {0.0002,0.0002,0.0003,0.0004,0.0005,0.0007,0.0010,0.0013,0.0018,0.0025,0.0035,0.0051,0.0070,0.0080};

	VEC GeV7_v2_20_30  = {0.0256,0.0424,0.0573,0.0704,0.0869,0.1013,0.1161,0.1295,0.1329,0.1526,0.1606,0.1782,0.1965,0.1783};//
	VEC GeV7_v2_20_30_err = {0.0005,0.0006,0.0008,0.0010,0.0013,0.0017,0.0023,0.0032,0.0044,0.0062,0.0088,0.0131,0.0183,0.0215};

	VEC GeV7_v2_30_40  = {0.0288,0.0473,0.0652,0.0803,0.0956,0.1150,0.1274,0.1470,0.1466,0.1656,0.1811,0.1945,0.1949,0.2287};
	VEC GeV7_v2_30_40_err = {0.0006,0.0008,0.0010,0.0013,0.0018,0.0024,0.0033,0.0045,0.0063,0.0090,0.0129,0.0192,0.0279,0.0326};

	VEC GeV7_v2_50_60  = {0.0280,0.0454,0.0615,0.0812,0.0976,0.1015,0.1118,0.1150,0.1858,0.1297,0.1948,0.0653,0.1988,0.0876};
	VEC GeV7_v2_50_60_err = {0.0017,0.0022,0.0030,0.0041,0.0057,0.0078,0.0111,0.0159,0.0234,0.0341,0.0505,0.0819,0.1197,0.1450};
	// ****************** beam Energy: 11.5 GeV ******************
	// Event plane method: eta-sub

	// ---------------- Particle species: Proton ----------------
	VEC GeV7_pt_bin_center_Pr_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.835};
	VEC GeV7_v2_values_Pr_0_80 = {0.00691689,0.0161549,0.0312339,0.0462411,0.0616304,0.0705497,0.0824948,0.0931209,0.0985242,0.0993178,0.107344,0.105471,0.0974711};
	VEC GeV7_v2_stat_error_Pr_0_80 = {0.00306278,0.00118166,0.00102077,0.00106674,0.00124804,0.00155582,0.00203117,0.00275382,0.00383099,0.00528279,0.00751169,0.0110623,0.0138618};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV7_pt_bin_center_Ant_0_80 = {0.465,0.795,1.155,1.635};
	VEC GeV7_v2_values_Ant_0_80 = {0.00686059,0.0109067,0.0131317,0.0637152};
	VEC GeV7_v2_stat_error_Ant_0_80 = {0.0141331,0.00971777,0.0138787,0.0205911};

	// ---------------- Particle species: PiM ----------------
	VEC GeV7_pt_bin_center_PiM_0_80 = {0.285,0.495,0.675,0.885,1.155,1.545,1.935,2.325};
	VEC GeV7_v2_values_PiM_0_80 = {0.0227645,0.0392692,0.0528099,0.0646252,0.0741718,0.0853517,0.100855,0.112864};
	VEC GeV7_v2_stat_error_PiM_0_80 = {0.000504488,0.000661156,0.000954283,0.00142493,0.00180491,0.00419101,0.0100161,0.0232811};

	// ---------------- Particle species: PiP ----------------
	VEC GeV7_pt_bin_center_PiP_0_80 = {0.285,0.495,0.675,0.885,1.155,1.545,1.935,2.325};
	VEC GeV7_v2_values_PiP_0_80 = {0.0190878,0.0361773,0.0478764,0.0578427,0.0667559,0.0828597,0.0895746,0.0515206};
	VEC GeV7_v2_stat_error_PiP_0_80 = {0.00052733,0.000676753,0.000966526,0.00143331,0.00180211,0.00413486,0.00979729,0.0216026};

	// ---------------- Particle species: KM ----------------
	VEC GeV7_pt_bin_center_KM_0_80 = {0.315,0.495,0.705,0.885,1.155,1.545,1.935,2.325};
	VEC GeV7_v2_values_KM_0_80 = {0.0109867,0.0225145,0.0310003,0.046957,0.0530698,0.0739074,0.0102221,-0.0178917};
	VEC GeV7_v2_stat_error_KM_0_80 = {0.00451129,0.00304574,0.00329745,0.00418776,0.0046704,0.00963526,0.0231336,0.0702516};

	// ---------------- Particle species: KP ----------------
	VEC GeV7_pt_bin_center_KP_0_80 = {0.315,0.495,0.705,0.885,1.155,1.545,1.935,2.325};
	VEC GeV7_v2_values_KP_0_80 = {0.010691,0.0263934,0.0397986,0.0537884,0.0622878,0.0765433,0.077197,0.0529762};
	VEC GeV7_v2_stat_error_KP_0_80 = {0.00280781,0.00182188,0.0019098,0.00236217,0.0025566,0.00510552,0.01125,0.0305002};

	// ---------------- Particle species: Proton ----------------
	VEC GeV7_pt_bin_center_Pr_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.145,2.625};
	VEC GeV7_v2_values_Pr_0_10 = {-0.000269201,0.00292256,0.0126519,0.0210576,0.0261742,0.0307903,0.0367998,0.0400284,0.0387608,0.0421431,0.0208087};
	VEC GeV7_v2_stat_error_Pr_0_10 = {0.00686715,0.00260173,0.00219494,0.00224465,0.00258058,0.00317244,0.0040983,0.00551364,0.00762922,0.0085211,0.017097};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV7_pt_bin_center_Ant_0_10 = {0.765,1.485};
	VEC GeV7_v2_values_Ant_0_10 = {-0.00707938,-0.00320015};
	VEC GeV7_v2_stat_error_Ant_0_10 = {0.0162858,0.0293924};

	// ---------------- Particle species: PiM ----------------
	VEC GeV7_pt_bin_center_PiM_0_10 = {0.375,0.735,1.125,1.545};
	VEC GeV7_v2_values_PiM_0_10 = {0.0177992,0.0314537,0.0367388,0.0391578};
	VEC GeV7_v2_stat_error_PiM_0_10 = {0.000866997,0.00168551,0.00375154,0.00859142};

	// ---------------- Particle species: PiP ----------------
	VEC GeV7_pt_bin_center_PiP_0_10 = {0.375,0.735,1.125,1.545};
	VEC GeV7_v2_values_PiP_0_10 = {0.0154442,0.0269776,0.0283281,0.037249};
	VEC GeV7_v2_stat_error_PiP_0_10 = {0.000900857,0.00169369,0.00371201,0.00835521};

	// ---------------- Particle species: KM ----------------
	VEC GeV7_pt_bin_center_KM_0_10 = {0.405,0.795,1.155,1.545};
	VEC GeV7_v2_values_KM_0_10 = {0.0152733,0.0136827,0.0278363,0.0242437};
	VEC GeV7_v2_stat_error_KM_0_10 = {0.00537366,0.00537836,0.00945312,0.0191164};

	// ---------------- Particle species: KP ----------------
	VEC GeV7_pt_bin_center_KP_0_10 = {0.405,0.795,1.155,1.545};
	VEC GeV7_v2_values_KP_0_10 = {0.0156202,0.022867,0.0325253,0.029094};
	VEC GeV7_v2_stat_error_KP_0_10 = {0.00324299,0.00308078,0.00518859,0.0102165};

	// ---------------- Particle species: Proton ----------------
	VEC GeV7_pt_bin_center_Pr_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.625};
	VEC GeV7_v2_values_Pr_10_40 = {0.00643856,0.0196543,0.0362226,0.0541206,0.0781074,0.090889,0.10682,0.121606,0.131171,0.13338,0.14221,0.152974};
	VEC GeV7_v2_stat_error_Pr_10_40 = {0.0028076,0.00108471,0.00093993,0.000987168,0.00115896,0.00144923,0.00189495,0.00257513,0.00358486,0.00494914,0.0070575,0.00820795};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV7_pt_bin_center_Ant_10_40 = {0.465,0.795,1.155,1.635};
	VEC GeV7_v2_values_Ant_10_40 = {-0.00393915,-0.0101251,0.0242419,0.0750553};
	VEC GeV7_v2_stat_error_Ant_10_40 = {0.0128592,0.00875286,0.0125874,0.0188931};

	// ---------------- Particle species: PiM ----------------
	VEC GeV7_pt_bin_center_PiM_10_40 = {0.285,0.495,0.675,0.885,1.155,1.545,1.935};
	VEC GeV7_v2_values_PiM_10_40 = {0.0264438,0.0453791,0.0613248,0.0752108,0.0886761,0.110444,0.110866};
	VEC GeV7_v2_stat_error_PiM_10_40 = {0.000459202,0.000602041,0.00087003,0.00129944,0.00164998,0.00383577,0.00917259};

	// ---------------- Particle species: PiP ----------------
	VEC GeV7_pt_bin_center_PiP_10_40 = {0.285,0.495,0.675,0.885,1.155,1.545,1.935};
	VEC GeV7_v2_values_PiP_10_40 = {0.0225062,0.0415791,0.0560643,0.0678896,0.0819491,0.100776,0.120508};
	VEC GeV7_v2_stat_error_PiP_10_40 = {0.000479722,0.000616763,0.000883004,0.00131197,0.00165343,0.00380281,0.00915605};

	// ---------------- Particle species: KM ----------------
	VEC GeV7_pt_bin_center_KM_10_40 = {0.315,0.495,0.705,0.885,1.155,1.545,1.935};
	VEC GeV7_v2_values_KM_10_40 = {0.00187643,0.027042,0.0401026,0.0520291,0.0636185,0.0922199,0.0851067};
	VEC GeV7_v2_stat_error_KM_10_40 = {0.00416041,0.00282332,0.00306204,0.00390723,0.0043618,0.00900285,0.0219624};

	// ---------------- Particle species: KP ----------------
	VEC GeV7_pt_bin_center_KP_10_40 = {0.315,0.495,0.705,0.885,1.155,1.545,1.965};
	VEC GeV7_v2_values_KP_10_40 = {0.0123349,0.0284567,0.0493235,0.0674997,0.080886,0.100798,0.0878592};
	VEC GeV7_v2_stat_error_KP_10_40 = {0.00259688,0.00169075,0.00177886,0.00220443,0.00239161,0.00478715,0.0106033};

	// ---------------- Particle species: Proton ----------------
	VEC GeV7_pt_bin_center_Pr_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.755,2.175};
	VEC GeV7_v2_values_Pr_40_80 = {0.0193443,0.0292042,0.0538348,0.0808153,0.0949687,0.107889,0.129731,0.156801,0.181197};
	VEC GeV7_v2_stat_error_Pr_40_80 = {0.0100635,0.00405841,0.00373599,0.00414954,0.00513363,0.00670999,0.00912322,0.0105659,0.0203899};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV7_pt_bin_center_Ant_40_80 = {0.645,1.245};
	VEC GeV7_v2_values_Ant_40_80 = {0.0747818,0.0396689};
	VEC GeV7_v2_stat_error_Ant_40_80 = {0.0246589,0.0437798};

	// ---------------- Particle species: PiM ----------------
	VEC GeV7_pt_bin_center_PiM_40_80 = {0.285,0.495,0.675,0.885,1.125,1.545};
	VEC GeV7_v2_values_PiM_40_80 = {0.0257936,0.0513509,0.0728293,0.0895198,0.110976,0.100925};
	VEC GeV7_v2_stat_error_PiM_40_80 = {0.00177321,0.00235618,0.00345922,0.00532486,0.00703612,0.0174535};

	// ---------------- Particle species: PiP ----------------
	VEC GeV7_pt_bin_center_PiP_40_80 = {0.285,0.495,0.675,0.885,1.125,1.545};
	VEC GeV7_v2_values_PiP_40_80 = {0.0224911,0.0479004,0.0653419,0.0906055,0.104877,0.142921};
	VEC GeV7_v2_stat_error_PiP_40_80 = {0.00184114,0.00242319,0.00353726,0.00542294,0.00715068,0.0176374};

	// ---------------- Particle species: KM ----------------
	VEC GeV7_pt_bin_center_KM_40_80 = {0.315,0.495,0.675,0.885,1.155,1.545};
	VEC GeV7_v2_values_KM_40_80 = {0.017888,0.0281836,0.0580301,0.0933571,0.0902883,0.185715};
	VEC GeV7_v2_stat_error_KM_40_80 = {0.0161862,0.0114554,0.0130659,0.0174445,0.0203021,0.0456036};

	// ---------------- Particle species: KP ----------------
	VEC GeV7_pt_bin_center_KP_40_80 = {0.315,0.495,0.705,0.885,1.155,1.545};
	VEC GeV7_v2_values_KP_40_80 = {0.00925424,0.0395317,0.0565827,0.0779383,0.0770185,0.142544};
	VEC GeV7_v2_stat_error_KP_40_80 = {0.0102705,0.00700009,0.00766815,0.00984592,0.0111964,0.0238864};

	if(Energy==7){
		DrawArticle( Energy, GeV7_pt_bin_center_Ant_0_80, GeV7_v2_values_Ant_0_80, GeV7_v2_stat_error_Ant_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Neg", 0.33,3.0);
		DrawArticle( Energy, GeV7_pt_bin_center_Pr_0_80, GeV7_v2_values_Pr_0_80, GeV7_v2_stat_error_Pr_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV7_pt_bin_center_KP_0_80, GeV7_v2_values_KP_0_80, GeV7_v2_stat_error_KP_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Pos", 0.28,3.0);
		DrawArticle( Energy, GeV7_pt_bin_center_KM_0_80, GeV7_v2_values_KM_0_80, GeV7_v2_stat_error_KM_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Neg", 0.28,3.0);
		DrawArticle( Energy, GeV7_pt_bin_center_PiP_0_80, GeV7_v2_values_PiP_0_80, GeV7_v2_stat_error_PiP_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Pos", 0.28,3.0);
		DrawArticle( Energy, GeV7_pt_bin_center_PiM_0_80, GeV7_v2_values_PiM_0_80, GeV7_v2_stat_error_PiM_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Neg", 0.28,3.0);
		DrawArticle( Energy, GeV7_pt_bin_center_Ant_0_10, GeV7_v2_values_Ant_0_10, GeV7_v2_stat_error_Ant_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Neg", 0.22,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_Pr_0_10, GeV7_v2_values_Pr_0_10, GeV7_v2_stat_error_Pr_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Pos", 0.22,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_KP_0_10, GeV7_v2_values_KP_0_10, GeV7_v2_stat_error_KP_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Pos", 0.16,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_KM_0_10, GeV7_v2_values_KM_0_10, GeV7_v2_stat_error_KM_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Neg", 0.16,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_PiP_0_10, GeV7_v2_values_PiP_0_10, GeV7_v2_stat_error_PiP_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Pos", 0.14,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_PiM_0_10, GeV7_v2_values_PiM_0_10, GeV7_v2_stat_error_PiM_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Neg", 0.14,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_Pr_10_40, GeV7_v2_values_Pr_10_40, GeV7_v2_stat_error_Pr_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Pos", 0.33,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_Ant_10_40, GeV7_v2_values_Ant_10_40, GeV7_v2_stat_error_Ant_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Neg", 0.33,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_KP_10_40, GeV7_v2_values_KP_10_40, GeV7_v2_stat_error_KP_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Pos", 0.28,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_KM_10_40, GeV7_v2_values_KM_10_40, GeV7_v2_stat_error_KM_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Neg", 0.28,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_PiP_10_40, GeV7_v2_values_PiP_10_40, GeV7_v2_stat_error_PiP_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Pos", 0.25,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_PiM_10_40, GeV7_v2_values_PiM_10_40, GeV7_v2_stat_error_PiM_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Neg", 0.25,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_Pr_40_80, GeV7_v2_values_Pr_40_80, GeV7_v2_stat_error_Pr_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Pos", 0.42,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_Ant_40_80, GeV7_v2_values_Ant_40_80, GeV7_v2_stat_error_Ant_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Neg", 0.42,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_KP_40_80, GeV7_v2_values_KP_40_80, GeV7_v2_stat_error_KP_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Pos", 0.33,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_KM_40_80, GeV7_v2_values_KM_40_80, GeV7_v2_stat_error_KM_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Neg", 0.33,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_PiP_40_80, GeV7_v2_values_PiP_40_80, GeV7_v2_stat_error_PiP_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Pos", 0.33,2.5);
		DrawArticle( Energy, GeV7_pt_bin_center_PiM_40_80, GeV7_v2_values_PiM_40_80, GeV7_v2_stat_error_PiM_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Neg", 0.33,2.5);

		DrawArticle( Energy, pt7GeV, GeV7_v2_10_20, GeV7_v2_10_20_err, ptGeV_err, 6, 6, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, pt7GeV, GeV7_v2_20_30, GeV7_v2_20_30_err, ptGeV_err, 5, 5, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, pt7GeV, GeV7_v2_30_40, GeV7_v2_30_40_err, ptGeV_err, 4, 4, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, pt7GeV, GeV7_v2_50_60, GeV7_v2_50_60_err, ptGeV_err, 2, 2, "Eta15","Hadrons", "", 0.33,3.6);

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////// 11.5 GeV //////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	VEC pt11GeV  = {0.3, 0.5, 0.7, 0.9, 1.1, 1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0};
	VEC pt11GeV_err = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	Double_t v2_cent_11GeV[9] = {0.01748, 0.02711, 0.03932, 0.05021, 0.05554, 0.05570, 0.05206, 0.04675, 0.04531};
	Double_t v2_cent_11GeV_err[9] = {0.00023, 0.00019, 0.00012, 0.00014, 0.00019, 0.00029, 0.00052, 0.00108, 0.00225};

	VEC GeV11_v2_10_20  = {0.0211,0.0356,0.0483,0.0612,0.0736,0.0855,0.0986,0.1061,0.1178,0.1297,0.1351,0.1396,0.1463,0.1523};
	VEC GeV11_v2_10_20_err = {0.0002,0.0002,0.0003,0.0004,0.0005,0.0007,0.0010,0.0013,0.0018,0.0025,0.0035,0.0051,0.0070,0.0080};

	VEC GeV11_v2_20_30  = {0.0263,0.0462,0.0632,0.0796,0.0967,0.1135,0.1280,0.1436,0.1528,0.1669,0.1753,0.1842,0.1890,0.2095};//
	VEC GeV11_v2_20_30_err = {0.0002,0.0002,0.0003,0.0004,0.0005,0.0006,0.0008,0.0011,0.0016,0.0022,0.0030,0.0042,0.0059,0.0067};

	VEC GeV11_v2_30_40  = {0.0295,0.0514,0.0717,0.0908,0.1097,0.1273,0.1450,0.1585,0.1740,0.1810,0.1907,0.1962,0.2224,0.2169};
	VEC GeV11_v2_30_40_err = {0.0003,0.0004,0.0005,0.0007,0.0009,0.0012,0.0016,0.0022,0.0031,0.0044,0.0061,0.0091,0.0128,0.0148};

	VEC GeV11_v2_50_60  = {0.0287,0.0508,0.0694,0.0916,0.1090,0.1223,0.1338,0.1471,0.1819,0.1567,0.2149,0.2111,0.1518,0.2272};
	VEC GeV11_v2_50_60_err = {0.0008,0.0010,0.0014,0.0019,0.0027,0.0037,0.0052,0.0074,0.0105,0.0151,0.0220,0.0333,0.0486,0.0575};
	// ****************** beam Energy: 11.5 GeV ******************
	// Event plane method: eta-sub

	// ---------------- Particle species: Proton ----------------
	VEC GeV11_pt_bin_center_Pr_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.805,3.435};
	VEC GeV11_v2_values_Pr_0_80 = {0.0120852,0.0155901,0.0334294,0.0499063,0.0690335,0.0836196,0.0989993,0.111464,0.119578,0.126679,0.13911,0.132912,0.153769,0.163517};
	VEC GeV11_v2_stat_error_Pr_0_80 = {0.0017164,0.000653056,0.000558847,0.00057748,0.00066537,0.000816543,0.00104704,0.00139407,0.00190311,0.00254527,0.0035386,0.00507522,0.00583886,0.0199232};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV11_pt_bin_center_Ant_0_80 = {0.315,0.495,0.705,0.885,1.095,1.365,1.755,2.145,2.655};
	VEC GeV11_v2_values_Ant_0_80 = {-0.0131254,-0.00122946,0.0181017,0.0268879,0.0428292,0.0656788,0.0875492,0.115173,0.123109};
	VEC GeV11_v2_stat_error_Ant_0_80 = {0.00937928,0.00394889,0.00337106,0.00354992,0.00419587,0.00428718,0.00705921,0.0138731,0.0305049};

	// ---------------- Particle species: PiM ----------------
	VEC GeV11_pt_bin_center_PiM_0_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.935,2.325,2.805};
	VEC GeV11_v2_values_PiM_0_80 = {0.0233576,0.0402824,0.055344,0.0714658,0.0808871,0.0905686,0.0977469,0.102013,0.105191,0.0963362,0.116218};
	VEC GeV11_v2_stat_error_PiM_0_80 = {0.000219189,0.000278528,0.000393855,0.000576135,0.000854327,0.00125873,0.00186999,0.00279223,0.0034078,0.00722939,0.0145114};

	// ---------------- Particle species: PiP ----------------
	VEC GeV11_pt_bin_center_PiP_0_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.935,2.355,2.805};
	VEC GeV11_v2_values_PiP_0_80 = {0.020229,0.0381335,0.0525583,0.0661144,0.0760887,0.0849049,0.0931519,0.0947235,0.107344,0.108577,0.141732};
	VEC GeV11_v2_stat_error_PiP_0_80 = {0.000225845,0.000283033,0.000397265,0.000578367,0.000852085,0.00125197,0.00185568,0.00277691,0.00342084,0.00702053,0.0139547};

	// ---------------- Particle species: KM ----------------
	VEC GeV11_pt_bin_center_KM_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.665,1.935,2.325,2.775};
	VEC GeV11_v2_values_KM_0_80 = {0.0113204,0.025042,0.0396243,0.0551058,0.0674386,0.0783379,0.0898145,0.0933213,0.0975876,0.0899405,0.104651};
	VEC GeV11_v2_stat_error_KM_0_80 = {0.0017241,0.00114052,0.00120486,0.00149062,0.001986,0.00269164,0.00379874,0.0055761,0.00704378,0.0176145,0.0408336};

	// ---------------- Particle species: KP ----------------
	VEC GeV11_pt_bin_center_KP_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.665,1.965,2.325,2.775};
	VEC GeV11_v2_values_KP_0_80 = {0.0138815,0.0284431,0.042206,0.0550338,0.0676825,0.0778151,0.0871085,0.0877634,0.0856474,0.0843797,0.0440532};
	VEC GeV11_v2_stat_error_KP_0_80 = {0.00123278,0.000793404,0.000820618,0.000995769,0.00130208,0.00174595,0.00242321,0.00343716,0.00406694,0.010179,0.033157};

	// ---------------- Particle species: Proton ----------------
	VEC GeV11_pt_bin_center_Pr_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805};
	VEC GeV11_v2_values_Pr_0_10 = {-0.000427349,0.0028464,0.0147849,0.0237616,0.0342703,0.0427704,0.0535675,0.0600192,0.0590241,0.0636355,0.0752996,0.0746523};
	VEC GeV11_v2_stat_error_Pr_0_10 = {0.00352098,0.00131852,0.00110456,0.00111843,0.00126779,0.00153667,0.00195472,0.00258985,0.00351736,0.00471402,0.0053951,0.0107001};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV11_pt_bin_center_Ant_0_10 = {0.465,0.795,1.245,1.845};
	VEC GeV11_v2_values_Ant_0_10 = {-0.0107589,0.00367793,0.0195834,0.0335654};
	VEC GeV11_v2_stat_error_Ant_0_10 = {0.00788978,0.00498118,0.00581923,0.0113896};

	// ---------------- Particle species: PiM ----------------
	VEC GeV11_pt_bin_center_PiM_0_10 = {0.285,0.495,0.675,0.885,1.155,1.545,1.935,2.355};
	VEC GeV11_v2_values_PiM_0_10 = {0.0162247,0.0254798,0.0314538,0.041526,0.0454085,0.0528672,0.0399127,0.0170665};
	VEC GeV11_v2_stat_error_PiM_0_10 = {0.00043004,0.000541769,0.000763113,0.00111113,0.00135423,0.00293635,0.00646275,0.0137012};

	// ---------------- Particle species: PiP ----------------
	VEC GeV11_pt_bin_center_PiP_0_10 = {0.285,0.495,0.675,0.885,1.155,1.545,1.935,2.355};
	VEC GeV11_v2_values_PiP_0_10 = {0.0123074,0.0221415,0.0293689,0.0332695,0.0383741,0.0452447,0.0500749,0.0426986};
	VEC GeV11_v2_stat_error_PiP_0_10 = {0.000444405,0.000549625,0.000766495,0.00110815,0.00134295,0.00288218,0.00638391,0.013107};

	// ---------------- Particle species: KM ----------------
	VEC GeV11_pt_bin_center_KM_0_10 = {0.315,0.495,0.705,0.885,1.155,1.545,1.935,2.325};
	VEC GeV11_v2_values_KM_0_10 = {0.00424446,0.015042,0.0252552,0.0307013,0.0401831,0.0540672,0.0588411,0.00606711};
	VEC GeV11_v2_stat_error_KM_0_10 = {0.00343309,0.00223039,0.00231722,0.00283951,0.00301211,0.00582309,0.0129096,0.0321052};

	// ---------------- Particle species: KP ----------------
	VEC GeV11_pt_bin_center_KP_0_10 = {0.315,0.495,0.705,0.885,1.155,1.545,1.935,2.325};
	VEC GeV11_v2_values_KP_0_10 = {0.00541312,0.0150757,0.0212609,0.0280978,0.0387445,0.0464258,0.0385034,0.0345413};
	VEC GeV11_v2_stat_error_KP_0_10 = {0.00244884,0.00154685,0.00157548,0.00189463,0.00197186,0.00370722,0.00765916,0.018113};

	// ---------------- Particle species: Proton ----------------
	VEC GeV11_pt_bin_center_Pr_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.685,3.045};
	VEC GeV11_v2_values_Pr_10_40 = {0.0120472,0.02022,0.0387996,0.0595139,0.0827635,0.102705,0.12077,0.135948,0.149776,0.157411,0.166608,0.174031,0.180611,0.199118};
	VEC GeV11_v2_stat_error_Pr_10_40 = {0.00170825,0.000652142,0.00056097,0.000582846,0.000674855,0.000830744,0.00106771,0.00142295,0.00194179,0.00259974,0.00362094,0.005181,0.00762002,0.00926275};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV11_pt_bin_center_Ant_10_40 = {0.315,0.495,0.705,0.885,1.095,1.365,1.755,2.205};
	VEC GeV11_v2_values_Ant_10_40 = {-0.00797661,-0.002644,0.0158792,0.0308917,0.0465937,0.0819031,0.101881,0.145431};
	VEC GeV11_v2_stat_error_Ant_10_40 = {0.00916032,0.00381289,0.00325416,0.003437,0.00408626,0.00420012,0.00704357,0.0129257};

	// ---------------- Particle species: PiM ----------------
	VEC GeV11_pt_bin_center_PiM_10_40 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.935,2.355};
	VEC GeV11_v2_values_PiM_10_40 = {0.0262688,0.0463362,0.0646118,0.0822958,0.0931763,0.106083,0.114697,0.123377,0.131838,0.126493};
	VEC GeV11_v2_stat_error_PiM_10_40 = {0.000216586,0.00027538,0.000389892,0.000570868,0.000847334,0.00124901,0.00185315,0.00276248,0.00338022,0.00715378};

	// ---------------- Particle species: PiP ----------------
	VEC GeV11_pt_bin_center_PiP_10_40 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.935,2.355};
	VEC GeV11_v2_values_PiP_10_40 = {0.0234571,0.0439218,0.0621634,0.0784782,0.0918088,0.103063,0.110708,0.120155,0.125528,0.127156};
	VEC GeV11_v2_stat_error_PiP_10_40 = {0.000223011,0.000280047,0.000393929,0.00057453,0.000847508,0.00124706,0.0018493,0.0027647,0.00341252,0.006978};

	// ---------------- Particle species: KM ----------------
	VEC GeV11_pt_bin_center_KM_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.665,1.935,2.325};
	VEC GeV11_v2_values_KM_10_40 = {0.013153,0.0283248,0.0439597,0.0640388,0.079496,0.0921375,0.107982,0.113846,0.121806,0.142012};
	VEC GeV11_v2_stat_error_KM_10_40 = {0.0017274,0.00114852,0.00121875,0.00151298,0.00201569,0.00274069,0.00387267,0.00570455,0.00718289,0.0192773};

	// ---------------- Particle species: KP ----------------
	VEC GeV11_pt_bin_center_KP_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.665,1.965,2.355};
	VEC GeV11_v2_values_KP_10_40 = {0.0147594,0.0330414,0.0518185,0.0691931,0.0822795,0.0936632,0.106092,0.110771,0.11221,0.126219};
	VEC GeV11_v2_stat_error_KP_10_40 = {0.00123896,0.000801614,0.000833044,0.00101324,0.00132714,0.0017836,0.00247474,0.00351944,0.0042141,0.00988277};

	// ---------------- Particle species: Proton ----------------
	VEC GeV11_pt_bin_center_Pr_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.535};
	VEC GeV11_v2_values_Pr_40_80 = {0.0328812,0.0253787,0.056364,0.0806939,0.114027,0.131065,0.151693,0.182512,0.202355,0.218501,0.236059,0.230715};
	VEC GeV11_v2_stat_error_Pr_40_80 = {0.00578687,0.00230746,0.00210492,0.00231838,0.00282459,0.00363044,0.00481897,0.00662249,0.00933959,0.0125801,0.017797,0.0219193};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV11_pt_bin_center_Ant_40_80 = {0.465,0.705,0.885,1.155,1.635};
	VEC GeV11_v2_values_Ant_40_80 = {0.00721496,0.0357429,0.0635109,0.0904557,0.187519};
	VEC GeV11_v2_stat_error_Ant_40_80 = {0.0104806,0.0106743,0.0122066,0.0125648,0.0205043};

	// ---------------- Particle species: PiM ----------------
	VEC GeV11_pt_bin_center_PiM_40_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.545,1.935};
	VEC GeV11_v2_values_PiM_40_80 = {0.0281522,0.0510958,0.0751902,0.09996,0.121391,0.133577,0.142167,0.166879};
	VEC GeV11_v2_stat_error_PiM_40_80 = {0.000792553,0.00102755,0.00147907,0.00222085,0.00335938,0.00507959,0.00644639,0.0146054};

	// ---------------- Particle species: PiP ----------------
	VEC GeV11_pt_bin_center_PiP_40_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.545,1.935};
	VEC GeV11_v2_values_PiP_40_80 = {0.0254122,0.0525961,0.0706521,0.0983145,0.10694,0.135477,0.139412,0.196748};
	VEC GeV11_v2_stat_error_PiP_40_80 = {0.000811388,0.00104739,0.00150464,0.00225616,0.00341003,0.00514678,0.00651868,0.0149829};

	// ---------------- Particle species: KM ----------------
	VEC GeV11_pt_bin_center_KM_40_80 = {0.315,0.495,0.675,0.885,1.065,1.275,1.545,1.935};
	VEC GeV11_v2_values_KM_40_80 = {0.0202858,0.0371685,0.0620032,0.0893671,0.103202,0.135015,0.133998,0.115089};
	VEC GeV11_v2_stat_error_KM_40_80 = {0.00629596,0.00438176,0.00485468,0.00626826,0.00862904,0.0120347,0.0144185,0.0350634};

	// ---------------- Particle species: KP ----------------
	VEC GeV11_pt_bin_center_KP_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.545,1.965};
	VEC GeV11_v2_values_KP_40_80 = {0.0297663,0.0449837,0.0620183,0.0761777,0.108157,0.118218,0.136687,0.132862};
	VEC GeV11_v2_stat_error_KP_40_80 = {0.0045977,0.00311375,0.00336253,0.00424135,0.0056802,0.00781436,0.00918004,0.0200752};

	if(Energy==11){
		DrawArticle( Energy, GeV11_pt_bin_center_Ant_0_80, GeV11_v2_values_Ant_0_80, GeV11_v2_stat_error_Ant_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Neg", 0.33,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_Pr_0_80, GeV11_v2_values_Pr_0_80, GeV11_v2_stat_error_Pr_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_KP_0_80, GeV11_v2_values_KP_0_80, GeV11_v2_stat_error_KP_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Pos", 0.28,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_KM_0_80, GeV11_v2_values_KM_0_80, GeV11_v2_stat_error_KM_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Neg", 0.28,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_PiP_0_80, GeV11_v2_values_PiP_0_80, GeV11_v2_stat_error_PiP_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Pos", 0.28,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_PiM_0_80, GeV11_v2_values_PiM_0_80, GeV11_v2_stat_error_PiM_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Neg", 0.28,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_Ant_0_10, GeV11_v2_values_Ant_0_10, GeV11_v2_stat_error_Ant_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Neg", 0.22,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_Pr_0_10, GeV11_v2_values_Pr_0_10, GeV11_v2_stat_error_Pr_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Pos", 0.22,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_KP_0_10, GeV11_v2_values_KP_0_10, GeV11_v2_stat_error_KP_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Pos", 0.16,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_KM_0_10, GeV11_v2_values_KM_0_10, GeV11_v2_stat_error_KM_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Neg", 0.16,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_PiP_0_10, GeV11_v2_values_PiP_0_10, GeV11_v2_stat_error_PiP_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Pos", 0.14,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_PiM_0_10, GeV11_v2_values_PiM_0_10, GeV11_v2_stat_error_PiM_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Neg", 0.14,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_Pr_10_40, GeV11_v2_values_Pr_10_40, GeV11_v2_stat_error_Pr_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_Ant_10_40, GeV11_v2_values_Ant_10_40, GeV11_v2_stat_error_Ant_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Neg", 0.33,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_KP_10_40, GeV11_v2_values_KP_10_40, GeV11_v2_stat_error_KP_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Pos", 0.28,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_KM_10_40, GeV11_v2_values_KM_10_40, GeV11_v2_stat_error_KM_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Neg", 0.28,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_PiP_10_40, GeV11_v2_values_PiP_10_40, GeV11_v2_stat_error_PiP_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Pos", 0.25,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_PiM_10_40, GeV11_v2_values_PiM_10_40, GeV11_v2_stat_error_PiM_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Neg", 0.25,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_Pr_40_80, GeV11_v2_values_Pr_40_80, GeV11_v2_stat_error_Pr_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Pos", 0.42,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_Ant_40_80, GeV11_v2_values_Ant_40_80, GeV11_v2_stat_error_Ant_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Neg", 0.42,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_KP_40_80, GeV11_v2_values_KP_40_80, GeV11_v2_stat_error_KP_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_KM_40_80, GeV11_v2_values_KM_40_80, GeV11_v2_stat_error_KM_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Neg", 0.33,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_PiP_40_80, GeV11_v2_values_PiP_40_80, GeV11_v2_stat_error_PiP_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV11_pt_bin_center_PiM_40_80, GeV11_v2_values_PiM_40_80, GeV11_v2_stat_error_PiM_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Neg", 0.33,3.0);
	
		DrawArticle( Energy, pt11GeV, GeV11_v2_10_20, GeV11_v2_10_20_err, ptGeV_err, 6, 6, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, pt11GeV, GeV11_v2_20_30, GeV11_v2_20_30_err, ptGeV_err, 5, 5, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, pt11GeV, GeV11_v2_30_40, GeV11_v2_30_40_err, ptGeV_err, 4, 4, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, pt11GeV, GeV11_v2_50_60, GeV11_v2_50_60_err, ptGeV_err, 2, 2, "Eta15","Hadrons", "", 0.33,3.6);

	}




	// ---------------- Particle species: Proton ----------------
	VEC GeV14_pt_bin_center_Pr_0_80 = {0.347289,0.515845,0.702633,0.89699,1.09391,1.292,1.49073,1.68987,1.88928,2.15596,2.59429};
	VEC GeV14_v2_values_Pr_0_80 = {0.00447695,0.0182758,0.0328822,0.0507055,0.0712785,0.0866777,0.102239,0.116333,0.127892,0.138503,0.148001};
	VEC GeV14_v2_stat_error_Pr_0_80 = {0.00180695,0.000600128,0.00048577,0.000491873,0.000560602,0.000679272,0.000861487,0.00113574,0.00157165,0.00194457,0.00504607};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV14_pt_bin_center_Ant_0_80 = {0.349959,0.516156,0.703589,0.897042,1.09333,1.29118,1.49003,1.68959,1.88966,2.16023,2.57829};
	VEC GeV14_v2_values_Ant_0_80 = {-0.0022704,0.00934183,0.0180102,0.0347868,0.056024,0.0690342,0.0798844,0.097556,0.113615,0.111802,0.147487};
	VEC GeV14_v2_stat_error_Ant_0_80 = {0.00828211,0.00272651,0.00215914,0.00217535,0.00248471,0.00302864,0.003883,0.00519901,0.00738286,0.00932063,0.0245266};

	// ---------------- Particle species: PiM ----------------
	VEC GeV14_pt_bin_center_PiM_0_80 = {0.295221,0.490592,0.688681,0.887672,1.08711,1.28685,1.48689,1.68722,1.88783,2.15616,2.59528};
	VEC GeV14_v2_values_PiM_0_80 = {0.0234208,0.0404856,0.0571415,0.0715853,0.0844605,0.0945863,0.101821,0.110021,0.11202,0.118464,0.121919};
	VEC GeV14_v2_stat_error_PiM_0_80 = {0.0001733,0.000202682,0.00028585,0.000424242,0.000660853,0.000926337,0.00135359,0.00198582,0.00291043,0.00348546,0.00687431};

	// ---------------- Particle species: PiP ----------------
	VEC GeV14_pt_bin_center_PiP_0_80 = {0.295297,0.490673,0.688767,0.88776,1.0872,1.28694,1.48696,1.68726,1.88782,2.15582,2.59647};
	VEC GeV14_v2_values_PiP_0_80 = {0.0211893,0.0390391,0.055217,0.0700864,0.0813231,0.092405,0.0967576,0.104074,0.109229,0.118892,0.120911};
	VEC GeV14_v2_stat_error_PiP_0_80 = {0.000178109,0.000205989,0.000287905,0.000425296,0.000660607,0.000929358,0.00136908,0.0020441,0.00306293,0.00369765,0.00769421};

	// ---------------- Particle species: KM ----------------
	VEC GeV14_pt_bin_center_KM_0_80 = {0.332933,0.505902,0.695388,0.890831,1.08869,1.28803,1.48851,1.68964,1.89064,2.16203,2.55108};
	VEC GeV14_v2_values_KM_0_80 = {0.0136032,0.0259191,0.0406916,0.0568871,0.0684453,0.0807841,0.091594,0.0993169,0.0977351,0.109567,0.126268};
	VEC GeV14_v2_stat_error_KM_0_80 = {0.00144301,0.000843419,0.000859905,0.00105577,0.0014177,0.00195579,0.00281425,0.00446621,0.00741984,0.0108014,0.0261304};

	// ---------------- Particle species: KP ----------------
	VEC GeV14_pt_bin_center_KP_0_80 = {0.332002,0.506593,0.69595,0.89134,1.08918,1.2885,1.48894,1.69004,1.89107,2.16457,2.56012};
	VEC GeV14_v2_values_KP_0_80 = {0.0148509,0.0280726,0.0440888,0.059656,0.0714214,0.0814646,0.0878935,0.0937463,0.0977549,0.106965,0.10216};
	VEC GeV14_v2_stat_error_KP_0_80 = {0.00107337,0.000618894,0.000620938,0.000749346,0.000993313,0.00134499,0.00184917,0.00269937,0.00415143,0.00715567,0.0252315};
	 
	// ---------------- Particle species: Proton ----------------
	VEC GeV14_pt_bin_center_Pr_0_10 = {0.350396,0.516707,0.704333,0.897901,1.09425,1.29207,1.49078,1.69011,1.88986,2.15985,2.58885};
	VEC GeV14_v2_values_Pr_0_10 = {-0.00465817,0.00761948,0.0143358,0.0229442,0.0369915,0.0448488,0.0556099,0.060298,0.07062,0.0725964,0.0801806};
	VEC GeV14_v2_stat_error_Pr_0_10 = {0.00378326,0.00124011,0.000983628,0.000976396,0.00109529,0.00131175,0.00165241,0.00216851,0.00299976,0.00368363,0.00956433};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV14_pt_bin_center_Ant_0_10 = {0.352354,0.517642,0.705114,0.898436,1.09448,1.29206,1.49062,1.68986,1.88961,2.1592,2.5867};
	VEC GeV14_v2_values_Ant_0_10 = {-0.0370647,-0.00559544,-0.00464836,0.00455593,0.0162261,0.0229623,0.0350269,0.0491747,0.0468177,0.0487421,0.00720832};
	VEC GeV14_v2_stat_error_Ant_0_10 = {0.0185394,0.00600283,0.00459227,0.00447922,0.00497676,0.00594595,0.00750407,0.00990898,0.0139721,0.0174732,0.0451476};

	// ---------------- Particle species: PiM ----------------
	VEC GeV14_pt_bin_center_PiM_0_10 = {0.298096,0.491153,0.687863,0.886051,1.08517,1.28521,1.48629,1.68827,1.89027,2.16544,2.58299};
	VEC GeV14_v2_values_PiM_0_10 = {0.0137955,0.0225186,0.0312685,0.0381084,0.0439899,0.0501234,0.0523408,0.0618582,0.0526906,0.0657157,0.0447023};
	VEC GeV14_v2_stat_error_PiM_0_10 = {0.00034887,0.000405204,0.000569817,0.000841573,0.0013034,0.00182778,0.00267542,0.00392206,0.00576864,0.00693652,0.0138484};

	// ---------------- Particle species: PiP ----------------
	VEC GeV14_pt_bin_center_PiP_0_10 = {0.298824,0.49155,0.688099,0.886197,1.08527,1.28525,1.4863,1.68825,1.89024,2.16527,2.58223};
	VEC GeV14_v2_values_PiP_0_10 = {0.0122915,0.0214318,0.0298809,0.038332,0.0424861,0.0487311,0.0481167,0.0528226,0.0516154,0.0549366,0.0502592};
	VEC GeV14_v2_stat_error_PiP_0_10 = {0.000359379,0.000411412,0.000571971,0.000840324,0.00129701,0.00182279,0.00268927,0.00402179,0.00602968,0.00733228,0.0155831};

	// ---------------- Particle species: KM ----------------
	VEC GeV14_pt_bin_center_KM_0_10 = {0.334965,0.506424,0.695909,0.891279,1.08898,1.28807,1.48827,1.68932,1.89064,2.16612,2.58389};
	VEC GeV14_v2_values_KM_0_10 = {0.00692282,0.0141286,0.0222386,0.0308528,0.0321586,0.0424212,0.0523672,0.0583173,0.0511273,0.048529,0.0935007};
	VEC GeV14_v2_stat_error_KM_0_10 = {0.00295206,0.00169558,0.00170577,0.00207347,0.00275989,0.00376054,0.00531266,0.00838491,0.0138107,0.0199117,0.0467237};

	// ---------------- Particle species: KP ----------------
	VEC GeV14_pt_bin_center_KP_0_10 = {0.335711,0.506913,0.696396,0.891737,1.08938,1.28838,1.48848,1.68951,1.89099,2.16958,2.60279};
	VEC GeV14_v2_values_KP_0_10 = {0.00952656,0.0145369,0.0236292,0.0325668,0.0393348,0.0505842,0.043207,0.0484112,0.0431522,0.057667,0.0229803};
	VEC GeV14_v2_stat_error_KP_0_10 = {0.00218677,0.00123798,0.00122619,0.00146688,0.00193236,0.00260169,0.00355389,0.00512107,0.00794758,0.0136407,0.0457996};

	VEC GeV14_pt_bin_center_Pr_10_40 = {0.34712,0.515889,0.702626,0.896969,1.09389,1.29198,1.49071,1.68985,1.88927,2.156,2.59464};
	VEC GeV14_v2_values_Pr_10_40 = {0.0076206,0.0203092,0.0386994,0.0609024,0.0854371,0.105293,0.123931,0.142986,0.155597,0.169599,0.185148};
	VEC GeV14_v2_stat_error_Pr_10_40 = {0.00181144,0.00060408,0.000491482,0.000500421,0.000573196,0.000697412,0.000886182,0.00116934,0.00162546,0.00200276,0.00520674};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV14_pt_bin_center_Ant_10_40 = {0.350397,0.516384,0.703786,0.897171,1.09338,1.29117,1.48998,1.68953,1.88961,2.16016,2.57842};
	VEC GeV14_v2_values_Ant_10_40 = {-0.00193207,0.0112685,0.0201522,0.0406254,0.0622747,0.0849408,0.0945954,0.108404,0.129379,0.139877,0.222264};
	VEC GeV14_v2_stat_error_Ant_10_40 = {0.00804927,0.00266784,0.00211773,0.00214997,0.00247358,0.0030351,0.00391108,0.00526031,0.00747988,0.00944059,0.0250765};

	// ---------------- Particle species: PiM ----------------
	VEC GeV14_pt_bin_center_PiM_10_40 = {0.297858,0.490942,0.687667,0.885859,1.08498,1.28502,1.48613,1.68816,1.89021,2.16528,2.5805};
	VEC GeV14_v2_values_PiM_10_40 = {0.0274649,0.0474804,0.0664517,0.0838382,0.0994407,0.109671,0.119446,0.127751,0.131659,0.137716,0.139558};
	VEC GeV14_v2_stat_error_PiM_10_40 = {0.000173267,0.00020288,0.000286251,0.000425194,0.000663382,0.000929292,0.00135624,0.00198716,0.00292075,0.00349173,0.00686911};

	// ---------------- Particle species: PiP ----------------
	VEC GeV14_pt_bin_center_PiP_10_40 = {0.298378,0.491226,0.687836,0.885967,1.08506,1.28508,1.48618,1.68821,1.89028,2.16563,2.58454};
	VEC GeV14_v2_values_PiP_10_40 = {0.0248813,0.0454857,0.0643,0.0813222,0.0959477,0.10721,0.114072,0.122654,0.13157,0.139833,0.157178};
	VEC GeV14_v2_stat_error_PiP_10_40 = {0.000177992,0.000206327,0.000288762,0.000427069,0.00066481,0.000934164,0.00137561,0.00204917,0.00302843,0.00364552,0.00782997};

	// ---------------- Particle species: KM ----------------
	VEC GeV14_pt_bin_center_KM_10_40 = {0.334361,0.505751,0.695383,0.890839,1.08865,1.28791,1.48837,1.68968,1.89108,2.16708,2.57463};
	VEC GeV14_v2_values_KM_10_40 = {0.0151304,0.0296265,0.0493236,0.066931,0.0819772,0.099713,0.109084,0.119107,0.126775,0.153734,0.147659};
	VEC GeV14_v2_stat_error_KM_10_40 = {0.00145661,0.000855301,0.000875341,0.00107788,0.00145085,0.00199688,0.00286309,0.00454215,0.00750365,0.0110345,0.0274931};

	// ---------------- Particle species: KP ----------------
	VEC GeV14_pt_bin_center_KP_10_40 = {0.333991,0.506346,0.695941,0.891361,1.08912,1.2883,1.48862,1.68979,1.89116,2.16813,2.58966};
	VEC GeV14_v2_values_KP_10_40 = {0.0162816,0.0334279,0.0533913,0.0709477,0.0854284,0.0962108,0.110287,0.117815,0.124677,0.129557,0.120061};
	VEC GeV14_v2_stat_error_KP_10_40 = {0.00108487,0.000629157,0.000634167,0.00076726,0.00101837,0.00137766,0.0018938,0.00275254,0.00423495,0.00741217,0.0257457};

	 // ---------------- Particle species: Proton ----------------
	VEC GeV14_pt_bin_center_Pr_40_80 = {0.347403,0.513985,0.700866,0.895297,1.09229,1.29047,1.48934,1.68871,1.88848,2.15576,2.60051};
	VEC GeV14_v2_values_Pr_40_80 = {0.0102291,0.0321039,0.0535407,0.0825052,0.112197,0.138716,0.161855,0.191263,0.207149,0.241324,0.220226};
	VEC GeV14_v2_stat_error_Pr_40_80 = {0.00599555,0.00207735,0.00178225,0.00192831,0.00231518,0.0029291,0.00383772,0.00519572,0.00737527,0.0092914,0.0242445};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV14_pt_bin_center_Ant_40_80 = {0.343809,0.512826,0.701221,0.894898,1.09119,1.28909,1.48818,1.68829,1.88911,2.16102,2.56645};
	VEC GeV14_v2_values_Ant_40_80 = {0.0252031,0.0207335,0.0426136,0.0686127,0.117497,0.120007,0.141241,0.201009,0.271719,0.207139,0.327501};
	VEC GeV14_v2_stat_error_Ant_40_80 = {0.0229306,0.00790462,0.00679639,0.00743444,0.00920554,0.0119811,0.0162085,0.0227628,0.033658,0.0443496,0.124101};

	// ---------------- Particle species: PiM ----------------
	VEC GeV14_pt_bin_center_PiM_40_80 = {0.297427,0.490246,0.686845,0.884965,1.08405,1.28413,1.48544,1.68784,1.89024,2.16575,2.57876};
	VEC GeV14_v2_values_PiM_40_80 = {0.0298475,0.0549447,0.0811457,0.103452,0.125139,0.144176,0.153663,0.160431,0.183169,0.169217,0.250727};
	VEC GeV14_v2_stat_error_PiM_40_80 = {0.000618149,0.000737065,0.0010604,0.00161567,0.00258533,0.00364416,0.00538727,0.0079768,0.0116669,0.0142668,0.0285448};

	// ---------------- Particle species: PiP ----------------
	VEC GeV14_pt_bin_center_PiP_40_80 = {0.297605,0.490369,0.686939,0.885035,1.08409,1.28411,1.48534,1.68768,1.89009,2.16559,2.5811};
	VEC GeV14_v2_values_PiP_40_80 = {0.0271215,0.0546943,0.0796372,0.103092,0.121423,0.144971,0.15315,0.163211,0.164864,0.20047,0.145828};
	VEC GeV14_v2_stat_error_PiP_40_80 = {0.000631853,0.000750223,0.00107438,0.00163335,0.00260818,0.00370183,0.00554412,0.0083759,0.0125396,0.015409,0.0327468};

	// ---------------- Particle species: KM ----------------
	VEC GeV14_pt_bin_center_KM_40_80 = {0.332206,0.504304,0.693879,0.889493,1.08765,1.28744,1.48835,1.68952,1.8899,2.15087,2.55914};
	VEC GeV14_v2_values_KM_40_80 = {0.022242,0.0400112,0.0546369,0.0892272,0.121609,0.114028,0.144942,0.153222,0.138086,0.121096,0.258792};
	VEC GeV14_v2_stat_error_KM_40_80 = {0.00513029,0.00314264,0.00335056,0.00428235,0.0059647,0.00859701,0.0134231,0.0235027,0.0398029,0.0535074,0.169258};

	// ---------------- Particle species: KP ----------------
	VEC GeV14_pt_bin_center_KP_40_80 = {0.332642,0.504919,0.694676,0.890262,1.08829,1.28789,1.48873,1.6902,1.89137,2.16509,2.55713};
	VEC GeV14_v2_values_KP_40_80 = {0.0216217,0.0412175,0.0623313,0.0913743,0.109836,0.112726,0.130896,0.130032,0.167509,0.176188,0.228367};
	VEC GeV14_v2_stat_error_KP_40_80 = {0.00390405,0.00235665,0.00246583,0.00307943,0.00420415,0.00582889,0.00828009,0.0125238,0.0196264,0.0311518,0.084073};
	
	VEC GeV14_v2_0_5 = {0.0112853,0.0169279,0.0225705,0.0272727,0.0319749,0.0347962,0.0394984,0.0442006,0.0460815,0.0489028,0.053605,0.0442006,0.0460815,0.0630094};
	VEC GeV14_v2_5_10 = {0.0159875,0.0263323,0.0357367,0.0423197,0.0507837,0.0583072,0.0658307,0.0761755,0.0799373,0.0846395,0.0855799,0.0931034,0.092163,0.101567};
	VEC GeV14_v2_10_20 = {0.0225705,0.0376176,0.0517241,0.0639498,0.0752351,0.0874608,0.100627,0.110031,0.121317,0.129781,0.133542,0.140125,0.141066,0.14953};
	VEC GeV14_v2_20_30 = {0.0291536,0.0489028,0.0667712,0.0836991,0.0996865,0.116614,0.129781,0.145768,0.157053,0.165517,0.175862,0.187147,0.19185,0.204075}; 
	VEC GeV14_v2_30_40 = {0.030094,0.053605,0.0742947,0.0940439,0.113793,0.130721,0.146708,0.160815,0.175862,0.187147,0.195611,0.207837,0.210658,0.228527};
	VEC GeV14_v2_40_50 = {0.030094,0.053605,0.0761755,0.0968652,0.116614,0.134483,0.148589,0.164577,0.176803,0.185266,0.19185, 0.195611,0.205016,0.239812};  

    VEC GeV14_err={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

	if(Energy==14){
		
		DrawArticle( Energy, GeV14_pt_bin_center_Ant_0_80, GeV14_v2_values_Ant_0_80, GeV14_v2_stat_error_Ant_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Neg", 0.33,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_Pr_0_80, GeV14_v2_values_Pr_0_80, GeV14_v2_stat_error_Pr_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_KP_0_80, GeV14_v2_values_KP_0_80, GeV14_v2_stat_error_KP_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Pos", 0.28,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_KM_0_80, GeV14_v2_values_KM_0_80, GeV14_v2_stat_error_KM_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Neg", 0.28,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_PiP_0_80, GeV14_v2_values_PiP_0_80, GeV14_v2_stat_error_PiP_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Pos", 0.28,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_PiM_0_80, GeV14_v2_values_PiM_0_80, GeV14_v2_stat_error_PiM_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Neg", 0.28,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_Ant_0_10, GeV14_v2_values_Ant_0_10, GeV14_v2_stat_error_Ant_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Neg", 0.22,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_Pr_0_10, GeV14_v2_values_Pr_0_10, GeV14_v2_stat_error_Pr_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Pos", 0.22,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_KP_0_10, GeV14_v2_values_KP_0_10, GeV14_v2_stat_error_KP_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Pos", 0.16,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_KM_0_10, GeV14_v2_values_KM_0_10, GeV14_v2_stat_error_KM_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Neg", 0.16,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_PiP_0_10, GeV14_v2_values_PiP_0_10, GeV14_v2_stat_error_PiP_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Pos", 0.14,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_PiM_0_10, GeV14_v2_values_PiM_0_10, GeV14_v2_stat_error_PiM_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Neg", 0.14,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_Pr_10_40, GeV14_v2_values_Pr_10_40, GeV14_v2_stat_error_Pr_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_Ant_10_40, GeV14_v2_values_Ant_10_40, GeV14_v2_stat_error_Ant_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Neg", 0.33,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_KP_10_40, GeV14_v2_values_KP_10_40, GeV14_v2_stat_error_KP_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Pos", 0.28,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_KM_10_40, GeV14_v2_values_KM_10_40, GeV14_v2_stat_error_KM_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Neg", 0.28,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_PiP_10_40, GeV14_v2_values_PiP_10_40, GeV14_v2_stat_error_PiP_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Pos", 0.25,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_PiM_10_40, GeV14_v2_values_PiM_10_40, GeV14_v2_stat_error_PiM_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Neg", 0.25,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_Pr_40_80, GeV14_v2_values_Pr_40_80, GeV14_v2_stat_error_Pr_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Pos", 0.42,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_Ant_40_80, GeV14_v2_values_Ant_40_80, GeV14_v2_stat_error_Ant_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Neg", 0.42,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_KP_40_80, GeV14_v2_values_KP_40_80, GeV14_v2_stat_error_KP_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_KM_40_80, GeV14_v2_values_KM_40_80, GeV14_v2_stat_error_KM_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Neg", 0.33,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_PiP_40_80, GeV14_v2_values_PiP_40_80, GeV14_v2_stat_error_PiP_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV14_pt_bin_center_PiM_40_80, GeV14_v2_values_PiM_40_80, GeV14_v2_stat_error_PiM_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Neg", 0.33,3.0);
	
		DrawArticle( Energy, ptHadrons, GeV14_v2_0_5, GeV14_err, GeV14_err, 8, 8, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, ptHadrons, GeV14_v2_5_10, GeV14_err, GeV14_err, 7, 7, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, ptHadrons, GeV14_v2_10_20, GeV14_err, GeV14_err, 6, 6, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, ptHadrons, GeV14_v2_20_30, GeV14_err, GeV14_err, 5, 5, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, ptHadrons, GeV14_v2_30_40, GeV14_err, GeV14_err, 4, 4, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, ptHadrons, GeV14_v2_40_50, GeV14_err, GeV14_err, 3, 3, "Eta15","Hadrons", "", 0.33,3.6);

	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////// beam energy: 19 GeV /////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Event plane method: eta-sub
	// ---------------- Particle species: Proton ----------------
	VEC GeV19_pt_bin_center_Pr_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV19_v2_values_Pr_0_80 = {0.00600749,0.0168896,0.0325285,0.052604,0.0717285,0.090242,0.106005,0.119548,0.131381,0.140943,0.148884,0.152717,0.167373,0.167586,0.174867};
	VEC GeV19_v2_stat_error_Pr_0_80 = {0.000948076,0.000369558,0.000315547,0.000322681,0.000366482,0.000444058,0.000562391,0.000741484,0.00101034,0.00128781,0.00175503,0.00246738,0.00351355,0.00412564,0.00863469};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV19_pt_bin_center_Ant_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.835};
	VEC GeV19_v2_values_Ant_0_80 = {0.00296996,0.00738706,0.0221654,0.0385976,0.0591853,0.074907,0.0919459,0.106479,0.116709,0.12735,0.134329,0.14933,0.149541};
	VEC GeV19_v2_stat_error_Ant_0_80 = {0.00280851,0.00117093,0.000988283,0.00101483,0.00116976,0.00145031,0.00188025,0.00256276,0.00312193,0.00431079,0.00595189,0.00859351,0.0103928};

	// ---------------- Particle species: PiM ----------------
	VEC GeV19_pt_bin_center_PiM_0_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.805,3.435};
	VEC GeV19_v2_values_PiM_0_80 = {0.0216746,0.0408128,0.0586074,0.0746435,0.0876909,0.0991869,0.107375,0.112568,0.117449,0.126703,0.127578,0.134613,0.143151};
	VEC GeV19_v2_stat_error_PiM_0_80 = {9.84287e-05,0.000120998,0.000168542,0.000243016,0.000353245,0.000512965,0.000749998,0.00109194,0.00155004,0.00217798,0.00249147,0.00449093,0.012078};

	// ---------------- Particle species: PiP ----------------
	VEC GeV19_pt_bin_center_PiP_0_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.805,3.465};
	VEC GeV19_v2_values_PiP_0_80 = {0.0200471,0.0393043,0.0571173,0.0723951,0.0856562,0.0955233,0.10471,0.111428,0.117715,0.119919,0.12482,0.12953,0.147494};
	VEC GeV19_v2_stat_error_PiP_0_80 = {9.98816e-05,0.000121897,0.000168728,0.00024214,0.00035,0.000505203,0.000737379,0.00107473,0.0015368,0.00215821,0.00246628,0.00440596,0.0115101};

	// ---------------- Particle species: KM ----------------
	VEC GeV19_pt_bin_center_KM_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.435};
	VEC GeV19_v2_values_KM_0_80 = {0.0125315,0.0266821,0.0434517,0.0591548,0.0720669,0.0829478,0.0938115,0.102492,0.107511,0.103098,0.114785,0.127232,0.123313};
	VEC GeV19_v2_stat_error_KM_0_80 = {0.000704145,0.000456717,0.000472272,0.000571875,0.00073736,0.000992619,0.00137958,0.00172655,0.00288389,0.00429175,0.00526242,0.0102767,0.0278666};

	// ---------------- Particle species: KP ----------------
	VEC GeV19_pt_bin_center_KP_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.435};
	VEC GeV19_v2_values_KP_0_80 = {0.0138283,0.0276401,0.0458078,0.0606341,0.0737375,0.0846323,0.0928695,0.0976066,0.109287,0.109269,0.113533,0.121311,0.118028};
	VEC GeV19_v2_stat_error_KP_0_80 = {0.000564368,0.000361223,0.000368676,0.000440327,0.000562521,0.000748538,0.00102465,0.0014212,0.00198029,0.00284172,0.00344385,0.00725556,0.0215644};
	// ---------------- Particle species: Proton ----------------
	VEC GeV19_pt_bin_center_Pr_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.535,3.045};
	VEC GeV19_v2_values_Pr_0_10 = {-0.00269425,0.00596139,0.0139036,0.0256201,0.0368745,0.0486643,0.0565901,0.0641655,0.0711467,0.0739312,0.079297,0.0861715,0.0865048};
	VEC GeV19_v2_stat_error_Pr_0_10 = {0.00189345,0.000728041,0.000609799,0.000611877,0.000684099,0.000819942,0.00103074,0.00135267,0.00183912,0.00234691,0.00320095,0.00375084,0.00726039};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV19_pt_bin_center_Ant_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.755,2.145,2.625};
	VEC GeV19_v2_values_Ant_0_10 = {-0.00237616,-0.00616014,0.00388062,0.0148079,0.0228002,0.0351128,0.0369546,0.0566465,0.0658163,0.0731895};
	VEC GeV19_v2_stat_error_Ant_0_10 = {0.00595461,0.00242394,0.00198431,0.00197878,0.00223006,0.00271744,0.00347042,0.00335779,0.00624111,0.0120133};

	// ---------------- Particle species: PiM ----------------
	VEC GeV19_pt_bin_center_PiM_0_10 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.965,2.355,2.805,3.435};
	VEC GeV19_v2_values_PiM_0_10 = {0.0132513,0.0236795,0.0330047,0.0409265,0.0488099,0.0544218,0.0566759,0.0560667,0.0570697,0.0610386,0.0665234,0.0265244};
	VEC GeV19_v2_stat_error_PiM_0_10 = {0.000188222,0.000229068,0.000318159,0.000457252,0.000662863,0.000962666,0.00140931,0.00205714,0.00239027,0.00470809,0.00834468,0.0238626};

	// ---------------- Particle species: PiP ----------------
	VEC GeV19_pt_bin_center_PiP_0_10 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.965,2.355,2.805,3.435};
	VEC GeV19_v2_values_PiP_0_10 = {0.0113998,0.0222471,0.0315409,0.0396261,0.0451527,0.0506575,0.0559122,0.0580258,0.0610479,0.0636403,0.061057,0.0720694};
	VEC GeV19_v2_stat_error_PiP_0_10 = {0.000191388,0.0002306,0.000317637,0.000453747,0.00065397,0.000942212,0.00137489,0.00200904,0.00236539,0.00460625,0.00825042,0.0198218};

	// ---------------- Particle species: KM ----------------
	VEC GeV19_pt_bin_center_KM_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.965,2.355,2.805,3.465};
	VEC GeV19_v2_values_KM_0_10 = {0.00731444,0.0158128,0.0243354,0.0334721,0.0403452,0.0460981,0.0525503,0.059217,0.0615636,0.0540598,0.0803659,0.0667919};
	VEC GeV19_v2_stat_error_KM_0_10 = {0.00137125,0.000874858,0.000892725,0.00107018,0.00137016,0.0018355,0.00254587,0.00359401,0.00429181,0.00937257,0.0200771,0.035088};

	// ---------------- Particle species: KP ----------------
	VEC GeV19_pt_bin_center_KP_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.965,2.355,2.805,3.465};
	VEC GeV19_v2_values_KP_0_10 = {0.00601667,0.0141115,0.0254069,0.0341006,0.0422171,0.0472486,0.05151,0.0495944,0.0607535,0.0603991,0.0409356,-0.0165935};
	VEC GeV19_v2_stat_error_KP_0_10 = {0.0010949,0.000688599,0.000694159,0.00082177,0.00104449,0.00138536,0.00189543,0.00261313,0.00298683,0.00644024,0.0137783,0.0306178};
	// ---------------- Particle species: Proton ----------------
	VEC GeV19_pt_bin_center_Pr_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV19_v2_values_Pr_10_40 = {0.0079273,0.0184922,0.0375052,0.0610131,0.0846079,0.107231,0.127185,0.144573,0.157226,0.172421,0.18198,0.188198,0.196454,0.199378,0.220326};
	VEC GeV19_v2_stat_error_Pr_10_40 = {0.000958179,0.000375148,0.000322329,0.000331998,0.000379322,0.000461302,0.000585403,0.000772619,0.00105291,0.00134161,0.00182673,0.00256816,0.00365334,0.00427941,0.00895054};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV19_pt_bin_center_Ant_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.535,3.015};
	VEC GeV19_v2_values_Ant_10_40 = {0.00378815,0.00628343,0.0240837,0.0435258,0.0662957,0.0844981,0.108297,0.125534,0.135189,0.147956,0.158745,0.172094,0.177293};
	VEC GeV19_v2_stat_error_Ant_10_40 = {0.00275103,0.00114717,0.000973641,0.00100899,0.00117106,0.00145787,0.00189569,0.00259915,0.00317992,0.004402,0.0060956,0.00742357,0.0152719};

	// ---------------- Particle species: PiM ----------------
	VEC GeV19_pt_bin_center_PiM_10_40 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.805,3.435};
	VEC GeV19_v2_values_PiM_10_40 = {0.0252193,0.0472496,0.0679023,0.0865443,0.101499,0.115212,0.125521,0.132089,0.139212,0.147458,0.148416,0.155312,0.175952};
	VEC GeV19_v2_stat_error_PiM_10_40 = {9.9107e-05,0.000122135,0.000170355,0.000246219,0.000357571,0.0005191,0.000757798,0.0011006,0.00156074,0.00218957,0.0025013,0.00448384,0.0120854};

	// ---------------- Particle species: PiP ----------------
	VEC GeV19_pt_bin_center_PiP_10_40 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.805,3.435};
	VEC GeV19_v2_values_PiP_10_40 = {0.023176,0.0455854,0.0663291,0.0844592,0.10061,0.112404,0.122129,0.130671,0.135372,0.139406,0.145359,0.155191,0.167817};
	VEC GeV19_v2_stat_error_PiP_10_40 = {0.000100528,0.00012314,0.000170853,0.000245534,0.000355214,0.000512996,0.000748384,0.00108813,0.00155717,0.00218161,0.00248504,0.00444503,0.0113176};

	// ---------------- Particle species: KM ----------------
	VEC GeV19_pt_bin_center_KM_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.435};
	VEC GeV19_v2_values_KM_10_40 = {0.0131962,0.0302134,0.0506047,0.0689414,0.0857591,0.0997649,0.111843,0.121981,0.126466,0.134673,0.140281,0.139083,0.243525};
	VEC GeV19_v2_stat_error_KM_10_40 = {0.000717721,0.000468306,0.000486654,0.000591505,0.000764339,0.00102893,0.0014294,0.00179417,0.00297409,0.00444372,0.0054256,0.0102166,0.0295786};

	// ---------------- Particle species: KP ----------------
	VEC GeV19_pt_bin_center_KP_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.465};
	VEC GeV19_v2_values_KP_10_40 = {0.0159419,0.032311,0.0541594,0.0721828,0.0884912,0.101852,0.112669,0.121144,0.130436,0.136791,0.13644,0.148666,0.155278};
	VEC GeV19_v2_stat_error_KP_10_40 = {0.000577535,0.000371908,0.000381561,0.000457142,0.0005847,0.000777754,0.00106251,0.0014733,0.00205421,0.00295317,0.00355388,0.00736402,0.0142418};
	// ---------------- Particle species: Proton ----------------
	VEC GeV19_pt_bin_center_Pr_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.535,3.015};
	VEC GeV19_v2_values_Pr_40_80 = {0.0146792,0.0322517,0.0561763,0.0894737,0.119115,0.147322,0.176066,0.196203,0.226235,0.232669,0.241099,0.259582,0.303711};
	VEC GeV19_v2_stat_error_Pr_40_80 = {0.0032456,0.00132755,0.00121138,0.00132148,0.00158945,0.00201447,0.00264225,0.00358479,0.00497677,0.0063807,0.00878495,0.0103135,0.0202461};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV19_pt_bin_center_Ant_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.665,1.875,2.145,2.625};
	VEC GeV19_v2_values_Ant_40_80 = {0.00651495,0.0249051,0.04333,0.0649464,0.10943,0.132852,0.166042,0.177224,0.199777,0.241737,0.295131};
	VEC GeV19_v2_stat_error_Ant_40_80 = {0.00811936,0.00356913,0.00326338,0.00364186,0.00453056,0.00595785,0.00814017,0.0116873,0.0148124,0.0170495,0.03512};

	// ---------------- Particle species: PiM ----------------
	VEC GeV19_pt_bin_center_PiM_40_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.965,2.355,2.805};
	VEC GeV19_v2_values_PiM_40_80 = {0.0273798,0.0559095,0.0831203,0.10911,0.128273,0.146238,0.160045,0.173936,0.189074,0.19978,0.172843};
	VEC GeV19_v2_stat_error_PiM_40_80 = {0.000367276,0.000464134,0.000660034,0.000970311,0.00143671,0.00211566,0.00312741,0.00456344,0.00533244,0.0106271,0.0194378};

	// ---------------- Particle species: PiP ----------------
	VEC GeV19_pt_bin_center_PiP_40_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.965,2.355,2.805};
	VEC GeV19_v2_values_PiP_40_80 = {0.0274083,0.0550259,0.0826116,0.105757,0.128192,0.142751,0.160306,0.170704,0.188821,0.202585,0.20945};
	VEC GeV19_v2_stat_error_PiP_40_80 = {0.000370986,0.000468384,0.000664699,0.000975501,0.00143788,0.0021199,0.00312621,0.004565,0.00539526,0.0105904,0.0192625};

	// ---------------- Particle species: KM ----------------
	VEC GeV19_pt_bin_center_KM_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.665,1.965,2.325,2.805};
	VEC GeV19_v2_values_KM_40_80 = {0.0210452,0.0396515,0.0656698,0.0920938,0.108824,0.123236,0.146519,0.154115,0.15688,0.228979,0.222329};
	VEC GeV19_v2_stat_error_KM_40_80 = {0.00261865,0.00178654,0.00193436,0.0024342,0.00322794,0.00444462,0.00626047,0.00912976,0.010568,0.0217863,0.0256622};

	// ---------------- Particle species: KP ----------------
	VEC GeV19_pt_bin_center_KP_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.965,2.355,2.805};
	VEC GeV19_v2_values_KP_40_80 = {0.0236532,0.0440603,0.0684008,0.0901399,0.106535,0.12711,0.137241,0.150715,0.168326,0.180259,0.221354};
	VEC GeV19_v2_stat_error_KP_40_80 = {0.00214311,0.00144386,0.0015403,0.00190202,0.00248372,0.00337005,0.00466216,0.00658947,0.0074734,0.0158791,0.0307244};
	
	VEC v2_cent_19GeV = {0.01895, 0.02904, 0.04138, 0.05304, 0.05857, 0.05916, 0.05525, 0.05042, 0.04727};
	VEC v2_cent_19GeV_err = {0.00013, 0.00011, 0.00007, 0.00008, 0.00011, 0.00017, 0.00030, 0.00059, 0.00120};

	VEC GeV19_v2_10_20 = {0.0208,0.0377,0.0526,0.0664,0.0799,0.0927,0.1046,0.1141,0.1244,0.1334,0.1379,0.1420,0.1477,0.1575};
	VEC GeV19_v2_10_20_err = {0.0001,0.0001,0.0002,0.0002,0.0003,0.0004,0.0006,0.0008,0.0010,0.0014,0.0019,0.0028,0.0038,0.0043};

	VEC GeV19_v2_20_30 = {0.0268,0.0487,0.0684,0.0862,0.1035,0.1202,0.1370,0.1498,0.1622,0.1705,0.1821,0.1870,0.1963,0.2003};
	VEC GeV19_v2_20_30_err = {0.0001,0.0002,0.0002,0.0003,0.0004,0.0005,0.0007,0.0009,0.0012,0.0017,0.0023,0.0033,0.0045,0.0051};

	VEC GeV19_v2_30_40 = {0.0294,0.0544,0.0772,0.0981,0.1183,0.1357,0.1531,0.1687,0.1789,0.1921,0.2065,0.2036,0.2229,0.2289};
	VEC GeV19_v2_30_40_err = {0.0002,0.0002,0.0003,0.0004,0.0005,0.0007,0.0009,0.0013,0.0017,0.0024,0.0033,0.0048,0.0066,0.0074};

	VEC GeV19_v2_50_60 = {0.0272,0.0530,0.0779,0.0991,0.1201,0.1359,0.1518,0.1665,0.1870,0.1900,0.1952,0.2054,0.2176,0.2026};
	VEC GeV19_v2_50_60_err = {0.0004,0.0006,0.0008,0.0011,0.0015,0.0021,0.0029,0.0040,0.0056,0.0079,0.0110,0.0162,0.0224,0.0255};

	if(Energy==19){
		DrawArticle( Energy, GeV19_pt_bin_center_Ant_0_80, GeV19_v2_values_Ant_0_80, GeV19_v2_stat_error_Ant_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Neg", 0.33,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_Pr_0_80, GeV19_v2_values_Pr_0_80, GeV19_v2_stat_error_Pr_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_KP_0_80, GeV19_v2_values_KP_0_80, GeV19_v2_stat_error_KP_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Pos", 0.28,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_KM_0_80, GeV19_v2_values_KM_0_80, GeV19_v2_stat_error_KM_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Neg", 0.28,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_PiP_0_80, GeV19_v2_values_PiP_0_80, GeV19_v2_stat_error_PiP_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Pos", 0.28,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_PiM_0_80, GeV19_v2_values_PiM_0_80, GeV19_v2_stat_error_PiM_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Neg", 0.28,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_Ant_0_10, GeV19_v2_values_Ant_0_10, GeV19_v2_stat_error_Ant_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Neg", 0.22,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_Pr_0_10, GeV19_v2_values_Pr_0_10, GeV19_v2_stat_error_Pr_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Pos", 0.22,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_KP_0_10, GeV19_v2_values_KP_0_10, GeV19_v2_stat_error_KP_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Pos", 0.16,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_KM_0_10, GeV19_v2_values_KM_0_10, GeV19_v2_stat_error_KM_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Neg", 0.16,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_PiP_0_10, GeV19_v2_values_PiP_0_10, GeV19_v2_stat_error_PiP_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Pos", 0.14,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_PiM_0_10, GeV19_v2_values_PiM_0_10, GeV19_v2_stat_error_PiM_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Neg", 0.14,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_Pr_10_40, GeV19_v2_values_Pr_10_40, GeV19_v2_stat_error_Pr_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_Ant_10_40, GeV19_v2_values_Ant_10_40, GeV19_v2_stat_error_Ant_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Neg", 0.33,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_KP_10_40, GeV19_v2_values_KP_10_40, GeV19_v2_stat_error_KP_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Pos", 0.28,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_KM_10_40, GeV19_v2_values_KM_10_40, GeV19_v2_stat_error_KM_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Neg", 0.28,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_PiP_10_40, GeV19_v2_values_PiP_10_40, GeV19_v2_stat_error_PiP_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Pos", 0.25,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_PiM_10_40, GeV19_v2_values_PiM_10_40, GeV19_v2_stat_error_PiM_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Neg", 0.25,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_Pr_40_80, GeV19_v2_values_Pr_40_80, GeV19_v2_stat_error_Pr_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Pos", 0.42,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_Ant_40_80, GeV19_v2_values_Ant_40_80, GeV19_v2_stat_error_Ant_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Neg", 0.42,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_KP_40_80, GeV19_v2_values_KP_40_80, GeV19_v2_stat_error_KP_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_KM_40_80, GeV19_v2_values_KM_40_80, GeV19_v2_stat_error_KM_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Neg", 0.33,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_PiP_40_80, GeV19_v2_values_PiP_40_80, GeV19_v2_stat_error_PiP_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Pos", 0.33,3.0);
		DrawArticle( Energy, GeV19_pt_bin_center_PiM_40_80, GeV19_v2_values_PiM_40_80, GeV19_v2_stat_error_PiM_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Neg", 0.33,3.0);
	

		DrawArticle( Energy, ptHadrons, GeV19_v2_10_20, GeV19_v2_10_20_err, ptGeV_err, 6, 6, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, ptHadrons, GeV19_v2_20_30, GeV19_v2_20_30_err, ptGeV_err, 5, 5, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, ptHadrons, GeV19_v2_30_40, GeV19_v2_30_40_err, ptGeV_err, 4, 4, "Eta15","Hadrons", "", 0.33,3.6);
		DrawArticle( Energy, ptHadrons, GeV19_v2_50_60, GeV19_v2_50_60_err, ptGeV_err, 3, 3, "Eta15","Hadrons", "", 0.33,3.6);

	}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////// 62.4 GeV //////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	VEC pt62GeV  = {0.3, 0.5, 0.7, 0.9, 1.1, 1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0};
	VEC pt62GeV_err = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	VEC v2_cent_62GeV = {0.02089, 0.03102, 0.04407, 0.05658, 0.06318, 0.06436, 0.06122, 0.05585, 0.05284};
	VEC v2_cent_62GeV_err = {0.00004, 0.00003, 0.00002, 0.00003, 0.00003, 0.00005, 0.00009, 0.00018, 0.00036};

	VEC GeV62_v2_10_20  = {0.0211,0.0356,0.0483,0.0612,0.0736,0.0855,0.0986,0.1061,0.1178,0.1297,0.1351,0.1396,0.1463,0.1523};
	VEC GeV62_v2_10_20_err = {0.0002,0.0002,0.0003,0.0004,0.0005,0.0007,0.0010,0.0013,0.0018,0.0025,0.0035,0.0051,0.0070,0.0080};

	VEC GeV62_v2_20_30  = {0.0263,0.0462,0.0632,0.0796,0.0967,0.1135,0.1280,0.1436,0.1528,0.1669,0.1753,0.1842,0.1890,0.2095};//
	VEC GeV62_v2_20_30_err = {0.0002,0.0002,0.0003,0.0004,0.0005,0.0006,0.0008,0.0011,0.0016,0.0022,0.0030,0.0042,0.0059,0.0067};

	VEC GeV62_v2_30_40  = {0.0295,0.0514,0.0717,0.0908,0.1097,0.1273,0.1450,0.1585,0.1740,0.1810,0.1907,0.1962,0.2224,0.2169};
	VEC GeV62_v2_30_40_err = {0.0003,0.0004,0.0005,0.0007,0.0009,0.0012,0.0016,0.0022,0.0031,0.0044,0.0061,0.0091,0.0128,0.0148};

	VEC GeV62_v2_50_60  = {0.0287,0.0508,0.0694,0.0916,0.1090,0.1223,0.1338,0.1471,0.1819,0.1567,0.2149,0.2111,0.1518,0.2272};
	VEC GeV62_v2_50_60_err = {0.0008,0.0010,0.0014,0.0019,0.0027,0.0037,0.0052,0.0074,0.0105,0.0151,0.0220,0.0333,0.0486,0.0575};
	// ****************** beam Energy: 11.5 GeV ******************
	// Event plane method: eta-sub

	// ---------------- Particle species: Proton ----------------
	VEC GeV62_pt_bin_center_Pr_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV62_v2_values_Pr_0_80 = {0.00683153,0.0158907,0.0315793,0.0516096,0.0727302,0.0925249,0.110851,0.126435,0.140123,0.150104,0.160199,0.165746,0.172611,0.179436,0.184446};
	VEC GeV62_v2_stat_error_Pr_0_80 = {0.000520759,0.000195399,0.000163827,0.000164381,0.00018273,0.000216225,0.000267061,0.000342627,0.000453176,0.000557974,0.000736763,0.00100473,0.00138938,0.00155757,0.00303019};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV62_pt_bin_center_Ant_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV62_v2_values_Ant_0_80 = {0.00494108,0.0114479,0.0261197,0.0449316,0.0663671,0.0860504,0.103827,0.120229,0.132341,0.142772,0.151479,0.165265,0.16601,0.17546,0.178903};
	VEC GeV62_v2_stat_error_Ant_0_80 = {0.000897342,0.000371594,0.000306912,0.000307998,0.000345265,0.000415456,0.000522919,0.000689985,0.000821512,0.00110148,0.00148064,0.00206483,0.00291789,0.00336361,0.00684634};

	// ---------------- Particle species: PiM ----------------
	VEC GeV62_pt_bin_center_PiM_0_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV62_v2_values_PiM_0_80 = {0.0223896,0.042742,0.0624241,0.0798173,0.0944707,0.106976,0.11663,0.125184,0.131548,0.135949,0.138445,0.14286,0.148789,0.143203,0.151072};
	VEC GeV62_v2_stat_error_PiM_0_80 = {4.02884e-05,4.85174e-05,6.65846e-05,9.43872e-05,0.000134239,0.000189959,0.000271414,0.000385102,0.00053152,0.000717481,0.000968061,0.00129079,0.00172756,0.00185626,0.00332027};

	// ---------------- Particle species: PiP ----------------
	VEC GeV62_pt_bin_center_PiP_0_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV62_v2_values_PiP_0_80 = {0.0208012,0.041395,0.0611387,0.0783936,0.0932306,0.105414,0.115887,0.123129,0.129545,0.134079,0.136608,0.141293,0.145453,0.144127,0.147837};
	VEC GeV62_v2_stat_error_PiP_0_80 = {4.06899e-05,4.87531e-05,6.64991e-05,9.38308e-05,0.000132824,0.000187112,0.000265895,0.000375892,0.000515683,0.000692813,0.000918648,0.00122249,0.00162044,0.00172633,0.00303448};

	// ---------------- Particle species: KM ----------------
	VEC GeV62_pt_bin_center_KM_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.685,2.955,3.345};
	VEC GeV62_v2_values_KM_0_80 = {0.0129623,0.02816,0.0468495,0.0647472,0.0797767,0.0927023,0.102755,0.109907,0.117535,0.120576,0.121677,0.124428,0.126282,0.144784,0.139207};
	VEC GeV62_v2_stat_error_KM_0_80 = {0.000275854,0.000175113,0.00017688,0.000208595,0.000262253,0.000341829,0.000461195,0.000635014,0.000904939,0.00131148,0.00196627,0.00289826,0.00425974,0.0049056,0.0110844};

	// ---------------- Particle species: KP ----------------
	VEC GeV62_pt_bin_center_KP_0_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.685,2.955,3.345};
	VEC GeV62_v2_values_KP_0_80 = {0.0131131,0.0280332,0.0473805,0.0650401,0.0797504,0.0919515,0.101452,0.108697,0.114705,0.120665,0.121731,0.128768,0.123222,0.127912,0.134311};
	VEC GeV62_v2_stat_error_KP_0_80 = {0.000244437,0.000153924,0.000154119,0.000180181,0.000224963,0.000290565,0.000386755,0.000523833,0.000719747,0.0010333,0.00151358,0.002247,0.00334946,0.00396006,0.0063682};

	// ---------------- Particle species: Proton ----------------
	VEC GeV62_pt_bin_center_Pr_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,3.045};
	VEC GeV62_v2_values_Pr_0_10 = {0.00130562,0.00184556,0.00907538,0.0195195,0.0327654,0.0467786,0.0594504,0.0698091,0.0784107,0.0876365,0.0923305,0.100487,0.101097,0.108177};
	VEC GeV62_v2_stat_error_Pr_0_10 = {0.00177184,0.000604111,0.000480588,0.000465061,0.00050042,0.000579048,0.000699772,0.000887725,0.00117942,0.00140895,0.0018797,0.00264607,0.00379723,0.00443438};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV62_pt_bin_center_Ant_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.535,3.045};
	VEC GeV62_v2_values_Ant_0_10 = {0.00332361,-0.00271889,0.00474151,0.0168899,0.0287387,0.0409905,0.0521202,0.0625626,0.0725457,0.0767136,0.0813841,0.092351,0.0830351};
	VEC GeV62_v2_stat_error_Ant_0_10 = {0.0025198,0.000969667,0.000754741,0.000727005,0.000783946,0.000913447,0.00111004,0.00142714,0.00165818,0.00222136,0.00300224,0.00368081,0.00782532};

	// ---------------- Particle species: PiM ----------------
	VEC GeV62_pt_bin_center_PiM_0_10 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.805,3.495};
	VEC GeV62_v2_values_PiM_0_10 = {0.0100948,0.0213515,0.0333816,0.0435058,0.0517501,0.0593174,0.063837,0.0699047,0.0729183,0.0730076,0.0760939,0.0730191,0.0692482};
	VEC GeV62_v2_stat_error_PiM_0_10 = {0.000104981,0.000122213,0.000165392,0.00023036,0.000322747,0.000452209,0.000639048,0.000901591,0.00124789,0.00168059,0.00178641,0.00283539,0.00606821};

	// ---------------- Particle species: PiP ----------------
	VEC GeV62_pt_bin_center_PiP_0_10 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.355,2.835,3.495};
	VEC GeV62_v2_values_PiP_0_10 = {0.0105363,0.0219187,0.0341168,0.0437394,0.0526304,0.0597022,0.0662818,0.0699636,0.0716711,0.0734173,0.0782672,0.0791601,0.0842406};
	VEC GeV62_v2_stat_error_PiP_0_10 = {0.000105537,0.000122407,0.000164809,0.000228831,0.000319575,0.000446836,0.000629297,0.000882756,0.00122416,0.00165379,0.00173078,0.00262822,0.00572041};

	// ---------------- Particle species: KM ----------------
	VEC GeV62_pt_bin_center_KM_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.465};
	VEC GeV62_v2_values_KM_0_10 = {0.00283136,0.0114375,0.0240432,0.035912,0.0451036,0.0538874,0.0607677,0.0660069,0.0743655,0.0721627,0.0730863,0.0696812,0.0949473};
	VEC GeV62_v2_stat_error_KM_0_10 = {0.000727205,0.000438211,0.000428096,0.000490031,0.000599288,0.000766747,0.00101962,0.00137493,0.00190889,0.00269757,0.00316511,0.00583505,0.0060618};

	// ---------------- Particle species: KP ----------------
	VEC GeV62_pt_bin_center_KP_0_10 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.355,2.805,3.465};
	VEC GeV62_v2_values_KP_0_10 = {0.00455492,0.0123065,0.0246732,0.0357745,0.046377,0.0525499,0.0599756,0.066733,0.0694622,0.0716595,0.0781367,0.0813609,0.0566004};
	VEC GeV62_v2_stat_error_KP_0_10 = {0.000663629,0.000399441,0.000388713,0.000442513,0.000539863,0.000689017,0.000909302,0.00120959,0.00164816,0.00230713,0.00261841,0.00489972,0.0126054};

	// ---------------- Particle species: Proton ----------------
	VEC GeV62_pt_bin_center_Pr_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV62_v2_values_Pr_10_40 = {0.00684454,0.0151163,0.0325327,0.0558758,0.0818952,0.106877,0.129578,0.148797,0.16564,0.181987,0.194254,0.199881,0.211623,0.2255,0.220917};
	VEC GeV62_v2_stat_error_Pr_10_40 = {0.000883377,0.000304589,0.000248505,0.000246802,0.000271826,0.0003192,0.000390438,0.000499595,0.000668752,0.000798948,0.00106583,0.00150297,0.00217511,0.00255858,0.00493119};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV62_pt_bin_center_Ant_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV62_v2_values_Ant_10_40 = {0.00655641,0.0108455,0.0268988,0.0497209,0.0742441,0.0979394,0.121136,0.14217,0.160819,0.171942,0.183975,0.193206,0.198831,0.223087,0.214498};
	VEC GeV62_v2_stat_error_Ant_10_40 = {0.00120701,0.000469091,0.000375989,0.000372987,0.000413322,0.000490455,0.000606019,0.000787414,0.00093507,0.00125693,0.00170006,0.00246151,0.00380559,0.004524,0.00870049};

	// ---------------- Particle species: PiM ----------------
	VEC GeV62_pt_bin_center_PiM_10_40 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.685,2.955,3.495};
	VEC GeV62_v2_values_PiM_10_40 = {0.0256089,0.049558,0.0729376,0.0940734,0.112055,0.126778,0.138555,0.146616,0.154424,0.158706,0.163852,0.164639,0.1638,0.16395,0.158083};
	VEC GeV62_v2_stat_error_PiM_10_40 = {5.53111e-05,6.55598e-05,8.90389e-05,0.000124527,0.000174859,0.000243997,0.000342874,0.000412283,0.000659089,0.00087937,0.00116265,0.00151799,0.00199145,0.00209862,0.00303673};

	// ---------------- Particle species: PiP ----------------
	VEC GeV62_pt_bin_center_PiP_10_40 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.495};
	VEC GeV62_v2_values_PiP_10_40 = {0.0246568,0.0488236,0.0724498,0.0933431,0.111199,0.125826,0.137535,0.145903,0.153263,0.159387,0.161891,0.1649,0.164148,0.164993,0.163504};
	VEC GeV62_v2_stat_error_PiP_10_40 = {5.5493e-05,6.57157e-05,8.89935e-05,0.000124016,0.000173554,0.000241259,0.000338897,0.000471602,0.000648294,0.000870067,0.00114163,0.00149154,0.00194562,0.00201045,0.00292764};

	// ---------------- Particle species: KM ----------------
	VEC GeV62_pt_bin_center_KM_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.465};
	VEC GeV62_v2_values_KM_10_40 = {0.0137527,0.032028,0.0553626,0.0779299,0.0969846,0.113642,0.127522,0.136508,0.145492,0.15104,0.15465,0.151516,0.158103,0.155772,0.149941};
	VEC GeV62_v2_stat_error_KM_10_40 = {0.00038068,0.000234438,0.000233434,0.000270744,0.000333662,0.000427493,0.000564609,0.00075668,0.00103883,0.00147094,0.00191321,0.00293076,0.00406229,0.00323258,0.00604343};

	// ---------------- Particle species: KP ----------------
	VEC GeV62_pt_bin_center_KP_10_40 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.495};
	VEC GeV62_v2_values_KP_10_40 = {0.0142023,0.0323987,0.0559368,0.0782457,0.0976515,0.112853,0.125878,0.137557,0.144958,0.1499,0.152787,0.156018,0.152189,0.149711,0.143515};
	VEC GeV62_v2_stat_error_KP_10_40 = {0.000348329,0.000214471,0.000212767,0.00024557,0.000301606,0.000385246,0.000504597,0.000668958,0.000906547,0.00124885,0.00174665,0.00243852,0.00329338,0.00395936,0.00568895};

	// ---------------- Particle species: Proton ----------------
	VEC GeV62_pt_bin_center_Pr_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
	VEC GeV62_v2_values_Pr_40_80 = {0.0191101,0.0323186,0.0582973,0.0891716,0.124343,0.155709,0.185097,0.211798,0.229861,0.248513,0.256268,0.274427,0.281107,0.295072,0.286854};
	VEC GeV62_v2_stat_error_Pr_40_80 = {0.00258841,0.000930903,0.000807939,0.000854547,0.000995694,0.00122282,0.00154932,0.00203734,0.00276412,0.00328633,0.00436492,0.00607112,0.00854475,0.0101247,0.0197444};

	// ---------------- Particle species: antiProton ----------------
	VEC GeV62_pt_bin_center_Ant_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.535,3.045};
	VEC GeV62_v2_values_Ant_40_80 = {0.0166509,0.0287498,0.0519787,0.0854519,0.117946,0.149035,0.179829,0.206576,0.221758,0.248716,0.258421,0.277343,0.285936};
	VEC GeV62_v2_stat_error_Ant_40_80 = {0.00331704,0.00135326,0.00116262,0.00123692,0.00145955,0.00182167,0.00235109,0.00316641,0.00380031,0.00513822,0.00697468,0.00862211,0.0165521};

	// ---------------- Particle species: PiM ----------------
	VEC GeV62_pt_bin_center_PiM_40_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.665,1.875,2.085,2.295,2.475,2.835};
	VEC GeV62_v2_values_PiM_40_80 = {0.031187,0.0631613,0.0945519,0.122517,0.146542,0.166996,0.18207,0.192386,0.204243,0.208179,0.217274,0.229866,0.221771};
	VEC GeV62_v2_stat_error_PiM_40_80 = {0.00018352,0.000224601,0.000312148,0.000445396,0.000634818,0.000898365,0.0012675,0.00153467,0.00239527,0.00316037,0.00412301,0.00537077,0.00505901};

	// ---------------- Particle species: PiP ----------------
	VEC GeV62_pt_bin_center_PiP_40_80 = {0.285,0.495,0.675,0.885,1.065,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.835};
	VEC GeV62_v2_values_PiP_40_80 = {0.0302612,0.0618884,0.0931795,0.120907,0.14497,0.16514,0.181052,0.191616,0.197537,0.212764,0.210707,0.217859,0.23746};
	VEC GeV62_v2_stat_error_PiP_40_80 = {0.0001836,0.000225071,0.000312447,0.000444503,0.000630998,0.000892833,0.00125877,0.00155694,0.00237389,0.00312628,0.00405695,0.00530303,0.00493076};

	// ---------------- Particle species: KM ----------------
	VEC GeV62_pt_bin_center_KM_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.835};
	VEC GeV62_v2_values_KM_40_80 = {0.0229626,0.0462706,0.0753394,0.102919,0.126531,0.144455,0.156266,0.17241,0.180039,0.182719,0.198842,0.19519,0.209907};
	VEC GeV62_v2_stat_error_KM_40_80 = {0.00124005,0.000796974,0.000825902,0.000989555,0.00125127,0.00163809,0.00218805,0.00295023,0.00407447,0.00572561,0.00720873,0.0111296,0.0103606};

	// ---------------- Particle species: KP ----------------
	VEC GeV62_pt_bin_center_KP_40_80 = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.835};
	VEC GeV62_v2_values_KP_40_80 = {0.023864,0.0477082,0.0741166,0.102174,0.125623,0.141866,0.157419,0.169195,0.177781,0.183851,0.194883,0.197105,0.200496};
	VEC GeV62_v2_stat_error_KP_40_80 = {0.0011443,0.000734718,0.000759666,0.000903888,0.00113393,0.00147602,0.0019548,0.00259636,0.00349887,0.00470513,0.00645972,0.0087786,0.00857531};

	if(Energy==62){
		DrawArticle( Energy, GeV62_pt_bin_center_Ant_0_80, GeV62_v2_values_Ant_0_80, GeV62_v2_stat_error_Ant_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Neg", 0.33,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_Pr_0_80, GeV62_v2_values_Pr_0_80, GeV62_v2_stat_error_Pr_0_80, pt_bin_err, 0, 8, "Eta01","Proton", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_KP_0_80, GeV62_v2_values_KP_0_80, GeV62_v2_stat_error_KP_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Pos", 0.28,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_KM_0_80, GeV62_v2_values_KM_0_80, GeV62_v2_stat_error_KM_0_80, pt_bin_err, 0, 8, "Eta01","Kaon", "Neg", 0.28,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_PiP_0_80, GeV62_v2_values_PiP_0_80, GeV62_v2_stat_error_PiP_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Pos", 0.28,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_PiM_0_80, GeV62_v2_values_PiM_0_80, GeV62_v2_stat_error_PiM_0_80, pt_bin_err, 0, 8, "Eta01","Pion", "Neg", 0.28,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_Ant_0_10, GeV62_v2_values_Ant_0_10, GeV62_v2_stat_error_Ant_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Neg", 0.22,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_Pr_0_10, GeV62_v2_values_Pr_0_10, GeV62_v2_stat_error_Pr_0_10, pt_bin_err, 7, 8, "Eta01","Proton", "Pos", 0.22,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_KP_0_10, GeV62_v2_values_KP_0_10, GeV62_v2_stat_error_KP_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Pos", 0.16,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_KM_0_10, GeV62_v2_values_KM_0_10, GeV62_v2_stat_error_KM_0_10, pt_bin_err, 7, 8, "Eta01","Kaon", "Neg", 0.16,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_PiP_0_10, GeV62_v2_values_PiP_0_10, GeV62_v2_stat_error_PiP_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Pos", 0.14,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_PiM_0_10, GeV62_v2_values_PiM_0_10, GeV62_v2_stat_error_PiM_0_10, pt_bin_err, 7, 8, "Eta01","Pion", "Neg", 0.14,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_Pr_10_40, GeV62_v2_values_Pr_10_40, GeV62_v2_stat_error_Pr_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_Ant_10_40, GeV62_v2_values_Ant_10_40, GeV62_v2_stat_error_Ant_10_40, pt_bin_err, 4, 6, "Eta01","Proton", "Neg", 0.33,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_KP_10_40, GeV62_v2_values_KP_10_40, GeV62_v2_stat_error_KP_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Pos", 0.28,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_KM_10_40, GeV62_v2_values_KM_10_40, GeV62_v2_stat_error_KM_10_40, pt_bin_err, 4, 6, "Eta01","Kaon", "Neg", 0.28,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_PiP_10_40, GeV62_v2_values_PiP_10_40, GeV62_v2_stat_error_PiP_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Pos", 0.25,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_PiM_10_40, GeV62_v2_values_PiM_10_40, GeV62_v2_stat_error_PiM_10_40, pt_bin_err, 4, 6, "Eta01","Pion", "Neg", 0.25,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_Pr_40_80, GeV62_v2_values_Pr_40_80, GeV62_v2_stat_error_Pr_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Pos", 0.42,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_Ant_40_80, GeV62_v2_values_Ant_40_80, GeV62_v2_stat_error_Ant_40_80, pt_bin_err, 0, 3, "Eta01","Proton", "Neg", 0.42,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_KP_40_80, GeV62_v2_values_KP_40_80, GeV62_v2_stat_error_KP_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_KM_40_80, GeV62_v2_values_KM_40_80, GeV62_v2_stat_error_KM_40_80, pt_bin_err, 0, 3, "Eta01","Kaon", "Neg", 0.33,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_PiP_40_80, GeV62_v2_values_PiP_40_80, GeV62_v2_stat_error_PiP_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Pos", 0.33,3.6);
		DrawArticle( Energy, GeV62_pt_bin_center_PiM_40_80, GeV62_v2_values_PiM_40_80, GeV62_v2_stat_error_PiM_40_80, pt_bin_err, 0, 3, "Eta01","Pion", "Neg", 0.33,3.6);
	
		//DrawArticle( Energy, pt62GeV, GeV62_v2_10_20, GeV62_v2_10_20_err, ptGeV_err, 6, 6, "Eta15","Hadrons", "", 0.33,3.6);
		//DrawArticle( Energy, pt62GeV, GeV62_v2_20_30, GeV62_v2_20_30_err, ptGeV_err, 5, 5, "Eta15","Hadrons", "", 0.33,3.6);
		//DrawArticle( Energy, pt62GeV, GeV62_v2_30_40, GeV62_v2_30_40_err, ptGeV_err, 4, 4, "Eta15","Hadrons", "", 0.33,3.6);
		//DrawArticle( Energy, pt62GeV, GeV62_v2_50_60, GeV62_v2_50_60_err, ptGeV_err, 2, 2, "Eta15","Hadrons", "", 0.33,3.6);

	}

}


void makeplotstyle(){
  
    TStyle *mystyle = new TStyle("PlottingInStyle", "Style for Summary Plots");
    mystyle->SetPalette(1);
    mystyle->SetCanvasColor(10);
    mystyle->SetHistFillColor(10);
    mystyle->SetHistFillStyle(0);
    mystyle->SetOptTitle(0);
    mystyle->SetOptStat(0);
    mystyle->SetCanvasBorderMode(0);//removes the yellow frame around the canvas
    mystyle->SetPadLeftMargin(0.16);
    mystyle->SetPadBottomMargin(0.15);
    mystyle->SetPadTickX(1);
    mystyle->SetPadTickY(1);
    mystyle->SetAxisColor(1, "X");
    mystyle->SetAxisColor(1, "Y");
    mystyle->SetLabelColor(1, "X");
    mystyle->SetLabelColor(1, "Y");
    mystyle->SetTickLength(0.03, "X");
    mystyle->SetTickLength(0.03, "Y");
    mystyle->SetTitleXSize(0.05);
    mystyle->SetTitleYSize(0.05);
    mystyle->SetNdivisions(508, "X");
    mystyle->SetNdivisions(505, "Y");
    mystyle->SetTitleXOffset(1.2);
    mystyle->SetTitleYOffset(1.4);
    mystyle->SetLabelOffset(0.02, "X");
    mystyle->SetLabelOffset(0.02, "Y");
    mystyle->SetLabelSize(0.05, "X");
    mystyle->SetLabelSize(0.05, "Y");
    //mystyle->SetGridx();

  	TFile f("style.root", "RECREATE");
  	f.cd();
  	mystyle->Write();
  	f.Close();
}
/*
void ParAntFlowVsPtForBES( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif){

	VEC rebin_vec;
	Int_t point = 1; 
	TString prefix = "TPCandTOF";
	
	Int_t Arr_energy[6]={11,14,19,27,39,62};
	const Int_t color[]={2, 4, 1};
	const Int_t style_par[]={33, 21, 22};
	const Int_t style_ant[]={27, 25, 26};
	Float_t size[3]={1.4,1.4,1.4};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasFlow = new Draw_Picture_new;
	Draw_Picture_new *canvasDifference = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow;
	std::vector<std::vector<Draw_Object *>> ratio;
	std::vector<std::vector<Draw_Object *>> difference;
	std::vector<Draw_Object *> object;
	
	Double_t coord_y_draw_ratio = min_y_axis_ratio + (max_y_axis_ratio - min_y_axis_ratio)*0.8;
	Double_t coord_y_draw_dif = min_y_axis_dif + (max_y_axis_dif - min_y_axis_dif)*0.8;
	
	Int_t b=0;

	for(Int_t j=0; j<6; j++){
		
		canvasFlow->SetTextInfo(Form("#font[42]{%s GeV}", Energy_scan_text.at(Arr_energy[j])), "Pad", 0.25, 0.75*max_y_axis, 0.13, 0.0);
		canvasDifference->SetTextInfo(Form("#font[42]{%s GeV}", Energy_scan_text.at(Arr_energy[j])), "Pad", 0.25, coord_y_draw_dif, 0.13, 0.0);
		canvasRatio->SetTextInfo(Form("#font[42]{%s GeV}", Energy_scan_text.at(Arr_energy[j])), "Pad", 0.25, coord_y_draw_ratio, 0.13, 0.0);
		
		for(Int_t harm=0; harm<2; harm++){

			rebin_vec = Rebin2( harm+2, Arr_energy[j], PID, "", CentBinMin, CentBinMax);
			
			object.push_back(new Draw_Object);
			Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, PID, "Pos", EtaGap, prefix);
			Analysis -> FlowVsPt_EtaSub(object.back(), CentBinMin, CentBinMax, rebin_vec);
    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.2]{%s}}", particle[GetNumberParticle(PID,"Pos")].Data()), color[0], style_par[0], 1.6);
			
			object.push_back(new Draw_Object);
			Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, PID, "Neg", EtaGap, prefix);
			Analysis -> FlowVsPt_EtaSub(object.back(), CentBinMin, CentBinMax, rebin_vec);
    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.2]{%s}}", particle[GetNumberParticle(PID,"Neg")].Data()), color[1], style_ant[0], 1.6);
			
			flow.push_back(object);
			object.clear();

			object.push_back(new Draw_Object);
			Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.2]{%s}}", particle[GetNumberParticle(PID,"Pos")].Data()), color[0], style_par[0], size[0]);
    		difference.push_back(object);
    		object.clear();

			object.push_back(new Draw_Object);
			Analysis -> RatioParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.2]{%s}}", particle[GetNumberParticle(PID,"Pos")].Data()), color[0], style_par[0], size[0]);
    		ratio.push_back(object);
    		object.clear();

			if(j==5){
				canvasFlow->SetTextInfo(Form("#font[42]{v_{%i}}", harm+2), Form("V%i",harm+2), 0.9*max_x_axis, 0.75*max_y_axis, 0.13, 0.0);
				canvasDifference->SetTextInfo(Form("#font[42]{v_{%i}}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis,coord_y_draw_dif, 0.13, 0.0);
				canvasRatio->SetTextInfo(Form("#font[42]{v_{%i}}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis, coord_y_draw_ratio, 0.13, 0.0);
			}

			canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.75*max_y_axis, 0.13, 0.0);
			canvasDifference->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, coord_y_draw_dif, 0.13, 0.0);
			canvasRatio->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, coord_y_draw_ratio, 0.13, 0.0);
			b++;

		}
	}
	canvasFlow->SetTextInfo("#font[42]{v_{n}}","Y",0.25,0.6,0.5,0.0);
	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",1.1,0.75*max_y_axis,0.13,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.13);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.15,0.5,0.35,0.92,1);
	TCanvas *result = canvasFlow->CanvasNxM(1000,1280, 6, 2, 0.4, 0.68, flow,1,0,1);

	result->SaveAs(Form("%s/picture/3_PID/BES_vn_%s_%s_cent_%i_%i.pdf",outPath, PID.Data(), EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result;

	canvasRatio->SetTextInfo(Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle(PID,"Pos")].Data(), particle[GetNumberParticle(PID,"Neg")].Data()),"Y",0.07,0.7,0.38,0);
	canvasRatio->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",1.1,coord_y_draw_ratio,0.13,0.0);
	canvasRatio->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio, max_y_axis_ratio, 0.13);
	canvasRatio->SetDrawObject(ratio[0],"OnePad");
	//canvasRatio->SetLegend(0.1,0.4,0.3,0.92,1,0,1);
	TCanvas *result2 = canvasRatio->CanvasNxM(1000,1280, 6, 2, 0.4, 0.68, ratio,1,0,1);

	result2->SaveAs(Form("%s/picture/3_PID/BES_vn_ratio_%s_%s_cent_%i_%i.pdf",outPath, PID.Data(), EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result2;

	canvasDifference->SetTextInfo(Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle(PID,"Pos")].Data(),particle[GetNumberParticle(PID,"Neg")].Data() ),"Y",0.4, 0.5,0.4,90);
	canvasDifference->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasDifference->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",1.1,coord_y_draw_dif,0.13,0.0);
	canvasDifference->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasDifference->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif, max_y_axis_dif, 0.13);
	canvasDifference->SetDrawObject(difference[0],"OnePad");
	//canvasDifference->SetLegend(0.1,0.7,0.3,0.92,1,0,1);
	TCanvas *result3 = canvasDifference->CanvasNxM(1150,1280, 6, 2, 0.4, 0.68, difference,1,0,1);

	result3->SaveAs(Form("%s/picture/3_PID/BES_vn_diff_%s_%s_cent_%i_%i.pdf",outPath, PID.Data(), EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result3;

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}
*/