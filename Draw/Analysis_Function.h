/* picture.h */
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
#include <algorithm>

#include "./Draw_Object.h"
#include "./Draw_Picture_new.h"

typedef vector<Double_t> Vec;

class Analysis_Function {
    public:
        // Деструктор
        ~Analysis_Function();

        Bool_t Sys_mod=false;
        Bool_t Sys_mod_DifParAnr=false;
        Bool_t Sys_mod_RatParAnr=false;
        Bool_t Sys_mod_Integral=false;

        TString NameInputFile;
        TString NameOutputFile;
        TString NameTProfileFlow;
        TString NameTProfileResolution;
        TString NameTProfileMeanPt;
        
        TString PID;
        TString charge;
        TString name_sys;
        TString method_PID;
		TString NameGraph;

		Double_t x_min_axis;
		Double_t x_max_axis;
		Double_t y_min_axis;
		Double_t y_max_axis;

		Int_t energy;
        Int_t harmonic;

        TFile *fileIn;
        TFile *fileOut;

        std::map<TString, std::vector<TString>> Sys_string = { {"Hadrons",{"Eta15","Eta03","Eta05","Eta07"}}, {"Pion",{"Eta01","Eta03","Eta05","Eta07"}}
    															, {"Kaon",{"Eta01","Eta03","Eta05","Eta07"}}, {"Proton",{"Eta01","Eta03","Eta05","Eta07"}}};
        std::map<TString, std::vector<Double_t>> Sys_double = {{"Hadrons",{0.075, 0.15, 0.25, 0.35}},         {"PID",{0.05, 0.15, 0.25, 0.35}}};

        void SetParametrs(Int_t _energy, TFile *file_in, TFile *file_out, int _harmonic, TString _PID, TString _charge, TString _name_sys, TString _method_PID  ){

       		PID = _PID;
        	charge = _charge;
        	name_sys = _name_sys;
        	method_PID = _method_PID;
        	harmonic = _harmonic;
        	fileIn = file_in;
        	fileOut = file_out;
        	energy = _energy;
        
        }

        void SetNameTProfile2DFlow(){

        	NameTProfileFlow = "tp_v";
        	NameTProfileFlow += harmonic;
        	NameTProfileFlow += "ewTPC";
        	NameTProfileFlow += PID;
        	NameTProfileFlow += charge;
        	NameTProfileFlow += name_sys;
        	NameTProfileFlow += method_PID;
        	//NameTProfileFlow += "_r";

        }

        void SetNameTProfile2DFlowEta(Int_t CentBin){

        	NameTProfileFlow = "tp_v";
        	NameTProfileFlow += harmonic;
        	NameTProfileFlow += "ewTPC";
        	NameTProfileFlow += PID;
        	NameTProfileFlow += charge;
        	NameTProfileFlow += name_sys;
        	NameTProfileFlow += "Eta";
        	NameTProfileFlow += method_PID;
        	NameTProfileFlow += "CentBin";
        	NameTProfileFlow += CentBin;

        }

        void SetNameTProfile2DMeanPt(){
        	
        	NameTProfileMeanPt = "tp_meanPt";
        	NameTProfileMeanPt += PID;
        	NameTProfileMeanPt += charge;
        	NameTProfileMeanPt += name_sys;
        	NameTProfileMeanPt += "TPCandTOF";
        	//NameTProfileMeanPt += method_PID;
        }

        void SetNameTProfileResolution(){

        	NameTProfileResolution = "tp_SqRes";
        	NameTProfileResolution += harmonic;
        	NameTProfileResolution += "TPC";
        	NameTProfileResolution += name_sys;
        }

        void SetNameTGraph(){

        	NameGraph = "gr_";
			NameGraph += energy;
			NameGraph += "GeV_v";
			NameGraph += harmonic;
			NameGraph += PID;
        	NameGraph += charge;
        	NameGraph += name_sys;
        	NameGraph += method_PID;

        }

        void ResolutionEP_EtaSub(){

			SetNameTProfileResolution();

			std::vector<Double_t> x = {2.5,7.5,15,25,35,45,55,65,75};
			std::vector<Double_t> y;
			std::vector<Double_t> x_err;
			std::vector<Double_t> y_err;

			TProfile *tp_Res = (TProfile*)fileIn -> Get(NameTProfileResolution.Data());

			for(Int_t i=9; i>0; i--){
				x_err.push_back(0.0);
				y.push_back( sqrt(abs(tp_Res -> GetBinContent(i))));
				y_err.push_back(  tp_Res -> GetBinError(i) / sqrt(abs(tp_Res -> GetBinContent(i))) );
			}

			//result->SetVectorsGraph(x, y, x_err, y_err);

			for(Int_t i=0; i<9; i++){
				std::cout<<y[i]<<",\t";
			}
			std::cout<<"\n\n";

			for(Int_t i=0; i<9; i++){
				std::cout<<y_err[i]<<",\t";
			}
			std::cout<<"\n";
			//return SetTGraph(NameTProfileResolution,TitleGraph,x,y,x_err,y_err, color, style, size_marker);

		}

		void FlowVsPt_EtaSub(Draw_Object *result, Int_t CentBinMin, Int_t CentBinMax,  Vec rebin_vector ){

			//SetInputFile();
			SetNameTProfile2DFlow();
			SetNameTProfile2DMeanPt();
			SetNameTProfileResolution();
			SetNameTGraph();

			//std::cout<<"\n"<<NameTProfileFlow<<"\n"<<NameTProfileMeanPt<<"\n"<<NameTProfileResolution<<"\n";
			std::vector<Double_t> pt;
			std::vector<Double_t> pt_err;
			std::vector<Double_t> flow;
			std::vector<Double_t> flow_err;
			
			Int_t n_rebin;
			Double_t rebin[100];
			ConvertVectorToArrey(rebin_vector, n_rebin, rebin);

			TProfile *tp_vn_bin[9];

			TProfile2D *tp_vn = (TProfile2D*)fileIn -> Get(NameTProfileFlow.Data());
			TProfile2D *tp_pt = (TProfile2D*)fileIn->Get(NameTProfileMeanPt.Data());
			TProfile *tp_Res = (TProfile*)fileIn -> Get(NameTProfileResolution.Data());
			
			// set pt axis
			TProfile *tp_pt2 = (TProfile*)tp_pt->ProfileX(Form("prof_%s",NameTProfileMeanPt.Data()),CentBinMin+1,CentBinMax+1);
			TH1 *tp_meanPt = tp_pt2->Rebin(n_rebin,Form("h2_%s",NameTProfileMeanPt.Data()),rebin);
			//set v_n axis
			for(Int_t j=CentBinMin;  j<=CentBinMax; j++){
				tp_vn_bin[j] = (TProfile*)tp_vn -> ProfileX(Form("profX_%s_cent%i%i", NameTProfileFlow.Data(), j, j),j+1,j+1);
				tp_vn_bin[j] -> Scale( 1.0 / TMath::Sqrt( abs(tp_Res->GetBinContent(tp_Res->FindBin((Double_t)j)))) );
			}

			if(CentBinMin != CentBinMax){
				for(Int_t j=CentBinMin+1; j<=CentBinMax; j++){
					tp_vn_bin[CentBinMin] -> Add(tp_vn_bin[j]);
				}
			}
			TH1 *tp_flow = tp_vn_bin[CentBinMin]->Rebin(n_rebin, NameTProfileFlow.Data(), rebin);
			
			for(Int_t i=0; i < (int)tp_meanPt->GetNbinsX()-1; i++){
				pt.push_back(tp_meanPt->GetBinContent(i+1));
				pt_err.push_back(0.);
				flow.push_back(tp_flow->GetBinContent(tp_flow->FindBin(pt[i])));
				flow_err.push_back(tp_flow->GetBinError(tp_flow->FindBin(pt[i])));
			}
			
			result->SetVectorsGraph(pt, flow, pt_err, flow_err);
			
			if(Sys_mod==true){
				for(Int_t i=0; i < (int)Sys_string.at(PID).size(); i++){
					FlowVsPt_EtaSubSysErrors(result, CentBinMin, CentBinMax, rebin_vector, Sys_string.at(PID)[i]);
				}
				result->SetSysErrors();
			}
			
			for(Int_t i=CentBinMin; i<=CentBinMax; i++){
				delete tp_vn_bin[i];
			}
			delete tp_meanPt;
			delete tp_pt2;

		}


		void FlowVsPt_EtaSubSysErrors(Draw_Object *result, Int_t CentBinMin, Int_t CentBinMax,  Vec rebin_vector , TString Sys_name){

			SetNameTProfile2DMeanPt();
			name_sys = Sys_name;
			SetNameTProfile2DFlow();
			SetNameTProfileResolution();
			SetNameTGraph();

			//std::cout<<"\n"<<NameTProfileFlow<<"\n"<<NameTProfileMeanPt<<"\n"<<NameTProfileResolution<<"\n";
			std::vector<Double_t> pt;
			std::vector<Double_t> flow;
			
			Int_t n_rebin;
			Double_t rebin[100];
			ConvertVectorToArrey(rebin_vector, n_rebin, rebin);

			TProfile *tp_vn_bin[9];

			TProfile2D *tp_vn = (TProfile2D*)fileIn -> Get(NameTProfileFlow.Data());
			TProfile2D *tp_pt = (TProfile2D*)fileIn->Get(NameTProfileMeanPt.Data());
			TProfile *tp_Res = (TProfile*)fileIn -> Get(NameTProfileResolution.Data());
			
			// set pt axis
			TProfile *tp_pt2 = (TProfile*)tp_pt->ProfileX(Form("prof_%s_sys",NameTProfileMeanPt.Data()),CentBinMin+1,CentBinMax+1);
			TH1 *tp_meanPt = tp_pt2->Rebin(n_rebin,Form("h2_%s_sys",NameTProfileMeanPt.Data()),rebin);
			//set v_n axis
			for(Int_t j=CentBinMin;  j<=CentBinMax; j++){
				tp_vn_bin[j] = (TProfile*)tp_vn -> ProfileX(Form("profX_%s_cent%i%i_sys", NameTProfileFlow.Data(), j, j),j+1,j+1);
				tp_vn_bin[j] -> Scale( 1.0 / TMath::Sqrt( abs(tp_Res->GetBinContent(tp_Res->FindBin((Double_t)j)))) );
			}

			if(CentBinMin != CentBinMax){
				for(Int_t j=CentBinMin+1; j<=CentBinMax; j++){
					tp_vn_bin[CentBinMin] -> Add(tp_vn_bin[j]);
				}
			}
			TH1 *tp_flow = tp_vn_bin[CentBinMin]->Rebin(n_rebin, NameTProfileFlow.Data(), rebin);
			
			for(Int_t i=0; i < (int)tp_meanPt->GetNbinsX()-1; i++){
				pt.push_back(tp_meanPt->GetBinContent(i+1));
				flow.push_back(tp_flow->GetBinContent(tp_flow->FindBin(pt[i])));
			}
			
			result->SetSysArrey(flow);
			
			for(Int_t i=CentBinMin; i<=CentBinMax; i++){
				delete tp_vn_bin[i];
			}
			delete tp_meanPt;
			delete tp_pt2;

		}

		void FlowVsPt_SP(Draw_Object *result, Int_t CentBinMin, Int_t CentBinMax,  Vec rebin_vector, TString axis ){

			//SetInputFile();
			SetNameTProfile2DFlow();
			SetNameTProfile2DMeanPt();
			SetNameTProfileResolution();
			SetNameTGraph();

			NameTProfileFlow += axis;
			NameTProfileResolution += axis;

			std::cout<<"\n"<<NameTProfileFlow<<"\n"<<NameTProfileMeanPt<<"\n"<<NameTProfileResolution<<"\n";
			std::vector<Double_t> pt;
			std::vector<Double_t> pt_err;
			std::vector<Double_t> flow;
			std::vector<Double_t> flow_err;
			
			Int_t n_rebin;
			Double_t rebin[100];
			ConvertVectorToArrey(rebin_vector, n_rebin, rebin);

			TProfile *tp_vn_bin[9];

			TProfile2D *tp_vn = (TProfile2D*)fileIn -> Get(NameTProfileFlow.Data());
			TProfile2D *tp_pt = (TProfile2D*)fileIn->Get(NameTProfileMeanPt.Data());
			TProfile *tp_Res = (TProfile*)fileIn -> Get(NameTProfileResolution.Data());
			
			// set pt axis
			TProfile *tp_pt2 = (TProfile*)tp_pt->ProfileX(Form("prof_%s",NameTProfileMeanPt.Data()),CentBinMin+1,CentBinMax+1);
			TH1 *tp_meanPt = tp_pt2->Rebin(n_rebin,Form("h2_%s",NameTProfileMeanPt.Data()),rebin);
			//set v_n axis
			for(Int_t j=CentBinMin;  j<=CentBinMax; j++){
				tp_vn_bin[j] = (TProfile*)tp_vn -> ProfileX(Form("profX_%s_cent%i%i", NameTProfileFlow.Data(), j, j),j+1,j+1);
				if(strncmp(axis, "X",1)==0 || strncmp(axis, "Y",1)==0){
					tp_vn_bin[j] -> Scale( 1.0 / TMath::Sqrt( 2.0 * abs(tp_Res->GetBinContent(tp_Res->FindBin((Double_t)j)))) );
				}else{
					tp_vn_bin[j] -> Scale( 1.0 / (2.0 * TMath::Sqrt( abs(tp_Res->GetBinContent(tp_Res->FindBin((Double_t)j))))) );
				}
			}	

			if(CentBinMin != CentBinMax){
				for(Int_t j=CentBinMin+1; j<=CentBinMax; j++){
					tp_vn_bin[CentBinMin] -> Add(tp_vn_bin[j]);
				}
			}
			TH1 *tp_flow = tp_vn_bin[CentBinMin]->Rebin(n_rebin, NameTProfileFlow.Data(), rebin);
			
			for(Int_t i=0; i < (int)tp_meanPt->GetNbinsX()-1; i++){
				pt.push_back(tp_meanPt->GetBinContent(i+1));
				pt_err.push_back(0.);
				flow.push_back(tp_flow->GetBinContent(tp_flow->FindBin(pt[i])));
				flow_err.push_back(tp_flow->GetBinError(tp_flow->FindBin(pt[i])));
			}
			
			result->SetVectorsGraph(pt, flow, pt_err, flow_err);
			
			for(Int_t i=CentBinMin; i<=CentBinMax; i++){
				delete tp_vn_bin[i];
			}
			delete tp_meanPt;
			delete tp_pt2;

		}

/*
		void FlowVsCent_EtaSub(Draw_Object *result, Double_t ptMin, Double_t ptMax,  Vec rebin_vector_bin, Vec rebin_vector_percent ){

			//SetInputFile();
			SetNameTProfile2DFlow();
			SetNameTProfileResolution();
			SetNameTGraph();
			NameTProfileFlow+="noWeight";

			//std::cout<<"\n"<<NameTProfileFlow<<"\n"<<NameTProfileMeanPt<<"\n"<<NameTProfileResolution<<"\n";
			std::vector<Double_t> cent;
			std::vector<Double_t> cent_err;
			std::vector<Double_t> flow;
			std::vector<Double_t> flow_err;

			Int_t nBinsPt=100;
			Double_t pt_min = 0.15;
			Double_t pt_max = 5.15;
			Double_t CenaDeleniy = (pt_max - pt_min)/100.;
			Int_t ptBinMin = (int)((ptMin-0.15) / CenaDeleniy) + 1;
			Int_t ptBinMax = (int)((ptMax-0.15) / CenaDeleniy);
			
			Int_t n_rebin;
			Double_t rebin[100];
			ConvertVectorToArrey(rebin_vector, n_rebin, rebin);

			TProfile2D *tp_vn = (TProfile2D*)fileIn -> Get(NameTProfileFlow.Data());
			TProfile *tp_Res = (TProfile*)fileIn -> Get(NameTProfileResolution.Data());
			TProfile *tp_flow = (TProfile*)tp_vn -> ProfileY(Form("profX_%s_cent%i%i", NameTProfileFlow.Data(), ptBinMin, ptBinMax),ptBinMin,ptBinMax);
			
			TH1 *h_FlowCentBin = tp_flow->Rebin(n_rebin,Form("h2_%s",NameTProfileFlow.Data()),rebin);
			TH1 *h_Res = tp_Res->Rebin(n_rebin,Form("h2_%s",NameTProfileFlow.Data()),rebin); 

			for(Int_t i=(int)tp_flow->GetNbinsX(); i>0; i--){
				cent.push_back(tp_Res->GetBinCenter(i));
				cent_err.push_back(0.);
				flow.push_back(tp_flow->GetBinContent(tp_Res->FindBin(cent[i])) / TMath::Sqrt( abs(tp_Res->GetBinContent(tp_Res->FindBin((Double_t)i)))));
				flow_err.push_back(tp_flow->GetBinError(tp_flow->FindBin(cent[i])));
			}

			//set v_n axis
			for(Int_t j=CentBinMin;  j<=CentBinMax; j++){
				tp_vn_bin[j] = (TProfile*)tp_vn -> ProfileY(Form("profX_%s_pt%i%i", NameTProfileFlow.Data(), ptBinMin, ptBinMax),ptBinMin,j+1);
				tp_vn_bin[j] -> Scale( 1.0 / TMath::Sqrt( abs(tp_Res->GetBinContent(tp_Res->FindBin((Double_t)j)))) );
			}

		}
*/
		void FlowVsEta_EtaSub(Draw_Object *result, Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, Vec vec_rebin ){

			//SetInputFile();
			SetNameTProfileResolution();
			SetNameTGraph();

			//std::cout<<"\n"<<NameTProfileFlow<<"\n"<<NameTProfileMeanPt<<"\n"<<NameTProfileResolution<<"\n";
			std::vector<Double_t> eta;
			std::vector<Double_t> eta_err;
			std::vector<Double_t> flow;
			std::vector<Double_t> flow_err;

			Int_t nBinsPt=100;
			Double_t pt_min = 0.15;
			Double_t pt_max = 5.15;
			Double_t CenaDeleniy = (pt_max - pt_min)/100.;
			Int_t ptBinMin = (int)((ptMin-0.15) / CenaDeleniy) + 1;
			Int_t ptBinMax = (int)((ptMax-0.15) / CenaDeleniy);
			
			Int_t n_rebin;
			Double_t rebin[100];
			ConvertVectorToArrey(vec_rebin, n_rebin, rebin);

			TProfile *tp_vn_bin[9];
			TProfile2D *tp_vn[9];
			TProfile *tp_Res = (TProfile*)fileIn -> Get(NameTProfileResolution.Data());
			
			//set v_n axis
			for(Int_t j=CentBinMin;  j<=CentBinMax; j++){
				SetNameTProfile2DFlowEta(j);
			//std::cout<<"\n"<<NameTProfileFlow<<"\n"<<NameTProfileMeanPt<<"\n"<<NameTProfileResolution<<"\n";

				tp_vn[j] = (TProfile2D*)fileIn -> Get(NameTProfileFlow.Data());
				tp_vn_bin[j] = (TProfile*)tp_vn[j] -> ProfileX(Form("profX_%s_cent%i%i", NameTProfileFlow.Data(), j, j),ptBinMin,ptBinMax);
				tp_vn_bin[j] -> Scale( 1.0 / TMath::Sqrt( abs(tp_Res->GetBinContent(tp_Res->FindBin((Double_t)j)))) );
			}

			if(CentBinMin != CentBinMax){
				for(Int_t j=CentBinMin+1; j<=CentBinMax; j++){
					tp_vn_bin[CentBinMin] -> Add(tp_vn_bin[j]);
				}
			}
			TH1 *tp_flow = tp_vn_bin[CentBinMin]->Rebin(n_rebin, NameTProfileFlow.Data(), rebin);
			
			for(Int_t i=0; i < (int)vec_rebin.size()-1; i++){
				eta.push_back(tp_flow->GetBinCenter(i+1));
				eta_err.push_back(0.);
				flow.push_back(tp_flow->GetBinContent(tp_flow->FindBin(eta[i])));
				flow_err.push_back(tp_flow->GetBinError(tp_flow->FindBin(eta[i])));
			}
			
			result->SetVectorsGraph(eta, flow, eta_err, flow_err);
			
			for(Int_t i=CentBinMin; i<=CentBinMax; i++){
				delete tp_vn_bin[i];
				delete tp_vn[i];
			}
			delete tp_Res;

		}

		void FlowVsPtDivideIntegralFLow_EtaSub(Draw_Object *result, Int_t CentBinMin, Int_t CentBinMax,  Vec rebin_vector, Vec rebin_integral , Int_t point_rebin_integral){
			
			Double_t x_int, y_int, ex_int, ey_int;
			Draw_Object *object = new Draw_Object;
			std::vector<Double_t> arrey_sys_err;

			FlowVsPt_EtaSub(result, CentBinMin,CentBinMax,rebin_vector );

			name_sys="Eta15";
			FlowVsPt_EtaSub(object, CentBinMin,CentBinMax,rebin_integral);

			object->GetPointDrawObject(point_rebin_integral, x_int, y_int, ex_int, ey_int);
			
			std::vector<std::vector<Double_t>> ey_int_sys_arrey = object->GetSysArrey();
			std::vector<std::vector<Double_t>> ey_sys_result = result->GetSysArrey();
			
			Double_t mean_res;
			Double_t mean_int;
			Double_t Err_res_int;

			Double_t res_err;
			Double_t res_vol;
			Double_t int_vol;
			Double_t int_err;

			std::cout<<"\n";
			for(Int_t i=0; i < result->GetSizeVector(); i++){
				Double_t x, y, ex, ey;
				result->GetPointDrawObject(i, x, y, ex, ey);
				
				//std::cout<<"v_n: "<<y<<"\tv_n_err: "<< result->GetPointYSysErrorGraph(i)<<"\tv_int: "<<y_int<<"\tv_int_err: "<< object->GetPointYSysErrorGraph(point_rebin_integral)<<"\n";
				//std::cout<<"v_n/v_int: "<<y/y_int<<"\tSys_err:"<<ErrorYRatio(y,y_int,result->GetPointYSysErrorGraph(i),object->GetPointYSysErrorGraph(point_rebin_integral))<<"\n";
				//std::cout<<y<<"\t"<<y_int<<"\t"<<y/y_int<<"\tSys:\t"<< result->GetPointYSysErrorGraph(i)<<"\t"<<object->GetPointYSysErrorGraph(point_rebin_integral)<<"\n";
				
				if(Sys_mod_Integral==true){
				
					mean_res=0.;
					mean_int=0.;
					res_err=0.;
					int_err=0.;
					res_vol=0.;
					int_vol=0.;
					Err_res_int=0.;
					for(Int_t s=0; s<(int)ey_int_sys_arrey.size(); s++){
						std::cout<<ey_int_sys_arrey[s][point_rebin_integral]<<"\n";
					}
					for(Int_t s=0; s<(int)ey_sys_result.size(); s++){
						mean_res=mean_res+ey_sys_result[s][i];
						mean_int=mean_int+ey_int_sys_arrey[s][point_rebin_integral];
					}
					mean_res=mean_res/(int)ey_sys_result.size();
					mean_int=mean_int/(int)ey_sys_result.size();
					
					for(Int_t s=0; s<(int)ey_sys_result.size(); s++){
						Err_res_int = Err_res_int + (ey_sys_result[s][i]-mean_res)*(ey_int_sys_arrey[s][point_rebin_integral]-mean_int);
					}
					Err_res_int = Err_res_int / ey_sys_result.size();
					res_vol = result->GetPointYGraph(i);
					res_err = result->GetPointYSysErrorGraph(i);
					int_vol = object->GetPointYGraph(point_rebin_integral);
					int_err = object->GetPointYSysErrorGraph(point_rebin_integral);
					arrey_sys_err.push_back( sqrt( pow(res_err/int_vol, 2.0) + pow( (int_err*res_vol/pow(mean_int,2.0)) ,2.0) - 2.0*Err_res_int*(res_vol/(pow(mean_int,2.)*int_vol) )) );
	
					//arrey_sys_err.push_back(ErrorYRatio(y,y_int,result->GetPointYSysErrorGraph(i),object->GetPointYSysErrorGraph(point_rebin_integral)));
				}

				result->СhangePointGraph(i, x, y/y_int, ex, ErrorYRatio(y,y_int,ey,ey_int),0.);
				result->ClearSysArrey();
				result->SetYSysErrors(arrey_sys_err);
			}
			std::cout<<"\n";

			/*
			if(Sys_mod_Integral==true){
				for(Int_t i=0; i < (int)Sys_string.at(PID).size(); i++){
					FlowVsPtDivideIntegralFLow_EtaSub_SysError(result, CentBinMin, CentBinMax, rebin_vector, rebin_integral, point_rebin_integral, Sys_string.at(PID)[i]);
				}
				result->SetSysErrors();
			}
			*/
		}

		void FlowVsPtDivideIntegralFLow_EtaSub_SysError(Draw_Object *result, Int_t CentBinMin, Int_t CentBinMax,  Vec rebin_vector, Vec rebin_integral , Int_t point_rebin_integral, TString sys_err){
			
			Double_t x_int, y_int, ex_int, ey_int;
			Draw_Object *object = new Draw_Object;
			Draw_Object *object_res = new Draw_Object;
			
			name_sys = sys_err;
			
			FlowVsPt_EtaSub(object_res, CentBinMin,CentBinMax,rebin_vector );
			FlowVsPt_EtaSub(object, CentBinMin,CentBinMax,rebin_integral);

			object->GetPointDrawObject(point_rebin_integral, x_int, y_int, ex_int, ey_int);
			/*
			Vec gggg;
			gggg = object->GetAxisY();
			for(Int_t i=0; i<gggg.size();i++){
				std::cout<<gggg[i]<<"\t";
			}
			std::cout<<"\n\n";
			*/
			std::vector<Double_t> flow;
			for(Int_t i=0; i < object_res->GetSizeVector(); i++){
				Double_t x, y, ex, ey;
				object_res->GetPointDrawObject(i, x, y, ex, ey);
				object_res->СhangePointGraph(i, x, y/y_int, ex, ErrorYRatio(y,y_int,ey,ey_int),0.);
				object_res->GetPointDrawObject(i, x, y, ex, ey);
				flow.push_back(y);
				//std::cout<<flow.back()<<"\t";
			}
			//std::cout<<"\n\n";
			result->SetSysArrey(flow);
			flow.clear();

		}



		void FlowVsPtDivideIntegralFLow_EtaSub_PID(Draw_Object *result, Int_t CentBinMin, Int_t CentBinMax,  Vec rebin_vector, Vec rebin_integral , Int_t point_rebin_integral){
			
			Double_t x_int, y_int, ex_int, ey_int, ey_int_sys;
			Draw_Object *object_hadrons = new Draw_Object;

			FlowVsPt_EtaSub(result, CentBinMin,CentBinMax,rebin_vector );
			
			PID="Hadrons";
			charge="";
			method_PID="";
			name_sys="Eta15";
			fileIn = new TFile(Form("/home/demanov/New_work/MyAnalysis/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",energy),"READ");
			FlowVsPt_EtaSub(object_hadrons, CentBinMin,CentBinMax,rebin_integral);

			object_hadrons->GetPointDrawObject(point_rebin_integral, x_int, y_int, ex_int, ey_int);
			
			std::vector<std::vector<Double_t>> ey_int_sys_arrey = object_hadrons->GetSysArrey();
			std::vector<std::vector<Double_t>> ey_sys_result = result->GetSysArrey();
			
			Double_t mean_res;
			Double_t mean_int;
			Double_t Err_res_int;

			Double_t res_err;
			Double_t res_vol;
			Double_t int_vol;
			Double_t int_err;

			std::vector<Double_t> arrey_sys_err;

			/*
			Vec gggg;
			gggg = object->GetAxisY();
			for(Int_t i=0; i<gggg.size();i++){
				std::cout<<gggg[i]<<"\t";
			}
			std::cout<<"\n\n";
			*/

			for(Int_t i=0; i < result->GetSizeVector(); i++){
				Double_t x, y, ex, ey;
				result->GetPointDrawObject(i, x, y, ex, ey);

				if(Sys_mod_Integral==true){
					mean_res=0.;
					mean_int=0.;
					res_err=0.;
					int_err=0.;
					res_vol=0.;
					int_vol=0.;
					Err_res_int=0.;
					for(Int_t s=0; s<(int)ey_int_sys_arrey.size(); s++){
						std::cout<<ey_int_sys_arrey[s][point_rebin_integral]<<"\n";
					}
					for(Int_t s=0; s<(int)ey_sys_result.size(); s++){
						mean_res=mean_res+ey_sys_result[s][i];
						mean_int=mean_int+ey_int_sys_arrey[s][point_rebin_integral];
					}
					mean_res=mean_res/(int)ey_sys_result.size();
					mean_int=mean_int/(int)ey_sys_result.size();
					
					for(Int_t s=0; s<(int)ey_sys_result.size(); s++){
						Err_res_int = Err_res_int + (ey_sys_result[s][i]-mean_res)*(ey_int_sys_arrey[s][point_rebin_integral]-mean_int);
					}
					Err_res_int = Err_res_int / ey_sys_result.size();
					//std::cout<<"sigmt\t"<<Err_res_int<<"\n";
					//std::cout<<par->GetPointYGraph(i)<<"\t"<<ant->GetPointYGraph(i)<<"\t"<<r<<"\t"<<par->GetPointYSysErrorGraph(i)<<"\t"<<ant->GetPointYSysErrorGraph(i)<<"\n";
					res_vol = result->GetPointYGraph(i);
					res_err = result->GetPointYSysErrorGraph(i);
					int_vol = object_hadrons->GetPointYGraph(point_rebin_integral);
					int_err = object_hadrons->GetPointYSysErrorGraph(point_rebin_integral);
					//std::cout<<res_vol<<"\t"<<res_err<<"\t"<<int_vol<<"\t"<<int_err<<"\t"<<res_vol/int_vol<<"\t"<<sqrt( pow(res_err/int_vol, 2.0) + pow( (int_err*res_vol/pow(mean_int,2.0)) ,2.0))<<"\n";
					arrey_sys_err.push_back( sqrt( pow(res_err/int_vol, 2.0) + pow( (int_err*res_vol/pow(mean_int,2.0)) ,2.0) - 2.0*Err_res_int*(res_vol/(pow(mean_int,2.)*int_vol) )) );
					//std::cout<<"\t"<<sqrt( pow(res_err/int_vol, 2.0) + pow( (int_err*res_vol/pow(mean_int,2.0)) ,2.0))<<"\n";
					//std::cout<<"\t"<<sqrt( pow(res_err/int_vol, 2.0) + pow( (int_err*res_vol/pow(mean_int,2.0)) ,2.0) - 2.0*Err_res_int*(res_vol/(pow(mean_int,2.)*int_vol) ))<<"\n";
					
					//
				}
				result->СhangePointGraph(i, x, y/y_int, ex, ErrorYRatio(y,y_int,ey,ey_int),0.);
				result->ClearSysArrey();
				result->SetYSysErrors(arrey_sys_err);
				/*
				if(Sys_mod_Integral==true){
					result->SetPointSysErrorGraph(ErrorYRatio(y,y_int,result->GetPointYSysErrorGraph(i),ey_int_sys));
				}
				*/
			}
			/*
			if(Sys_mod_Integral==true){

				for(Int_t i=0; i < (int)Sys_string.at(PID).size(); i++){
					FlowVsPtDivideIntegralFLow_EtaSub_SysError(result, CentBinMin, CentBinMax, rebin_vector, rebin_integral, point_rebin_integral, Sys_string.at(PID)[i]);
				}
				result->SetSysErrors();
			
			}
			*/

		}



		void DifferentParAnt(Draw_Object *result, Int_t CentBinMin, Int_t CentBinMax,  Vec rebin_pt){

			Draw_Object *par = new Draw_Object;
			Draw_Object *ant = new Draw_Object;

			charge="Pos";
			FlowVsPt_EtaSub(par, CentBinMin,CentBinMax,rebin_pt);
			charge="Neg";
			name_sys="Eta01";
			FlowVsPt_EtaSub(ant, CentBinMin,CentBinMax,rebin_pt);

			Double_t r, r_err, r_err_sys;
			std::vector<std::vector<Double_t>> sys_par = par->GetSysArrey();
			std::vector<std::vector<Double_t>> sys_ant = ant->GetSysArrey();
			Double_t mean_par;
			Double_t mean_ant;
			Double_t Err_par_ant;

			for(Int_t i=0; i < par->GetSizeVector(); i++){
				Different(par->GetPointYGraph(i),par->GetPointYErrorGraph(i),ant->GetPointYGraph(i),ant->GetPointYErrorGraph(i),r,r_err);
				result->SetPointGraph(par->GetPointXGraph(i), r, 0., r_err);
				if(Sys_mod_DifParAnr==true){
					mean_par=0.;
					mean_ant=0.;
					Err_par_ant=0.;
					for(Int_t s=0; s<(int)sys_par.size(); s++){
						mean_par=mean_par+sys_par[s][i];
						mean_ant=mean_ant+sys_ant[s][i];
					}
					mean_par=mean_par/(int)sys_par.size();
					mean_ant=mean_ant/(int)sys_par.size();
					
					for(Int_t s=0; s<(int)sys_par.size(); s++){
						Err_par_ant = Err_par_ant + (sys_par[s][i]-mean_par)*(sys_ant[s][i]-mean_ant);
					}
					Err_par_ant = Err_par_ant / sys_par.size();
					result->SetPointSysErrorGraph( sqrt( pow(par->GetPointYSysErrorGraph(i),2.0) + pow(ant->GetPointYSysErrorGraph(i),2.0) - 2.0*Err_par_ant));
				}
			}


			/*

			Double_t r, r_err, r_err_sys;
			for(Int_t i=0; i < par->GetSizeVector(); i++){
				Different(par->GetPointYGraph(i),par->GetPointYErrorGraph(i),ant->GetPointYGraph(i),ant->GetPointYErrorGraph(i),r,r_err);
				//std::cout<<par->GetPointYGraph(i)<<"\t\t"<<par->GetPointYErrorGraph(i)<<"\t\t"<<par->GetPointYSysErrorGraph(i)<<"\n";
				//std::cout<<ant->GetPointYGraph(i)<<"\t\t"<<ant->GetPointYErrorGraph(i)<<"\t\t"<<ant->GetPointYSysErrorGraph(i)<<"\n";
				//result->SetPointGraph(par->GetPointXGraph(i), r, 0., r_err);
				//result->SetPointSysErrorGraph( sqrt( pow(par->GetPointYSysErrorGraph(i),2.0) + pow(ant->GetPointYSysErrorGraph(i),2.0)));
			}
			*/
			/*
			std::cout<<"stop\n";
			r_err=0.;
			if(Sys_mod_DifParAnr==true){
				for(Int_t i=0; i<par->GetSizeVector(); i++){
					//Different(par->GetPointYGraph(i),par->GetPointYSysErrorGraph(i),ant->GetPointYGraph(i),ant->GetPointYSysErrorGraph(i),r,r_err);
					//result->SetPointSysErrorGraph(r_err);
					//std::cout<<"\t"<<par->GetPointYSysErrorGraph(i)/par->GetPointYGraph(i)<<"\t"<<ant->GetPointYSysErrorGraph(i)/ant->GetPointYGraph(i)<<"\t"<<r_err<<"\n";
					//std::cout<<result->GetSizeVectorSys()<<"\t"<<result->GetPointYSysErrorGraph(i)<<"\n";
				}
				
				//
				for(Int_t i=0; i < (int)Sys_string.at(PID).size(); i++){
					DifferentParAnt_SysErrors(result, CentBinMin, CentBinMax, rebin_pt, Sys_string.at(PID)[i]);
				}
				result->SetSysErrors();
				//
			}
			*/

			par->ClearObject();
			ant->ClearObject();
		}

		void DifferentParAnt_SysErrors(Draw_Object *result, Int_t CentBinMin, Int_t CentBinMax,  Vec rebin_pt, TString sys_err){

			Draw_Object *par = new Draw_Object;
			Draw_Object *ant = new Draw_Object;

			name_sys = sys_err;

			charge="Pos";
			//name_sys="Eta01";
			FlowVsPt_EtaSub(par, CentBinMin,CentBinMax,rebin_pt);
			charge="Neg";
			//name_sys="Eta01";
			FlowVsPt_EtaSub(ant, CentBinMin,CentBinMax,rebin_pt);

			Double_t r, r_err;
			std::vector<Double_t> flow;
			for(Int_t i=0; i < par->GetSizeVector(); i++){
				Different(par->GetPointYGraph(i),par->GetPointYErrorGraph(i),ant->GetPointYGraph(i),ant->GetPointYErrorGraph(i),r,r_err);
				flow.push_back(r);
			}
			result->SetSysArrey(flow);
		}

		void Different(Double_t a, Double_t a_err , Double_t b, Double_t b_err, Double_t &result, Double_t &result_err){
			result = a - b;
			result_err = sqrt( a_err*a_err + b_err*b_err );
		}

		void RatioParAnt(Draw_Object *result, Int_t CentBinMin, Int_t CentBinMax,  Vec rebin_pt){

			Draw_Object *par = new Draw_Object;
			Draw_Object *ant = new Draw_Object;

			std::cout<<"StartRatio\n";

			charge="Pos";
			FlowVsPt_EtaSub(par, CentBinMin,CentBinMax,rebin_pt);
			charge="Neg";
			name_sys="Eta01";
			FlowVsPt_EtaSub(ant, CentBinMin,CentBinMax,rebin_pt);

			Double_t r, r_err, r_err_sys;
			std::vector<std::vector<Double_t>> sys_par = par->GetSysArrey();
			std::vector<std::vector<Double_t>> sys_ant = ant->GetSysArrey();
			Double_t mean_par;
			Double_t mean_ant;
			Double_t Err_par_ant;

			Double_t par_err;
			Double_t par_vol;
			Double_t ant_vol;
			Double_t ant_err;

			/*
			for(Int_t i=0; i < par->GetSizeVector(); i++){
				Ratio(par->GetPointYGraph(i),par->GetPointYErrorGraph(i),ant->GetPointYGraph(i),ant->GetPointYErrorGraph(i),r,r_err);
				result->SetPointGraph(par->GetPointXGraph(i), r, 0., r_err);
			}
			*/
			
			for(Int_t i=0; i < par->GetSizeVector(); i++){
				Ratio(par->GetPointYGraph(i),par->GetPointYErrorGraph(i),ant->GetPointYGraph(i),ant->GetPointYErrorGraph(i),r,r_err);
				result->SetPointGraph(par->GetPointXGraph(i), r, 0., r_err);
				if(Sys_mod_RatParAnr==true){
					mean_par=0.;
					mean_ant=0.;
					par_err=0.;
					ant_err=0.;
					par_vol=0.;
					ant_vol=0.;
					Err_par_ant=0.;
					for(Int_t s=0; s<(int)sys_par.size(); s++){
						mean_par=mean_par+sys_par[s][i];
						mean_ant=mean_ant+sys_ant[s][i];
					}
					mean_par=mean_par/(int)sys_par.size();
					mean_ant=mean_ant/(int)sys_par.size();
					
					for(Int_t s=0; s<(int)sys_par.size(); s++){
						Err_par_ant = Err_par_ant + (sys_par[s][i]-mean_par)*(sys_ant[s][i]-mean_ant);
					}
					Err_par_ant = Err_par_ant / sys_par.size();
					std::cout<<Err_par_ant<<"\n";
					std::cout<<par->GetPointYGraph(i)<<"\t"<<ant->GetPointYGraph(i)<<"\t"<<r<<"\t"<<par->GetPointYSysErrorGraph(i)<<"\t"<<ant->GetPointYSysErrorGraph(i)<<"\n";
					std::cout<<Err_par_ant<<"\n";
					std::cout<<r<<"\t"<<ErrorYRatio( par->GetPointYGraph(i), ant->GetPointYGraph(i), par->GetPointYSysErrorGraph(i), ant->GetPointYSysErrorGraph(i))<<"\n";
					std::cout<<r<<"\t"<<sqrt( pow(par->GetPointYSysErrorGraph(i)/ant->GetPointYGraph(i), 2.0) + pow( (ant->GetPointYSysErrorGraph(i)* par->GetPointYGraph(i)/pow(ant->GetPointYGraph(i),2.0)) ,2.0))<<"\n";
					//result->SetPointGraph(par->GetPointXGraph(i), r, 0., r_err);
					par_vol = par->GetPointYGraph(i);
					par_err = par->GetPointYSysErrorGraph(i);
					ant_vol = ant->GetPointYGraph(i);
					ant_err = ant->GetPointYSysErrorGraph(i);
					//std::cout<<par_vol<<"\t"<<par_err<<"\t"<<ant_vol<<"\t"<<ant_err<<"\t"<<sqrt( pow(par_err/ant_vol, 2.0) + pow( (ant_err*par_vol/pow(ant_vol,2.0)) ,2.0))<<"\n";
					result->SetPointSysErrorGraph( sqrt( pow(par_err/ant_vol, 2.0) + pow( (ant_err*par_vol/pow(mean_ant,2.0)) ,2.0) - 2.0*Err_par_ant*(par_vol/(pow(mean_ant,2.)*ant_vol) )) );
					//std::cout<<r<<"\t"<<sqrt( pow(par_err/ant_vol, 2.0) + pow( (ant_err*par_vol/pow(mean_ant,2.0)) ,2.0))<<"\n";
					//std::cout<<r<<"\t"<<           sqrt( pow(par_err/ant_vol, 2.0) + pow( (ant_err*par_vol/pow(mean_ant,2.0)) ,2.0) - 2.0*Err_par_ant*(par_vol/(pow(mean_ant,2.)*ant_vol) ))<<"\n";
					
					//
				}
			}
			
			par->ClearObject();
			ant->ClearObject();
			std::cout<<"endRatio\n";
			/*
			r_err=0.;
			
			if(Sys_mod_RatParAnr==true){
				
				for(Int_t i=0; i<par->GetSizeVector(); i++){
					Ratio(par->GetPointYGraph(i),par->GetPointYSysErrorGraph(i),ant->GetPointYGraph(i),ant->GetPointYSysErrorGraph(i),r,r_err);
					result->SetPointSysErrorGraph(r_err);
				}
				
				
				for(Int_t i=0; i < (int)Sys_string.at(PID).size(); i++){
					RatioParAnt_SysErrors(result, CentBinMin, CentBinMax, rebin_pt, Sys_string.at(PID)[i]);
				}
				result->SetSysErrors();
				
			}
			*/
			
		}
		
		void RatioParAnt_SysErrors(Draw_Object *result, Int_t CentBinMin, Int_t CentBinMax,  Vec rebin_pt, TString sys_err){

			Draw_Object *par = new Draw_Object;
			Draw_Object *ant = new Draw_Object;

			name_sys = sys_err;

			charge="Pos";
			FlowVsPt_EtaSub(par, CentBinMin,CentBinMax,rebin_pt);
			charge="Neg";
			FlowVsPt_EtaSub(ant, CentBinMin,CentBinMax,rebin_pt);

			Double_t r, r_err;
			std::vector<Double_t > flow;
			for(Int_t i=0; i < par->GetSizeVector(); i++){
				Ratio(par->GetPointYGraph(i),par->GetPointYErrorGraph(i),ant->GetPointYGraph(i),ant->GetPointYErrorGraph(i),r,r_err);
				flow.push_back(r);
			}
			result->SetSysArrey(flow);
		}

		void Ratio(Double_t a, Double_t a_err , Double_t b, Double_t b_err, Double_t &result, Double_t &result_err){
			result = a/b;
			result_err = ErrorYRatio(a, b, a_err, b_err);
		}
		

		Double_t ErrorYRatio(Double_t num, Double_t den, Double_t num_err, Double_t den_err){

			return sqrt( pow( num_err/den, 2) + pow( (num*den_err) / (den*den) , 2) );
		
		}

		void RatioGraphToOnePoint(Draw_Object *result, Double_t den, Double_t den_err){
			
			Double_t x,y,ex,ey;
			
			for(Int_t i=0; i < result->GetSizeVector(); i++){
				result->GetPointDrawObject(i,x,y,ex,ey);
				result->СhangePointGraph(i,x,y/den,ex, ErrorYRatio(y,den,ey,den_err),0); ///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			}
		}

		void ConvertVectorToArrey(Vec vector, Int_t &n_arrey, Double_t array[]){
			
			for(Int_t i=0; i < (int)vector.size(); i++){
				array[i] = vector[i];
			}
			n_arrey = (int)vector.size()-1;
			
			if( (vector[0]==0.2 || vector[0]==0.15) && vector[1]==0.5 ){
				for(Int_t i=1; i < (int)vector.size(); i++){
					array[i-1] = vector[i];
				}
				n_arrey = (int)vector.size()-2;
				array[(int)vector.size()]=0.0;
			}
			
		}

		void RatioGraphPointToPoint(Draw_Object *result, Draw_Object *numerator, Draw_Object *denominator, TString legend){

			Double_t num_x;
			Double_t num_y;
			Double_t num_x_err;
			Double_t num_y_err;

			Double_t den_x;
			Double_t den_y;
			Double_t den_x_err;
			Double_t den_y_err;

			TString Name ="Ratio_" + numerator -> GetNameGraph()+ "_" + denominator -> GetNameGraph();

			for(Int_t i = 0; i < numerator->GetSizeVector(); i++){
				
				numerator -> GetPointDrawObject(i,num_x,num_y,num_x_err,num_y_err);
				denominator -> GetPointDrawObject(i,den_x,den_y, den_x_err, den_y_err);

				result->SetPointGraph(num_x, num_y / den_y, 0., ErrorYRatio(num_y,den_y,num_y_err,den_y_err) );
				//std::cout<<ErrorYRatio(num_y,den_y,num_y_err,den_y_err)<<"\t"<<den_y_err<<"\t"<<num_y_err<<"\n";

			}
			//std::cout<<"\n";
			result->SetParametrsGraph(Name, legend, numerator->GetColor(), numerator->GetMarker(), numerator->GetMarkerSize());

		}


		void RatioGraphEVAL(Draw_Object *result, Draw_Object *numerator, Draw_Object *denominator,TString legend){

			std::vector<Double_t> x;
			std::vector<Double_t> y;
			std::vector<Double_t> x_err;
			std::vector<Double_t> y_err;

			Double_t num_x;
			Double_t num_y;
			Double_t num_y_err;

			Double_t den_x;
			Double_t den_y;
			Double_t den_y_err;

			TGraphErrors *graph_numerator = numerator->GetTGraph();
			TGraphErrors *graph_denominator = denominator->GetTGraph();

			for(Int_t i = 0; i < (int)graph_numerator->GetN(); i++){
				
				graph_numerator -> GetPoint(i,num_x,num_y);
				graph_denominator -> GetPoint(i,den_x,den_y);
				
				//std::cout<<den_x<<"\t"<<den_y<<"\t"<<num_x<<"\t"<<num_y<<"\n";
				
				num_y_err = graph_numerator->GetErrorY(i);
				//den_y_err = graph_denominator->GetErrorY(i);

				x.push_back( num_x );
				y.push_back( abs(num_y / graph_denominator->Eval(x[i],0,"S")) );
				x_err.push_back( 0.0 );

				if((int)graph_denominator->GetN()!=(int)graph_numerator->GetN()){
					for(Int_t j=0; j<(int)graph_denominator->GetN()-1; j++){
						
						Double_t x0,y0;
						Double_t x1,y1;
						graph_denominator -> GetPoint(j,x0,y0);
						graph_denominator -> GetPoint(j+1,x1,y1);
						
						if(j==0 && x[i] < x0){
							den_y_err = graph_denominator->Eval(x[i],0,"S") * graph_denominator->GetErrorY(j) / y0 ;
							j = (int)graph_denominator->GetN();
						}
						if(x[i]>x0 && x[i]<x1){
							den_y_err = (graph_denominator->GetErrorY(j) + graph_denominator->GetErrorY(j+1))/2.0;
							j = (int)graph_denominator->GetN();
						}
						if(j==(int)graph_denominator->GetN() && x[i]>x1){
							den_y_err = graph_denominator->Eval(x[i],0,"S") * graph_denominator->GetErrorY(j+1) / y1 ;
							j = (int)graph_denominator->GetN();
						}	
					}
				}else{
					den_y_err = graph_denominator->GetErrorY(i);
				}
				den_y = graph_denominator->Eval(x[i],0,"S");
				y_err.push_back( ErrorYRatio(num_y,den_y,num_y_err,den_y_err) );
				//std::cout<<ErrorYRatio(num_y,den_y,num_y_err,den_y_err)<<"\t"<<den_y_err<<"\t"<<num_y_err<<"\n";

			}
			//std::cout<<"\n";

			result->SetVectorsGraph(x,y,x_err,y_err);
			result->SetParametrsGraph(Form("Ratio_%s_%s", graph_numerator -> GetName(), graph_denominator -> GetName()),legend, graph_numerator->GetMarkerColor(), graph_numerator->GetMarkerStyle(), graph_numerator->GetMarkerSize());

		}



		void DifferentGraphEVAL(Draw_Object *result, Draw_Object *numerator, Draw_Object *denominator,TString legend){

			std::vector<Double_t> x;
			std::vector<Double_t> y;
			std::vector<Double_t> x_err;
			std::vector<Double_t> y_err;

			Double_t num_x;
			Double_t num_y;
			Double_t num_y_err;

			Double_t den_x;
			Double_t den_y;
			Double_t den_y_err;

			TGraphErrors *graph_numerator = numerator->GetTGraph();
			TGraphErrors *graph_denominator = denominator->GetTGraph();

			for(Int_t i = 0; i < (int)graph_numerator->GetN(); i++){
				
				graph_numerator -> GetPoint(i,num_x,num_y);
				graph_denominator -> GetPoint(i,den_x,den_y);
				
				num_y_err = graph_numerator->GetErrorY(i);
				//den_y_err = graph_denominator->GetErrorY(i);

				x.push_back( num_x );
				//y.push_back( abs(num_y / graph_denominator->Eval(x[i],0,"S")) );
				y.push_back( num_y - graph_denominator->Eval(x[i],0,"S"));
				x_err.push_back( 0.0 );

				if((int)graph_denominator->GetN()!=(int)graph_numerator->GetN()){
					for(Int_t j=0; j<(int)graph_denominator->GetN()-1; j++){
						
						Double_t x0,y0;
						Double_t x1,y1;
						graph_denominator -> GetPoint(j,x0,y0);
						graph_denominator -> GetPoint(j+1,x1,y1);
						
						if(j==0 && x[i] < x0){
							den_y_err = graph_denominator->Eval(x[i],0,"S") * graph_denominator->GetErrorY(j) / y0 ;
							j = (int)graph_denominator->GetN();
						}
						if(x[i]>x0 && x[i]<x1){
							den_y_err = (graph_denominator->GetErrorY(j) + graph_denominator->GetErrorY(j+1))/2.0;
							j = (int)graph_denominator->GetN();
						}
						if(j==(int)graph_denominator->GetN() && x[i]>x1){
							den_y_err = graph_denominator->Eval(x[i],0,"S") * graph_denominator->GetErrorY(j+1) / y1 ;
							j = (int)graph_denominator->GetN();
						}	
					}
				}else{
					den_y_err = graph_denominator->GetErrorY(i);
				}
				den_y = graph_denominator->Eval(x[i],0,"S");
				//y_err.push_back( ErrorYRatio(num_y,den_y,num_y_err,den_y_err) );
				y_err.push_back( sqrt( pow(den_y_err ,2.) + pow(num_y_err ,2.) ) );
				//std::cout<<ErrorYRatio(num_y,den_y,num_y_err,den_y_err)<<"\t"<<den_y_err<<"\t"<<num_y_err<<"\n";


			}
			//std::cout<<"\n";

			result->SetVectorsGraph(x,y,x_err,y_err);
			result->SetParametrsGraph(Form("Ratio_%s_%s", graph_numerator -> GetName(), graph_denominator -> GetName()),legend, graph_numerator->GetMarkerColor(), graph_numerator->GetMarkerStyle(), graph_numerator->GetMarkerSize());

		}



	
};