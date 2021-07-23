#include <map>

class Draw_Picture_new {
    public:
        // Деструктор
        ~Draw_Picture_new();
        
        struct TextInfo
        {
        	TString text;
        	TString place; //X , Y , Pad, Info
        	Float_t x;
        	Float_t y;
        	Float_t size;
        	Float_t angle;
        };

        std::vector<TextInfo> labels;
        std::vector<TextInfo> labels_b;

        Bool_t Sys_mod=false;
        Bool_t legend_flag=false;

        TFile *fileIn;
        TFile *fileOut;

        TCanvas *canvas;

        TH2F *hPadUp;
        TH2F *hPadDown;
        TH2F *hPadOne;
        std::vector<TH2F *> hPadNM;

        TLegend *legend;

		TLatex _textInfo;

		Int_t type_canvas;

		std::vector<TLine *> Line_vector;
		std::vector<TString> Line_vector_DrawOptions;

		std::vector<TGraphErrors *> graph_up;
		std::vector<TGraphErrors *> graph_up_sys;
		std::vector<TGraphErrors *> graph_down;

		void SetTLine(Float_t mXmin, Float_t mYmin, Float_t mXmax, Float_t mYmax, Int_t mColour, Int_t mLineWidth, Int_t mLineStyle, TString options) {

			Line_vector.push_back(new TLine( mXmin, mYmin, mXmax, mYmax ));
			Line_vector_DrawOptions.push_back(options);
			Line_vector[Line_vector.size()-1] -> SetLineColor( mColour );
			Line_vector[Line_vector.size()-1] -> SetLineWidth( mLineWidth );
			Line_vector[Line_vector.size()-1] -> SetLineStyle( mLineStyle );

		}// 

		void SetAxisToCanvsWithBottomPanel(Double_t Xmin, Double_t Xmax,Double_t YminUp, Double_t YmaxUp, Double_t YminDown, Double_t YmaxDown){

			hPadUp = new TH2F(Form("hPadUp%.4f", YmaxUp),"", 2,Xmin,Xmax,2,YminUp,YmaxUp);
            hPadUp->GetYaxis()->SetLabelSize(0.08);

			hPadDown = new TH2F(Form("hPadDown%.4f", YmaxUp),"", 2,Xmin,Xmax,2,YminDown,YmaxDown);
            hPadDown->GetXaxis()->SetLabelSize(0.18);
            hPadDown->GetYaxis()->SetLabelSize(0.17);

		}


		void SetLegend(Double_t xMin, Double_t yMin, Double_t xMax, Double_t yMax, Int_t N){
			



		TGraphErrors *grl = new TGraphErrors;
        grl->SetMarkerColor(1);
		grl->SetMarkerStyle(25);
		grl->SetMarkerSize(1.4);
		TGraphErrors *grlpubl = new TGraphErrors;
		grlpubl->SetMarkerStyle(21);
        grlpubl->SetMarkerColor(1);
        grlpubl->SetMarkerSize(1.4);



			legend = new TLegend(xMin,yMin,xMax,yMax);
			legend->SetNColumns(N);
			
			for(Int_t i = 0; i < graph_up.size(); i++){
				legend->AddEntry(graph_up[i],graph_up[i]->GetTitle(),"p");
			}
			
			//	legend->AddEntry(grlpubl,"#font[42]{#scale[1.4]{STAR, Phys.Rev.C 93 (2016)}}","p");
			//	legend->AddEntry(grl,"#font[42]{#scale[1.4]{vHLLE+UrQMD}}","p");

			legend_flag=true;
			
		}

		void SetAxisToCanvsOne( Double_t Xmin, Double_t Xmax ,Double_t Ymin, Double_t Ymax){

			hPadOne = new TH2F("hPadOne","", 2,Xmin,Xmax,2,Ymin,Ymax);
            hPadOne->GetYaxis()->SetLabelSize(0.05);

		}

		void SetAxisToCanvsNM( Double_t Xmin, Double_t Xmax ,Double_t Ymin, Double_t Ymax, Double_t size){

			hPadNM.push_back( new TH2F("","", 2,Xmin,Xmax,2,Ymin,Ymax));
            hPadNM.back()->GetYaxis()->SetLabelSize(size);
            hPadNM.back()->GetXaxis()->SetLabelSize(size);

		}

		void SetTextInfo(TString text, TString place, Float_t x, Float_t y, Float_t size, Float_t angle){
			labels.push_back({text,place,x,y,size, angle});
		}

		void SetNumberPad(TString text, TString place, Float_t x, Float_t y, Float_t size, Float_t angle){
			labels_b.push_back({text,place,x,y,size, angle});
		}

		void SetDrawObject( std::vector<Draw_Object *> object_vec, TString options){
			for(Int_t i=0; i < object_vec.size(); i++){
				if(strncmp(options, "Up",2)==0 || strncmp(options, "OnePad",6)==0){
					graph_up.push_back(object_vec[i]->GetTGraph());
					if(Sys_mod==true) graph_up_sys.push_back(object_vec[i]->GetTGraphSys());
				}
				if(strncmp(options, "Down",4)==0){
					graph_down.push_back(object_vec[i]->GetTGraph());
				}
			}
		}

		TCanvas *CanvasWithBottomPanel(TString text, Double_t xText, Double_t yText, TString Xname, TString YnameUp, TString YnameDown, TString PhysRef){

			canvas = new TCanvas("canvas","plot",1080,960);
			//canvas = new TCanvas("canvas","plot",1280,960);
            canvas->GetFrame()->SetFillColor(21);
            canvas->GetFrame()->SetBorderSize(115);
            canvas->cd();
            
            TPad *pad_up = new TPad("pad_up", "", 0.1, 0.38, 1.0, 1.0, 0, 0, 0);
            TPad *pad_down = new TPad("pad_down", "", 0.1, 0.06, 1.0, 0.38, 0, 0, 0);
            TPad *pad_AxisY = new TPad("pad_AxisY", "", 0.0, 0.06, 0.1, 1.0, 0, 0, 0);
			TPad *pad_AxisX = new TPad("pad_AxisX", "", 0.0, 0.0, 1.0, 0.06, 0, 0, 0);

            pad_up->SetFrameBorderMode(0);
            pad_up->SetBottomMargin(0);
            pad_up->SetLeftMargin(0.09);
            
            pad_down->SetFrameBorderMode(0);
            pad_down->SetTopMargin(0.005);
            pad_down->SetBottomMargin(0.26);
            pad_down->SetLeftMargin(0.09);
            
            pad_AxisY->SetFrameBorderMode(0);
            pad_AxisY->SetRightMargin(0.23);
            pad_AxisY->SetBottomMargin(0);
            
            pad_AxisX->SetTopMargin(0);
            pad_AxisX->SetFrameBorderMode(0);
            pad_AxisX->SetBottomMargin(0);

            pad_up->Draw();
            pad_down->Draw();
            pad_AxisY->Draw();
            pad_AxisX->Draw();

			pad_up->cd();
			pad_up->Draw();
			hPadUp->Draw();
			
			_textInfo.SetTextSize(0.07);   	
			_textInfo.DrawLatex(xText, yText, text.Data() );	

			//_textInfo.SetTextSize(0.07);   	
			//_textInfo.DrawLatex(1.8, 0.005, "#font[42]{ #scale[1.0]{STAR Preliminary}}" ); 
			
			legend->Draw("same");
			for(Int_t i = 0; i < graph_up.size(); i++){
				graph_up[i]->Draw("Psame");
			}
			for(Int_t i = 0; i < Line_vector.size(); i++){
				if(strncmp(Line_vector_DrawOptions[i], "Up",2)==0){
					Line_vector[i]->Draw("same");
				}
			}
 	
			pad_down-> cd();
			hPadDown->Draw();

			for(Int_t i = 0; i < graph_down.size(); i++){
				graph_down[i]->Draw("Psame");
			}
			for(Int_t i = 0; i < Line_vector.size(); i++){
				if(strncmp(Line_vector_DrawOptions[i], "Down",4)==0){
					Line_vector[i]->Draw("same");
				}
			}

			pad_AxisY->cd();
			TLatex textY;
			textY.SetTextSize(0.4);
			textY.SetTextAngle(90);
			textY.DrawLatex(0.5,0.55,YnameUp.Data());
			//textY.DrawLatex(0.3,0.85,YnameUp.Data());

			textY.SetTextSize(0.4);
			textY.SetTextAngle(90);
			textY.DrawLatex(0.6,0.09,YnameDown.Data() );

			pad_AxisX->cd();
			TLatex textArticle;
			textArticle.SetTextSize(0.8);
			textArticle.DrawLatex(0.75,0.35,Xname.Data());
			textArticle.DrawLatex(0.05,0.4,PhysRef.Data());
		
			return canvas;
		
		}


		TCanvas *CanvasOne(TString text, Double_t xText, Double_t yText, TString Xname, TString YnameUp, TString PhysRef){

			canvas = new TCanvas("canvas","plot",1280,960);
            canvas->GetFrame()->SetFillColor(21);
            canvas->GetFrame()->SetBorderSize(115);
            canvas->cd();
            std::cout<<"good\n";
	
            TPad *pad_up = new TPad("pad_up", "", 0.1, 0.06, 1.0, 1.0, 0, 0, 0);
            TPad *pad_AxisY = new TPad("pad_AxisY", "", 0.0, 0.06, 0.1, 1.0, 0, 0, 0);
			TPad *pad_AxisX = new TPad("pad_AxisX", "", 0.0, 0.0, 1.0, 0.06, 0, 0, 0);

            pad_up->SetFrameBorderMode(0);
            pad_up->SetBottomMargin(0.1);
            pad_up->SetLeftMargin(0.09);
            pad_up->SetTopMargin(0.05);

            
            pad_AxisY->SetFrameBorderMode(0);
            pad_AxisY->SetRightMargin(0.23);
            pad_AxisY->SetBottomMargin(0);
            
            pad_AxisX->SetTopMargin(0);
            pad_AxisX->SetFrameBorderMode(0);
            pad_AxisX->SetBottomMargin(0);

            pad_up->Draw();
            pad_AxisY->Draw();
            pad_AxisX->Draw();

            pad_up->cd();
			hPadOne->Draw();

			_textInfo.SetTextSize(0.07);   	
			_textInfo.DrawLatex(xText, yText, text.Data() );	 

			legend->Draw("same");
			for(Int_t i = 0; i < graph_up.size(); i++){
				graph_up[i]->Draw("Psame");
			}
			for(Int_t i = 0; i < Line_vector.size(); i++){
				if(strncmp(Line_vector_DrawOptions[i], "Up",2)==0){
					Line_vector[i]->Draw("same");
				}
			}

			pad_AxisY->cd();
			TLatex textY;
			textY.SetTextSize(0.5);
			//textY.SetTextAngle(90);
			textY.DrawLatex(0.3,0.85,YnameUp.Data());

			pad_AxisX->cd();
			TLatex textArticle;
			textArticle.SetTextSize(1.0);
			textArticle.DrawLatex(0.75,0.4,Xname.Data());
			textArticle.DrawLatex(0.05,0.4,PhysRef.Data());

			return canvas;
		
		}

		/*
		TCanvas *CanvasNxM(Int_t xCanvas, Int_t yCanvas, Int_t line, Int_t column, Float_t left, Float_t down, std::vector<std::vector<Draw_Object*>> DrawResult){

			//Float_t
			Float_t otstupUp=0.02;
			Float_t otstupDown=0.1;
			Float_t otstupLeft=0.1;
			Float_t otstupRight=0.02;
			Float_t h = (1.0 - otstupUp - otstupDown )/line;
			Float_t l = (1.0 - otstupRight - otstupLeft )/column;

			canvas = new TCanvas("canvas","plot",xCanvas,yCanvas);
            canvas->GetFrame()->SetFillColor(21);
            canvas->GetFrame()->SetBorderSize(115);
            canvas->cd();
            
            std::vector<TPad *> pad;
            std::vector<TPad *> pad_axisDown;
            std::vector<TPad *> pad_axisLeft;

            for(Int_t y=line; y>0; y--){
            	for(Int_t x=0; x<column; x++){
            		//std::cout<<0.1+l*x<<"\t"<< 1.0-h*y<<"\t"<< 0.1+l*(x+1)<<"\t"<< 1.0-h*(y-1)<<"\n\n";
	            	pad.push_back(new TPad(Form("Pad_%i_%i",x,y),"",(Float_t)(otstupLeft + l*x),     (Float_t)( otstupDown + h*(y-1)), 
	            													(Float_t)(otstupLeft + l*(x+1)), (Float_t)( otstupDown + h*y), 0,0,0));
	            	pad.back()->SetFrameBorderMode(1);
	            	pad.back()->SetLeftMargin(0);
	            	pad.back()->SetRightMargin(0);
	            	pad.back()->SetBottomMargin(0);
	            	pad.back()->SetTopMargin(0);
		            
		            if(x==0){
						pad_axisLeft.push_back(new TPad(Form("PadAx_%i_%i",x,y),"",(Float_t)( l*x),     (Float_t)( otstupDown + h*(y-1)), 
	            													               (Float_t)( l*(x+1)), (Float_t)( otstupDown + h*y), 0,0,0));
						pad_axisLeft.back()->SetFrameBorderMode(1);
	            		pad_axisLeft.back()->SetLeftMargin(left);
	            		pad_axisLeft.back()->SetRightMargin(0);
	            		pad_axisLeft.back()->SetBottomMargin(0);
	            		pad_axisLeft.back()->SetTopMargin(0);
		            }
		            if(y==1){
            			//std::cout<<(Float_t)(otstupLeft + l*x)<<"\t"<< (Float_t)( 1.0 - otstupUp - otstupDown - h*y)<<"\t"<< (Float_t)(otstupLeft + l*(x+1))<<"\t"<< (Float_t)( 1.0 - otstupUp - otstupDown - h*(y-1))<<"\n\n";
						pad_axisDown.push_back(new TPad(Form("Pad_%i_%i",x,y),"",(Float_t)(otstupLeft + l*x),     (Float_t)( h*(y-1)), 
	            															     (Float_t)(otstupLeft + l*(x+1)), (Float_t)( h*y), 0,0,0));
						pad_axisDown.back()->SetFrameBorderMode(1);
	            		pad_axisDown.back()->SetLeftMargin(0.0);
	            		pad_axisDown.back()->SetRightMargin(0);
	            		pad_axisDown.back()->SetBottomMargin(down);
	            		pad_axisDown.back()->SetTopMargin(0);
		            }

            	}
            	
            }

            TPad *pad_AxisY = new TPad("pad_AxisY", "", 0.0, 0.05, otstupLeft, 1.0, 0, 0, 0);
			TPad *pad_AxisX = new TPad("pad_AxisX", "", 0.0, 0.0, 1.0, otstupLeft, 0, 0, 0);

			pad_AxisY->SetFrameBorderMode(0);
            pad_AxisY->SetRightMargin(0.23);
            pad_AxisY->SetBottomMargin(0);
            
            pad_AxisX->SetTopMargin(0);
            pad_AxisX->SetFrameBorderMode(0);
            pad_AxisX->SetBottomMargin(0);

            pad_AxisY->Draw();
            pad_AxisX->Draw();

			for(Int_t i=0; i<pad_axisDown.size(); i++){
            	pad_axisDown[i]->Draw();
            }
            for(Int_t i=0; i<pad_axisLeft.size(); i++){
            	pad_axisLeft[i]->Draw();
            }
            for(Int_t i=0; i<pad.size(); i++){
            	pad[i]->Draw();
            }

            for(Int_t i=0; i<pad_axisDown.size(); i++){
            	pad_axisDown[i]->cd();
            	hPadOne->Draw();
            }
            for(Int_t i=0; i<pad_axisLeft.size(); i++){
            	pad_axisLeft[i]->cd();
            	hPadOne->Draw();
            }
            for(Int_t i=0; i<pad.size(); i++){
            	pad[i]->cd();
            	hPadOne->Draw();
            }

            for(Int_t i=0; i<pad.size(); i++){
            	pad[i]->cd();
            	
            }
            
            pad[0]->cd();
            legend->Draw("same");
            
            for(Int_t i=0; i<DrawResult.size(); i++){
            	pad[i]->cd();
            	_textInfo.SetTextSize(labels[i].size);
            	_textInfo.SetTextAngle(labels[i].angle);
				_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );

            	SetDrawObject(DrawResult[i],"OnePad");
            	for(Int_t j=0; j<graph_up.size();j++){
            		graph_up[j]->Draw("Psame");
            	}
            	for(Int_t l = 0; l < Line_vector.size(); l++){
					if(strncmp(Line_vector_DrawOptions[l], "Up",2)==0){
						Line_vector[l]->Draw("same");
					}
				}
            	graph_up.clear();
            }

            pad_AxisY->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "Y",1)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
				}
			}

			pad_AxisX->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "X",1)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
				}
			}

			pad[1]->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "Info",4)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
					std::cout<<"Good\n";
				}
			}

            return canvas;

		}
		*/

		TCanvas *CanvasNxM(Int_t xCanvas, Int_t yCanvas, Int_t line, Int_t column, Float_t left, Float_t down, std::vector<std::vector<Draw_Object*>> DrawResult, Int_t fiew, Int_t plInfo, Int_t plLeg){

			//Float_t
			Float_t otstupUp=0.02;
			//Float_t otstupDown=0.1;
			Float_t otstupDown=0.15;
			//Float_t otstupLeft=0.13;
			//Float_t otstupLeft=0.16;
			Float_t otstupLeft=0.15;
			Float_t otstupRight=0.1;
			Float_t h = (1.0 - otstupUp - otstupDown )/line;
			Float_t l = (1.0 - otstupRight - otstupLeft )/column;

			canvas = new TCanvas("canvas","plot",xCanvas,yCanvas);
            canvas->GetFrame()->SetFillColor(21);
            canvas->GetFrame()->SetBorderSize(115);
            canvas->cd();
            
            std::vector<TPad *> pad;
            std::vector<TPad *> pad_axisDown;
            std::vector<TPad *> pad_axisLeft;
            Int_t n=0;

            for(Int_t y=line; y>0; y--){
            	for(Int_t x=0; x<column; x++){
            		//std::cout<<0.1+l*x<<"\t"<< 1.0-h*y<<"\t"<< 0.1+l*(x+1)<<"\t"<< 1.0-h*(y-1)<<"\n\n";
	            	pad.push_back(new TPad(Form("Pad_%i_%i",x,y),"",(Float_t)(otstupLeft + l*x),     (Float_t)( otstupDown + h*(y-1)), 
	            													(Float_t)(otstupLeft + l*(x+1)), (Float_t)( otstupDown + h*y), 0,0,0));
	            	pad.back()->SetFrameBorderMode(1);
	            	pad.back()->SetLeftMargin(0);
	            	pad.back()->SetRightMargin(0);
	            	pad.back()->SetBottomMargin(0);
	            	pad.back()->SetTopMargin(0);
		            
		            if(x==0){
						pad_axisLeft.push_back(new TPad(Form("PadAxisLeft_%i_%i",x,y),"",(Float_t)( l*x),     (Float_t)( otstupDown + h*(y-1)), 
	            													               (Float_t)( l*(x+1)), (Float_t)( otstupDown + h*y), 0,0,0));
						pad_axisLeft.back()->SetFrameBorderMode(1);
	            		pad_axisLeft.back()->SetLeftMargin(left);
	            		pad_axisLeft.back()->SetRightMargin(0);
	            		pad_axisLeft.back()->SetBottomMargin(0);
	            		pad_axisLeft.back()->SetTopMargin(0);
		            }
		            if(y==1){
            			//std::cout<<(Float_t)(otstupLeft + l*x)<<"\t"<< (Float_t)( 1.0 - otstupUp - otstupDown - h*y)<<"\t"<< (Float_t)(otstupLeft + l*(x+1))<<"\t"<< (Float_t)( 1.0 - otstupUp - otstupDown - h*(y-1))<<"\n\n";
						pad_axisDown.push_back(new TPad(Form("PadAxisDown_%i_%i",x,y),"",(Float_t)(otstupLeft + l*x),     (Float_t)( h*(y-1)), 
	            															     (Float_t)(otstupLeft + l*(x+1)), (Float_t)( h*y), 0,0,0));
						pad_axisDown.back()->SetFrameBorderMode(1);
	            		pad_axisDown.back()->SetLeftMargin(0.0);
	            		pad_axisDown.back()->SetRightMargin(0);
	            		pad_axisDown.back()->SetBottomMargin(down);
	            		pad_axisDown.back()->SetTopMargin(0);
		            }

            	}
            	
            }

            TPad *pad_AxisY = new TPad("pad_AxisY", "", 0.0, 0.1, 0.1, 1.0, 0, 0, 0);
			TPad *pad_AxisX = new TPad("pad_AxisX", "", 0.0, 0.0, 1.0, 0.07, 0, 0, 0);

			pad_AxisY->SetFrameBorderMode(0);
            pad_AxisY->SetRightMargin(0.23);
            pad_AxisY->SetBottomMargin(0);
            
            pad_AxisX->SetTopMargin(0);
            pad_AxisX->SetFrameBorderMode(0);
            pad_AxisX->SetBottomMargin(0);

			for(Int_t i=0; i<pad_axisDown.size(); i++){
            	pad_axisDown[i]->Draw();
            }
            for(Int_t i=0; i<pad_axisLeft.size(); i++){
            	pad_axisLeft[i]->Draw();
            }
            for(Int_t i=0; i<pad.size(); i++){
            	pad[i]->Draw();
            }
            pad_AxisY->Draw();
            pad_AxisX->Draw();

            for(Int_t i=0; i<pad_axisDown.size(); i++){
            	pad_axisDown[i]->cd();
            	hPadNM[0]->Draw();
            }
            n=0;
            for(Int_t i=0; i<pad_axisLeft.size(); i++){
            	pad_axisLeft[i]->cd();
            	hPadNM[n]->Draw();
            	if(n<fiew-1){
            		n++;
            	}
            }

            n=0;
            for(Int_t i=0; i<pad.size(); i++){
            	pad[i]->cd();
            	if(i%column==0 && i!=0 && n<fiew-1){
            		n++;
            	}
            	hPadNM[n]->Draw();
            }
            
            pad[plLeg]->cd();
            if(legend_flag==true){
    	       legend->Draw("same");
            }

            Int_t p=0;
            std::cout<<"goog\n";
            for(Int_t i=0; i<DrawResult.size(); i++){
            	pad[i]->cd();
				_textInfo.SetTextSize(labels_b[p].size);
	            _textInfo.SetTextAngle(labels_b[p].angle);
            	_textInfo.DrawLatex(labels_b[i].x, labels_b[i].y, (labels_b[i].text).Data() );
            	if(strncmp(labels[p].place, "Pad",3)==0 /*&& i%2==0*/ ){
	            	_textInfo.SetTextSize(labels[p].size);
	            	_textInfo.SetTextAngle(labels[p].angle);
					_textInfo.DrawLatex(labels[p].x, labels[p].y, (labels[p].text).Data() );
					p++;
				}
            	SetDrawObject(DrawResult[i],"OnePad");
	            if(Sys_mod==true){
	            	for(Int_t j=0; j<graph_up_sys.size();j++){
	            		graph_up_sys[j]->Draw("e5same");
	            	}
	            }
            	for(Int_t j=0; j<graph_up.size();j++){
            		graph_up[j]->Draw("Psame");
            	}
            	for(Int_t l = 0; l < Line_vector.size(); l++){
					if(strncmp(Line_vector_DrawOptions[l], "Up",2)==0){
						Line_vector[l]->Draw("same");
					}
				}
            	graph_up.clear();
            	graph_up_sys.clear();
            }
            std::cout<<"goog\n";

            pad_AxisY->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "Y",1)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
				}
			}

			pad_AxisX->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "X",1)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
				}
			}
			
			pad[plInfo]->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "Info",4)==0 || strncmp(labels[i].place, "V2",2)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
					//std::cout<<"Good\n";
				}
			}
			
			pad[1]->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "V3",2)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
					//std::cout<<"Good\n";
				}
			}
			
			//pad[line*column-1]->cd();
			//pad[line*column-2]->cd();
			pad[0]->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "Prel",4)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
					//std::cout<<"Good\n";
				}
			}

			pad[1]->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "2Prel",5)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
					//std::cout<<"Good\n";
				}
			}

			/*
			pad_AxisX->cd();
			TLatex textArticle;
			textArticle.SetTextSize(0.3);
			textArticle.DrawLatex(0.05,0.4,"#font[42]{ *DOI: 10.1103/PhysRevC.93.014907 }");
			*/
			/*
			pad[5]->cd();
			TLatex textArticle;
			textArticle.SetTextSize(0.12);
			textArticle.DrawLatex(0.2,0.21,"#font[42]{ #scale[1.0]{STAR Preliminary}}");
			*/	
////////////////////////////////////////////////
			/*
			pad[0]->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "C0",2)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
					//std::cout<<"Good\n";
				}
			}

			pad[1]->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "C1",2)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
					//std::cout<<"Good\n";
				}
			}

			pad[2]->cd();
			for(Int_t i=0; i<labels.size(); i++){
				if(strncmp(labels[i].place, "C2",2)==0 ){
					_textInfo.SetTextSize(labels[i].size);
            		_textInfo.SetTextAngle(labels[i].angle);
					_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
					//std::cout<<"Good\n";
				}
			}
			*/
            return canvas;

		}

		TCanvas *CanvasNxMRatio(Int_t column, Float_t left, Float_t down, std::vector<std::vector<Draw_Object*>> DrawResult, std::vector<std::vector<Draw_Object*>> DrawRatio, TString text, Double_t xText, Double_t yText, std::vector<TString> PadText, Double_t xTextPad, Double_t yTextPad, TString Xname, TString YnameUp, TString YnameDown, TString PhysRef){

			Int_t xCanvas=1600;
			Int_t yCanvas=1000;

			//Float_t
			Int_t line=2;
			Float_t otstupUp=0.02;
			Float_t otstupDown=0.1;
			Float_t otstupLeft=0.1;
			Float_t otstupRight=0.02;
			Float_t h = (1.0 - otstupUp - otstupDown )/line;
			Float_t l = (1.0 - otstupRight - otstupLeft )/column;

			TH2F *hPad = new TH2F("hPad","", 2,-1.05,1.05,2,0.6,1.15);
            hPad->GetYaxis()->SetLabelSize(0.07);
            hPad->GetXaxis()->SetLabelSize(0.07);

			canvas = new TCanvas("canvas","plot",xCanvas,yCanvas);
            canvas->GetFrame()->SetFillColor(21);
            canvas->GetFrame()->SetBorderSize(115);
            canvas->cd();
            
            std::vector<TPad *> pad;
            std::vector<TPad *> pad_axisDown;
            std::vector<TPad *> pad_axisLeft;

            for(Int_t y=line; y>0; y--){
            	for(Int_t x=0; x<column; x++){
            		//std::cout<<0.1+l*x<<"\t"<< 1.0-h*y<<"\t"<< 0.1+l*(x+1)<<"\t"<< 1.0-h*(y-1)<<"\n\n";
	            	pad.push_back(new TPad(Form("Pad_%i_%i",x,y),"",(Float_t)(otstupLeft + l*x),     (Float_t)( otstupDown + h*(y-1)), 
	            													(Float_t)(otstupLeft + l*(x+1)), (Float_t)( otstupDown + h*y), 0,0,0));
	            	pad.back()->SetFrameBorderMode(1);
	            	pad.back()->SetLeftMargin(0);
	            	pad.back()->SetRightMargin(0);
	            	pad.back()->SetBottomMargin(0);
	            	pad.back()->SetTopMargin(0);
		            
		            if(x==0){
						pad_axisLeft.push_back(new TPad(Form("PadAx_%i_%i",x,y),"",(Float_t)( l*x),     (Float_t)( otstupDown + h*(y-1)), 
	            													               (Float_t)( l*(x+1)), (Float_t)( otstupDown + h*y), 0,0,0));
						pad_axisLeft.back()->SetFrameBorderMode(1);
	            		pad_axisLeft.back()->SetLeftMargin(left);
	            		pad_axisLeft.back()->SetRightMargin(0);
	            		pad_axisLeft.back()->SetBottomMargin(0);
	            		pad_axisLeft.back()->SetTopMargin(0);
		            }
		            if(y==1){
            			//std::cout<<(Float_t)(otstupLeft + l*x)<<"\t"<< (Float_t)( 1.0 - otstupUp - otstupDown - h*y)<<"\t"<< (Float_t)(otstupLeft + l*(x+1))<<"\t"<< (Float_t)( 1.0 - otstupUp - otstupDown - h*(y-1))<<"\n\n";
						pad_axisDown.push_back(new TPad(Form("Pad_%i_%i",x,y),"",(Float_t)(otstupLeft + l*x),     (Float_t)( h*(y-1)), 
	            															     (Float_t)(otstupLeft + l*(x+1)), (Float_t)( h*y), 0,0,0));
						pad_axisDown.back()->SetFrameBorderMode(1);
	            		pad_axisDown.back()->SetLeftMargin(0.0);
	            		pad_axisDown.back()->SetRightMargin(0);
	            		pad_axisDown.back()->SetBottomMargin(down);
	            		pad_axisDown.back()->SetTopMargin(0);
		            }

		            /*
		            if(y==1){
	            		pad.back()->SetTopMargin(0.1);
	            	}
	            	if(y==line){
	            		pad.back()->SetBottomMargin(0.1);
	            	}
	            	if(x==0){
	            		pad.back()->SetLeftMargin(0.1);
	            	}
	            	if(x==column){
	            		pad.back()->SetRightMargin(0.1);
	            	}
	            	*/
            	}
            	
            }

            TPad *pad_AxisY = new TPad("pad_AxisY", "", 0.0, 0.05, 0.05, 1.0, 0, 0, 0);
			TPad *pad_AxisX = new TPad("pad_AxisX", "", 0.0, 0.0, 1.0, 0.05, 0, 0, 0);

			pad_AxisY->SetFrameBorderMode(0);
            pad_AxisY->SetRightMargin(0.23);
            pad_AxisY->SetBottomMargin(0);
            
            pad_AxisX->SetTopMargin(0);
            pad_AxisX->SetFrameBorderMode(0);
            pad_AxisX->SetBottomMargin(0);

			for(Int_t i=0; i<pad_axisDown.size(); i++){
            	pad_axisDown[i]->Draw();
            }
            for(Int_t i=0; i<pad_axisLeft.size(); i++){
            	pad_axisLeft[i]->Draw();
            }
            for(Int_t i=0; i<pad.size(); i++){
            	pad[i]->Draw();
            }
            pad_AxisY->Draw();
            pad_AxisX->Draw();


            for(Int_t i=0; i<pad_axisDown.size(); i++){
            	pad_axisDown[i]->cd();
            	hPadOne->Draw();
            }
            pad_axisLeft[0]->cd();
            hPadOne->Draw();
            pad_axisLeft[1]->cd();
            hPad->Draw();

            for(Int_t i=0; i<column; i++){
            	pad[i]->cd();
            	hPadOne->Draw();
            }            
            for(Int_t i=column; i<pad.size(); i++){
            	pad[i]->cd();
            	hPad->Draw();
            }
            pad[1]->cd();
            legend->Draw("same");
  			
			pad[0]->cd();
  			_textInfo.SetTextSize(0.07);   	
			_textInfo.DrawLatex(xText, yText, text.Data() );	
  			
  			TLatex textCent;
            textCent.SetTextSize(0.07);
            
            for(Int_t i=0; i<column; i++){
            	pad[i]->cd();
            	textCent.DrawLatex(0.55*1.05, yText*1.05,PadText[i].Data());
            	SetDrawObject(DrawResult[i],"OnePad");
            	for(Int_t j=0; j<graph_up.size();j++){
            		graph_up[j]->Draw("Psame");
            	}
            	graph_up.clear();
            }
            
            for(Int_t i=0; i<column; i++){
            	pad[i+3]->cd();
            	SetDrawObject(DrawRatio[i],"OnePad");
            	for(Int_t j=0; j<graph_up.size();j++){
            		graph_up[j]->Draw("Psame");
            	}
            	for(Int_t l = 0; l < Line_vector.size(); l++){
					if(strncmp(Line_vector_DrawOptions[l], "Down",2)==0){
						Line_vector[l]->Draw("same");
					}
				}
            	graph_up.clear();
            }

            pad_AxisY->cd();
			TLatex textY;
			textY.SetTextSize(0.5);
			//textY.SetTextAngle(90);
			textY.DrawLatex(0.3,0.85,YnameUp.Data());

			textY.SetTextSize(0.4);
			textY.SetTextAngle(90);
			textY.DrawLatex(0.6,0.09,YnameDown.Data() );

			pad_AxisX->cd();
			TLatex textArticle;
			textArticle.SetTextSize(1.0);
			textArticle.DrawLatex(0.75,0.4,Xname.Data());
			textArticle.DrawLatex(0.05,0.4,PhysRef.Data());
            
            return canvas;

		}

};