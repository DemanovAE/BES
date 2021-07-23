typedef vector<Double_t> Vec;

class Draw_Object {
    public:
        // Деструктор
        ~Draw_Object();

        TH1D *histo;

        TGraphErrors *grFromFile;

        Vec axisX;
        Vec axisY;
        Vec axisXerror;
        Vec axisYerror;
        Vec axisSysYerror;

        std::vector<std::vector<Double_t>> sys_arrey;

        TString NameGraph;
        TString Legend;

        Int_t marker_color;
        Int_t marker_style;
       	Float_t marker_size;

       	void SetParametrsGraph(TString Name, TString TitleLegend, Int_t color, Int_t style, Float_t size ){
       		NameGraph = Name;
       		Legend = TitleLegend;
       		marker_color = color;
       		marker_size = size;
       		marker_style = style;
       	}

       	void SetVectorsGraph(Vec x, Vec y, Vec ex, Vec ey){
       		axisX = x;
       		axisY = y;
       		axisXerror = ex;
       		axisYerror = ey;
       	}

       	void SetPointGraph(Double_t x, Double_t y, Double_t ex, Double_t ey){
       		axisX.push_back(x);
       		axisY.push_back(y);
       		axisXerror.push_back(ex);
       		axisYerror.push_back(ey);
       	}

        void SetPointSysErrorGraph(Double_t yse){
          axisSysYerror.push_back(yse);
        }

        void СhangePointGraph(Int_t i, Double_t x, Double_t y, Double_t ex, Double_t ey, Double_t eys){
          axisX[i]=x;
          axisY[i]=y;
          axisXerror[i]=ex;
          axisYerror[i]=ey;
          if(axisSysYerror.size()!=0){ 
            axisSysYerror[i]=eys;
          }
        }

        void SetGraphFromFile(TFile *f){

          grFromFile = (TGraphErrors*) f->Get(NameGraph);
          //grFromFile->SetName(NameGraph.Data());
          grFromFile->SetTitle(Legend.Data());
          grFromFile->SetMarkerColor(marker_color);
          grFromFile->SetLineColor(marker_color);
          grFromFile->SetLineWidth(2.);
          grFromFile->SetMarkerStyle(marker_style);
          grFromFile->SetMarkerSize(marker_size);
          //gr->SetMarkerWidth(2.0);
        }

        void SetSysArrey(std::vector<Double_t> vec){
          sys_arrey.push_back(vec);
        }

        void SetYSysErrors(std::vector<Double_t> y){
          axisSysYerror=y;
        }

        void SetSysErrors(){

        Double_t ey;
        Double_t mean;

        if(sys_arrey.size()!=0){
          
          for(Int_t p = 0; p < sys_arrey[0].size(); p++){
            ey=0.;
            mean=0.;
            for(Int_t i=0; i<(int)sys_arrey.size(); i++){
              mean=mean+sys_arrey[i][p];
            }
            //std::cout<<"meanE\t"<<mean<<"\n";
            mean = mean / sys_arrey.size();
            //std::cout<<"meanE\t"<<mean<<"\n";

            for(Int_t i=0; i<(int)sys_arrey.size(); i++){
              ey = ey + pow( mean - sys_arrey[i][p] ,2.);
            }
            //std::cout<<"E\t"<<ey<<"\n";
            ey = sqrt( ey / sys_arrey.size() );
            axisSysYerror.push_back( ey );
            //std::cout<<"Objecterr\t"<<(int)sys_arrey.size()<<"\t"<<(int)sys_arrey[0].size()<<"\t"<<axisSysYerror.back()<<"\t"<<"\n";
          }
          //std::cout<<"next\n";
          /*
          for(Int_t p = 0; p < sys_arrey[0].size(); p++){
            ey=0.;
            for(Int_t i=1; i<(int)sys_arrey.size(); i++){
              ey = ey + pow( (sys_arrey[0][p] - sys_arrey[i][p]) , 2 );
            }
            axisSysYerror.push_back( sqrt( ey / (int)sys_arrey.size() ));
            std::cout<<(int)sys_arrey.size()<<"\t"<<ey<<"\t"<<ey/sys_arrey[0][p]<<"\t"<<sqrt( ey / (int)sys_arrey.size())<<"\n";
          }
          */
          /*
          std::vector<Double_t> err;
          for(Int_t p = 0; p < sys_arrey[0].size(); p++){
            for(Int_t i=0; i<(int)sys_arrey.size(); i++){
              err.push_back(sys_arrey[i][p]);
            }
            auto minmax = std::minmax_element(err.begin(), err.end());
            axisSysYerror.push_back( abs(*minmax.first-*minmax.second)/2. );

            std::cout<<sys_arrey[0][p]<<"\t"<<*minmax.first<<"\t"<<*minmax.second<<"\t"<<abs(*minmax.first-*minmax.second)<<"\t"<<abs(*minmax.first-*minmax.second)/(2.*sys_arrey[0][p])<<"\n";
            err.clear();
          }
          */
          }
        }

        void ClearSysArrey(){
          sys_arrey.clear();
        }

        Double_t GetPointXGraph(Int_t i){
          return axisX[i]; 
        }
        Double_t GetPointXErrorGraph(Int_t i){
          return axisXerror[i]; 
        }
        Double_t GetPointYGraph(Int_t i){
          return axisY[i]; 
        }
        Double_t GetPointYErrorGraph(Int_t i){
          return axisYerror[i]; 
        }
        Double_t GetPointYSysErrorGraph(Int_t i){
          if(axisSysYerror.size()!=0) return axisSysYerror[i];
          return 0.; 
        }

        std::vector<std::vector<Double_t>> GetSysArrey(){
          return sys_arrey;
        }

       	// get functions

    		TGraphErrors *GetTGraph(){

    			Double_t x_arrey[100];
    			Double_t y_arrey[100];
    			Double_t x_arrey_err[100];
    			Double_t y_arrey_err[100];

    			for(Int_t i=0; i < Draw_Object::GetSizeVector(); i++){
    				x_arrey[i] = axisX[i];
    				y_arrey[i] = axisY[i];
    				x_arrey_err[i] = axisXerror[i];
    				y_arrey_err[i] = axisYerror[i];
    				//std::cout<<y_arrey[i]<<"\t";
    			}
    			//std::cout<<"stop\n";			
    			TGraphErrors *gr = new TGraphErrors(Draw_Object::GetSizeVector(), x_arrey, y_arrey, x_arrey_err, y_arrey_err);
    			gr->SetName(NameGraph.Data());
    			gr->SetTitle(Legend.Data());
    		  gr->SetMarkerColor(marker_color);
    		  gr->SetLineColor(marker_color);
    		  gr->SetMarkerStyle(marker_style);
    		  gr->SetMarkerSize(marker_size);
          //gr->SetMarkerWidth(2.0);


          return gr;
    		  //return grFromFile;
    		
    		}

        TGraphErrors *GetTGraphSys(){

          Double_t x_arrey[100];
          Double_t y_arrey[100];
          Double_t x_arrey_err[100];
          Double_t y_arrey_err[100];

          for(Int_t i=0; i < Draw_Object::GetSizeVector(); i++){
            x_arrey[i] = axisX[i];
            y_arrey[i] = axisY[i];
            x_arrey_err[i] = 0.03;
            //x_arrey_err[i] = 1.5;
            y_arrey_err[i] = axisSysYerror[i];
            //std::cout<<y_arrey_err[i]<<"\t";
          }
          //std::cout<<"stop\n";      
          TGraphErrors *gr = new TGraphErrors(Draw_Object::GetSizeVector(), x_arrey, y_arrey, x_arrey_err, y_arrey_err);
          gr-> SetFillColorAlpha(GetColorRoot(marker_color), 0.1);
          gr->SetName( Form("%s_systematic",NameGraph.Data()) );
          gr->SetTitle(Legend.Data());
          gr->SetMarkerColor(marker_color);
          gr->SetLineColor(marker_color);
          //gr->SetLineWidth(2.0);
          //gr->SetMarkerStyle(marker_style);
          //gr->SetMarkerSize(marker_size);

          return gr;
        
        }

        Color_t GetColorRoot(Int_t c){
          if(c == 6) return kMagenta; 
          if(c == 4) return kBlue; 
          if(c == 3) return kGreen; 
          if(c == 2) return kRed; 
          if(c == 1) return kBlack;
          return kBlack;
        }

        void GetPointDrawObject(Int_t i, Double_t &x, Double_t &y, Double_t &ex, Double_t &ey){
          x = axisX[i];
          y = axisY[i];
          ex = axisXerror[i];
          ey = axisYerror[i];
        }

        void ClearObject(){
          axisX.clear();
          axisY.clear();
          axisXerror.clear();
          axisYerror.clear();
          axisSysYerror.clear();
          sys_arrey.clear();
        }

    		Vec GetAxisX(){ return axisX; }
        Vec GetAxisY(){ return axisY; }
        Vec GetAxisXerrors(){ return axisXerror; }
        Vec GetAxisYerrors(){ return axisYerror; }
           	
        Int_t GetColor(){ return marker_color; }
        Int_t GetMarker(){ return marker_style; }
        Float_t GetMarkerSize(){ return marker_size; }
        
        TString GetNameGraph(){ return NameGraph; }
        TString GetLegend(){ return Legend; }

        Int_t GetSizeVector(){
          std::vector<int> VecSize = {(int)axisX.size(),(int)axisY.size(),(int)axisXerror.size(),(int)axisYerror.size()};
          std::sort(VecSize.begin(),VecSize.end());
          return VecSize[0];
        }

        Int_t GetSizeVectorSys(){
          return axisSysYerror.size();
        }
};