#include <TMath.h>

const Char_t *path="/scratch2/demanov/STAR/BES/";

const int nBinCent = 9;
const int nBinVtxZ_Hadrons = 2;
const int nBinVtxZ_PID = 10;

// eta-gap значения
const int nEtaGap=1;
const int nEtaGapHadrons=1;

const Char_t *particles[] = {"Pion","Kaon","Proton"};
const Char_t *particlesSign[] = {"Pos","Neg"};
const Char_t *partLateX[] = {"#pi^{+}","pi^{-}","K^{+}","K^{-}","p","#bar{p}"};

//Double_t pt_min[3] = {0.2,0.2,0.5};
//Double_t pt_max[3] = {1.5,1.5,2.0};

const std::vector<const Char_t *> NameEtaHadrons = {"Eta15","Eta03","Eta05","Eta07"};
const std::vector<Double_t> EtaVecHadrons = {0.075, 0.15, 0.25, 0.35};

const std::vector<const Char_t *> NameEtaPID = {"Eta01","Eta03","Eta05","Eta07"};
const std::vector<Double_t> EtaVecPID = {0.05, 0.15, 0.25, 0.35};

//RunId range from energy = 7.7 | 11.5 | 14.5 | 19.6 | 27.0 | 39.0 | 62.4 |  GeV
const std::map<const Int_t, Int_t> RunIdMin = {{7, 11110000}, {11, 11145000}, {14, 15045000}, {19, 12110000}, {27, 19130000}, {39, 11095000}, {62, 11077023}};
const std::map<const Int_t, Int_t> RunIdMax = {{7, 11150000}, {11, 11165000}, {14, 15075000}, {19, 12123000}, {27, 19180000}, {39, 11115000}, {62, 11098060}};

//Event cut
const std::map<const Int_t, Double_t> CutVtxZ = {{7, 70.}, {11, 50.}, {14, 70.}, {19, 40.}, {27, 70.}, {39, 40.}, {62, 40.}};
const std::map<const Int_t, Double_t> CutVtxR = {{7, 2.}, {11, 2.}, {14, 1.}, {19, 2.}, {27, 2.}, {39, 2.}, {62, 2.}};
const std::map<const Int_t, Double_t> CutDeltaVtxY = {{7, 0.}, {11, 0.}, {14, 0.8847}, {19, 0.}, {27, 0.}, {39, 0.}, {62, 0.}};

//Track cut
const std::map<const Int_t, Double_t> CutDCAHadronsEP = {{7, 2.}, {11, 2.}, {14, 2.}, {19, 2.}, {27, 2.}, {39, 2.}, {62, 2.}};
const std::map<const Int_t, Double_t> CutDCAHadronsFlow = {{7, 2.}, {11, 2.}, {14, 3.}, {19, 2.}, {27, 2.}, {39, 2.}, {62, 2.}};

//Track cut
const std::map<const Int_t, Double_t> CutDCApidEP = {{7, 2.}, {11, 2.}, {14, 2.}, {19, 2.}, {27, 1.}, {39, 2.}, {62, 2.}};
const std::map<const Int_t, Double_t> CutDCApidFlow = {{7, 2.}, {11, 2.}, {14, 3.}, {19, 2.}, {27, 1.}, {39, 2.}, {62, 2.}};

const int CutnHits = 15;
const double CutnHitsRatio = 0.52;

const double CutPtotMin = 0.15;
const double CutPtotMax = 10.0;

const double CutPtEPMin_Hadrons = 0.2;
const double CutPtEPMax_Hadrons = 2.0;
const double CutPtFlowMax_Hadrons = 5.0;

const double CutPtotEPMin_PID = 0.15;
const double CutPtotEPMax_PID = 5.0;
const double CutPtotFlowMin_PID = 0.15;
//const double CutPtotFlowMax_PID = 10.0;

//const double CutPtotEPMin_PID = 0.2;
//const double CutPtotEPMax_PID = 2.0;
//const double CutPtotFlowMin_PID = 0.2;
const double CutPtotFlowMax_PID = 5.0;


const Char_t *resol[] = {"ew","e","w"};
const Char_t *direction[] = {"East","West"};

//Масса и квадрат массы частиц
const Float_t electron_mass = 0.0005485799;
const Float_t pion_mass = 0.13957061;
const Float_t kaon_mass = 0.493677;
const Float_t proton_mass = 0.9382720813;

const Float_t electron_mass_sqr = 0.000000301;
const Float_t pion_mass_sqr = 0.019479955;
const Float_t kaon_mass_sqr = 0.24371698;
const Float_t proton_mass_sqr = 0.880354499;

const std::map<const Int_t, std::vector<Int_t> > BadRuns = {
	{7,  {11110000, 11110001}}, 
	{11, {11148001, 11148008, 11148009, 11148010, 11148036, 
		  11148055, 11149017, 11149018, 11149040, 11149043, 
		  11150017, 11150029, 11151051, 11151057, 11153045, 
		  11154026, 11154040, 11154059, 11156036, 11156043, 
		  11156044, 11156045, 11157039}}, 
	{14, {15053027, 15053028, 15053029, 15053034, 15053035, 
		  15053048, 15053052, 15053053, 15053054, 15053055, 
		  15053056, 15053057, 15054019, 15054053, 15054054, 
		  15055018, 15055131, 15055133, 15055134, 15055135, 
		  15055136, 15055137, 15055138, 15055139, 15055140, 
		  15055141, 15056001, 15056002, 15056003, 15056004, 
		  15056005, 15056006, 15056007, 15056008, 15056009, 
		  15056113, 15056114, 15056116, 15056117, 15056124, 
		  15056125, 15057001, 15057003, 15057004, 15057006, 
		  15057007, 15057010, 15057011, 15057013, 15057014, 
		  15057018, 15057055, 15057059, 15058006, 15058011, 
		  15058018, 15060061, 15060069, 15061001, 15061002, 
		  15062006, 15065012, 15065014, 15066008, 15066013, 
		  15066017, 15068013, 15068014, 15068016, 15069034, 
		  15069036, 15070009, 15070010}}, 
	{19, {12113081, 12113084, 12114077, 12114079, 12114085, 
		  12114088, 12114089, 12114091, 12114094, 12114095, 
		  12114097, 12114098, 12114099, 12114100, 12114101, 
		  12114102, 12114103, 12114104, 12114110, 12115025, 
		  12115026, 12116015, 12116016, 12116063, 12116084, 
		  12117009, 12117047, 12118010, 12118045, 12119008, 
		  12119009, 12119011, 12119015, 12119016, 12119017, 
		  12119019, 12119020, 12119021, 12119022, 12119023, 
		  12119024, 12119025, 12119027, 12119028, 12119029, 
		  12119030, 12119032, 12119035, 12119036, 12119039, 
		  12120018, 12120073}}, 
	 {27, {19130078,19130079,19130084,19130086,19131001,
		19131003,19131004,19131005,19131006,19131007,
		19131013,19131014,19131015,19131016,19131019,
		19131020,19131026,19131027,19131030,19131031,
		19131039,19131042,19131048,19131049,19131050,
		19132016,19132029,19132031,19132037,19132038,
		19132039,19134008,19134009,19134025,19134044,
		19134050,19135001,19135012,19135039,19135040,
		19136014,19136044,19137008,19137013,19137014,
		19137025,19138004,19138010,19138019,19138025,
		19140040,19141004,
		19142027,19142034,19142035,19142038,19142041,
		19142045,19142049,19142055,19142057,19142059,
		19142062,19142064,19142065,19142068,19143001,
		19143003,19143006,19143007,19144012,19144013,
		19144014,19144018,19144019,19144020,19144024,
		19144025,19144026,19144031,19144032,19144033,
		19144038,19144043,19144044,19145004,19145005,
		19145009,19145010,19145011,19145013,19145014,
		19145015,19145017,19145020,19146007,19149035,
		19156004,19156005,19156006,19156008,19156033,
		19156047,19159003,19159024,19159025,19165006,
		19165008,19165011,19165013,19165016,19165018,
		19167026 }},  
	{39, {11099102, 11099103, 11099104, 11099106, 11099107, 
    	  11099125, 11100004, 11100005, 11100008, 11100010, 
    	  11100011, 11100016, 11100020, 11100071, 11101014, 
    	  11101104, 11102098, 11103009, 11103047, 11103065, 
    	  11105011, 11105029, 11106026, 11106027, 11106028, 
    	  11106029, 11106030, 11106040, 11106041, 11107008, 
    	  11107046, 11107083, 11108040, 11108053, 11108065, 
    	  11108075, 11109092, 11109102, 11109105, 11109104, 
    	  11110005, 11110041, 11110042, 11110086}}, 
	{62, {11080057, 11080060, 11080061, 11080062, 11080069, 
          11081002, 11081026, 11081034, 11081037, 11081052, 
          11082017, 11084010, 11085013, 11085050, 11085051, 
          11085053, 11085054, 11086001, 11086002, 11086006, 
          11086019, 11087057, 11087058, 11087059, 11088003, 
          11089048, 11089070, 11091045, 11091067, 11091068, 
          11091068, 11092013, 11092014, 11092015, 11092016, 
          11092024, 11092025, 11092031, 11092055, 11092065, 
          11092066, 11092068, 11092069, 11092070, 11092071, 
          11092073, 11092075, 11092077, 11092086, 11092087, 
          11092089, 11092090, 11093024, 11093042, 11093080, 
          11095029, 11095029, 11095041, 11095049, 11095065, 
          11095066, 11095069, 11095070, 11095077, 11095078, 
          11095089, 11096005, 11096006, 11096007, 11096008, 
          11096010, 11096011, 11096012, 11096014, 11096015, 
          11096016, 11096017, 11096018, 11096023, 11096045, 
          11096057, 11096093, 11096110}}
};


const std::map<const Int_t, std::vector<Int_t> > Trigger = {
	{7,  {290001,290004,
          290003,         
          290501,         
          290053, 290054, 
          290002,         
          290070,         
          290060}},
	{11, {310014,310004}},
	{14, {440005, 440015}},
	{19, {340001,340011,340021}},
	{27, {610001, 610011, 610016, 610021, 610026, 610031, 610041, 610051}},
	{39, {280001,280002,280501}},
	{64, {270001,270011,270021}},
};

/// PID ////
//const std::vector<Double_t> ptBinRange = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.2,3.6,4.0,4.4,4.8,5.2};
//const std::vector<Double_t> ptBinRange = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.2,3.6,4.0,4.4,4.8,5.2};
const std::vector<Double_t> ptBinRange = {0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4,2.6,2.8,3.0,3.2};

const Double_t ptUpTPC[3] = {0.55,0.45,1.0};
const Double_t ptDownTPC[3] = {0.2,0.2,0.5};
//const Double_t SqMdown[3] = {-0.04, 0.16, 0.49};
//const Double_t SqMup[3] = {0.09, 0.36, 1.44};
const Double_t SqMdown[3] = {-0.15, 0.2, 0.74};
const Double_t SqMup[3] = {0.1, 0.32, 1.20};

const Double_t SqMdown1[3] = {-0.15, 0.2, 0.74};
const Double_t SqMup1[3] = {0.1, 0.27, 1.20};

const Double_t SqMdown2[3] = {-0.15, 0.23, 0.8};
const Double_t SqMup2[3] = {0.05, 0.25, 1.10};

//const Double_t SqMdown2[3] = {-0.15, 0.23, 0.74};
//const Double_t SqMup2[3] = {0.1, 0.27, 1.20};

const std::map<const Int_t, Double_t> nSigmaTofTpc = {{7, 3.}, {11, 3.}, {14, 3.}, {19, 3.}, {27, 2.}, {39, 3.}, {62, 1.5}};
const std::map<const Int_t, Double_t> nSigmaTpc = {{7, 1.5}, {11, 1.5}, {14, 1.5}, {19, 1.5}, {27, 2.}, {39, 1.5}, {62, 0.75}};
const std::map<const Int_t, Double_t> nSigmaTpcAntSym = {{7, 2.}, {11, 2.}, {14, 2.}, {19, 2.}, {27, 2.}, {39, 2.}, {62, 1.}};

/*
{27, {12172044, 12172045, 12172046, 12172048, 12172049, 
		  12172056, 12173009, 12173018, 12173026, 12173053, 
		  12173054, 12173055, 12173056, 12173057, 12173072, 
		  12174077, 12174096, 12174109, 12175089, 12176046, 
		  12176047, 12176067, 12176069, 12176104, 12177097, 
		  12178034, 12178051, 12178093, 12179068, 12179083, 
		  12179084, 12179085, 12179086}}

		  360001
		  



{19130078, 19130079, 19130085, 19130086, 19131001, 19131003, 19131006, 19131007, 19131009, 19131010, 19131012,
19131014, 19131015, 19131016, 19131019, 19131020, 19131039, 19131045, 19131049, 19131050, 19131052, 19131057,
19131062, 19132016, 19132029, 19132031, 19132035, 19132046, 19132047, 19132063, 19132083, 19133009, 19133012,
19133013, 19133014, 19133018, 19133021, 19133023, 19133032, 19133033, 19133040, 19133041, 19133050, 19133061,
19134005, 19134008, 19134010, 19134013, 19134017, 19134019, 19134025, 19134045, 19134046, 19134047, 19135012,
19135013, 19135014, 19135022, 19135029, 19135039, 19136001, 19136042, 19137001, 19137003, 19137008, 19137009,
19137010, 19137011, 19137013, 19137015, 19137022, 19137027, 19137050, 19137051, 19137052, 19137056, 19138004,
19138008, 19138010, 19138014, 19138025, 19139023, 19139026, 19139027, 19139028, 19139032, 19139033, 19139037,
19139038, 19139041, 19140014, 19140043, 19141004, 19141008, 19141018, 19141019, 19143008, 19143009, 19143010,
19143011, 19143012, 19143013, 19143014, 19143015, 19143016, 19143017, 19144012, 19144013, 19144014, 19144018,
19144019, 19144020, 19144024, 19144025, 19144026, 19144031, 19144032, 19144033, 19144036, 19144042, 19144044,
19144046, 19144047, 19145001, 19145004, 19145005, 19145006, 19145008, 19145009, 19145010, 19145011, 19145013,
19145014, 19145015, 19145017, 19145019, 19145020, 19145031, 19145035, 19145036, 19145040, 19145042, 19145044,
19145050, 19146002, 19146003, 19146007, 19146008, 19146009, 19146012, 19146013, 19146014, 19146016, 19146017,
19146019, 19146020, 19146024, 19146025, 19146026, 19147007,
19133010,19134011,19135011,19136016,19137047,19137053,19137057,19138009,19139022,19139024,19139034,19140009,19142005,19142048,
19147008,19147009,19147010,19147014,19147015,19147016,19156002,19156032,19156044,19156045,19156046,19157013,19158003,19158007,
19158009,19158010,19158011,19158013,19158014,19158015,19158017,19158018,19158019,19160018,19162002,19162005,19165015,19165020,19165021,19167042,
19130078,19130079,19130084,19130086,19131001,
19131003,19131004,19131005,19131007,19131013,
19131014,19131015,19131016,19131019,19131020,
19131026,19131027,19131031,19131042,19131048,
19131050,19132016,19132038,19132039,19134008,
19134009,19134013,19134025,19134044,19134050,
19135001,19135012,19135039,19136014,19136044,
19137008,19137013,19137014,19138010,19138025,
19140040,19141020,19141048,19141049,19141051}



		  */

const Double_t width_nsigma[2][14] = {
{0.501428, 0.490315, 0.486264, 0.479067, 0.464057, 0.469122, 0.465152, 0.465419, 0.469524, 0.474803, 0.455775, 0.448508, 0.438007, 0.42},
{0.500687, 0.490418, 0.483683, 0.48006, 0.474214, 0.469039, 0.467315, 0.464093, 0.467129, 0.466397, 0.465381, 0.444488, 0.430990, 0.42}
};
const Double_t width_m2[2][14] = {
{0.00193145, 0.00477391, 0.0104395, 0.0141619, 0.0204434, 0.0333461, 0.0441183, 0.0526289, 0.0632522, 0.0747594, 0.09149, 0.117983, 0.135, 0.16},
{0.0019165, 0.00472773, 0.00887696, 0.0141702, 0.0205458, 0.0333006, 0.0431941, 0.0519919, 0.0633056, 0.0770517, 0.0902082, 0.112218, 0.135, 0.16}
};
const Double_t mean_nsigma_Pion[2][14] = {
{0.050654, 0.101588, 0.104956, 0.107116, 0.0777653, 0.0802931, 0.0651657, 0.0576706, 0.0466698, 0.0438218, 0.0304439, 0.0181008, 0.00909384, 0.00456414},
{0.048906, 0.100404, 0.105523, 0.105961, 0.0890426, 0.0767333, 0.0604225, 0.0580056, 0.0456894, 0.0316674, 0.0316152, 0.0200551, 0.00941279, 0.00450711}
};
const Double_t mean_nsigma_Kaon[2][14] = {
{4.68427, 2.43037, 1.00972, 0.270689, -0.152001, -0.654553, -0.598583, -0.716279, -0.798108, -0.862709, -0.804169, -0.832204, -0.9,-1.0489},
{4.69686, 2.59199, 1.01539, 0.277866, -0.147357, -0.64979, -0.626797, -0.71425, -0.799115, -0.838327, -0.895996, -0.844301, -0.9, -1.003802}
};
const Double_t mean_m2_Pion[2][14] = {
{0.0195335, 0.0186167, 0.0174857, 0.0151291, 0.0124021, 0.0102237, 0.00647416, 0.00105405, -0.00535572, -0.0129372, -0.0136014, -0.019053, -0.0002141, -0.0199976},
{0.0195749, 0.0186741, 0.0172804, 0.0152415, 0.0125077, 0.0104346, 0.00650542, 0.00165492, -0.00458698, -0.00991593, -0.0179754, -0.0219985, -0.0002141, -0.0199976}
};
const Double_t mean_m2_Kaon[2][14] = {
{0.244438, 0.248274, 0.242241, 0.24068, 0.238808, 0.23953, 0.233085, 0.230375, 0.226031, 0.223484, 0.218067, 0.211638, 0.23149, 0.024046},
{0.244923, 0.243345, 0.242588, 0.241325, 0.239509, 0.241016, 0.235063, 0.229098, 0.228016, 0.222352, 0.220963, 0.219194, 0.0231382, 0.19808}
};


const Double_t NewXPeakPion[2][14]={
{-6.4114e-06, 5.10704e-05, 9.17645e-06, 0.000773449, 0.00139434, 0.000992906, 0.000966396, 0.00163846, 0.00153716, 0.00107657, -0.00578114, -0.0267947, -0.0504839, -0.00605438},
{-1.09202e-05, 4.95479e-05, 0.000288614, 0.000756128, 0.0014363, 0.000917606, 0.00104634, 0.00131817, 0.000967495, -0.00214529, 1.11686, -0.00904713, -0.00816162, -0.00605438}
};
const Double_t NewYPeakPion[2][14]={
{-4.60885e-06, 1.02398e-05, -5.02405e-07, 0.000144706, -0.000562995, 0.000232557, 0.000733689, 0.00203394, 0.00212163, 0.00301092, 0.00424645, 0.000594141, 0.00911696, -0.00265998},
{-8.42874e-06, 1.34208e-05, 6.87048e-05, 0.00020262, 8.08946e-05, 0.000240139, 0.000384291, 0.00197433, 0.00155731, 0.0013366, 0.0500682, 0.00274295, 0.00249722, 0.00263695}
};
const Double_t NewXWidthPion[2][14]={
{0.00224698, 0.0053578, 0.0104495, 0.0172082, 0.0256597, 0.0354797, 0.0458705, 0.0562039, 0.0672196, 0.0782728, 0.0937541, 0.0844937, 0.0937289, 0.0943056},
{0.00222042, 0.00532091, 0.0103844, 0.0171394, 0.0255121, 0.0352279, 0.045411, 0.0557478, 0.0669434, 0.0791651, 0.08143, 0.0902066, 0.101536, 0.110107}
};
const Double_t NewYWidthPion[2][14]={
{0.00206844, 0.00483793, 0.010412, 0.0141794, 0.0209679, 0.0337795, 0.0449596, 0.0538429, 0.064819, 0.0764339, 0.0953892, 0.102646, 0.1146614, 0.1200864},
{0.002056, 0.00479186, 0.00891549, 0.0141533, 0.0205757, 0.0335974, 0.0436557, 0.0530118, 0.064323, 0.0781557, 0.0531249, 0.0993474, 0.113968, 0.123146}
};
const Double_t NewXPeakKaon[2][14]={
{0.225835, 0.225577, 0.22541, 0.22551, 0.226439, 0.227063, 0.232626, 0.240396, 0.249667, 0.259343, 0.276361, 0.261093, 0.27295, 0.273045},
{0.459969, 0.225918, 0.225853, 0.226025, 0.227002, 0.22699, 0.231263, 0.237511, 0.245107, 0.251781, 0.257, 0.251942, 0.251942, 0.02533439}
};
const Double_t NewYPeakKaon[2][14]={
{-0.00120109, -0.000400447, 0.00053529, 0.000321254, -0.000147668, -0.0175279, -0.00517167, -0.0079061, -0.0129533, -0.0190561, -0.0075946, -0.00287417, -0.01072725, -0.0082671},
{-0.00120109, 0.00160636, 0.000393884, 0.000414966, -0.000117658, -0.0175496, -0.00853181, -0.0107189, -0.0163173, -0.0218954, -0.0200054, -0.00904274, -0.0162725, -0.0194811}
};
const Double_t NewXWidthKaon[2][14]={
{0.00737309, 0.00896097, 0.0134489, 0.0202414, 0.0293864, 0.0411273, 0.0537365, 0.0663426, 0.0803673, 0.0966468, 0.116847, 0.129325, 0.134838, 0.130516},
{0.00737309, 0.00909016, 0.0135239, 0.0204198, 0.0297797, 0.0421915, 0.0553183, 0.0683152, 0.0834815, 0.101966, 0.116847, 0.124976, 0.134838, 0.1387359}
};
const Double_t NewYWidthKaon[2][14]={
{0.00500194, 0.00913308, 0.014487, 0.0174426, 0.0244737, 0.0387578, 0.0510497, 0.060453, 0.0723871, 0.0858504, 0.108043, 0.115812, 0.1200804, 0.1376166},
{0.901324, 0.0091286, 0.0124455, 0.0173681, 0.0239469, 0.0385383, 0.0497944, 0.0601819, 0.0735474, 0.0905886, 0.108043, 0.111531, 0.1211449, 0.133486}
};
const Double_t NewXPeakProton[2][14]={
{0.864598, 0.862577, 0.863655, 0.862281, 0.859688, 0.834124, 0.833898, 0.825772, 0.817884, 0.812939, 0.813759, 0.812092,  0.812092, 0.852204},
{18.7488, 0.863808, 0.864488, 0.864343, 0.862138, 0.837297, 0.836481, 0.827905, 0.819204, 0.810399, 0.813759, 0.822194, 0.812092, 0.852204}
};
const Double_t NewYPeakProton[2][14]={
{0.0304799, 0.0187874, -0.0119754, -0.0500748, -0.093592, -0.22743, -0.23357, -0.27469, -0.315704, -0.351302, -0.3692, -0.427598, -0.437598, -0.393866},
{-13.163, 0.026383, -0.0106411, -0.0501371, -0.0942948, -0.227598, -0.238149, -0.278685, -0.320774, -0.361854, -0.401488, -0.422701, -0.437598, -0.393866}
};
const Double_t NewXWidthProton[2][14]={
{0.0294637, 0.0287788, 0.0297273, 0.0346567, 0.0418724, 0.0504145, 0.0616918, 0.0725177, 0.0842019, 0.0962552, 0.115308, 0.138908, 0.144908, 0.184778},
{0.0294637, 0.0299935, 0.0315782, 0.0371736, 0.0447013, 0.053033, 0.0639467, 0.0746115, 0.0862362, 0.100805, 0.115308, 0.16621, 0.144908, 0.184778}
};
const Double_t NewYWidthProton[2][14]={
{0.00498244, 0.0147606, 0.0243394, 0.0264955, 0.0333857, 0.0490984, 0.0611153, 0.0691049, 0.0798645, 0.0920837, 0.11515, 0.142586, 0.15884, 0.189029},
{0.00498244, 0.0147363, 0.0209833, 0.0267325, 0.0327508, 0.0485005, 0.0587, 0.0674838, 0.0792831, 0.0955951, 0.11515, 0.143592, 0.15884, 0.189029}
};


const Double_t grPion90[2][20] = {{0.054, 0.084, 0.0985, 0.0965, 0.104, 0.1, 0.097, 0.096, 0.092, 0.089, 0.0855, 0.082, 0.078, 0.0705, 0.0715, 0.086, 0.1065, 0.1125, 0.136, 0.1895 },
		{ 0.055, 0.0835, 0.098, 0.101, 0.1, 0.1015, 0.1005, 0.0965, 0.0945, 0.092, 0.0905, 0.088, 0.0885, 0.0865, 0.0945, 0.1125, 0.1365, 0.165, 0.2365, 0.2365 }
};

const Double_t grKaon90[2][20] = {{0.211469, 0.20796, 0.199287, 0.187528, 0.174166, 0.153414, 0.141093, 0.1355, 0.1425, 0.1515, 0.1615, 0.1735, 0.186, 0.199, 0.216, 0.248, 0.3, 0.343, 0.5215, 0.417907},
		{ 0.212096, 0.208294, 0.199907, 0.188166, 0.173077, 0.153887, 0.142139, 0.137, 0.146, 0.1555, 0.167, 0.18, 0.1975, 0.214, 0.235, 0.2685, 0.325, 0.411, 0.4495, 0.056 }
};


/*
const Double_t CutPion27GeV_95[2][20][4] = {
{	{ -0.00451606, 0.05175, -0.0067, 0.00448802}, { -0.0107417, 0.0805, -0.0057, 0.010754}, { -0.021113, 0.097, -0.1813, 0.0210006}, { -0.0345039, 0.098, -0.0069, 0.0345607}, { -0.0514383,
0.098, -0.02209, 0.05125}, { -0.0671147, 0.09425, -0.03712, 0.0667657}, { -0.0775403, 0.09075, -0.04396, 0.0770301}, { -0.0880087, 0.08725, -0.04984, 0.0872053}, { -0.0989049, 0.0835, -0.05522,
0.0978002}, { -0.110667, 0.08025, -0.06068, 0.109326}, { -0.123239, 0.07725, -0.06528, 0.12185}, { -0.137085, 0.07625, -0.06912, 0.136884}, { -0.151293, 0.07625, -0.06976, 0.152414}, { -0.167113,
0.0775, -0.06928, 0.169945}, { -0.184375, 0.079, -0.06597, 0.188281}, { -0.213875, 0.0865, -0.05305, 0.21697}, { -0.257401, 0.095, -0.03443, 0.253996}, { -0.290503, 0.195011, -0.01482, 0.195011},
{ -0.359016, 0.009, 0.0382, 0.279447}, { -0.402484, 0.1045, 0.02192, 0.394717}},
{	{ -0.00447781, 0.0505, -0.0069, 0.00445074}, { -0.0106796, 0.0805, -0.008, 0.0106859}, { -0.0210118, 0.09675, -0.1421, 0.020903}, { -0.0343298, 0.09925, -0.0172, 0.034344}, { -0.0511555,
0.099, -0.03441, 0.0510691}, { -0.0666948, 0.09525, -0.0508, 0.0665445}, { -0.0769689, 0.09175, -0.05932, 0.0766972}, { -0.0873592, 0.0885, -0.06712, 0.0868434}, { -0.0982525, 0.0855, -0.07412,
0.0973525}, { -0.109793, 0.08275, -0.08, 0.10832}, { -0.122401, 0.0805, -0.08656, 0.120409}, { -0.13597, 0.079, -0.09088, 0.133792}, { -0.151458, 0.08275, -0.09824, 0.151263}, { -0.166432, 0.086,
-0.09808, 0.167659}, { -0.183697, 0.091, -0.09447, 0.185505}, { -0.21108, 0.101, -0.08592, 0.21083}, { -0.252578, 0.111, -0.07072, 0.244579}, { -0.302314, 0.1385, -0.01482, 0.298736}, { -0.347207,
0.177, 0.0382, 0.389529}, { -0.373315, 0.04475, 0.02192, 0.289488}}
};

const Double_t CutKaon27GeV_95[2][20][4] = {
{	{ 0.210656, 0.240986, -0.0067, 0.0151652}, { 0.205815, 0.244246, -0.0057, 0.0192158}, { 0.197657, 0.252493, -0.1813, 0.0274179}, { 0.108, 0.266914, -0.0069, 0.0425865}, { 0.11825,
0.284274, -0.02209, 0.0595727}, { 0.12775, 0.302542, -0.03712, 0.0767123}, { 0.135, 0.315871, -0.04396, 0.0886677}, { 0.14325, 0.32994, -0.04984, 0.100314}, { 0.1525, 0.345382, -0.05522, 0.11221},
{ 0.164, 0.362786, -0.06068, 0.12433}, { 0.17675, 0.381902, -0.06528, 0.136785}, { 0.19325, 0.403807, -0.06912, 0.148911}, { 0.2105, 0.424301, -0.06976, 0.159354}, { 0.2305, 0.445364, -0.06928,
0.169345}, { 0.25175, 0.44354, -0.06597, 0.17835}, { 0.28575, 0.42112, -0.05305, 0.186055}, { 0.33325, 0.38198, -0.03443, 0.194728}, { 0.2565, 0.34664, -0.01482, 0.627847}, { 0.317, 0.355, 0.0382,
0.315102}, { 0.317, 0.32042, 0.02192, 0.162774}},
{	{ 0.210576, 0.241728, -0.0069, 0.0155761}, { 0.206315, 0.244789, -0.008, 0.0192373}, { 0.19798, 0.253134, -0.1421, 0.0275773}, { 0.109, 0.267226, -0.0172, 0.0421412}, { 0.11925, 0.285452,
-0.03441, 0.0600471}, { 0.12925, 0.304879, -0.0508, 0.0786886}, { 0.13675, 0.318735, -0.05932, 0.0916993}, { 0.14575, 0.333308, -0.06712, 0.104315}, { 0.15625, 0.348548, -0.07412, 0.116487}, {
0.168, 0.365112, -0.08, 0.128845}, { 0.182, 0.383157, -0.08656, 0.141253}, { 0.19825, 0.403006, -0.09088, 0.153658}, { 0.2205, 0.428123, -0.09824, 0.165125}, { 0.2415, 0.448595, -0.09808,
0.174406}, { 0.26475, 0.46596, -0.09447, 0.182499}, { 0.2995, 0.4519, -0.08592, 0.190526}, { 0.34875, 0.46216, -0.07072, 0.204424}, { 0.4075, 0.42644, -0.01482, 0.195877}, { 0.204, 0.2562, 0.0382,
0.104393}, { 0.317, 0.43024, 0.02192, 0.130647}}
};
*/


/////////////////////////////////////////////////run 18

///////////////////////////////

const Double_t CutPion27GeV_95[2][20][4] = {
{	{ -0.00535041, 0.00533786, -0.0085, 0.00533786}, { -0.0130432, 0.0130702, -0.0098, 0.0130702}, { -0.0258856, 0.02573, -0.01366, 0.02573}, { -0.0421851, 0.0422683, -0.0167, 0.0422683}, { 
-0.0631293, 0.0629584, -0.03342, 0.0629584}, { -0.0814936, 0.0811245, -0.05092, 0.0811245}, { -0.0932503, 0.08425, -0.05944, 0.0925651}, { -0.105796, 0.081, -0.06676, 0.104714}, { -0.119616, 
0.07775, -0.07202, 0.118118}, { -0.134931, 0.0745, -0.07594, 0.133409}, { -0.151536, 0.07325, -0.07824, 0.151129}, { -0.169284, 0.07575, -0.07968, 0.172012}, { -0.187776, 0.081, -0.07984, 
0.194873}, { -0.206274, 0.0875, -0.07968, 0.218566}, { -0.228004, 0.09775, -0.07984, 0.246834}, { -0.25722, 0.12025, -0.07984, 0.290456}, { -0.29269, 0.14925, -0.03728, 0.356429}, { -0.34112, 
0.157, 0.01378, 0.42982}, { -0.401467, 0.1625, 0.09694, 0.530991}, { -0.452745, -0.297, -0.15012, 0.293636}},
{	{ -0.00526192, 0.00526535, -0.0045, 0.00526535}, { -0.0131702, 0.0130934, -0.0095, 0.0130934}, { -0.0257196, 0.0255733, -0.01366, 0.0255733}, { -0.042044, 0.0420459, -0.0281, 0.0420459}, { 
-0.0626952, 0.0626165, -0.04695, 0.0626165}, { -0.0806608, 0.0804208, -0.06496, 0.0804208}, { -0.0922116, 0.084, -0.07468, 0.0916802}, { -0.104602, 0.08175, -0.08368, 0.103663}, { -0.118165, 
0.07925, -0.09092, 0.116779}, { -0.133071, 0.07775, -0.09582, 0.131506}, { -0.149281, 0.0785, -0.09968, 0.148504}, { -0.166575, 0.08225, -0.10256, 0.16787}, { -0.184884, 0.09175, -0.10528, 
0.19044}, { -0.203405, 0.10375, -0.10544, 0.214132}, { -0.224489, 0.11825, -0.10454, 0.240576}, { -0.206113, 0.173872, -0.09865, 0.173872}, { -0.267964, 0.274103, -0.03728, 0.274103}, { -0.30157, 
0.303959, 0.01378, 0.303959}, { -0.411076, 0.25375, 0.09694, 0.522529}, { -0.364756, 0.19144, -0.15012, 0.19144}}
};

const Double_t CutKaon27GeV_95[2][20][4] = {
{	{ 0.209699, 0.243204, -0.0085, 0.0167522}, { 0.20397, 0.24845, -0.0098, 0.0222398}, { 0.193138, 0.259048, -0.0071366, 0.0329549}, { 0.174726, 0.276456, -0.0167, 0.050865}, { 0.151408, 
0.299673, -0.03342, 0.0741329}, { 0.13775, 0.321965, -0.05092, 0.0953449}, { 0.14775, 0.336384, -0.05944, 0.107235}, { 0.15925, 0.352373, -0.06676, 0.119081}, { 0.17275, 0.370095, -0.07202, 
0.131288}, { 0.18925, 0.38973, -0.07594, 0.143843}, { 0.20875, 0.410302, -0.07824, 0.15536}, { 0.2325, 0.430889, -0.07968, 0.164386}, { 0.2595, 0.451095, -0.07984, 0.171698}, { 0.2885, 0.45874, 
-0.07968, 0.178962}, { 0.32525, 0.44734, -0.07984, 0.183668}, { 0.387, 0.4177, -0.07984, 0.188813}, { 0.179428, 0.53854, -0.03728, 0.193729}, { 0.185501, 0.4025, 0.01378, 0.204639}, { 0.188647, 
0.37362, 0.09694, 0.208428}, { 6.89967e-310, 0.413723, -0.15012, 0.297241}},
{	{ 0.209075, 0.24344, -0.0045, 0.0171822}, { 0.204396, 0.248253, -0.0099, 0.0219285}, { 0.19251, 0.259289, -0.0071366, 0.0333898}, { 0.174329, 0.276534, -0.0281, 0.0511024}, { 0.148039, 
0.301651, -0.04695, 0.076806}, { 0.13875, 0.324957, -0.06496, 0.100172}, { 0.14975, 0.339284, -0.07468, 0.112469}, { 0.1625, 0.355031, -0.08368, 0.124472}, { 0.17725, 0.372838, -0.09092, 
0.136999}, { 0.19475, 0.391877, -0.09582, 0.14922}, { 0.216, 0.412318, -0.09968, 0.160375}, { 0.24075, 0.433473, -0.10256, 0.169644}, { 0.27075, 0.456735, -0.10528, 0.176662}, { 0.30325, 0.47964, 
-0.10544, 0.182683}, { 0.34125, 0.47698, -0.10454, 0.186333}, { 0.2305, 0.49636, -0.09865, 0.431174}, { 0.3705, 0.32156, -0.03728, 0.896985}, { 1.18649, -0.398296, 0.01378, -0.792394}, { 0.242936, 
0.30104, 0.09694, 0.0316124}, { 0.175119, 0.250214, -0.15012, 0.0375475}}
};
