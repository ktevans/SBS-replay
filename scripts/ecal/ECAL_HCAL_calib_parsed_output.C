#include "TTree.h"
#include "TChain.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "gep_tree_parsed.C"
#include <vector>
#include <set>
#include <map>
#include <vector>
#include "TCut.h"
#include "TTreeFormula.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"

using namespace std;

const int nchan_ECAL = 1656;
const int nchan_HCAL = 288;

const int nrows = 69;
const int ncols[nrows] =
  { 12, 12, 12,
    18, 18, 18,
    21, 21, 21, 21, 21, 21,
    24, 24, 24, 24, 24, 24,
    27, 27, 27, 27, 27, 27,
    27, 27, 27, 27, 27, 27,
    27, 27, 27, 27, 27, 27,
    27, 27, 27, 27, 27, 27,
    27, 27, 27, 27, 27, 27,
    27, 27, 27, 27, 27, 27,
    24, 24, 24, 24, 24, 24,
    21, 21, 21, 21, 21, 21,
    18, 18, 18 };
    
    
    

void ECAL_HCAL_calib_parsed_output(const char *configfilename, const char *outfilename="GEP_ECAL_HCAL_calib_temp.root"){
  // What inputs do we need?
  //Main thing is we want a list of "bad/dead" channels. To start with we'll simply ignore the effects of dead channels until we have a proof-of-concept:
  // We need to define the list of files, the list of global cut and/or cut(s), minimum number of events per channel, etc.   

  TFile *fout = new TFile(outfilename,"RECREATE");
  
  ifstream configfile(configfilename);

  if( !configfile ) return;

  TString currentline;

  TChain *C = new TChain("Tout");
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  TCut globalcut = "";
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }

  int minhits_ECAL = 100;
  int minhits_HCAL = 500; 

  int niter = 3;
  double minMoverB = 0.01;
  double mincoeff = 0.0005;
  double maxcoeff = 0.05;
    
  //We'll also need to load the "xpos" and "ypos" values of the
  //bad channels (might as well define for all the channels) to exclude them effectively

  vector<double> xposECAL(nchan_ECAL);
  vector<double> yposECAL(nchan_ECAL);
  vector<double> defaultGainECAL(nchan_ECAL,0.01); //to be loaded from a file

  bool gotoldgain = false;
  
  TString oldgainfname="";
  TString fnamexpos="";
  TString fnameypos="";
  
  set<int> badchanlist_ECAL;
  set<int> badchanlist_HCAL; 

  double dtcut_ECAL = 10.0; //Max time difference from primary
  double dtcut_HCAL = 10.0; //Max time difference from primary

  double Rmax_ECAL = 0.15; //cm (max distance from seed to include in the fit)
  double Rmax_HCAL = 0.5; //cm (max distance from seed to include in the fit)

  double hit_thresh_ECAL = 0.01; //20 MeV
  double seed_thresh_ECAL = 0.1; // 100 MeV
  double clust_thresh_ECAL = 0.5; // 500 MeV

  // These are in GeV
  double hit_thresh_HCAL = 0.02;
  double seed_thresh_HCAL = 0.1;
  double clust_thresh_HCAL = 0.22; 
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");

      if( tokens->GetEntries() >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	TString sval = ( (TObjString*) (*tokens)[1] )->GetString();

	if( skey == "niter" ){
	  niter = sval.Atoi();
	}

	if( skey == "mincoeff" ){
	  mincoeff = sval.Atof();
	}

	if( skey == "maxcoeff" ){
	  maxcoeff = sval.Atof();
	}
	
	if( skey == "minhits_ECAL" ){
	  minhits_ECAL = sval.Atoi();
	}

	if( skey == "minhits_HCAL" ){
	  minhits_HCAL = sval.Atoi();
	}

	if( skey == "dtcut_ECAL" ){
	  dtcut_ECAL = sval.Atof();
	}
	
	if( skey == "dtcut_HCAL" ){
	  dtcut_HCAL = sval.Atof();
	}
	if( skey == "Rmax_ECAL" ){
	  Rmax_ECAL = sval.Atof();
	}
	if( skey == "Rmax_HCAL" ){
	  Rmax_HCAL = sval.Atof();
	}

	if( skey == "badchanlist_ECAL" ){
	  for( int i=1; i<tokens->GetEntries(); i++ ){
	    TString schan = ( (TObjString*) (*tokens)[i] )->GetString(); 
	    badchanlist_ECAL.insert( schan.Atoi() );
	  }
	}

	if( skey == "badchanlist_HCAL" ){
	  for( int i=1; i<tokens->GetEntries(); i++ ){
	    TString schan = ( (TObjString*) (*tokens)[i] )->GetString(); 
	    badchanlist_HCAL.insert( schan.Atoi() );
	  }
	}
	//TODO: make hit/seed/cluster thresholds configurable
	if( skey == "oldgaincoeffs" ){
	  oldgainfname = sval;

	  ifstream oldgainfile(oldgainfname.Data());

	  if( oldgainfile ){
	    int ncoeffsread = 0;
	    double coefftemp;
	    while( oldgainfile >> coefftemp && ncoeffsread < nchan_ECAL ){
	      defaultGainECAL[ncoeffsread] = coefftemp;
	      ncoeffsread++;
	    }
	    if( ncoeffsread == nchan_ECAL ) {
	      gotoldgain = true;
	      cout << "successfully loaded old gain coeffs from file " << oldgainfname << endl;
	    }
	  }
	  
	}
	
      }
    }
  }

  vector<int> Nhit_chan_ECAL(nchan_ECAL,0);
  vector<int> Nhit_chan_HCAL(nchan_HCAL,0);
  
  
  TMatrixD MECAL(nchan_ECAL,nchan_ECAL);
  TVectorD bECAL(nchan_ECAL);

  vector<bool> goodcalibECAL(nchan_ECAL,false);
  vector<double> PreviousSolutionECAL(nchan_ECAL);
  vector<double> sumEmiss_ratio(nchan_ECAL,0.0);
  vector<double> neventsEmiss_ratio(nchan_ECAL,0.0);

  set<int> list_uncalibrated; 
  
  for( int i=0; i<nchan_ECAL; i++ ){
    bECAL(i) = 0.0;
    //goodcalibECAL[i] = false;
    for( int j=0; j<nchan_ECAL; j++ ){
      MECAL(i,j) = 0.0;
    }
  }
  
  TMatrixD MHCAL(nchan_HCAL,nchan_HCAL);
  TVectorD bHCAL(nchan_HCAL);

  for( int i=0; i<nchan_HCAL; i++ ){
    bHCAL(i) = 0.0;
    for( int j=0; j<nchan_HCAL; j++ ){
      MHCAL(i,j) = 0.0;
    }
  }
  
  
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  //What branches do we need to activate? Let's copy-paste from the "ElasticParse" macro:

  gep_tree_parsed *T = new gep_tree_parsed(C);

  long nevent=0;

  int treenum=-1, oldtreenum = -1;

  fout->cd();
  
  //Plots to make: OLD E/p
  TH1D *hEoverP_old_ECAL = new TH1D("hEoverP_old_ECAL","OLD;ECAL E/p; ", 200,0,2.0);
  TH2D *hEoverP_x_old_ECAL = new TH2D("hEoverP_x_old_ECAL", "OLD;xECAL (m);E/p",160,-1.6,1.6,200,0,2);
  TH2D *hEoverP_y_old_ECAL = new TH2D("hEoverP_y_old_ECAL", "OLD;yECAL (m);E/p",160,-0.8,0.8,200,0,2);

  TH1D *hEoverP_temp_ECAL = new TH1D("hEoverP_temp_ECAL","",500,0,5);

  double meanEoverP = 1.0;
  double sigmaEoverP = 0.13;
  
  for( int iter=0; iter<niter; iter++ ){
    cout << "ECAL calib iteration " << iter << ": " <<  endl;

    nevent = 0;
    treenum = -1;
    oldtreenum = -1;

    hEoverP_temp_ECAL->SetTitle(Form("ECAL iteration %d; E/p;",iter));
    hEoverP_temp_ECAL->Reset();
    
    for( int i=0; i<nchan_ECAL; i++ ){
      bECAL(i) = 0.0;
      Nhit_chan_ECAL[i] = 0;
      for( int j=0; j<nchan_ECAL; j++ ){
	MECAL(i,j) = 0.0;
      }
      sumEmiss_ratio[i] = 0.0;
      neventsEmiss_ratio[i] = 0.0;
    }
    
    while( C->GetEntry(nevent++) ){
      treenum = C->GetTreeNumber();
      
      if( treenum != oldtreenum ){
	oldtreenum = treenum;
	GlobalCut->UpdateFormulaLeaves();
      }

      bool passedcut = GlobalCut->EvalInstance(0) != 0;

      //chi^2 = sum_i [ (E_i - sum_j c_j A_j)^2 ]
      // dchi2/dc_k = -2 sum_i [ E_i - sum_j c_j A_j ] A_k = 0
      // coeff multiplying c_j in dchi2/dc_k = sum_j A_j A_k
      // const term = E_i A_k
      
      if( passedcut ){
	//loop on ECAL and HCAL hits, increment matrices:
      
	double Eprime = T->eprime_eth;
	double Tp = T->Q2_eth/(2.0*0.938272);
	
	double xseedECAL, yseedECAL, eseedECAL, tseedECAL; int idseedECAL;
	
	double Esum_good_old = 0.0;
	// double xclust_old = T->xECAL;
	// double yclust_old = T->yECAL;
	double Esum_good_new = 0.0;
	
	bool goodcluster = true;

	int nbadcalib_thiscluster = 0;
	
	if( iter > 0 ){ //on calibration iterations after the first one (if applicable), calculate new EoverP value using the result from previous iteration and check it 
	  for( int i=0; i<T->nblkECAL; i++ ){
	    double A_i = T->ECAL_clusblk_e[i]/T->ECAL_clusblk_again[i];
	    double t_i = T->ECAL_clusblk_atime[i];
	    double x_i = T->ECAL_clusblk_x[i];
	    double y_i = T->ECAL_clusblk_y[i];
	    double e_i = A_i * PreviousSolutionECAL[T->ECAL_clusblk_id[i]];
	    if( i == 0 ){
	      tseedECAL = t_i;
	      xseedECAL = x_i;
	      yseedECAL = y_i;
	      eseedECAL = e_i;
	    }

	    if( fabs( t_i - tseedECAL ) <= dtcut_ECAL &&
		sqrt(pow(x_i - xseedECAL,2) + pow( y_i-yseedECAL,2)) <= Rmax_ECAL &&
		e_i > hit_thresh_ECAL ){
	      Esum_good_new += e_i;
	    }

	    if( !goodcalibECAL[T->ECAL_clusblk_id[i]] ){
	      nbadcalib_thiscluster++;
	    }
	  }

	  if( eseedECAL < seed_thresh_ECAL ) goodcluster = false;
	  
	  if ( fabs( Esum_good_new/Eprime - meanEoverP ) > 3.5*sigmaEoverP ){
	    goodcluster = false;
	  }
	}
	
	if( goodcluster ){
	  for( int i=0; i<T->nblkECAL; i++ ){
	    int id_i = T->ECAL_clusblk_id[i];
	    
	    double x_i = T->ECAL_clusblk_x[i];
	    double y_i = T->ECAL_clusblk_y[i];
	    double e_i = T->ECAL_clusblk_e[i];
	    double t_i = T->ECAL_clusblk_atime[i];
	    double oldgain_i = T->ECAL_clusblk_again[i];
	    double A_i = e_i / oldgain_i;
	    
	    if( iter > 0 ){
	      e_i = A_i * PreviousSolutionECAL[id_i];
	    }
	    
	    if( i == 0 ){
	      xseedECAL = T->ECAL_clusblk_x[i];
	      yseedECAL = T->ECAL_clusblk_y[i];
	      eseedECAL = e_i;
	      tseedECAL = T->ECAL_clusblk_atime[i];
	      idseedECAL = T->ECAL_clusblk_id[i];
	    }	
	    
	    if( eseedECAL < seed_thresh_ECAL ) continue;
	    
	    bool goodblock_i = true;
	    
	    //	if( i > 0 ){
	    double distfromseed = sqrt(pow(x_i-xseedECAL,2)+pow(y_i-yseedECAL,2));
	    if( fabs( t_i - tseedECAL ) > dtcut_ECAL ||
		distfromseed > Rmax_ECAL || e_i < hit_thresh_ECAL ){
	      goodblock_i = false;
	    }
	    
	    if( !goodblock_i ) continue;

	    // If this block DOESN'T have a good calibration and it is the only one
	    // in an otherwise good cluster that doesn't, calculate the average of the
	    // ratio of the calibrated energy of this channel to the "missing" energy from the cluster:
	    
	    if( iter > 0 && !goodcalibECAL[id_i] && nbadcalib_thiscluster == 1 ){
	      // Esum_good_new - e_i is the sum of all the other energies in the cluster
	      // e_i / (Eprime - Eother) is the ratio of the calibrated energy of this block
	      // to the "missing energy" from the cluster
	      sumEmiss_ratio[id_i] += e_i/(Eprime - (Esum_good_new-e_i) );
	      neventsEmiss_ratio[id_i] += 1.0;
	    }
	      
	    if( iter == 0 ) Esum_good_old += e_i;
	    
	    for( int j=0; j<T->nblkECAL; j++ ){
	      int id_j = T->ECAL_clusblk_id[j];
	      
	      double x_j = T->ECAL_clusblk_x[j];
	      double y_j = T->ECAL_clusblk_y[j];
	      double e_j = T->ECAL_clusblk_e[j];
	      double t_j = T->ECAL_clusblk_atime[j];
	      double oldgain_j = T->ECAL_clusblk_again[j];
	      double A_j = e_j / oldgain_j;

	      if( iter > 0 ) e_j = A_j * PreviousSolutionECAL[id_j];
	      
	      bool goodblock_j = true;
	      
	      distfromseed = sqrt(pow(x_j-xseedECAL,2)+pow(y_j-yseedECAL,2));
	      
	      if( fabs( t_j - tseedECAL ) > dtcut_ECAL ||
		  distfromseed > Rmax_ECAL || e_j < hit_thresh_ECAL ){
		goodblock_j = false;
	      }
	      
	      if( !goodblock_j ) continue;
	      
	      MECAL(id_i, id_j) += A_i * A_j;
	      
	    }
	  
	    bECAL(id_i) += A_i * Eprime;
	  
	    Nhit_chan_ECAL[id_i] += 1;
	  
	  }
	}

	if( iter == 0 ){
	  hEoverP_old_ECAL->Fill( Esum_good_old / Eprime );
	  hEoverP_x_old_ECAL->Fill( T->xECAL, Esum_good_old / Eprime );
	  hEoverP_y_old_ECAL->Fill( T->yECAL, Esum_good_old / Eprime );
	  Esum_good_new = Esum_good_old;
	}

	hEoverP_temp_ECAL->Fill( Esum_good_new/Eprime );
	
      }
    
      if( nevent % 1000 == 0 ) cout << "Event " << nevent << endl;
    }
 
    list_uncalibrated.clear();
    //list of uncalibrated channels:

    for( int i=0; i<nchan_ECAL; i++ ){
      if( Nhit_chan_ECAL[i] < minhits_ECAL ){
	bECAL(i) = 1.0;
	MECAL(i,i) = 1.0;
	
	list_uncalibrated.insert( i );
	for( int j=0; j<nchan_ECAL; j++ ){
	  if( j != i ) {
	    MECAL(i,j) = 0.0;
	    MECAL(j,i) = 0.0; //Not necessary, but also not harmful;;
	  }
	}
      }
    }
    
    cout << "constructing SVD matrix... ";
    
    TDecompSVD AECAL(MECAL);
    
    cout << "done. Solving for coefficients... ";
    AECAL.Solve( bECAL );
    cout << "done" << endl;

    if( iter > 0 ){
      hEoverP_temp_ECAL->Fit("gaus","SQ","",meanEoverP-1.5*sigmaEoverP,meanEoverP+1.5*sigmaEoverP);

      meanEoverP = ( (TF1*) (hEoverP_temp_ECAL->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      sigmaEoverP = ( (TF1*) (hEoverP_temp_ECAL->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
      
    }

    double ngoodcalib = 0;
    double sum_calib = 0.0;
    for( int i=0; i<nchan_ECAL; i++ ){
      if( list_uncalibrated.find( i ) == list_uncalibrated.end() && bECAL(i) >= mincoeff && bECAL(i) <= maxcoeff ){
	ngoodcalib += 1.0;
	sum_calib += bECAL(i);
	PreviousSolutionECAL[i] = bECAL(i);
	goodcalibECAL[i] = true;
      } else {
	goodcalibECAL[i] = false;
	if( iter > 0 ){
	  if( neventsEmiss_ratio[i] >= 10. ){
	    double avgEmiss_ratio = sumEmiss_ratio[i]/neventsEmiss_ratio[i]; 
	    // ratio of e_i / (Etrue - Eother) ~= e_i / e_i(true)
	    // NOTE: we will have used "PreviousSolution" to do the energy reconstruction in the preceding loop.
	    double newcoeff = PreviousSolutionECAL[i] / avgEmiss_ratio; // old coeff / (e_i/e_i(true)) = old coeff * e_i("true")/e_i (calib)
	    if( newcoeff >= mincoeff && newcoeff <= maxcoeff ){
	      cout << "NOTE: channel " << i << " gain being calibrated based on missing energy ratio, (ratio, newcoeff)=(" << avgEmiss_ratio << ", " << newcoeff << ")" << endl;
	      goodcalibECAL[i] = true;
	      PreviousSolutionECAL[i] = newcoeff; //unclear how many additional channels we might be able to calibrate this way.
	      bECAL(i) = newcoeff; //Just replace bECAL with new coeff here. This should fix the possibility to have factors outside the user-defined range
	    }
	  }
	}
      }
    }

    double avgcoeff = sum_calib/ngoodcalib;

    for( int i=0; i<nchan_ECAL; i++ ){
      if( !goodcalibECAL[i] ){ //replace this channel's gain with the average gain of all "good" channels or with the default gain from the file if provided by the user
	if( gotoldgain ){
	  bECAL(i) = defaultGainECAL[i];
	  PreviousSolutionECAL[i] = bECAL(i);
	} else { // if the user didn't provide anything, use the average "good" calibration result:
	  bECAL(i) = avgcoeff;
	  PreviousSolutionECAL[i] = avgcoeff;
	}
      }
    }

    cout << "After iteration " << iter << ", (ngood, meanEoverP, sigmaEoverP)=("
	 << ngoodcalib << ", " << meanEoverP << ", " << sigmaEoverP << ")" << endl;
  } //end loop over iterations:
    
  
  //Now loop over the data again and see how we did:

  nevent = 0;

  oldtreenum = -1;

  TH1D *hEoverP_new_ECAL = new TH1D("hEoverP_new_ECAL","NEW;ECAL E/p; ", 250,0,5.0);
  TH2D *hEoverP_x_new_ECAL = new TH2D("hEoverP_x_new_ECAL", "NEW;xECAL (m);E/p",160,-1.6,1.6,250,0,5);
  TH2D *hEoverP_y_new_ECAL = new TH2D("hEoverP_y_new_ECAL", "NEW;yECAL (m);E/p",160,-0.8,0.8,250,0,5);

  TH1D *hGainCoeffNew = new TH1D("hGainCoeffNew", "NEW Gain Coeffs; Gain (GeV/pC);", 200, 0.0, 0.2 ); 

  double chan[nchan_ECAL], newcoeff[nchan_ECAL], rationewold[nchan_ECAL];

  double sumcoeff = 0.0, ngoodcalib = 0.0;

  cout << "(mincoeff,maxcoeff)=(" << mincoeff << ", " << maxcoeff << ")" << endl;
  
  for( int i=0; i<nchan_ECAL; i++ ){
    //if( list_uncalibrated.find( i ) == list_uncalibrated.end() ){
    if( goodcalibECAL[i] && list_uncalibrated.find(i) == list_uncalibrated.end() ){
      hGainCoeffNew->Fill( bECAL(i) );
      sumcoeff += bECAL(i);
      ngoodcalib += 1.0;
      newcoeff[i] = bECAL(i);
    
      
      if( gotoldgain ){
	rationewold[i] = newcoeff[i]/defaultGainECAL[i]; 
      }
    }
    
    chan[i] = i;
  }

  double avgcoeff = sumcoeff/ngoodcalib;

  for( int i=0; i<nchan_ECAL; i++ ){
    if( list_uncalibrated.find(i) != list_uncalibrated.end() || !goodcalibECAL[i] ){
      if( gotoldgain ){
	newcoeff[i] = defaultGainECAL[i];
	rationewold[i] = newcoeff[i] / defaultGainECAL[i]; //1
      } else {
	newcoeff[i] = avgcoeff;
      }
      
    }
  }

  TString dbfilename = outfilename;

  dbfilename.ReplaceAll(".root",".dat");
  dbfilename.Prepend("db_");
  
  ofstream dbfile( dbfilename.Data() );

  dbfile << "earm.ecal.adc.gain = " << endl;
  int rowtemp=0, coltemp=0;
  for( int i=0; i<nchan_ECAL; i++ ){
    if( coltemp == ncols[rowtemp] ){
      rowtemp++;
      coltemp = 0;
      dbfile << endl;
    }
    dbfile << Form( "%10.6g ", newcoeff[i] );
    coltemp++;
  }

  TGraph *gNewCoeff = new TGraph( nchan_ECAL, chan, newcoeff );
  gNewCoeff->SetMarkerStyle(20);
  fout->cd();
  gNewCoeff->Write("gNewCoeff");

  if( gotoldgain ){
    TGraph *gNewOldRatio = new TGraph( nchan_ECAL, chan, rationewold );
    gNewOldRatio->SetMarkerStyle(20);

    TH1D *hNewOldRatio = new TH1D("hNewOldRatio", "NEW/OLD; Ratio;", 500,0,5);

    for( int i=0; i<nchan_ECAL; i++ ){
      hNewOldRatio->Fill( rationewold[i] );
    }
    
    fout->cd();
    gNewOldRatio->Write("gNewOldRatio");
  }
  
  while( C->GetEntry(nevent++) ){
    treenum = C->GetTreeNumber();
    
    if( treenum != oldtreenum ){
      oldtreenum = treenum;
      GlobalCut->UpdateFormulaLeaves();
    }
    
    bool passedcut = GlobalCut->EvalInstance(0) != 0;
    
    //chi^2 = sum_i [ (E_i - sum_j c_j A_j)^2 ]
    // dchi2/dc_k = -2 sum_i [ E_i - sum_j c_j A_j ] A_k = 0
    // coeff multiplying c_j in dchi2/dc_k = sum_j A_j A_k
    // const term = E_i A_k

    if( passedcut ){
      //loop on ECAL and HCAL hits, increment matrices:
      

      double Eprime = T->eprime_eth;
      double Tp = T->Q2_eth/(2.0*0.938272);
      
      double xseedECAL, yseedECAL, eseedECAL, tseedECAL; int idseedECAL;

      double Esum_good_new = 0.0;
      double xclust_old = T->xECAL;
      double yclust_old = T->yECAL;
      
      for( int i=0; i<T->nblkECAL; i++ ){
	int id_i = T->ECAL_clusblk_id[i];
	double x_i = T->ECAL_clusblk_x[i];
	double y_i = T->ECAL_clusblk_y[i];
	double e_i = T->ECAL_clusblk_e[i];
	double t_i = T->ECAL_clusblk_atime[i];
	double oldgain_i = T->ECAL_clusblk_again[i];
	double A_i = e_i / oldgain_i;

	if( i == 0 ){
	  xseedECAL = T->ECAL_clusblk_x[i];
	  yseedECAL = T->ECAL_clusblk_y[i];
	  eseedECAL = T->ECAL_clusblk_e[i];
	  tseedECAL = T->ECAL_clusblk_atime[i];
	  idseedECAL = T->ECAL_clusblk_id[i];
	}

	if( eseedECAL < seed_thresh_ECAL ) continue;

	bool goodblock_i = true;

	//	if( i > 0 ){
	double distfromseed = sqrt(pow(x_i-xseedECAL,2)+pow(y_i-yseedECAL,2));
	if( fabs( t_i - tseedECAL ) > dtcut_ECAL ||
	    distfromseed > Rmax_ECAL || e_i < hit_thresh_ECAL ){
	  goodblock_i = false;
	}
	
	if( !goodblock_i ) continue;
	
	Esum_good_new += bECAL(id_i) * A_i;
      }

      hEoverP_new_ECAL->Fill( Esum_good_new/ Eprime );
      hEoverP_x_new_ECAL->Fill( T->xECAL, Esum_good_new / Eprime );
      hEoverP_y_new_ECAL->Fill( T->yECAL, Esum_good_new / Eprime );
    }
    if( nevent % 1000 == 0 ) cout << "Final loop, checking results, event " << nevent << endl;
  }
	
  
  fout->Write();
  

}
