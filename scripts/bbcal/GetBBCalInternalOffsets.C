#include "gmn_tree.C"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCut.h"
#include <iostream>
#include "TString.h"
#include "TCanvas.h"
#include <vector>
#include "TPad.h"
#include "TSystem.h"
#include "TLine.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TMinuit.h"
#include "TSpectrum.h"
#include "TClonesArray.h"

void GetBBCalInternalOffsets(const char *configfilename, const std::string det, const std::string offsetBefore="bb_sh_meantime_gen3.dat", const char *outputfilename="BBCalInternalOffsets_temp.root")
{
  gErrorIgnoreLevel = kError; // Ignore all ROOT warnings

  TFile *fout = new TFile(outputfilename,"RECREATE");

  ifstream configfile(configfilename);
  if( !configfile ) return;

  TString currentline;

  TChain *C = new TChain("T");

  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") )
  {
    if( !currentline.BeginsWith("#") )
    {
      C->Add( currentline.Data() );
    }
  } // end while

  gmn_tree *T = new gmn_tree(C);

  long nevent = 0;

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("bb.tr.*",1);
  C->SetBranchStatus("bb.hodotdc.*",1);
  C->SetBranchStatus("bb.tdctrig.*",1);
  C->SetBranchStatus("bb.ps.*",1);
  C->SetBranchStatus("bb.sh.*",1);
  C->SetBranchStatus("sbs.hcal.*",1);
  C->SetBranchStatus("e.kine.*",1);
  C->SetBranchStatus("Ndata.sbs.hcal.*",1);
  C->SetBranchStatus("Ndata.bb.sh.*",1);
  C->SetBranchStatus("Ndata.bb.ps.*",1);

  double c = 0.299792458;
  double etof0 = (1.96+3.0)/c;

  double coinmin = -20.0;
  double coinmax = 20.0;
  double tdiff_max = 20.0;
  int nBlocks = 0;
  int refID = 0;
  double reft0 = 0.0;

  double tr_p = T->bb_tr_p[0];

  if(det == "SH")
  {
    nBlocks = 189;
    refID = 73;
    //refID = 0;
  }
  if(det == "PS")
  {
    nBlocks = 52;
    refID = 25;
  }
  if(det == "HCal")
  {
    nBlocks = 289;
    refID = 150;
  }

  //double walk_default = (0.2113+0.2794)/2.0; // I think this is just for hodo
  vector<double> t0(nBlocks,0.0);

  //Read parameters from initial fit:
  ifstream t0file(offsetBefore);

  if(t0file)
  {
    //cout << "There is a mean offset data file!" << endl;
    for(int i=0; i<nBlocks; i++)
    {
      t0file >> t0[i];
    }
  } // end file if

  TMatrixD M(nBlocks,nBlocks);
  TVectorD b(nBlocks);
  vector<double> nevent_by_block(nBlocks,0.0);

  cout << endl;

  int bins = 0;
  int hMin = 0;
  int hMax = 0;

  if(det == "SH" || det=="PS")
  {
    bins = 600;
    hMin = -30;
    hMax = 30;
  }

  if(det == "HCal")
  {
    bins = 1800;
    hMin = 60;
    hMax = 240;
  }

  TH1D *hdt_ji_before = new TH1D("hdt_ji_before", "All hits uncorrected time difference; t_{j}-t_{i} [ns];",400,-20,20);
  TH1D *hcoin_before = new TH1D("hcoin_before", "Coincidence Time (HCal - BBCal) BEFORE Alignment [ns];",800,50,250);

  TH2D *htime_blockID_before = new TH2D("htime_blockID_before", "All hits; Block ID; ADC Time [ns];",nBlocks,-0.5,nBlocks-0.5,bins,hMin,hMax);

  while(C->GetEntry(nevent++))
  {

    if( nevent % 100000 == 0 ) cout << "event " << nevent << endl;

    if(T->bb_tr_n >= 1 && T->bb_ps_e>0.2 && T->sbs_hcal_e>0.025 && T->e_kine_W2<1.5 && T->e_kine_W2>0.5 && abs(T->bb_tr_vz[0])<0.27)// && abs(((T->bb_ps_e + T->bb_sh_e)/tr_p) - 1.0) <= 0.2)// && T->bb_tdctrig_trigtime < 1.e10 )
    {

      //      Double_t sh_clus_id = T->bb_sh_clus_id;
      //Double_t hcal_clus_id = T->sbs_hcal_clus_id;

      //grab needed track parameters:
      double vz = T->bb_tr_vz[0];
      double pathl = T->bb_tr_pathl[0];
      double etof = pathl/c-etof0; //electron TOF from vertex to hodo.
      double yBBCal = 0.0; 
      double xBBCal = 0.0;
      double maxBlocks = 0.0;
      double hcalTime = T->sbs_hcal_adctime;
      Double_t bbcalTime = T->bb_sh_atimeblk;

      if(det == "SH")
      {
        yBBCal = T->bb_sh_clus_blk_y[0];
        xBBCal = T->bb_sh_clus_blk_x[0];
        //maxBlocks = T->bb_sh_clus_nblk[0];
        maxBlocks = 15;
	      //maxBlocks = T->Ndata_bb_sh_clus_nblk;
      }
      if(det == "PS") 
      {
        yBBCal = T->bb_ps_clus_blk_y[0];
        xBBCal = T->bb_ps_clus_blk_x[0];
        //maxBlocks = T->bb_ps_clus_nblk[0];
        maxBlocks = 10;
	      //maxBlocks = T->Ndata_bb_ps_clus_nblk;
      }
      if(det == "HCal") 
      {
        yBBCal = T->sbs_hcal_clus_blk_y[0];
        xBBCal = T->sbs_hcal_clus_blk_x[0];
        //maxBlocks = T->sbs_hcal_clus_nblk[0];
        maxBlocks = 15;
	      //maxBlocks = T->Ndata_sbs_hcal_clus_nblk;
      }

      double tmean0 = 0.0;
      double blkTime = 0.0;
      double IDblock = 0.0;
      double eBlock = 0.0;
      double blkTime_j = 0.0;
      double IDblock_j = 0.0;
      double eBlock_j = 0.0;

      for(int iblock=0; iblock<maxBlocks; iblock++)
      {
        if(det == "SH")
        {
          IDblock = T->bb_sh_clus_blk_id[iblock];
          blkTime = T->bb_sh_clus_blk_atime[iblock];
          //double xpos = T->bb_sh_clus_blk_x[iblock];
          eBlock = T->bb_sh_clus_blk_e[iblock];
        }
        if(det == "PS")
        {
          IDblock = T->bb_ps_clus_blk_id[iblock];
          blkTime = T->bb_ps_clus_blk_atime[iblock];
          eBlock = T->bb_ps_clus_blk_e[iblock];
        }
        if(det == "HCal")
        {
          IDblock = T->sbs_hcal_clus_blk_id[iblock];
          blkTime = T->sbs_hcal_clus_blk_atime[iblock];
          eBlock = T->sbs_hcal_clus_blk_e[iblock];
        }

        if(iblock == 0) tmean0 = blkTime;
        if(abs(blkTime-tmean0)>tdiff_max) continue;
        if(iblock>=maxBlocks) continue;

        for(int jblock=iblock+1; jblock<maxBlocks; jblock++)
        {
          if(det == "SH")
          {
            IDblock_j = T->bb_sh_clus_blk_id[jblock];
            blkTime_j = T->bb_sh_clus_blk_atime[jblock];
            eBlock_j = T->bb_sh_clus_blk_e[jblock];
          }
          if(det == "PS")
          {
            IDblock_j = T->bb_ps_clus_blk_id[jblock];
            blkTime_j = T->bb_ps_clus_blk_atime[jblock];
            eBlock_j = T->bb_ps_clus_blk_e[jblock];
          }
          if(det == "HCal")
          {
            IDblock_j = T->sbs_hcal_clus_blk_id[jblock];
            blkTime_j = T->sbs_hcal_clus_blk_atime[jblock];
            eBlock_j = T->sbs_hcal_clus_blk_e[jblock];
          }

          if(jblock>=maxBlocks) continue; //redundant?
          if(abs(blkTime-blkTime_j)>tdiff_max) continue;

          nevent_by_block[IDblock] += 1.0;
          nevent_by_block[IDblock_j] += 1.0;

          b(IDblock) += (blkTime-blkTime_j);
          b(IDblock_j) += (blkTime_j-blkTime);
          M(IDblock,IDblock) += 1.0;
          M(IDblock_j,IDblock_j) += 1.0;
          M(IDblock,IDblock_j) += -1.0;
          M(IDblock_j,IDblock) += -1.0;

          hdt_ji_before->Fill(blkTime_j-blkTime);
          if(det=="PS") hcoin_before->Fill(hcalTime-bbcalTime);
          if(det=="SH" && iblock==0) hcoin_before->Fill(T->sbs_hcal_clus_blk_atime[iblock] - blkTime);
          if(det=="HCal" && iblock==0) hcoin_before->Fill(blkTime - T->bb_sh_clus_blk_atime[iblock]);
          htime_blockID_before->Fill(IDblock,blkTime);

        } // end loop through j blocks

      } // end loop through i blocks

    } // end if of basic timing cuts

  } // end while loop through events

  //cout << "I got through the event loop." << endl;

  for(int i=0; i<nBlocks; i++)
  {
    //cout << "I made it into the event block loop." << endl;
    if(nevent_by_block[i] <= 100.0)
    {
      cout << "Warning: block " << i << " fewer than 100 events... zeroing out offsets" << endl;
      b(i) = 0.0;
      M(i,i) = 1.0;
      for(int j=0; j<nBlocks; j++)
      {
	      if(j != i)
        {
	        M(i,j) = M(j,i) = 0.0;
	      }
      }
    } // end if too few events
  } // end for loop through blocks

  TDecompSVD A(M);
  A.Solve(b);
  b.Print();

  double corr = reft0 - b(refID); // reft0 = 0

  for(int i=0; i<nBlocks; i++ )
  {
    b(i) += corr; 
    t0[i] += b(i); 
  }

  b.Print();

  TString dbfilename = outputfilename;
  dbfilename.ReplaceAll(".root",".dat");
  
  ofstream dbfile(dbfilename.Data());

  if (det == "SH") dbfile << "bb.sh.adc.timeoffset = " << endl;
  if (det == "PS") dbfile << "bb.ps.adc.timeoffset = " << endl;

  for(int i=0; i<nBlocks; i++)
  {
    if(nBlocks==189)
    {
      if(i>0 && i%7==0)
      {
        dbfile << endl;
      }
    }
    if(nBlocks==52)
    {
      if(i>0 && i%2==0)
      {
        dbfile << endl;
      }
    }
    if(nBlocks==288)
    {
      if(i>0 && i%12==0)
      {
        dbfile << endl;
      }
    }
    TString entry;
    dbfile << entry.Format("%14.6g", t0[i]) << " ";
  } // end loop over blocks to fill output txt file

  // We will loop through the events again and make some histograms
  nevent = 0;

  TH1D *hdt_ji = new TH1D("hdt_ji", "All hits corrected time difference; t_{j}-t_{i} [ns];",400,-20,20);
  TH1D *hcoin_after = new TH1D("hcoin_after", "Coincidence Time (HCal - BBCal) AFTER Alignment [ns];",800,50,250);

  TH2D *hdt_ji_IDj = new TH2D("hdt_ji_IDj", "All hits; ID j; t_{j}-t_{i} [ns];",nBlocks,-0.5,nBlocks-0.5,400,-20,20);
  TH2D *hdt_ji_IDi = new TH2D("hdt_ji_IDi", "All hits; ID i; t_{j}-t_{i} [ns];",nBlocks,-0.5,nBlocks-0.5,400,-20,20);

  TH2D *h_offsets_ID = new TH2D("h_offsets_ID", "Alignment Offsets vs Block ID; Block ID; Alignment Offsets [ns]",nBlocks,-0.5,nBlocks-0.5,200,-10,10);
  h_offsets_ID->SetMarkerSize(0.5);
  h_offsets_ID->SetMarkerColor(kRed);
  h_offsets_ID->SetMarkerStyle(kFullCircle);

  TH2D *htime_blockID_after = new TH2D("htime_blockID_after", "All hits; Block ID; ADC Time [ns];",nBlocks,-0.5,nBlocks-0.5,bins,hMin,hMax);

  while(C->GetEntry(nevent++))
  {

    if( nevent % 100000 == 0 ) cout << "event " << nevent << endl;

    if(T->bb_tr_n >= 1 && T->bb_ps_e>0.2 && T->sbs_hcal_e>0.025 && T->e_kine_W2<1.5 && T->e_kine_W2>0.5 && abs(T->bb_tr_vz[0])<0.27)// && abs(((T->bb_ps_e + T->bb_sh_e)/tr_p) - 1.0) <= 0.2)// && T->bb_tdctrig_trigtime < 1.e10 )
    {

      //Double_t sh_clus_id = T->bb_sh_clus_id;
      //Double_t hcal_clus_id = T->sbs_hcal_clus_id;

      //grab needed track parameters:
      double vz = T->bb_tr_vz[0];
      double pathl = T->bb_tr_pathl[0];
      double etof = pathl/c-etof0; //electron TOF from vertex to hodo.
      double yBBCal = 0.0; 
      double xBBCal = 0.0;
      double maxBlocks = 0.0;

      double hcalTime = T->sbs_hcal_adctime;
      Double_t bbcalTime = T->bb_sh_atimeblk;

      double idSH = T->bb_sh_idblk;
      double idPS = T->bb_ps_idblk;
      double idHC = T->sbs_hcal_idblk;

      if(det == "SH")
      {
        yBBCal = T->bb_sh_clus_blk_y[0];
        xBBCal = T->bb_sh_clus_blk_x[0];
        //maxBlocks = T->bb_sh_clus_nblk[0];
        maxBlocks = 15;
      }
      if(det == "PS") 
      {
        yBBCal = T->bb_ps_clus_blk_y[0];
        xBBCal = T->bb_ps_clus_blk_x[0];
        //maxBlocks = T->bb_ps_clus_nblk[0];
        maxBlocks = 10;
      }
      if(det == "HCal") 
      {
        yBBCal = T->sbs_hcal_clus_blk_y[0];
        xBBCal = T->sbs_hcal_clus_blk_x[0];
        //maxBlocks = T->bb_ps_clus_nblk[0];
        maxBlocks = 15;
      }

      double tmean0 = 0.0;
      double blkTime = 0.0;
      double IDblock = 0.0;
      double eBlock = 0.0;
      double blkTime_j = 0.0;
      double IDblock_j = 0.0;
      double eBlock_j = 0.0;

      for(int iblock=0; iblock<maxBlocks; iblock++)
      {
        if(det == "SH")
        {
          IDblock = T->bb_sh_clus_blk_id[iblock];
          blkTime = T->bb_sh_clus_blk_atime[iblock]-t0[IDblock];
          eBlock = T->bb_sh_clus_blk_e[iblock];
        }
        if(det == "PS")
        {
          IDblock = T->bb_ps_clus_blk_id[iblock];
          blkTime = T->bb_ps_clus_blk_atime[iblock]-t0[IDblock];
          eBlock = T->bb_ps_clus_blk_e[iblock];
        }
        if(det == "HCal")
        {
          IDblock = T->sbs_hcal_clus_blk_id[iblock];
          blkTime = T->sbs_hcal_clus_blk_atime[iblock]-t0[IDblock];
          eBlock = T->sbs_hcal_clus_blk_e[iblock];
        }

        if(iblock == 0) tmean0 = blkTime;
        if(abs(blkTime-tmean0)>tdiff_max) continue;
        if(iblock>=maxBlocks) continue;

        for(int jblock=iblock+1; jblock<maxBlocks; jblock++)
        {
          if(det == "SH")
          {
            IDblock_j = T->bb_sh_clus_blk_id[jblock];
            blkTime_j = T->bb_sh_clus_blk_atime[jblock]-t0[IDblock_j];
            eBlock_j = T->bb_sh_clus_blk_e[jblock];
          }
          if(det == "PS")
          {
            IDblock_j = T->bb_ps_clus_blk_id[jblock];
            blkTime_j = T->bb_ps_clus_blk_atime[jblock]-t0[IDblock_j];
            eBlock_j = T->bb_ps_clus_blk_e[jblock];
          }
          if(det == "HCal")
          {
            IDblock_j = T->sbs_hcal_clus_blk_id[jblock];
            blkTime_j = T->sbs_hcal_clus_blk_atime[jblock]-t0[IDblock_j];
            eBlock_j = T->sbs_hcal_clus_blk_e[jblock];
          }

          if(jblock>=maxBlocks) continue;
          if(abs(blkTime-blkTime_j)>tdiff_max) continue;

          hdt_ji->Fill(blkTime_j-blkTime);
          hdt_ji_IDi->Fill(IDblock,blkTime_j-blkTime);
          hdt_ji_IDj->Fill(IDblock_j,blkTime_j-blkTime);
          if(det=="PS") hcoin_after->Fill(hcalTime-bbcalTime);
          if(det=="SH" && iblock==0) hcoin_after->Fill(T->sbs_hcal_clus_blk_atime[iblock] - blkTime);
          if(det=="HCal" && iblock==0) hcoin_after->Fill(blkTime - T->bb_sh_clus_blk_atime[iblock]);
          htime_blockID_after->Fill(IDblock,blkTime);

        } // end loop through j blocks

      } // end for loop through i blocks

      for (int i=0; i<nBlocks; i++)
      {
        h_offsets_ID->Fill(i,t0[i]);
      } // end for loop through blocks

    } // end if cuts

  } // end second loop through events

  hdt_ji->Draw();
  hdt_ji_before->Draw();
  hdt_ji_IDi->Draw("colz");
  hdt_ji_IDj->Draw("colz");
  hcoin_before->Draw();
  hcoin_after->Draw();
  h_offsets_ID->Draw("SCAT");

  TCanvas *T1DcointimeAlignment = new TCanvas("T1DcointimeAlignment","ADC Coincidence Time Alignment", 1000, 1000, 1000, 1000);
  T1DcointimeAlignment -> Divide(2,1);
  T1DcointimeAlignment -> cd(1);
  hcoin_before -> Draw();
  T1DcointimeAlignment -> cd(2);
  hcoin_after -> Draw();
  T1DcointimeAlignment -> Write();

  TCanvas *T2DblockTimeAlignment = new TCanvas("T2DblockTimeAlignment","ADC Block Time Alignment", 1000, 1000, 1000, 1000);
  T2DblockTimeAlignment -> Divide(1,2);
  T2DblockTimeAlignment -> cd(1);
  htime_blockID_before -> Draw("colz");
  T2DblockTimeAlignment -> cd(2);
  htime_blockID_after -> Draw("colz");
  T2DblockTimeAlignment -> Write();

  fout->Write();
  C->Reset();

} // end main 
