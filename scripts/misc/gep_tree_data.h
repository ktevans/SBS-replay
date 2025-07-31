//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 21 17:03:48 2025 by ROOT version 6.30/04
// from TChain T/
//////////////////////////////////////////////////////////

#ifndef gep_tree_data_h
#define gep_tree_data_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "THaEvent.h"
#include "THaEvent.h"

class gep_tree_data {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        ax_new;
   Double_t        bx_new;
   Double_t        ay_new;
   Double_t        by_new;
   Double_t        targx;
   Double_t        targy;
   Double_t        scalhel_true_hel;
   Double_t        scalhel_hel_stable;
   Double_t        dEcalo;
   Double_t        detheta;
   Double_t        dptheta;
   Int_t           Ndata_earm_cdet_hit_col;
   Double_t        earm_cdet_hit_col[166];   //[Ndata.earm.cdet.hit.col]
   Int_t           Ndata_earm_cdet_hit_layer;
   Double_t        earm_cdet_hit_layer[166];   //[Ndata.earm.cdet.hit.layer]
   Int_t           Ndata_earm_cdet_hit_pmtnum;
   Double_t        earm_cdet_hit_pmtnum[166];   //[Ndata.earm.cdet.hit.pmtnum]
   Int_t           Ndata_earm_cdet_hit_row;
   Double_t        earm_cdet_hit_row[166];   //[Ndata.earm.cdet.hit.row]
   Int_t           Ndata_earm_cdet_hit_tdc_le;
   Double_t        earm_cdet_hit_tdc_le[166];   //[Ndata.earm.cdet.hit.tdc_le]
   Int_t           Ndata_earm_cdet_hit_tdc_te;
   Double_t        earm_cdet_hit_tdc_te[166];   //[Ndata.earm.cdet.hit.tdc_te]
   Int_t           Ndata_earm_cdet_hit_tdc_tot;
   Double_t        earm_cdet_hit_tdc_tot[166];   //[Ndata.earm.cdet.hit.tdc_tot]
   Int_t           Ndata_earm_cdet_hit_xhit;
   Double_t        earm_cdet_hit_xhit[166];   //[Ndata.earm.cdet.hit.xhit]
   Int_t           Ndata_earm_cdet_hit_yhit;
   Double_t        earm_cdet_hit_yhit[166];   //[Ndata.earm.cdet.hit.yhit]
   Int_t           Ndata_earm_cdet_hit_zhit;
   Double_t        earm_cdet_hit_zhit[166];   //[Ndata.earm.cdet.hit.zhit]
   Int_t           Ndata_earm_cdet_hits_TDCelemID;
   Double_t        earm_cdet_hits_TDCelemID[1312];   //[Ndata.earm.cdet.hits.TDCelemID]
   Int_t           Ndata_earm_cdet_hits_t;
   Double_t        earm_cdet_hits_t[1312];   //[Ndata.earm.cdet.hits.t]
   Int_t           Ndata_earm_cdet_hits_t_te;
   Double_t        earm_cdet_hits_t_te[1312];   //[Ndata.earm.cdet.hits.t_te]
   Int_t           Ndata_earm_cdet_hits_t_tot;
   Double_t        earm_cdet_hits_t_tot[1312];   //[Ndata.earm.cdet.hits.t_tot]
   Int_t           Ndata_earm_cdet_tdc;
   Double_t        earm_cdet_tdc[166];   //[Ndata.earm.cdet.tdc]
   Int_t           Ndata_earm_cdet_tdc_mult;
   Double_t        earm_cdet_tdc_mult[166];   //[Ndata.earm.cdet.tdc_mult]
   Int_t           Ndata_earm_cdet_tdc_te;
   Double_t        earm_cdet_tdc_te[166];   //[Ndata.earm.cdet.tdc_te]
   Int_t           Ndata_earm_cdet_tdc_tot;
   Double_t        earm_cdet_tdc_tot[166];   //[Ndata.earm.cdet.tdc_tot]
   Int_t           Ndata_earm_cdet_tdccol;
   Double_t        earm_cdet_tdccol[166];   //[Ndata.earm.cdet.tdccol]
   Int_t           Ndata_earm_cdet_tdcelemID;
   Double_t        earm_cdet_tdcelemID[166];   //[Ndata.earm.cdet.tdcelemID]
   Int_t           Ndata_earm_cdet_tdclayer;
   Double_t        earm_cdet_tdclayer[166];   //[Ndata.earm.cdet.tdclayer]
   Int_t           Ndata_earm_cdet_tdcrow;
   Double_t        earm_cdet_tdcrow[166];   //[Ndata.earm.cdet.tdcrow]
   Int_t           Ndata_earm_ecal_a_c;
   Double_t        earm_ecal_a_c[92];   //[Ndata.earm.ecal.a_c]
   Int_t           Ndata_earm_ecal_a_p;
   Double_t        earm_ecal_a_p[92];   //[Ndata.earm.ecal.a_p]
   Int_t           Ndata_earm_ecal_a_time;
   Double_t        earm_ecal_a_time[92];   //[Ndata.earm.ecal.a_time]
   Int_t           Ndata_earm_ecal_adccol;
   Double_t        earm_ecal_adccol[92];   //[Ndata.earm.ecal.adccol]
   Int_t           Ndata_earm_ecal_adcelemID;
   Double_t        earm_ecal_adcelemID[92];   //[Ndata.earm.ecal.adcelemID]
   Int_t           Ndata_earm_ecal_adcrow;
   Double_t        earm_ecal_adcrow[92];   //[Ndata.earm.ecal.adcrow]
   Int_t           Ndata_earm_ecal_adcxpos;
   Double_t        earm_ecal_adcxpos[92];   //[Ndata.earm.ecal.adcxpos]
   Int_t           Ndata_earm_ecal_adcypos;
   Double_t        earm_ecal_adcypos[92];   //[Ndata.earm.ecal.adcypos]
   Int_t           Ndata_earm_ecal_clus_atimeblk;
   Double_t        earm_ecal_clus_atimeblk[5];   //[Ndata.earm.ecal.clus.atimeblk]
   Int_t           Ndata_earm_ecal_clus_col;
   Double_t        earm_ecal_clus_col[5];   //[Ndata.earm.ecal.clus.col]
   Int_t           Ndata_earm_ecal_clus_e;
   Double_t        earm_ecal_clus_e[5];   //[Ndata.earm.ecal.clus.e]
   Int_t           Ndata_earm_ecal_clus_eblk;
   Double_t        earm_ecal_clus_eblk[5];   //[Ndata.earm.ecal.clus.eblk]
   Int_t           Ndata_earm_ecal_clus_id;
   Double_t        earm_ecal_clus_id[5];   //[Ndata.earm.ecal.clus.id]
   Int_t           Ndata_earm_ecal_clus_nblk;
   Double_t        earm_ecal_clus_nblk[5];   //[Ndata.earm.ecal.clus.nblk]
   Int_t           Ndata_earm_ecal_clus_row;
   Double_t        earm_ecal_clus_row[5];   //[Ndata.earm.ecal.clus.row]
   Int_t           Ndata_earm_ecal_clus_x;
   Double_t        earm_ecal_clus_x[5];   //[Ndata.earm.ecal.clus.x]
   Int_t           Ndata_earm_ecal_clus_y;
   Double_t        earm_ecal_clus_y[5];   //[Ndata.earm.ecal.clus.y]
   Int_t           Ndata_earm_ecal_clus_blk_atime;
   Double_t        earm_ecal_clus_blk_atime[13];   //[Ndata.earm.ecal.clus_blk.atime]
   Int_t           Ndata_earm_ecal_clus_blk_col;
   Double_t        earm_ecal_clus_blk_col[13];   //[Ndata.earm.ecal.clus_blk.col]
   Int_t           Ndata_earm_ecal_clus_blk_e;
   Double_t        earm_ecal_clus_blk_e[13];   //[Ndata.earm.ecal.clus_blk.e]
   Int_t           Ndata_earm_ecal_clus_blk_id;
   Double_t        earm_ecal_clus_blk_id[13];   //[Ndata.earm.ecal.clus_blk.id]
   Int_t           Ndata_earm_ecal_clus_blk_row;
   Double_t        earm_ecal_clus_blk_row[13];   //[Ndata.earm.ecal.clus_blk.row]
   Int_t           Ndata_earm_ecal_clus_blk_x;
   Double_t        earm_ecal_clus_blk_x[13];   //[Ndata.earm.ecal.clus_blk.x]
   Int_t           Ndata_earm_ecal_clus_blk_y;
   Double_t        earm_ecal_clus_blk_y[13];   //[Ndata.earm.ecal.clus_blk.y]
   Int_t           Ndata_earm_ecal_ped;
   Double_t        earm_ecal_ped[92];   //[Ndata.earm.ecal.ped]
   Int_t           Ndata_earm_ecal_vtp_clus_e;
   Double_t        earm_ecal_vtp_clus_e[8];   //[Ndata.earm.ecal.vtp.clus.e]
   Int_t           Ndata_earm_ecal_vtp_clus_size;
   Double_t        earm_ecal_vtp_clus_size[8];   //[Ndata.earm.ecal.vtp.clus.size]
   Int_t           Ndata_earm_ecal_vtp_clus_time;
   Double_t        earm_ecal_vtp_clus_time[8];   //[Ndata.earm.ecal.vtp.clus.time]
   Int_t           Ndata_earm_ecal_vtp_clus_x;
   Double_t        earm_ecal_vtp_clus_x[8];   //[Ndata.earm.ecal.vtp.clus.x]
   Int_t           Ndata_earm_ecal_vtp_clus_y;
   Double_t        earm_ecal_vtp_clus_y[8];   //[Ndata.earm.ecal.vtp.clus.y]
   Int_t           Ndata_earm_tr_beta;
   Double_t        earm_tr_beta[1];   //[Ndata.earm.tr.beta]
   Int_t           Ndata_earm_tr_chi2;
   Double_t        earm_tr_chi2[1];   //[Ndata.earm.tr.chi2]
   Int_t           Ndata_earm_tr_d_ph;
   Double_t        earm_tr_d_ph[1];   //[Ndata.earm.tr.d_ph]
   Int_t           Ndata_earm_tr_d_th;
   Double_t        earm_tr_d_th[1];   //[Ndata.earm.tr.d_th]
   Int_t           Ndata_earm_tr_d_x;
   Double_t        earm_tr_d_x[1];   //[Ndata.earm.tr.d_x]
   Int_t           Ndata_earm_tr_d_y;
   Double_t        earm_tr_d_y[1];   //[Ndata.earm.tr.d_y]
   Int_t           Ndata_earm_tr_dbeta;
   Double_t        earm_tr_dbeta[1];   //[Ndata.earm.tr.dbeta]
   Int_t           Ndata_earm_tr_dtime;
   Double_t        earm_tr_dtime[1];   //[Ndata.earm.tr.dtime]
   Int_t           Ndata_earm_tr_flag;
   Double_t        earm_tr_flag[1];   //[Ndata.earm.tr.flag]
   Int_t           Ndata_earm_tr_ndof;
   Double_t        earm_tr_ndof[1];   //[Ndata.earm.tr.ndof]
   Int_t           Ndata_earm_tr_p;
   Double_t        earm_tr_p[1];   //[Ndata.earm.tr.p]
   Int_t           Ndata_earm_tr_pathl;
   Double_t        earm_tr_pathl[1];   //[Ndata.earm.tr.pathl]
   Int_t           Ndata_earm_tr_ph;
   Double_t        earm_tr_ph[1];   //[Ndata.earm.tr.ph]
   Int_t           Ndata_earm_tr_px;
   Double_t        earm_tr_px[1];   //[Ndata.earm.tr.px]
   Int_t           Ndata_earm_tr_py;
   Double_t        earm_tr_py[1];   //[Ndata.earm.tr.py]
   Int_t           Ndata_earm_tr_pz;
   Double_t        earm_tr_pz[1];   //[Ndata.earm.tr.pz]
   Int_t           Ndata_earm_tr_r_ph;
   Double_t        earm_tr_r_ph[1];   //[Ndata.earm.tr.r_ph]
   Int_t           Ndata_earm_tr_r_th;
   Double_t        earm_tr_r_th[1];   //[Ndata.earm.tr.r_th]
   Int_t           Ndata_earm_tr_r_x;
   Double_t        earm_tr_r_x[1];   //[Ndata.earm.tr.r_x]
   Int_t           Ndata_earm_tr_r_y;
   Double_t        earm_tr_r_y[1];   //[Ndata.earm.tr.r_y]
   Int_t           Ndata_earm_tr_tg_dp;
   Double_t        earm_tr_tg_dp[1];   //[Ndata.earm.tr.tg_dp]
   Int_t           Ndata_earm_tr_tg_ph;
   Double_t        earm_tr_tg_ph[1];   //[Ndata.earm.tr.tg_ph]
   Int_t           Ndata_earm_tr_tg_th;
   Double_t        earm_tr_tg_th[1];   //[Ndata.earm.tr.tg_th]
   Int_t           Ndata_earm_tr_tg_y;
   Double_t        earm_tr_tg_y[1];   //[Ndata.earm.tr.tg_y]
   Int_t           Ndata_earm_tr_th;
   Double_t        earm_tr_th[1];   //[Ndata.earm.tr.th]
   Int_t           Ndata_earm_tr_time;
   Double_t        earm_tr_time[1];   //[Ndata.earm.tr.time]
   Int_t           Ndata_earm_tr_vx;
   Double_t        earm_tr_vx[1];   //[Ndata.earm.tr.vx]
   Int_t           Ndata_earm_tr_vy;
   Double_t        earm_tr_vy[1];   //[Ndata.earm.tr.vy]
   Int_t           Ndata_earm_tr_vz;
   Double_t        earm_tr_vz[1];   //[Ndata.earm.tr.vz]
   Int_t           Ndata_earm_tr_x;
   Double_t        earm_tr_x[1];   //[Ndata.earm.tr.x]
   Int_t           Ndata_earm_tr_y;
   Double_t        earm_tr_y[1];   //[Ndata.earm.tr.y]
   Int_t           Ndata_sbs_gemFPP_hit_ADCU;
   Double_t        sbs_gemFPP_hit_ADCU[32];   //[Ndata.sbs.gemFPP.hit.ADCU]
   Int_t           Ndata_sbs_gemFPP_hit_ADCU_deconv;
   Double_t        sbs_gemFPP_hit_ADCU_deconv[32];   //[Ndata.sbs.gemFPP.hit.ADCU_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCV;
   Double_t        sbs_gemFPP_hit_ADCV[32];   //[Ndata.sbs.gemFPP.hit.ADCV]
   Int_t           Ndata_sbs_gemFPP_hit_ADCV_deconv;
   Double_t        sbs_gemFPP_hit_ADCV_deconv[32];   //[Ndata.sbs.gemFPP.hit.ADCV_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCasym;
   Double_t        sbs_gemFPP_hit_ADCasym[32];   //[Ndata.sbs.gemFPP.hit.ADCasym]
   Int_t           Ndata_sbs_gemFPP_hit_ADCasym_deconv;
   Double_t        sbs_gemFPP_hit_ADCasym_deconv[32];   //[Ndata.sbs.gemFPP.hit.ADCasym_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCavg;
   Double_t        sbs_gemFPP_hit_ADCavg[32];   //[Ndata.sbs.gemFPP.hit.ADCavg]
   Int_t           Ndata_sbs_gemFPP_hit_ADCavg_deconv;
   Double_t        sbs_gemFPP_hit_ADCavg_deconv[32];   //[Ndata.sbs.gemFPP.hit.ADCavg_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac0_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac0_Umax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac0_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac0_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac0_Vmax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac0_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac1_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac1_Umax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac1_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac1_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac1_Vmax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac1_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac2_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac2_Umax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac2_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac2_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac2_Vmax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac2_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac3_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac3_Umax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac3_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac3_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac3_Vmax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac3_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac4_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac4_Umax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac4_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac4_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac4_Vmax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac4_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac5_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac5_Umax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac5_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac5_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac5_Vmax[32];   //[Ndata.sbs.gemFPP.hit.ADCfrac5_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv;
   Double_t        sbs_gemFPP_hit_ADCmaxcomboUclust_deconv[32];   //[Ndata.sbs.gemFPP.hit.ADCmaxcomboUclust_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv;
   Double_t        sbs_gemFPP_hit_ADCmaxcomboVclust_deconv[32];   //[Ndata.sbs.gemFPP.hit.ADCmaxcomboVclust_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampU;
   Double_t        sbs_gemFPP_hit_ADCmaxsampU[32];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampU]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampUclust;
   Double_t        sbs_gemFPP_hit_ADCmaxsampUclust[32];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampUclust]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampUclust_deconv;
   Double_t        sbs_gemFPP_hit_ADCmaxsampUclust_deconv[32];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampUclust_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampV;
   Double_t        sbs_gemFPP_hit_ADCmaxsampV[32];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampV]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampVclust;
   Double_t        sbs_gemFPP_hit_ADCmaxsampVclust[32];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampVclust]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampVclust_deconv;
   Double_t        sbs_gemFPP_hit_ADCmaxsampVclust_deconv[32];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampVclust_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxstripU;
   Double_t        sbs_gemFPP_hit_ADCmaxstripU[32];   //[Ndata.sbs.gemFPP.hit.ADCmaxstripU]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxstripV;
   Double_t        sbs_gemFPP_hit_ADCmaxstripV[32];   //[Ndata.sbs.gemFPP.hit.ADCmaxstripV]
   Int_t           Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U;
   Double_t        sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U[32];   //[Ndata.sbs.gemFPP.hit.BUILD_ALL_SAMPLES_U]
   Int_t           Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V;
   Double_t        sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V[32];   //[Ndata.sbs.gemFPP.hit.BUILD_ALL_SAMPLES_V]
   Int_t           Ndata_sbs_gemFPP_hit_CM_GOOD_U;
   Double_t        sbs_gemFPP_hit_CM_GOOD_U[32];   //[Ndata.sbs.gemFPP.hit.CM_GOOD_U]
   Int_t           Ndata_sbs_gemFPP_hit_CM_GOOD_V;
   Double_t        sbs_gemFPP_hit_CM_GOOD_V[32];   //[Ndata.sbs.gemFPP.hit.CM_GOOD_V]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC0_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC0_Umax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC0_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC0_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC0_Vmax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC0_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC1_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC1_Umax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC1_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC1_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC1_Vmax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC1_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC2_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC2_Umax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC2_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC2_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC2_Vmax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC2_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC3_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC3_Umax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC3_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC3_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC3_Vmax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC3_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC4_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC4_Umax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC4_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC4_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC4_Vmax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC4_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC5_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC5_Umax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC5_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC5_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC5_Vmax[32];   //[Ndata.sbs.gemFPP.hit.DeconvADC5_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboU;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxcomboU[32];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxcomboU]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboV;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxcomboV[32];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxcomboV]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxsampU;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxsampU[32];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxsampU]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxsampV;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxsampV[32];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxsampV]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxstripU;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxstripU[32];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxstripU]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxstripV;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxstripV[32];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxstripV]
   Int_t           Ndata_sbs_gemFPP_hit_ENABLE_CM_U;
   Double_t        sbs_gemFPP_hit_ENABLE_CM_U[32];   //[Ndata.sbs.gemFPP.hit.ENABLE_CM_U]
   Int_t           Ndata_sbs_gemFPP_hit_ENABLE_CM_V;
   Double_t        sbs_gemFPP_hit_ENABLE_CM_V[32];   //[Ndata.sbs.gemFPP.hit.ENABLE_CM_V]
   Int_t           Ndata_sbs_gemFPP_hit_TSchi2_Umax;
   Double_t        sbs_gemFPP_hit_TSchi2_Umax[32];   //[Ndata.sbs.gemFPP.hit.TSchi2_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_TSchi2_Vmax;
   Double_t        sbs_gemFPP_hit_TSchi2_Vmax[32];   //[Ndata.sbs.gemFPP.hit.TSchi2_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_TSprob_Umax;
   Double_t        sbs_gemFPP_hit_TSprob_Umax[32];   //[Ndata.sbs.gemFPP.hit.TSprob_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_TSprob_Vmax;
   Double_t        sbs_gemFPP_hit_TSprob_Vmax[32];   //[Ndata.sbs.gemFPP.hit.TSprob_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_Tavg;
   Double_t        sbs_gemFPP_hit_Tavg[32];   //[Ndata.sbs.gemFPP.hit.Tavg]
   Int_t           Ndata_sbs_gemFPP_hit_Tavg_corr;
   Double_t        sbs_gemFPP_hit_Tavg_corr[32];   //[Ndata.sbs.gemFPP.hit.Tavg_corr]
   Int_t           Ndata_sbs_gemFPP_hit_Tavg_deconv;
   Double_t        sbs_gemFPP_hit_Tavg_deconv[32];   //[Ndata.sbs.gemFPP.hit.Tavg_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_Tavg_fit;
   Double_t        sbs_gemFPP_hit_Tavg_fit[32];   //[Ndata.sbs.gemFPP.hit.Tavg_fit]
   Int_t           Ndata_sbs_gemFPP_hit_Ugain;
   Double_t        sbs_gemFPP_hit_Ugain[32];   //[Ndata.sbs.gemFPP.hit.Ugain]
   Int_t           Ndata_sbs_gemFPP_hit_Utime;
   Double_t        sbs_gemFPP_hit_Utime[32];   //[Ndata.sbs.gemFPP.hit.Utime]
   Int_t           Ndata_sbs_gemFPP_hit_UtimeDeconv;
   Double_t        sbs_gemFPP_hit_UtimeDeconv[32];   //[Ndata.sbs.gemFPP.hit.UtimeDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_UtimeFit;
   Double_t        sbs_gemFPP_hit_UtimeFit[32];   //[Ndata.sbs.gemFPP.hit.UtimeFit]
   Int_t           Ndata_sbs_gemFPP_hit_UtimeMaxStrip;
   Double_t        sbs_gemFPP_hit_UtimeMaxStrip[32];   //[Ndata.sbs.gemFPP.hit.UtimeMaxStrip]
   Int_t           Ndata_sbs_gemFPP_hit_UtimeMaxStripDeconv;
   Double_t        sbs_gemFPP_hit_UtimeMaxStripDeconv[32];   //[Ndata.sbs.gemFPP.hit.UtimeMaxStripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_UtimeMaxStripFit;
   Double_t        sbs_gemFPP_hit_UtimeMaxStripFit[32];   //[Ndata.sbs.gemFPP.hit.UtimeMaxStripFit]
   Int_t           Ndata_sbs_gemFPP_hit_Vgain;
   Double_t        sbs_gemFPP_hit_Vgain[32];   //[Ndata.sbs.gemFPP.hit.Vgain]
   Int_t           Ndata_sbs_gemFPP_hit_Vtime;
   Double_t        sbs_gemFPP_hit_Vtime[32];   //[Ndata.sbs.gemFPP.hit.Vtime]
   Int_t           Ndata_sbs_gemFPP_hit_VtimeDeconv;
   Double_t        sbs_gemFPP_hit_VtimeDeconv[32];   //[Ndata.sbs.gemFPP.hit.VtimeDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_VtimeFit;
   Double_t        sbs_gemFPP_hit_VtimeFit[32];   //[Ndata.sbs.gemFPP.hit.VtimeFit]
   Int_t           Ndata_sbs_gemFPP_hit_VtimeMaxStrip;
   Double_t        sbs_gemFPP_hit_VtimeMaxStrip[32];   //[Ndata.sbs.gemFPP.hit.VtimeMaxStrip]
   Int_t           Ndata_sbs_gemFPP_hit_VtimeMaxStripDeconv;
   Double_t        sbs_gemFPP_hit_VtimeMaxStripDeconv[32];   //[Ndata.sbs.gemFPP.hit.VtimeMaxStripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_VtimeMaxStripFit;
   Double_t        sbs_gemFPP_hit_VtimeMaxStripFit[32];   //[Ndata.sbs.gemFPP.hit.VtimeMaxStripFit]
   Int_t           Ndata_sbs_gemFPP_hit_adc_id_U;
   Double_t        sbs_gemFPP_hit_adc_id_U[32];   //[Ndata.sbs.gemFPP.hit.adc_id_U]
   Int_t           Ndata_sbs_gemFPP_hit_adc_id_V;
   Double_t        sbs_gemFPP_hit_adc_id_V[32];   //[Ndata.sbs.gemFPP.hit.adc_id_V]
   Int_t           Ndata_sbs_gemFPP_hit_ccor_clust;
   Double_t        sbs_gemFPP_hit_ccor_clust[32];   //[Ndata.sbs.gemFPP.hit.ccor_clust]
   Int_t           Ndata_sbs_gemFPP_hit_ccor_clust_deconv;
   Double_t        sbs_gemFPP_hit_ccor_clust_deconv[32];   //[Ndata.sbs.gemFPP.hit.ccor_clust_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ccor_strip;
   Double_t        sbs_gemFPP_hit_ccor_strip[32];   //[Ndata.sbs.gemFPP.hit.ccor_strip]
   Int_t           Ndata_sbs_gemFPP_hit_ccor_strip_deconv;
   Double_t        sbs_gemFPP_hit_ccor_strip_deconv[32];   //[Ndata.sbs.gemFPP.hit.ccor_strip_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_crate_U;
   Double_t        sbs_gemFPP_hit_crate_U[32];   //[Ndata.sbs.gemFPP.hit.crate_U]
   Int_t           Ndata_sbs_gemFPP_hit_crate_V;
   Double_t        sbs_gemFPP_hit_crate_V[32];   //[Ndata.sbs.gemFPP.hit.crate_V]
   Int_t           Ndata_sbs_gemFPP_hit_deltat;
   Double_t        sbs_gemFPP_hit_deltat[32];   //[Ndata.sbs.gemFPP.hit.deltat]
   Int_t           Ndata_sbs_gemFPP_hit_deltat_deconv;
   Double_t        sbs_gemFPP_hit_deltat_deconv[32];   //[Ndata.sbs.gemFPP.hit.deltat_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_deltat_fit;
   Double_t        sbs_gemFPP_hit_deltat_fit[32];   //[Ndata.sbs.gemFPP.hit.deltat_fit]
   Int_t           Ndata_sbs_gemFPP_hit_eresidu;
   Double_t        sbs_gemFPP_hit_eresidu[32];   //[Ndata.sbs.gemFPP.hit.eresidu]
   Int_t           Ndata_sbs_gemFPP_hit_eresidv;
   Double_t        sbs_gemFPP_hit_eresidv[32];   //[Ndata.sbs.gemFPP.hit.eresidv]
   Int_t           Ndata_sbs_gemFPP_hit_icombomaxUclustDeconv;
   Double_t        sbs_gemFPP_hit_icombomaxUclustDeconv[32];   //[Ndata.sbs.gemFPP.hit.icombomaxUclustDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_icombomaxUstripDeconv;
   Double_t        sbs_gemFPP_hit_icombomaxUstripDeconv[32];   //[Ndata.sbs.gemFPP.hit.icombomaxUstripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_icombomaxVclustDeconv;
   Double_t        sbs_gemFPP_hit_icombomaxVclustDeconv[32];   //[Ndata.sbs.gemFPP.hit.icombomaxVclustDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_icombomaxVstripDeconv;
   Double_t        sbs_gemFPP_hit_icombomaxVstripDeconv[32];   //[Ndata.sbs.gemFPP.hit.icombomaxVstripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxUclust;
   Double_t        sbs_gemFPP_hit_isampmaxUclust[32];   //[Ndata.sbs.gemFPP.hit.isampmaxUclust]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxUclustDeconv;
   Double_t        sbs_gemFPP_hit_isampmaxUclustDeconv[32];   //[Ndata.sbs.gemFPP.hit.isampmaxUclustDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxUstrip;
   Double_t        sbs_gemFPP_hit_isampmaxUstrip[32];   //[Ndata.sbs.gemFPP.hit.isampmaxUstrip]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxUstripDeconv;
   Double_t        sbs_gemFPP_hit_isampmaxUstripDeconv[32];   //[Ndata.sbs.gemFPP.hit.isampmaxUstripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxVclust;
   Double_t        sbs_gemFPP_hit_isampmaxVclust[32];   //[Ndata.sbs.gemFPP.hit.isampmaxVclust]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxVclustDeconv;
   Double_t        sbs_gemFPP_hit_isampmaxVclustDeconv[32];   //[Ndata.sbs.gemFPP.hit.isampmaxVclustDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxVstrip;
   Double_t        sbs_gemFPP_hit_isampmaxVstrip[32];   //[Ndata.sbs.gemFPP.hit.isampmaxVstrip]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxVstripDeconv;
   Double_t        sbs_gemFPP_hit_isampmaxVstripDeconv[32];   //[Ndata.sbs.gemFPP.hit.isampmaxVstripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_layer;
   Double_t        sbs_gemFPP_hit_layer[32];   //[Ndata.sbs.gemFPP.hit.layer]
   Int_t           Ndata_sbs_gemFPP_hit_module;
   Double_t        sbs_gemFPP_hit_module[32];   //[Ndata.sbs.gemFPP.hit.module]
   Int_t           Ndata_sbs_gemFPP_hit_mpd_U;
   Double_t        sbs_gemFPP_hit_mpd_U[32];   //[Ndata.sbs.gemFPP.hit.mpd_U]
   Int_t           Ndata_sbs_gemFPP_hit_mpd_V;
   Double_t        sbs_gemFPP_hit_mpd_V[32];   //[Ndata.sbs.gemFPP.hit.mpd_V]
   Int_t           Ndata_sbs_gemFPP_hit_nstripu;
   Double_t        sbs_gemFPP_hit_nstripu[32];   //[Ndata.sbs.gemFPP.hit.nstripu]
   Int_t           Ndata_sbs_gemFPP_hit_nstripv;
   Double_t        sbs_gemFPP_hit_nstripv[32];   //[Ndata.sbs.gemFPP.hit.nstripv]
   Int_t           Ndata_sbs_gemFPP_hit_residu;
   Double_t        sbs_gemFPP_hit_residu[32];   //[Ndata.sbs.gemFPP.hit.residu]
   Int_t           Ndata_sbs_gemFPP_hit_residv;
   Double_t        sbs_gemFPP_hit_residv[32];   //[Ndata.sbs.gemFPP.hit.residv]
   Int_t           Ndata_sbs_gemFPP_hit_trackindex;
   Double_t        sbs_gemFPP_hit_trackindex[32];   //[Ndata.sbs.gemFPP.hit.trackindex]
   Int_t           Ndata_sbs_gemFPP_hit_u;
   Double_t        sbs_gemFPP_hit_u[32];   //[Ndata.sbs.gemFPP.hit.u]
   Int_t           Ndata_sbs_gemFPP_hit_umoment;
   Double_t        sbs_gemFPP_hit_umoment[32];   //[Ndata.sbs.gemFPP.hit.umoment]
   Int_t           Ndata_sbs_gemFPP_hit_usigma;
   Double_t        sbs_gemFPP_hit_usigma[32];   //[Ndata.sbs.gemFPP.hit.usigma]
   Int_t           Ndata_sbs_gemFPP_hit_ustriphi;
   Double_t        sbs_gemFPP_hit_ustriphi[32];   //[Ndata.sbs.gemFPP.hit.ustriphi]
   Int_t           Ndata_sbs_gemFPP_hit_ustriplo;
   Double_t        sbs_gemFPP_hit_ustriplo[32];   //[Ndata.sbs.gemFPP.hit.ustriplo]
   Int_t           Ndata_sbs_gemFPP_hit_ustripmax;
   Double_t        sbs_gemFPP_hit_ustripmax[32];   //[Ndata.sbs.gemFPP.hit.ustripmax]
   Int_t           Ndata_sbs_gemFPP_hit_v;
   Double_t        sbs_gemFPP_hit_v[32];   //[Ndata.sbs.gemFPP.hit.v]
   Int_t           Ndata_sbs_gemFPP_hit_vmoment;
   Double_t        sbs_gemFPP_hit_vmoment[32];   //[Ndata.sbs.gemFPP.hit.vmoment]
   Int_t           Ndata_sbs_gemFPP_hit_vsigma;
   Double_t        sbs_gemFPP_hit_vsigma[32];   //[Ndata.sbs.gemFPP.hit.vsigma]
   Int_t           Ndata_sbs_gemFPP_hit_vstriphi;
   Double_t        sbs_gemFPP_hit_vstriphi[32];   //[Ndata.sbs.gemFPP.hit.vstriphi]
   Int_t           Ndata_sbs_gemFPP_hit_vstriplo;
   Double_t        sbs_gemFPP_hit_vstriplo[32];   //[Ndata.sbs.gemFPP.hit.vstriplo]
   Int_t           Ndata_sbs_gemFPP_hit_vstripmax;
   Double_t        sbs_gemFPP_hit_vstripmax[32];   //[Ndata.sbs.gemFPP.hit.vstripmax]
   Int_t           Ndata_sbs_gemFPP_hit_xglobal;
   Double_t        sbs_gemFPP_hit_xglobal[32];   //[Ndata.sbs.gemFPP.hit.xglobal]
   Int_t           Ndata_sbs_gemFPP_hit_xlocal;
   Double_t        sbs_gemFPP_hit_xlocal[32];   //[Ndata.sbs.gemFPP.hit.xlocal]
   Int_t           Ndata_sbs_gemFPP_hit_yglobal;
   Double_t        sbs_gemFPP_hit_yglobal[32];   //[Ndata.sbs.gemFPP.hit.yglobal]
   Int_t           Ndata_sbs_gemFPP_hit_ylocal;
   Double_t        sbs_gemFPP_hit_ylocal[32];   //[Ndata.sbs.gemFPP.hit.ylocal]
   Int_t           Ndata_sbs_gemFPP_hit_zglobal;
   Double_t        sbs_gemFPP_hit_zglobal[32];   //[Ndata.sbs.gemFPP.hit.zglobal]
   Int_t           Ndata_sbs_gemFPP_n2Dhit_layer;
   Double_t        sbs_gemFPP_n2Dhit_layer[8];   //[Ndata.sbs.gemFPP.n2Dhit_layer]
   Int_t           Ndata_sbs_gemFPP_nclustu_layer;
   Double_t        sbs_gemFPP_nclustu_layer[8];   //[Ndata.sbs.gemFPP.nclustu_layer]
   Int_t           Ndata_sbs_gemFPP_nclustv_layer;
   Double_t        sbs_gemFPP_nclustv_layer[8];   //[Ndata.sbs.gemFPP.nclustv_layer]
   Int_t           Ndata_sbs_gemFPP_nstripsu_layer;
   Double_t        sbs_gemFPP_nstripsu_layer[8];   //[Ndata.sbs.gemFPP.nstripsu_layer]
   Int_t           Ndata_sbs_gemFPP_nstripsv_layer;
   Double_t        sbs_gemFPP_nstripsv_layer[8];   //[Ndata.sbs.gemFPP.nstripsv_layer]
   Int_t           Ndata_sbs_gemFPP_track_chi2ndf;
   Double_t        sbs_gemFPP_track_chi2ndf[8];   //[Ndata.sbs.gemFPP.track.chi2ndf]
   Int_t           Ndata_sbs_gemFPP_track_chi2ndf_hitquality;
   Double_t        sbs_gemFPP_track_chi2ndf_hitquality[8];   //[Ndata.sbs.gemFPP.track.chi2ndf_hitquality]
   Int_t           Ndata_sbs_gemFPP_track_ngoodhits;
   Double_t        sbs_gemFPP_track_ngoodhits[8];   //[Ndata.sbs.gemFPP.track.ngoodhits]
   Int_t           Ndata_sbs_gemFPP_track_nhits;
   Double_t        sbs_gemFPP_track_nhits[8];   //[Ndata.sbs.gemFPP.track.nhits]
   Int_t           Ndata_sbs_gemFPP_track_phi;
   Double_t        sbs_gemFPP_track_phi[8];   //[Ndata.sbs.gemFPP.track.phi]
   Int_t           Ndata_sbs_gemFPP_track_sclose;
   Double_t        sbs_gemFPP_track_sclose[8];   //[Ndata.sbs.gemFPP.track.sclose]
   Int_t           Ndata_sbs_gemFPP_track_t0;
   Double_t        sbs_gemFPP_track_t0[8];   //[Ndata.sbs.gemFPP.track.t0]
   Int_t           Ndata_sbs_gemFPP_track_theta;
   Double_t        sbs_gemFPP_track_theta[8];   //[Ndata.sbs.gemFPP.track.theta]
   Int_t           Ndata_sbs_gemFPP_track_x;
   Double_t        sbs_gemFPP_track_x[8];   //[Ndata.sbs.gemFPP.track.x]
   Int_t           Ndata_sbs_gemFPP_track_xp;
   Double_t        sbs_gemFPP_track_xp[8];   //[Ndata.sbs.gemFPP.track.xp]
   Int_t           Ndata_sbs_gemFPP_track_y;
   Double_t        sbs_gemFPP_track_y[8];   //[Ndata.sbs.gemFPP.track.y]
   Int_t           Ndata_sbs_gemFPP_track_yp;
   Double_t        sbs_gemFPP_track_yp[8];   //[Ndata.sbs.gemFPP.track.yp]
   Int_t           Ndata_sbs_gemFPP_track_zclose;
   Double_t        sbs_gemFPP_track_zclose[8];   //[Ndata.sbs.gemFPP.track.zclose]
   Int_t           Ndata_sbs_gemFT_hit_ADCU;
   Double_t        sbs_gemFT_hit_ADCU[29];   //[Ndata.sbs.gemFT.hit.ADCU]
   Int_t           Ndata_sbs_gemFT_hit_ADCU_deconv;
   Double_t        sbs_gemFT_hit_ADCU_deconv[29];   //[Ndata.sbs.gemFT.hit.ADCU_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCV;
   Double_t        sbs_gemFT_hit_ADCV[29];   //[Ndata.sbs.gemFT.hit.ADCV]
   Int_t           Ndata_sbs_gemFT_hit_ADCV_deconv;
   Double_t        sbs_gemFT_hit_ADCV_deconv[29];   //[Ndata.sbs.gemFT.hit.ADCV_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCasym;
   Double_t        sbs_gemFT_hit_ADCasym[29];   //[Ndata.sbs.gemFT.hit.ADCasym]
   Int_t           Ndata_sbs_gemFT_hit_ADCasym_deconv;
   Double_t        sbs_gemFT_hit_ADCasym_deconv[29];   //[Ndata.sbs.gemFT.hit.ADCasym_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCavg;
   Double_t        sbs_gemFT_hit_ADCavg[29];   //[Ndata.sbs.gemFT.hit.ADCavg]
   Int_t           Ndata_sbs_gemFT_hit_ADCavg_deconv;
   Double_t        sbs_gemFT_hit_ADCavg_deconv[29];   //[Ndata.sbs.gemFT.hit.ADCavg_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac0_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac0_Umax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac0_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac0_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac0_Vmax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac0_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac1_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac1_Umax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac1_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac1_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac1_Vmax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac1_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac2_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac2_Umax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac2_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac2_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac2_Vmax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac2_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac3_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac3_Umax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac3_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac3_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac3_Vmax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac3_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac4_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac4_Umax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac4_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac4_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac4_Vmax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac4_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac5_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac5_Umax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac5_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac5_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac5_Vmax[29];   //[Ndata.sbs.gemFT.hit.ADCfrac5_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxcomboUclust_deconv;
   Double_t        sbs_gemFT_hit_ADCmaxcomboUclust_deconv[29];   //[Ndata.sbs.gemFT.hit.ADCmaxcomboUclust_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxcomboVclust_deconv;
   Double_t        sbs_gemFT_hit_ADCmaxcomboVclust_deconv[29];   //[Ndata.sbs.gemFT.hit.ADCmaxcomboVclust_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampU;
   Double_t        sbs_gemFT_hit_ADCmaxsampU[29];   //[Ndata.sbs.gemFT.hit.ADCmaxsampU]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampUclust;
   Double_t        sbs_gemFT_hit_ADCmaxsampUclust[29];   //[Ndata.sbs.gemFT.hit.ADCmaxsampUclust]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampUclust_deconv;
   Double_t        sbs_gemFT_hit_ADCmaxsampUclust_deconv[29];   //[Ndata.sbs.gemFT.hit.ADCmaxsampUclust_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampV;
   Double_t        sbs_gemFT_hit_ADCmaxsampV[29];   //[Ndata.sbs.gemFT.hit.ADCmaxsampV]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampVclust;
   Double_t        sbs_gemFT_hit_ADCmaxsampVclust[29];   //[Ndata.sbs.gemFT.hit.ADCmaxsampVclust]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampVclust_deconv;
   Double_t        sbs_gemFT_hit_ADCmaxsampVclust_deconv[29];   //[Ndata.sbs.gemFT.hit.ADCmaxsampVclust_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxstripU;
   Double_t        sbs_gemFT_hit_ADCmaxstripU[29];   //[Ndata.sbs.gemFT.hit.ADCmaxstripU]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxstripV;
   Double_t        sbs_gemFT_hit_ADCmaxstripV[29];   //[Ndata.sbs.gemFT.hit.ADCmaxstripV]
   Int_t           Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U;
   Double_t        sbs_gemFT_hit_BUILD_ALL_SAMPLES_U[29];   //[Ndata.sbs.gemFT.hit.BUILD_ALL_SAMPLES_U]
   Int_t           Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V;
   Double_t        sbs_gemFT_hit_BUILD_ALL_SAMPLES_V[29];   //[Ndata.sbs.gemFT.hit.BUILD_ALL_SAMPLES_V]
   Int_t           Ndata_sbs_gemFT_hit_CM_GOOD_U;
   Double_t        sbs_gemFT_hit_CM_GOOD_U[29];   //[Ndata.sbs.gemFT.hit.CM_GOOD_U]
   Int_t           Ndata_sbs_gemFT_hit_CM_GOOD_V;
   Double_t        sbs_gemFT_hit_CM_GOOD_V[29];   //[Ndata.sbs.gemFT.hit.CM_GOOD_V]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC0_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC0_Umax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC0_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC0_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC0_Vmax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC0_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC1_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC1_Umax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC1_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC1_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC1_Vmax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC1_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC2_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC2_Umax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC2_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC2_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC2_Vmax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC2_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC3_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC3_Umax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC3_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC3_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC3_Vmax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC3_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC4_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC4_Umax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC4_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC4_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC4_Vmax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC4_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC5_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC5_Umax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC5_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC5_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC5_Vmax[29];   //[Ndata.sbs.gemFT.hit.DeconvADC5_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxcomboU;
   Double_t        sbs_gemFT_hit_DeconvADCmaxcomboU[29];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxcomboU]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxcomboV;
   Double_t        sbs_gemFT_hit_DeconvADCmaxcomboV[29];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxcomboV]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxsampU;
   Double_t        sbs_gemFT_hit_DeconvADCmaxsampU[29];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxsampU]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxsampV;
   Double_t        sbs_gemFT_hit_DeconvADCmaxsampV[29];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxsampV]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxstripU;
   Double_t        sbs_gemFT_hit_DeconvADCmaxstripU[29];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxstripU]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxstripV;
   Double_t        sbs_gemFT_hit_DeconvADCmaxstripV[29];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxstripV]
   Int_t           Ndata_sbs_gemFT_hit_ENABLE_CM_U;
   Double_t        sbs_gemFT_hit_ENABLE_CM_U[29];   //[Ndata.sbs.gemFT.hit.ENABLE_CM_U]
   Int_t           Ndata_sbs_gemFT_hit_ENABLE_CM_V;
   Double_t        sbs_gemFT_hit_ENABLE_CM_V[29];   //[Ndata.sbs.gemFT.hit.ENABLE_CM_V]
   Int_t           Ndata_sbs_gemFT_hit_TSchi2_Umax;
   Double_t        sbs_gemFT_hit_TSchi2_Umax[29];   //[Ndata.sbs.gemFT.hit.TSchi2_Umax]
   Int_t           Ndata_sbs_gemFT_hit_TSchi2_Vmax;
   Double_t        sbs_gemFT_hit_TSchi2_Vmax[29];   //[Ndata.sbs.gemFT.hit.TSchi2_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_TSprob_Umax;
   Double_t        sbs_gemFT_hit_TSprob_Umax[29];   //[Ndata.sbs.gemFT.hit.TSprob_Umax]
   Int_t           Ndata_sbs_gemFT_hit_TSprob_Vmax;
   Double_t        sbs_gemFT_hit_TSprob_Vmax[29];   //[Ndata.sbs.gemFT.hit.TSprob_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_Tavg;
   Double_t        sbs_gemFT_hit_Tavg[29];   //[Ndata.sbs.gemFT.hit.Tavg]
   Int_t           Ndata_sbs_gemFT_hit_Tavg_corr;
   Double_t        sbs_gemFT_hit_Tavg_corr[29];   //[Ndata.sbs.gemFT.hit.Tavg_corr]
   Int_t           Ndata_sbs_gemFT_hit_Tavg_deconv;
   Double_t        sbs_gemFT_hit_Tavg_deconv[29];   //[Ndata.sbs.gemFT.hit.Tavg_deconv]
   Int_t           Ndata_sbs_gemFT_hit_Tavg_fit;
   Double_t        sbs_gemFT_hit_Tavg_fit[29];   //[Ndata.sbs.gemFT.hit.Tavg_fit]
   Int_t           Ndata_sbs_gemFT_hit_Ugain;
   Double_t        sbs_gemFT_hit_Ugain[29];   //[Ndata.sbs.gemFT.hit.Ugain]
   Int_t           Ndata_sbs_gemFT_hit_Utime;
   Double_t        sbs_gemFT_hit_Utime[29];   //[Ndata.sbs.gemFT.hit.Utime]
   Int_t           Ndata_sbs_gemFT_hit_UtimeDeconv;
   Double_t        sbs_gemFT_hit_UtimeDeconv[29];   //[Ndata.sbs.gemFT.hit.UtimeDeconv]
   Int_t           Ndata_sbs_gemFT_hit_UtimeFit;
   Double_t        sbs_gemFT_hit_UtimeFit[29];   //[Ndata.sbs.gemFT.hit.UtimeFit]
   Int_t           Ndata_sbs_gemFT_hit_UtimeMaxStrip;
   Double_t        sbs_gemFT_hit_UtimeMaxStrip[29];   //[Ndata.sbs.gemFT.hit.UtimeMaxStrip]
   Int_t           Ndata_sbs_gemFT_hit_UtimeMaxStripDeconv;
   Double_t        sbs_gemFT_hit_UtimeMaxStripDeconv[29];   //[Ndata.sbs.gemFT.hit.UtimeMaxStripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_UtimeMaxStripFit;
   Double_t        sbs_gemFT_hit_UtimeMaxStripFit[29];   //[Ndata.sbs.gemFT.hit.UtimeMaxStripFit]
   Int_t           Ndata_sbs_gemFT_hit_Vgain;
   Double_t        sbs_gemFT_hit_Vgain[29];   //[Ndata.sbs.gemFT.hit.Vgain]
   Int_t           Ndata_sbs_gemFT_hit_Vtime;
   Double_t        sbs_gemFT_hit_Vtime[29];   //[Ndata.sbs.gemFT.hit.Vtime]
   Int_t           Ndata_sbs_gemFT_hit_VtimeDeconv;
   Double_t        sbs_gemFT_hit_VtimeDeconv[29];   //[Ndata.sbs.gemFT.hit.VtimeDeconv]
   Int_t           Ndata_sbs_gemFT_hit_VtimeFit;
   Double_t        sbs_gemFT_hit_VtimeFit[29];   //[Ndata.sbs.gemFT.hit.VtimeFit]
   Int_t           Ndata_sbs_gemFT_hit_VtimeMaxStrip;
   Double_t        sbs_gemFT_hit_VtimeMaxStrip[29];   //[Ndata.sbs.gemFT.hit.VtimeMaxStrip]
   Int_t           Ndata_sbs_gemFT_hit_VtimeMaxStripDeconv;
   Double_t        sbs_gemFT_hit_VtimeMaxStripDeconv[29];   //[Ndata.sbs.gemFT.hit.VtimeMaxStripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_VtimeMaxStripFit;
   Double_t        sbs_gemFT_hit_VtimeMaxStripFit[29];   //[Ndata.sbs.gemFT.hit.VtimeMaxStripFit]
   Int_t           Ndata_sbs_gemFT_hit_adc_id_U;
   Double_t        sbs_gemFT_hit_adc_id_U[29];   //[Ndata.sbs.gemFT.hit.adc_id_U]
   Int_t           Ndata_sbs_gemFT_hit_adc_id_V;
   Double_t        sbs_gemFT_hit_adc_id_V[29];   //[Ndata.sbs.gemFT.hit.adc_id_V]
   Int_t           Ndata_sbs_gemFT_hit_ccor_clust;
   Double_t        sbs_gemFT_hit_ccor_clust[29];   //[Ndata.sbs.gemFT.hit.ccor_clust]
   Int_t           Ndata_sbs_gemFT_hit_ccor_clust_deconv;
   Double_t        sbs_gemFT_hit_ccor_clust_deconv[29];   //[Ndata.sbs.gemFT.hit.ccor_clust_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ccor_strip;
   Double_t        sbs_gemFT_hit_ccor_strip[29];   //[Ndata.sbs.gemFT.hit.ccor_strip]
   Int_t           Ndata_sbs_gemFT_hit_ccor_strip_deconv;
   Double_t        sbs_gemFT_hit_ccor_strip_deconv[29];   //[Ndata.sbs.gemFT.hit.ccor_strip_deconv]
   Int_t           Ndata_sbs_gemFT_hit_crate_U;
   Double_t        sbs_gemFT_hit_crate_U[29];   //[Ndata.sbs.gemFT.hit.crate_U]
   Int_t           Ndata_sbs_gemFT_hit_crate_V;
   Double_t        sbs_gemFT_hit_crate_V[29];   //[Ndata.sbs.gemFT.hit.crate_V]
   Int_t           Ndata_sbs_gemFT_hit_deltat;
   Double_t        sbs_gemFT_hit_deltat[29];   //[Ndata.sbs.gemFT.hit.deltat]
   Int_t           Ndata_sbs_gemFT_hit_deltat_deconv;
   Double_t        sbs_gemFT_hit_deltat_deconv[29];   //[Ndata.sbs.gemFT.hit.deltat_deconv]
   Int_t           Ndata_sbs_gemFT_hit_deltat_fit;
   Double_t        sbs_gemFT_hit_deltat_fit[29];   //[Ndata.sbs.gemFT.hit.deltat_fit]
   Int_t           Ndata_sbs_gemFT_hit_eresidu;
   Double_t        sbs_gemFT_hit_eresidu[29];   //[Ndata.sbs.gemFT.hit.eresidu]
   Int_t           Ndata_sbs_gemFT_hit_eresidv;
   Double_t        sbs_gemFT_hit_eresidv[29];   //[Ndata.sbs.gemFT.hit.eresidv]
   Int_t           Ndata_sbs_gemFT_hit_icombomaxUclustDeconv;
   Double_t        sbs_gemFT_hit_icombomaxUclustDeconv[29];   //[Ndata.sbs.gemFT.hit.icombomaxUclustDeconv]
   Int_t           Ndata_sbs_gemFT_hit_icombomaxUstripDeconv;
   Double_t        sbs_gemFT_hit_icombomaxUstripDeconv[29];   //[Ndata.sbs.gemFT.hit.icombomaxUstripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_icombomaxVclustDeconv;
   Double_t        sbs_gemFT_hit_icombomaxVclustDeconv[29];   //[Ndata.sbs.gemFT.hit.icombomaxVclustDeconv]
   Int_t           Ndata_sbs_gemFT_hit_icombomaxVstripDeconv;
   Double_t        sbs_gemFT_hit_icombomaxVstripDeconv[29];   //[Ndata.sbs.gemFT.hit.icombomaxVstripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxUclust;
   Double_t        sbs_gemFT_hit_isampmaxUclust[29];   //[Ndata.sbs.gemFT.hit.isampmaxUclust]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxUclustDeconv;
   Double_t        sbs_gemFT_hit_isampmaxUclustDeconv[29];   //[Ndata.sbs.gemFT.hit.isampmaxUclustDeconv]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxUstrip;
   Double_t        sbs_gemFT_hit_isampmaxUstrip[29];   //[Ndata.sbs.gemFT.hit.isampmaxUstrip]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxUstripDeconv;
   Double_t        sbs_gemFT_hit_isampmaxUstripDeconv[29];   //[Ndata.sbs.gemFT.hit.isampmaxUstripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxVclust;
   Double_t        sbs_gemFT_hit_isampmaxVclust[29];   //[Ndata.sbs.gemFT.hit.isampmaxVclust]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxVclustDeconv;
   Double_t        sbs_gemFT_hit_isampmaxVclustDeconv[29];   //[Ndata.sbs.gemFT.hit.isampmaxVclustDeconv]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxVstrip;
   Double_t        sbs_gemFT_hit_isampmaxVstrip[29];   //[Ndata.sbs.gemFT.hit.isampmaxVstrip]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxVstripDeconv;
   Double_t        sbs_gemFT_hit_isampmaxVstripDeconv[29];   //[Ndata.sbs.gemFT.hit.isampmaxVstripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_layer;
   Double_t        sbs_gemFT_hit_layer[29];   //[Ndata.sbs.gemFT.hit.layer]
   Int_t           Ndata_sbs_gemFT_hit_module;
   Double_t        sbs_gemFT_hit_module[29];   //[Ndata.sbs.gemFT.hit.module]
   Int_t           Ndata_sbs_gemFT_hit_mpd_U;
   Double_t        sbs_gemFT_hit_mpd_U[29];   //[Ndata.sbs.gemFT.hit.mpd_U]
   Int_t           Ndata_sbs_gemFT_hit_mpd_V;
   Double_t        sbs_gemFT_hit_mpd_V[29];   //[Ndata.sbs.gemFT.hit.mpd_V]
   Int_t           Ndata_sbs_gemFT_hit_nstripu;
   Double_t        sbs_gemFT_hit_nstripu[29];   //[Ndata.sbs.gemFT.hit.nstripu]
   Int_t           Ndata_sbs_gemFT_hit_nstripv;
   Double_t        sbs_gemFT_hit_nstripv[29];   //[Ndata.sbs.gemFT.hit.nstripv]
   Int_t           Ndata_sbs_gemFT_hit_residu;
   Double_t        sbs_gemFT_hit_residu[29];   //[Ndata.sbs.gemFT.hit.residu]
   Int_t           Ndata_sbs_gemFT_hit_residv;
   Double_t        sbs_gemFT_hit_residv[29];   //[Ndata.sbs.gemFT.hit.residv]
   Int_t           Ndata_sbs_gemFT_hit_trackindex;
   Double_t        sbs_gemFT_hit_trackindex[29];   //[Ndata.sbs.gemFT.hit.trackindex]
   Int_t           Ndata_sbs_gemFT_hit_u;
   Double_t        sbs_gemFT_hit_u[29];   //[Ndata.sbs.gemFT.hit.u]
   Int_t           Ndata_sbs_gemFT_hit_umoment;
   Double_t        sbs_gemFT_hit_umoment[29];   //[Ndata.sbs.gemFT.hit.umoment]
   Int_t           Ndata_sbs_gemFT_hit_usigma;
   Double_t        sbs_gemFT_hit_usigma[29];   //[Ndata.sbs.gemFT.hit.usigma]
   Int_t           Ndata_sbs_gemFT_hit_ustriphi;
   Double_t        sbs_gemFT_hit_ustriphi[29];   //[Ndata.sbs.gemFT.hit.ustriphi]
   Int_t           Ndata_sbs_gemFT_hit_ustriplo;
   Double_t        sbs_gemFT_hit_ustriplo[29];   //[Ndata.sbs.gemFT.hit.ustriplo]
   Int_t           Ndata_sbs_gemFT_hit_ustripmax;
   Double_t        sbs_gemFT_hit_ustripmax[29];   //[Ndata.sbs.gemFT.hit.ustripmax]
   Int_t           Ndata_sbs_gemFT_hit_v;
   Double_t        sbs_gemFT_hit_v[29];   //[Ndata.sbs.gemFT.hit.v]
   Int_t           Ndata_sbs_gemFT_hit_vmoment;
   Double_t        sbs_gemFT_hit_vmoment[29];   //[Ndata.sbs.gemFT.hit.vmoment]
   Int_t           Ndata_sbs_gemFT_hit_vsigma;
   Double_t        sbs_gemFT_hit_vsigma[29];   //[Ndata.sbs.gemFT.hit.vsigma]
   Int_t           Ndata_sbs_gemFT_hit_vstriphi;
   Double_t        sbs_gemFT_hit_vstriphi[29];   //[Ndata.sbs.gemFT.hit.vstriphi]
   Int_t           Ndata_sbs_gemFT_hit_vstriplo;
   Double_t        sbs_gemFT_hit_vstriplo[29];   //[Ndata.sbs.gemFT.hit.vstriplo]
   Int_t           Ndata_sbs_gemFT_hit_vstripmax;
   Double_t        sbs_gemFT_hit_vstripmax[29];   //[Ndata.sbs.gemFT.hit.vstripmax]
   Int_t           Ndata_sbs_gemFT_hit_xglobal;
   Double_t        sbs_gemFT_hit_xglobal[29];   //[Ndata.sbs.gemFT.hit.xglobal]
   Int_t           Ndata_sbs_gemFT_hit_xlocal;
   Double_t        sbs_gemFT_hit_xlocal[29];   //[Ndata.sbs.gemFT.hit.xlocal]
   Int_t           Ndata_sbs_gemFT_hit_yglobal;
   Double_t        sbs_gemFT_hit_yglobal[29];   //[Ndata.sbs.gemFT.hit.yglobal]
   Int_t           Ndata_sbs_gemFT_hit_ylocal;
   Double_t        sbs_gemFT_hit_ylocal[29];   //[Ndata.sbs.gemFT.hit.ylocal]
   Int_t           Ndata_sbs_gemFT_hit_zglobal;
   Double_t        sbs_gemFT_hit_zglobal[29];   //[Ndata.sbs.gemFT.hit.zglobal]
   Int_t           Ndata_sbs_gemFT_n2Dhit_layer;
   Double_t        sbs_gemFT_n2Dhit_layer[8];   //[Ndata.sbs.gemFT.n2Dhit_layer]
   Int_t           Ndata_sbs_gemFT_nclustu_layer;
   Double_t        sbs_gemFT_nclustu_layer[8];   //[Ndata.sbs.gemFT.nclustu_layer]
   Int_t           Ndata_sbs_gemFT_nclustv_layer;
   Double_t        sbs_gemFT_nclustv_layer[8];   //[Ndata.sbs.gemFT.nclustv_layer]
   Int_t           Ndata_sbs_gemFT_nstripsu_layer;
   Double_t        sbs_gemFT_nstripsu_layer[8];   //[Ndata.sbs.gemFT.nstripsu_layer]
   Int_t           Ndata_sbs_gemFT_nstripsv_layer;
   Double_t        sbs_gemFT_nstripsv_layer[8];   //[Ndata.sbs.gemFT.nstripsv_layer]
   Int_t           Ndata_sbs_gemFT_track_chi2ndf;
   Double_t        sbs_gemFT_track_chi2ndf[7];   //[Ndata.sbs.gemFT.track.chi2ndf]
   Int_t           Ndata_sbs_gemFT_track_chi2ndf_hitquality;
   Double_t        sbs_gemFT_track_chi2ndf_hitquality[7];   //[Ndata.sbs.gemFT.track.chi2ndf_hitquality]
   Int_t           Ndata_sbs_gemFT_track_ngoodhits;
   Double_t        sbs_gemFT_track_ngoodhits[7];   //[Ndata.sbs.gemFT.track.ngoodhits]
   Int_t           Ndata_sbs_gemFT_track_nhits;
   Double_t        sbs_gemFT_track_nhits[7];   //[Ndata.sbs.gemFT.track.nhits]
   Int_t           Ndata_sbs_gemFT_track_t0;
   Double_t        sbs_gemFT_track_t0[7];   //[Ndata.sbs.gemFT.track.t0]
   Int_t           Ndata_sbs_gemFT_track_x;
   Double_t        sbs_gemFT_track_x[7];   //[Ndata.sbs.gemFT.track.x]
   Int_t           Ndata_sbs_gemFT_track_xp;
   Double_t        sbs_gemFT_track_xp[7];   //[Ndata.sbs.gemFT.track.xp]
   Int_t           Ndata_sbs_gemFT_track_y;
   Double_t        sbs_gemFT_track_y[7];   //[Ndata.sbs.gemFT.track.y]
   Int_t           Ndata_sbs_gemFT_track_yp;
   Double_t        sbs_gemFT_track_yp[7];   //[Ndata.sbs.gemFT.track.yp]
   Int_t           Ndata_sbs_hcal_Ref_tdc;
   Double_t        sbs_hcal_Ref_tdc[4];   //[Ndata.sbs.hcal.Ref.tdc]
   Int_t           Ndata_sbs_hcal_Ref_tdc_mult;
   Double_t        sbs_hcal_Ref_tdc_mult[4];   //[Ndata.sbs.hcal.Ref.tdc_mult]
   Int_t           Ndata_sbs_hcal_Ref_tdcelemID;
   Double_t        sbs_hcal_Ref_tdcelemID[4];   //[Ndata.sbs.hcal.Ref.tdcelemID]
   Int_t           Ndata_sbs_hcal_clus_adctime;
   Double_t        sbs_hcal_clus_adctime[8];   //[Ndata.sbs.hcal.clus.adctime]
   Int_t           Ndata_sbs_hcal_clus_again;
   Double_t        sbs_hcal_clus_again[8];   //[Ndata.sbs.hcal.clus.again]
   Int_t           Ndata_sbs_hcal_clus_atimeblk;
   Double_t        sbs_hcal_clus_atimeblk[8];   //[Ndata.sbs.hcal.clus.atimeblk]
   Int_t           Ndata_sbs_hcal_clus_col;
   Double_t        sbs_hcal_clus_col[8];   //[Ndata.sbs.hcal.clus.col]
   Int_t           Ndata_sbs_hcal_clus_col_goodtdc;
   Double_t        sbs_hcal_clus_col_goodtdc[8];   //[Ndata.sbs.hcal.clus.col_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_e;
   Double_t        sbs_hcal_clus_e[8];   //[Ndata.sbs.hcal.clus.e]
   Int_t           Ndata_sbs_hcal_clus_e_goodtdc;
   Double_t        sbs_hcal_clus_e_goodtdc[8];   //[Ndata.sbs.hcal.clus.e_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_eblk;
   Double_t        sbs_hcal_clus_eblk[8];   //[Ndata.sbs.hcal.clus.eblk]
   Int_t           Ndata_sbs_hcal_clus_eblk_goodtdc;
   Double_t        sbs_hcal_clus_eblk_goodtdc[8];   //[Ndata.sbs.hcal.clus.eblk_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_id;
   Double_t        sbs_hcal_clus_id[8];   //[Ndata.sbs.hcal.clus.id]
   Int_t           Ndata_sbs_hcal_clus_id_goodtdc;
   Double_t        sbs_hcal_clus_id_goodtdc[8];   //[Ndata.sbs.hcal.clus.id_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_nblk;
   Double_t        sbs_hcal_clus_nblk[8];   //[Ndata.sbs.hcal.clus.nblk]
   Int_t           Ndata_sbs_hcal_clus_nblk_goodtdc;
   Double_t        sbs_hcal_clus_nblk_goodtdc[8];   //[Ndata.sbs.hcal.clus.nblk_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_row;
   Double_t        sbs_hcal_clus_row[8];   //[Ndata.sbs.hcal.clus.row]
   Int_t           Ndata_sbs_hcal_clus_row_goodtdc;
   Double_t        sbs_hcal_clus_row_goodtdc[8];   //[Ndata.sbs.hcal.clus.row_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_tdctime;
   Double_t        sbs_hcal_clus_tdctime[8];   //[Ndata.sbs.hcal.clus.tdctime]
   Int_t           Ndata_sbs_hcal_clus_tdctime_tw;
   Double_t        sbs_hcal_clus_tdctime_tw[8];   //[Ndata.sbs.hcal.clus.tdctime_tw]
   Int_t           Ndata_sbs_hcal_clus_tdctimeblk;
   Double_t        sbs_hcal_clus_tdctimeblk[8];   //[Ndata.sbs.hcal.clus.tdctimeblk]
   Int_t           Ndata_sbs_hcal_clus_tdctimeblk_tw;
   Double_t        sbs_hcal_clus_tdctimeblk_tw[8];   //[Ndata.sbs.hcal.clus.tdctimeblk_tw]
   Int_t           Ndata_sbs_hcal_clus_x;
   Double_t        sbs_hcal_clus_x[8];   //[Ndata.sbs.hcal.clus.x]
   Int_t           Ndata_sbs_hcal_clus_y;
   Double_t        sbs_hcal_clus_y[8];   //[Ndata.sbs.hcal.clus.y]
   Int_t           Ndata_sbs_hcal_clus_blk_again;
   Double_t        sbs_hcal_clus_blk_again[19];   //[Ndata.sbs.hcal.clus_blk.again]
   Int_t           Ndata_sbs_hcal_clus_blk_atime;
   Double_t        sbs_hcal_clus_blk_atime[19];   //[Ndata.sbs.hcal.clus_blk.atime]
   Int_t           Ndata_sbs_hcal_clus_blk_col;
   Double_t        sbs_hcal_clus_blk_col[19];   //[Ndata.sbs.hcal.clus_blk.col]
   Int_t           Ndata_sbs_hcal_clus_blk_e;
   Double_t        sbs_hcal_clus_blk_e[19];   //[Ndata.sbs.hcal.clus_blk.e]
   Int_t           Ndata_sbs_hcal_clus_blk_id;
   Double_t        sbs_hcal_clus_blk_id[19];   //[Ndata.sbs.hcal.clus_blk.id]
   Int_t           Ndata_sbs_hcal_clus_blk_row;
   Double_t        sbs_hcal_clus_blk_row[19];   //[Ndata.sbs.hcal.clus_blk.row]
   Int_t           Ndata_sbs_hcal_clus_blk_tdctime;
   Double_t        sbs_hcal_clus_blk_tdctime[19];   //[Ndata.sbs.hcal.clus_blk.tdctime]
   Int_t           Ndata_sbs_hcal_clus_blk_tdctime_tw;
   Double_t        sbs_hcal_clus_blk_tdctime_tw[19];   //[Ndata.sbs.hcal.clus_blk.tdctime_tw]
   Int_t           Ndata_sbs_hcal_clus_blk_x;
   Double_t        sbs_hcal_clus_blk_x[19];   //[Ndata.sbs.hcal.clus_blk.x]
   Int_t           Ndata_sbs_hcal_clus_blk_y;
   Double_t        sbs_hcal_clus_blk_y[19];   //[Ndata.sbs.hcal.clus_blk.y]
   Int_t           Ndata_sbs_hcal_goodblock_atime;
   Double_t        sbs_hcal_goodblock_atime[38];   //[Ndata.sbs.hcal.goodblock.atime]
   Int_t           Ndata_sbs_hcal_goodblock_cid;
   Double_t        sbs_hcal_goodblock_cid[38];   //[Ndata.sbs.hcal.goodblock.cid]
   Int_t           Ndata_sbs_hcal_goodblock_col;
   Double_t        sbs_hcal_goodblock_col[38];   //[Ndata.sbs.hcal.goodblock.col]
   Int_t           Ndata_sbs_hcal_goodblock_e;
   Double_t        sbs_hcal_goodblock_e[38];   //[Ndata.sbs.hcal.goodblock.e]
   Int_t           Ndata_sbs_hcal_goodblock_id;
   Double_t        sbs_hcal_goodblock_id[38];   //[Ndata.sbs.hcal.goodblock.id]
   Int_t           Ndata_sbs_hcal_goodblock_row;
   Double_t        sbs_hcal_goodblock_row[38];   //[Ndata.sbs.hcal.goodblock.row]
   Int_t           Ndata_sbs_hcal_goodblock_tdctime;
   Double_t        sbs_hcal_goodblock_tdctime[38];   //[Ndata.sbs.hcal.goodblock.tdctime]
   Int_t           Ndata_sbs_hcal_goodblock_tdctime_tw;
   Double_t        sbs_hcal_goodblock_tdctime_tw[38];   //[Ndata.sbs.hcal.goodblock.tdctime_tw]
   Int_t           Ndata_sbs_hcal_goodblock_x;
   Double_t        sbs_hcal_goodblock_x[38];   //[Ndata.sbs.hcal.goodblock.x]
   Int_t           Ndata_sbs_hcal_goodblock_y;
   Double_t        sbs_hcal_goodblock_y[38];   //[Ndata.sbs.hcal.goodblock.y]
   Int_t           Ndata_sbs_hcal_vtp_clus_e;
   Double_t        sbs_hcal_vtp_clus_e[8];   //[Ndata.sbs.hcal.vtp.clus.e]
   Int_t           Ndata_sbs_hcal_vtp_clus_size;
   Double_t        sbs_hcal_vtp_clus_size[8];   //[Ndata.sbs.hcal.vtp.clus.size]
   Int_t           Ndata_sbs_hcal_vtp_clus_time;
   Double_t        sbs_hcal_vtp_clus_time[8];   //[Ndata.sbs.hcal.vtp.clus.time]
   Int_t           Ndata_sbs_hcal_vtp_clus_x;
   Double_t        sbs_hcal_vtp_clus_x[8];   //[Ndata.sbs.hcal.vtp.clus.x]
   Int_t           Ndata_sbs_hcal_vtp_clus_y;
   Double_t        sbs_hcal_vtp_clus_y[8];   //[Ndata.sbs.hcal.vtp.clus.y]
   Int_t           Ndata_sbs_tr_beta;
   Double_t        sbs_tr_beta[7];   //[Ndata.sbs.tr.beta]
   Int_t           Ndata_sbs_tr_chi2;
   Double_t        sbs_tr_chi2[7];   //[Ndata.sbs.tr.chi2]
   Int_t           Ndata_sbs_tr_d_ph;
   Double_t        sbs_tr_d_ph[7];   //[Ndata.sbs.tr.d_ph]
   Int_t           Ndata_sbs_tr_d_th;
   Double_t        sbs_tr_d_th[7];   //[Ndata.sbs.tr.d_th]
   Int_t           Ndata_sbs_tr_d_x;
   Double_t        sbs_tr_d_x[7];   //[Ndata.sbs.tr.d_x]
   Int_t           Ndata_sbs_tr_d_y;
   Double_t        sbs_tr_d_y[7];   //[Ndata.sbs.tr.d_y]
   Int_t           Ndata_sbs_tr_dbeta;
   Double_t        sbs_tr_dbeta[7];   //[Ndata.sbs.tr.dbeta]
   Int_t           Ndata_sbs_tr_dtime;
   Double_t        sbs_tr_dtime[7];   //[Ndata.sbs.tr.dtime]
   Int_t           Ndata_sbs_tr_flag;
   Double_t        sbs_tr_flag[7];   //[Ndata.sbs.tr.flag]
   Int_t           Ndata_sbs_tr_ndof;
   Double_t        sbs_tr_ndof[7];   //[Ndata.sbs.tr.ndof]
   Int_t           Ndata_sbs_tr_p;
   Double_t        sbs_tr_p[7];   //[Ndata.sbs.tr.p]
   Int_t           Ndata_sbs_tr_pathl;
   Double_t        sbs_tr_pathl[7];   //[Ndata.sbs.tr.pathl]
   Int_t           Ndata_sbs_tr_ph;
   Double_t        sbs_tr_ph[7];   //[Ndata.sbs.tr.ph]
   Int_t           Ndata_sbs_tr_px;
   Double_t        sbs_tr_px[7];   //[Ndata.sbs.tr.px]
   Int_t           Ndata_sbs_tr_py;
   Double_t        sbs_tr_py[7];   //[Ndata.sbs.tr.py]
   Int_t           Ndata_sbs_tr_pz;
   Double_t        sbs_tr_pz[7];   //[Ndata.sbs.tr.pz]
   Int_t           Ndata_sbs_tr_r_ph;
   Double_t        sbs_tr_r_ph[7];   //[Ndata.sbs.tr.r_ph]
   Int_t           Ndata_sbs_tr_r_th;
   Double_t        sbs_tr_r_th[7];   //[Ndata.sbs.tr.r_th]
   Int_t           Ndata_sbs_tr_r_x;
   Double_t        sbs_tr_r_x[7];   //[Ndata.sbs.tr.r_x]
   Int_t           Ndata_sbs_tr_r_y;
   Double_t        sbs_tr_r_y[7];   //[Ndata.sbs.tr.r_y]
   Int_t           Ndata_sbs_tr_tg_dp;
   Double_t        sbs_tr_tg_dp[7];   //[Ndata.sbs.tr.tg_dp]
   Int_t           Ndata_sbs_tr_tg_ph;
   Double_t        sbs_tr_tg_ph[7];   //[Ndata.sbs.tr.tg_ph]
   Int_t           Ndata_sbs_tr_tg_th;
   Double_t        sbs_tr_tg_th[7];   //[Ndata.sbs.tr.tg_th]
   Int_t           Ndata_sbs_tr_tg_y;
   Double_t        sbs_tr_tg_y[7];   //[Ndata.sbs.tr.tg_y]
   Int_t           Ndata_sbs_tr_th;
   Double_t        sbs_tr_th[7];   //[Ndata.sbs.tr.th]
   Int_t           Ndata_sbs_tr_time;
   Double_t        sbs_tr_time[7];   //[Ndata.sbs.tr.time]
   Int_t           Ndata_sbs_tr_vx;
   Double_t        sbs_tr_vx[7];   //[Ndata.sbs.tr.vx]
   Int_t           Ndata_sbs_tr_vy;
   Double_t        sbs_tr_vy[7];   //[Ndata.sbs.tr.vy]
   Int_t           Ndata_sbs_tr_vz;
   Double_t        sbs_tr_vz[7];   //[Ndata.sbs.tr.vz]
   Int_t           Ndata_sbs_tr_x;
   Double_t        sbs_tr_x[7];   //[Ndata.sbs.tr.x]
   Int_t           Ndata_sbs_tr_y;
   Double_t        sbs_tr_y[7];   //[Ndata.sbs.tr.y]
   Int_t           Ndata_sbs_trig_a_amp_p;
   Double_t        sbs_trig_a_amp_p[11];   //[Ndata.sbs.trig.a_amp_p]
   Int_t           Ndata_sbs_trig_a_p;
   Double_t        sbs_trig_a_p[11];   //[Ndata.sbs.trig.a_p]
   Int_t           Ndata_sbs_trig_a_time;
   Double_t        sbs_trig_a_time[11];   //[Ndata.sbs.trig.a_time]
   Int_t           Ndata_sbs_trig_adcelemID;
   Double_t        sbs_trig_adcelemID[11];   //[Ndata.sbs.trig.adcelemID]
   Int_t           Ndata_sbs_x_bcp;
   Double_t        sbs_x_bcp[35];   //[Ndata.sbs.x_bcp]
   Int_t           Ndata_sbs_x_fcp;
   Double_t        sbs_x_fcp[35];   //[Ndata.sbs.x_fcp]
   Int_t           Ndata_sbs_y_bcp;
   Double_t        sbs_y_bcp[35];   //[Ndata.sbs.y_bcp]
   Int_t           Ndata_sbs_y_fcp;
   Double_t        sbs_y_fcp[35];   //[Ndata.sbs.y_fcp]
   Int_t           Ndata_sbs_z_bcp;
   Double_t        sbs_z_bcp[35];   //[Ndata.sbs.z_bcp]
   Int_t           Ndata_sbs_z_fcp;
   Double_t        sbs_z_fcp[35];   //[Ndata.sbs.z_fcp]
   Double_t        FTROI_ECAL_energy;
   Double_t        FTROI_Eprime;
   Double_t        FTROI_ephi;
   Double_t        FTROI_etheta;
   Double_t        FTROI_good;
   Double_t        FTROI_pp;
   Double_t        FTROI_pphi;
   Double_t        FTROI_ptheta;
   Double_t        FTROI_xECAL_global;
   Double_t        FTROI_yECAL_global;
   Double_t        FTROI_zECAL_global;
   Double_t        Lrb_BPMA_rawcur_1;
   Double_t        Lrb_BPMA_rawcur_2;
   Double_t        Lrb_BPMA_rawcur_3;
   Double_t        Lrb_BPMA_rawcur_4;
   Double_t        Lrb_BPMA_rotpos1;
   Double_t        Lrb_BPMA_rotpos2;
   Double_t        Lrb_BPMA_x;
   Double_t        Lrb_BPMA_y;
   Double_t        Lrb_BPMA_z;
   Double_t        Lrb_BPMB_rawcur_1;
   Double_t        Lrb_BPMB_rawcur_2;
   Double_t        Lrb_BPMB_rawcur_3;
   Double_t        Lrb_BPMB_rawcur_4;
   Double_t        Lrb_BPMB_rotpos1;
   Double_t        Lrb_BPMB_rotpos2;
   Double_t        Lrb_BPMB_x;
   Double_t        Lrb_BPMB_y;
   Double_t        Lrb_BPMB_z;
   Double_t        Lrb_Raster_bpma_x;
   Double_t        Lrb_Raster_bpma_y;
   Double_t        Lrb_Raster_bpma_z;
   Double_t        Lrb_Raster_bpmb_x;
   Double_t        Lrb_Raster_bpmb_y;
   Double_t        Lrb_Raster_bpmb_z;
   Double_t        Lrb_Raster_rawcur_x;
   Double_t        Lrb_Raster_rawcur_y;
   Double_t        Lrb_Raster_rawslope_x;
   Double_t        Lrb_Raster_rawslope_y;
   Double_t        Lrb_Raster_target_dir_x;
   Double_t        Lrb_Raster_target_dir_y;
   Double_t        Lrb_Raster_target_dir_z;
   Double_t        Lrb_Raster_target_x;
   Double_t        Lrb_Raster_target_y;
   Double_t        Lrb_Raster_target_z;
   Double_t        Lrb_Raster2_bpma_x;
   Double_t        Lrb_Raster2_bpma_y;
   Double_t        Lrb_Raster2_bpma_z;
   Double_t        Lrb_Raster2_bpmb_x;
   Double_t        Lrb_Raster2_bpmb_y;
   Double_t        Lrb_Raster2_bpmb_z;
   Double_t        Lrb_Raster2_rawcur_x;
   Double_t        Lrb_Raster2_rawcur_y;
   Double_t        Lrb_Raster2_rawslope_x;
   Double_t        Lrb_Raster2_rawslope_y;
   Double_t        Lrb_Raster2_target_dir_x;
   Double_t        Lrb_Raster2_target_dir_y;
   Double_t        Lrb_Raster2_target_dir_z;
   Double_t        Lrb_Raster2_target_x;
   Double_t        Lrb_Raster2_target_y;
   Double_t        Lrb_Raster2_target_z;
   Double_t        SBS_gold_beta;
   Double_t        SBS_gold_dp;
   Double_t        SBS_gold_index;
   Double_t        SBS_gold_ok;
   Double_t        SBS_gold_p;
   Double_t        SBS_gold_ph;
   Double_t        SBS_gold_px;
   Double_t        SBS_gold_py;
   Double_t        SBS_gold_pz;
   Double_t        SBS_gold_th;
   Double_t        SBS_gold_x;
   Double_t        SBS_gold_y;
   Double_t        SBSrb_BPMA_rawcur_1;
   Double_t        SBSrb_BPMA_rawcur_2;
   Double_t        SBSrb_BPMA_rawcur_3;
   Double_t        SBSrb_BPMA_rawcur_4;
   Double_t        SBSrb_BPMA_rotpos1;
   Double_t        SBSrb_BPMA_rotpos2;
   Double_t        SBSrb_BPMA_x;
   Double_t        SBSrb_BPMA_y;
   Double_t        SBSrb_BPMA_z;
   Double_t        SBSrb_BPMB_rawcur_1;
   Double_t        SBSrb_BPMB_rawcur_2;
   Double_t        SBSrb_BPMB_rawcur_3;
   Double_t        SBSrb_BPMB_rawcur_4;
   Double_t        SBSrb_BPMB_rotpos1;
   Double_t        SBSrb_BPMB_rotpos2;
   Double_t        SBSrb_BPMB_x;
   Double_t        SBSrb_BPMB_y;
   Double_t        SBSrb_BPMB_z;
   Double_t        SBSrb_Raster_bpma_x;
   Double_t        SBSrb_Raster_bpma_y;
   Double_t        SBSrb_Raster_bpma_z;
   Double_t        SBSrb_Raster_bpmb_x;
   Double_t        SBSrb_Raster_bpmb_y;
   Double_t        SBSrb_Raster_bpmb_z;
   Double_t        SBSrb_Raster_rawcur_x;
   Double_t        SBSrb_Raster_rawcur_y;
   Double_t        SBSrb_Raster_rawslope_x;
   Double_t        SBSrb_Raster_rawslope_y;
   Double_t        SBSrb_Raster_target_dir_x;
   Double_t        SBSrb_Raster_target_dir_y;
   Double_t        SBSrb_Raster_target_dir_z;
   Double_t        SBSrb_Raster_target_x;
   Double_t        SBSrb_Raster_target_y;
   Double_t        SBSrb_Raster_target_z;
   Double_t        SBSrb_Raster2_bpma_x;
   Double_t        SBSrb_Raster2_bpma_y;
   Double_t        SBSrb_Raster2_bpma_z;
   Double_t        SBSrb_Raster2_bpmb_x;
   Double_t        SBSrb_Raster2_bpmb_y;
   Double_t        SBSrb_Raster2_bpmb_z;
   Double_t        SBSrb_Raster2_rawcur_x;
   Double_t        SBSrb_Raster2_rawcur_y;
   Double_t        SBSrb_Raster2_rawslope_x;
   Double_t        SBSrb_Raster2_rawslope_y;
   Double_t        SBSrb_Raster2_target_dir_x;
   Double_t        SBSrb_Raster2_target_dir_y;
   Double_t        SBSrb_Raster2_target_dir_z;
   Double_t        SBSrb_Raster2_target_x;
   Double_t        SBSrb_Raster2_target_y;
   Double_t        SBSrb_Raster2_target_z;
   Double_t        earm_cdet_ngoodADChits;
   Double_t        earm_cdet_ngoodTDChits;
   Double_t        earm_cdet_ngoodhits;
   Double_t        earm_cdet_nhits;
   Double_t        earm_cdet_nrefhits;
   Double_t        earm_ecal_adctime;
   Double_t        earm_ecal_againblk;
   Double_t        earm_ecal_atimeblk;
   Double_t        earm_ecal_colblk;
   Double_t        earm_ecal_e;
   Double_t        earm_ecal_eblk;
   Double_t        earm_ecal_idblk;
   Double_t        earm_ecal_index;
   Double_t        earm_ecal_nblk;
   Double_t        earm_ecal_nclus;
   Double_t        earm_ecal_ngoodADChits;
   Double_t        earm_ecal_rowblk;
   Double_t        earm_ecal_vtp_detid;
   Double_t        earm_ecal_vtp_errflag;
   Double_t        earm_ecal_x;
   Double_t        earm_ecal_y;
   Double_t        earm_tr_n;
   Double_t        g_datatype;
   Double_t        g_evlen;
   Double_t        g_evnum;
   Double_t        g_evtime;
   Double_t        g_evtyp;
   Double_t        g_runnum;
   Double_t        g_runtime;
   Double_t        g_runtype;
   Double_t        g_trigbits;
   Double_t        heep_K_eth;
   Double_t        heep_K_p4vect;
   Double_t        heep_K_pp;
   Double_t        heep_K_pth;
   Double_t        heep_Q2_e4vect;
   Double_t        heep_Q2_eth;
   Double_t        heep_Q2_p4vect;
   Double_t        heep_Q2_pp;
   Double_t        heep_Q2_pth;
   Double_t        heep_acoplanarity;
   Double_t        heep_datavalid;
   Double_t        heep_dpe;
   Double_t        heep_dphi;
   Double_t        heep_dpp;
   Double_t        heep_dt_ADC;
   Double_t        heep_dt_TDC;
   Double_t        heep_dxECAL;
   Double_t        heep_dxECAL4vect;
   Double_t        heep_dyECAL;
   Double_t        heep_dyECAL4vect;
   Double_t        heep_ecalo;
   Double_t        heep_ephi;
   Double_t        heep_eprime_eth;
   Double_t        heep_eprime_pp;
   Double_t        heep_eprime_pth;
   Double_t        heep_eps_4vect;
   Double_t        heep_eps_eth;
   Double_t        heep_eps_pp;
   Double_t        heep_eps_pth;
   Double_t        heep_eth_pp;
   Double_t        heep_eth_pth;
   Double_t        heep_etheta;
   Double_t        heep_inelasticity;
   Double_t        heep_pp;
   Double_t        heep_pp_eth;
   Double_t        heep_pp_pth;
   Double_t        heep_pphi;
   Double_t        heep_pth_eth;
   Double_t        heep_ptheta;
   Double_t        sbs_HCALdir_x;
   Double_t        sbs_HCALdir_y;
   Double_t        sbs_HCALdir_z;
   Double_t        sbs_HCALph_n;
   Double_t        sbs_HCALth_n;
   Double_t        sbs_gemFPP_hit_ngoodhits;
   Double_t        sbs_gemFPP_nlayershit;
   Double_t        sbs_gemFPP_nlayershitu;
   Double_t        sbs_gemFPP_nlayershituv;
   Double_t        sbs_gemFPP_nlayershitv;
   Double_t        sbs_gemFPP_track_besttrack;
   Double_t        sbs_gemFPP_track_ntrack;
   Double_t        sbs_gemFT_hit_ngoodhits;
   Double_t        sbs_gemFT_nlayershit;
   Double_t        sbs_gemFT_nlayershitu;
   Double_t        sbs_gemFT_nlayershituv;
   Double_t        sbs_gemFT_nlayershitv;
   Double_t        sbs_gemFT_track_besttrack;
   Double_t        sbs_gemFT_track_ntrack;
   Double_t        sbs_hcal_adctime;
   Double_t        sbs_hcal_againblk;
   Double_t        sbs_hcal_atimeblk;
   Double_t        sbs_hcal_colblk;
   Double_t        sbs_hcal_colblk_goodtdc;
   Double_t        sbs_hcal_e;
   Double_t        sbs_hcal_e_goodtdc;
   Double_t        sbs_hcal_eblk;
   Double_t        sbs_hcal_eblk_goodtdc;
   Double_t        sbs_hcal_idblk;
   Double_t        sbs_hcal_idblk_goodtdc;
   Double_t        sbs_hcal_index;
   Double_t        sbs_hcal_nblk;
   Double_t        sbs_hcal_nblk_goodtdc;
   Double_t        sbs_hcal_nclus;
   Double_t        sbs_hcal_ngoodADChits;
   Double_t        sbs_hcal_rowblk;
   Double_t        sbs_hcal_rowblk_goodtdc;
   Double_t        sbs_hcal_tdctime;
   Double_t        sbs_hcal_tdctime_tw;
   Double_t        sbs_hcal_tdctimeblk;
   Double_t        sbs_hcal_tdctimeblk_tw;
   Double_t        sbs_hcal_trigtime;
   Double_t        sbs_hcal_vtp_detid;
   Double_t        sbs_hcal_vtp_errflag;
   Double_t        sbs_hcal_x;
   Double_t        sbs_hcal_y;
   Double_t        sbs_tr_n;
   Double_t        scalhel_errcode;
   Double_t        scalhel_evtcount;
   Double_t        scalhel_hel;
   Double_t        scalhel_lhrs_fadc_hel;
   Double_t        scalhel_lhrs_fadc_pat;
   Double_t        scalhel_lhrs_fadc_tsettle;
   Double_t        scalhel_patcount;
   Double_t        scalhel_patphase;
   Double_t        scalhel_seed;
   Double_t        singletrack;
   Double_t        anytrack;
   Double_t        singletrack;
   Double_t        anytrack;
   Double_t        goodheep;
   Double_t        gooddxdy;
   Double_t        gooddpp;
   Double_t        goodctime;
   Double_t        HALLA_p;
   Double_t        hac_bcm_average;
   Double_t        IPM1H04A_XPOS;
   Double_t        IPM1H04A_YPOS;
   Double_t        IPM1H04E_XPOS;
   Double_t        IPM1H04E_YPOS;
   Double_t        IGL1I00OD16_16;
 //THaEvent        *Event_Branch;
   ULong64_t       fEvtHdr_fEvtTime;
   UInt_t          fEvtHdr_fEvtNum;
   UInt_t          fEvtHdr_fEvtType;
   UInt_t          fEvtHdr_fEvtLen;
   Int_t           fEvtHdr_fHelicity;
   UInt_t          fEvtHdr_fTrigBits;
   UInt_t          fEvtHdr_fRun;

   // List of branches
   TBranch        *b_ax_new;   //!
   TBranch        *b_bx_new;   //!
   TBranch        *b_ay_new;   //!
   TBranch        *b_by_new;   //!
   TBranch        *b_targx;   //!
   TBranch        *b_targy;   //!
   TBranch        *b_scalhel_true_hel;   //!
   TBranch        *b_scalhel_hel_stable;   //!
   TBranch        *b_dEcalo;   //!
   TBranch        *b_detheta;   //!
   TBranch        *b_dptheta;   //!
   TBranch        *b_Ndata_earm_cdet_hit_col;   //!
   TBranch        *b_earm_cdet_hit_col;   //!
   TBranch        *b_Ndata_earm_cdet_hit_layer;   //!
   TBranch        *b_earm_cdet_hit_layer;   //!
   TBranch        *b_Ndata_earm_cdet_hit_pmtnum;   //!
   TBranch        *b_earm_cdet_hit_pmtnum;   //!
   TBranch        *b_Ndata_earm_cdet_hit_row;   //!
   TBranch        *b_earm_cdet_hit_row;   //!
   TBranch        *b_Ndata_earm_cdet_hit_tdc_le;   //!
   TBranch        *b_earm_cdet_hit_tdc_le;   //!
   TBranch        *b_Ndata_earm_cdet_hit_tdc_te;   //!
   TBranch        *b_earm_cdet_hit_tdc_te;   //!
   TBranch        *b_Ndata_earm_cdet_hit_tdc_tot;   //!
   TBranch        *b_earm_cdet_hit_tdc_tot;   //!
   TBranch        *b_Ndata_earm_cdet_hit_xhit;   //!
   TBranch        *b_earm_cdet_hit_xhit;   //!
   TBranch        *b_Ndata_earm_cdet_hit_yhit;   //!
   TBranch        *b_earm_cdet_hit_yhit;   //!
   TBranch        *b_Ndata_earm_cdet_hit_zhit;   //!
   TBranch        *b_earm_cdet_hit_zhit;   //!
   TBranch        *b_Ndata_earm_cdet_hits_TDCelemID;   //!
   TBranch        *b_earm_cdet_hits_TDCelemID;   //!
   TBranch        *b_Ndata_earm_cdet_hits_t;   //!
   TBranch        *b_earm_cdet_hits_t;   //!
   TBranch        *b_Ndata_earm_cdet_hits_t_te;   //!
   TBranch        *b_earm_cdet_hits_t_te;   //!
   TBranch        *b_Ndata_earm_cdet_hits_t_tot;   //!
   TBranch        *b_earm_cdet_hits_t_tot;   //!
   TBranch        *b_Ndata_earm_cdet_tdc;   //!
   TBranch        *b_earm_cdet_tdc;   //!
   TBranch        *b_Ndata_earm_cdet_tdc_mult;   //!
   TBranch        *b_earm_cdet_tdc_mult;   //!
   TBranch        *b_Ndata_earm_cdet_tdc_te;   //!
   TBranch        *b_earm_cdet_tdc_te;   //!
   TBranch        *b_Ndata_earm_cdet_tdc_tot;   //!
   TBranch        *b_earm_cdet_tdc_tot;   //!
   TBranch        *b_Ndata_earm_cdet_tdccol;   //!
   TBranch        *b_earm_cdet_tdccol;   //!
   TBranch        *b_Ndata_earm_cdet_tdcelemID;   //!
   TBranch        *b_earm_cdet_tdcelemID;   //!
   TBranch        *b_Ndata_earm_cdet_tdclayer;   //!
   TBranch        *b_earm_cdet_tdclayer;   //!
   TBranch        *b_Ndata_earm_cdet_tdcrow;   //!
   TBranch        *b_earm_cdet_tdcrow;   //!
   TBranch        *b_Ndata_earm_ecal_a_c;   //!
   TBranch        *b_earm_ecal_a_c;   //!
   TBranch        *b_Ndata_earm_ecal_a_p;   //!
   TBranch        *b_earm_ecal_a_p;   //!
   TBranch        *b_Ndata_earm_ecal_a_time;   //!
   TBranch        *b_earm_ecal_a_time;   //!
   TBranch        *b_Ndata_earm_ecal_adccol;   //!
   TBranch        *b_earm_ecal_adccol;   //!
   TBranch        *b_Ndata_earm_ecal_adcelemID;   //!
   TBranch        *b_earm_ecal_adcelemID;   //!
   TBranch        *b_Ndata_earm_ecal_adcrow;   //!
   TBranch        *b_earm_ecal_adcrow;   //!
   TBranch        *b_Ndata_earm_ecal_adcxpos;   //!
   TBranch        *b_earm_ecal_adcxpos;   //!
   TBranch        *b_Ndata_earm_ecal_adcypos;   //!
   TBranch        *b_earm_ecal_adcypos;   //!
   TBranch        *b_Ndata_earm_ecal_clus_atimeblk;   //!
   TBranch        *b_earm_ecal_clus_atimeblk;   //!
   TBranch        *b_Ndata_earm_ecal_clus_col;   //!
   TBranch        *b_earm_ecal_clus_col;   //!
   TBranch        *b_Ndata_earm_ecal_clus_e;   //!
   TBranch        *b_earm_ecal_clus_e;   //!
   TBranch        *b_Ndata_earm_ecal_clus_eblk;   //!
   TBranch        *b_earm_ecal_clus_eblk;   //!
   TBranch        *b_Ndata_earm_ecal_clus_id;   //!
   TBranch        *b_earm_ecal_clus_id;   //!
   TBranch        *b_Ndata_earm_ecal_clus_nblk;   //!
   TBranch        *b_earm_ecal_clus_nblk;   //!
   TBranch        *b_Ndata_earm_ecal_clus_row;   //!
   TBranch        *b_earm_ecal_clus_row;   //!
   TBranch        *b_Ndata_earm_ecal_clus_x;   //!
   TBranch        *b_earm_ecal_clus_x;   //!
   TBranch        *b_Ndata_earm_ecal_clus_y;   //!
   TBranch        *b_earm_ecal_clus_y;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_atime;   //!
   TBranch        *b_earm_ecal_clus_blk_atime;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_col;   //!
   TBranch        *b_earm_ecal_clus_blk_col;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_e;   //!
   TBranch        *b_earm_ecal_clus_blk_e;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_id;   //!
   TBranch        *b_earm_ecal_clus_blk_id;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_row;   //!
   TBranch        *b_earm_ecal_clus_blk_row;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_x;   //!
   TBranch        *b_earm_ecal_clus_blk_x;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_y;   //!
   TBranch        *b_earm_ecal_clus_blk_y;   //!
   TBranch        *b_Ndata_earm_ecal_ped;   //!
   TBranch        *b_earm_ecal_ped;   //!
   TBranch        *b_Ndata_earm_ecal_vtp_clus_e;   //!
   TBranch        *b_earm_ecal_vtp_clus_e;   //!
   TBranch        *b_Ndata_earm_ecal_vtp_clus_size;   //!
   TBranch        *b_earm_ecal_vtp_clus_size;   //!
   TBranch        *b_Ndata_earm_ecal_vtp_clus_time;   //!
   TBranch        *b_earm_ecal_vtp_clus_time;   //!
   TBranch        *b_Ndata_earm_ecal_vtp_clus_x;   //!
   TBranch        *b_earm_ecal_vtp_clus_x;   //!
   TBranch        *b_Ndata_earm_ecal_vtp_clus_y;   //!
   TBranch        *b_earm_ecal_vtp_clus_y;   //!
   TBranch        *b_Ndata_earm_tr_beta;   //!
   TBranch        *b_earm_tr_beta;   //!
   TBranch        *b_Ndata_earm_tr_chi2;   //!
   TBranch        *b_earm_tr_chi2;   //!
   TBranch        *b_Ndata_earm_tr_d_ph;   //!
   TBranch        *b_earm_tr_d_ph;   //!
   TBranch        *b_Ndata_earm_tr_d_th;   //!
   TBranch        *b_earm_tr_d_th;   //!
   TBranch        *b_Ndata_earm_tr_d_x;   //!
   TBranch        *b_earm_tr_d_x;   //!
   TBranch        *b_Ndata_earm_tr_d_y;   //!
   TBranch        *b_earm_tr_d_y;   //!
   TBranch        *b_Ndata_earm_tr_dbeta;   //!
   TBranch        *b_earm_tr_dbeta;   //!
   TBranch        *b_Ndata_earm_tr_dtime;   //!
   TBranch        *b_earm_tr_dtime;   //!
   TBranch        *b_Ndata_earm_tr_flag;   //!
   TBranch        *b_earm_tr_flag;   //!
   TBranch        *b_Ndata_earm_tr_ndof;   //!
   TBranch        *b_earm_tr_ndof;   //!
   TBranch        *b_Ndata_earm_tr_p;   //!
   TBranch        *b_earm_tr_p;   //!
   TBranch        *b_Ndata_earm_tr_pathl;   //!
   TBranch        *b_earm_tr_pathl;   //!
   TBranch        *b_Ndata_earm_tr_ph;   //!
   TBranch        *b_earm_tr_ph;   //!
   TBranch        *b_Ndata_earm_tr_px;   //!
   TBranch        *b_earm_tr_px;   //!
   TBranch        *b_Ndata_earm_tr_py;   //!
   TBranch        *b_earm_tr_py;   //!
   TBranch        *b_Ndata_earm_tr_pz;   //!
   TBranch        *b_earm_tr_pz;   //!
   TBranch        *b_Ndata_earm_tr_r_ph;   //!
   TBranch        *b_earm_tr_r_ph;   //!
   TBranch        *b_Ndata_earm_tr_r_th;   //!
   TBranch        *b_earm_tr_r_th;   //!
   TBranch        *b_Ndata_earm_tr_r_x;   //!
   TBranch        *b_earm_tr_r_x;   //!
   TBranch        *b_Ndata_earm_tr_r_y;   //!
   TBranch        *b_earm_tr_r_y;   //!
   TBranch        *b_Ndata_earm_tr_tg_dp;   //!
   TBranch        *b_earm_tr_tg_dp;   //!
   TBranch        *b_Ndata_earm_tr_tg_ph;   //!
   TBranch        *b_earm_tr_tg_ph;   //!
   TBranch        *b_Ndata_earm_tr_tg_th;   //!
   TBranch        *b_earm_tr_tg_th;   //!
   TBranch        *b_Ndata_earm_tr_tg_y;   //!
   TBranch        *b_earm_tr_tg_y;   //!
   TBranch        *b_Ndata_earm_tr_th;   //!
   TBranch        *b_earm_tr_th;   //!
   TBranch        *b_Ndata_earm_tr_time;   //!
   TBranch        *b_earm_tr_time;   //!
   TBranch        *b_Ndata_earm_tr_vx;   //!
   TBranch        *b_earm_tr_vx;   //!
   TBranch        *b_Ndata_earm_tr_vy;   //!
   TBranch        *b_earm_tr_vy;   //!
   TBranch        *b_Ndata_earm_tr_vz;   //!
   TBranch        *b_earm_tr_vz;   //!
   TBranch        *b_Ndata_earm_tr_x;   //!
   TBranch        *b_earm_tr_x;   //!
   TBranch        *b_Ndata_earm_tr_y;   //!
   TBranch        *b_earm_tr_y;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCU;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCU_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCU_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCV;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCV_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCV_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCasym;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCasym;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCasym_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCasym_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCavg;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCavg;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCavg_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCavg_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac0_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac0_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac0_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac0_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac1_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac1_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac1_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac1_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac2_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac2_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac3_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac3_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac3_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac3_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac4_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac4_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac4_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac4_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac5_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac5_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac5_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac5_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampU;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampUclust;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampUclust;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampUclust_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampUclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampV;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampVclust;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampVclust;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampVclust_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampVclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxstripU;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxstripU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxstripV;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxstripV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_CM_GOOD_U;   //!
   TBranch        *b_sbs_gemFPP_hit_CM_GOOD_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_CM_GOOD_V;   //!
   TBranch        *b_sbs_gemFPP_hit_CM_GOOD_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC0_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC0_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC0_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC0_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC1_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC1_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC1_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC1_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC2_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC2_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC3_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC3_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC3_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC3_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC4_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC4_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC4_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC4_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC5_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC5_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC5_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC5_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboU;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxcomboU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboV;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxcomboV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxsampU;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxsampU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxsampV;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxsampV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxstripU;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxstripU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxstripV;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxstripV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ENABLE_CM_U;   //!
   TBranch        *b_sbs_gemFPP_hit_ENABLE_CM_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ENABLE_CM_V;   //!
   TBranch        *b_sbs_gemFPP_hit_ENABLE_CM_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_TSchi2_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_TSchi2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_TSchi2_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_TSchi2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_TSprob_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_TSprob_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_TSprob_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_TSprob_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Tavg;   //!
   TBranch        *b_sbs_gemFPP_hit_Tavg;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Tavg_corr;   //!
   TBranch        *b_sbs_gemFPP_hit_Tavg_corr;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Tavg_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_Tavg_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Tavg_fit;   //!
   TBranch        *b_sbs_gemFPP_hit_Tavg_fit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Ugain;   //!
   TBranch        *b_sbs_gemFPP_hit_Ugain;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Utime;   //!
   TBranch        *b_sbs_gemFPP_hit_Utime;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_UtimeDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_UtimeDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_UtimeFit;   //!
   TBranch        *b_sbs_gemFPP_hit_UtimeFit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_UtimeMaxStrip;   //!
   TBranch        *b_sbs_gemFPP_hit_UtimeMaxStrip;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_UtimeMaxStripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_UtimeMaxStripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_UtimeMaxStripFit;   //!
   TBranch        *b_sbs_gemFPP_hit_UtimeMaxStripFit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Vgain;   //!
   TBranch        *b_sbs_gemFPP_hit_Vgain;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Vtime;   //!
   TBranch        *b_sbs_gemFPP_hit_Vtime;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_VtimeDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_VtimeDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_VtimeFit;   //!
   TBranch        *b_sbs_gemFPP_hit_VtimeFit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_VtimeMaxStrip;   //!
   TBranch        *b_sbs_gemFPP_hit_VtimeMaxStrip;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_VtimeMaxStripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_VtimeMaxStripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_VtimeMaxStripFit;   //!
   TBranch        *b_sbs_gemFPP_hit_VtimeMaxStripFit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_adc_id_U;   //!
   TBranch        *b_sbs_gemFPP_hit_adc_id_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_adc_id_V;   //!
   TBranch        *b_sbs_gemFPP_hit_adc_id_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ccor_clust;   //!
   TBranch        *b_sbs_gemFPP_hit_ccor_clust;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ccor_clust_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ccor_clust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ccor_strip;   //!
   TBranch        *b_sbs_gemFPP_hit_ccor_strip;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ccor_strip_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ccor_strip_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_crate_U;   //!
   TBranch        *b_sbs_gemFPP_hit_crate_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_crate_V;   //!
   TBranch        *b_sbs_gemFPP_hit_crate_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_deltat;   //!
   TBranch        *b_sbs_gemFPP_hit_deltat;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_deltat_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_deltat_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_deltat_fit;   //!
   TBranch        *b_sbs_gemFPP_hit_deltat_fit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_eresidu;   //!
   TBranch        *b_sbs_gemFPP_hit_eresidu;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_eresidv;   //!
   TBranch        *b_sbs_gemFPP_hit_eresidv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_icombomaxUclustDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_icombomaxUclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_icombomaxUstripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_icombomaxUstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_icombomaxVclustDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_icombomaxVclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_icombomaxVstripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_icombomaxVstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxUclust;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxUclust;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxUclustDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxUclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxUstrip;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxUstrip;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxUstripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxUstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxVclust;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxVclust;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxVclustDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxVclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxVstrip;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxVstrip;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxVstripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxVstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_layer;   //!
   TBranch        *b_sbs_gemFPP_hit_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_module;   //!
   TBranch        *b_sbs_gemFPP_hit_module;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_mpd_U;   //!
   TBranch        *b_sbs_gemFPP_hit_mpd_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_mpd_V;   //!
   TBranch        *b_sbs_gemFPP_hit_mpd_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_nstripu;   //!
   TBranch        *b_sbs_gemFPP_hit_nstripu;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_nstripv;   //!
   TBranch        *b_sbs_gemFPP_hit_nstripv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_residu;   //!
   TBranch        *b_sbs_gemFPP_hit_residu;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_residv;   //!
   TBranch        *b_sbs_gemFPP_hit_residv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_trackindex;   //!
   TBranch        *b_sbs_gemFPP_hit_trackindex;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_u;   //!
   TBranch        *b_sbs_gemFPP_hit_u;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_umoment;   //!
   TBranch        *b_sbs_gemFPP_hit_umoment;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_usigma;   //!
   TBranch        *b_sbs_gemFPP_hit_usigma;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ustriphi;   //!
   TBranch        *b_sbs_gemFPP_hit_ustriphi;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ustriplo;   //!
   TBranch        *b_sbs_gemFPP_hit_ustriplo;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ustripmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ustripmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_v;   //!
   TBranch        *b_sbs_gemFPP_hit_v;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_vmoment;   //!
   TBranch        *b_sbs_gemFPP_hit_vmoment;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_vsigma;   //!
   TBranch        *b_sbs_gemFPP_hit_vsigma;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_vstriphi;   //!
   TBranch        *b_sbs_gemFPP_hit_vstriphi;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_vstriplo;   //!
   TBranch        *b_sbs_gemFPP_hit_vstriplo;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_vstripmax;   //!
   TBranch        *b_sbs_gemFPP_hit_vstripmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_xglobal;   //!
   TBranch        *b_sbs_gemFPP_hit_xglobal;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_xlocal;   //!
   TBranch        *b_sbs_gemFPP_hit_xlocal;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_yglobal;   //!
   TBranch        *b_sbs_gemFPP_hit_yglobal;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ylocal;   //!
   TBranch        *b_sbs_gemFPP_hit_ylocal;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_zglobal;   //!
   TBranch        *b_sbs_gemFPP_hit_zglobal;   //!
   TBranch        *b_Ndata_sbs_gemFPP_n2Dhit_layer;   //!
   TBranch        *b_sbs_gemFPP_n2Dhit_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_nclustu_layer;   //!
   TBranch        *b_sbs_gemFPP_nclustu_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_nclustv_layer;   //!
   TBranch        *b_sbs_gemFPP_nclustv_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_nstripsu_layer;   //!
   TBranch        *b_sbs_gemFPP_nstripsu_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_nstripsv_layer;   //!
   TBranch        *b_sbs_gemFPP_nstripsv_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_chi2ndf;   //!
   TBranch        *b_sbs_gemFPP_track_chi2ndf;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_chi2ndf_hitquality;   //!
   TBranch        *b_sbs_gemFPP_track_chi2ndf_hitquality;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_ngoodhits;   //!
   TBranch        *b_sbs_gemFPP_track_ngoodhits;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_nhits;   //!
   TBranch        *b_sbs_gemFPP_track_nhits;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_phi;   //!
   TBranch        *b_sbs_gemFPP_track_phi;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_sclose;   //!
   TBranch        *b_sbs_gemFPP_track_sclose;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_t0;   //!
   TBranch        *b_sbs_gemFPP_track_t0;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_theta;   //!
   TBranch        *b_sbs_gemFPP_track_theta;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_x;   //!
   TBranch        *b_sbs_gemFPP_track_x;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_xp;   //!
   TBranch        *b_sbs_gemFPP_track_xp;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_y;   //!
   TBranch        *b_sbs_gemFPP_track_y;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_yp;   //!
   TBranch        *b_sbs_gemFPP_track_yp;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_zclose;   //!
   TBranch        *b_sbs_gemFPP_track_zclose;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCU;   //!
   TBranch        *b_sbs_gemFT_hit_ADCU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCU_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCU_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCV;   //!
   TBranch        *b_sbs_gemFT_hit_ADCV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCV_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCV_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCasym;   //!
   TBranch        *b_sbs_gemFT_hit_ADCasym;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCasym_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCasym_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCavg;   //!
   TBranch        *b_sbs_gemFT_hit_ADCavg;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCavg_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCavg_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac0_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac0_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac0_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac0_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac1_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac1_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac1_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac1_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac2_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac2_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac3_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac3_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac3_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac3_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac4_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac4_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac4_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac4_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac5_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac5_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac5_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac5_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxcomboUclust_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxcomboUclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxcomboVclust_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxcomboVclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampU;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampUclust;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampUclust;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampUclust_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampUclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampV;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampVclust;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampVclust;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampVclust_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampVclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxstripU;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxstripU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxstripV;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxstripV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_CM_GOOD_U;   //!
   TBranch        *b_sbs_gemFT_hit_CM_GOOD_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_CM_GOOD_V;   //!
   TBranch        *b_sbs_gemFT_hit_CM_GOOD_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC0_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC0_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC0_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC0_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC1_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC1_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC1_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC1_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC2_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC2_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC3_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC3_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC3_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC3_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC4_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC4_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC4_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC4_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC5_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC5_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC5_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC5_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxcomboU;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxcomboU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxcomboV;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxcomboV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxsampU;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxsampU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxsampV;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxsampV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxstripU;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxstripU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxstripV;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxstripV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ENABLE_CM_U;   //!
   TBranch        *b_sbs_gemFT_hit_ENABLE_CM_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ENABLE_CM_V;   //!
   TBranch        *b_sbs_gemFT_hit_ENABLE_CM_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_TSchi2_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_TSchi2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_TSchi2_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_TSchi2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_TSprob_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_TSprob_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_TSprob_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_TSprob_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Tavg;   //!
   TBranch        *b_sbs_gemFT_hit_Tavg;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Tavg_corr;   //!
   TBranch        *b_sbs_gemFT_hit_Tavg_corr;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Tavg_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_Tavg_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Tavg_fit;   //!
   TBranch        *b_sbs_gemFT_hit_Tavg_fit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Ugain;   //!
   TBranch        *b_sbs_gemFT_hit_Ugain;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Utime;   //!
   TBranch        *b_sbs_gemFT_hit_Utime;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_UtimeDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_UtimeDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_UtimeFit;   //!
   TBranch        *b_sbs_gemFT_hit_UtimeFit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_UtimeMaxStrip;   //!
   TBranch        *b_sbs_gemFT_hit_UtimeMaxStrip;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_UtimeMaxStripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_UtimeMaxStripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_UtimeMaxStripFit;   //!
   TBranch        *b_sbs_gemFT_hit_UtimeMaxStripFit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Vgain;   //!
   TBranch        *b_sbs_gemFT_hit_Vgain;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Vtime;   //!
   TBranch        *b_sbs_gemFT_hit_Vtime;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_VtimeDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_VtimeDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_VtimeFit;   //!
   TBranch        *b_sbs_gemFT_hit_VtimeFit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_VtimeMaxStrip;   //!
   TBranch        *b_sbs_gemFT_hit_VtimeMaxStrip;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_VtimeMaxStripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_VtimeMaxStripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_VtimeMaxStripFit;   //!
   TBranch        *b_sbs_gemFT_hit_VtimeMaxStripFit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_adc_id_U;   //!
   TBranch        *b_sbs_gemFT_hit_adc_id_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_adc_id_V;   //!
   TBranch        *b_sbs_gemFT_hit_adc_id_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ccor_clust;   //!
   TBranch        *b_sbs_gemFT_hit_ccor_clust;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ccor_clust_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ccor_clust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ccor_strip;   //!
   TBranch        *b_sbs_gemFT_hit_ccor_strip;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ccor_strip_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ccor_strip_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_crate_U;   //!
   TBranch        *b_sbs_gemFT_hit_crate_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_crate_V;   //!
   TBranch        *b_sbs_gemFT_hit_crate_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_deltat;   //!
   TBranch        *b_sbs_gemFT_hit_deltat;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_deltat_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_deltat_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_deltat_fit;   //!
   TBranch        *b_sbs_gemFT_hit_deltat_fit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_eresidu;   //!
   TBranch        *b_sbs_gemFT_hit_eresidu;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_eresidv;   //!
   TBranch        *b_sbs_gemFT_hit_eresidv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_icombomaxUclustDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_icombomaxUclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_icombomaxUstripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_icombomaxUstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_icombomaxVclustDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_icombomaxVclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_icombomaxVstripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_icombomaxVstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxUclust;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxUclust;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxUclustDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxUclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxUstrip;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxUstrip;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxUstripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxUstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxVclust;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxVclust;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxVclustDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxVclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxVstrip;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxVstrip;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxVstripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxVstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_layer;   //!
   TBranch        *b_sbs_gemFT_hit_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_module;   //!
   TBranch        *b_sbs_gemFT_hit_module;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_mpd_U;   //!
   TBranch        *b_sbs_gemFT_hit_mpd_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_mpd_V;   //!
   TBranch        *b_sbs_gemFT_hit_mpd_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_nstripu;   //!
   TBranch        *b_sbs_gemFT_hit_nstripu;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_nstripv;   //!
   TBranch        *b_sbs_gemFT_hit_nstripv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_residu;   //!
   TBranch        *b_sbs_gemFT_hit_residu;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_residv;   //!
   TBranch        *b_sbs_gemFT_hit_residv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_trackindex;   //!
   TBranch        *b_sbs_gemFT_hit_trackindex;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_u;   //!
   TBranch        *b_sbs_gemFT_hit_u;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_umoment;   //!
   TBranch        *b_sbs_gemFT_hit_umoment;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_usigma;   //!
   TBranch        *b_sbs_gemFT_hit_usigma;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ustriphi;   //!
   TBranch        *b_sbs_gemFT_hit_ustriphi;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ustriplo;   //!
   TBranch        *b_sbs_gemFT_hit_ustriplo;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ustripmax;   //!
   TBranch        *b_sbs_gemFT_hit_ustripmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_v;   //!
   TBranch        *b_sbs_gemFT_hit_v;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_vmoment;   //!
   TBranch        *b_sbs_gemFT_hit_vmoment;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_vsigma;   //!
   TBranch        *b_sbs_gemFT_hit_vsigma;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_vstriphi;   //!
   TBranch        *b_sbs_gemFT_hit_vstriphi;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_vstriplo;   //!
   TBranch        *b_sbs_gemFT_hit_vstriplo;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_vstripmax;   //!
   TBranch        *b_sbs_gemFT_hit_vstripmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_xglobal;   //!
   TBranch        *b_sbs_gemFT_hit_xglobal;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_xlocal;   //!
   TBranch        *b_sbs_gemFT_hit_xlocal;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_yglobal;   //!
   TBranch        *b_sbs_gemFT_hit_yglobal;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ylocal;   //!
   TBranch        *b_sbs_gemFT_hit_ylocal;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_zglobal;   //!
   TBranch        *b_sbs_gemFT_hit_zglobal;   //!
   TBranch        *b_Ndata_sbs_gemFT_n2Dhit_layer;   //!
   TBranch        *b_sbs_gemFT_n2Dhit_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_nclustu_layer;   //!
   TBranch        *b_sbs_gemFT_nclustu_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_nclustv_layer;   //!
   TBranch        *b_sbs_gemFT_nclustv_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_nstripsu_layer;   //!
   TBranch        *b_sbs_gemFT_nstripsu_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_nstripsv_layer;   //!
   TBranch        *b_sbs_gemFT_nstripsv_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_chi2ndf;   //!
   TBranch        *b_sbs_gemFT_track_chi2ndf;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_chi2ndf_hitquality;   //!
   TBranch        *b_sbs_gemFT_track_chi2ndf_hitquality;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_ngoodhits;   //!
   TBranch        *b_sbs_gemFT_track_ngoodhits;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_nhits;   //!
   TBranch        *b_sbs_gemFT_track_nhits;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_t0;   //!
   TBranch        *b_sbs_gemFT_track_t0;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_x;   //!
   TBranch        *b_sbs_gemFT_track_x;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_xp;   //!
   TBranch        *b_sbs_gemFT_track_xp;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_y;   //!
   TBranch        *b_sbs_gemFT_track_y;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_yp;   //!
   TBranch        *b_sbs_gemFT_track_yp;   //!
   TBranch        *b_Ndata_sbs_hcal_Ref_tdc;   //!
   TBranch        *b_sbs_hcal_Ref_tdc;   //!
   TBranch        *b_Ndata_sbs_hcal_Ref_tdc_mult;   //!
   TBranch        *b_sbs_hcal_Ref_tdc_mult;   //!
   TBranch        *b_Ndata_sbs_hcal_Ref_tdcelemID;   //!
   TBranch        *b_sbs_hcal_Ref_tdcelemID;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_adctime;   //!
   TBranch        *b_sbs_hcal_clus_adctime;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_again;   //!
   TBranch        *b_sbs_hcal_clus_again;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_atimeblk;   //!
   TBranch        *b_sbs_hcal_clus_atimeblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_col;   //!
   TBranch        *b_sbs_hcal_clus_col;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_col_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_col_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_e;   //!
   TBranch        *b_sbs_hcal_clus_e;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_e_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_e_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_eblk;   //!
   TBranch        *b_sbs_hcal_clus_eblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_eblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_eblk_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_id;   //!
   TBranch        *b_sbs_hcal_clus_id;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_id_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_id_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_nblk;   //!
   TBranch        *b_sbs_hcal_clus_nblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_nblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_nblk_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_row;   //!
   TBranch        *b_sbs_hcal_clus_row;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_row_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_row_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_tdctime;   //!
   TBranch        *b_sbs_hcal_clus_tdctime;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_tdctime_tw;   //!
   TBranch        *b_sbs_hcal_clus_tdctime_tw;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_tdctimeblk;   //!
   TBranch        *b_sbs_hcal_clus_tdctimeblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_tdctimeblk_tw;   //!
   TBranch        *b_sbs_hcal_clus_tdctimeblk_tw;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_x;   //!
   TBranch        *b_sbs_hcal_clus_x;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_y;   //!
   TBranch        *b_sbs_hcal_clus_y;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_again;   //!
   TBranch        *b_sbs_hcal_clus_blk_again;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_atime;   //!
   TBranch        *b_sbs_hcal_clus_blk_atime;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_col;   //!
   TBranch        *b_sbs_hcal_clus_blk_col;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_e;   //!
   TBranch        *b_sbs_hcal_clus_blk_e;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_id;   //!
   TBranch        *b_sbs_hcal_clus_blk_id;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_row;   //!
   TBranch        *b_sbs_hcal_clus_blk_row;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_tdctime;   //!
   TBranch        *b_sbs_hcal_clus_blk_tdctime;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_tdctime_tw;   //!
   TBranch        *b_sbs_hcal_clus_blk_tdctime_tw;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_x;   //!
   TBranch        *b_sbs_hcal_clus_blk_x;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_y;   //!
   TBranch        *b_sbs_hcal_clus_blk_y;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_atime;   //!
   TBranch        *b_sbs_hcal_goodblock_atime;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_cid;   //!
   TBranch        *b_sbs_hcal_goodblock_cid;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_col;   //!
   TBranch        *b_sbs_hcal_goodblock_col;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_e;   //!
   TBranch        *b_sbs_hcal_goodblock_e;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_id;   //!
   TBranch        *b_sbs_hcal_goodblock_id;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_row;   //!
   TBranch        *b_sbs_hcal_goodblock_row;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_tdctime;   //!
   TBranch        *b_sbs_hcal_goodblock_tdctime;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_tdctime_tw;   //!
   TBranch        *b_sbs_hcal_goodblock_tdctime_tw;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_x;   //!
   TBranch        *b_sbs_hcal_goodblock_x;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_y;   //!
   TBranch        *b_sbs_hcal_goodblock_y;   //!
   TBranch        *b_Ndata_sbs_hcal_vtp_clus_e;   //!
   TBranch        *b_sbs_hcal_vtp_clus_e;   //!
   TBranch        *b_Ndata_sbs_hcal_vtp_clus_size;   //!
   TBranch        *b_sbs_hcal_vtp_clus_size;   //!
   TBranch        *b_Ndata_sbs_hcal_vtp_clus_time;   //!
   TBranch        *b_sbs_hcal_vtp_clus_time;   //!
   TBranch        *b_Ndata_sbs_hcal_vtp_clus_x;   //!
   TBranch        *b_sbs_hcal_vtp_clus_x;   //!
   TBranch        *b_Ndata_sbs_hcal_vtp_clus_y;   //!
   TBranch        *b_sbs_hcal_vtp_clus_y;   //!
   TBranch        *b_Ndata_sbs_tr_beta;   //!
   TBranch        *b_sbs_tr_beta;   //!
   TBranch        *b_Ndata_sbs_tr_chi2;   //!
   TBranch        *b_sbs_tr_chi2;   //!
   TBranch        *b_Ndata_sbs_tr_d_ph;   //!
   TBranch        *b_sbs_tr_d_ph;   //!
   TBranch        *b_Ndata_sbs_tr_d_th;   //!
   TBranch        *b_sbs_tr_d_th;   //!
   TBranch        *b_Ndata_sbs_tr_d_x;   //!
   TBranch        *b_sbs_tr_d_x;   //!
   TBranch        *b_Ndata_sbs_tr_d_y;   //!
   TBranch        *b_sbs_tr_d_y;   //!
   TBranch        *b_Ndata_sbs_tr_dbeta;   //!
   TBranch        *b_sbs_tr_dbeta;   //!
   TBranch        *b_Ndata_sbs_tr_dtime;   //!
   TBranch        *b_sbs_tr_dtime;   //!
   TBranch        *b_Ndata_sbs_tr_flag;   //!
   TBranch        *b_sbs_tr_flag;   //!
   TBranch        *b_Ndata_sbs_tr_ndof;   //!
   TBranch        *b_sbs_tr_ndof;   //!
   TBranch        *b_Ndata_sbs_tr_p;   //!
   TBranch        *b_sbs_tr_p;   //!
   TBranch        *b_Ndata_sbs_tr_pathl;   //!
   TBranch        *b_sbs_tr_pathl;   //!
   TBranch        *b_Ndata_sbs_tr_ph;   //!
   TBranch        *b_sbs_tr_ph;   //!
   TBranch        *b_Ndata_sbs_tr_px;   //!
   TBranch        *b_sbs_tr_px;   //!
   TBranch        *b_Ndata_sbs_tr_py;   //!
   TBranch        *b_sbs_tr_py;   //!
   TBranch        *b_Ndata_sbs_tr_pz;   //!
   TBranch        *b_sbs_tr_pz;   //!
   TBranch        *b_Ndata_sbs_tr_r_ph;   //!
   TBranch        *b_sbs_tr_r_ph;   //!
   TBranch        *b_Ndata_sbs_tr_r_th;   //!
   TBranch        *b_sbs_tr_r_th;   //!
   TBranch        *b_Ndata_sbs_tr_r_x;   //!
   TBranch        *b_sbs_tr_r_x;   //!
   TBranch        *b_Ndata_sbs_tr_r_y;   //!
   TBranch        *b_sbs_tr_r_y;   //!
   TBranch        *b_Ndata_sbs_tr_tg_dp;   //!
   TBranch        *b_sbs_tr_tg_dp;   //!
   TBranch        *b_Ndata_sbs_tr_tg_ph;   //!
   TBranch        *b_sbs_tr_tg_ph;   //!
   TBranch        *b_Ndata_sbs_tr_tg_th;   //!
   TBranch        *b_sbs_tr_tg_th;   //!
   TBranch        *b_Ndata_sbs_tr_tg_y;   //!
   TBranch        *b_sbs_tr_tg_y;   //!
   TBranch        *b_Ndata_sbs_tr_th;   //!
   TBranch        *b_sbs_tr_th;   //!
   TBranch        *b_Ndata_sbs_tr_time;   //!
   TBranch        *b_sbs_tr_time;   //!
   TBranch        *b_Ndata_sbs_tr_vx;   //!
   TBranch        *b_sbs_tr_vx;   //!
   TBranch        *b_Ndata_sbs_tr_vy;   //!
   TBranch        *b_sbs_tr_vy;   //!
   TBranch        *b_Ndata_sbs_tr_vz;   //!
   TBranch        *b_sbs_tr_vz;   //!
   TBranch        *b_Ndata_sbs_tr_x;   //!
   TBranch        *b_sbs_tr_x;   //!
   TBranch        *b_Ndata_sbs_tr_y;   //!
   TBranch        *b_sbs_tr_y;   //!
   TBranch        *b_Ndata_sbs_trig_a_amp_p;   //!
   TBranch        *b_sbs_trig_a_amp_p;   //!
   TBranch        *b_Ndata_sbs_trig_a_p;   //!
   TBranch        *b_sbs_trig_a_p;   //!
   TBranch        *b_Ndata_sbs_trig_a_time;   //!
   TBranch        *b_sbs_trig_a_time;   //!
   TBranch        *b_Ndata_sbs_trig_adcelemID;   //!
   TBranch        *b_sbs_trig_adcelemID;   //!
   TBranch        *b_Ndata_sbs_x_bcp;   //!
   TBranch        *b_sbs_x_bcp;   //!
   TBranch        *b_Ndata_sbs_x_fcp;   //!
   TBranch        *b_sbs_x_fcp;   //!
   TBranch        *b_Ndata_sbs_y_bcp;   //!
   TBranch        *b_sbs_y_bcp;   //!
   TBranch        *b_Ndata_sbs_y_fcp;   //!
   TBranch        *b_sbs_y_fcp;   //!
   TBranch        *b_Ndata_sbs_z_bcp;   //!
   TBranch        *b_sbs_z_bcp;   //!
   TBranch        *b_Ndata_sbs_z_fcp;   //!
   TBranch        *b_sbs_z_fcp;   //!
   TBranch        *b_FTROI_ECAL_energy;   //!
   TBranch        *b_FTROI_Eprime;   //!
   TBranch        *b_FTROI_ephi;   //!
   TBranch        *b_FTROI_etheta;   //!
   TBranch        *b_FTROI_good;   //!
   TBranch        *b_FTROI_pp;   //!
   TBranch        *b_FTROI_pphi;   //!
   TBranch        *b_FTROI_ptheta;   //!
   TBranch        *b_FTROI_xECAL_global;   //!
   TBranch        *b_FTROI_yECAL_global;   //!
   TBranch        *b_FTROI_zECAL_global;   //!
   TBranch        *b_Lrb_BPMA_rawcur_1;   //!
   TBranch        *b_Lrb_BPMA_rawcur_2;   //!
   TBranch        *b_Lrb_BPMA_rawcur_3;   //!
   TBranch        *b_Lrb_BPMA_rawcur_4;   //!
   TBranch        *b_Lrb_BPMA_rotpos1;   //!
   TBranch        *b_Lrb_BPMA_rotpos2;   //!
   TBranch        *b_Lrb_BPMA_x;   //!
   TBranch        *b_Lrb_BPMA_y;   //!
   TBranch        *b_Lrb_BPMA_z;   //!
   TBranch        *b_Lrb_BPMB_rawcur_1;   //!
   TBranch        *b_Lrb_BPMB_rawcur_2;   //!
   TBranch        *b_Lrb_BPMB_rawcur_3;   //!
   TBranch        *b_Lrb_BPMB_rawcur_4;   //!
   TBranch        *b_Lrb_BPMB_rotpos1;   //!
   TBranch        *b_Lrb_BPMB_rotpos2;   //!
   TBranch        *b_Lrb_BPMB_x;   //!
   TBranch        *b_Lrb_BPMB_y;   //!
   TBranch        *b_Lrb_BPMB_z;   //!
   TBranch        *b_Lrb_Raster_bpma_x;   //!
   TBranch        *b_Lrb_Raster_bpma_y;   //!
   TBranch        *b_Lrb_Raster_bpma_z;   //!
   TBranch        *b_Lrb_Raster_bpmb_x;   //!
   TBranch        *b_Lrb_Raster_bpmb_y;   //!
   TBranch        *b_Lrb_Raster_bpmb_z;   //!
   TBranch        *b_Lrb_Raster_rawcur_x;   //!
   TBranch        *b_Lrb_Raster_rawcur_y;   //!
   TBranch        *b_Lrb_Raster_rawslope_x;   //!
   TBranch        *b_Lrb_Raster_rawslope_y;   //!
   TBranch        *b_Lrb_Raster_target_dir_x;   //!
   TBranch        *b_Lrb_Raster_target_dir_y;   //!
   TBranch        *b_Lrb_Raster_target_dir_z;   //!
   TBranch        *b_Lrb_Raster_target_x;   //!
   TBranch        *b_Lrb_Raster_target_y;   //!
   TBranch        *b_Lrb_Raster_target_z;   //!
   TBranch        *b_Lrb_Raster2_bpma_x;   //!
   TBranch        *b_Lrb_Raster2_bpma_y;   //!
   TBranch        *b_Lrb_Raster2_bpma_z;   //!
   TBranch        *b_Lrb_Raster2_bpmb_x;   //!
   TBranch        *b_Lrb_Raster2_bpmb_y;   //!
   TBranch        *b_Lrb_Raster2_bpmb_z;   //!
   TBranch        *b_Lrb_Raster2_rawcur_x;   //!
   TBranch        *b_Lrb_Raster2_rawcur_y;   //!
   TBranch        *b_Lrb_Raster2_rawslope_x;   //!
   TBranch        *b_Lrb_Raster2_rawslope_y;   //!
   TBranch        *b_Lrb_Raster2_target_dir_x;   //!
   TBranch        *b_Lrb_Raster2_target_dir_y;   //!
   TBranch        *b_Lrb_Raster2_target_dir_z;   //!
   TBranch        *b_Lrb_Raster2_target_x;   //!
   TBranch        *b_Lrb_Raster2_target_y;   //!
   TBranch        *b_Lrb_Raster2_target_z;   //!
   TBranch        *b_SBS_gold_beta;   //!
   TBranch        *b_SBS_gold_dp;   //!
   TBranch        *b_SBS_gold_index;   //!
   TBranch        *b_SBS_gold_ok;   //!
   TBranch        *b_SBS_gold_p;   //!
   TBranch        *b_SBS_gold_ph;   //!
   TBranch        *b_SBS_gold_px;   //!
   TBranch        *b_SBS_gold_py;   //!
   TBranch        *b_SBS_gold_pz;   //!
   TBranch        *b_SBS_gold_th;   //!
   TBranch        *b_SBS_gold_x;   //!
   TBranch        *b_SBS_gold_y;   //!
   TBranch        *b_SBSrb_BPMA_rawcur_1;   //!
   TBranch        *b_SBSrb_BPMA_rawcur_2;   //!
   TBranch        *b_SBSrb_BPMA_rawcur_3;   //!
   TBranch        *b_SBSrb_BPMA_rawcur_4;   //!
   TBranch        *b_SBSrb_BPMA_rotpos1;   //!
   TBranch        *b_SBSrb_BPMA_rotpos2;   //!
   TBranch        *b_SBSrb_BPMA_x;   //!
   TBranch        *b_SBSrb_BPMA_y;   //!
   TBranch        *b_SBSrb_BPMA_z;   //!
   TBranch        *b_SBSrb_BPMB_rawcur_1;   //!
   TBranch        *b_SBSrb_BPMB_rawcur_2;   //!
   TBranch        *b_SBSrb_BPMB_rawcur_3;   //!
   TBranch        *b_SBSrb_BPMB_rawcur_4;   //!
   TBranch        *b_SBSrb_BPMB_rotpos1;   //!
   TBranch        *b_SBSrb_BPMB_rotpos2;   //!
   TBranch        *b_SBSrb_BPMB_x;   //!
   TBranch        *b_SBSrb_BPMB_y;   //!
   TBranch        *b_SBSrb_BPMB_z;   //!
   TBranch        *b_SBSrb_Raster_bpma_x;   //!
   TBranch        *b_SBSrb_Raster_bpma_y;   //!
   TBranch        *b_SBSrb_Raster_bpma_z;   //!
   TBranch        *b_SBSrb_Raster_bpmb_x;   //!
   TBranch        *b_SBSrb_Raster_bpmb_y;   //!
   TBranch        *b_SBSrb_Raster_bpmb_z;   //!
   TBranch        *b_SBSrb_Raster_rawcur_x;   //!
   TBranch        *b_SBSrb_Raster_rawcur_y;   //!
   TBranch        *b_SBSrb_Raster_rawslope_x;   //!
   TBranch        *b_SBSrb_Raster_rawslope_y;   //!
   TBranch        *b_SBSrb_Raster_target_dir_x;   //!
   TBranch        *b_SBSrb_Raster_target_dir_y;   //!
   TBranch        *b_SBSrb_Raster_target_dir_z;   //!
   TBranch        *b_SBSrb_Raster_target_x;   //!
   TBranch        *b_SBSrb_Raster_target_y;   //!
   TBranch        *b_SBSrb_Raster_target_z;   //!
   TBranch        *b_SBSrb_Raster2_bpma_x;   //!
   TBranch        *b_SBSrb_Raster2_bpma_y;   //!
   TBranch        *b_SBSrb_Raster2_bpma_z;   //!
   TBranch        *b_SBSrb_Raster2_bpmb_x;   //!
   TBranch        *b_SBSrb_Raster2_bpmb_y;   //!
   TBranch        *b_SBSrb_Raster2_bpmb_z;   //!
   TBranch        *b_SBSrb_Raster2_rawcur_x;   //!
   TBranch        *b_SBSrb_Raster2_rawcur_y;   //!
   TBranch        *b_SBSrb_Raster2_rawslope_x;   //!
   TBranch        *b_SBSrb_Raster2_rawslope_y;   //!
   TBranch        *b_SBSrb_Raster2_target_dir_x;   //!
   TBranch        *b_SBSrb_Raster2_target_dir_y;   //!
   TBranch        *b_SBSrb_Raster2_target_dir_z;   //!
   TBranch        *b_SBSrb_Raster2_target_x;   //!
   TBranch        *b_SBSrb_Raster2_target_y;   //!
   TBranch        *b_SBSrb_Raster2_target_z;   //!
   TBranch        *b_earm_cdet_ngoodADChits;   //!
   TBranch        *b_earm_cdet_ngoodTDChits;   //!
   TBranch        *b_earm_cdet_ngoodhits;   //!
   TBranch        *b_earm_cdet_nhits;   //!
   TBranch        *b_earm_cdet_nrefhits;   //!
   TBranch        *b_earm_ecal_adctime;   //!
   TBranch        *b_earm_ecal_againblk;   //!
   TBranch        *b_earm_ecal_atimeblk;   //!
   TBranch        *b_earm_ecal_colblk;   //!
   TBranch        *b_earm_ecal_e;   //!
   TBranch        *b_earm_ecal_eblk;   //!
   TBranch        *b_earm_ecal_idblk;   //!
   TBranch        *b_earm_ecal_index;   //!
   TBranch        *b_earm_ecal_nblk;   //!
   TBranch        *b_earm_ecal_nclus;   //!
   TBranch        *b_earm_ecal_ngoodADChits;   //!
   TBranch        *b_earm_ecal_rowblk;   //!
   TBranch        *b_earm_ecal_vtp_detid;   //!
   TBranch        *b_earm_ecal_vtp_errflag;   //!
   TBranch        *b_earm_ecal_x;   //!
   TBranch        *b_earm_ecal_y;   //!
   TBranch        *b_earm_tr_n;   //!
   TBranch        *b_g_datatype;   //!
   TBranch        *b_g_evlen;   //!
   TBranch        *b_g_evnum;   //!
   TBranch        *b_g_evtime;   //!
   TBranch        *b_g_evtyp;   //!
   TBranch        *b_g_runnum;   //!
   TBranch        *b_g_runtime;   //!
   TBranch        *b_g_runtype;   //!
   TBranch        *b_g_trigbits;   //!
   TBranch        *b_heep_K_eth;   //!
   TBranch        *b_heep_K_p4vect;   //!
   TBranch        *b_heep_K_pp;   //!
   TBranch        *b_heep_K_pth;   //!
   TBranch        *b_heep_Q2_e4vect;   //!
   TBranch        *b_heep_Q2_eth;   //!
   TBranch        *b_heep_Q2_p4vect;   //!
   TBranch        *b_heep_Q2_pp;   //!
   TBranch        *b_heep_Q2_pth;   //!
   TBranch        *b_heep_acoplanarity;   //!
   TBranch        *b_heep_datavalid;   //!
   TBranch        *b_heep_dpe;   //!
   TBranch        *b_heep_dphi;   //!
   TBranch        *b_heep_dpp;   //!
   TBranch        *b_heep_dt_ADC;   //!
   TBranch        *b_heep_dt_TDC;   //!
   TBranch        *b_heep_dxECAL;   //!
   TBranch        *b_heep_dxECAL4vect;   //!
   TBranch        *b_heep_dyECAL;   //!
   TBranch        *b_heep_dyECAL4vect;   //!
   TBranch        *b_heep_ecalo;   //!
   TBranch        *b_heep_ephi;   //!
   TBranch        *b_heep_eprime_eth;   //!
   TBranch        *b_heep_eprime_pp;   //!
   TBranch        *b_heep_eprime_pth;   //!
   TBranch        *b_heep_eps_4vect;   //!
   TBranch        *b_heep_eps_eth;   //!
   TBranch        *b_heep_eps_pp;   //!
   TBranch        *b_heep_eps_pth;   //!
   TBranch        *b_heep_eth_pp;   //!
   TBranch        *b_heep_eth_pth;   //!
   TBranch        *b_heep_etheta;   //!
   TBranch        *b_heep_inelasticity;   //!
   TBranch        *b_heep_pp;   //!
   TBranch        *b_heep_pp_eth;   //!
   TBranch        *b_heep_pp_pth;   //!
   TBranch        *b_heep_pphi;   //!
   TBranch        *b_heep_pth_eth;   //!
   TBranch        *b_heep_ptheta;   //!
   TBranch        *b_sbs_HCALdir_x;   //!
   TBranch        *b_sbs_HCALdir_y;   //!
   TBranch        *b_sbs_HCALdir_z;   //!
   TBranch        *b_sbs_HCALph_n;   //!
   TBranch        *b_sbs_HCALth_n;   //!
   TBranch        *b_sbs_gemFPP_hit_ngoodhits;   //!
   TBranch        *b_sbs_gemFPP_nlayershit;   //!
   TBranch        *b_sbs_gemFPP_nlayershitu;   //!
   TBranch        *b_sbs_gemFPP_nlayershituv;   //!
   TBranch        *b_sbs_gemFPP_nlayershitv;   //!
   TBranch        *b_sbs_gemFPP_track_besttrack;   //!
   TBranch        *b_sbs_gemFPP_track_ntrack;   //!
   TBranch        *b_sbs_gemFT_hit_ngoodhits;   //!
   TBranch        *b_sbs_gemFT_nlayershit;   //!
   TBranch        *b_sbs_gemFT_nlayershitu;   //!
   TBranch        *b_sbs_gemFT_nlayershituv;   //!
   TBranch        *b_sbs_gemFT_nlayershitv;   //!
   TBranch        *b_sbs_gemFT_track_besttrack;   //!
   TBranch        *b_sbs_gemFT_track_ntrack;   //!
   TBranch        *b_sbs_hcal_adctime;   //!
   TBranch        *b_sbs_hcal_againblk;   //!
   TBranch        *b_sbs_hcal_atimeblk;   //!
   TBranch        *b_sbs_hcal_colblk;   //!
   TBranch        *b_sbs_hcal_colblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_e;   //!
   TBranch        *b_sbs_hcal_e_goodtdc;   //!
   TBranch        *b_sbs_hcal_eblk;   //!
   TBranch        *b_sbs_hcal_eblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_idblk;   //!
   TBranch        *b_sbs_hcal_idblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_index;   //!
   TBranch        *b_sbs_hcal_nblk;   //!
   TBranch        *b_sbs_hcal_nblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_nclus;   //!
   TBranch        *b_sbs_hcal_ngoodADChits;   //!
   TBranch        *b_sbs_hcal_rowblk;   //!
   TBranch        *b_sbs_hcal_rowblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_tdctime;   //!
   TBranch        *b_sbs_hcal_tdctime_tw;   //!
   TBranch        *b_sbs_hcal_tdctimeblk;   //!
   TBranch        *b_sbs_hcal_tdctimeblk_tw;   //!
   TBranch        *b_sbs_hcal_trigtime;   //!
   TBranch        *b_sbs_hcal_vtp_detid;   //!
   TBranch        *b_sbs_hcal_vtp_errflag;   //!
   TBranch        *b_sbs_hcal_x;   //!
   TBranch        *b_sbs_hcal_y;   //!
   TBranch        *b_sbs_tr_n;   //!
   TBranch        *b_scalhel_errcode;   //!
   TBranch        *b_scalhel_evtcount;   //!
   TBranch        *b_scalhel_hel;   //!
   TBranch        *b_scalhel_lhrs_fadc_hel;   //!
   TBranch        *b_scalhel_lhrs_fadc_pat;   //!
   TBranch        *b_scalhel_lhrs_fadc_tsettle;   //!
   TBranch        *b_scalhel_patcount;   //!
   TBranch        *b_scalhel_patphase;   //!
   TBranch        *b_scalhel_seed;   //!
   TBranch        *b_singletrack;   //!
   TBranch        *b_anytrack;   //!
   TBranch        *b_singletrack;   //!
   TBranch        *b_anytrack;   //!
   TBranch        *b_goodheep;   //!
   TBranch        *b_gooddxdy;   //!
   TBranch        *b_gooddpp;   //!
   TBranch        *b_goodctime;   //!
   TBranch        *b_HALLA_p;   //!
   TBranch        *b_hac_bcm_average;   //!
   TBranch        *b_IPM1H04A_XPOS;   //!
   TBranch        *b_IPM1H04A_YPOS;   //!
   TBranch        *b_IPM1H04E_XPOS;   //!
   TBranch        *b_IPM1H04E_YPOS;   //!
   TBranch        *b_IGL1I00OD16_16;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtTime;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtNum;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtType;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtLen;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fHelicity;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fTrigBits;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fRun;   //!

   gep_tree_data(TTree *tree=0);
   virtual ~gep_tree_data();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gep_tree_data_cxx
gep_tree_data::gep_tree_data(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("T",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("T","");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg0_0.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg10_10.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg11_11.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg1_1.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg2_2.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg3_3.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg4_4.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg5_5.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg6_6.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg7_7.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg8_8.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream0_0_seg9_9.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg0_0.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg10_10.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg11_11.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg1_1.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg2_2.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg3_3.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg4_4.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg5_5.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg6_6.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg7_7.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg8_8.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream1_1_seg9_9.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg0_0.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg10_10.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg11_11.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg1_1.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg2_2.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg3_3.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg4_4.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg5_5.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg6_6.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg7_7.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg8_8.root/T");
      chain->Add("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign/rootfiles/gep5_fullreplay_2947_stream2_2_seg9_9.root/T");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

gep_tree_data::~gep_tree_data()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gep_tree_data::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gep_tree_data::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gep_tree_data::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ax_new", &ax_new, &b_ax_new);
   fChain->SetBranchAddress("bx_new", &bx_new, &b_bx_new);
   fChain->SetBranchAddress("ay_new", &ay_new, &b_ay_new);
   fChain->SetBranchAddress("by_new", &by_new, &b_by_new);
   fChain->SetBranchAddress("targx", &targx, &b_targx);
   fChain->SetBranchAddress("targy", &targy, &b_targy);
   fChain->SetBranchAddress("scalhel.true_hel", &scalhel_true_hel, &b_scalhel_true_hel);
   fChain->SetBranchAddress("scalhel.hel_stable", &scalhel_hel_stable, &b_scalhel_hel_stable);
   fChain->SetBranchAddress("dEcalo", &dEcalo, &b_dEcalo);
   fChain->SetBranchAddress("detheta", &detheta, &b_detheta);
   fChain->SetBranchAddress("dptheta", &dptheta, &b_dptheta);
   fChain->SetBranchAddress("Ndata.earm.cdet.hit.col", &Ndata_earm_cdet_hit_col, &b_Ndata_earm_cdet_hit_col);
   fChain->SetBranchAddress("earm.cdet.hit.col", earm_cdet_hit_col, &b_earm_cdet_hit_col);
   fChain->SetBranchAddress("Ndata.earm.cdet.hit.layer", &Ndata_earm_cdet_hit_layer, &b_Ndata_earm_cdet_hit_layer);
   fChain->SetBranchAddress("earm.cdet.hit.layer", earm_cdet_hit_layer, &b_earm_cdet_hit_layer);
   fChain->SetBranchAddress("Ndata.earm.cdet.hit.pmtnum", &Ndata_earm_cdet_hit_pmtnum, &b_Ndata_earm_cdet_hit_pmtnum);
   fChain->SetBranchAddress("earm.cdet.hit.pmtnum", earm_cdet_hit_pmtnum, &b_earm_cdet_hit_pmtnum);
   fChain->SetBranchAddress("Ndata.earm.cdet.hit.row", &Ndata_earm_cdet_hit_row, &b_Ndata_earm_cdet_hit_row);
   fChain->SetBranchAddress("earm.cdet.hit.row", earm_cdet_hit_row, &b_earm_cdet_hit_row);
   fChain->SetBranchAddress("Ndata.earm.cdet.hit.tdc_le", &Ndata_earm_cdet_hit_tdc_le, &b_Ndata_earm_cdet_hit_tdc_le);
   fChain->SetBranchAddress("earm.cdet.hit.tdc_le", earm_cdet_hit_tdc_le, &b_earm_cdet_hit_tdc_le);
   fChain->SetBranchAddress("Ndata.earm.cdet.hit.tdc_te", &Ndata_earm_cdet_hit_tdc_te, &b_Ndata_earm_cdet_hit_tdc_te);
   fChain->SetBranchAddress("earm.cdet.hit.tdc_te", earm_cdet_hit_tdc_te, &b_earm_cdet_hit_tdc_te);
   fChain->SetBranchAddress("Ndata.earm.cdet.hit.tdc_tot", &Ndata_earm_cdet_hit_tdc_tot, &b_Ndata_earm_cdet_hit_tdc_tot);
   fChain->SetBranchAddress("earm.cdet.hit.tdc_tot", earm_cdet_hit_tdc_tot, &b_earm_cdet_hit_tdc_tot);
   fChain->SetBranchAddress("Ndata.earm.cdet.hit.xhit", &Ndata_earm_cdet_hit_xhit, &b_Ndata_earm_cdet_hit_xhit);
   fChain->SetBranchAddress("earm.cdet.hit.xhit", earm_cdet_hit_xhit, &b_earm_cdet_hit_xhit);
   fChain->SetBranchAddress("Ndata.earm.cdet.hit.yhit", &Ndata_earm_cdet_hit_yhit, &b_Ndata_earm_cdet_hit_yhit);
   fChain->SetBranchAddress("earm.cdet.hit.yhit", earm_cdet_hit_yhit, &b_earm_cdet_hit_yhit);
   fChain->SetBranchAddress("Ndata.earm.cdet.hit.zhit", &Ndata_earm_cdet_hit_zhit, &b_Ndata_earm_cdet_hit_zhit);
   fChain->SetBranchAddress("earm.cdet.hit.zhit", earm_cdet_hit_zhit, &b_earm_cdet_hit_zhit);
   fChain->SetBranchAddress("Ndata.earm.cdet.hits.TDCelemID", &Ndata_earm_cdet_hits_TDCelemID, &b_Ndata_earm_cdet_hits_TDCelemID);
   fChain->SetBranchAddress("earm.cdet.hits.TDCelemID", earm_cdet_hits_TDCelemID, &b_earm_cdet_hits_TDCelemID);
   fChain->SetBranchAddress("Ndata.earm.cdet.hits.t", &Ndata_earm_cdet_hits_t, &b_Ndata_earm_cdet_hits_t);
   fChain->SetBranchAddress("earm.cdet.hits.t", earm_cdet_hits_t, &b_earm_cdet_hits_t);
   fChain->SetBranchAddress("Ndata.earm.cdet.hits.t_te", &Ndata_earm_cdet_hits_t_te, &b_Ndata_earm_cdet_hits_t_te);
   fChain->SetBranchAddress("earm.cdet.hits.t_te", earm_cdet_hits_t_te, &b_earm_cdet_hits_t_te);
   fChain->SetBranchAddress("Ndata.earm.cdet.hits.t_tot", &Ndata_earm_cdet_hits_t_tot, &b_Ndata_earm_cdet_hits_t_tot);
   fChain->SetBranchAddress("earm.cdet.hits.t_tot", earm_cdet_hits_t_tot, &b_earm_cdet_hits_t_tot);
   fChain->SetBranchAddress("Ndata.earm.cdet.tdc", &Ndata_earm_cdet_tdc, &b_Ndata_earm_cdet_tdc);
   fChain->SetBranchAddress("earm.cdet.tdc", earm_cdet_tdc, &b_earm_cdet_tdc);
   fChain->SetBranchAddress("Ndata.earm.cdet.tdc_mult", &Ndata_earm_cdet_tdc_mult, &b_Ndata_earm_cdet_tdc_mult);
   fChain->SetBranchAddress("earm.cdet.tdc_mult", earm_cdet_tdc_mult, &b_earm_cdet_tdc_mult);
   fChain->SetBranchAddress("Ndata.earm.cdet.tdc_te", &Ndata_earm_cdet_tdc_te, &b_Ndata_earm_cdet_tdc_te);
   fChain->SetBranchAddress("earm.cdet.tdc_te", earm_cdet_tdc_te, &b_earm_cdet_tdc_te);
   fChain->SetBranchAddress("Ndata.earm.cdet.tdc_tot", &Ndata_earm_cdet_tdc_tot, &b_Ndata_earm_cdet_tdc_tot);
   fChain->SetBranchAddress("earm.cdet.tdc_tot", earm_cdet_tdc_tot, &b_earm_cdet_tdc_tot);
   fChain->SetBranchAddress("Ndata.earm.cdet.tdccol", &Ndata_earm_cdet_tdccol, &b_Ndata_earm_cdet_tdccol);
   fChain->SetBranchAddress("earm.cdet.tdccol", earm_cdet_tdccol, &b_earm_cdet_tdccol);
   fChain->SetBranchAddress("Ndata.earm.cdet.tdcelemID", &Ndata_earm_cdet_tdcelemID, &b_Ndata_earm_cdet_tdcelemID);
   fChain->SetBranchAddress("earm.cdet.tdcelemID", earm_cdet_tdcelemID, &b_earm_cdet_tdcelemID);
   fChain->SetBranchAddress("Ndata.earm.cdet.tdclayer", &Ndata_earm_cdet_tdclayer, &b_Ndata_earm_cdet_tdclayer);
   fChain->SetBranchAddress("earm.cdet.tdclayer", earm_cdet_tdclayer, &b_earm_cdet_tdclayer);
   fChain->SetBranchAddress("Ndata.earm.cdet.tdcrow", &Ndata_earm_cdet_tdcrow, &b_Ndata_earm_cdet_tdcrow);
   fChain->SetBranchAddress("earm.cdet.tdcrow", earm_cdet_tdcrow, &b_earm_cdet_tdcrow);
   fChain->SetBranchAddress("Ndata.earm.ecal.a_c", &Ndata_earm_ecal_a_c, &b_Ndata_earm_ecal_a_c);
   fChain->SetBranchAddress("earm.ecal.a_c", earm_ecal_a_c, &b_earm_ecal_a_c);
   fChain->SetBranchAddress("Ndata.earm.ecal.a_p", &Ndata_earm_ecal_a_p, &b_Ndata_earm_ecal_a_p);
   fChain->SetBranchAddress("earm.ecal.a_p", earm_ecal_a_p, &b_earm_ecal_a_p);
   fChain->SetBranchAddress("Ndata.earm.ecal.a_time", &Ndata_earm_ecal_a_time, &b_Ndata_earm_ecal_a_time);
   fChain->SetBranchAddress("earm.ecal.a_time", earm_ecal_a_time, &b_earm_ecal_a_time);
   fChain->SetBranchAddress("Ndata.earm.ecal.adccol", &Ndata_earm_ecal_adccol, &b_Ndata_earm_ecal_adccol);
   fChain->SetBranchAddress("earm.ecal.adccol", earm_ecal_adccol, &b_earm_ecal_adccol);
   fChain->SetBranchAddress("Ndata.earm.ecal.adcelemID", &Ndata_earm_ecal_adcelemID, &b_Ndata_earm_ecal_adcelemID);
   fChain->SetBranchAddress("earm.ecal.adcelemID", earm_ecal_adcelemID, &b_earm_ecal_adcelemID);
   fChain->SetBranchAddress("Ndata.earm.ecal.adcrow", &Ndata_earm_ecal_adcrow, &b_Ndata_earm_ecal_adcrow);
   fChain->SetBranchAddress("earm.ecal.adcrow", earm_ecal_adcrow, &b_earm_ecal_adcrow);
   fChain->SetBranchAddress("Ndata.earm.ecal.adcxpos", &Ndata_earm_ecal_adcxpos, &b_Ndata_earm_ecal_adcxpos);
   fChain->SetBranchAddress("earm.ecal.adcxpos", earm_ecal_adcxpos, &b_earm_ecal_adcxpos);
   fChain->SetBranchAddress("Ndata.earm.ecal.adcypos", &Ndata_earm_ecal_adcypos, &b_Ndata_earm_ecal_adcypos);
   fChain->SetBranchAddress("earm.ecal.adcypos", earm_ecal_adcypos, &b_earm_ecal_adcypos);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.atimeblk", &Ndata_earm_ecal_clus_atimeblk, &b_Ndata_earm_ecal_clus_atimeblk);
   fChain->SetBranchAddress("earm.ecal.clus.atimeblk", earm_ecal_clus_atimeblk, &b_earm_ecal_clus_atimeblk);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.col", &Ndata_earm_ecal_clus_col, &b_Ndata_earm_ecal_clus_col);
   fChain->SetBranchAddress("earm.ecal.clus.col", earm_ecal_clus_col, &b_earm_ecal_clus_col);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.e", &Ndata_earm_ecal_clus_e, &b_Ndata_earm_ecal_clus_e);
   fChain->SetBranchAddress("earm.ecal.clus.e", earm_ecal_clus_e, &b_earm_ecal_clus_e);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.eblk", &Ndata_earm_ecal_clus_eblk, &b_Ndata_earm_ecal_clus_eblk);
   fChain->SetBranchAddress("earm.ecal.clus.eblk", earm_ecal_clus_eblk, &b_earm_ecal_clus_eblk);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.id", &Ndata_earm_ecal_clus_id, &b_Ndata_earm_ecal_clus_id);
   fChain->SetBranchAddress("earm.ecal.clus.id", earm_ecal_clus_id, &b_earm_ecal_clus_id);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.nblk", &Ndata_earm_ecal_clus_nblk, &b_Ndata_earm_ecal_clus_nblk);
   fChain->SetBranchAddress("earm.ecal.clus.nblk", earm_ecal_clus_nblk, &b_earm_ecal_clus_nblk);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.row", &Ndata_earm_ecal_clus_row, &b_Ndata_earm_ecal_clus_row);
   fChain->SetBranchAddress("earm.ecal.clus.row", earm_ecal_clus_row, &b_earm_ecal_clus_row);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.x", &Ndata_earm_ecal_clus_x, &b_Ndata_earm_ecal_clus_x);
   fChain->SetBranchAddress("earm.ecal.clus.x", earm_ecal_clus_x, &b_earm_ecal_clus_x);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.y", &Ndata_earm_ecal_clus_y, &b_Ndata_earm_ecal_clus_y);
   fChain->SetBranchAddress("earm.ecal.clus.y", earm_ecal_clus_y, &b_earm_ecal_clus_y);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.atime", &Ndata_earm_ecal_clus_blk_atime, &b_Ndata_earm_ecal_clus_blk_atime);
   fChain->SetBranchAddress("earm.ecal.clus_blk.atime", earm_ecal_clus_blk_atime, &b_earm_ecal_clus_blk_atime);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.col", &Ndata_earm_ecal_clus_blk_col, &b_Ndata_earm_ecal_clus_blk_col);
   fChain->SetBranchAddress("earm.ecal.clus_blk.col", earm_ecal_clus_blk_col, &b_earm_ecal_clus_blk_col);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.e", &Ndata_earm_ecal_clus_blk_e, &b_Ndata_earm_ecal_clus_blk_e);
   fChain->SetBranchAddress("earm.ecal.clus_blk.e", earm_ecal_clus_blk_e, &b_earm_ecal_clus_blk_e);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.id", &Ndata_earm_ecal_clus_blk_id, &b_Ndata_earm_ecal_clus_blk_id);
   fChain->SetBranchAddress("earm.ecal.clus_blk.id", earm_ecal_clus_blk_id, &b_earm_ecal_clus_blk_id);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.row", &Ndata_earm_ecal_clus_blk_row, &b_Ndata_earm_ecal_clus_blk_row);
   fChain->SetBranchAddress("earm.ecal.clus_blk.row", earm_ecal_clus_blk_row, &b_earm_ecal_clus_blk_row);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.x", &Ndata_earm_ecal_clus_blk_x, &b_Ndata_earm_ecal_clus_blk_x);
   fChain->SetBranchAddress("earm.ecal.clus_blk.x", earm_ecal_clus_blk_x, &b_earm_ecal_clus_blk_x);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.y", &Ndata_earm_ecal_clus_blk_y, &b_Ndata_earm_ecal_clus_blk_y);
   fChain->SetBranchAddress("earm.ecal.clus_blk.y", earm_ecal_clus_blk_y, &b_earm_ecal_clus_blk_y);
   fChain->SetBranchAddress("Ndata.earm.ecal.ped", &Ndata_earm_ecal_ped, &b_Ndata_earm_ecal_ped);
   fChain->SetBranchAddress("earm.ecal.ped", earm_ecal_ped, &b_earm_ecal_ped);
   fChain->SetBranchAddress("Ndata.earm.ecal.vtp.clus.e", &Ndata_earm_ecal_vtp_clus_e, &b_Ndata_earm_ecal_vtp_clus_e);
   fChain->SetBranchAddress("earm.ecal.vtp.clus.e", earm_ecal_vtp_clus_e, &b_earm_ecal_vtp_clus_e);
   fChain->SetBranchAddress("Ndata.earm.ecal.vtp.clus.size", &Ndata_earm_ecal_vtp_clus_size, &b_Ndata_earm_ecal_vtp_clus_size);
   fChain->SetBranchAddress("earm.ecal.vtp.clus.size", earm_ecal_vtp_clus_size, &b_earm_ecal_vtp_clus_size);
   fChain->SetBranchAddress("Ndata.earm.ecal.vtp.clus.time", &Ndata_earm_ecal_vtp_clus_time, &b_Ndata_earm_ecal_vtp_clus_time);
   fChain->SetBranchAddress("earm.ecal.vtp.clus.time", earm_ecal_vtp_clus_time, &b_earm_ecal_vtp_clus_time);
   fChain->SetBranchAddress("Ndata.earm.ecal.vtp.clus.x", &Ndata_earm_ecal_vtp_clus_x, &b_Ndata_earm_ecal_vtp_clus_x);
   fChain->SetBranchAddress("earm.ecal.vtp.clus.x", earm_ecal_vtp_clus_x, &b_earm_ecal_vtp_clus_x);
   fChain->SetBranchAddress("Ndata.earm.ecal.vtp.clus.y", &Ndata_earm_ecal_vtp_clus_y, &b_Ndata_earm_ecal_vtp_clus_y);
   fChain->SetBranchAddress("earm.ecal.vtp.clus.y", earm_ecal_vtp_clus_y, &b_earm_ecal_vtp_clus_y);
   fChain->SetBranchAddress("Ndata.earm.tr.beta", &Ndata_earm_tr_beta, &b_Ndata_earm_tr_beta);
   fChain->SetBranchAddress("earm.tr.beta", earm_tr_beta, &b_earm_tr_beta);
   fChain->SetBranchAddress("Ndata.earm.tr.chi2", &Ndata_earm_tr_chi2, &b_Ndata_earm_tr_chi2);
   fChain->SetBranchAddress("earm.tr.chi2", earm_tr_chi2, &b_earm_tr_chi2);
   fChain->SetBranchAddress("Ndata.earm.tr.d_ph", &Ndata_earm_tr_d_ph, &b_Ndata_earm_tr_d_ph);
   fChain->SetBranchAddress("earm.tr.d_ph", earm_tr_d_ph, &b_earm_tr_d_ph);
   fChain->SetBranchAddress("Ndata.earm.tr.d_th", &Ndata_earm_tr_d_th, &b_Ndata_earm_tr_d_th);
   fChain->SetBranchAddress("earm.tr.d_th", earm_tr_d_th, &b_earm_tr_d_th);
   fChain->SetBranchAddress("Ndata.earm.tr.d_x", &Ndata_earm_tr_d_x, &b_Ndata_earm_tr_d_x);
   fChain->SetBranchAddress("earm.tr.d_x", earm_tr_d_x, &b_earm_tr_d_x);
   fChain->SetBranchAddress("Ndata.earm.tr.d_y", &Ndata_earm_tr_d_y, &b_Ndata_earm_tr_d_y);
   fChain->SetBranchAddress("earm.tr.d_y", earm_tr_d_y, &b_earm_tr_d_y);
   fChain->SetBranchAddress("Ndata.earm.tr.dbeta", &Ndata_earm_tr_dbeta, &b_Ndata_earm_tr_dbeta);
   fChain->SetBranchAddress("earm.tr.dbeta", earm_tr_dbeta, &b_earm_tr_dbeta);
   fChain->SetBranchAddress("Ndata.earm.tr.dtime", &Ndata_earm_tr_dtime, &b_Ndata_earm_tr_dtime);
   fChain->SetBranchAddress("earm.tr.dtime", earm_tr_dtime, &b_earm_tr_dtime);
   fChain->SetBranchAddress("Ndata.earm.tr.flag", &Ndata_earm_tr_flag, &b_Ndata_earm_tr_flag);
   fChain->SetBranchAddress("earm.tr.flag", earm_tr_flag, &b_earm_tr_flag);
   fChain->SetBranchAddress("Ndata.earm.tr.ndof", &Ndata_earm_tr_ndof, &b_Ndata_earm_tr_ndof);
   fChain->SetBranchAddress("earm.tr.ndof", earm_tr_ndof, &b_earm_tr_ndof);
   fChain->SetBranchAddress("Ndata.earm.tr.p", &Ndata_earm_tr_p, &b_Ndata_earm_tr_p);
   fChain->SetBranchAddress("earm.tr.p", earm_tr_p, &b_earm_tr_p);
   fChain->SetBranchAddress("Ndata.earm.tr.pathl", &Ndata_earm_tr_pathl, &b_Ndata_earm_tr_pathl);
   fChain->SetBranchAddress("earm.tr.pathl", earm_tr_pathl, &b_earm_tr_pathl);
   fChain->SetBranchAddress("Ndata.earm.tr.ph", &Ndata_earm_tr_ph, &b_Ndata_earm_tr_ph);
   fChain->SetBranchAddress("earm.tr.ph", earm_tr_ph, &b_earm_tr_ph);
   fChain->SetBranchAddress("Ndata.earm.tr.px", &Ndata_earm_tr_px, &b_Ndata_earm_tr_px);
   fChain->SetBranchAddress("earm.tr.px", earm_tr_px, &b_earm_tr_px);
   fChain->SetBranchAddress("Ndata.earm.tr.py", &Ndata_earm_tr_py, &b_Ndata_earm_tr_py);
   fChain->SetBranchAddress("earm.tr.py", earm_tr_py, &b_earm_tr_py);
   fChain->SetBranchAddress("Ndata.earm.tr.pz", &Ndata_earm_tr_pz, &b_Ndata_earm_tr_pz);
   fChain->SetBranchAddress("earm.tr.pz", earm_tr_pz, &b_earm_tr_pz);
   fChain->SetBranchAddress("Ndata.earm.tr.r_ph", &Ndata_earm_tr_r_ph, &b_Ndata_earm_tr_r_ph);
   fChain->SetBranchAddress("earm.tr.r_ph", earm_tr_r_ph, &b_earm_tr_r_ph);
   fChain->SetBranchAddress("Ndata.earm.tr.r_th", &Ndata_earm_tr_r_th, &b_Ndata_earm_tr_r_th);
   fChain->SetBranchAddress("earm.tr.r_th", earm_tr_r_th, &b_earm_tr_r_th);
   fChain->SetBranchAddress("Ndata.earm.tr.r_x", &Ndata_earm_tr_r_x, &b_Ndata_earm_tr_r_x);
   fChain->SetBranchAddress("earm.tr.r_x", earm_tr_r_x, &b_earm_tr_r_x);
   fChain->SetBranchAddress("Ndata.earm.tr.r_y", &Ndata_earm_tr_r_y, &b_Ndata_earm_tr_r_y);
   fChain->SetBranchAddress("earm.tr.r_y", earm_tr_r_y, &b_earm_tr_r_y);
   fChain->SetBranchAddress("Ndata.earm.tr.tg_dp", &Ndata_earm_tr_tg_dp, &b_Ndata_earm_tr_tg_dp);
   fChain->SetBranchAddress("earm.tr.tg_dp", earm_tr_tg_dp, &b_earm_tr_tg_dp);
   fChain->SetBranchAddress("Ndata.earm.tr.tg_ph", &Ndata_earm_tr_tg_ph, &b_Ndata_earm_tr_tg_ph);
   fChain->SetBranchAddress("earm.tr.tg_ph", earm_tr_tg_ph, &b_earm_tr_tg_ph);
   fChain->SetBranchAddress("Ndata.earm.tr.tg_th", &Ndata_earm_tr_tg_th, &b_Ndata_earm_tr_tg_th);
   fChain->SetBranchAddress("earm.tr.tg_th", earm_tr_tg_th, &b_earm_tr_tg_th);
   fChain->SetBranchAddress("Ndata.earm.tr.tg_y", &Ndata_earm_tr_tg_y, &b_Ndata_earm_tr_tg_y);
   fChain->SetBranchAddress("earm.tr.tg_y", earm_tr_tg_y, &b_earm_tr_tg_y);
   fChain->SetBranchAddress("Ndata.earm.tr.th", &Ndata_earm_tr_th, &b_Ndata_earm_tr_th);
   fChain->SetBranchAddress("earm.tr.th", earm_tr_th, &b_earm_tr_th);
   fChain->SetBranchAddress("Ndata.earm.tr.time", &Ndata_earm_tr_time, &b_Ndata_earm_tr_time);
   fChain->SetBranchAddress("earm.tr.time", earm_tr_time, &b_earm_tr_time);
   fChain->SetBranchAddress("Ndata.earm.tr.vx", &Ndata_earm_tr_vx, &b_Ndata_earm_tr_vx);
   fChain->SetBranchAddress("earm.tr.vx", earm_tr_vx, &b_earm_tr_vx);
   fChain->SetBranchAddress("Ndata.earm.tr.vy", &Ndata_earm_tr_vy, &b_Ndata_earm_tr_vy);
   fChain->SetBranchAddress("earm.tr.vy", earm_tr_vy, &b_earm_tr_vy);
   fChain->SetBranchAddress("Ndata.earm.tr.vz", &Ndata_earm_tr_vz, &b_Ndata_earm_tr_vz);
   fChain->SetBranchAddress("earm.tr.vz", earm_tr_vz, &b_earm_tr_vz);
   fChain->SetBranchAddress("Ndata.earm.tr.x", &Ndata_earm_tr_x, &b_Ndata_earm_tr_x);
   fChain->SetBranchAddress("earm.tr.x", earm_tr_x, &b_earm_tr_x);
   fChain->SetBranchAddress("Ndata.earm.tr.y", &Ndata_earm_tr_y, &b_Ndata_earm_tr_y);
   fChain->SetBranchAddress("earm.tr.y", earm_tr_y, &b_earm_tr_y);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCU", &Ndata_sbs_gemFPP_hit_ADCU, &b_Ndata_sbs_gemFPP_hit_ADCU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCU", sbs_gemFPP_hit_ADCU, &b_sbs_gemFPP_hit_ADCU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCU_deconv", &Ndata_sbs_gemFPP_hit_ADCU_deconv, &b_Ndata_sbs_gemFPP_hit_ADCU_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCU_deconv", sbs_gemFPP_hit_ADCU_deconv, &b_sbs_gemFPP_hit_ADCU_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCV", &Ndata_sbs_gemFPP_hit_ADCV, &b_Ndata_sbs_gemFPP_hit_ADCV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCV", sbs_gemFPP_hit_ADCV, &b_sbs_gemFPP_hit_ADCV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCV_deconv", &Ndata_sbs_gemFPP_hit_ADCV_deconv, &b_Ndata_sbs_gemFPP_hit_ADCV_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCV_deconv", sbs_gemFPP_hit_ADCV_deconv, &b_sbs_gemFPP_hit_ADCV_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCasym", &Ndata_sbs_gemFPP_hit_ADCasym, &b_Ndata_sbs_gemFPP_hit_ADCasym);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCasym", sbs_gemFPP_hit_ADCasym, &b_sbs_gemFPP_hit_ADCasym);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCasym_deconv", &Ndata_sbs_gemFPP_hit_ADCasym_deconv, &b_Ndata_sbs_gemFPP_hit_ADCasym_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCasym_deconv", sbs_gemFPP_hit_ADCasym_deconv, &b_sbs_gemFPP_hit_ADCasym_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCavg", &Ndata_sbs_gemFPP_hit_ADCavg, &b_Ndata_sbs_gemFPP_hit_ADCavg);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCavg", sbs_gemFPP_hit_ADCavg, &b_sbs_gemFPP_hit_ADCavg);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCavg_deconv", &Ndata_sbs_gemFPP_hit_ADCavg_deconv, &b_Ndata_sbs_gemFPP_hit_ADCavg_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCavg_deconv", sbs_gemFPP_hit_ADCavg_deconv, &b_sbs_gemFPP_hit_ADCavg_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac0_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac0_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac0_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac0_Umax", sbs_gemFPP_hit_ADCfrac0_Umax, &b_sbs_gemFPP_hit_ADCfrac0_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac0_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac0_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac0_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac0_Vmax", sbs_gemFPP_hit_ADCfrac0_Vmax, &b_sbs_gemFPP_hit_ADCfrac0_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac1_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac1_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac1_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac1_Umax", sbs_gemFPP_hit_ADCfrac1_Umax, &b_sbs_gemFPP_hit_ADCfrac1_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac1_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac1_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac1_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac1_Vmax", sbs_gemFPP_hit_ADCfrac1_Vmax, &b_sbs_gemFPP_hit_ADCfrac1_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac2_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac2_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac2_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac2_Umax", sbs_gemFPP_hit_ADCfrac2_Umax, &b_sbs_gemFPP_hit_ADCfrac2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac2_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac2_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac2_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac2_Vmax", sbs_gemFPP_hit_ADCfrac2_Vmax, &b_sbs_gemFPP_hit_ADCfrac2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac3_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac3_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac3_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac3_Umax", sbs_gemFPP_hit_ADCfrac3_Umax, &b_sbs_gemFPP_hit_ADCfrac3_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac3_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac3_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac3_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac3_Vmax", sbs_gemFPP_hit_ADCfrac3_Vmax, &b_sbs_gemFPP_hit_ADCfrac3_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac4_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac4_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac4_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac4_Umax", sbs_gemFPP_hit_ADCfrac4_Umax, &b_sbs_gemFPP_hit_ADCfrac4_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac4_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac4_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac4_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac4_Vmax", sbs_gemFPP_hit_ADCfrac4_Vmax, &b_sbs_gemFPP_hit_ADCfrac4_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac5_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac5_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac5_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac5_Umax", sbs_gemFPP_hit_ADCfrac5_Umax, &b_sbs_gemFPP_hit_ADCfrac5_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac5_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac5_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac5_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac5_Vmax", sbs_gemFPP_hit_ADCfrac5_Vmax, &b_sbs_gemFPP_hit_ADCfrac5_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxcomboUclust_deconv", &Ndata_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv, &b_Ndata_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxcomboUclust_deconv", sbs_gemFPP_hit_ADCmaxcomboUclust_deconv, &b_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxcomboVclust_deconv", &Ndata_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv, &b_Ndata_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxcomboVclust_deconv", sbs_gemFPP_hit_ADCmaxcomboVclust_deconv, &b_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampU", &Ndata_sbs_gemFPP_hit_ADCmaxsampU, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampU", sbs_gemFPP_hit_ADCmaxsampU, &b_sbs_gemFPP_hit_ADCmaxsampU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampUclust", &Ndata_sbs_gemFPP_hit_ADCmaxsampUclust, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampUclust", sbs_gemFPP_hit_ADCmaxsampUclust, &b_sbs_gemFPP_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampUclust_deconv", &Ndata_sbs_gemFPP_hit_ADCmaxsampUclust_deconv, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampUclust_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampUclust_deconv", sbs_gemFPP_hit_ADCmaxsampUclust_deconv, &b_sbs_gemFPP_hit_ADCmaxsampUclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampV", &Ndata_sbs_gemFPP_hit_ADCmaxsampV, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampV", sbs_gemFPP_hit_ADCmaxsampV, &b_sbs_gemFPP_hit_ADCmaxsampV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampVclust", &Ndata_sbs_gemFPP_hit_ADCmaxsampVclust, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampVclust", sbs_gemFPP_hit_ADCmaxsampVclust, &b_sbs_gemFPP_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampVclust_deconv", &Ndata_sbs_gemFPP_hit_ADCmaxsampVclust_deconv, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampVclust_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampVclust_deconv", sbs_gemFPP_hit_ADCmaxsampVclust_deconv, &b_sbs_gemFPP_hit_ADCmaxsampVclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxstripU", &Ndata_sbs_gemFPP_hit_ADCmaxstripU, &b_Ndata_sbs_gemFPP_hit_ADCmaxstripU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxstripU", sbs_gemFPP_hit_ADCmaxstripU, &b_sbs_gemFPP_hit_ADCmaxstripU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxstripV", &Ndata_sbs_gemFPP_hit_ADCmaxstripV, &b_Ndata_sbs_gemFPP_hit_ADCmaxstripV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxstripV", sbs_gemFPP_hit_ADCmaxstripV, &b_sbs_gemFPP_hit_ADCmaxstripV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.BUILD_ALL_SAMPLES_U", &Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U, &b_Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.BUILD_ALL_SAMPLES_U", sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U, &b_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.BUILD_ALL_SAMPLES_V", &Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V, &b_Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.BUILD_ALL_SAMPLES_V", sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V, &b_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.CM_GOOD_U", &Ndata_sbs_gemFPP_hit_CM_GOOD_U, &b_Ndata_sbs_gemFPP_hit_CM_GOOD_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.CM_GOOD_U", sbs_gemFPP_hit_CM_GOOD_U, &b_sbs_gemFPP_hit_CM_GOOD_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.CM_GOOD_V", &Ndata_sbs_gemFPP_hit_CM_GOOD_V, &b_Ndata_sbs_gemFPP_hit_CM_GOOD_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.CM_GOOD_V", sbs_gemFPP_hit_CM_GOOD_V, &b_sbs_gemFPP_hit_CM_GOOD_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC0_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC0_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC0_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC0_Umax", sbs_gemFPP_hit_DeconvADC0_Umax, &b_sbs_gemFPP_hit_DeconvADC0_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC0_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC0_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC0_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC0_Vmax", sbs_gemFPP_hit_DeconvADC0_Vmax, &b_sbs_gemFPP_hit_DeconvADC0_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC1_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC1_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC1_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC1_Umax", sbs_gemFPP_hit_DeconvADC1_Umax, &b_sbs_gemFPP_hit_DeconvADC1_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC1_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC1_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC1_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC1_Vmax", sbs_gemFPP_hit_DeconvADC1_Vmax, &b_sbs_gemFPP_hit_DeconvADC1_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC2_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC2_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC2_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC2_Umax", sbs_gemFPP_hit_DeconvADC2_Umax, &b_sbs_gemFPP_hit_DeconvADC2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC2_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC2_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC2_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC2_Vmax", sbs_gemFPP_hit_DeconvADC2_Vmax, &b_sbs_gemFPP_hit_DeconvADC2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC3_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC3_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC3_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC3_Umax", sbs_gemFPP_hit_DeconvADC3_Umax, &b_sbs_gemFPP_hit_DeconvADC3_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC3_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC3_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC3_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC3_Vmax", sbs_gemFPP_hit_DeconvADC3_Vmax, &b_sbs_gemFPP_hit_DeconvADC3_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC4_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC4_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC4_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC4_Umax", sbs_gemFPP_hit_DeconvADC4_Umax, &b_sbs_gemFPP_hit_DeconvADC4_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC4_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC4_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC4_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC4_Vmax", sbs_gemFPP_hit_DeconvADC4_Vmax, &b_sbs_gemFPP_hit_DeconvADC4_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC5_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC5_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC5_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC5_Umax", sbs_gemFPP_hit_DeconvADC5_Umax, &b_sbs_gemFPP_hit_DeconvADC5_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC5_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC5_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC5_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC5_Vmax", sbs_gemFPP_hit_DeconvADC5_Vmax, &b_sbs_gemFPP_hit_DeconvADC5_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxcomboU", &Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboU, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxcomboU", sbs_gemFPP_hit_DeconvADCmaxcomboU, &b_sbs_gemFPP_hit_DeconvADCmaxcomboU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxcomboV", &Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboV, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxcomboV", sbs_gemFPP_hit_DeconvADCmaxcomboV, &b_sbs_gemFPP_hit_DeconvADCmaxcomboV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxsampU", &Ndata_sbs_gemFPP_hit_DeconvADCmaxsampU, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxsampU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxsampU", sbs_gemFPP_hit_DeconvADCmaxsampU, &b_sbs_gemFPP_hit_DeconvADCmaxsampU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxsampV", &Ndata_sbs_gemFPP_hit_DeconvADCmaxsampV, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxsampV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxsampV", sbs_gemFPP_hit_DeconvADCmaxsampV, &b_sbs_gemFPP_hit_DeconvADCmaxsampV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxstripU", &Ndata_sbs_gemFPP_hit_DeconvADCmaxstripU, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxstripU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxstripU", sbs_gemFPP_hit_DeconvADCmaxstripU, &b_sbs_gemFPP_hit_DeconvADCmaxstripU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxstripV", &Ndata_sbs_gemFPP_hit_DeconvADCmaxstripV, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxstripV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxstripV", sbs_gemFPP_hit_DeconvADCmaxstripV, &b_sbs_gemFPP_hit_DeconvADCmaxstripV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ENABLE_CM_U", &Ndata_sbs_gemFPP_hit_ENABLE_CM_U, &b_Ndata_sbs_gemFPP_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ENABLE_CM_U", sbs_gemFPP_hit_ENABLE_CM_U, &b_sbs_gemFPP_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ENABLE_CM_V", &Ndata_sbs_gemFPP_hit_ENABLE_CM_V, &b_Ndata_sbs_gemFPP_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ENABLE_CM_V", sbs_gemFPP_hit_ENABLE_CM_V, &b_sbs_gemFPP_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.TSchi2_Umax", &Ndata_sbs_gemFPP_hit_TSchi2_Umax, &b_Ndata_sbs_gemFPP_hit_TSchi2_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.TSchi2_Umax", sbs_gemFPP_hit_TSchi2_Umax, &b_sbs_gemFPP_hit_TSchi2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.TSchi2_Vmax", &Ndata_sbs_gemFPP_hit_TSchi2_Vmax, &b_Ndata_sbs_gemFPP_hit_TSchi2_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.TSchi2_Vmax", sbs_gemFPP_hit_TSchi2_Vmax, &b_sbs_gemFPP_hit_TSchi2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.TSprob_Umax", &Ndata_sbs_gemFPP_hit_TSprob_Umax, &b_Ndata_sbs_gemFPP_hit_TSprob_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.TSprob_Umax", sbs_gemFPP_hit_TSprob_Umax, &b_sbs_gemFPP_hit_TSprob_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.TSprob_Vmax", &Ndata_sbs_gemFPP_hit_TSprob_Vmax, &b_Ndata_sbs_gemFPP_hit_TSprob_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.TSprob_Vmax", sbs_gemFPP_hit_TSprob_Vmax, &b_sbs_gemFPP_hit_TSprob_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Tavg", &Ndata_sbs_gemFPP_hit_Tavg, &b_Ndata_sbs_gemFPP_hit_Tavg);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Tavg", sbs_gemFPP_hit_Tavg, &b_sbs_gemFPP_hit_Tavg);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Tavg_corr", &Ndata_sbs_gemFPP_hit_Tavg_corr, &b_Ndata_sbs_gemFPP_hit_Tavg_corr);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Tavg_corr", sbs_gemFPP_hit_Tavg_corr, &b_sbs_gemFPP_hit_Tavg_corr);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Tavg_deconv", &Ndata_sbs_gemFPP_hit_Tavg_deconv, &b_Ndata_sbs_gemFPP_hit_Tavg_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Tavg_deconv", sbs_gemFPP_hit_Tavg_deconv, &b_sbs_gemFPP_hit_Tavg_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Tavg_fit", &Ndata_sbs_gemFPP_hit_Tavg_fit, &b_Ndata_sbs_gemFPP_hit_Tavg_fit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Tavg_fit", sbs_gemFPP_hit_Tavg_fit, &b_sbs_gemFPP_hit_Tavg_fit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Ugain", &Ndata_sbs_gemFPP_hit_Ugain, &b_Ndata_sbs_gemFPP_hit_Ugain);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Ugain", sbs_gemFPP_hit_Ugain, &b_sbs_gemFPP_hit_Ugain);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Utime", &Ndata_sbs_gemFPP_hit_Utime, &b_Ndata_sbs_gemFPP_hit_Utime);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Utime", sbs_gemFPP_hit_Utime, &b_sbs_gemFPP_hit_Utime);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.UtimeDeconv", &Ndata_sbs_gemFPP_hit_UtimeDeconv, &b_Ndata_sbs_gemFPP_hit_UtimeDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.UtimeDeconv", sbs_gemFPP_hit_UtimeDeconv, &b_sbs_gemFPP_hit_UtimeDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.UtimeFit", &Ndata_sbs_gemFPP_hit_UtimeFit, &b_Ndata_sbs_gemFPP_hit_UtimeFit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.UtimeFit", sbs_gemFPP_hit_UtimeFit, &b_sbs_gemFPP_hit_UtimeFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.UtimeMaxStrip", &Ndata_sbs_gemFPP_hit_UtimeMaxStrip, &b_Ndata_sbs_gemFPP_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("sbs.gemFPP.hit.UtimeMaxStrip", sbs_gemFPP_hit_UtimeMaxStrip, &b_sbs_gemFPP_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.UtimeMaxStripDeconv", &Ndata_sbs_gemFPP_hit_UtimeMaxStripDeconv, &b_Ndata_sbs_gemFPP_hit_UtimeMaxStripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.UtimeMaxStripDeconv", sbs_gemFPP_hit_UtimeMaxStripDeconv, &b_sbs_gemFPP_hit_UtimeMaxStripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.UtimeMaxStripFit", &Ndata_sbs_gemFPP_hit_UtimeMaxStripFit, &b_Ndata_sbs_gemFPP_hit_UtimeMaxStripFit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.UtimeMaxStripFit", sbs_gemFPP_hit_UtimeMaxStripFit, &b_sbs_gemFPP_hit_UtimeMaxStripFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Vgain", &Ndata_sbs_gemFPP_hit_Vgain, &b_Ndata_sbs_gemFPP_hit_Vgain);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Vgain", sbs_gemFPP_hit_Vgain, &b_sbs_gemFPP_hit_Vgain);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Vtime", &Ndata_sbs_gemFPP_hit_Vtime, &b_Ndata_sbs_gemFPP_hit_Vtime);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Vtime", sbs_gemFPP_hit_Vtime, &b_sbs_gemFPP_hit_Vtime);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.VtimeDeconv", &Ndata_sbs_gemFPP_hit_VtimeDeconv, &b_Ndata_sbs_gemFPP_hit_VtimeDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.VtimeDeconv", sbs_gemFPP_hit_VtimeDeconv, &b_sbs_gemFPP_hit_VtimeDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.VtimeFit", &Ndata_sbs_gemFPP_hit_VtimeFit, &b_Ndata_sbs_gemFPP_hit_VtimeFit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.VtimeFit", sbs_gemFPP_hit_VtimeFit, &b_sbs_gemFPP_hit_VtimeFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.VtimeMaxStrip", &Ndata_sbs_gemFPP_hit_VtimeMaxStrip, &b_Ndata_sbs_gemFPP_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("sbs.gemFPP.hit.VtimeMaxStrip", sbs_gemFPP_hit_VtimeMaxStrip, &b_sbs_gemFPP_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.VtimeMaxStripDeconv", &Ndata_sbs_gemFPP_hit_VtimeMaxStripDeconv, &b_Ndata_sbs_gemFPP_hit_VtimeMaxStripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.VtimeMaxStripDeconv", sbs_gemFPP_hit_VtimeMaxStripDeconv, &b_sbs_gemFPP_hit_VtimeMaxStripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.VtimeMaxStripFit", &Ndata_sbs_gemFPP_hit_VtimeMaxStripFit, &b_Ndata_sbs_gemFPP_hit_VtimeMaxStripFit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.VtimeMaxStripFit", sbs_gemFPP_hit_VtimeMaxStripFit, &b_sbs_gemFPP_hit_VtimeMaxStripFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.adc_id_U", &Ndata_sbs_gemFPP_hit_adc_id_U, &b_Ndata_sbs_gemFPP_hit_adc_id_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.adc_id_U", sbs_gemFPP_hit_adc_id_U, &b_sbs_gemFPP_hit_adc_id_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.adc_id_V", &Ndata_sbs_gemFPP_hit_adc_id_V, &b_Ndata_sbs_gemFPP_hit_adc_id_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.adc_id_V", sbs_gemFPP_hit_adc_id_V, &b_sbs_gemFPP_hit_adc_id_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ccor_clust", &Ndata_sbs_gemFPP_hit_ccor_clust, &b_Ndata_sbs_gemFPP_hit_ccor_clust);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ccor_clust", sbs_gemFPP_hit_ccor_clust, &b_sbs_gemFPP_hit_ccor_clust);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ccor_clust_deconv", &Ndata_sbs_gemFPP_hit_ccor_clust_deconv, &b_Ndata_sbs_gemFPP_hit_ccor_clust_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ccor_clust_deconv", sbs_gemFPP_hit_ccor_clust_deconv, &b_sbs_gemFPP_hit_ccor_clust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ccor_strip", &Ndata_sbs_gemFPP_hit_ccor_strip, &b_Ndata_sbs_gemFPP_hit_ccor_strip);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ccor_strip", sbs_gemFPP_hit_ccor_strip, &b_sbs_gemFPP_hit_ccor_strip);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ccor_strip_deconv", &Ndata_sbs_gemFPP_hit_ccor_strip_deconv, &b_Ndata_sbs_gemFPP_hit_ccor_strip_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ccor_strip_deconv", sbs_gemFPP_hit_ccor_strip_deconv, &b_sbs_gemFPP_hit_ccor_strip_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.crate_U", &Ndata_sbs_gemFPP_hit_crate_U, &b_Ndata_sbs_gemFPP_hit_crate_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.crate_U", sbs_gemFPP_hit_crate_U, &b_sbs_gemFPP_hit_crate_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.crate_V", &Ndata_sbs_gemFPP_hit_crate_V, &b_Ndata_sbs_gemFPP_hit_crate_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.crate_V", sbs_gemFPP_hit_crate_V, &b_sbs_gemFPP_hit_crate_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.deltat", &Ndata_sbs_gemFPP_hit_deltat, &b_Ndata_sbs_gemFPP_hit_deltat);
   fChain->SetBranchAddress("sbs.gemFPP.hit.deltat", sbs_gemFPP_hit_deltat, &b_sbs_gemFPP_hit_deltat);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.deltat_deconv", &Ndata_sbs_gemFPP_hit_deltat_deconv, &b_Ndata_sbs_gemFPP_hit_deltat_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.deltat_deconv", sbs_gemFPP_hit_deltat_deconv, &b_sbs_gemFPP_hit_deltat_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.deltat_fit", &Ndata_sbs_gemFPP_hit_deltat_fit, &b_Ndata_sbs_gemFPP_hit_deltat_fit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.deltat_fit", sbs_gemFPP_hit_deltat_fit, &b_sbs_gemFPP_hit_deltat_fit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.eresidu", &Ndata_sbs_gemFPP_hit_eresidu, &b_Ndata_sbs_gemFPP_hit_eresidu);
   fChain->SetBranchAddress("sbs.gemFPP.hit.eresidu", sbs_gemFPP_hit_eresidu, &b_sbs_gemFPP_hit_eresidu);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.eresidv", &Ndata_sbs_gemFPP_hit_eresidv, &b_Ndata_sbs_gemFPP_hit_eresidv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.eresidv", sbs_gemFPP_hit_eresidv, &b_sbs_gemFPP_hit_eresidv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.icombomaxUclustDeconv", &Ndata_sbs_gemFPP_hit_icombomaxUclustDeconv, &b_Ndata_sbs_gemFPP_hit_icombomaxUclustDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.icombomaxUclustDeconv", sbs_gemFPP_hit_icombomaxUclustDeconv, &b_sbs_gemFPP_hit_icombomaxUclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.icombomaxUstripDeconv", &Ndata_sbs_gemFPP_hit_icombomaxUstripDeconv, &b_Ndata_sbs_gemFPP_hit_icombomaxUstripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.icombomaxUstripDeconv", sbs_gemFPP_hit_icombomaxUstripDeconv, &b_sbs_gemFPP_hit_icombomaxUstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.icombomaxVclustDeconv", &Ndata_sbs_gemFPP_hit_icombomaxVclustDeconv, &b_Ndata_sbs_gemFPP_hit_icombomaxVclustDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.icombomaxVclustDeconv", sbs_gemFPP_hit_icombomaxVclustDeconv, &b_sbs_gemFPP_hit_icombomaxVclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.icombomaxVstripDeconv", &Ndata_sbs_gemFPP_hit_icombomaxVstripDeconv, &b_Ndata_sbs_gemFPP_hit_icombomaxVstripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.icombomaxVstripDeconv", sbs_gemFPP_hit_icombomaxVstripDeconv, &b_sbs_gemFPP_hit_icombomaxVstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxUclust", &Ndata_sbs_gemFPP_hit_isampmaxUclust, &b_Ndata_sbs_gemFPP_hit_isampmaxUclust);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxUclust", sbs_gemFPP_hit_isampmaxUclust, &b_sbs_gemFPP_hit_isampmaxUclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxUclustDeconv", &Ndata_sbs_gemFPP_hit_isampmaxUclustDeconv, &b_Ndata_sbs_gemFPP_hit_isampmaxUclustDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxUclustDeconv", sbs_gemFPP_hit_isampmaxUclustDeconv, &b_sbs_gemFPP_hit_isampmaxUclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxUstrip", &Ndata_sbs_gemFPP_hit_isampmaxUstrip, &b_Ndata_sbs_gemFPP_hit_isampmaxUstrip);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxUstrip", sbs_gemFPP_hit_isampmaxUstrip, &b_sbs_gemFPP_hit_isampmaxUstrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxUstripDeconv", &Ndata_sbs_gemFPP_hit_isampmaxUstripDeconv, &b_Ndata_sbs_gemFPP_hit_isampmaxUstripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxUstripDeconv", sbs_gemFPP_hit_isampmaxUstripDeconv, &b_sbs_gemFPP_hit_isampmaxUstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxVclust", &Ndata_sbs_gemFPP_hit_isampmaxVclust, &b_Ndata_sbs_gemFPP_hit_isampmaxVclust);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxVclust", sbs_gemFPP_hit_isampmaxVclust, &b_sbs_gemFPP_hit_isampmaxVclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxVclustDeconv", &Ndata_sbs_gemFPP_hit_isampmaxVclustDeconv, &b_Ndata_sbs_gemFPP_hit_isampmaxVclustDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxVclustDeconv", sbs_gemFPP_hit_isampmaxVclustDeconv, &b_sbs_gemFPP_hit_isampmaxVclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxVstrip", &Ndata_sbs_gemFPP_hit_isampmaxVstrip, &b_Ndata_sbs_gemFPP_hit_isampmaxVstrip);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxVstrip", sbs_gemFPP_hit_isampmaxVstrip, &b_sbs_gemFPP_hit_isampmaxVstrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxVstripDeconv", &Ndata_sbs_gemFPP_hit_isampmaxVstripDeconv, &b_Ndata_sbs_gemFPP_hit_isampmaxVstripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxVstripDeconv", sbs_gemFPP_hit_isampmaxVstripDeconv, &b_sbs_gemFPP_hit_isampmaxVstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.layer", &Ndata_sbs_gemFPP_hit_layer, &b_Ndata_sbs_gemFPP_hit_layer);
   fChain->SetBranchAddress("sbs.gemFPP.hit.layer", sbs_gemFPP_hit_layer, &b_sbs_gemFPP_hit_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.module", &Ndata_sbs_gemFPP_hit_module, &b_Ndata_sbs_gemFPP_hit_module);
   fChain->SetBranchAddress("sbs.gemFPP.hit.module", sbs_gemFPP_hit_module, &b_sbs_gemFPP_hit_module);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.mpd_U", &Ndata_sbs_gemFPP_hit_mpd_U, &b_Ndata_sbs_gemFPP_hit_mpd_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.mpd_U", sbs_gemFPP_hit_mpd_U, &b_sbs_gemFPP_hit_mpd_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.mpd_V", &Ndata_sbs_gemFPP_hit_mpd_V, &b_Ndata_sbs_gemFPP_hit_mpd_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.mpd_V", sbs_gemFPP_hit_mpd_V, &b_sbs_gemFPP_hit_mpd_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.nstripu", &Ndata_sbs_gemFPP_hit_nstripu, &b_Ndata_sbs_gemFPP_hit_nstripu);
   fChain->SetBranchAddress("sbs.gemFPP.hit.nstripu", sbs_gemFPP_hit_nstripu, &b_sbs_gemFPP_hit_nstripu);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.nstripv", &Ndata_sbs_gemFPP_hit_nstripv, &b_Ndata_sbs_gemFPP_hit_nstripv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.nstripv", sbs_gemFPP_hit_nstripv, &b_sbs_gemFPP_hit_nstripv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.residu", &Ndata_sbs_gemFPP_hit_residu, &b_Ndata_sbs_gemFPP_hit_residu);
   fChain->SetBranchAddress("sbs.gemFPP.hit.residu", sbs_gemFPP_hit_residu, &b_sbs_gemFPP_hit_residu);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.residv", &Ndata_sbs_gemFPP_hit_residv, &b_Ndata_sbs_gemFPP_hit_residv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.residv", sbs_gemFPP_hit_residv, &b_sbs_gemFPP_hit_residv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.trackindex", &Ndata_sbs_gemFPP_hit_trackindex, &b_Ndata_sbs_gemFPP_hit_trackindex);
   fChain->SetBranchAddress("sbs.gemFPP.hit.trackindex", sbs_gemFPP_hit_trackindex, &b_sbs_gemFPP_hit_trackindex);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.u", &Ndata_sbs_gemFPP_hit_u, &b_Ndata_sbs_gemFPP_hit_u);
   fChain->SetBranchAddress("sbs.gemFPP.hit.u", sbs_gemFPP_hit_u, &b_sbs_gemFPP_hit_u);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.umoment", &Ndata_sbs_gemFPP_hit_umoment, &b_Ndata_sbs_gemFPP_hit_umoment);
   fChain->SetBranchAddress("sbs.gemFPP.hit.umoment", sbs_gemFPP_hit_umoment, &b_sbs_gemFPP_hit_umoment);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.usigma", &Ndata_sbs_gemFPP_hit_usigma, &b_Ndata_sbs_gemFPP_hit_usigma);
   fChain->SetBranchAddress("sbs.gemFPP.hit.usigma", sbs_gemFPP_hit_usigma, &b_sbs_gemFPP_hit_usigma);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ustriphi", &Ndata_sbs_gemFPP_hit_ustriphi, &b_Ndata_sbs_gemFPP_hit_ustriphi);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ustriphi", sbs_gemFPP_hit_ustriphi, &b_sbs_gemFPP_hit_ustriphi);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ustriplo", &Ndata_sbs_gemFPP_hit_ustriplo, &b_Ndata_sbs_gemFPP_hit_ustriplo);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ustriplo", sbs_gemFPP_hit_ustriplo, &b_sbs_gemFPP_hit_ustriplo);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ustripmax", &Ndata_sbs_gemFPP_hit_ustripmax, &b_Ndata_sbs_gemFPP_hit_ustripmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ustripmax", sbs_gemFPP_hit_ustripmax, &b_sbs_gemFPP_hit_ustripmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.v", &Ndata_sbs_gemFPP_hit_v, &b_Ndata_sbs_gemFPP_hit_v);
   fChain->SetBranchAddress("sbs.gemFPP.hit.v", sbs_gemFPP_hit_v, &b_sbs_gemFPP_hit_v);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.vmoment", &Ndata_sbs_gemFPP_hit_vmoment, &b_Ndata_sbs_gemFPP_hit_vmoment);
   fChain->SetBranchAddress("sbs.gemFPP.hit.vmoment", sbs_gemFPP_hit_vmoment, &b_sbs_gemFPP_hit_vmoment);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.vsigma", &Ndata_sbs_gemFPP_hit_vsigma, &b_Ndata_sbs_gemFPP_hit_vsigma);
   fChain->SetBranchAddress("sbs.gemFPP.hit.vsigma", sbs_gemFPP_hit_vsigma, &b_sbs_gemFPP_hit_vsigma);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.vstriphi", &Ndata_sbs_gemFPP_hit_vstriphi, &b_Ndata_sbs_gemFPP_hit_vstriphi);
   fChain->SetBranchAddress("sbs.gemFPP.hit.vstriphi", sbs_gemFPP_hit_vstriphi, &b_sbs_gemFPP_hit_vstriphi);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.vstriplo", &Ndata_sbs_gemFPP_hit_vstriplo, &b_Ndata_sbs_gemFPP_hit_vstriplo);
   fChain->SetBranchAddress("sbs.gemFPP.hit.vstriplo", sbs_gemFPP_hit_vstriplo, &b_sbs_gemFPP_hit_vstriplo);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.vstripmax", &Ndata_sbs_gemFPP_hit_vstripmax, &b_Ndata_sbs_gemFPP_hit_vstripmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.vstripmax", sbs_gemFPP_hit_vstripmax, &b_sbs_gemFPP_hit_vstripmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.xglobal", &Ndata_sbs_gemFPP_hit_xglobal, &b_Ndata_sbs_gemFPP_hit_xglobal);
   fChain->SetBranchAddress("sbs.gemFPP.hit.xglobal", sbs_gemFPP_hit_xglobal, &b_sbs_gemFPP_hit_xglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.xlocal", &Ndata_sbs_gemFPP_hit_xlocal, &b_Ndata_sbs_gemFPP_hit_xlocal);
   fChain->SetBranchAddress("sbs.gemFPP.hit.xlocal", sbs_gemFPP_hit_xlocal, &b_sbs_gemFPP_hit_xlocal);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.yglobal", &Ndata_sbs_gemFPP_hit_yglobal, &b_Ndata_sbs_gemFPP_hit_yglobal);
   fChain->SetBranchAddress("sbs.gemFPP.hit.yglobal", sbs_gemFPP_hit_yglobal, &b_sbs_gemFPP_hit_yglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ylocal", &Ndata_sbs_gemFPP_hit_ylocal, &b_Ndata_sbs_gemFPP_hit_ylocal);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ylocal", sbs_gemFPP_hit_ylocal, &b_sbs_gemFPP_hit_ylocal);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.zglobal", &Ndata_sbs_gemFPP_hit_zglobal, &b_Ndata_sbs_gemFPP_hit_zglobal);
   fChain->SetBranchAddress("sbs.gemFPP.hit.zglobal", sbs_gemFPP_hit_zglobal, &b_sbs_gemFPP_hit_zglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.n2Dhit_layer", &Ndata_sbs_gemFPP_n2Dhit_layer, &b_Ndata_sbs_gemFPP_n2Dhit_layer);
   fChain->SetBranchAddress("sbs.gemFPP.n2Dhit_layer", sbs_gemFPP_n2Dhit_layer, &b_sbs_gemFPP_n2Dhit_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.nclustu_layer", &Ndata_sbs_gemFPP_nclustu_layer, &b_Ndata_sbs_gemFPP_nclustu_layer);
   fChain->SetBranchAddress("sbs.gemFPP.nclustu_layer", sbs_gemFPP_nclustu_layer, &b_sbs_gemFPP_nclustu_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.nclustv_layer", &Ndata_sbs_gemFPP_nclustv_layer, &b_Ndata_sbs_gemFPP_nclustv_layer);
   fChain->SetBranchAddress("sbs.gemFPP.nclustv_layer", sbs_gemFPP_nclustv_layer, &b_sbs_gemFPP_nclustv_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.nstripsu_layer", &Ndata_sbs_gemFPP_nstripsu_layer, &b_Ndata_sbs_gemFPP_nstripsu_layer);
   fChain->SetBranchAddress("sbs.gemFPP.nstripsu_layer", sbs_gemFPP_nstripsu_layer, &b_sbs_gemFPP_nstripsu_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.nstripsv_layer", &Ndata_sbs_gemFPP_nstripsv_layer, &b_Ndata_sbs_gemFPP_nstripsv_layer);
   fChain->SetBranchAddress("sbs.gemFPP.nstripsv_layer", sbs_gemFPP_nstripsv_layer, &b_sbs_gemFPP_nstripsv_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.chi2ndf", &Ndata_sbs_gemFPP_track_chi2ndf, &b_Ndata_sbs_gemFPP_track_chi2ndf);
   fChain->SetBranchAddress("sbs.gemFPP.track.chi2ndf", sbs_gemFPP_track_chi2ndf, &b_sbs_gemFPP_track_chi2ndf);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.chi2ndf_hitquality", &Ndata_sbs_gemFPP_track_chi2ndf_hitquality, &b_Ndata_sbs_gemFPP_track_chi2ndf_hitquality);
   fChain->SetBranchAddress("sbs.gemFPP.track.chi2ndf_hitquality", sbs_gemFPP_track_chi2ndf_hitquality, &b_sbs_gemFPP_track_chi2ndf_hitquality);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.ngoodhits", &Ndata_sbs_gemFPP_track_ngoodhits, &b_Ndata_sbs_gemFPP_track_ngoodhits);
   fChain->SetBranchAddress("sbs.gemFPP.track.ngoodhits", sbs_gemFPP_track_ngoodhits, &b_sbs_gemFPP_track_ngoodhits);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.nhits", &Ndata_sbs_gemFPP_track_nhits, &b_Ndata_sbs_gemFPP_track_nhits);
   fChain->SetBranchAddress("sbs.gemFPP.track.nhits", sbs_gemFPP_track_nhits, &b_sbs_gemFPP_track_nhits);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.phi", &Ndata_sbs_gemFPP_track_phi, &b_Ndata_sbs_gemFPP_track_phi);
   fChain->SetBranchAddress("sbs.gemFPP.track.phi", sbs_gemFPP_track_phi, &b_sbs_gemFPP_track_phi);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.sclose", &Ndata_sbs_gemFPP_track_sclose, &b_Ndata_sbs_gemFPP_track_sclose);
   fChain->SetBranchAddress("sbs.gemFPP.track.sclose", sbs_gemFPP_track_sclose, &b_sbs_gemFPP_track_sclose);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.t0", &Ndata_sbs_gemFPP_track_t0, &b_Ndata_sbs_gemFPP_track_t0);
   fChain->SetBranchAddress("sbs.gemFPP.track.t0", sbs_gemFPP_track_t0, &b_sbs_gemFPP_track_t0);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.theta", &Ndata_sbs_gemFPP_track_theta, &b_Ndata_sbs_gemFPP_track_theta);
   fChain->SetBranchAddress("sbs.gemFPP.track.theta", sbs_gemFPP_track_theta, &b_sbs_gemFPP_track_theta);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.x", &Ndata_sbs_gemFPP_track_x, &b_Ndata_sbs_gemFPP_track_x);
   fChain->SetBranchAddress("sbs.gemFPP.track.x", sbs_gemFPP_track_x, &b_sbs_gemFPP_track_x);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.xp", &Ndata_sbs_gemFPP_track_xp, &b_Ndata_sbs_gemFPP_track_xp);
   fChain->SetBranchAddress("sbs.gemFPP.track.xp", sbs_gemFPP_track_xp, &b_sbs_gemFPP_track_xp);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.y", &Ndata_sbs_gemFPP_track_y, &b_Ndata_sbs_gemFPP_track_y);
   fChain->SetBranchAddress("sbs.gemFPP.track.y", sbs_gemFPP_track_y, &b_sbs_gemFPP_track_y);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.yp", &Ndata_sbs_gemFPP_track_yp, &b_Ndata_sbs_gemFPP_track_yp);
   fChain->SetBranchAddress("sbs.gemFPP.track.yp", sbs_gemFPP_track_yp, &b_sbs_gemFPP_track_yp);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.zclose", &Ndata_sbs_gemFPP_track_zclose, &b_Ndata_sbs_gemFPP_track_zclose);
   fChain->SetBranchAddress("sbs.gemFPP.track.zclose", sbs_gemFPP_track_zclose, &b_sbs_gemFPP_track_zclose);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCU", &Ndata_sbs_gemFT_hit_ADCU, &b_Ndata_sbs_gemFT_hit_ADCU);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCU", sbs_gemFT_hit_ADCU, &b_sbs_gemFT_hit_ADCU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCU_deconv", &Ndata_sbs_gemFT_hit_ADCU_deconv, &b_Ndata_sbs_gemFT_hit_ADCU_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCU_deconv", sbs_gemFT_hit_ADCU_deconv, &b_sbs_gemFT_hit_ADCU_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCV", &Ndata_sbs_gemFT_hit_ADCV, &b_Ndata_sbs_gemFT_hit_ADCV);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCV", sbs_gemFT_hit_ADCV, &b_sbs_gemFT_hit_ADCV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCV_deconv", &Ndata_sbs_gemFT_hit_ADCV_deconv, &b_Ndata_sbs_gemFT_hit_ADCV_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCV_deconv", sbs_gemFT_hit_ADCV_deconv, &b_sbs_gemFT_hit_ADCV_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCasym", &Ndata_sbs_gemFT_hit_ADCasym, &b_Ndata_sbs_gemFT_hit_ADCasym);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCasym", sbs_gemFT_hit_ADCasym, &b_sbs_gemFT_hit_ADCasym);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCasym_deconv", &Ndata_sbs_gemFT_hit_ADCasym_deconv, &b_Ndata_sbs_gemFT_hit_ADCasym_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCasym_deconv", sbs_gemFT_hit_ADCasym_deconv, &b_sbs_gemFT_hit_ADCasym_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCavg", &Ndata_sbs_gemFT_hit_ADCavg, &b_Ndata_sbs_gemFT_hit_ADCavg);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCavg", sbs_gemFT_hit_ADCavg, &b_sbs_gemFT_hit_ADCavg);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCavg_deconv", &Ndata_sbs_gemFT_hit_ADCavg_deconv, &b_Ndata_sbs_gemFT_hit_ADCavg_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCavg_deconv", sbs_gemFT_hit_ADCavg_deconv, &b_sbs_gemFT_hit_ADCavg_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac0_Umax", &Ndata_sbs_gemFT_hit_ADCfrac0_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac0_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac0_Umax", sbs_gemFT_hit_ADCfrac0_Umax, &b_sbs_gemFT_hit_ADCfrac0_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac0_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac0_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac0_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac0_Vmax", sbs_gemFT_hit_ADCfrac0_Vmax, &b_sbs_gemFT_hit_ADCfrac0_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac1_Umax", &Ndata_sbs_gemFT_hit_ADCfrac1_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac1_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac1_Umax", sbs_gemFT_hit_ADCfrac1_Umax, &b_sbs_gemFT_hit_ADCfrac1_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac1_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac1_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac1_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac1_Vmax", sbs_gemFT_hit_ADCfrac1_Vmax, &b_sbs_gemFT_hit_ADCfrac1_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac2_Umax", &Ndata_sbs_gemFT_hit_ADCfrac2_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac2_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac2_Umax", sbs_gemFT_hit_ADCfrac2_Umax, &b_sbs_gemFT_hit_ADCfrac2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac2_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac2_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac2_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac2_Vmax", sbs_gemFT_hit_ADCfrac2_Vmax, &b_sbs_gemFT_hit_ADCfrac2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac3_Umax", &Ndata_sbs_gemFT_hit_ADCfrac3_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac3_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac3_Umax", sbs_gemFT_hit_ADCfrac3_Umax, &b_sbs_gemFT_hit_ADCfrac3_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac3_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac3_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac3_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac3_Vmax", sbs_gemFT_hit_ADCfrac3_Vmax, &b_sbs_gemFT_hit_ADCfrac3_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac4_Umax", &Ndata_sbs_gemFT_hit_ADCfrac4_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac4_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac4_Umax", sbs_gemFT_hit_ADCfrac4_Umax, &b_sbs_gemFT_hit_ADCfrac4_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac4_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac4_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac4_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac4_Vmax", sbs_gemFT_hit_ADCfrac4_Vmax, &b_sbs_gemFT_hit_ADCfrac4_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac5_Umax", &Ndata_sbs_gemFT_hit_ADCfrac5_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac5_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac5_Umax", sbs_gemFT_hit_ADCfrac5_Umax, &b_sbs_gemFT_hit_ADCfrac5_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac5_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac5_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac5_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac5_Vmax", sbs_gemFT_hit_ADCfrac5_Vmax, &b_sbs_gemFT_hit_ADCfrac5_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxcomboUclust_deconv", &Ndata_sbs_gemFT_hit_ADCmaxcomboUclust_deconv, &b_Ndata_sbs_gemFT_hit_ADCmaxcomboUclust_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxcomboUclust_deconv", sbs_gemFT_hit_ADCmaxcomboUclust_deconv, &b_sbs_gemFT_hit_ADCmaxcomboUclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxcomboVclust_deconv", &Ndata_sbs_gemFT_hit_ADCmaxcomboVclust_deconv, &b_Ndata_sbs_gemFT_hit_ADCmaxcomboVclust_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxcomboVclust_deconv", sbs_gemFT_hit_ADCmaxcomboVclust_deconv, &b_sbs_gemFT_hit_ADCmaxcomboVclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampU", &Ndata_sbs_gemFT_hit_ADCmaxsampU, &b_Ndata_sbs_gemFT_hit_ADCmaxsampU);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampU", sbs_gemFT_hit_ADCmaxsampU, &b_sbs_gemFT_hit_ADCmaxsampU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampUclust", &Ndata_sbs_gemFT_hit_ADCmaxsampUclust, &b_Ndata_sbs_gemFT_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampUclust", sbs_gemFT_hit_ADCmaxsampUclust, &b_sbs_gemFT_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampUclust_deconv", &Ndata_sbs_gemFT_hit_ADCmaxsampUclust_deconv, &b_Ndata_sbs_gemFT_hit_ADCmaxsampUclust_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampUclust_deconv", sbs_gemFT_hit_ADCmaxsampUclust_deconv, &b_sbs_gemFT_hit_ADCmaxsampUclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampV", &Ndata_sbs_gemFT_hit_ADCmaxsampV, &b_Ndata_sbs_gemFT_hit_ADCmaxsampV);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampV", sbs_gemFT_hit_ADCmaxsampV, &b_sbs_gemFT_hit_ADCmaxsampV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampVclust", &Ndata_sbs_gemFT_hit_ADCmaxsampVclust, &b_Ndata_sbs_gemFT_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampVclust", sbs_gemFT_hit_ADCmaxsampVclust, &b_sbs_gemFT_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampVclust_deconv", &Ndata_sbs_gemFT_hit_ADCmaxsampVclust_deconv, &b_Ndata_sbs_gemFT_hit_ADCmaxsampVclust_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampVclust_deconv", sbs_gemFT_hit_ADCmaxsampVclust_deconv, &b_sbs_gemFT_hit_ADCmaxsampVclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxstripU", &Ndata_sbs_gemFT_hit_ADCmaxstripU, &b_Ndata_sbs_gemFT_hit_ADCmaxstripU);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxstripU", sbs_gemFT_hit_ADCmaxstripU, &b_sbs_gemFT_hit_ADCmaxstripU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxstripV", &Ndata_sbs_gemFT_hit_ADCmaxstripV, &b_Ndata_sbs_gemFT_hit_ADCmaxstripV);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxstripV", sbs_gemFT_hit_ADCmaxstripV, &b_sbs_gemFT_hit_ADCmaxstripV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.BUILD_ALL_SAMPLES_U", &Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U, &b_Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.BUILD_ALL_SAMPLES_U", sbs_gemFT_hit_BUILD_ALL_SAMPLES_U, &b_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.BUILD_ALL_SAMPLES_V", &Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V, &b_Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.BUILD_ALL_SAMPLES_V", sbs_gemFT_hit_BUILD_ALL_SAMPLES_V, &b_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.CM_GOOD_U", &Ndata_sbs_gemFT_hit_CM_GOOD_U, &b_Ndata_sbs_gemFT_hit_CM_GOOD_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.CM_GOOD_U", sbs_gemFT_hit_CM_GOOD_U, &b_sbs_gemFT_hit_CM_GOOD_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.CM_GOOD_V", &Ndata_sbs_gemFT_hit_CM_GOOD_V, &b_Ndata_sbs_gemFT_hit_CM_GOOD_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.CM_GOOD_V", sbs_gemFT_hit_CM_GOOD_V, &b_sbs_gemFT_hit_CM_GOOD_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC0_Umax", &Ndata_sbs_gemFT_hit_DeconvADC0_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC0_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC0_Umax", sbs_gemFT_hit_DeconvADC0_Umax, &b_sbs_gemFT_hit_DeconvADC0_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC0_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC0_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC0_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC0_Vmax", sbs_gemFT_hit_DeconvADC0_Vmax, &b_sbs_gemFT_hit_DeconvADC0_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC1_Umax", &Ndata_sbs_gemFT_hit_DeconvADC1_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC1_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC1_Umax", sbs_gemFT_hit_DeconvADC1_Umax, &b_sbs_gemFT_hit_DeconvADC1_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC1_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC1_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC1_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC1_Vmax", sbs_gemFT_hit_DeconvADC1_Vmax, &b_sbs_gemFT_hit_DeconvADC1_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC2_Umax", &Ndata_sbs_gemFT_hit_DeconvADC2_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC2_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC2_Umax", sbs_gemFT_hit_DeconvADC2_Umax, &b_sbs_gemFT_hit_DeconvADC2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC2_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC2_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC2_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC2_Vmax", sbs_gemFT_hit_DeconvADC2_Vmax, &b_sbs_gemFT_hit_DeconvADC2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC3_Umax", &Ndata_sbs_gemFT_hit_DeconvADC3_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC3_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC3_Umax", sbs_gemFT_hit_DeconvADC3_Umax, &b_sbs_gemFT_hit_DeconvADC3_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC3_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC3_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC3_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC3_Vmax", sbs_gemFT_hit_DeconvADC3_Vmax, &b_sbs_gemFT_hit_DeconvADC3_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC4_Umax", &Ndata_sbs_gemFT_hit_DeconvADC4_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC4_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC4_Umax", sbs_gemFT_hit_DeconvADC4_Umax, &b_sbs_gemFT_hit_DeconvADC4_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC4_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC4_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC4_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC4_Vmax", sbs_gemFT_hit_DeconvADC4_Vmax, &b_sbs_gemFT_hit_DeconvADC4_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC5_Umax", &Ndata_sbs_gemFT_hit_DeconvADC5_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC5_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC5_Umax", sbs_gemFT_hit_DeconvADC5_Umax, &b_sbs_gemFT_hit_DeconvADC5_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC5_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC5_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC5_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC5_Vmax", sbs_gemFT_hit_DeconvADC5_Vmax, &b_sbs_gemFT_hit_DeconvADC5_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxcomboU", &Ndata_sbs_gemFT_hit_DeconvADCmaxcomboU, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxcomboU);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxcomboU", sbs_gemFT_hit_DeconvADCmaxcomboU, &b_sbs_gemFT_hit_DeconvADCmaxcomboU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxcomboV", &Ndata_sbs_gemFT_hit_DeconvADCmaxcomboV, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxcomboV);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxcomboV", sbs_gemFT_hit_DeconvADCmaxcomboV, &b_sbs_gemFT_hit_DeconvADCmaxcomboV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxsampU", &Ndata_sbs_gemFT_hit_DeconvADCmaxsampU, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxsampU);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxsampU", sbs_gemFT_hit_DeconvADCmaxsampU, &b_sbs_gemFT_hit_DeconvADCmaxsampU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxsampV", &Ndata_sbs_gemFT_hit_DeconvADCmaxsampV, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxsampV);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxsampV", sbs_gemFT_hit_DeconvADCmaxsampV, &b_sbs_gemFT_hit_DeconvADCmaxsampV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxstripU", &Ndata_sbs_gemFT_hit_DeconvADCmaxstripU, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxstripU);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxstripU", sbs_gemFT_hit_DeconvADCmaxstripU, &b_sbs_gemFT_hit_DeconvADCmaxstripU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxstripV", &Ndata_sbs_gemFT_hit_DeconvADCmaxstripV, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxstripV);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxstripV", sbs_gemFT_hit_DeconvADCmaxstripV, &b_sbs_gemFT_hit_DeconvADCmaxstripV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ENABLE_CM_U", &Ndata_sbs_gemFT_hit_ENABLE_CM_U, &b_Ndata_sbs_gemFT_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.ENABLE_CM_U", sbs_gemFT_hit_ENABLE_CM_U, &b_sbs_gemFT_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ENABLE_CM_V", &Ndata_sbs_gemFT_hit_ENABLE_CM_V, &b_Ndata_sbs_gemFT_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.ENABLE_CM_V", sbs_gemFT_hit_ENABLE_CM_V, &b_sbs_gemFT_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.TSchi2_Umax", &Ndata_sbs_gemFT_hit_TSchi2_Umax, &b_Ndata_sbs_gemFT_hit_TSchi2_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.TSchi2_Umax", sbs_gemFT_hit_TSchi2_Umax, &b_sbs_gemFT_hit_TSchi2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.TSchi2_Vmax", &Ndata_sbs_gemFT_hit_TSchi2_Vmax, &b_Ndata_sbs_gemFT_hit_TSchi2_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.TSchi2_Vmax", sbs_gemFT_hit_TSchi2_Vmax, &b_sbs_gemFT_hit_TSchi2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.TSprob_Umax", &Ndata_sbs_gemFT_hit_TSprob_Umax, &b_Ndata_sbs_gemFT_hit_TSprob_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.TSprob_Umax", sbs_gemFT_hit_TSprob_Umax, &b_sbs_gemFT_hit_TSprob_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.TSprob_Vmax", &Ndata_sbs_gemFT_hit_TSprob_Vmax, &b_Ndata_sbs_gemFT_hit_TSprob_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.TSprob_Vmax", sbs_gemFT_hit_TSprob_Vmax, &b_sbs_gemFT_hit_TSprob_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Tavg", &Ndata_sbs_gemFT_hit_Tavg, &b_Ndata_sbs_gemFT_hit_Tavg);
   fChain->SetBranchAddress("sbs.gemFT.hit.Tavg", sbs_gemFT_hit_Tavg, &b_sbs_gemFT_hit_Tavg);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Tavg_corr", &Ndata_sbs_gemFT_hit_Tavg_corr, &b_Ndata_sbs_gemFT_hit_Tavg_corr);
   fChain->SetBranchAddress("sbs.gemFT.hit.Tavg_corr", sbs_gemFT_hit_Tavg_corr, &b_sbs_gemFT_hit_Tavg_corr);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Tavg_deconv", &Ndata_sbs_gemFT_hit_Tavg_deconv, &b_Ndata_sbs_gemFT_hit_Tavg_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.Tavg_deconv", sbs_gemFT_hit_Tavg_deconv, &b_sbs_gemFT_hit_Tavg_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Tavg_fit", &Ndata_sbs_gemFT_hit_Tavg_fit, &b_Ndata_sbs_gemFT_hit_Tavg_fit);
   fChain->SetBranchAddress("sbs.gemFT.hit.Tavg_fit", sbs_gemFT_hit_Tavg_fit, &b_sbs_gemFT_hit_Tavg_fit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Ugain", &Ndata_sbs_gemFT_hit_Ugain, &b_Ndata_sbs_gemFT_hit_Ugain);
   fChain->SetBranchAddress("sbs.gemFT.hit.Ugain", sbs_gemFT_hit_Ugain, &b_sbs_gemFT_hit_Ugain);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Utime", &Ndata_sbs_gemFT_hit_Utime, &b_Ndata_sbs_gemFT_hit_Utime);
   fChain->SetBranchAddress("sbs.gemFT.hit.Utime", sbs_gemFT_hit_Utime, &b_sbs_gemFT_hit_Utime);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.UtimeDeconv", &Ndata_sbs_gemFT_hit_UtimeDeconv, &b_Ndata_sbs_gemFT_hit_UtimeDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.UtimeDeconv", sbs_gemFT_hit_UtimeDeconv, &b_sbs_gemFT_hit_UtimeDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.UtimeFit", &Ndata_sbs_gemFT_hit_UtimeFit, &b_Ndata_sbs_gemFT_hit_UtimeFit);
   fChain->SetBranchAddress("sbs.gemFT.hit.UtimeFit", sbs_gemFT_hit_UtimeFit, &b_sbs_gemFT_hit_UtimeFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.UtimeMaxStrip", &Ndata_sbs_gemFT_hit_UtimeMaxStrip, &b_Ndata_sbs_gemFT_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("sbs.gemFT.hit.UtimeMaxStrip", sbs_gemFT_hit_UtimeMaxStrip, &b_sbs_gemFT_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.UtimeMaxStripDeconv", &Ndata_sbs_gemFT_hit_UtimeMaxStripDeconv, &b_Ndata_sbs_gemFT_hit_UtimeMaxStripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.UtimeMaxStripDeconv", sbs_gemFT_hit_UtimeMaxStripDeconv, &b_sbs_gemFT_hit_UtimeMaxStripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.UtimeMaxStripFit", &Ndata_sbs_gemFT_hit_UtimeMaxStripFit, &b_Ndata_sbs_gemFT_hit_UtimeMaxStripFit);
   fChain->SetBranchAddress("sbs.gemFT.hit.UtimeMaxStripFit", sbs_gemFT_hit_UtimeMaxStripFit, &b_sbs_gemFT_hit_UtimeMaxStripFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Vgain", &Ndata_sbs_gemFT_hit_Vgain, &b_Ndata_sbs_gemFT_hit_Vgain);
   fChain->SetBranchAddress("sbs.gemFT.hit.Vgain", sbs_gemFT_hit_Vgain, &b_sbs_gemFT_hit_Vgain);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Vtime", &Ndata_sbs_gemFT_hit_Vtime, &b_Ndata_sbs_gemFT_hit_Vtime);
   fChain->SetBranchAddress("sbs.gemFT.hit.Vtime", sbs_gemFT_hit_Vtime, &b_sbs_gemFT_hit_Vtime);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.VtimeDeconv", &Ndata_sbs_gemFT_hit_VtimeDeconv, &b_Ndata_sbs_gemFT_hit_VtimeDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.VtimeDeconv", sbs_gemFT_hit_VtimeDeconv, &b_sbs_gemFT_hit_VtimeDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.VtimeFit", &Ndata_sbs_gemFT_hit_VtimeFit, &b_Ndata_sbs_gemFT_hit_VtimeFit);
   fChain->SetBranchAddress("sbs.gemFT.hit.VtimeFit", sbs_gemFT_hit_VtimeFit, &b_sbs_gemFT_hit_VtimeFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.VtimeMaxStrip", &Ndata_sbs_gemFT_hit_VtimeMaxStrip, &b_Ndata_sbs_gemFT_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("sbs.gemFT.hit.VtimeMaxStrip", sbs_gemFT_hit_VtimeMaxStrip, &b_sbs_gemFT_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.VtimeMaxStripDeconv", &Ndata_sbs_gemFT_hit_VtimeMaxStripDeconv, &b_Ndata_sbs_gemFT_hit_VtimeMaxStripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.VtimeMaxStripDeconv", sbs_gemFT_hit_VtimeMaxStripDeconv, &b_sbs_gemFT_hit_VtimeMaxStripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.VtimeMaxStripFit", &Ndata_sbs_gemFT_hit_VtimeMaxStripFit, &b_Ndata_sbs_gemFT_hit_VtimeMaxStripFit);
   fChain->SetBranchAddress("sbs.gemFT.hit.VtimeMaxStripFit", sbs_gemFT_hit_VtimeMaxStripFit, &b_sbs_gemFT_hit_VtimeMaxStripFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.adc_id_U", &Ndata_sbs_gemFT_hit_adc_id_U, &b_Ndata_sbs_gemFT_hit_adc_id_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.adc_id_U", sbs_gemFT_hit_adc_id_U, &b_sbs_gemFT_hit_adc_id_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.adc_id_V", &Ndata_sbs_gemFT_hit_adc_id_V, &b_Ndata_sbs_gemFT_hit_adc_id_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.adc_id_V", sbs_gemFT_hit_adc_id_V, &b_sbs_gemFT_hit_adc_id_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ccor_clust", &Ndata_sbs_gemFT_hit_ccor_clust, &b_Ndata_sbs_gemFT_hit_ccor_clust);
   fChain->SetBranchAddress("sbs.gemFT.hit.ccor_clust", sbs_gemFT_hit_ccor_clust, &b_sbs_gemFT_hit_ccor_clust);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ccor_clust_deconv", &Ndata_sbs_gemFT_hit_ccor_clust_deconv, &b_Ndata_sbs_gemFT_hit_ccor_clust_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ccor_clust_deconv", sbs_gemFT_hit_ccor_clust_deconv, &b_sbs_gemFT_hit_ccor_clust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ccor_strip", &Ndata_sbs_gemFT_hit_ccor_strip, &b_Ndata_sbs_gemFT_hit_ccor_strip);
   fChain->SetBranchAddress("sbs.gemFT.hit.ccor_strip", sbs_gemFT_hit_ccor_strip, &b_sbs_gemFT_hit_ccor_strip);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ccor_strip_deconv", &Ndata_sbs_gemFT_hit_ccor_strip_deconv, &b_Ndata_sbs_gemFT_hit_ccor_strip_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ccor_strip_deconv", sbs_gemFT_hit_ccor_strip_deconv, &b_sbs_gemFT_hit_ccor_strip_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.crate_U", &Ndata_sbs_gemFT_hit_crate_U, &b_Ndata_sbs_gemFT_hit_crate_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.crate_U", sbs_gemFT_hit_crate_U, &b_sbs_gemFT_hit_crate_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.crate_V", &Ndata_sbs_gemFT_hit_crate_V, &b_Ndata_sbs_gemFT_hit_crate_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.crate_V", sbs_gemFT_hit_crate_V, &b_sbs_gemFT_hit_crate_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.deltat", &Ndata_sbs_gemFT_hit_deltat, &b_Ndata_sbs_gemFT_hit_deltat);
   fChain->SetBranchAddress("sbs.gemFT.hit.deltat", sbs_gemFT_hit_deltat, &b_sbs_gemFT_hit_deltat);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.deltat_deconv", &Ndata_sbs_gemFT_hit_deltat_deconv, &b_Ndata_sbs_gemFT_hit_deltat_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.deltat_deconv", sbs_gemFT_hit_deltat_deconv, &b_sbs_gemFT_hit_deltat_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.deltat_fit", &Ndata_sbs_gemFT_hit_deltat_fit, &b_Ndata_sbs_gemFT_hit_deltat_fit);
   fChain->SetBranchAddress("sbs.gemFT.hit.deltat_fit", sbs_gemFT_hit_deltat_fit, &b_sbs_gemFT_hit_deltat_fit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.eresidu", &Ndata_sbs_gemFT_hit_eresidu, &b_Ndata_sbs_gemFT_hit_eresidu);
   fChain->SetBranchAddress("sbs.gemFT.hit.eresidu", sbs_gemFT_hit_eresidu, &b_sbs_gemFT_hit_eresidu);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.eresidv", &Ndata_sbs_gemFT_hit_eresidv, &b_Ndata_sbs_gemFT_hit_eresidv);
   fChain->SetBranchAddress("sbs.gemFT.hit.eresidv", sbs_gemFT_hit_eresidv, &b_sbs_gemFT_hit_eresidv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.icombomaxUclustDeconv", &Ndata_sbs_gemFT_hit_icombomaxUclustDeconv, &b_Ndata_sbs_gemFT_hit_icombomaxUclustDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.icombomaxUclustDeconv", sbs_gemFT_hit_icombomaxUclustDeconv, &b_sbs_gemFT_hit_icombomaxUclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.icombomaxUstripDeconv", &Ndata_sbs_gemFT_hit_icombomaxUstripDeconv, &b_Ndata_sbs_gemFT_hit_icombomaxUstripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.icombomaxUstripDeconv", sbs_gemFT_hit_icombomaxUstripDeconv, &b_sbs_gemFT_hit_icombomaxUstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.icombomaxVclustDeconv", &Ndata_sbs_gemFT_hit_icombomaxVclustDeconv, &b_Ndata_sbs_gemFT_hit_icombomaxVclustDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.icombomaxVclustDeconv", sbs_gemFT_hit_icombomaxVclustDeconv, &b_sbs_gemFT_hit_icombomaxVclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.icombomaxVstripDeconv", &Ndata_sbs_gemFT_hit_icombomaxVstripDeconv, &b_Ndata_sbs_gemFT_hit_icombomaxVstripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.icombomaxVstripDeconv", sbs_gemFT_hit_icombomaxVstripDeconv, &b_sbs_gemFT_hit_icombomaxVstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxUclust", &Ndata_sbs_gemFT_hit_isampmaxUclust, &b_Ndata_sbs_gemFT_hit_isampmaxUclust);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxUclust", sbs_gemFT_hit_isampmaxUclust, &b_sbs_gemFT_hit_isampmaxUclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxUclustDeconv", &Ndata_sbs_gemFT_hit_isampmaxUclustDeconv, &b_Ndata_sbs_gemFT_hit_isampmaxUclustDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxUclustDeconv", sbs_gemFT_hit_isampmaxUclustDeconv, &b_sbs_gemFT_hit_isampmaxUclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxUstrip", &Ndata_sbs_gemFT_hit_isampmaxUstrip, &b_Ndata_sbs_gemFT_hit_isampmaxUstrip);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxUstrip", sbs_gemFT_hit_isampmaxUstrip, &b_sbs_gemFT_hit_isampmaxUstrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxUstripDeconv", &Ndata_sbs_gemFT_hit_isampmaxUstripDeconv, &b_Ndata_sbs_gemFT_hit_isampmaxUstripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxUstripDeconv", sbs_gemFT_hit_isampmaxUstripDeconv, &b_sbs_gemFT_hit_isampmaxUstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxVclust", &Ndata_sbs_gemFT_hit_isampmaxVclust, &b_Ndata_sbs_gemFT_hit_isampmaxVclust);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxVclust", sbs_gemFT_hit_isampmaxVclust, &b_sbs_gemFT_hit_isampmaxVclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxVclustDeconv", &Ndata_sbs_gemFT_hit_isampmaxVclustDeconv, &b_Ndata_sbs_gemFT_hit_isampmaxVclustDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxVclustDeconv", sbs_gemFT_hit_isampmaxVclustDeconv, &b_sbs_gemFT_hit_isampmaxVclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxVstrip", &Ndata_sbs_gemFT_hit_isampmaxVstrip, &b_Ndata_sbs_gemFT_hit_isampmaxVstrip);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxVstrip", sbs_gemFT_hit_isampmaxVstrip, &b_sbs_gemFT_hit_isampmaxVstrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxVstripDeconv", &Ndata_sbs_gemFT_hit_isampmaxVstripDeconv, &b_Ndata_sbs_gemFT_hit_isampmaxVstripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxVstripDeconv", sbs_gemFT_hit_isampmaxVstripDeconv, &b_sbs_gemFT_hit_isampmaxVstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.layer", &Ndata_sbs_gemFT_hit_layer, &b_Ndata_sbs_gemFT_hit_layer);
   fChain->SetBranchAddress("sbs.gemFT.hit.layer", sbs_gemFT_hit_layer, &b_sbs_gemFT_hit_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.module", &Ndata_sbs_gemFT_hit_module, &b_Ndata_sbs_gemFT_hit_module);
   fChain->SetBranchAddress("sbs.gemFT.hit.module", sbs_gemFT_hit_module, &b_sbs_gemFT_hit_module);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.mpd_U", &Ndata_sbs_gemFT_hit_mpd_U, &b_Ndata_sbs_gemFT_hit_mpd_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.mpd_U", sbs_gemFT_hit_mpd_U, &b_sbs_gemFT_hit_mpd_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.mpd_V", &Ndata_sbs_gemFT_hit_mpd_V, &b_Ndata_sbs_gemFT_hit_mpd_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.mpd_V", sbs_gemFT_hit_mpd_V, &b_sbs_gemFT_hit_mpd_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.nstripu", &Ndata_sbs_gemFT_hit_nstripu, &b_Ndata_sbs_gemFT_hit_nstripu);
   fChain->SetBranchAddress("sbs.gemFT.hit.nstripu", sbs_gemFT_hit_nstripu, &b_sbs_gemFT_hit_nstripu);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.nstripv", &Ndata_sbs_gemFT_hit_nstripv, &b_Ndata_sbs_gemFT_hit_nstripv);
   fChain->SetBranchAddress("sbs.gemFT.hit.nstripv", sbs_gemFT_hit_nstripv, &b_sbs_gemFT_hit_nstripv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.residu", &Ndata_sbs_gemFT_hit_residu, &b_Ndata_sbs_gemFT_hit_residu);
   fChain->SetBranchAddress("sbs.gemFT.hit.residu", sbs_gemFT_hit_residu, &b_sbs_gemFT_hit_residu);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.residv", &Ndata_sbs_gemFT_hit_residv, &b_Ndata_sbs_gemFT_hit_residv);
   fChain->SetBranchAddress("sbs.gemFT.hit.residv", sbs_gemFT_hit_residv, &b_sbs_gemFT_hit_residv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.trackindex", &Ndata_sbs_gemFT_hit_trackindex, &b_Ndata_sbs_gemFT_hit_trackindex);
   fChain->SetBranchAddress("sbs.gemFT.hit.trackindex", sbs_gemFT_hit_trackindex, &b_sbs_gemFT_hit_trackindex);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.u", &Ndata_sbs_gemFT_hit_u, &b_Ndata_sbs_gemFT_hit_u);
   fChain->SetBranchAddress("sbs.gemFT.hit.u", sbs_gemFT_hit_u, &b_sbs_gemFT_hit_u);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.umoment", &Ndata_sbs_gemFT_hit_umoment, &b_Ndata_sbs_gemFT_hit_umoment);
   fChain->SetBranchAddress("sbs.gemFT.hit.umoment", sbs_gemFT_hit_umoment, &b_sbs_gemFT_hit_umoment);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.usigma", &Ndata_sbs_gemFT_hit_usigma, &b_Ndata_sbs_gemFT_hit_usigma);
   fChain->SetBranchAddress("sbs.gemFT.hit.usigma", sbs_gemFT_hit_usigma, &b_sbs_gemFT_hit_usigma);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ustriphi", &Ndata_sbs_gemFT_hit_ustriphi, &b_Ndata_sbs_gemFT_hit_ustriphi);
   fChain->SetBranchAddress("sbs.gemFT.hit.ustriphi", sbs_gemFT_hit_ustriphi, &b_sbs_gemFT_hit_ustriphi);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ustriplo", &Ndata_sbs_gemFT_hit_ustriplo, &b_Ndata_sbs_gemFT_hit_ustriplo);
   fChain->SetBranchAddress("sbs.gemFT.hit.ustriplo", sbs_gemFT_hit_ustriplo, &b_sbs_gemFT_hit_ustriplo);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ustripmax", &Ndata_sbs_gemFT_hit_ustripmax, &b_Ndata_sbs_gemFT_hit_ustripmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ustripmax", sbs_gemFT_hit_ustripmax, &b_sbs_gemFT_hit_ustripmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.v", &Ndata_sbs_gemFT_hit_v, &b_Ndata_sbs_gemFT_hit_v);
   fChain->SetBranchAddress("sbs.gemFT.hit.v", sbs_gemFT_hit_v, &b_sbs_gemFT_hit_v);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.vmoment", &Ndata_sbs_gemFT_hit_vmoment, &b_Ndata_sbs_gemFT_hit_vmoment);
   fChain->SetBranchAddress("sbs.gemFT.hit.vmoment", sbs_gemFT_hit_vmoment, &b_sbs_gemFT_hit_vmoment);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.vsigma", &Ndata_sbs_gemFT_hit_vsigma, &b_Ndata_sbs_gemFT_hit_vsigma);
   fChain->SetBranchAddress("sbs.gemFT.hit.vsigma", sbs_gemFT_hit_vsigma, &b_sbs_gemFT_hit_vsigma);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.vstriphi", &Ndata_sbs_gemFT_hit_vstriphi, &b_Ndata_sbs_gemFT_hit_vstriphi);
   fChain->SetBranchAddress("sbs.gemFT.hit.vstriphi", sbs_gemFT_hit_vstriphi, &b_sbs_gemFT_hit_vstriphi);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.vstriplo", &Ndata_sbs_gemFT_hit_vstriplo, &b_Ndata_sbs_gemFT_hit_vstriplo);
   fChain->SetBranchAddress("sbs.gemFT.hit.vstriplo", sbs_gemFT_hit_vstriplo, &b_sbs_gemFT_hit_vstriplo);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.vstripmax", &Ndata_sbs_gemFT_hit_vstripmax, &b_Ndata_sbs_gemFT_hit_vstripmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.vstripmax", sbs_gemFT_hit_vstripmax, &b_sbs_gemFT_hit_vstripmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.xglobal", &Ndata_sbs_gemFT_hit_xglobal, &b_Ndata_sbs_gemFT_hit_xglobal);
   fChain->SetBranchAddress("sbs.gemFT.hit.xglobal", sbs_gemFT_hit_xglobal, &b_sbs_gemFT_hit_xglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.xlocal", &Ndata_sbs_gemFT_hit_xlocal, &b_Ndata_sbs_gemFT_hit_xlocal);
   fChain->SetBranchAddress("sbs.gemFT.hit.xlocal", sbs_gemFT_hit_xlocal, &b_sbs_gemFT_hit_xlocal);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.yglobal", &Ndata_sbs_gemFT_hit_yglobal, &b_Ndata_sbs_gemFT_hit_yglobal);
   fChain->SetBranchAddress("sbs.gemFT.hit.yglobal", sbs_gemFT_hit_yglobal, &b_sbs_gemFT_hit_yglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ylocal", &Ndata_sbs_gemFT_hit_ylocal, &b_Ndata_sbs_gemFT_hit_ylocal);
   fChain->SetBranchAddress("sbs.gemFT.hit.ylocal", sbs_gemFT_hit_ylocal, &b_sbs_gemFT_hit_ylocal);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.zglobal", &Ndata_sbs_gemFT_hit_zglobal, &b_Ndata_sbs_gemFT_hit_zglobal);
   fChain->SetBranchAddress("sbs.gemFT.hit.zglobal", sbs_gemFT_hit_zglobal, &b_sbs_gemFT_hit_zglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.n2Dhit_layer", &Ndata_sbs_gemFT_n2Dhit_layer, &b_Ndata_sbs_gemFT_n2Dhit_layer);
   fChain->SetBranchAddress("sbs.gemFT.n2Dhit_layer", sbs_gemFT_n2Dhit_layer, &b_sbs_gemFT_n2Dhit_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.nclustu_layer", &Ndata_sbs_gemFT_nclustu_layer, &b_Ndata_sbs_gemFT_nclustu_layer);
   fChain->SetBranchAddress("sbs.gemFT.nclustu_layer", sbs_gemFT_nclustu_layer, &b_sbs_gemFT_nclustu_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.nclustv_layer", &Ndata_sbs_gemFT_nclustv_layer, &b_Ndata_sbs_gemFT_nclustv_layer);
   fChain->SetBranchAddress("sbs.gemFT.nclustv_layer", sbs_gemFT_nclustv_layer, &b_sbs_gemFT_nclustv_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.nstripsu_layer", &Ndata_sbs_gemFT_nstripsu_layer, &b_Ndata_sbs_gemFT_nstripsu_layer);
   fChain->SetBranchAddress("sbs.gemFT.nstripsu_layer", sbs_gemFT_nstripsu_layer, &b_sbs_gemFT_nstripsu_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.nstripsv_layer", &Ndata_sbs_gemFT_nstripsv_layer, &b_Ndata_sbs_gemFT_nstripsv_layer);
   fChain->SetBranchAddress("sbs.gemFT.nstripsv_layer", sbs_gemFT_nstripsv_layer, &b_sbs_gemFT_nstripsv_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.chi2ndf", &Ndata_sbs_gemFT_track_chi2ndf, &b_Ndata_sbs_gemFT_track_chi2ndf);
   fChain->SetBranchAddress("sbs.gemFT.track.chi2ndf", sbs_gemFT_track_chi2ndf, &b_sbs_gemFT_track_chi2ndf);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.chi2ndf_hitquality", &Ndata_sbs_gemFT_track_chi2ndf_hitquality, &b_Ndata_sbs_gemFT_track_chi2ndf_hitquality);
   fChain->SetBranchAddress("sbs.gemFT.track.chi2ndf_hitquality", sbs_gemFT_track_chi2ndf_hitquality, &b_sbs_gemFT_track_chi2ndf_hitquality);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.ngoodhits", &Ndata_sbs_gemFT_track_ngoodhits, &b_Ndata_sbs_gemFT_track_ngoodhits);
   fChain->SetBranchAddress("sbs.gemFT.track.ngoodhits", sbs_gemFT_track_ngoodhits, &b_sbs_gemFT_track_ngoodhits);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.nhits", &Ndata_sbs_gemFT_track_nhits, &b_Ndata_sbs_gemFT_track_nhits);
   fChain->SetBranchAddress("sbs.gemFT.track.nhits", sbs_gemFT_track_nhits, &b_sbs_gemFT_track_nhits);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.t0", &Ndata_sbs_gemFT_track_t0, &b_Ndata_sbs_gemFT_track_t0);
   fChain->SetBranchAddress("sbs.gemFT.track.t0", sbs_gemFT_track_t0, &b_sbs_gemFT_track_t0);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.x", &Ndata_sbs_gemFT_track_x, &b_Ndata_sbs_gemFT_track_x);
   fChain->SetBranchAddress("sbs.gemFT.track.x", sbs_gemFT_track_x, &b_sbs_gemFT_track_x);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.xp", &Ndata_sbs_gemFT_track_xp, &b_Ndata_sbs_gemFT_track_xp);
   fChain->SetBranchAddress("sbs.gemFT.track.xp", sbs_gemFT_track_xp, &b_sbs_gemFT_track_xp);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.y", &Ndata_sbs_gemFT_track_y, &b_Ndata_sbs_gemFT_track_y);
   fChain->SetBranchAddress("sbs.gemFT.track.y", sbs_gemFT_track_y, &b_sbs_gemFT_track_y);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.yp", &Ndata_sbs_gemFT_track_yp, &b_Ndata_sbs_gemFT_track_yp);
   fChain->SetBranchAddress("sbs.gemFT.track.yp", sbs_gemFT_track_yp, &b_sbs_gemFT_track_yp);
   fChain->SetBranchAddress("Ndata.sbs.hcal.Ref.tdc", &Ndata_sbs_hcal_Ref_tdc, &b_Ndata_sbs_hcal_Ref_tdc);
   fChain->SetBranchAddress("sbs.hcal.Ref.tdc", sbs_hcal_Ref_tdc, &b_sbs_hcal_Ref_tdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.Ref.tdc_mult", &Ndata_sbs_hcal_Ref_tdc_mult, &b_Ndata_sbs_hcal_Ref_tdc_mult);
   fChain->SetBranchAddress("sbs.hcal.Ref.tdc_mult", sbs_hcal_Ref_tdc_mult, &b_sbs_hcal_Ref_tdc_mult);
   fChain->SetBranchAddress("Ndata.sbs.hcal.Ref.tdcelemID", &Ndata_sbs_hcal_Ref_tdcelemID, &b_Ndata_sbs_hcal_Ref_tdcelemID);
   fChain->SetBranchAddress("sbs.hcal.Ref.tdcelemID", sbs_hcal_Ref_tdcelemID, &b_sbs_hcal_Ref_tdcelemID);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.adctime", &Ndata_sbs_hcal_clus_adctime, &b_Ndata_sbs_hcal_clus_adctime);
   fChain->SetBranchAddress("sbs.hcal.clus.adctime", sbs_hcal_clus_adctime, &b_sbs_hcal_clus_adctime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.again", &Ndata_sbs_hcal_clus_again, &b_Ndata_sbs_hcal_clus_again);
   fChain->SetBranchAddress("sbs.hcal.clus.again", sbs_hcal_clus_again, &b_sbs_hcal_clus_again);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.atimeblk", &Ndata_sbs_hcal_clus_atimeblk, &b_Ndata_sbs_hcal_clus_atimeblk);
   fChain->SetBranchAddress("sbs.hcal.clus.atimeblk", sbs_hcal_clus_atimeblk, &b_sbs_hcal_clus_atimeblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.col", &Ndata_sbs_hcal_clus_col, &b_Ndata_sbs_hcal_clus_col);
   fChain->SetBranchAddress("sbs.hcal.clus.col", sbs_hcal_clus_col, &b_sbs_hcal_clus_col);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.col_goodtdc", &Ndata_sbs_hcal_clus_col_goodtdc, &b_Ndata_sbs_hcal_clus_col_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.col_goodtdc", sbs_hcal_clus_col_goodtdc, &b_sbs_hcal_clus_col_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.e", &Ndata_sbs_hcal_clus_e, &b_Ndata_sbs_hcal_clus_e);
   fChain->SetBranchAddress("sbs.hcal.clus.e", sbs_hcal_clus_e, &b_sbs_hcal_clus_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.e_goodtdc", &Ndata_sbs_hcal_clus_e_goodtdc, &b_Ndata_sbs_hcal_clus_e_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.e_goodtdc", sbs_hcal_clus_e_goodtdc, &b_sbs_hcal_clus_e_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.eblk", &Ndata_sbs_hcal_clus_eblk, &b_Ndata_sbs_hcal_clus_eblk);
   fChain->SetBranchAddress("sbs.hcal.clus.eblk", sbs_hcal_clus_eblk, &b_sbs_hcal_clus_eblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.eblk_goodtdc", &Ndata_sbs_hcal_clus_eblk_goodtdc, &b_Ndata_sbs_hcal_clus_eblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.eblk_goodtdc", sbs_hcal_clus_eblk_goodtdc, &b_sbs_hcal_clus_eblk_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.id", &Ndata_sbs_hcal_clus_id, &b_Ndata_sbs_hcal_clus_id);
   fChain->SetBranchAddress("sbs.hcal.clus.id", sbs_hcal_clus_id, &b_sbs_hcal_clus_id);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.id_goodtdc", &Ndata_sbs_hcal_clus_id_goodtdc, &b_Ndata_sbs_hcal_clus_id_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.id_goodtdc", sbs_hcal_clus_id_goodtdc, &b_sbs_hcal_clus_id_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.nblk", &Ndata_sbs_hcal_clus_nblk, &b_Ndata_sbs_hcal_clus_nblk);
   fChain->SetBranchAddress("sbs.hcal.clus.nblk", sbs_hcal_clus_nblk, &b_sbs_hcal_clus_nblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.nblk_goodtdc", &Ndata_sbs_hcal_clus_nblk_goodtdc, &b_Ndata_sbs_hcal_clus_nblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.nblk_goodtdc", sbs_hcal_clus_nblk_goodtdc, &b_sbs_hcal_clus_nblk_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.row", &Ndata_sbs_hcal_clus_row, &b_Ndata_sbs_hcal_clus_row);
   fChain->SetBranchAddress("sbs.hcal.clus.row", sbs_hcal_clus_row, &b_sbs_hcal_clus_row);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.row_goodtdc", &Ndata_sbs_hcal_clus_row_goodtdc, &b_Ndata_sbs_hcal_clus_row_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.row_goodtdc", sbs_hcal_clus_row_goodtdc, &b_sbs_hcal_clus_row_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.tdctime", &Ndata_sbs_hcal_clus_tdctime, &b_Ndata_sbs_hcal_clus_tdctime);
   fChain->SetBranchAddress("sbs.hcal.clus.tdctime", sbs_hcal_clus_tdctime, &b_sbs_hcal_clus_tdctime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.tdctime_tw", &Ndata_sbs_hcal_clus_tdctime_tw, &b_Ndata_sbs_hcal_clus_tdctime_tw);
   fChain->SetBranchAddress("sbs.hcal.clus.tdctime_tw", sbs_hcal_clus_tdctime_tw, &b_sbs_hcal_clus_tdctime_tw);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.tdctimeblk", &Ndata_sbs_hcal_clus_tdctimeblk, &b_Ndata_sbs_hcal_clus_tdctimeblk);
   fChain->SetBranchAddress("sbs.hcal.clus.tdctimeblk", sbs_hcal_clus_tdctimeblk, &b_sbs_hcal_clus_tdctimeblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.tdctimeblk_tw", &Ndata_sbs_hcal_clus_tdctimeblk_tw, &b_Ndata_sbs_hcal_clus_tdctimeblk_tw);
   fChain->SetBranchAddress("sbs.hcal.clus.tdctimeblk_tw", sbs_hcal_clus_tdctimeblk_tw, &b_sbs_hcal_clus_tdctimeblk_tw);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.x", &Ndata_sbs_hcal_clus_x, &b_Ndata_sbs_hcal_clus_x);
   fChain->SetBranchAddress("sbs.hcal.clus.x", sbs_hcal_clus_x, &b_sbs_hcal_clus_x);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.y", &Ndata_sbs_hcal_clus_y, &b_Ndata_sbs_hcal_clus_y);
   fChain->SetBranchAddress("sbs.hcal.clus.y", sbs_hcal_clus_y, &b_sbs_hcal_clus_y);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.again", &Ndata_sbs_hcal_clus_blk_again, &b_Ndata_sbs_hcal_clus_blk_again);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.again", sbs_hcal_clus_blk_again, &b_sbs_hcal_clus_blk_again);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.atime", &Ndata_sbs_hcal_clus_blk_atime, &b_Ndata_sbs_hcal_clus_blk_atime);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.atime", sbs_hcal_clus_blk_atime, &b_sbs_hcal_clus_blk_atime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.col", &Ndata_sbs_hcal_clus_blk_col, &b_Ndata_sbs_hcal_clus_blk_col);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.col", sbs_hcal_clus_blk_col, &b_sbs_hcal_clus_blk_col);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.e", &Ndata_sbs_hcal_clus_blk_e, &b_Ndata_sbs_hcal_clus_blk_e);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.e", sbs_hcal_clus_blk_e, &b_sbs_hcal_clus_blk_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.id", &Ndata_sbs_hcal_clus_blk_id, &b_Ndata_sbs_hcal_clus_blk_id);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.id", sbs_hcal_clus_blk_id, &b_sbs_hcal_clus_blk_id);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.row", &Ndata_sbs_hcal_clus_blk_row, &b_Ndata_sbs_hcal_clus_blk_row);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.row", sbs_hcal_clus_blk_row, &b_sbs_hcal_clus_blk_row);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.tdctime", &Ndata_sbs_hcal_clus_blk_tdctime, &b_Ndata_sbs_hcal_clus_blk_tdctime);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.tdctime", sbs_hcal_clus_blk_tdctime, &b_sbs_hcal_clus_blk_tdctime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.tdctime_tw", &Ndata_sbs_hcal_clus_blk_tdctime_tw, &b_Ndata_sbs_hcal_clus_blk_tdctime_tw);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.tdctime_tw", sbs_hcal_clus_blk_tdctime_tw, &b_sbs_hcal_clus_blk_tdctime_tw);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.x", &Ndata_sbs_hcal_clus_blk_x, &b_Ndata_sbs_hcal_clus_blk_x);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.x", sbs_hcal_clus_blk_x, &b_sbs_hcal_clus_blk_x);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.y", &Ndata_sbs_hcal_clus_blk_y, &b_Ndata_sbs_hcal_clus_blk_y);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.y", sbs_hcal_clus_blk_y, &b_sbs_hcal_clus_blk_y);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.atime", &Ndata_sbs_hcal_goodblock_atime, &b_Ndata_sbs_hcal_goodblock_atime);
   fChain->SetBranchAddress("sbs.hcal.goodblock.atime", sbs_hcal_goodblock_atime, &b_sbs_hcal_goodblock_atime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.cid", &Ndata_sbs_hcal_goodblock_cid, &b_Ndata_sbs_hcal_goodblock_cid);
   fChain->SetBranchAddress("sbs.hcal.goodblock.cid", sbs_hcal_goodblock_cid, &b_sbs_hcal_goodblock_cid);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.col", &Ndata_sbs_hcal_goodblock_col, &b_Ndata_sbs_hcal_goodblock_col);
   fChain->SetBranchAddress("sbs.hcal.goodblock.col", sbs_hcal_goodblock_col, &b_sbs_hcal_goodblock_col);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.e", &Ndata_sbs_hcal_goodblock_e, &b_Ndata_sbs_hcal_goodblock_e);
   fChain->SetBranchAddress("sbs.hcal.goodblock.e", sbs_hcal_goodblock_e, &b_sbs_hcal_goodblock_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.id", &Ndata_sbs_hcal_goodblock_id, &b_Ndata_sbs_hcal_goodblock_id);
   fChain->SetBranchAddress("sbs.hcal.goodblock.id", sbs_hcal_goodblock_id, &b_sbs_hcal_goodblock_id);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.row", &Ndata_sbs_hcal_goodblock_row, &b_Ndata_sbs_hcal_goodblock_row);
   fChain->SetBranchAddress("sbs.hcal.goodblock.row", sbs_hcal_goodblock_row, &b_sbs_hcal_goodblock_row);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.tdctime", &Ndata_sbs_hcal_goodblock_tdctime, &b_Ndata_sbs_hcal_goodblock_tdctime);
   fChain->SetBranchAddress("sbs.hcal.goodblock.tdctime", sbs_hcal_goodblock_tdctime, &b_sbs_hcal_goodblock_tdctime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.tdctime_tw", &Ndata_sbs_hcal_goodblock_tdctime_tw, &b_Ndata_sbs_hcal_goodblock_tdctime_tw);
   fChain->SetBranchAddress("sbs.hcal.goodblock.tdctime_tw", sbs_hcal_goodblock_tdctime_tw, &b_sbs_hcal_goodblock_tdctime_tw);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.x", &Ndata_sbs_hcal_goodblock_x, &b_Ndata_sbs_hcal_goodblock_x);
   fChain->SetBranchAddress("sbs.hcal.goodblock.x", sbs_hcal_goodblock_x, &b_sbs_hcal_goodblock_x);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.y", &Ndata_sbs_hcal_goodblock_y, &b_Ndata_sbs_hcal_goodblock_y);
   fChain->SetBranchAddress("sbs.hcal.goodblock.y", sbs_hcal_goodblock_y, &b_sbs_hcal_goodblock_y);
   fChain->SetBranchAddress("Ndata.sbs.hcal.vtp.clus.e", &Ndata_sbs_hcal_vtp_clus_e, &b_Ndata_sbs_hcal_vtp_clus_e);
   fChain->SetBranchAddress("sbs.hcal.vtp.clus.e", sbs_hcal_vtp_clus_e, &b_sbs_hcal_vtp_clus_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.vtp.clus.size", &Ndata_sbs_hcal_vtp_clus_size, &b_Ndata_sbs_hcal_vtp_clus_size);
   fChain->SetBranchAddress("sbs.hcal.vtp.clus.size", sbs_hcal_vtp_clus_size, &b_sbs_hcal_vtp_clus_size);
   fChain->SetBranchAddress("Ndata.sbs.hcal.vtp.clus.time", &Ndata_sbs_hcal_vtp_clus_time, &b_Ndata_sbs_hcal_vtp_clus_time);
   fChain->SetBranchAddress("sbs.hcal.vtp.clus.time", sbs_hcal_vtp_clus_time, &b_sbs_hcal_vtp_clus_time);
   fChain->SetBranchAddress("Ndata.sbs.hcal.vtp.clus.x", &Ndata_sbs_hcal_vtp_clus_x, &b_Ndata_sbs_hcal_vtp_clus_x);
   fChain->SetBranchAddress("sbs.hcal.vtp.clus.x", sbs_hcal_vtp_clus_x, &b_sbs_hcal_vtp_clus_x);
   fChain->SetBranchAddress("Ndata.sbs.hcal.vtp.clus.y", &Ndata_sbs_hcal_vtp_clus_y, &b_Ndata_sbs_hcal_vtp_clus_y);
   fChain->SetBranchAddress("sbs.hcal.vtp.clus.y", sbs_hcal_vtp_clus_y, &b_sbs_hcal_vtp_clus_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.beta", &Ndata_sbs_tr_beta, &b_Ndata_sbs_tr_beta);
   fChain->SetBranchAddress("sbs.tr.beta", sbs_tr_beta, &b_sbs_tr_beta);
   fChain->SetBranchAddress("Ndata.sbs.tr.chi2", &Ndata_sbs_tr_chi2, &b_Ndata_sbs_tr_chi2);
   fChain->SetBranchAddress("sbs.tr.chi2", sbs_tr_chi2, &b_sbs_tr_chi2);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_ph", &Ndata_sbs_tr_d_ph, &b_Ndata_sbs_tr_d_ph);
   fChain->SetBranchAddress("sbs.tr.d_ph", sbs_tr_d_ph, &b_sbs_tr_d_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_th", &Ndata_sbs_tr_d_th, &b_Ndata_sbs_tr_d_th);
   fChain->SetBranchAddress("sbs.tr.d_th", sbs_tr_d_th, &b_sbs_tr_d_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_x", &Ndata_sbs_tr_d_x, &b_Ndata_sbs_tr_d_x);
   fChain->SetBranchAddress("sbs.tr.d_x", sbs_tr_d_x, &b_sbs_tr_d_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_y", &Ndata_sbs_tr_d_y, &b_Ndata_sbs_tr_d_y);
   fChain->SetBranchAddress("sbs.tr.d_y", sbs_tr_d_y, &b_sbs_tr_d_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.dbeta", &Ndata_sbs_tr_dbeta, &b_Ndata_sbs_tr_dbeta);
   fChain->SetBranchAddress("sbs.tr.dbeta", sbs_tr_dbeta, &b_sbs_tr_dbeta);
   fChain->SetBranchAddress("Ndata.sbs.tr.dtime", &Ndata_sbs_tr_dtime, &b_Ndata_sbs_tr_dtime);
   fChain->SetBranchAddress("sbs.tr.dtime", sbs_tr_dtime, &b_sbs_tr_dtime);
   fChain->SetBranchAddress("Ndata.sbs.tr.flag", &Ndata_sbs_tr_flag, &b_Ndata_sbs_tr_flag);
   fChain->SetBranchAddress("sbs.tr.flag", sbs_tr_flag, &b_sbs_tr_flag);
   fChain->SetBranchAddress("Ndata.sbs.tr.ndof", &Ndata_sbs_tr_ndof, &b_Ndata_sbs_tr_ndof);
   fChain->SetBranchAddress("sbs.tr.ndof", sbs_tr_ndof, &b_sbs_tr_ndof);
   fChain->SetBranchAddress("Ndata.sbs.tr.p", &Ndata_sbs_tr_p, &b_Ndata_sbs_tr_p);
   fChain->SetBranchAddress("sbs.tr.p", sbs_tr_p, &b_sbs_tr_p);
   fChain->SetBranchAddress("Ndata.sbs.tr.pathl", &Ndata_sbs_tr_pathl, &b_Ndata_sbs_tr_pathl);
   fChain->SetBranchAddress("sbs.tr.pathl", sbs_tr_pathl, &b_sbs_tr_pathl);
   fChain->SetBranchAddress("Ndata.sbs.tr.ph", &Ndata_sbs_tr_ph, &b_Ndata_sbs_tr_ph);
   fChain->SetBranchAddress("sbs.tr.ph", sbs_tr_ph, &b_sbs_tr_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.px", &Ndata_sbs_tr_px, &b_Ndata_sbs_tr_px);
   fChain->SetBranchAddress("sbs.tr.px", sbs_tr_px, &b_sbs_tr_px);
   fChain->SetBranchAddress("Ndata.sbs.tr.py", &Ndata_sbs_tr_py, &b_Ndata_sbs_tr_py);
   fChain->SetBranchAddress("sbs.tr.py", sbs_tr_py, &b_sbs_tr_py);
   fChain->SetBranchAddress("Ndata.sbs.tr.pz", &Ndata_sbs_tr_pz, &b_Ndata_sbs_tr_pz);
   fChain->SetBranchAddress("sbs.tr.pz", sbs_tr_pz, &b_sbs_tr_pz);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_ph", &Ndata_sbs_tr_r_ph, &b_Ndata_sbs_tr_r_ph);
   fChain->SetBranchAddress("sbs.tr.r_ph", sbs_tr_r_ph, &b_sbs_tr_r_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_th", &Ndata_sbs_tr_r_th, &b_Ndata_sbs_tr_r_th);
   fChain->SetBranchAddress("sbs.tr.r_th", sbs_tr_r_th, &b_sbs_tr_r_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_x", &Ndata_sbs_tr_r_x, &b_Ndata_sbs_tr_r_x);
   fChain->SetBranchAddress("sbs.tr.r_x", sbs_tr_r_x, &b_sbs_tr_r_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_y", &Ndata_sbs_tr_r_y, &b_Ndata_sbs_tr_r_y);
   fChain->SetBranchAddress("sbs.tr.r_y", sbs_tr_r_y, &b_sbs_tr_r_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_dp", &Ndata_sbs_tr_tg_dp, &b_Ndata_sbs_tr_tg_dp);
   fChain->SetBranchAddress("sbs.tr.tg_dp", sbs_tr_tg_dp, &b_sbs_tr_tg_dp);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_ph", &Ndata_sbs_tr_tg_ph, &b_Ndata_sbs_tr_tg_ph);
   fChain->SetBranchAddress("sbs.tr.tg_ph", sbs_tr_tg_ph, &b_sbs_tr_tg_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_th", &Ndata_sbs_tr_tg_th, &b_Ndata_sbs_tr_tg_th);
   fChain->SetBranchAddress("sbs.tr.tg_th", sbs_tr_tg_th, &b_sbs_tr_tg_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_y", &Ndata_sbs_tr_tg_y, &b_Ndata_sbs_tr_tg_y);
   fChain->SetBranchAddress("sbs.tr.tg_y", sbs_tr_tg_y, &b_sbs_tr_tg_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.th", &Ndata_sbs_tr_th, &b_Ndata_sbs_tr_th);
   fChain->SetBranchAddress("sbs.tr.th", sbs_tr_th, &b_sbs_tr_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.time", &Ndata_sbs_tr_time, &b_Ndata_sbs_tr_time);
   fChain->SetBranchAddress("sbs.tr.time", sbs_tr_time, &b_sbs_tr_time);
   fChain->SetBranchAddress("Ndata.sbs.tr.vx", &Ndata_sbs_tr_vx, &b_Ndata_sbs_tr_vx);
   fChain->SetBranchAddress("sbs.tr.vx", sbs_tr_vx, &b_sbs_tr_vx);
   fChain->SetBranchAddress("Ndata.sbs.tr.vy", &Ndata_sbs_tr_vy, &b_Ndata_sbs_tr_vy);
   fChain->SetBranchAddress("sbs.tr.vy", sbs_tr_vy, &b_sbs_tr_vy);
   fChain->SetBranchAddress("Ndata.sbs.tr.vz", &Ndata_sbs_tr_vz, &b_Ndata_sbs_tr_vz);
   fChain->SetBranchAddress("sbs.tr.vz", sbs_tr_vz, &b_sbs_tr_vz);
   fChain->SetBranchAddress("Ndata.sbs.tr.x", &Ndata_sbs_tr_x, &b_Ndata_sbs_tr_x);
   fChain->SetBranchAddress("sbs.tr.x", sbs_tr_x, &b_sbs_tr_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.y", &Ndata_sbs_tr_y, &b_Ndata_sbs_tr_y);
   fChain->SetBranchAddress("sbs.tr.y", sbs_tr_y, &b_sbs_tr_y);
   fChain->SetBranchAddress("Ndata.sbs.trig.a_amp_p", &Ndata_sbs_trig_a_amp_p, &b_Ndata_sbs_trig_a_amp_p);
   fChain->SetBranchAddress("sbs.trig.a_amp_p", sbs_trig_a_amp_p, &b_sbs_trig_a_amp_p);
   fChain->SetBranchAddress("Ndata.sbs.trig.a_p", &Ndata_sbs_trig_a_p, &b_Ndata_sbs_trig_a_p);
   fChain->SetBranchAddress("sbs.trig.a_p", sbs_trig_a_p, &b_sbs_trig_a_p);
   fChain->SetBranchAddress("Ndata.sbs.trig.a_time", &Ndata_sbs_trig_a_time, &b_Ndata_sbs_trig_a_time);
   fChain->SetBranchAddress("sbs.trig.a_time", sbs_trig_a_time, &b_sbs_trig_a_time);
   fChain->SetBranchAddress("Ndata.sbs.trig.adcelemID", &Ndata_sbs_trig_adcelemID, &b_Ndata_sbs_trig_adcelemID);
   fChain->SetBranchAddress("sbs.trig.adcelemID", sbs_trig_adcelemID, &b_sbs_trig_adcelemID);
   fChain->SetBranchAddress("Ndata.sbs.x_bcp", &Ndata_sbs_x_bcp, &b_Ndata_sbs_x_bcp);
   fChain->SetBranchAddress("sbs.x_bcp", sbs_x_bcp, &b_sbs_x_bcp);
   fChain->SetBranchAddress("Ndata.sbs.x_fcp", &Ndata_sbs_x_fcp, &b_Ndata_sbs_x_fcp);
   fChain->SetBranchAddress("sbs.x_fcp", sbs_x_fcp, &b_sbs_x_fcp);
   fChain->SetBranchAddress("Ndata.sbs.y_bcp", &Ndata_sbs_y_bcp, &b_Ndata_sbs_y_bcp);
   fChain->SetBranchAddress("sbs.y_bcp", sbs_y_bcp, &b_sbs_y_bcp);
   fChain->SetBranchAddress("Ndata.sbs.y_fcp", &Ndata_sbs_y_fcp, &b_Ndata_sbs_y_fcp);
   fChain->SetBranchAddress("sbs.y_fcp", sbs_y_fcp, &b_sbs_y_fcp);
   fChain->SetBranchAddress("Ndata.sbs.z_bcp", &Ndata_sbs_z_bcp, &b_Ndata_sbs_z_bcp);
   fChain->SetBranchAddress("sbs.z_bcp", sbs_z_bcp, &b_sbs_z_bcp);
   fChain->SetBranchAddress("Ndata.sbs.z_fcp", &Ndata_sbs_z_fcp, &b_Ndata_sbs_z_fcp);
   fChain->SetBranchAddress("sbs.z_fcp", sbs_z_fcp, &b_sbs_z_fcp);
   fChain->SetBranchAddress("FTROI.ECAL_energy", &FTROI_ECAL_energy, &b_FTROI_ECAL_energy);
   fChain->SetBranchAddress("FTROI.Eprime", &FTROI_Eprime, &b_FTROI_Eprime);
   fChain->SetBranchAddress("FTROI.ephi", &FTROI_ephi, &b_FTROI_ephi);
   fChain->SetBranchAddress("FTROI.etheta", &FTROI_etheta, &b_FTROI_etheta);
   fChain->SetBranchAddress("FTROI.good", &FTROI_good, &b_FTROI_good);
   fChain->SetBranchAddress("FTROI.pp", &FTROI_pp, &b_FTROI_pp);
   fChain->SetBranchAddress("FTROI.pphi", &FTROI_pphi, &b_FTROI_pphi);
   fChain->SetBranchAddress("FTROI.ptheta", &FTROI_ptheta, &b_FTROI_ptheta);
   fChain->SetBranchAddress("FTROI.xECAL_global", &FTROI_xECAL_global, &b_FTROI_xECAL_global);
   fChain->SetBranchAddress("FTROI.yECAL_global", &FTROI_yECAL_global, &b_FTROI_yECAL_global);
   fChain->SetBranchAddress("FTROI.zECAL_global", &FTROI_zECAL_global, &b_FTROI_zECAL_global);
   fChain->SetBranchAddress("Lrb.BPMA.rawcur.1", &Lrb_BPMA_rawcur_1, &b_Lrb_BPMA_rawcur_1);
   fChain->SetBranchAddress("Lrb.BPMA.rawcur.2", &Lrb_BPMA_rawcur_2, &b_Lrb_BPMA_rawcur_2);
   fChain->SetBranchAddress("Lrb.BPMA.rawcur.3", &Lrb_BPMA_rawcur_3, &b_Lrb_BPMA_rawcur_3);
   fChain->SetBranchAddress("Lrb.BPMA.rawcur.4", &Lrb_BPMA_rawcur_4, &b_Lrb_BPMA_rawcur_4);
   fChain->SetBranchAddress("Lrb.BPMA.rotpos1", &Lrb_BPMA_rotpos1, &b_Lrb_BPMA_rotpos1);
   fChain->SetBranchAddress("Lrb.BPMA.rotpos2", &Lrb_BPMA_rotpos2, &b_Lrb_BPMA_rotpos2);
   fChain->SetBranchAddress("Lrb.BPMA.x", &Lrb_BPMA_x, &b_Lrb_BPMA_x);
   fChain->SetBranchAddress("Lrb.BPMA.y", &Lrb_BPMA_y, &b_Lrb_BPMA_y);
   fChain->SetBranchAddress("Lrb.BPMA.z", &Lrb_BPMA_z, &b_Lrb_BPMA_z);
   fChain->SetBranchAddress("Lrb.BPMB.rawcur.1", &Lrb_BPMB_rawcur_1, &b_Lrb_BPMB_rawcur_1);
   fChain->SetBranchAddress("Lrb.BPMB.rawcur.2", &Lrb_BPMB_rawcur_2, &b_Lrb_BPMB_rawcur_2);
   fChain->SetBranchAddress("Lrb.BPMB.rawcur.3", &Lrb_BPMB_rawcur_3, &b_Lrb_BPMB_rawcur_3);
   fChain->SetBranchAddress("Lrb.BPMB.rawcur.4", &Lrb_BPMB_rawcur_4, &b_Lrb_BPMB_rawcur_4);
   fChain->SetBranchAddress("Lrb.BPMB.rotpos1", &Lrb_BPMB_rotpos1, &b_Lrb_BPMB_rotpos1);
   fChain->SetBranchAddress("Lrb.BPMB.rotpos2", &Lrb_BPMB_rotpos2, &b_Lrb_BPMB_rotpos2);
   fChain->SetBranchAddress("Lrb.BPMB.x", &Lrb_BPMB_x, &b_Lrb_BPMB_x);
   fChain->SetBranchAddress("Lrb.BPMB.y", &Lrb_BPMB_y, &b_Lrb_BPMB_y);
   fChain->SetBranchAddress("Lrb.BPMB.z", &Lrb_BPMB_z, &b_Lrb_BPMB_z);
   fChain->SetBranchAddress("Lrb.Raster.bpma.x", &Lrb_Raster_bpma_x, &b_Lrb_Raster_bpma_x);
   fChain->SetBranchAddress("Lrb.Raster.bpma.y", &Lrb_Raster_bpma_y, &b_Lrb_Raster_bpma_y);
   fChain->SetBranchAddress("Lrb.Raster.bpma.z", &Lrb_Raster_bpma_z, &b_Lrb_Raster_bpma_z);
   fChain->SetBranchAddress("Lrb.Raster.bpmb.x", &Lrb_Raster_bpmb_x, &b_Lrb_Raster_bpmb_x);
   fChain->SetBranchAddress("Lrb.Raster.bpmb.y", &Lrb_Raster_bpmb_y, &b_Lrb_Raster_bpmb_y);
   fChain->SetBranchAddress("Lrb.Raster.bpmb.z", &Lrb_Raster_bpmb_z, &b_Lrb_Raster_bpmb_z);
   fChain->SetBranchAddress("Lrb.Raster.rawcur.x", &Lrb_Raster_rawcur_x, &b_Lrb_Raster_rawcur_x);
   fChain->SetBranchAddress("Lrb.Raster.rawcur.y", &Lrb_Raster_rawcur_y, &b_Lrb_Raster_rawcur_y);
   fChain->SetBranchAddress("Lrb.Raster.rawslope.x", &Lrb_Raster_rawslope_x, &b_Lrb_Raster_rawslope_x);
   fChain->SetBranchAddress("Lrb.Raster.rawslope.y", &Lrb_Raster_rawslope_y, &b_Lrb_Raster_rawslope_y);
   fChain->SetBranchAddress("Lrb.Raster.target.dir.x", &Lrb_Raster_target_dir_x, &b_Lrb_Raster_target_dir_x);
   fChain->SetBranchAddress("Lrb.Raster.target.dir.y", &Lrb_Raster_target_dir_y, &b_Lrb_Raster_target_dir_y);
   fChain->SetBranchAddress("Lrb.Raster.target.dir.z", &Lrb_Raster_target_dir_z, &b_Lrb_Raster_target_dir_z);
   fChain->SetBranchAddress("Lrb.Raster.target.x", &Lrb_Raster_target_x, &b_Lrb_Raster_target_x);
   fChain->SetBranchAddress("Lrb.Raster.target.y", &Lrb_Raster_target_y, &b_Lrb_Raster_target_y);
   fChain->SetBranchAddress("Lrb.Raster.target.z", &Lrb_Raster_target_z, &b_Lrb_Raster_target_z);
   fChain->SetBranchAddress("Lrb.Raster2.bpma.x", &Lrb_Raster2_bpma_x, &b_Lrb_Raster2_bpma_x);
   fChain->SetBranchAddress("Lrb.Raster2.bpma.y", &Lrb_Raster2_bpma_y, &b_Lrb_Raster2_bpma_y);
   fChain->SetBranchAddress("Lrb.Raster2.bpma.z", &Lrb_Raster2_bpma_z, &b_Lrb_Raster2_bpma_z);
   fChain->SetBranchAddress("Lrb.Raster2.bpmb.x", &Lrb_Raster2_bpmb_x, &b_Lrb_Raster2_bpmb_x);
   fChain->SetBranchAddress("Lrb.Raster2.bpmb.y", &Lrb_Raster2_bpmb_y, &b_Lrb_Raster2_bpmb_y);
   fChain->SetBranchAddress("Lrb.Raster2.bpmb.z", &Lrb_Raster2_bpmb_z, &b_Lrb_Raster2_bpmb_z);
   fChain->SetBranchAddress("Lrb.Raster2.rawcur.x", &Lrb_Raster2_rawcur_x, &b_Lrb_Raster2_rawcur_x);
   fChain->SetBranchAddress("Lrb.Raster2.rawcur.y", &Lrb_Raster2_rawcur_y, &b_Lrb_Raster2_rawcur_y);
   fChain->SetBranchAddress("Lrb.Raster2.rawslope.x", &Lrb_Raster2_rawslope_x, &b_Lrb_Raster2_rawslope_x);
   fChain->SetBranchAddress("Lrb.Raster2.rawslope.y", &Lrb_Raster2_rawslope_y, &b_Lrb_Raster2_rawslope_y);
   fChain->SetBranchAddress("Lrb.Raster2.target.dir.x", &Lrb_Raster2_target_dir_x, &b_Lrb_Raster2_target_dir_x);
   fChain->SetBranchAddress("Lrb.Raster2.target.dir.y", &Lrb_Raster2_target_dir_y, &b_Lrb_Raster2_target_dir_y);
   fChain->SetBranchAddress("Lrb.Raster2.target.dir.z", &Lrb_Raster2_target_dir_z, &b_Lrb_Raster2_target_dir_z);
   fChain->SetBranchAddress("Lrb.Raster2.target.x", &Lrb_Raster2_target_x, &b_Lrb_Raster2_target_x);
   fChain->SetBranchAddress("Lrb.Raster2.target.y", &Lrb_Raster2_target_y, &b_Lrb_Raster2_target_y);
   fChain->SetBranchAddress("Lrb.Raster2.target.z", &Lrb_Raster2_target_z, &b_Lrb_Raster2_target_z);
   fChain->SetBranchAddress("SBS.gold.beta", &SBS_gold_beta, &b_SBS_gold_beta);
   fChain->SetBranchAddress("SBS.gold.dp", &SBS_gold_dp, &b_SBS_gold_dp);
   fChain->SetBranchAddress("SBS.gold.index", &SBS_gold_index, &b_SBS_gold_index);
   fChain->SetBranchAddress("SBS.gold.ok", &SBS_gold_ok, &b_SBS_gold_ok);
   fChain->SetBranchAddress("SBS.gold.p", &SBS_gold_p, &b_SBS_gold_p);
   fChain->SetBranchAddress("SBS.gold.ph", &SBS_gold_ph, &b_SBS_gold_ph);
   fChain->SetBranchAddress("SBS.gold.px", &SBS_gold_px, &b_SBS_gold_px);
   fChain->SetBranchAddress("SBS.gold.py", &SBS_gold_py, &b_SBS_gold_py);
   fChain->SetBranchAddress("SBS.gold.pz", &SBS_gold_pz, &b_SBS_gold_pz);
   fChain->SetBranchAddress("SBS.gold.th", &SBS_gold_th, &b_SBS_gold_th);
   fChain->SetBranchAddress("SBS.gold.x", &SBS_gold_x, &b_SBS_gold_x);
   fChain->SetBranchAddress("SBS.gold.y", &SBS_gold_y, &b_SBS_gold_y);
   fChain->SetBranchAddress("SBSrb.BPMA.rawcur.1", &SBSrb_BPMA_rawcur_1, &b_SBSrb_BPMA_rawcur_1);
   fChain->SetBranchAddress("SBSrb.BPMA.rawcur.2", &SBSrb_BPMA_rawcur_2, &b_SBSrb_BPMA_rawcur_2);
   fChain->SetBranchAddress("SBSrb.BPMA.rawcur.3", &SBSrb_BPMA_rawcur_3, &b_SBSrb_BPMA_rawcur_3);
   fChain->SetBranchAddress("SBSrb.BPMA.rawcur.4", &SBSrb_BPMA_rawcur_4, &b_SBSrb_BPMA_rawcur_4);
   fChain->SetBranchAddress("SBSrb.BPMA.rotpos1", &SBSrb_BPMA_rotpos1, &b_SBSrb_BPMA_rotpos1);
   fChain->SetBranchAddress("SBSrb.BPMA.rotpos2", &SBSrb_BPMA_rotpos2, &b_SBSrb_BPMA_rotpos2);
   fChain->SetBranchAddress("SBSrb.BPMA.x", &SBSrb_BPMA_x, &b_SBSrb_BPMA_x);
   fChain->SetBranchAddress("SBSrb.BPMA.y", &SBSrb_BPMA_y, &b_SBSrb_BPMA_y);
   fChain->SetBranchAddress("SBSrb.BPMA.z", &SBSrb_BPMA_z, &b_SBSrb_BPMA_z);
   fChain->SetBranchAddress("SBSrb.BPMB.rawcur.1", &SBSrb_BPMB_rawcur_1, &b_SBSrb_BPMB_rawcur_1);
   fChain->SetBranchAddress("SBSrb.BPMB.rawcur.2", &SBSrb_BPMB_rawcur_2, &b_SBSrb_BPMB_rawcur_2);
   fChain->SetBranchAddress("SBSrb.BPMB.rawcur.3", &SBSrb_BPMB_rawcur_3, &b_SBSrb_BPMB_rawcur_3);
   fChain->SetBranchAddress("SBSrb.BPMB.rawcur.4", &SBSrb_BPMB_rawcur_4, &b_SBSrb_BPMB_rawcur_4);
   fChain->SetBranchAddress("SBSrb.BPMB.rotpos1", &SBSrb_BPMB_rotpos1, &b_SBSrb_BPMB_rotpos1);
   fChain->SetBranchAddress("SBSrb.BPMB.rotpos2", &SBSrb_BPMB_rotpos2, &b_SBSrb_BPMB_rotpos2);
   fChain->SetBranchAddress("SBSrb.BPMB.x", &SBSrb_BPMB_x, &b_SBSrb_BPMB_x);
   fChain->SetBranchAddress("SBSrb.BPMB.y", &SBSrb_BPMB_y, &b_SBSrb_BPMB_y);
   fChain->SetBranchAddress("SBSrb.BPMB.z", &SBSrb_BPMB_z, &b_SBSrb_BPMB_z);
   fChain->SetBranchAddress("SBSrb.Raster.bpma.x", &SBSrb_Raster_bpma_x, &b_SBSrb_Raster_bpma_x);
   fChain->SetBranchAddress("SBSrb.Raster.bpma.y", &SBSrb_Raster_bpma_y, &b_SBSrb_Raster_bpma_y);
   fChain->SetBranchAddress("SBSrb.Raster.bpma.z", &SBSrb_Raster_bpma_z, &b_SBSrb_Raster_bpma_z);
   fChain->SetBranchAddress("SBSrb.Raster.bpmb.x", &SBSrb_Raster_bpmb_x, &b_SBSrb_Raster_bpmb_x);
   fChain->SetBranchAddress("SBSrb.Raster.bpmb.y", &SBSrb_Raster_bpmb_y, &b_SBSrb_Raster_bpmb_y);
   fChain->SetBranchAddress("SBSrb.Raster.bpmb.z", &SBSrb_Raster_bpmb_z, &b_SBSrb_Raster_bpmb_z);
   fChain->SetBranchAddress("SBSrb.Raster.rawcur.x", &SBSrb_Raster_rawcur_x, &b_SBSrb_Raster_rawcur_x);
   fChain->SetBranchAddress("SBSrb.Raster.rawcur.y", &SBSrb_Raster_rawcur_y, &b_SBSrb_Raster_rawcur_y);
   fChain->SetBranchAddress("SBSrb.Raster.rawslope.x", &SBSrb_Raster_rawslope_x, &b_SBSrb_Raster_rawslope_x);
   fChain->SetBranchAddress("SBSrb.Raster.rawslope.y", &SBSrb_Raster_rawslope_y, &b_SBSrb_Raster_rawslope_y);
   fChain->SetBranchAddress("SBSrb.Raster.target.dir.x", &SBSrb_Raster_target_dir_x, &b_SBSrb_Raster_target_dir_x);
   fChain->SetBranchAddress("SBSrb.Raster.target.dir.y", &SBSrb_Raster_target_dir_y, &b_SBSrb_Raster_target_dir_y);
   fChain->SetBranchAddress("SBSrb.Raster.target.dir.z", &SBSrb_Raster_target_dir_z, &b_SBSrb_Raster_target_dir_z);
   fChain->SetBranchAddress("SBSrb.Raster.target.x", &SBSrb_Raster_target_x, &b_SBSrb_Raster_target_x);
   fChain->SetBranchAddress("SBSrb.Raster.target.y", &SBSrb_Raster_target_y, &b_SBSrb_Raster_target_y);
   fChain->SetBranchAddress("SBSrb.Raster.target.z", &SBSrb_Raster_target_z, &b_SBSrb_Raster_target_z);
   fChain->SetBranchAddress("SBSrb.Raster2.bpma.x", &SBSrb_Raster2_bpma_x, &b_SBSrb_Raster2_bpma_x);
   fChain->SetBranchAddress("SBSrb.Raster2.bpma.y", &SBSrb_Raster2_bpma_y, &b_SBSrb_Raster2_bpma_y);
   fChain->SetBranchAddress("SBSrb.Raster2.bpma.z", &SBSrb_Raster2_bpma_z, &b_SBSrb_Raster2_bpma_z);
   fChain->SetBranchAddress("SBSrb.Raster2.bpmb.x", &SBSrb_Raster2_bpmb_x, &b_SBSrb_Raster2_bpmb_x);
   fChain->SetBranchAddress("SBSrb.Raster2.bpmb.y", &SBSrb_Raster2_bpmb_y, &b_SBSrb_Raster2_bpmb_y);
   fChain->SetBranchAddress("SBSrb.Raster2.bpmb.z", &SBSrb_Raster2_bpmb_z, &b_SBSrb_Raster2_bpmb_z);
   fChain->SetBranchAddress("SBSrb.Raster2.rawcur.x", &SBSrb_Raster2_rawcur_x, &b_SBSrb_Raster2_rawcur_x);
   fChain->SetBranchAddress("SBSrb.Raster2.rawcur.y", &SBSrb_Raster2_rawcur_y, &b_SBSrb_Raster2_rawcur_y);
   fChain->SetBranchAddress("SBSrb.Raster2.rawslope.x", &SBSrb_Raster2_rawslope_x, &b_SBSrb_Raster2_rawslope_x);
   fChain->SetBranchAddress("SBSrb.Raster2.rawslope.y", &SBSrb_Raster2_rawslope_y, &b_SBSrb_Raster2_rawslope_y);
   fChain->SetBranchAddress("SBSrb.Raster2.target.dir.x", &SBSrb_Raster2_target_dir_x, &b_SBSrb_Raster2_target_dir_x);
   fChain->SetBranchAddress("SBSrb.Raster2.target.dir.y", &SBSrb_Raster2_target_dir_y, &b_SBSrb_Raster2_target_dir_y);
   fChain->SetBranchAddress("SBSrb.Raster2.target.dir.z", &SBSrb_Raster2_target_dir_z, &b_SBSrb_Raster2_target_dir_z);
   fChain->SetBranchAddress("SBSrb.Raster2.target.x", &SBSrb_Raster2_target_x, &b_SBSrb_Raster2_target_x);
   fChain->SetBranchAddress("SBSrb.Raster2.target.y", &SBSrb_Raster2_target_y, &b_SBSrb_Raster2_target_y);
   fChain->SetBranchAddress("SBSrb.Raster2.target.z", &SBSrb_Raster2_target_z, &b_SBSrb_Raster2_target_z);
   fChain->SetBranchAddress("earm.cdet.ngoodADChits", &earm_cdet_ngoodADChits, &b_earm_cdet_ngoodADChits);
   fChain->SetBranchAddress("earm.cdet.ngoodTDChits", &earm_cdet_ngoodTDChits, &b_earm_cdet_ngoodTDChits);
   fChain->SetBranchAddress("earm.cdet.ngoodhits", &earm_cdet_ngoodhits, &b_earm_cdet_ngoodhits);
   fChain->SetBranchAddress("earm.cdet.nhits", &earm_cdet_nhits, &b_earm_cdet_nhits);
   fChain->SetBranchAddress("earm.cdet.nrefhits", &earm_cdet_nrefhits, &b_earm_cdet_nrefhits);
   fChain->SetBranchAddress("earm.ecal.adctime", &earm_ecal_adctime, &b_earm_ecal_adctime);
   fChain->SetBranchAddress("earm.ecal.againblk", &earm_ecal_againblk, &b_earm_ecal_againblk);
   fChain->SetBranchAddress("earm.ecal.atimeblk", &earm_ecal_atimeblk, &b_earm_ecal_atimeblk);
   fChain->SetBranchAddress("earm.ecal.colblk", &earm_ecal_colblk, &b_earm_ecal_colblk);
   fChain->SetBranchAddress("earm.ecal.e", &earm_ecal_e, &b_earm_ecal_e);
   fChain->SetBranchAddress("earm.ecal.eblk", &earm_ecal_eblk, &b_earm_ecal_eblk);
   fChain->SetBranchAddress("earm.ecal.idblk", &earm_ecal_idblk, &b_earm_ecal_idblk);
   fChain->SetBranchAddress("earm.ecal.index", &earm_ecal_index, &b_earm_ecal_index);
   fChain->SetBranchAddress("earm.ecal.nblk", &earm_ecal_nblk, &b_earm_ecal_nblk);
   fChain->SetBranchAddress("earm.ecal.nclus", &earm_ecal_nclus, &b_earm_ecal_nclus);
   fChain->SetBranchAddress("earm.ecal.ngoodADChits", &earm_ecal_ngoodADChits, &b_earm_ecal_ngoodADChits);
   fChain->SetBranchAddress("earm.ecal.rowblk", &earm_ecal_rowblk, &b_earm_ecal_rowblk);
   fChain->SetBranchAddress("earm.ecal.vtp.detid", &earm_ecal_vtp_detid, &b_earm_ecal_vtp_detid);
   fChain->SetBranchAddress("earm.ecal.vtp.errflag", &earm_ecal_vtp_errflag, &b_earm_ecal_vtp_errflag);
   fChain->SetBranchAddress("earm.ecal.x", &earm_ecal_x, &b_earm_ecal_x);
   fChain->SetBranchAddress("earm.ecal.y", &earm_ecal_y, &b_earm_ecal_y);
   fChain->SetBranchAddress("earm.tr.n", &earm_tr_n, &b_earm_tr_n);
   fChain->SetBranchAddress("g.datatype", &g_datatype, &b_g_datatype);
   fChain->SetBranchAddress("g.evlen", &g_evlen, &b_g_evlen);
   fChain->SetBranchAddress("g.evnum", &g_evnum, &b_g_evnum);
   fChain->SetBranchAddress("g.evtime", &g_evtime, &b_g_evtime);
   fChain->SetBranchAddress("g.evtyp", &g_evtyp, &b_g_evtyp);
   fChain->SetBranchAddress("g.runnum", &g_runnum, &b_g_runnum);
   fChain->SetBranchAddress("g.runtime", &g_runtime, &b_g_runtime);
   fChain->SetBranchAddress("g.runtype", &g_runtype, &b_g_runtype);
   fChain->SetBranchAddress("g.trigbits", &g_trigbits, &b_g_trigbits);
   fChain->SetBranchAddress("heep.K_eth", &heep_K_eth, &b_heep_K_eth);
   fChain->SetBranchAddress("heep.K_p4vect", &heep_K_p4vect, &b_heep_K_p4vect);
   fChain->SetBranchAddress("heep.K_pp", &heep_K_pp, &b_heep_K_pp);
   fChain->SetBranchAddress("heep.K_pth", &heep_K_pth, &b_heep_K_pth);
   fChain->SetBranchAddress("heep.Q2_e4vect", &heep_Q2_e4vect, &b_heep_Q2_e4vect);
   fChain->SetBranchAddress("heep.Q2_eth", &heep_Q2_eth, &b_heep_Q2_eth);
   fChain->SetBranchAddress("heep.Q2_p4vect", &heep_Q2_p4vect, &b_heep_Q2_p4vect);
   fChain->SetBranchAddress("heep.Q2_pp", &heep_Q2_pp, &b_heep_Q2_pp);
   fChain->SetBranchAddress("heep.Q2_pth", &heep_Q2_pth, &b_heep_Q2_pth);
   fChain->SetBranchAddress("heep.acoplanarity", &heep_acoplanarity, &b_heep_acoplanarity);
   fChain->SetBranchAddress("heep.datavalid", &heep_datavalid, &b_heep_datavalid);
   fChain->SetBranchAddress("heep.dpe", &heep_dpe, &b_heep_dpe);
   fChain->SetBranchAddress("heep.dphi", &heep_dphi, &b_heep_dphi);
   fChain->SetBranchAddress("heep.dpp", &heep_dpp, &b_heep_dpp);
   fChain->SetBranchAddress("heep.dt_ADC", &heep_dt_ADC, &b_heep_dt_ADC);
   fChain->SetBranchAddress("heep.dt_TDC", &heep_dt_TDC, &b_heep_dt_TDC);
   fChain->SetBranchAddress("heep.dxECAL", &heep_dxECAL, &b_heep_dxECAL);
   fChain->SetBranchAddress("heep.dxECAL4vect", &heep_dxECAL4vect, &b_heep_dxECAL4vect);
   fChain->SetBranchAddress("heep.dyECAL", &heep_dyECAL, &b_heep_dyECAL);
   fChain->SetBranchAddress("heep.dyECAL4vect", &heep_dyECAL4vect, &b_heep_dyECAL4vect);
   fChain->SetBranchAddress("heep.ecalo", &heep_ecalo, &b_heep_ecalo);
   fChain->SetBranchAddress("heep.ephi", &heep_ephi, &b_heep_ephi);
   fChain->SetBranchAddress("heep.eprime_eth", &heep_eprime_eth, &b_heep_eprime_eth);
   fChain->SetBranchAddress("heep.eprime_pp", &heep_eprime_pp, &b_heep_eprime_pp);
   fChain->SetBranchAddress("heep.eprime_pth", &heep_eprime_pth, &b_heep_eprime_pth);
   fChain->SetBranchAddress("heep.eps_4vect", &heep_eps_4vect, &b_heep_eps_4vect);
   fChain->SetBranchAddress("heep.eps_eth", &heep_eps_eth, &b_heep_eps_eth);
   fChain->SetBranchAddress("heep.eps_pp", &heep_eps_pp, &b_heep_eps_pp);
   fChain->SetBranchAddress("heep.eps_pth", &heep_eps_pth, &b_heep_eps_pth);
   fChain->SetBranchAddress("heep.eth_pp", &heep_eth_pp, &b_heep_eth_pp);
   fChain->SetBranchAddress("heep.eth_pth", &heep_eth_pth, &b_heep_eth_pth);
   fChain->SetBranchAddress("heep.etheta", &heep_etheta, &b_heep_etheta);
   fChain->SetBranchAddress("heep.inelasticity", &heep_inelasticity, &b_heep_inelasticity);
   fChain->SetBranchAddress("heep.pp", &heep_pp, &b_heep_pp);
   fChain->SetBranchAddress("heep.pp_eth", &heep_pp_eth, &b_heep_pp_eth);
   fChain->SetBranchAddress("heep.pp_pth", &heep_pp_pth, &b_heep_pp_pth);
   fChain->SetBranchAddress("heep.pphi", &heep_pphi, &b_heep_pphi);
   fChain->SetBranchAddress("heep.pth_eth", &heep_pth_eth, &b_heep_pth_eth);
   fChain->SetBranchAddress("heep.ptheta", &heep_ptheta, &b_heep_ptheta);
   fChain->SetBranchAddress("sbs.HCALdir_x", &sbs_HCALdir_x, &b_sbs_HCALdir_x);
   fChain->SetBranchAddress("sbs.HCALdir_y", &sbs_HCALdir_y, &b_sbs_HCALdir_y);
   fChain->SetBranchAddress("sbs.HCALdir_z", &sbs_HCALdir_z, &b_sbs_HCALdir_z);
   fChain->SetBranchAddress("sbs.HCALph_n", &sbs_HCALph_n, &b_sbs_HCALph_n);
   fChain->SetBranchAddress("sbs.HCALth_n", &sbs_HCALth_n, &b_sbs_HCALth_n);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ngoodhits", &sbs_gemFPP_hit_ngoodhits, &b_sbs_gemFPP_hit_ngoodhits);
   fChain->SetBranchAddress("sbs.gemFPP.nlayershit", &sbs_gemFPP_nlayershit, &b_sbs_gemFPP_nlayershit);
   fChain->SetBranchAddress("sbs.gemFPP.nlayershitu", &sbs_gemFPP_nlayershitu, &b_sbs_gemFPP_nlayershitu);
   fChain->SetBranchAddress("sbs.gemFPP.nlayershituv", &sbs_gemFPP_nlayershituv, &b_sbs_gemFPP_nlayershituv);
   fChain->SetBranchAddress("sbs.gemFPP.nlayershitv", &sbs_gemFPP_nlayershitv, &b_sbs_gemFPP_nlayershitv);
   fChain->SetBranchAddress("sbs.gemFPP.track.besttrack", &sbs_gemFPP_track_besttrack, &b_sbs_gemFPP_track_besttrack);
   fChain->SetBranchAddress("sbs.gemFPP.track.ntrack", &sbs_gemFPP_track_ntrack, &b_sbs_gemFPP_track_ntrack);
   fChain->SetBranchAddress("sbs.gemFT.hit.ngoodhits", &sbs_gemFT_hit_ngoodhits, &b_sbs_gemFT_hit_ngoodhits);
   fChain->SetBranchAddress("sbs.gemFT.nlayershit", &sbs_gemFT_nlayershit, &b_sbs_gemFT_nlayershit);
   fChain->SetBranchAddress("sbs.gemFT.nlayershitu", &sbs_gemFT_nlayershitu, &b_sbs_gemFT_nlayershitu);
   fChain->SetBranchAddress("sbs.gemFT.nlayershituv", &sbs_gemFT_nlayershituv, &b_sbs_gemFT_nlayershituv);
   fChain->SetBranchAddress("sbs.gemFT.nlayershitv", &sbs_gemFT_nlayershitv, &b_sbs_gemFT_nlayershitv);
   fChain->SetBranchAddress("sbs.gemFT.track.besttrack", &sbs_gemFT_track_besttrack, &b_sbs_gemFT_track_besttrack);
   fChain->SetBranchAddress("sbs.gemFT.track.ntrack", &sbs_gemFT_track_ntrack, &b_sbs_gemFT_track_ntrack);
   fChain->SetBranchAddress("sbs.hcal.adctime", &sbs_hcal_adctime, &b_sbs_hcal_adctime);
   fChain->SetBranchAddress("sbs.hcal.againblk", &sbs_hcal_againblk, &b_sbs_hcal_againblk);
   fChain->SetBranchAddress("sbs.hcal.atimeblk", &sbs_hcal_atimeblk, &b_sbs_hcal_atimeblk);
   fChain->SetBranchAddress("sbs.hcal.colblk", &sbs_hcal_colblk, &b_sbs_hcal_colblk);
   fChain->SetBranchAddress("sbs.hcal.colblk_goodtdc", &sbs_hcal_colblk_goodtdc, &b_sbs_hcal_colblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.e", &sbs_hcal_e, &b_sbs_hcal_e);
   fChain->SetBranchAddress("sbs.hcal.e_goodtdc", &sbs_hcal_e_goodtdc, &b_sbs_hcal_e_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.eblk", &sbs_hcal_eblk, &b_sbs_hcal_eblk);
   fChain->SetBranchAddress("sbs.hcal.eblk_goodtdc", &sbs_hcal_eblk_goodtdc, &b_sbs_hcal_eblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.idblk", &sbs_hcal_idblk, &b_sbs_hcal_idblk);
   fChain->SetBranchAddress("sbs.hcal.idblk_goodtdc", &sbs_hcal_idblk_goodtdc, &b_sbs_hcal_idblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.index", &sbs_hcal_index, &b_sbs_hcal_index);
   fChain->SetBranchAddress("sbs.hcal.nblk", &sbs_hcal_nblk, &b_sbs_hcal_nblk);
   fChain->SetBranchAddress("sbs.hcal.nblk_goodtdc", &sbs_hcal_nblk_goodtdc, &b_sbs_hcal_nblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus, &b_sbs_hcal_nclus);
   fChain->SetBranchAddress("sbs.hcal.ngoodADChits", &sbs_hcal_ngoodADChits, &b_sbs_hcal_ngoodADChits);
   fChain->SetBranchAddress("sbs.hcal.rowblk", &sbs_hcal_rowblk, &b_sbs_hcal_rowblk);
   fChain->SetBranchAddress("sbs.hcal.rowblk_goodtdc", &sbs_hcal_rowblk_goodtdc, &b_sbs_hcal_rowblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.tdctime", &sbs_hcal_tdctime, &b_sbs_hcal_tdctime);
   fChain->SetBranchAddress("sbs.hcal.tdctime_tw", &sbs_hcal_tdctime_tw, &b_sbs_hcal_tdctime_tw);
   fChain->SetBranchAddress("sbs.hcal.tdctimeblk", &sbs_hcal_tdctimeblk, &b_sbs_hcal_tdctimeblk);
   fChain->SetBranchAddress("sbs.hcal.tdctimeblk_tw", &sbs_hcal_tdctimeblk_tw, &b_sbs_hcal_tdctimeblk_tw);
   fChain->SetBranchAddress("sbs.hcal.trigtime", &sbs_hcal_trigtime, &b_sbs_hcal_trigtime);
   fChain->SetBranchAddress("sbs.hcal.vtp.detid", &sbs_hcal_vtp_detid, &b_sbs_hcal_vtp_detid);
   fChain->SetBranchAddress("sbs.hcal.vtp.errflag", &sbs_hcal_vtp_errflag, &b_sbs_hcal_vtp_errflag);
   fChain->SetBranchAddress("sbs.hcal.x", &sbs_hcal_x, &b_sbs_hcal_x);
   fChain->SetBranchAddress("sbs.hcal.y", &sbs_hcal_y, &b_sbs_hcal_y);
   fChain->SetBranchAddress("sbs.tr.n", &sbs_tr_n, &b_sbs_tr_n);
   fChain->SetBranchAddress("scalhel.errcode", &scalhel_errcode, &b_scalhel_errcode);
   fChain->SetBranchAddress("scalhel.evtcount", &scalhel_evtcount, &b_scalhel_evtcount);
   fChain->SetBranchAddress("scalhel.hel", &scalhel_hel, &b_scalhel_hel);
   fChain->SetBranchAddress("scalhel.lhrs.fadc.hel", &scalhel_lhrs_fadc_hel, &b_scalhel_lhrs_fadc_hel);
   fChain->SetBranchAddress("scalhel.lhrs.fadc.pat", &scalhel_lhrs_fadc_pat, &b_scalhel_lhrs_fadc_pat);
   fChain->SetBranchAddress("scalhel.lhrs.fadc.tsettle", &scalhel_lhrs_fadc_tsettle, &b_scalhel_lhrs_fadc_tsettle);
   fChain->SetBranchAddress("scalhel.patcount", &scalhel_patcount, &b_scalhel_patcount);
   fChain->SetBranchAddress("scalhel.patphase", &scalhel_patphase, &b_scalhel_patphase);
   fChain->SetBranchAddress("scalhel.seed", &scalhel_seed, &b_scalhel_seed);
   fChain->SetBranchAddress("singletrack", &singletrack, &b_singletrack);
   fChain->SetBranchAddress("anytrack", &anytrack, &b_anytrack);
//    fChain->SetBranchAddress("singletrack", &singletrack, &b_singletrack);
//    fChain->SetBranchAddress("anytrack", &anytrack, &b_anytrack);
   fChain->SetBranchAddress("goodheep", &goodheep, &b_goodheep);
   fChain->SetBranchAddress("gooddxdy", &gooddxdy, &b_gooddxdy);
   fChain->SetBranchAddress("gooddpp", &gooddpp, &b_gooddpp);
   fChain->SetBranchAddress("goodctime", &goodctime, &b_goodctime);
   fChain->SetBranchAddress("HALLA_p", &HALLA_p, &b_HALLA_p);
   fChain->SetBranchAddress("hac_bcm_average", &hac_bcm_average, &b_hac_bcm_average);
   fChain->SetBranchAddress("IPM1H04A.XPOS", &IPM1H04A_XPOS, &b_IPM1H04A_XPOS);
   fChain->SetBranchAddress("IPM1H04A.YPOS", &IPM1H04A_YPOS, &b_IPM1H04A_YPOS);
   fChain->SetBranchAddress("IPM1H04E.XPOS", &IPM1H04E_XPOS, &b_IPM1H04E_XPOS);
   fChain->SetBranchAddress("IPM1H04E.YPOS", &IPM1H04E_YPOS, &b_IPM1H04E_YPOS);
   fChain->SetBranchAddress("IGL1I00OD16_16", &IGL1I00OD16_16, &b_IGL1I00OD16_16);
   fChain->SetBranchAddress("fEvtHdr.fEvtTime", &fEvtHdr_fEvtTime, &b_Event_Branch_fEvtHdr_fEvtTime);
   fChain->SetBranchAddress("fEvtHdr.fEvtNum", &fEvtHdr_fEvtNum, &b_Event_Branch_fEvtHdr_fEvtNum);
   fChain->SetBranchAddress("fEvtHdr.fEvtType", &fEvtHdr_fEvtType, &b_Event_Branch_fEvtHdr_fEvtType);
   fChain->SetBranchAddress("fEvtHdr.fEvtLen", &fEvtHdr_fEvtLen, &b_Event_Branch_fEvtHdr_fEvtLen);
   fChain->SetBranchAddress("fEvtHdr.fHelicity", &fEvtHdr_fHelicity, &b_Event_Branch_fEvtHdr_fHelicity);
   fChain->SetBranchAddress("fEvtHdr.fTrigBits", &fEvtHdr_fTrigBits, &b_Event_Branch_fEvtHdr_fTrigBits);
   fChain->SetBranchAddress("fEvtHdr.fRun", &fEvtHdr_fRun, &b_Event_Branch_fEvtHdr_fRun);
   Notify();
}

Bool_t gep_tree_data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gep_tree_data::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gep_tree_data::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gep_tree_data_cxx
