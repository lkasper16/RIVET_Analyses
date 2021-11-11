// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
//#include "RHICCentrality.hh" //external header for Centrality calculation
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES

namespace Rivet {

  class CMS_2020_I064906 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2020_I064906);

    /// Book histograms and initialise projections before the run
    void init() {

      // ----Initialise and register projections----
      // The basic final-state projection:
      // all final-state particles within the given eta acceptance
      // Particles: K0S, Lambda, Xi-, Anti-Xi+, Omega-, Anti-Omega+ ((((Anti Lambdas?))))
    	std::initializer_list<int> pdgIds = {310, 3122, 3312, -3312, 3334, -3334};
    	const PrimaryParticles fs(pdgIds, Cuts::abscharge > 0 && Cuts::absrap < 1.8);
    	declare(fs, "fs");
    	const PrimaryParticles ns(pdgIds, Cuts::abscharge == 0 && Cuts::absrap < 1.8);
    	declare(ns, "ns");
    	
    	beamOpt = getOption<string>("beam", "NONE");

	if (beamOpt == "PP") collSys = pp;
	else if (beamOpt == "pPB") collSys = pPB;
	
	//Create various counters
	book(sow["sow_pp"], "sow_pp");
	book(sow["sow_pPB"], "sow_pPB");

        //--------Begin booking histograms-----------
    	
    	//Inv. pT spectra of K0s for various y_CM (fig 2.1)
	book(hInvariantPTK0S["pp_-1.8<yCM<0"], 1, 1, 1);
	book(hInvariantPTK0S["pp_-1.8<yCM<1.8"], 1, 1, 2);
	book(hInvariantPTK0S["pp_0<yCM<1.8"], 1, 1, 3);
	book(hInvariantPTK0S["pPB_-1.8<yCM<0"], 1, 1, 4);
	book(hInvariantPTK0S["pPB_-1.8<yCM<1.8"], 1, 1, 5);
	book(hInvariantPTK0S["pPB_0<yCM<1.8"], 1, 1, 6);
	
	//Inv. pT spectra of Lambda for various y_CM (fig 2.2)
	book(hInvariantPTLambda["pp_-1.8<yCM<0"], 2, 1, 1);
	book(hInvariantPTLambda["pp_-1.8<yCM<1.8"], 2, 1, 2);
	book(hInvariantPTLambda["pp_0<yCM<1.8"], 2, 1, 3);
	book(hInvariantPTLambda["pPB_-1.8<yCM<0"], 2, 1, 4);
	book(hInvariantPTLambda["pPB_-1.8<yCM<1.8"], 2, 1, 5);
	book(hInvariantPTLambda["pPB_0<yCM<1.8"], 2, 1, 6);
	
	//Inv. pT spectra of Xi for various y_CM (fig 2.3)
	book(hInvariantPTXi["pp_-1.8<yCM<0"], 3, 1, 1);
	book(hInvariantPTXi["pp_-1.8<yCM<1.8"], 3, 1, 2);
	book(hInvariantPTXi["pp_0<yCM<1.8"], 3, 1, 3);
	book(hInvariantPTXi["pPB_-1.8<yCM<0"], 3, 1, 4);
	book(hInvariantPTXi["pPB_-1.8<yCM<1.8"], 3, 1, 5);
	book(hInvariantPTXi["pPB_0<yCM<1.8"], 3, 1, 6);
	
	//Inv. pT spectra of Omega for various y_CM (fig 2.4)
	book(hInvariantPTOmega["pp_-1.8<yCM<1.8"], 4, 1, 1);
	book(hInvariantPTOmega["pPB_-1.8<yCM<1.8"], 4, 1, 2);
	
	//R_pPB for -1.8<y_CM<1.8 in pPB (fig 3)
	//book(RpPBFullK0s["RpPB_K0S_-1.8<yCM<1.8"], 5, 1, 1);
	//book(RpPBFullLambda["RpPB_Lambda_-1.8<yCM<1.8"], 6, 1, 1);
	//book(RpPBFullXi["RpPB_Xi_-1.8<yCM<1.8"], 7, 1, 1);
	//book(RpPBFullOmega["RpPB_Omega_-1.8<yCM<1.8"], 8, 1, 1);
	
	//R_pPB for -1.8<y_CM<0 in pPB (fig 4.1)
	//book(hRpPBLowK0s["RpPB_K0S_-1.8<yCM<0"], 9, 1, 1);
	//book(hRpPBLowLambda["RpPB_Lambda_-1.8<yCM<0"], 10, 1, 1);
	//book(hRpPBLowXi["RpPB_Xi_-1.8<yCM<0"], 11, 1, 1);
	
	//R_pPB for 0<y_CM<1.8 in pPB (fig 4.2)
	//book(hRpPBHighK0s["RpPB_K0S_0<yCM<1.8"], 12, 1, 1);
	//book(hRpPBHighLambda["RpPB_Lambda_0<yCM<1.8"], 13, 1, 1);
	//book(hRpPBHighXi["RpPB_Xi_0<yCM<1.8"], 14, 1, 1);
	
	//Inv. pT of K0S for various y_CM in pPB (fig 5.1)
	book(hInvariantPTK0SpPB["-1.8<yCM<-1.3"], 15, 1, 1);
	book(hInvariantPTK0SpPB["-1.3<yCM<-0.8"], 15, 1, 2);
	book(hInvariantPTK0SpPB["-0.8<yCM<-0.3"], 15, 1, 3);
	book(hInvariantPTK0SpPB["0.3<yCM<0.8"], 15, 1, 4);
	book(hInvariantPTK0SpPB["0.8<yCM<1.3"], 15, 1, 5);
	book(hInvariantPTK0SpPB["0.3<yCM<1.8"], 15, 1, 6);
	
	//Inv. pT of Lambda for various y_CM in pPB (fig 5.2)
	book(hInvariantPTLambdapPB["-1.8<yCM<-1.3"], 16, 1, 1);
	book(hInvariantPTLambdapPB["-1.3<yCM<-0.8"], 16, 1, 2);
	book(hInvariantPTLambdapPB["-0.8<yCM<-0.3"], 16, 1, 3);
	book(hInvariantPTLambdapPB["0.3<yCM<0.8"], 16, 1, 4);
	book(hInvariantPTLambdapPB["0.8<yCM<1.3"], 16, 1, 5);
	book(hInvariantPTLambdapPB["1.3<yCM<1.8"], 16, 1, 6);
	
	//-----Book Scatter Plots from division-------
	
	//Fig 6.1, Y_asym for 0.3 < |yCM| < 0.8
	//Y_asym Low for K0s
	string refname1 = mkAxisCode(17, 1, 1);
	const Scatter2D& refdata1 = refData(refname1);
	book(hInvariantPTK0SpPB["LowNeg"], refname1 + "_LowNeg", refdata1);
	book(hInvariantPTK0SpPB["LowPos"], refname1 + "_LowPos", refdata1);
	book(YasymLowK0s["-0.8<yCM<-0.3/0.3<yCM<0.8"], refname1);
	
	//Y_asym Low for Lambda
	string refname2 = mkAxisCode(18, 1, 1);
	const Scatter2D& refdata2 = refData(refname2);
	book(hInvariantPTLambdapPB["LowNeg"], refname2 + "_LowNeg", refdata2);
	book(hInvariantPTLambdapPB["LowPos"], refname2 + "_LowPos", refdata2);
	book(YasymLowLambda["-0.8<yCM<-0.3/0.3<yCM<0.8"], refname2);
	
	//Y_asym Low for h+/-
	string refname3 = mkAxisCode(19, 1, 1);
	const Scatter2D& refdata3 = refData(refname3);
	book(h["negative_charged_yCM_low"], refname3 + "_LowNeg", refdata3);
	book(h["positive_charged_yCM_low"], refname3 + "_LowPos", refdata3);
	book(YasymLowh+/-["-0.8<yCM<-0.3/0.3<yCM<0.8"], refname3);
	
	//Fig 6.2, Y_asym for 0.8 < |yCM| < 1.3
	//Y_asym mid for K0s
	string refname4 = mkAxisCode(20, 1, 1);
	const Scatter2D& refdata4 = refData(refname4);
	book(hInvariantPTK0SpPB["MidNeg"], refname4 + "_-1.3<yCM<-0.8", refdata4);
	book(hInvariantPTK0SpPB["MidPos"], refname4 + "_0.8<yCM<1.3", refdata4);
	book(YasymMidK0s["-1.3<yCM<-0.8/0.8<yCM<1.3"], refname4);
	
	//Y_asym mid for Lambda
	string refname5 = mkAxisCode(21, 1, 1);
	const Scatter2D& refdata5 = refData(refname5);
	book(hInvariantPTLambdapPB["MidNeg"], refname5 + "_-1.3<yCM<-0.8", refdata5);
	book(hInvariantPTLambdapPB["MidPos"], refname5 + "_0.8<yCM<1.3", refdata5);
	book(YasymMidLambda["-1.3<yCM<-0.8/0.8<yCM<1.3"], refname5);
	
	//Y_asym mid for h+/-
	string refname6 = mkAxisCode(22, 1, 1);
	const Scatter2D& refdata6 = refData(refname6);
	book(h["negative_charged_yCM_mid"], refname6 + "_-1.3<yCM<-0.8", refdata6);
	book(h["positive_charged_yCM_mid"], refname6 + "_0.8<yCM<1.3", refdata6);
	book(YasymMidh+/-["-1.3<yCM<-0.8/0.8<yCM<1.3"], refname6);
	
	//Fig 6.2, Y_asym for 1.3 < |yCM| < 1.8
	//Y_asym high for K0s
	string refname7 = mkAxisCode(23, 1, 1);
	const Scatter2D& refdata7 = refData(refname7);
	book(hInvariantPTK0SpPB["HighNeg"], refname7 + "_-1.8<yCM<-1.3", refdata7);
	book(hInvariantPTK0SpPB["HighPos"], refname7 + "_1.3<yCM<1.8", refdata7);
	book(YasymHighK0s["-1.8<yCM<-1.3/1.3<yCM<1.8"], refname7);
	
	//Y_asym high for Lambda
	string refname8 = mkAxisCode(24, 1, 1);
	const Scatter2D& refdata8 = refData(refname8);
	book(hInvariantPTLambdapPB["HighNeg"], refname8 + "_-1.8<yCM<-1.3", refdata8);
	book(hInvariantPTLambdapPB["HighPos"], refname8 + "_1.3<yCM<1.8", refdata8);
	book(YasymHighLambda["-1.8<yCM<-1.3/1.3<yCM<1.8"], refname8);
	
	//Y_asym high for h+/-
	string refname9 = mkAxisCode(25, 1, 1);
	const Scatter2D& refdata9 = refData(refname9);
	book(h["negative_charged_yCM_high"], refname9 + "_-1.8<yCM<-1.3", refdata9);
	book(h["positive_charged_yCM_high"], refname9 + "_1.3<yCM<1.8", refdata9);
	book(YasymHighh+/-["-1.8<yCM<-1.3/1.3<yCM<1.8"], refname9);
	
	//Figure 3, R_pPB for -1.8<y_CM<1.8
	//R_pPB full for K0s
	string refname10 = mkAxisCode(5, 1, 1);
	const Scatter2D& refdata10 = refData(refname10);
	book(hInvariantPTK0S["pT_pPB_full"], refname10 + "_full_pPB", refdata10);
	book(hInvariantPTK0S["pT_pp_full"], refname10 + "_full_pp", refdata10);
	book(RpPBFullK0s["pPB/pp_-1.8<yCM<1.8"], refname10);
	
	//R_pPB full for Lambda
	string refname11 = mkAxisCode(6, 1, 1);
	const Scatter2D& refdata11 = refData(refname11);
	book(hInvariantPTLambda["pT_pPB_full"], refname11 + "_full_pPB", refdata11);
	book(hInvariantPTLambda["pT_pp_full"], refname11 + "_full_pp", refdata11);
	book(RpPBFullLambda["pPB/pp_-1.8<yCM<1.8"], refname11);
	
	//R_pPB full for Xi
	string refname12 = mkAxisCode(7, 1, 1);
	const Scatter2D& refdata12 = refData(refname12);
	book(hInvariantPTXi["pT_pPB_full"], refname12 + "_full_pPB", refdata12);
	book(hInvariantPTXi["pT_pp_full"], refname12 + "_full_pp", refdata12);
	book(RpPBFullXi["pPB/pp_-1.8<yCM<1.8"], refname12);
	
	//R_pPB full for Omega
	string refname13 = mkAxisCode(8, 1, 1);
	const Scatter2D& refdata13 = refData(refname13);
	book(hInvariantPTOmega["pT_pPB_full"], refname13 + "_full_pPB", refdata13);
	book(hInvariantPTOmega["pT_pp_full"], refname13 + "_full_pp", refdata13);
	book(RpPBFullOmega["pPB/pp_-1.8<yCM<1.8"], refname13);
	
	//Figure 4.1, R_pPB for -1.8<y_CM<0
	//R_pPB low for K0s
	string refname14 = mkAxisCode(9, 1, 1);
	const Scatter2D& refdata14 = refData(refname14);
	book(hInvariantPTK0S["pT_pPB_low"], refname14 + "_low_pPB", refdata14);
	book(hInvariantPTK0S["pT_pp_low"], refname14 + "_low_pp", refdata14);
	book(RpPBLowK0s["pPB/pp_-1.8<yCM<0"], refname14);
	
	//R_pPB low for Lambda
	string refname15 = mkAxisCode(10, 1, 1);
	const Scatter2D& refdata15 = refData(refname15);
	book(hInvariantPTLambda["pT_pPB_low"], refname15 + "_low_pPB", refdata15);
	book(hInvariantPTLambda["pT_pp_low"], refname15 + "_low_pp", refdata15);
	book(RpPBLowLambda["pPB/pp_-1.8<yCM<0"], refname15);
	
	//R_pPB low for Xi
	string refname16 = mkAxisCode(11, 1, 1);
	const Scatter2D& refdata16 = refData(refname16);
	book(hInvariantPTXi["pT_pPB_low"], refname16 + "_low_pPB", refdata16);
	book(hInvariantPTXi["pT_pp_low"], refname16 + "_low_pp", refdata16);
	book(RpPBLowXi["pPB/pp_-1.8<yCM<0"], refname16);
	
	//Figure 4.2, R_pPB for 0<y_CM<1.8
	//R_pPB high for K0s
	string refname17 = mkAxisCode(12, 1, 1);
	const Scatter2D& refdata17 = refData(refname17);
	book(hInvariantPTK0S["pT_pPB_high"], refname17 + "_high_pPB", refdata17);
	book(hInvariantPTK0S["pT_pp_high"], refname17 + "_high_pp", refdata17);
	book(RpPBHighK0s["pPB/pp_0<yCM<1.8"], refname17);
	
	//R_pPB high for Lambda
	string refname18 = mkAxisCode(13, 1, 1);
	const Scatter2D& refdata18 = refData(refname18);
	book(hInvariantPTLambda["pT_pPB_high"], refname18 + "_high_pPB", refdata18);
	book(hInvariantPTLambda["pT_pp_high"], refname18 + "_high_pp", refdata18);
	book(RpPBHighLambda["pPB/pp_0<yCM<1.8"], refname18);
	
	//R_pPB high for Xi
	string refname19 = mkAxisCode(14, 1, 1);
	const Scatter2D& refdata19 = refData(refname19);
	book(hInvariantPTXi["pT_pPB_high"], refname19 + "_high_pPB", refdata19);
	book(hInvariantPTXi["pT_pp_high"], refname19 + "_high_pp", refdata19);
	book(RpPBHighXi["pPB/pp_0<yCM<1.8"], refname19);
	

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
    

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h["XXXX"]); // normalize to unity
      normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
      scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}


  };


  DECLARE_RIVET_PLUGIN(CMS_2020_I064906);

}
