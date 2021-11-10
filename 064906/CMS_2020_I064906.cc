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
	book(hInvariantPTK0S["pT_K0S_pp_-1.8<yCM<0"], 1, 1, 1);
	book(hInvariantPTK0S["pT_K0S_pp_-1.8<yCM<1.8"], 1, 1, 2);
	book(hInvariantPTK0S["pT_K0S_pp_0<yCM<1.8"], 1, 1, 3);
	book(hInvariantPTK0S["pT_K0S_pPB_-1.8<yCM<0"], 1, 1, 4);
	book(hInvariantPTK0S["pT_K0S_pPB_-1.8<yCM<1.8"], 1, 1, 5);
	book(hInvariantPTK0S["pT_K0S_pPB_0<yCM<1.8"], 1, 1, 6);
	
	//Inv. pT spectra of Lambda for various y_CM (fig 2.2)
	book(hInvariantPTLambda["pT_Lambda_pp_-1.8<yCM<0"], 2, 1, 1);
	book(hInvariantPTLambda["pT_Lambda_pp_-1.8<yCM<1.8"], 2, 1, 2);
	book(hInvariantPTLambda["pT_Lambda_pp_0<yCM<1.8"], 2, 1, 3);
	book(hInvariantPTLambda["pT_Lambda_pPB_-1.8<yCM<0"], 2, 1, 4);
	book(hInvariantPTLambda["pT_Lambda_pPB_-1.8<yCM<1.8"], 2, 1, 5);
	book(hInvariantPTLambda["pT_Lambda_pPB_0<yCM<1.8"], 2, 1, 6);
	
	//Inv. pT spectra of Xi for various y_CM (fig 2.3)
	book(hInvariantPTXi["pT_Xi_pp_-1.8<yCM<0"], 3, 1, 1);
	book(hInvariantPTXi["pT_Xi_pp_-1.8<yCM<1.8"], 3, 1, 2);
	book(hInvariantPTXi["pT_Xi_pp_0<yCM<1.8"], 3, 1, 3);
	book(hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<0"], 3, 1, 4);
	book(hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<1.8"], 3, 1, 5);
	book(hInvariantPTXi["pT_Xi_pPB_0<yCM<1.8"], 3, 1, 6);
	
	//Inv. pT spectra of Omega for various y_CM (fig 2.4)
	book(hInvariantPTOmega["pT_Omega_pp_-1.8<yCM<1.8"], 4, 1, 1);
	book(hInvariantPTOmega["pT_Omega_pPB_-1.8<yCM<1.8"], 4, 1, 2);
	
	//R_pPB for -1.8<y_CM<1.8 in pPB (fig 3)
	book(RpPBFullK0s["RpPB_K0S_-1.8<yCM<1.8"], 5, 1, 1);
	book(RpPBFullLambda["RpPB_Lambda_-1.8<yCM<1.8"], 6, 1, 1);
	book(RpPBFullXi["RpPB_Xi_-1.8<yCM<1.8"], 7, 1, 1);
	book(RpPBFullOmega["RpPB_Omega_-1.8<yCM<1.8"], 8, 1, 1);
	
	//R_pPB for -1.8<y_CM<0 in pPB (fig 4.1)
	book(hRpPBLowK0s["RpPB_K0S_-1.8<yCM<0"], 9, 1, 1);
	book(hRpPBLowLambda["RpPB_Lambda_-1.8<yCM<0"], 10, 1, 1);
	book(hRpPBLowXi["RpPB_Xi_-1.8<yCM<0"], 11, 1, 1);
	
	//R_pPB for 0<y_CM<1.8 in pPB (fig 4.2)
	book(hRpPBHighK0s["RpPB_K0S_0<yCM<1.8"], 12, 1, 1);
	book(hRpPBHighLambda["RpPB_Lambda_0<yCM<1.8"], 13, 1, 1);
	book(hRpPBHighXi["RpPB_Xi_0<yCM<1.8"], 14, 1, 1);
	
	//Inv. pT of K0S for various y_CM in pPB (fig 5.1)
	book(hInvariantPTK0SpPB["pT_K0S_pPB_-1.8<yCM<-1.3"], 15, 1, 1);
	book(hInvariantPTK0SpPB["pT_K0S_pPB_-1.3<yCM<-0.8"], 15, 1, 2);
	book(hInvariantPTK0SpPB["pT_K0S_pPB_-0.8<yCM<-0.3"], 15, 1, 3);
	book(hInvariantPTK0SpPB["pT_K0S_pPB_0.3<yCM<0.8"], 15, 1, 4);
	book(hInvariantPTK0SpPB["pT_K0S_pPB_0.8<yCM<1.3"], 15, 1, 5);
	book(hInvariantPTK0SpPB["pT_K0S_pPB_1.3<yCM<1.8"], 15, 1, 6);
	
	//Inv. pT of Lambda for various y_CM in pPB (fig 5.2)
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_-1.8<yCM<-1.3"], 16, 1, 1);
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_-1.3<yCM<-0.8"], 16, 1, 2);
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_-0.8<yCM<-0.3"], 16, 1, 3);
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_0.3<yCM<0.8"], 16, 1, 4);
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_0.8<yCM<1.3"], 16, 1, 5);
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_1.3<yCM<1.8"], 16, 1, 6);
	
	//-----Book Scatter Plots from division-------
	
	//Fig 6.1, Y_asym for 0.3 < |yCM| < 0.8
	//Y_asym Low for K0s
	string refname1 = mkAxisCode(17, 1, 1);
	const Scatter2D& refdata1 = refData(refname1);
	book(hInvariantPTK0SpPB["K0S_pT_pPB_LowNeg"], refname1 + "_LowNeg", refdata1);
	book(hInvariantPTK0SpPB["K0S_pT_pPB_LowPos"], refname1 + "_LowPos", refdata1);
	book(YasymLowK0s["-0.8<yCM<-0.3/0.3<yCM<0.8"], refname1);
	
	//Y_asym Low for Lambda
	string refname2 = mkAxisCode(18, 1, 1);
	const Scatter2D& refdata2 = refData(refname2);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_LowNeg"], refname2 + "_LowNeg", refdata2);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_LowPos"], refname2 + "_LowPos", refdata2);
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
	book(hInvariantPTK0SpPB["K0S_pT_pPB_MidNeg"], refname4 + "_-1.3<yCM<-0.8", refdata4);
	book(hInvariantPTK0SpPB["K0S_pT_pPB_MidPos"], refname4 + "_0.8<yCM<1.3", refdata4);
	book(YasymMidK0s["-1.3<yCM<-0.8/0.8<yCM<1.3"], refname4);
	
	//Y_asym mid for Lambda
	string refname5 = mkAxisCode(21, 1, 1);
	const Scatter2D& refdata5 = refData(refname5);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_MidNeg"], refname5 + "_-1.3<yCM<-0.8", refdata5);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_MidPos"], refname5 + "_0.8<yCM<1.3", refdata5);
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
	book(hInvariantPTK0SpPB["K0S_pT_pPB_HighNeg"], refname7 + "_-1.8<yCM<-1.3", refdata7);
	book(hInvariantPTK0SpPB["K0S_pT_pPB_HighPos"], refname7 + "_1.3<yCM<1.8", refdata7);
	book(YasymHighK0s["-1.8<yCM<-1.3/1.3<yCM<1.8"], refname7);
	
	//Y_asym high for Lambda
	string refname8 = mkAxisCode(24, 1, 1);
	const Scatter2D& refdata8 = refData(refname8);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_HighNeg"], refname8 + "_-1.8<yCM<-1.3", refdata8);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_HighPos"], refname8 + "_1.3<yCM<1.8", refdata8);
	book(YasymHighLambda["-1.8<yCM<-1.3/1.3<yCM<1.8"], refname8);
	
	//Y_asym high for h+/-
	string refname9 = mkAxisCode(25, 1, 1);
	const Scatter2D& refdata9 = refData(refname9);
	book(h["negative_charged_yCM_high"], refname9 + "_-1.8<yCM<-1.3", refdata9);
	book(h["positive_charged_yCM_high"], refname9 + "_1.3<yCM<1.8", refdata9);
	book(YasymHighh+/-["-1.8<yCM<-1.3/1.3<yCM<1.8"], refname9);
	
	//------Divisions for Figures 3 and 4-------------

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
