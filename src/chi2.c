#include "include.h"

//#define DEBUG

/*--------------------------------------------------------------------*/

void set_nuisance(struct nuisance* nuisparam)
/* Sets the structure nuisparam with the default values of its parameters */
{
	nuisparam->alphas_MZ=(struct indnuis){0.1180,0.0009,1,"alphas_MZ"};

	nuisparam->mass_b=(struct indnuis){4.183,0.007,1,"mass_b"};
	nuisparam->mass_c=(struct indnuis){1.273,0.0046,1,"mass_c"};
	nuisparam->mass_s=(struct indnuis){0.0935,0.008,1,"mass_s"};
	nuisparam->mass_top_pole=(struct indnuis){172.57,0.29,1,"mass_top_pole"};

	nuisparam->mass_h0=(struct indnuis){125.20,0.11,1,"mass_h0"};

	nuisparam->CKM_lambda=(struct indnuis){0.22501,0.00068,1,"CKM_lambda"};
	nuisparam->CKM_A=(struct indnuis){0.826,0.016,1,"CKM_A"};
	nuisparam->CKM_rhobar=(struct indnuis){0.1591,0.0094,1,"CKM_rhobar"};
	nuisparam->CKM_etabar=(struct indnuis){0.3523,0.0073,1,"CKM_etabar"};

	nuisparam->log_mu_W_mass_W=(struct indnuis){0.,1.,2,"log_mu_W_mass_W"}; /* log here is log2 */
	nuisparam->log_mu_b_mass_b=(struct indnuis){0.,1.,2,"log_mu_b_mass_b"}; /* log here is log2 */

/* inclusive b -> s */
	nuisparam->BR_BXclnu_exp=(struct indnuis){0.1065,0.0016,1,"BR_BXclnu_exp"};
/* b -> s gamma */ 
	nuisparam->mu_G2_bsg=(struct indnuis){0.336,0.064,1,"mu_G2_bsg"};
	nuisparam->rho_D3_bsg=(struct indnuis){0.153,0.045,1,"rho_D3_bsg"};
	nuisparam->rho_LS3_bsg=(struct indnuis){-0.145,0.098,1,"rho_LS3_bsg"};
	nuisparam->bsgamma_rand=(struct indnuis){0.,0.04,1,"bsgamma_rand"};
	nuisparam->mu_c_bsg=(struct indnuis){2.45,1.55,2,"mu_c_bsg"};
/* b -> s mu mu */
	nuisparam->BRBXsmumu_lowq2_rand=(struct indnuis){0.,0.05,1,"BRBXsmumu_lowq2_rand"};
	nuisparam->BRBXsmumu_highq2_rand=(struct indnuis){0.,0.05,1,"BRBXsmumu_highq2_rand"};
	nuisparam->BRBXsmumu_full_rand=(struct indnuis){0.,0.05,1,"BRBXsmumu_full_rand"};
/* b -> s e e */
	nuisparam->BRBXsee_lowq2_rand=(struct indnuis){0.,0.05,1,"BRBXsee_lowq2_rand"};
	nuisparam->BRBXsee_highq2_rand=(struct indnuis){0.,0.05,1,"BRBXsee_highq2_rand"};
	nuisparam->BRBXsee_full_rand=(struct indnuis){0.,0.05,1,"BRBXsee_full_rand"};
/* b -> s tau tau */
	nuisparam->BRBXstautau_lowq2_rand=(struct indnuis){0.,0.05,1,"BRBXstautau_lowq2_rand"};
	nuisparam->BRBXstautau_highq2_rand=(struct indnuis){0.,0.05,1,"BRBXstautau_highq2_rand"};
	nuisparam->BRBXstautau_full_rand=(struct indnuis){0.,0.05,1,"BRBXstautau_full_rand"};

/* B */	
	nuisparam->f_B=(struct indnuis){0.1900,0.0013,1,"f_B"};
	nuisparam->lambda_Bp=(struct indnuis){0.46,0.11,1,"lambda_Bp"};
/* B -> K* */
	nuisparam->f_Kstar_par=(struct indnuis){0.204,0.007,1,"f_Kstar_par"};
	nuisparam->f_Kstar_perp=(struct indnuis){0.159,0.006,1,"f_Kstar_perp"};
	nuisparam->a1perp=(struct indnuis){0.04,0.03,1,"a1perp"};
	nuisparam->a2perp=(struct indnuis){0.10,0.08,1,"a2perp"};
	nuisparam->a1par=(struct indnuis){0.06,0.04,1,"a1par"};
	nuisparam->a2par=(struct indnuis){0.16,0.09,1,"a2par"};
/* B -> K* gamma */
	nuisparam->T1_BKstar=(struct indnuis){0.312,0.027,1,"T1_BKstar"};
	nuisparam->log_mu_spec_lambda_h_mass_b=(struct indnuis){0.,1.,2,"log_mu_spec_lambda_h_mass_b"}; /* log here is log2 */
	/* low */
	nuisparam->BtoKstarlow_ALperp_err_noq2=(struct indnuis){0.,0.1,1,"BtoKstarlow_ALperp_err_noq2"};
	nuisparam->BtoKstarlow_ARperp_err_noq2=(struct indnuis){0.,0.1,1,"BtoKstarlow_ARperp_err_noq2"};
	nuisparam->BtoKstarlow_ALpar_err_noq2=(struct indnuis){0.,0.1,1,"BtoKstarlow_ALpar_err_noq2"};
	nuisparam->BtoKstarlow_ARpar_err_noq2=(struct indnuis){0.,0.1,1,"BtoKstarlow_ARpar_err_noq2"};
	nuisparam->BtoKstarlow_AL0_err_noq2=(struct indnuis){0.,0.1,1,"BtoKstarlow_AL0_err_noq2"};
	nuisparam->BtoKstarlow_AR0_err_noq2=(struct indnuis){0.,0.1,1,"BtoKstarlow_AR0_err_noq2"};
	nuisparam->BtoKstarlow_At_err_noq2=(struct indnuis){0.,0.1,1,"BtoKstarlow_At_err_noq2"};
	nuisparam->BtoKstarlow_AS_err_noq2=(struct indnuis){0.,0.1,1,"BtoKstarlow_AS_err_noq2"};

	nuisparam->BtoKstarlow_ALperp_err_q2=(struct indnuis){0.,0.25,1,"BtoKstarlow_ALperp_err_q2"};
	nuisparam->BtoKstarlow_ARperp_err_q2=(struct indnuis){0.,0.25,1,"BtoKstarlow_ARperp_err_q2"};
	nuisparam->BtoKstarlow_ALpar_err_q2=(struct indnuis){0.,0.25,1,"BtoKstarlow_ALpar_err_q2"};
	nuisparam->BtoKstarlow_ARpar_err_q2=(struct indnuis){0.,0.25,1,"BtoKstarlow_ARpar_err_q2"};
	nuisparam->BtoKstarlow_AL0_err_q2=(struct indnuis){0.,0.25,1,"BtoKstarlow_AL0_err_q2"};
	nuisparam->BtoKstarlow_AR0_err_q2=(struct indnuis){0.,0.25,1,"BtoKstarlow_AR0_err_q2"};
	nuisparam->BtoKstarlow_At_err_q2=(struct indnuis){0.,0.25,1,"BtoKstarlow_At_err_q2"};
	nuisparam->BtoKstarlow_AS_err_q2=(struct indnuis){0.,0.25,1,"BtoKstarlow_AS_err_q2"};

	/* SN: Result from B->K* Hadronic Fit (to be used to check how it impacts predictions) - 13/04/2025
	 * the default we put to zero*/
	nuisparam->hplus0r=(struct indnuis){0.0, 0.0, 1,"hplus0r"};
	nuisparam->hminus0r=(struct indnuis){0.0, 0.0, 1,"hminus0r"};
	nuisparam->hplus1r=(struct indnuis){0.0, 0.0, 1,"hplus1r"};
	nuisparam->hminus1r=(struct indnuis){0.0, 0.0, 1,"hminus1r"};
	nuisparam->hplus2r=(struct indnuis){0.0, 0.0, 1,"hplus2r"};
	nuisparam->hminus2r=(struct indnuis){0.0, 0.0, 1,"hminus2r"};
	nuisparam->hzero0r=(struct indnuis){0.0, 0.0, 1,"hzero0r"};
	nuisparam->hzero1r=(struct indnuis){0.0, 0.0, 1,"hzero1r"};
	nuisparam->hzero2r=(struct indnuis){0.0, 0.0, 1,"hzero2r"};
	nuisparam->hplus0i=(struct indnuis){0.0, 0.0, 1,"hplus0i"};
	nuisparam->hminus0i=(struct indnuis){0.0, 0.0, 1,"hminus0i"};
	nuisparam->hplus1i=(struct indnuis){0.0, 0.0, 1,"hplus1i"};
	nuisparam->hminus1i=(struct indnuis){0.0, 0.0, 1,"hminus1i"};
	nuisparam->hplus2i=(struct indnuis){0.0, 0.0, 1,"hplus2i"};
	nuisparam->hminus2i=(struct indnuis){0.0, 0.0, 1,"hminus2i"};
	nuisparam->hzero0i=(struct indnuis){0.0, 0.0, 1,"hzero0i"};
	nuisparam->hzero1i=(struct indnuis){0.0, 0.0, 1,"hzero1i"};
	nuisparam->hzero2i=(struct indnuis){0.0, 0.0, 1,"hzero2i"};

	nuisparam->real_alpha_perp0=(struct indnuis){-0.06e-4,0.21e-4,1,"real_alpha_perp0"};
	nuisparam->real_alpha_perp1=(struct indnuis){-6.77e-4,0.27e-4,1,"real_alpha_perp1"};
	nuisparam->real_alpha_perp2=(struct indnuis){18.96e-4,0.59e-4,1,"real_alpha_perp2"};
	nuisparam->real_alpha_par0=(struct indnuis){-0.35e-4,0.62e-4,1,"real_alpha_par0"};
	nuisparam->real_alpha_par1=(struct indnuis){-3.13e-4,0.41e-4,1,"real_alpha_par1"};
	nuisparam->real_alpha_par2=(struct indnuis){12.20e-4,1.34e-4,1,"real_alpha_par2"};
	nuisparam->real_alpha_zero0=(struct indnuis){0.05e-4,1.52e-4,1,"real_alpha_zero0"};
	nuisparam->real_alpha_zero1=(struct indnuis){17.26e-4,1.64e-4,1,"real_alpha_zero1"};
	nuisparam->imag_alpha_perp0=(struct indnuis){-0.21e-4,2.25e-4,1,"imag_alpha_perp0"};
	nuisparam->imag_alpha_perp1=(struct indnuis){1.17e-4,3.58e-4,1,"imag_alpha_perp1"};
	nuisparam->imag_alpha_perp2=(struct indnuis){-0.08e-4,2.24e-4,1,"imag_alpha_perp2"};
	nuisparam->imag_alpha_par0=(struct indnuis){-0.04e-4,3.67e-4,1,"imag_alpha_par0"};
	nuisparam->imag_alpha_par1=(struct indnuis){-2.14e-4,2.46e-4,1,"imag_alpha_par1"};
	nuisparam->imag_alpha_par2=(struct indnuis){6.03e-4,2.50e-4,1,"imag_alpha_par2"};
	nuisparam->imag_alpha_zero0=(struct indnuis){-0.05e-4,4.99e-4,1,"imag_alpha_zero0"};
	nuisparam->imag_alpha_zero1=(struct indnuis){4.29e-4,3.14e-4,1,"imag_alpha_zero1"};
	
	nuisparam->DeltaC9_M1_q2bar=(struct indnuis){0.72,0.47,1,"DeltaC9_M1_q2bar"};
	nuisparam->r1_M1=(struct indnuis){0.10,0.01,1,"r1_M1"};
	nuisparam->r2_M1=(struct indnuis){1.13,0.01,1,"r2_M1"};
	nuisparam->DeltaC9_M2_q2bar=(struct indnuis){0.76,0.56,1,"DeltaC9_M2_q2bar"};
	nuisparam->r1_M2=(struct indnuis){0.09,0.01,1,"r1_M2"};
	nuisparam->r2_M2=(struct indnuis){1.12,0.01,1,"r2_M2"};
	nuisparam->DeltaC9_M3_q2bar=(struct indnuis){1.11,0.92,1,"DeltaC9_M3_q2bar"};
	nuisparam->r1_M3=(struct indnuis){0.06,0.07,1,"r1_M3"};
	nuisparam->r2_M3=(struct indnuis){1.05,0.05,1,"r2_M3"};
	
	/* high */
	nuisparam->BtoKstarhigh_ALperp_err=(struct indnuis){0.,0.1,1,"BtoKstarhigh_ALperp_err"};
	nuisparam->BtoKstarhigh_ARperp_err=(struct indnuis){0.,0.1,1,"BtoKstarhigh_ARperp_err"};
	nuisparam->BtoKstarhigh_ALpar_err=(struct indnuis){0.,0.1,1,"BtoKstarhigh_ALpar_err"};
	nuisparam->BtoKstarhigh_ARpar_err=(struct indnuis){0.,0.1,1,"BtoKstarhigh_ARpar_err"};
	nuisparam->BtoKstarhigh_AL0_err=(struct indnuis){0.,0.1,1,"BtoKstarhigh_AL0_err"};
	nuisparam->BtoKstarhigh_AR0_err=(struct indnuis){0.,0.1,1,"BtoKstarhigh_AR0_err"};
	nuisparam->BtoKstarhigh_At_err=(struct indnuis){0.,0.1,1,"BtoKstarhigh_At_err"};
	nuisparam->BtoKstarhigh_AS_err=(struct indnuis){0.,0.1,1,"BtoKstarhigh_AS_err"};

	/* Form factors B->K* ll */	
	nuisparam->a0A0_BKstar=(struct indnuis){0.369196,0.0289419,1,"a0A0_BKstar"};
	nuisparam->a1A0_BKstar=(struct indnuis){-1.36584,0.256849,1,"a1A0_BKstar"};
	nuisparam->a2A0_BKstar=(struct indnuis){0.128191,1.63438,1,"a2A0_BKstar"};
	nuisparam->a0A1_BKstar=(struct indnuis){0.29725,0.026356,1,"a0A1_BKstar"};
	nuisparam->a1A1_BKstar=(struct indnuis){0.392378,0.187894,1,"a1A1_BKstar"};
	nuisparam->a2A1_BKstar=(struct indnuis){1.18916,1.02531,1,"a2A1_BKstar"};
	nuisparam->a0A12_BKstar=(struct indnuis){0.265375,0.0208033,1,"a0A12_BKstar"};
	nuisparam->a1A12_BKstar=(struct indnuis){0.533638,0.128772,1,"a1A12_BKstar"};
	nuisparam->a2A12_BKstar=(struct indnuis){0.483166,0.656273,1,"a2A12_BKstar"};
	nuisparam->a0V_BKstar=(struct indnuis){0.376313,0.0332944,1,"a0V_BKstar"};
	nuisparam->a1V_BKstar=(struct indnuis){-1.16597,0.261268,1,"a1V_BKstar"};
	nuisparam->a2V_BKstar=(struct indnuis){2.42443,1.53102,1,"a2V_BKstar"};
	nuisparam->a0T1_BKstar=(struct indnuis){0.312055,0.027496,1,"a0T1_BKstar"};
	nuisparam->a1T1_BKstar=(struct indnuis){-1.00893,0.189575,1,"a1T1_BKstar"};
	nuisparam->a2T1_BKstar=(struct indnuis){1.5272,1.63965,1,"a2T1_BKstar"};
	nuisparam->a0T2_BKstar=(struct indnuis){0.312055,0.027496,1,"a0T2_BKstar"};
	nuisparam->a1T2_BKstar=(struct indnuis){0.496846,0.166219,1,"a1T2_BKstar"};
	nuisparam->a2T2_BKstar=(struct indnuis){1.61431,0.803783,1,"a2T2_BKstar"};
	nuisparam->a0T23_BKstar=(struct indnuis){0.667412,0.0633333,1,"a0T23_BKstar"};
	nuisparam->a1T23_BKstar=(struct indnuis){1.31812,0.222129,1,"a1T23_BKstar"};
	nuisparam->a2T23_BKstar=(struct indnuis){3.82334,2.20338,1,"a2T23_BKstar"};

	/* SN: B -> Kstar FF van Dyk et al. 2305.06301 (BKstar_FormFactor_choice:2)*/
	nuisparam->a0A0_BKstar_GRvDV_BSZ=(struct indnuis){0.34213768,0.02653931227443545,1, "a0A0_BKstar_GRvDV_BSZ"};
	nuisparam->a1A0_BKstar_GRvDV_BSZ=(struct indnuis){-1.14740932,0.19807532658057136,1, "a1A0_BKstar_GRvDV_BSZ"};
	nuisparam->a2A0_BKstar_GRvDV_BSZ=(struct indnuis){2.37275841,1.7807500891478287,1, "a2A0_BKstar_GRvDV_BSZ"};
	nuisparam->a0A1_BKstar_GRvDV_BSZ=(struct indnuis){0.28926565,0.019194694657639125,1, "a0A1_BKstar_GRvDV_BSZ"};
	nuisparam->a1A1_BKstar_GRvDV_BSZ=(struct indnuis){0.45994786,0.12940883779711493,1, "a1A1_BKstar_GRvDV_BSZ"};
	nuisparam->a2A1_BKstar_GRvDV_BSZ=(struct indnuis){1.22417351,0.741928818014235,1, "a2A1_BKstar_GRvDV_BSZ"};
	nuisparam->a1A12_BKstar_GRvDV_BSZ=(struct indnuis){0.5511594,0.3403667507263305,1, "a1A12_BKstar_GRvDV_BSZ"};
	nuisparam->a2A12_BKstar_GRvDV_BSZ=(struct indnuis){0.57813167,2.0588237272773013,1, "a2A12_BKstar_GRvDV_BSZ"};
	nuisparam->a0V_BKstar_GRvDV_BSZ=(struct indnuis){0.36280141,0.02869665212180682,1, "a0V_BKstar_GRvDV_BSZ"};
	nuisparam->a1V_BKstar_GRvDV_BSZ=(struct indnuis){-1.08931588,0.16983992934525144,1, "a1V_BKstar_GRvDV_BSZ"};
	nuisparam->a2V_BKstar_GRvDV_BSZ=(struct indnuis){2.80292123,2.0353314029906775,1, "a2V_BKstar_GRvDV_BSZ"};
	nuisparam->a0T1_BKstar_GRvDV_BSZ=(struct indnuis){0.31924221,0.020953170690852493,1, "a0T1_BKstar_GRvDV_BSZ"};
	nuisparam->a1T1_BKstar_GRvDV_BSZ=(struct indnuis){-0.96188597,0.13737650599720463,1, "a1T1_BKstar_GRvDV_BSZ"};
	nuisparam->a2T1_BKstar_GRvDV_BSZ=(struct indnuis){2.09305925,1.270197252398225,1, "a2T1_BKstar_GRvDV_BSZ"};
	nuisparam->a1T2_BKstar_GRvDV_BSZ=(struct indnuis){0.59323779,0.1802750043683261,1, "a1T2_BKstar_GRvDV_BSZ"};
	nuisparam->a2T2_BKstar_GRvDV_BSZ=(struct indnuis){1.6795425,0.9860532445055896,1, "a2T2_BKstar_GRvDV_BSZ"};
	nuisparam->a0T23_BKstar_GRvDV_BSZ=(struct indnuis){0.62123502,0.034463685380411654,1, "a0T23_BKstar_GRvDV_BSZ"};
	nuisparam->a1T23_BKstar_GRvDV_BSZ=(struct indnuis){0.97330702,0.3118172612284317,1, "a1T23_BKstar_GRvDV_BSZ"};
	nuisparam->a2T23_BKstar_GRvDV_BSZ=(struct indnuis){1.82128227,2.420460491311519,1, "a2T23_BKstar_GRvDV_BSZ"};

	/* SN: LCSR+Lattice B -> Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:3)*/
	nuisparam->a0A0_BKstar_GKvD_LCSR_Lattice=(struct indnuis){0.3433920997134539,0.0294990807136424,1,"a0A0_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a1A0_BKstar_GKvD_LCSR_Lattice=(struct indnuis){-1.1296312022549930,0.1991670817248993,1,"a1A0_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a2A0_BKstar_GKvD_LCSR_Lattice=(struct indnuis){2.3262919488225120,1.7880455795392056,1,"a2A0_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a0A1_BKstar_GKvD_LCSR_Lattice=(struct indnuis){0.3012383696993843,0.0267260683919918,1,"a0A1_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a1A1_BKstar_GKvD_LCSR_Lattice=(struct indnuis){0.5046280709625639,0.1495448534660424,1,"a1A1_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a2A1_BKstar_GKvD_LCSR_Lattice=(struct indnuis){1.1417500321900249,0.7525585788011253,1,"a2A1_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a1A12_BKstar_GKvD_LCSR_Lattice=(struct indnuis){0.5553962796553867,0.3482916361416120,1,"a1A12_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a2A12_BKstar_GKvD_LCSR_Lattice=(struct indnuis){0.6647053900993466,2.0808778959979359,1,"a2A12_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a0V_BKstar_GKvD_LCSR_Lattice=(struct indnuis){0.3815471168511042,0.0389344776806712,1,"a0V_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a1V_BKstar_GKvD_LCSR_Lattice=(struct indnuis){-1.1035111471855810,0.1697989659892044,1,"a1V_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a2V_BKstar_GKvD_LCSR_Lattice=(struct indnuis){2.1066350644810390,2.1513996843903724,1,"a2V_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a0T1_BKstar_GKvD_LCSR_Lattice=(struct indnuis){0.3378196251357424,0.0335821013401675,1,"a0T1_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a1T1_BKstar_GKvD_LCSR_Lattice=(struct indnuis){-0.8677845009935159,0.3895946174684902,1,"a1T1_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a2T1_BKstar_GKvD_LCSR_Lattice=(struct indnuis){-0.1169443466566312,5.8256054171228353,1,"a2T1_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a1T2_BKstar_GKvD_LCSR_Lattice=(struct indnuis){0.6959696346749498,0.3604650719791208,1,"a1T2_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a2T2_BKstar_GKvD_LCSR_Lattice=(struct indnuis){-1.1136735066833230,6.2808611737816298,1,"a2T2_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a0T23_BKstar_GKvD_LCSR_Lattice=(struct indnuis){0.6454874220120908,0.0518036855188885,1,"a0T23_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a1T23_BKstar_GKvD_LCSR_Lattice=(struct indnuis){0.7287741888333001,0.7420233695564566,1,"a1T23_BKstar_GKvD_LCSR_Lattice"};
	nuisparam->a2T23_BKstar_GKvD_LCSR_Lattice=(struct indnuis){-0.3112553766962693,13.2808707930017214,1,"a2T23_BKstar_GKvD_LCSR_Lattice"};


	/* SN: LCSR-only B -> Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:4)*/
	nuisparam->a0A0_BKstar_GKvD_LCSRonly=(struct indnuis){0.3052008490740308,0.0684995456216797,1,"a0A0_BKstar_GKvD_LCSRonly"};
	nuisparam->a1A0_BKstar_GKvD_LCSRonly=(struct indnuis){-0.9062138334213113,0.4558019393454409,1,"a1A0_BKstar_GKvD_LCSRonly"};
	nuisparam->a2A0_BKstar_GKvD_LCSRonly=(struct indnuis){0.6180529797336485,6.5937546287679067,1,"a2A0_BKstar_GKvD_LCSRonly"};
	nuisparam->a0A1_BKstar_GKvD_LCSRonly=(struct indnuis){0.2593945344444926,0.0818906969253259,1,"a0A1_BKstar_GKvD_LCSRonly"};
	nuisparam->a1A1_BKstar_GKvD_LCSRonly=(struct indnuis){0.4373081045870070,0.3552149777814105,1,"a1A1_BKstar_GKvD_LCSRonly"};
	nuisparam->a2A1_BKstar_GKvD_LCSRonly=(struct indnuis){-0.5746207172478975,4.9537305108459408,1,"a2A1_BKstar_GKvD_LCSRonly"};
	nuisparam->a1A12_BKstar_GKvD_LCSRonly=(struct indnuis){0.2562656111706637,0.8166863097948055,1,"a1A12_BKstar_GKvD_LCSRonly"};
	nuisparam->a2A12_BKstar_GKvD_LCSRonly=(struct indnuis){-0.5503409807129199,14.1653668308197052,1,"a2A12_BKstar_GKvD_LCSRonly"};
	nuisparam->a0V_BKstar_GKvD_LCSRonly=(struct indnuis){ 0.3245160186497879,0.1103948113715595,1,"a0V_BKstar_GKvD_LCSRonly"};
	nuisparam->a1V_BKstar_GKvD_LCSRonly=(struct indnuis){ -0.7759254346044123,0.5451351881898320,1,"a1V_BKstar_GKvD_LCSRonly"};
	nuisparam->a2V_BKstar_GKvD_LCSRonly=(struct indnuis){ -0.5615478365432042,6.7782283380918873,1,"a2V_BKstar_GKvD_LCSRonly"};
	nuisparam->a0T1_BKstar_GKvD_LCSRonly=(struct indnuis){0.2874948700463931,0.0939529049713691,1,"a0T1_BKstar_GKvD_LCSRonly"};
	nuisparam->a1T1_BKstar_GKvD_LCSRonly=(struct indnuis){-0.7209004764950822,0.4640113111771070,1,"a1T1_BKstar_GKvD_LCSRonly"};
	nuisparam->a2T1_BKstar_GKvD_LCSRonly=(struct indnuis){-0.2457768503570863,5.8235899428506981,1,"a2T1_BKstar_GKvD_LCSRonly"};
	nuisparam->a1T2_BKstar_GKvD_LCSRonly=(struct indnuis){0.5947831288310546,0.4026827682824656,1,"a1T2_BKstar_GKvD_LCSRonly"};
	nuisparam->a2T2_BKstar_GKvD_LCSRonly=(struct indnuis){-0.9282192907339308,6.2824253518395023,1,"a2T2_BKstar_GKvD_LCSRonly"};
	nuisparam->a0T23_BKstar_GKvD_LCSRonly=(struct indnuis){0.5798572662094974,0.1240371504520887,1,"a0T23_BKstar_GKvD_LCSRonly"};
	nuisparam->a1T23_BKstar_GKvD_LCSRonly=(struct indnuis){0.6536261917256043,0.7535129610959626,1,"a1T23_BKstar_GKvD_LCSRonly"};
	nuisparam->a2T23_BKstar_GKvD_LCSRonly=(struct indnuis){-0.4504317678393136,13.3064988735766967,1,"a2T23_BKstar_GKvD_LCSRonly"};

	/* SN: LCSR-only B->Kstar FF Zwicky et al. 1503.05534 (BKstar_FormFactor_choice:5)*/
	nuisparam->a0V_BKstar_BSZ_LCSRonly=(struct indnuis){0.341428,0.0360552,1, "a0V_BKstar_BSZ_LCSRonly"};
	nuisparam->a1V_BKstar_BSZ_LCSRonly=(struct indnuis){-1.04834,0.240855,1, "a1V_BKstar_BSZ_LCSRonly"};
	nuisparam->a2V_BKstar_BSZ_LCSRonly=(struct indnuis){2.37143,1.38825,1, "a2V_BKstar_BSZ_LCSRonly"};
	nuisparam->a0T1_BKstar_BSZ_LCSRonly=(struct indnuis){0.28235,0.0306646,1, "a0T1_BKstar_BSZ_LCSRonly"};
	nuisparam->a1T1_BKstar_BSZ_LCSRonly=(struct indnuis){-0.888396,0.185458,1, "a1T1_BKstar_BSZ_LCSRonly"};
	nuisparam->a2T1_BKstar_BSZ_LCSRonly=(struct indnuis){1.94823,1.09865,1, "a2T1_BKstar_BSZ_LCSRonly"};
	nuisparam->a0T2_BKstar_BSZ_LCSRonly=(struct indnuis){0.28235,0.0306646,1, "a0T2_BKstar_BSZ_LCSRonly"};
	nuisparam->a1T2_BKstar_BSZ_LCSRonly=(struct indnuis){0.398974,0.183893,1, "a1T2_BKstar_BSZ_LCSRonly"};
	nuisparam->a2T2_BKstar_BSZ_LCSRonly=(struct indnuis){0.36137,0.514046,1, "a2T2_BKstar_BSZ_LCSRonly"};
	nuisparam->a0T23_BKstar_BSZ_LCSRonly=(struct indnuis){0.667768,0.0831775,1, "a0T23_BKstar_BSZ_LCSRonly"};
	nuisparam->a1T23_BKstar_BSZ_LCSRonly=(struct indnuis){1.47676,0.49123,1, "a1T23_BKstar_BSZ_LCSRonly"};
	nuisparam->a2T23_BKstar_BSZ_LCSRonly=(struct indnuis){1.92352,1.95757,1, "a2T23_BKstar_BSZ_LCSRonly"};
	nuisparam->a0A0_BKstar_BSZ_LCSRonly=(struct indnuis){0.355851,0.0462313,1, "a0A0_BKstar_BSZ_LCSRonly"};
	nuisparam->a1A0_BKstar_BSZ_LCSRonly=(struct indnuis){-1.04363,0.268857,1, "a1A0_BKstar_BSZ_LCSRonly"};
	nuisparam->a2A0_BKstar_BSZ_LCSRonly=(struct indnuis){1.12403,1.35295,1, "a2A0_BKstar_BSZ_LCSRonly"};
	nuisparam->a0A1_BKstar_BSZ_LCSRonly=(struct indnuis){0.269264,0.0293981,1, "a0A1_BKstar_BSZ_LCSRonly"};
	nuisparam->a1A1_BKstar_BSZ_LCSRonly=(struct indnuis){0.304578,0.189926,1, "a1A1_BKstar_BSZ_LCSRonly"};
	nuisparam->a2A1_BKstar_BSZ_LCSRonly=(struct indnuis){-0.10662,0.479212,1, "a2A1_BKstar_BSZ_LCSRonly"};
	nuisparam->a0A12_BKstar_BSZ_LCSRonly=(struct indnuis){0.255783,0.0332308,1, "a0A12_BKstar_BSZ_LCSRonly"};
	nuisparam->a1A12_BKstar_BSZ_LCSRonly=(struct indnuis){0.601902,0.203756,1, "a1A12_BKstar_BSZ_LCSRonly"};
	nuisparam->a2A12_BKstar_BSZ_LCSRonly=(struct indnuis){0.117626,0.839054,1, "a2A12_BKstar_BSZ_LCSRonly"};

	/* SN: Lattice-only B->Kstar FF Wingate et al. 1310.3722 + 1501.00367 (BKstar_FormFactor_choice:6*/
	nuisparam->a0V_BKstar_HLMW=(struct indnuis){0.497528420813,0.0666823829052322,1, "a0V_BKstar_HLMW"};
	nuisparam->a1V_BKstar_HLMW=(struct indnuis){-2.0150648425,0.916490839633981,1, "a1V_BKstar_HLMW"};
	nuisparam->a0A0_BKstar_HLMW=(struct indnuis){0.502302643141,0.0370063274557203,1, "a0A0_BKstar_HLMW"};
	nuisparam->a1A0_BKstar_HLMW=(struct indnuis){-1.60847600112,0.447220134489716,1, "a1A0_BKstar_HLMW"};
	nuisparam->a0A1_BKstar_HLMW=(struct indnuis){0.284803612469,0.0233211047821067,1, "a0A1_BKstar_HLMW"};
	nuisparam->a1A1_BKstar_HLMW=(struct indnuis){0.191411584999,0.280445814343876,1, "a1A1_BKstar_HLMW"};
	nuisparam->a0A12_BKstar_HLMW=(struct indnuis){0.219527522577,0.0237567768756622,1, "a0A12_BKstar_HLMW"};
	nuisparam->a1A12_BKstar_HLMW=(struct indnuis){0.332374417965,0.30017771284774,1, "a1A12_BKstar_HLMW"};
	nuisparam->a0T1_BKstar_HLMW=(struct indnuis){0.419712711234,0.0240583451819738,1, "a0T1_BKstar_HLMW"};
	nuisparam->a1T1_BKstar_HLMW=(struct indnuis){-1.36328797303,0.258574427178714,1, "a1T1_BKstar_HLMW"};
	nuisparam->a0T2_BKstar_HLMW=(struct indnuis){0.279965680646,0.0194805665334456,1, "a0T2_BKstar_HLMW"};
	nuisparam->a1T2_BKstar_HLMW=(struct indnuis){0.11706607999,0.236379782057603,1, "a1T2_BKstar_HLMW"};
	nuisparam->a0T23_BKstar_HLMW=(struct indnuis){0.523462757534,0.0451285804023348,1, "a0T23_BKstar_HLMW"};
	nuisparam->a1T23_BKstar_HLMW=(struct indnuis){-0.27135970224,0.579122648412234,1, "a1T23_BKstar_HLMW"};
	nuisparam->systErr_BKstar_HLMW=(struct indnuis){0.,0.05,1, "systErr_BKstar_HLMW"};


/* B -> K */
	nuisparam->f_K=(struct indnuis){0.156,0.005,1,"f_K"};
	nuisparam->a1K=(struct indnuis){0.06,0.03,1,"a1K"};
	nuisparam->a2K=(struct indnuis){0.25,0.15,1,"a2K"};
	
	/* Form factors B->K ll */
	/* SN: Commenting out, replacing with more explicit name, adding missing a3f0 - 28/03/2025 */
	// nuisparam->a00_BK=(struct indnuis){0.3290899285463681,0.02772821914511607,1,"a00_BK"};
	// nuisparam->a10_BK=(struct indnuis){0.195117199570376,0.1680967198344544,1,"a10_BK"};
	// nuisparam->a20_BK=(struct indnuis){-0.4461264576740431,0.40894628457968657,1,"a20_BK"};
	// nuisparam->a0p_BK=(struct indnuis){0.3290899285463681,0.02772821914511607,1,"a0p_BK"};
	// nuisparam->a1p_BK=(struct indnuis){-0.8669465867427438,0.13814361973628353,1,"a1p_BK"};
	// nuisparam->a2p_BK=(struct indnuis){0.006095669673341246,0.7508262983344038,1,"a2p_BK"};
	// nuisparam->a0T_BK=(struct indnuis){0.2993831492500477,0.025872030646127343,1,"a0T_BK"};
	// nuisparam->a1T_BK=(struct indnuis){-0.7735456824025474,0.1501283039139932,1,"a1T_BK"};
	// nuisparam->a2T_BK=(struct indnuis){0.009554375102511888,0.8715760989774238,1,"a2T_BK"};
	/* SN: B->K FF Altmannshofer, Straub 1411.3161 (BK_FormFactor_choice:1)*/
	nuisparam->a0f0_BK_AS_LCSR_Lattice=(struct indnuis){0.54,0.03,1,"a0f0_BK_AS_LCSR_Lattice"};
	nuisparam->a1f0_BK_AS_LCSR_Lattice=(struct indnuis){-1.91,0.10,1,"a1f0_BK_AS_LCSR_Lattice"};
	nuisparam->a2f0_BK_AS_LCSR_Lattice=(struct indnuis){1.83,1.07,1,"a2f0_BK_AS_LCSR_Lattice"};
	nuisparam->a3f0_BK_AS_LCSR_Lattice=(struct indnuis){-0.02,2.74,1,"a3f0_BK_AS_LCSR_Lattice"};
	nuisparam->a0fp_BK_AS_LCSR_Lattice=(struct indnuis){0.43,0.02,1,"a0fp_BK_AS_LCSR_Lattice"};
	nuisparam->a1fp_BK_AS_LCSR_Lattice=(struct indnuis){-0.67,0.09,1,"a1fp_BK_AS_LCSR_Lattice"};
	nuisparam->a2fp_BK_AS_LCSR_Lattice=(struct indnuis){-1.12,0.76,1,"a2fp_BK_AS_LCSR_Lattice"};
	nuisparam->a0fT_BK_AS_LCSR_Lattice=(struct indnuis){0.40,0.02,1,"a0fT_BK_AS_LCSR_Lattice"};
	nuisparam->a1fT_BK_AS_LCSR_Lattice=(struct indnuis){-0.53,0.13,1,"a1fT_BK_AS_LCSR_Lattice"};
	nuisparam->a2fT_BK_AS_LCSR_Lattice=(struct indnuis){-0.29,1.00,1,"a2fT_BK_AS_LCSR_Lattice"};

	/* SN: B->K FF van Dyk et al. 2305.06301 (BK_FormFactor_choice:2)*/
	nuisparam->a0fp_BK_GRvDV_BSZ=(struct indnuis){0.33196382,0.010306670170331444,1, "a0fp_BK_GRvDV_BSZ"};
	nuisparam->a1fp_BK_GRvDV_BSZ=(struct indnuis){-0.96403506,0.06670450636951,1, "a1fp_BK_GRvDV_BSZ"};
	nuisparam->a2fp_BK_GRvDV_BSZ=(struct indnuis){-0.30273752,0.08985931086982583,1, "a2fp_BK_GRvDV_BSZ"};
	nuisparam->a1f0_BK_GRvDV_BSZ=(struct indnuis){0.33478906,0.0775275266598903,1, "a1f0_BK_GRvDV_BSZ"};
	nuisparam->a2f0_BK_GRvDV_BSZ=(struct indnuis){0.00372445,0.14655033606239187,1, "a2f0_BK_GRvDV_BSZ"};
	nuisparam->a0fT_BK_GRvDV_BSZ=(struct indnuis){0.32361471,0.021538583959954283,1, "a0fT_BK_GRvDV_BSZ"};
	nuisparam->a1fT_BK_GRvDV_BSZ=(struct indnuis){-1.02911783,0.17444811348936967,1, "a1fT_BK_GRvDV_BSZ"};
	nuisparam->a2fT_BK_GRvDV_BSZ=(struct indnuis){-0.39954681,0.34388168314116413,1, "a2fT_BK_GRvDV_BSZ"};

	/* SN: LCSR+Lattice B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:3)*/
	nuisparam->a0fp_BK_GKvD_LCSR_Lattice=(struct indnuis){0.3290899285463681,0.0277282191451161,1,"a0fp_BK_GKvD_LCSR_Lattice"};
	nuisparam->a1fp_BK_GKvD_LCSR_Lattice=(struct indnuis){-0.8669465867427438,0.1381436197362835,1,"a1fp_BK_GKvD_LCSR_Lattice"};
	nuisparam->a2fp_BK_GKvD_LCSR_Lattice=(struct indnuis){0.0060956696733412,0.7508262983344038,1,"a2fp_BK_GKvD_LCSR_Lattice"};
	nuisparam->a1f0_BK_GKvD_LCSR_Lattice=(struct indnuis){0.1951171995703760,0.1680967198344544,1,"a1f0_BK_GKvD_LCSR_Lattice"};
	nuisparam->a2f0_BK_GKvD_LCSR_Lattice=(struct indnuis){-0.4461264576740431,0.4089462845796866,1,"a2f0_BK_GKvD_LCSR_Lattice"};
	nuisparam->a0fT_BK_GKvD_LCSR_Lattice=(struct indnuis){0.2993831492500477,0.0258720306461273,1,"a0fT_BK_GKvD_LCSR_Lattice"};
	nuisparam->a1fT_BK_GKvD_LCSR_Lattice=(struct indnuis){-0.7735456824025474,0.1501283039139932,1,"a1fT_BK_GKvD_LCSR_Lattice"};
	nuisparam->a2fT_BK_GKvD_LCSR_Lattice=(struct indnuis){0.0095543751025119,0.8715760989774238,1,"a2fT_BK_GKvD_LCSR_Lattice"};

	/* SN: LCSR-only B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:4)*/
	nuisparam->a0fp_BK_GKvD_LCSRonly=(struct indnuis){0.2655528728950212,0.0800990560984796,1,"a0fp_BK_GKvD_LCSRonly"};
	nuisparam->a1fp_BK_GKvD_LCSRonly=(struct indnuis){-0.6466140804171657,0.4882976078914857,1,"a1fp_BK_GKvD_LCSRonly"};
	nuisparam->a2fp_BK_GKvD_LCSRonly=(struct indnuis){-0.1337677825631754,5.7534987880040864,1,"a2fp_BK_GKvD_LCSRonly"};
	nuisparam->a1f0_BK_GKvD_LCSRonly=(struct indnuis){0.3841222758978433,0.3828140046234049,1,"a1f0_BK_GKvD_LCSRonly"};
	nuisparam->a2f0_BK_GKvD_LCSRonly=(struct indnuis){-0.6628825163091753,5.9955080919574426,1,"a2f0_BK_GKvD_LCSRonly"};
	nuisparam->a0fT_BK_GKvD_LCSRonly=(struct indnuis){0.2510192324158927,0.0675867314925144,1,"a0fT_BK_GKvD_LCSRonly"};
	nuisparam->a1fT_BK_GKvD_LCSRonly=(struct indnuis){-0.6508680050905388,0.4608724471636461,1,"a1fT_BK_GKvD_LCSRonly"};
	nuisparam->a2fT_BK_GKvD_LCSRonly=(struct indnuis){0.0999901466869552,5.6460204741706210,1,"a2fT_BK_GKvD_LCSRonly"};

	/* SN: Lattice FLAG24 average 2411.04268 of HPQCD'13 and FNAL/MILC'15 (BK_FormFactor_choice:6) - 26/04/2025*/
	nuisparam->a0fp_BK_FLAG24=(struct indnuis){0.471,0.014,1,"a0fp_BK_FLAG24"};
	nuisparam->a1fp_BK_FLAG24=(struct indnuis){-0.74,0.16,1,"a1fp_BK_FLAG24"};
	nuisparam->a2fp_BK_FLAG24=(struct indnuis){0.32,0.71,1,"a2fp_BK_FLAG24"};
	nuisparam->a0f0_BK_FLAG24=(struct indnuis){0.301,0.010,1,"a0f0_BK_FLAG24"};
	nuisparam->a1f0_BK_FLAG24=(struct indnuis){0.40,0.15,1,"a1f0_BK_FLAG24"};
	nuisparam->a0fT_BK_FLAG24=(struct indnuis){0.455,0.021,1,"a0fT_BK_FLAG24"};
	nuisparam->a1fT_BK_FLAG24=(struct indnuis){-1.00,0.31,1,"a1fT_BK_FLAG24"};
	nuisparam->a2fT_BK_FLAG24=(struct indnuis){-0.9,1.3,1,"a2fT_BK_FLAG24"};

	/* SN: Lattice full q2-region HPQCD22 2207.12468 (BK_FormFactor_choice:7) - 26/04/2025*/
	nuisparam->a0f0_BK_HPQCD22=(struct indnuis){0.2545,0.090,1,"a0f0_BK_HPQCD22"};
	nuisparam->a1f0_BK_HPQCD22=(struct indnuis){0.21,0.076,1,"a1f0_BK_HPQCD22"};
	nuisparam->a2f0_BK_HPQCD22=(struct indnuis){0.02,0.17,1,"a2f0_BK_HPQCD22"};
	nuisparam->a1fp_BK_HPQCD22=(struct indnuis){-0.71,0.14,1,"a1fp_BK_HPQCD22"};
	nuisparam->a2fp_BK_HPQCD22=(struct indnuis){0.32,0.59,1,"a2fp_BK_HPQCD22"};
	nuisparam->a0fT_BK_HPQCD22=(struct indnuis){0.255,0.018,1,"a0fT_BK_HPQCD22"};
	nuisparam->a1fT_BK_HPQCD22=(struct indnuis){-0.66,0.23,1,"a1fT_BK_HPQCD22"};
	nuisparam->a2fT_BK_HPQCD22=(struct indnuis){0.36,0.84,1,"a2fT_BK_HPQCD22"};
	nuisparam->L_BK_HPQCD22=(struct indnuis){1.304,0.010,1,"L_BK_HPQCD22"};

	/* SN: B->K FF Khodjamirian & Rusov 1703.04765 (BK_FormFactor_choice:8) - 02/08/2025*/
	nuisparam->f0fp_BK_KR=(struct indnuis){0.395,0.033,1,"f0fp_BK_KR"};
	nuisparam->b1fp_BK_KR=(struct indnuis){-1.42,1.52,1,"b1fp_BK_KR"};
	nuisparam->f0fT_BK_KR=(struct indnuis){0.381,0.027,1,"f0fT_BK_KR"};
	nuisparam->b1fT_BK_KR=(struct indnuis){-0.87,1.72,1,"b1fT_BK_KR"};


	
	/* low */
	nuisparam->BtoKlow_FV_err_noq2=(struct indnuis){0.,0.1,1,"BtoKlow_FV_err_noq2"};
	nuisparam->BtoKlow_FA_err_noq2=(struct indnuis){0.,0.1,1,"BtoKlow_FA_err_noq2"};
	nuisparam->BtoKlow_FS_err_noq2=(struct indnuis){0.,0.1,1,"BtoKlow_FS_err_noq2"};
	nuisparam->BtoKlow_FP_err_noq2=(struct indnuis){0.,0.1,1,"BtoKlow_FP_err_noq2"};

	nuisparam->BtoKlow_FV_err_q2=(struct indnuis){0.,0.25,1,"BtoKlow_FV_err_q2"};
	nuisparam->BtoKlow_FA_err_q2=(struct indnuis){0.,0.25,1,"BtoKlow_FA_err_q2"};
	nuisparam->BtoKlow_FS_err_q2=(struct indnuis){0.,0.25,1,"BtoKlow_FS_err_q2"};
	nuisparam->BtoKlow_FP_err_q2=(struct indnuis){0.,0.25,1,"BtoKlow_FP_err_q2"};

	/* high */
	nuisparam->BtoKhigh_FV_err=(struct indnuis){0.,0.1,1,"BtoKhigh_FV_err"};
	nuisparam->BtoKhigh_FA_err=(struct indnuis){0.,0.1,1,"BtoKhigh_FA_err"};
	nuisparam->BtoKhigh_FS_err=(struct indnuis){0.,0.1,1,"BtoKhigh_FS_err"};
	nuisparam->BtoKhigh_FP_err=(struct indnuis){0.,0.1,1,"BtoKhigh_FP_err"};

/* Bs */
	nuisparam->life_Bs=(struct indnuis){1.527e-12,0.011e-12,1,"life_Bs"};
	nuisparam->f_Bs=(struct indnuis){0.2303,0.0013,1,"f_Bs"};
	nuisparam->lambda_Bsp=(struct indnuis){0.46,0.11,1,"lambda_Bsp"};
	nuisparam->ys_Bs=(struct indnuis){0.0635,0.0035,1,"lambda_Bsp"};
	
/* Bs -> phi */
	nuisparam->f_phi_par=(struct indnuis){0.233,0.004,1,"f_phi_par"};
	nuisparam->f_phi_perp=(struct indnuis){0.191,0.004,1,"f_phi_perp"};
	nuisparam->a1phi_perp=(struct indnuis){0.,0.,1,"a1phi_perp"};
	nuisparam->a1phi_par=(struct indnuis){0.,0.,1,"a1phi_par"};
	nuisparam->a2phi_perp=(struct indnuis){0.14,0.07,1,"a2phi_perp"};
	nuisparam->a2phi_par=(struct indnuis){0.23,0.08,1,"a2phi_par"};

	/* low */
	nuisparam->Bstophilow_ALperp_err_noq2=(struct indnuis){0.,0.1,1,"Bstophilow_ALperp_err_noq2"};
	nuisparam->Bstophilow_ARperp_err_noq2=(struct indnuis){0.,0.1,1,"Bstophilow_ARperp_err_noq2"};
	nuisparam->Bstophilow_ALpar_err_noq2=(struct indnuis){0.,0.1,1,"Bstophilow_ALpar_err_noq2"};
	nuisparam->Bstophilow_ARpar_err_noq2=(struct indnuis){0.,0.1,1,"Bstophilow_ARpar_err_noq2"};
	nuisparam->Bstophilow_AL0_err_noq2=(struct indnuis){0.,0.1,1,"Bstophilow_AL0_err_noq2"};
	nuisparam->Bstophilow_AR0_err_noq2=(struct indnuis){0.,0.1,1,"Bstophilow_AR0_err_noq2"};
	nuisparam->Bstophilow_At_err_noq2=(struct indnuis){0.,0.1,1,"Bstophilow_At_err_noq2"};
	nuisparam->Bstophilow_AS_err_noq2=(struct indnuis){0.,0.1,1,"Bstophilow_AS_err_noq2"};

	nuisparam->Bstophilow_ALperp_err_q2=(struct indnuis){0.,0.25,1,"Bstophilow_ALperp_err_q2"};
	nuisparam->Bstophilow_ARperp_err_q2=(struct indnuis){0.,0.25,1,"Bstophilow_ARperp_err_q2"};
	nuisparam->Bstophilow_ALpar_err_q2=(struct indnuis){0.,0.25,1,"Bstophilow_ALpar_err_q2"};
	nuisparam->Bstophilow_ARpar_err_q2=(struct indnuis){0.,0.25,1,"Bstophilow_ARpar_err_q2"};
	nuisparam->Bstophilow_AL0_err_q2=(struct indnuis){0.,0.25,1,"Bstophilow_AL0_err_q2"};
	nuisparam->Bstophilow_AR0_err_q2=(struct indnuis){0.,0.25,1,"Bstophilow_AR0_err_q2"};
	nuisparam->Bstophilow_At_err_q2=(struct indnuis){0.,0.25,1,"Bstophilow_At_err_q2"};
	nuisparam->Bstophilow_AS_err_q2=(struct indnuis){0.,0.25,1,"Bstophilow_AS_err_q2"};
	
	/* high */
	nuisparam->Bstophihigh_ALperp_err=(struct indnuis){0.,0.1,1,"Bstophihigh_ALperp_err"};
	nuisparam->Bstophihigh_ARperp_err=(struct indnuis){0.,0.1,1,"Bstophihigh_ARperp_err"};
	nuisparam->Bstophihigh_ALpar_err=(struct indnuis){0.,0.1,1,"Bstophihigh_ALpar_err"};
	nuisparam->Bstophihigh_ARpar_err=(struct indnuis){0.,0.1,1,"Bstophihigh_ARpar_err"};
	nuisparam->Bstophihigh_AL0_err=(struct indnuis){0.,0.1,1,"Bstophihigh_AL0_err"};
	nuisparam->Bstophihigh_AR0_err=(struct indnuis){0.,0.1,1,"Bstophihigh_AR0_err"};
	nuisparam->Bstophihigh_At_err=(struct indnuis){0.,0.1,1,"Bstophihigh_At_err"};
	nuisparam->Bstophihigh_AS_err=(struct indnuis){0.,0.1,1,"Bstophihigh_AS_err"};
			
	/* Form factors Bs->phi ll */	
	nuisparam->a0A0_Bsphi=(struct indnuis){0.421328,0.0240514,1,"a0A0_Bsphi"};
	nuisparam->a1A0_Bsphi=(struct indnuis){-0.976454,0.23757,1,"a1A0_Bsphi"};
	nuisparam->a2A0_Bsphi=(struct indnuis){3.2714,1.35909,1,"a2A0_Bsphi"};
	nuisparam->a0A1_Bsphi=(struct indnuis){0.288007,0.0105759,1,"a0A1_Bsphi"};
	nuisparam->a1A1_Bsphi=(struct indnuis){0.350826,0.103813,1,"a1A1_Bsphi"};
	nuisparam->a2A1_Bsphi=(struct indnuis){1.69688,0.790297,1,"a2A1_Bsphi"};
	nuisparam->a0A12_Bsphi=(struct indnuis){0.267053,0.0152447,1,"a0A12_Bsphi"};
	nuisparam->a1A12_Bsphi=(struct indnuis){0.954402,0.125879,1,"a1A12_Bsphi"};
	nuisparam->a2A12_Bsphi=(struct indnuis){2.15263,0.47881,1,"a2A12_Bsphi"};
	nuisparam->a0V_Bsphi=(struct indnuis){0.364478,0.0141353,1,"a0V_Bsphi"};
	nuisparam->a1V_Bsphi=(struct indnuis){-1.22389,0.164176,1,"a1V_Bsphi"};
	nuisparam->a2V_Bsphi=(struct indnuis){3.74061,1.72682,1,"a2V_Bsphi"};
	nuisparam->a0T1_Bsphi=(struct indnuis){0.299475,0.0120688,1,"a0T1_Bsphi"};
	nuisparam->a1T1_Bsphi=(struct indnuis){-1.1013,0.0835241,1,"a1T1_Bsphi"};
	nuisparam->a2T1_Bsphi=(struct indnuis){0.58459,1.00324,1,"a2T1_Bsphi"};
	nuisparam->a0T2_Bsphi=(struct indnuis){0.299475,0.0120688,1,"a0T2_Bsphi"};
	nuisparam->a1T2_Bsphi=(struct indnuis){0.403564,0.0803563,1,"a1T2_Bsphi"};
	nuisparam->a2T2_Bsphi=(struct indnuis){1.03987,0.608872,1,"a2T2_Bsphi"};
	nuisparam->a0T23_Bsphi=(struct indnuis){0.65233,0.0357469,1,"a0T23_Bsphi"};
	nuisparam->a1T23_Bsphi=(struct indnuis){2.09622,0.330082,1,"a1T23_Bsphi"};
	nuisparam->a2T23_Bsphi=(struct indnuis){6.73572,1.79526,1,"a2T23_Bsphi"};

	/* SN: Bs->phi FF van Dyk et al. 2305.06301 (Bsphi_FormFactor_choice:2)*/
	nuisparam->a0A0_Bsphi_GRvDV_BSZ=(struct indnuis){0.38000311,0.041091335704744374,1, "a0A0_Bsphi_GRvDV_BSZ"};
	nuisparam->a1A0_Bsphi_GRvDV_BSZ=(struct indnuis){-1.27092227,0.39364220810273887,1, "a1A0_Bsphi_GRvDV_BSZ"};
	nuisparam->a2A0_Bsphi_GRvDV_BSZ=(struct indnuis){2.75074752,3.7986696355434755,1, "a2A0_Bsphi_GRvDV_BSZ"};
	nuisparam->a0A1_Bsphi_GRvDV_BSZ=(struct indnuis){0.30467219,0.0322816124132609,1, "a0A1_Bsphi_GRvDV_BSZ"};
	nuisparam->a1A1_Bsphi_GRvDV_BSZ=(struct indnuis){0.47300921,0.2949495439901544,1, "a1A1_Bsphi_GRvDV_BSZ"};
	nuisparam->a2A1_Bsphi_GRvDV_BSZ=(struct indnuis){1.96851627,1.8263304465512258,1, "a2A1_Bsphi_GRvDV_BSZ"};
	nuisparam->a1A12_Bsphi_GRvDV_BSZ=(struct indnuis){0.32277817,0.5477838487943945,1, "a1A12_Bsphi_GRvDV_BSZ"};
	nuisparam->a2A12_Bsphi_GRvDV_BSZ=(struct indnuis){-0.64655359,3.1437102363290417,1, "a2A12_Bsphi_GRvDV_BSZ"};
	nuisparam->a0V_Bsphi_GRvDV_BSZ=(struct indnuis){0.37782957,0.0506990521607653,1, "a0V_Bsphi_GRvDV_BSZ"};
	nuisparam->a1V_Bsphi_GRvDV_BSZ=(struct indnuis){-0.89784964,0.38314267185997436,1, "a1V_Bsphi_GRvDV_BSZ"};
	nuisparam->a2V_Bsphi_GRvDV_BSZ=(struct indnuis){3.89603141,4.018822066227864,1, "a2V_Bsphi_GRvDV_BSZ"};
	nuisparam->a0T1_Bsphi_GRvDV_BSZ=(struct indnuis){0.33630568,0.03560373547817701,1, "a0T1_Bsphi_GRvDV_BSZ"};
	nuisparam->a1T1_Bsphi_GRvDV_BSZ=(struct indnuis){-0.76226543,0.32113585754318996,1, "a1T1_Bsphi_GRvDV_BSZ"};
	nuisparam->a2T1_Bsphi_GRvDV_BSZ=(struct indnuis){0.99707755,2.8145841593386405,1, "a2T1_Bsphi_GRvDV_BSZ"};
	nuisparam->a1T2_Bsphi_GRvDV_BSZ=(struct indnuis){0.74223041,0.4237922911049704,1, "a1T2_Bsphi_GRvDV_BSZ"};
	nuisparam->a2T2_Bsphi_GRvDV_BSZ=(struct indnuis){2.00139229,2.576371580343177,1, "a2T2_Bsphi_GRvDV_BSZ"};
	nuisparam->a0T23_Bsphi_GRvDV_BSZ=(struct indnuis){0.63570005,0.05749057409697698,1, "a0T23_Bsphi_GRvDV_BSZ"};
	nuisparam->a1T23_Bsphi_GRvDV_BSZ=(struct indnuis){1.1784008,0.6973286083332592,1, "a1T23_Bsphi_GRvDV_BSZ"};
	nuisparam->a2T23_Bsphi_GRvDV_BSZ=(struct indnuis){1.79990344,5.038399378771,1, "a2T23_Bsphi_GRvDV_BSZ"};

	/* SN: Bs->phi FF Zwicky et al. 1503.05534 (Bsphi_FormFactor_choice:5)*/
	nuisparam->a0V_Bsphi_BSZ_LCSRonly=(struct indnuis){0.386556,0.033258,1, "a0V_Bsphi_BSZ_LCSRonly"};
	nuisparam->a1V_Bsphi_BSZ_LCSRonly=(struct indnuis){-1.02957,0.254779,1, "a1V_Bsphi_BSZ_LCSRonly"};
	nuisparam->a2V_Bsphi_BSZ_LCSRonly=(struct indnuis){3.50499,1.55247,1, "a2V_Bsphi_BSZ_LCSRonly"};
	nuisparam->a0T1_Bsphi_BSZ_LCSRonly=(struct indnuis){0.309106,0.0271593,1, "a0T1_Bsphi_BSZ_LCSRonly"};
	nuisparam->a1T1_Bsphi_BSZ_LCSRonly=(struct indnuis){-0.867808,0.190096,1, "a1T1_Bsphi_BSZ_LCSRonly"};
	nuisparam->a2T1_Bsphi_BSZ_LCSRonly=(struct indnuis){2.74535,1.19396,1, "a2T1_Bsphi_BSZ_LCSRonly"};
	nuisparam->a0T2_Bsphi_BSZ_LCSRonly=(struct indnuis){0.309106,0.0271593,1, "a0T2_Bsphi_BSZ_LCSRonly"};
	nuisparam->a1T2_Bsphi_BSZ_LCSRonly=(struct indnuis){0.58328,0.18578,1, "a1T2_Bsphi_BSZ_LCSRonly"};
	nuisparam->a2T2_Bsphi_BSZ_LCSRonly=(struct indnuis){0.893065,0.713992,1, "a2T2_Bsphi_BSZ_LCSRonly"};
	nuisparam->a0T23_Bsphi_BSZ_LCSRonly=(struct indnuis){0.675618,0.0707825,1, "a0T23_Bsphi_BSZ_LCSRonly"};
	nuisparam->a1T23_Bsphi_BSZ_LCSRonly=(struct indnuis){2.11435,0.46018,1, "a1T23_Bsphi_BSZ_LCSRonly"};
	nuisparam->a2T23_Bsphi_BSZ_LCSRonly=(struct indnuis){4.93525,2.24931,1, "a2T23_Bsphi_BSZ_LCSRonly"};
	nuisparam->a0A0_Bsphi_BSZ_LCSRonly=(struct indnuis){0.388675,0.0454248,1, "a0A0_Bsphi_BSZ_LCSRonly"};
	nuisparam->a1A0_Bsphi_BSZ_LCSRonly=(struct indnuis){-0.783874,0.262704,1, "a1A0_Bsphi_BSZ_LCSRonly"};
	nuisparam->a2A0_Bsphi_BSZ_LCSRonly=(struct indnuis){2.41473,1.4762,1, "a2A0_Bsphi_BSZ_LCSRonly"};
	nuisparam->a0A1_Bsphi_BSZ_LCSRonly=(struct indnuis){0.295607,0.026534,1, "a0A1_Bsphi_BSZ_LCSRonly"};
	nuisparam->a1A1_Bsphi_BSZ_LCSRonly=(struct indnuis){0.483914,0.194885,1, "a1A1_Bsphi_BSZ_LCSRonly"};
	nuisparam->a2A1_Bsphi_BSZ_LCSRonly=(struct indnuis){0.291817,0.646679,1, "a2A1_Bsphi_BSZ_LCSRonly"};
	nuisparam->a0A12_Bsphi_BSZ_LCSRonly=(struct indnuis){0.246357,0.0287919,1, "a0A12_Bsphi_BSZ_LCSRonly"};
	nuisparam->a1A12_Bsphi_BSZ_LCSRonly=(struct indnuis){0.764562,0.195545,1, "a1A12_Bsphi_BSZ_LCSRonly"};
	nuisparam->a2A12_Bsphi_BSZ_LCSRonly=(struct indnuis){0.711083,0.957192,1, "a2A12_Bsphi_BSZ_LCSRonly"};

	/* SN: Lattice-only Bs->phi FF Wingate et al. 1310.3722 + 1501.00367 (Bsphi_FormFactor_choice:6*/
	nuisparam->a0V_Bsphi_HLMW=(struct indnuis){0.452535336822,0.0303218620590161,1, "a0V_Bsphi_HLMW"};
	nuisparam->a1V_Bsphi_HLMW=(struct indnuis){-2.39909309211,0.496307597798986,1, "a1V_Bsphi_HLMW"};
	nuisparam->a0A0_Bsphi_HLMW=(struct indnuis){0.529352133811,0.0163321988745851,1, "a0A0_Bsphi_HLMW"};
	nuisparam->a1A0_Bsphi_HLMW=(struct indnuis){-1.59248596874,0.242539978251834,1, "a1A0_Bsphi_HLMW"};
	nuisparam->a0A1_Bsphi_HLMW=(struct indnuis){0.282831524425,0.011073567487445,1, "a0A1_Bsphi_HLMW"};
	nuisparam->a1A1_Bsphi_HLMW=(struct indnuis){0.123793996986,0.15010983866789,1, "a1A1_Bsphi_HLMW"};
	nuisparam->a0A12_Bsphi_HLMW=(struct indnuis){0.20660689485,0.00977633713173293,1, "a0A12_Bsphi_HLMW"};
	nuisparam->a1A12_Bsphi_HLMW=(struct indnuis){0.476257121291,0.147272721746765,1, "a1A12_Bsphi_HLMW"};
	nuisparam->a0T1_Bsphi_HLMW=(struct indnuis){0.401597756404,0.0105678841189237,1, "a0T1_Bsphi_HLMW"};
	nuisparam->a1T1_Bsphi_HLMW=(struct indnuis){-1.13396307551,0.132528712189472,1, "a1T1_Bsphi_HLMW"};
	nuisparam->a0T2_Bsphi_HLMW=(struct indnuis){0.282969214555,0.0075170904050703,1, "a0T2_Bsphi_HLMW"};
	nuisparam->a1T2_Bsphi_HLMW=(struct indnuis){0.248733391065,0.0974614020059223,1, "a1T2_Bsphi_HLMW"};
	nuisparam->a0T23_Bsphi_HLMW=(struct indnuis){0.521926043325,0.0156643526510992,1, "a0T23_Bsphi_HLMW"};
	nuisparam->a1T23_Bsphi_HLMW=(struct indnuis){0.384393588684,0.236134214881918,1, "a1T23_Bsphi_HLMW"};
	nuisparam->systErr_Bsphi_HLMW=(struct indnuis){0.,0.05,1, "systErr_Bsphi_HLMW"};


	/* Lambda_b -> Lambda l+l- */
	nuisparam->life_Lb=(struct indnuis){1.466e-12,0.010e-12,1,"life_Lb"};
    nuisparam->alphaL_LbLll=(struct indnuis){0.642,0.013,1,"alphaL_LbLll"};
    nuisparam->a0_HO_fplus_LbLll=(struct indnuis){0.422899,0.0273787,1,"a0_HO_fplus_LbLll"};
    nuisparam->a1_HO_fplus_LbLll=(struct indnuis){-1.37281,0.30678,1,"a1_HO_fplus_LbLll"};
    nuisparam->a2_HO_fplus_LbLll=(struct indnuis){1.79720,1.15063,1,"a2_HO_fplus_LbLll"};
    nuisparam->a0_HO_fperp_LbLll=(struct indnuis){0.514815,0.0352842,1,"a0_HO_fperp_LbLll"};
    nuisparam->a1_HO_fperp_LbLll=(struct indnuis){-1.47811,0.40296,1,"a1_HO_fperp_LbLll"};
    nuisparam->a2_HO_fperp_LbLll=(struct indnuis){1.24964,1.63958,1,"a2_HO_fperp_LbLll"};
    nuisparam->a0_HO_gpp_LbLll=(struct indnuis){0.352156,0.0205424,1,"a0_HO_gpp_LbLll"};
    nuisparam->a1_HO_gplus_LbLll=(struct indnuis){-1.29677,0.27319,1,"a1_HO_gplus_LbLll"};
    nuisparam->a2_HO_gplus_LbLll=(struct indnuis){2.71059,1.06655,1,"a2_HO_gplus_LbLll"};
    nuisparam->a1_HO_gperp_LbLll=(struct indnuis){-1.36073,0.29489,1,"a1_HO_gperp_LbLll"};
    nuisparam->a2_HO_gperp_LbLll=(struct indnuis){2.46213,1.37114,1,"a2_HO_gperp_LbLll"};
    nuisparam->a0_HO_hplus_LbLll=(struct indnuis){0.475314,0.0422776,1,"a0_HO_hplus_LbLll"};
    nuisparam->a1_HO_hplus_LbLll=(struct indnuis){-0.883985,0.3996916,1,"a1_HO_hplus_LbLll"};
    nuisparam->a2_HO_hplus_LbLll=(struct indnuis){-0.819016,1.6759566,1,"a2_HO_hplus_LbLll"};
    nuisparam->a0_HO_hperp_LbLll=(struct indnuis){0.374496,0.0312506,1,"a0_HO_hperp_LbLll"};
    nuisparam->a1_HO_hperp_LbLll=(struct indnuis){-0.943936,0.2766479,1,"a1_HO_hperp_LbLll"};
    nuisparam->a2_HO_hperp_LbLll=(struct indnuis){1.16058,1.07574,1,"a2_HO_hperp_LbLll"};
    nuisparam->a0_HO_htildepp_LbLll=(struct indnuis){0.325580,0.0247773,1,"a0_HO_htildepp_LbLll"};
    nuisparam->a1_HO_htildeplus_LbLll=(struct indnuis){-0.960313,0.2303462,1,"a1_HO_htildeplus_LbLll"};
    nuisparam->a2_HO_htildeplus_LbLll=(struct indnuis){2.97798,1.00407,1,"a2_HO_htildeplus_LbLll"};
    nuisparam->a1_HO_htildeperp_LbLll=(struct indnuis){-0.963386,0.2267987,1,"a1_HO_htildeperp_LbLll"};
    nuisparam->a2_HO_htildeperp_LbLll=(struct indnuis){2.47823,0.95487,1,"a2_HO_htildeperp_LbLll"};

	/* kaons */
	nuisparam->deltaPcu_Kppipnunu=(struct indnuis){0.04,0.02,1,"deltaPcu_Kppipnunu"};
	nuisparam->err_Pc_Xlambda_Kppipnunu=(struct indnuis){0.,1.,1,"err_Pc_Xlambda_Kppipnunu"};
	nuisparam->BR_KLgammagamma_exp=(struct indnuis){5.47e-4,0.04e-4,1,"BR_KLgammagamma_exp"};
	nuisparam->Aterm_mu_KLmumu=(struct indnuis){0.71,1.01,1,"Aterm_mu_KLmumu"};
	nuisparam->chi_gg_Mrho=(struct indnuis){5.83,1.01,1,"chi_gg_Mrho"};
	nuisparam->BR_KSgammagamma_exp=(struct indnuis){2.63e-6,0.81e-6,1,"BR_KSgammagamma_exp"};
	nuisparam->Iterm_mu_KSmumu=(struct indnuis){-2.821,0.001,1,"Iterm_mu_KSmumu"};
	nuisparam->log_muK_1GeV=(struct indnuis){0.,1.,2,"log_muK_1GeV"};
	nuisparam->KLpill_Ce_Cdir=(struct indnuis){4.62,0.24,1,"KLpill_Ce_Cdir"};
    nuisparam->KLpill_Cmu_Cdir=(struct indnuis){1.09,0.05,1,"KLpill_Cmu_Cdir"};
	nuisparam->KLpill_Ce_Cint=(struct indnuis){11.3,0.3,1,"KLpill_Ce_Cint"};
	nuisparam->KLpill_Cmu_Cint=(struct indnuis){2.63,0.06,1,"KLpill_Cmu_Cint"};
	nuisparam->KLpill_Ce_Cmix=(struct indnuis){14.5,0.5,1,"KLpill_Ce_Cmix"};
	nuisparam->KLpill_Cmu_Cmix=(struct indnuis){3.36,0.20,1,"KLpill_Cmu_Cmix"};
	nuisparam->KLpill_Cmu_CPC=(struct indnuis){5.2,1.6,1,"KLpill_Cmu_CPC"};
	nuisparam->KLpill_abs_aS=(struct indnuis){1.2,0.2,1,"KLpill_abs_aS"};
	/* Scalar and Pseudoscalar contributions included - 08/10/2025 */
	nuisparam->KLpill_Ce_Sgg=(struct indnuis){1.5e-4,0.3e-4,1,"KLpill_Ce_Sgg"};
	nuisparam->KLpill_Cmu_Sgg=(struct indnuis){4.0e-2,1.0e-2,1,"KLpill_Cmu_Sgg"};
	/* AFBmu_KLpi0mumu  - 09/10/2025 */
    nuisparam->KLpill_AFB_a1pi_gg_7V=(struct indnuis){1.3,0.1,1,"KLpill_AFB_a1pi_gg_7V"};
	nuisparam->KLpill_AFB_a1pi_gg_aS=(struct indnuis){1.7,0.2,1,"KLpill_AFB_a1pi_gg_aS"};

	return;
}

/*--------------------------------------------------------------------*/

void set_nuisance_deviation_to_zero(struct nuisance* nuisparam)
/* Sets the standard deviations of the parameters of the structure nuisparam to zero */
{
	nuisparam->alphas_MZ.dev=0.;
	nuisparam->mass_b.dev=0.;
	nuisparam->mass_c.dev=0.;
	nuisparam->mass_s.dev=0.;
	nuisparam->mass_top_pole.dev=0.;

	nuisparam->mass_h0.dev=0.;

	nuisparam->CKM_lambda.dev=0.;
	nuisparam->CKM_A.dev=0.;
	nuisparam->CKM_rhobar.dev=0.;
	nuisparam->CKM_etabar.dev=0.;

	nuisparam->log_mu_W_mass_W.dev=0.;
	nuisparam->log_mu_b_mass_b.dev=0.;

/* inclusive b -> s */
	nuisparam->BR_BXclnu_exp.dev=0.;
/* b -> s gamma */ 
	nuisparam->mu_G2_bsg.dev=0.;
	nuisparam->rho_D3_bsg.dev=0.;
	nuisparam->rho_LS3_bsg.dev=0.;
	nuisparam->bsgamma_rand.dev=0.;
	nuisparam->mu_c_bsg.dev=0.;
/* b -> s mu mu */
	nuisparam->BRBXsmumu_lowq2_rand.dev=0.;
	nuisparam->BRBXsmumu_highq2_rand.dev=0.;
	nuisparam->BRBXsmumu_full_rand.dev=0.;
/* b -> s e e */
	nuisparam->BRBXsee_lowq2_rand.dev=0.;
	nuisparam->BRBXsee_highq2_rand.dev=0.;
	nuisparam->BRBXsee_full_rand.dev=0.;
/* b -> s tau tau */
	nuisparam->BRBXstautau_lowq2_rand.dev=0.;
	nuisparam->BRBXstautau_highq2_rand.dev=0.;
	nuisparam->BRBXstautau_full_rand.dev=0.;

/* B */	
	nuisparam->f_B.dev=0.;
	nuisparam->lambda_Bp.dev=0.;
/* B -> K* */
	nuisparam->f_Kstar_par.dev=0.;
	nuisparam->f_Kstar_perp.dev=0.;
	nuisparam->a1perp.dev=0.;
	nuisparam->a2perp.dev=0.;
	nuisparam->a1par.dev=0.;
	nuisparam->a2par.dev=0.;
/* B -> K* gamma */
	nuisparam->T1_BKstar.dev=0.;
	nuisparam->log_mu_spec_lambda_h_mass_b.dev=0.;
	/* low */
	nuisparam->BtoKstarlow_ALperp_err_noq2.dev=0.;
	nuisparam->BtoKstarlow_ARperp_err_noq2.dev=0.;
	nuisparam->BtoKstarlow_ALpar_err_noq2.dev=0.;
	nuisparam->BtoKstarlow_ARpar_err_noq2.dev=0.;
	nuisparam->BtoKstarlow_AL0_err_noq2.dev=0.;
	nuisparam->BtoKstarlow_AR0_err_noq2.dev=0.;
	nuisparam->BtoKstarlow_At_err_noq2.dev=0.;
	nuisparam->BtoKstarlow_AS_err_noq2.dev=0.;

	nuisparam->BtoKstarlow_ALperp_err_q2.dev=0.;
	nuisparam->BtoKstarlow_ARperp_err_q2.dev=0.;
	nuisparam->BtoKstarlow_ALpar_err_q2.dev=0.;
	nuisparam->BtoKstarlow_ARpar_err_q2.dev=0.;
	nuisparam->BtoKstarlow_AL0_err_q2.dev=0.;
	nuisparam->BtoKstarlow_AR0_err_q2.dev=0.;
	nuisparam->BtoKstarlow_At_err_q2.dev=0.;
	nuisparam->BtoKstarlow_AS_err_q2.dev=0.;

	/* SN: Result from B->K* Hadronic Fit (to be used to check how it impacts predictions) - 13/04/2025 */
	nuisparam->hplus0r.dev=0.;
	nuisparam->hminus0r.dev=0.;
	nuisparam->hplus1r.dev=0.;
	nuisparam->hminus1r.dev=0.;
	nuisparam->hplus2r.dev=0.;
	nuisparam->hminus2r.dev=0.;
	nuisparam->hzero0r.dev=0.;
	nuisparam->hzero1r.dev=0.;
	nuisparam->hzero2r.dev=0.;
	nuisparam->hplus0i.dev=0.;
	nuisparam->hminus0i.dev=0.;
	nuisparam->hplus1i.dev=0.;
	nuisparam->hminus1i.dev=0.;
	nuisparam->hplus2i.dev=0.;
	nuisparam->hminus2i.dev=0.;
	nuisparam->hzero0i.dev=0.;
	nuisparam->hzero1i.dev=0.;
	nuisparam->hzero2i.dev=0.;


	nuisparam->real_alpha_perp0.dev=0.;
	nuisparam->real_alpha_perp1.dev=0.;
	nuisparam->real_alpha_perp2.dev=0.;
	nuisparam->real_alpha_par0.dev=0.;
	nuisparam->real_alpha_par1.dev=0.;
	nuisparam->real_alpha_par2.dev=0.;
	nuisparam->real_alpha_zero0.dev=0.;
	nuisparam->real_alpha_zero1.dev=0.;
	nuisparam->imag_alpha_perp0.dev=0.;
	nuisparam->imag_alpha_perp1.dev=0.;
	nuisparam->imag_alpha_perp2.dev=0.;
	nuisparam->imag_alpha_par0.dev=0.;
	nuisparam->imag_alpha_par1.dev=0.;
	nuisparam->imag_alpha_par2.dev=0.;
	nuisparam->imag_alpha_zero0.dev=0.;
	nuisparam->imag_alpha_zero1.dev=0.;
	
	nuisparam->DeltaC9_M1_q2bar.dev=0.;
	nuisparam->r1_M1.dev=0.;
	nuisparam->r2_M1.dev=0.;
	nuisparam->DeltaC9_M2_q2bar.dev=0.;
	nuisparam->r1_M2.dev=0.;
	nuisparam->r2_M2.dev=0.;
	nuisparam->DeltaC9_M3_q2bar.dev=0.;
	nuisparam->r1_M3.dev=0.;
	nuisparam->r2_M3.dev=0.;
	
	/* high */
	nuisparam->BtoKstarhigh_ALperp_err.dev=0.;
	nuisparam->BtoKstarhigh_ARperp_err.dev=0.;
	nuisparam->BtoKstarhigh_ALpar_err.dev=0.;
	nuisparam->BtoKstarhigh_ARpar_err.dev=0.;
	nuisparam->BtoKstarhigh_AL0_err.dev=0.;
	nuisparam->BtoKstarhigh_AR0_err.dev=0.;
	nuisparam->BtoKstarhigh_At_err.dev=0.;
	nuisparam->BtoKstarhigh_AS_err.dev=0.;

	/* Form factors B->K* ll */	
	nuisparam->a0A0_BKstar.dev=0.;
	nuisparam->a1A0_BKstar.dev=0.;
	nuisparam->a2A0_BKstar.dev=0.;
	nuisparam->a0A1_BKstar.dev=0.;
	nuisparam->a1A1_BKstar.dev=0.;
	nuisparam->a2A1_BKstar.dev=0.;
	nuisparam->a0A12_BKstar.dev=0.;
	nuisparam->a1A12_BKstar.dev=0.;
	nuisparam->a2A12_BKstar.dev=0.;
	nuisparam->a0V_BKstar.dev=0.;
	nuisparam->a1V_BKstar.dev=0.;
	nuisparam->a2V_BKstar.dev=0.;
	nuisparam->a0T1_BKstar.dev=0.;
	nuisparam->a1T1_BKstar.dev=0.;
	nuisparam->a2T1_BKstar.dev=0.;
	nuisparam->a0T2_BKstar.dev=0.;
	nuisparam->a1T2_BKstar.dev=0.;
	nuisparam->a2T2_BKstar.dev=0.;
	nuisparam->a0T23_BKstar.dev=0.;
	nuisparam->a1T23_BKstar.dev=0.;
	nuisparam->a2T23_BKstar.dev=0.;

	/* SN: B -> Kstar FF van Dyk et al. 2305.06301 (BKstar_FormFactor_choice:2)*/
	nuisparam->a0A0_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a1A0_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a2A0_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a0A1_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a1A1_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a2A1_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a1A12_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a2A12_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a0V_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a1V_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a2V_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a0T1_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a1T1_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a2T1_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a1T2_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a2T2_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a0T23_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a1T23_BKstar_GRvDV_BSZ.dev=0.;
	nuisparam->a2T23_BKstar_GRvDV_BSZ.dev=0.;

	/* SN: LCSR+Lattice B -> Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:3)*/
	nuisparam->a0A0_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a1A0_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a2A0_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a0A1_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a1A1_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a2A1_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a1A12_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a2A12_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a0V_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a1V_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a2V_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a0T1_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a1T1_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a2T1_BKstar_GKvD_LCSR_Lattice.dev=0.;
	// nuisparam->a0T2_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a1T2_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a2T2_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a0T23_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a1T23_BKstar_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a2T23_BKstar_GKvD_LCSR_Lattice.dev=0.;

	/* SN: LCSR-only B -> Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:4)*/
	nuisparam->a0A0_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a1A0_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a2A0_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a0A1_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a1A1_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a2A1_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a1A12_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a2A12_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a0V_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a1V_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a2V_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a0T1_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a1T1_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a2T1_BKstar_GKvD_LCSRonly.dev=0.;
	// nuisparam->a0T2_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a1T2_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a2T2_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a0T23_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a1T23_BKstar_GKvD_LCSRonly.dev=0.;
	nuisparam->a2T23_BKstar_GKvD_LCSRonly.dev=0.;

	/* SN: LCSR-only B->Kstar FF Zwicky et al. 1503.05534 (BKstar_FormFactor_choice:5)*/
	nuisparam->a0A0_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a1A0_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a2A0_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a0A1_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a1A1_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a2A1_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a0A12_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a1A12_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a2A12_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a0V_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a1V_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a2V_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a0T1_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a1T1_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a2T1_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a0T2_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a1T2_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a2T2_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a0T23_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a1T23_BKstar_BSZ_LCSRonly.dev=0.;
	nuisparam->a2T23_BKstar_BSZ_LCSRonly.dev=0.;

	/* SN: Lattice-only B->Kstar FF Wingate et al. 1310.3722 + 1501.00367 (BKstar_FormFactor_choice:6*/
	nuisparam->a0A0_BKstar_HLMW.dev=0.;
	nuisparam->a1A0_BKstar_HLMW.dev=0.;
	nuisparam->a0A1_BKstar_HLMW.dev=0.;
	nuisparam->a1A1_BKstar_HLMW.dev=0.;
	nuisparam->a0A12_BKstar_HLMW.dev=0.;
	nuisparam->a1A12_BKstar_HLMW.dev=0.;
	nuisparam->a0V_BKstar_HLMW.dev=0.;
	nuisparam->a1V_BKstar_HLMW.dev=0.;
	nuisparam->a0T1_BKstar_HLMW.dev=0.;
	nuisparam->a1T1_BKstar_HLMW.dev=0.;
	nuisparam->a0T2_BKstar_HLMW.dev=0.;
	nuisparam->a1T2_BKstar_HLMW.dev=0.;
	nuisparam->a0T23_BKstar_HLMW.dev=0.;
	nuisparam->a1T23_BKstar_HLMW.dev=0.;
	nuisparam->systErr_BKstar_HLMW.dev=0.;


/* B -> K */
	nuisparam->f_K.dev=0.;
	nuisparam->a1K.dev=0.;
	nuisparam->a2K.dev=0.;
	/* Form factors B->K ll */
	/* SN: Renaming form factor parameters to be more explicit - 28/03/2025 */
	// nuisparam->a00_BK.dev=0.;
	// nuisparam->a10_BK.dev=0.;
	// nuisparam->a20_BK.dev=0.;
	// nuisparam->a0p_BK.dev=0.;
	// nuisparam->a1p_BK.dev=0.;
	// nuisparam->a2p_BK.dev=0.;
	// nuisparam->a0T_BK.dev=0.;
	// nuisparam->a1T_BK.dev=0.;
	// nuisparam->a2T_BK.dev=0.;
	/* SN: B->K FF Altmannshofer, Straub 1411.3161 (BK_FormFactor_choice:1)*/
	nuisparam->a0f0_BK_AS_LCSR_Lattice.dev=0.;
	nuisparam->a1f0_BK_AS_LCSR_Lattice.dev=0.;
	nuisparam->a2f0_BK_AS_LCSR_Lattice.dev=0.;
	nuisparam->a3f0_BK_AS_LCSR_Lattice.dev=0.;
	nuisparam->a0fp_BK_AS_LCSR_Lattice.dev=0.;
	nuisparam->a1fp_BK_AS_LCSR_Lattice.dev=0.;
	nuisparam->a2fp_BK_AS_LCSR_Lattice.dev=0.;
	nuisparam->a0fT_BK_AS_LCSR_Lattice.dev=0.;
	nuisparam->a1fT_BK_AS_LCSR_Lattice.dev=0.;
	nuisparam->a2fT_BK_AS_LCSR_Lattice.dev=0.;
	/* SN: B->K FF van Dyk et al. 2305.06301 (BK_FormFactor_choice:2)*/
	nuisparam->a0fp_BK_GRvDV_BSZ.dev=0.;
	nuisparam->a1fp_BK_GRvDV_BSZ.dev=0.;
	nuisparam->a2fp_BK_GRvDV_BSZ.dev=0.;
	nuisparam->a1f0_BK_GRvDV_BSZ.dev=0.;
	nuisparam->a2f0_BK_GRvDV_BSZ.dev=0.;
	nuisparam->a0fT_BK_GRvDV_BSZ.dev=0.;
	nuisparam->a1fT_BK_GRvDV_BSZ.dev=0.;
	nuisparam->a2fT_BK_GRvDV_BSZ.dev=0.;
	/* SN: LCSR+Lattice B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:3)*/
	nuisparam->a0fp_BK_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a1fp_BK_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a2fp_BK_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a1f0_BK_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a2f0_BK_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a0fT_BK_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a1fT_BK_GKvD_LCSR_Lattice.dev=0.;
	nuisparam->a2fT_BK_GKvD_LCSR_Lattice.dev=0.;
	/* SN: LCSR-only B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:4)*/
	nuisparam->a0fp_BK_GKvD_LCSRonly.dev=0.;
	nuisparam->a1fp_BK_GKvD_LCSRonly.dev=0.;
	nuisparam->a2fp_BK_GKvD_LCSRonly.dev=0.;
	nuisparam->a1f0_BK_GKvD_LCSRonly.dev=0.;
	nuisparam->a2f0_BK_GKvD_LCSRonly.dev=0.;
	nuisparam->a0fT_BK_GKvD_LCSRonly.dev=0.;
	nuisparam->a1fT_BK_GKvD_LCSRonly.dev=0.;
	nuisparam->a2fT_BK_GKvD_LCSRonly.dev=0.;
	/* SN: Lattice FLAG24 average 2411.04268 of HPQCD'13 and FNAL/MILC'15 (BK_FormFactor_choice:6) - 26/04/2025*/
	nuisparam->a0fp_BK_FLAG24.dev=0.;
	nuisparam->a1fp_BK_FLAG24.dev=0.;
	nuisparam->a2fp_BK_FLAG24.dev=0.;
	nuisparam->a0f0_BK_FLAG24.dev=0.;
	nuisparam->a1f0_BK_FLAG24.dev=0.;
	nuisparam->a0fT_BK_FLAG24.dev=0.;
	nuisparam->a1fT_BK_FLAG24.dev=0.;
	nuisparam->a2fT_BK_FLAG24.dev=0.;
	/* SN: Lattice full q2-region HPQCD22 2207.12468 (BK_FormFactor_choice:7) - 26/04/2025*/
	nuisparam->a0f0_BK_HPQCD22.dev=0.;
	nuisparam->a1f0_BK_HPQCD22.dev=0.;
	nuisparam->a2f0_BK_HPQCD22.dev=0.;
	nuisparam->a1fp_BK_HPQCD22.dev=0.;
	nuisparam->a2fp_BK_HPQCD22.dev=0.;
	nuisparam->a0fT_BK_HPQCD22.dev=0.;
	nuisparam->a1fT_BK_HPQCD22.dev=0.;
	nuisparam->a2fT_BK_HPQCD22.dev=0.;
	nuisparam->L_BK_HPQCD22.dev=0.;
	/* SN: LCSR-only Khodjamirian & Rusov 1703.04765 (BK_FormFactor_choice:8) - 02/08/2025*/
	nuisparam->f0fp_BK_KR.dev=0.;
	nuisparam->b1fp_BK_KR.dev=0.;
	nuisparam->f0fT_BK_KR.dev=0.;
	nuisparam->b1fT_BK_KR.dev=0.;


	/* low */
	nuisparam->BtoKlow_FV_err_noq2.dev=0.;
	nuisparam->BtoKlow_FA_err_noq2.dev=0.;
	nuisparam->BtoKlow_FS_err_noq2.dev=0.;
	nuisparam->BtoKlow_FP_err_noq2.dev=0.;

	nuisparam->BtoKlow_FV_err_q2.dev=0.;
	nuisparam->BtoKlow_FA_err_q2.dev=0.;
	nuisparam->BtoKlow_FS_err_q2.dev=0.;
	nuisparam->BtoKlow_FP_err_q2.dev=0.;

	/* high */
	nuisparam->BtoKhigh_FV_err.dev=0.;
	nuisparam->BtoKhigh_FA_err.dev=0.;
	nuisparam->BtoKhigh_FS_err.dev=0.;
	nuisparam->BtoKhigh_FP_err.dev=0.;

/* Bs */
	nuisparam->life_Bs.dev=0.;
	nuisparam->f_Bs.dev=0.;
	nuisparam->lambda_Bsp.dev=0.;
	nuisparam->ys_Bs.dev=0.;
	
/* Bs -> phi */
	nuisparam->f_phi_par.dev=0.;
	nuisparam->f_phi_perp.dev=0.;
	nuisparam->a1phi_perp.dev=0.;
	nuisparam->a1phi_par.dev=0.;
	nuisparam->a2phi_perp.dev=0.;
	nuisparam->a2phi_par.dev=0.;

	/* low */
	nuisparam->Bstophilow_ALperp_err_noq2.dev=0.;
	nuisparam->Bstophilow_ARperp_err_noq2.dev=0.;
	nuisparam->Bstophilow_ALpar_err_noq2.dev=0.;
	nuisparam->Bstophilow_ARpar_err_noq2.dev=0.;
	nuisparam->Bstophilow_AL0_err_noq2.dev=0.;
	nuisparam->Bstophilow_AR0_err_noq2.dev=0.;
	nuisparam->Bstophilow_At_err_noq2.dev=0.;
	nuisparam->Bstophilow_AS_err_noq2.dev=0.;

	nuisparam->Bstophilow_ALperp_err_q2.dev=0.;
	nuisparam->Bstophilow_ARperp_err_q2.dev=0.;
	nuisparam->Bstophilow_ALpar_err_q2.dev=0.;
	nuisparam->Bstophilow_ARpar_err_q2.dev=0.;
	nuisparam->Bstophilow_AL0_err_q2.dev=0.;
	nuisparam->Bstophilow_AR0_err_q2.dev=0.;
	nuisparam->Bstophilow_At_err_q2.dev=0.;
	nuisparam->Bstophilow_AS_err_q2.dev=0.;
	
	/* high */
	nuisparam->Bstophihigh_ALperp_err.dev=0.;
	nuisparam->Bstophihigh_ARperp_err.dev=0.;
	nuisparam->Bstophihigh_ALpar_err.dev=0.;
	nuisparam->Bstophihigh_ARpar_err.dev=0.;
	nuisparam->Bstophihigh_AL0_err.dev=0.;
	nuisparam->Bstophihigh_AR0_err.dev=0.;
	nuisparam->Bstophihigh_At_err.dev=0.;
	nuisparam->Bstophihigh_AS_err.dev=0.;
			
	/* Form factors Bs->phi ll */	
	nuisparam->a0A0_Bsphi.dev=0.;
	nuisparam->a1A0_Bsphi.dev=0.;
	nuisparam->a2A0_Bsphi.dev=0.;
	nuisparam->a0A1_Bsphi.dev=0.;
	nuisparam->a1A1_Bsphi.dev=0.;
	nuisparam->a2A1_Bsphi.dev=0.;
	nuisparam->a0A12_Bsphi.dev=0.;
	nuisparam->a1A12_Bsphi.dev=0.;
	nuisparam->a2A12_Bsphi.dev=0.;
	nuisparam->a0V_Bsphi.dev=0.;
	nuisparam->a1V_Bsphi.dev=0.;
	nuisparam->a2V_Bsphi.dev=0.;
	nuisparam->a0T1_Bsphi.dev=0.;
	nuisparam->a1T1_Bsphi.dev=0.;
	nuisparam->a2T1_Bsphi.dev=0.;
	nuisparam->a0T2_Bsphi.dev=0.;
	nuisparam->a1T2_Bsphi.dev=0.;
	nuisparam->a2T2_Bsphi.dev=0.;
	nuisparam->a0T23_Bsphi.dev=0.;
	nuisparam->a1T23_Bsphi.dev=0.;
	nuisparam->a2T23_Bsphi.dev=0.;

	/* SN: Bs->phi FF van Dyk et al. 2305.06301 (Bsphi_FormFactor_choice:2)*/
	nuisparam->a0A0_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a1A0_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a2A0_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a0A1_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a1A1_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a2A1_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a1A12_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a2A12_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a0V_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a1V_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a2V_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a0T1_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a1T1_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a2T1_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a1T2_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a2T2_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a0T23_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a1T23_Bsphi_GRvDV_BSZ.dev=0.;
	nuisparam->a2T23_Bsphi_GRvDV_BSZ.dev=0.;

	/* SN: Bs->phi FF Zwicky et al. 1503.05534 (Bsphi_FormFactor_choice:5)*/
	nuisparam->a0A0_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a1A0_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a2A0_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a0A1_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a1A1_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a2A1_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a0A12_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a1A12_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a2A12_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a0V_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a1V_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a2V_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a0T1_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a1T1_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a2T1_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a0T2_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a1T2_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a2T2_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a0T23_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a1T23_Bsphi_BSZ_LCSRonly.dev=0.;
	nuisparam->a2T23_Bsphi_BSZ_LCSRonly.dev=0.;

	/* SN: Lattice-only Bs->phi FF Wingate et al. 1310.3722 + 1501.00367 (Bsphi_FormFactor_choice:6*/
	nuisparam->a0A0_Bsphi_HLMW.dev=0.;
	nuisparam->a1A0_Bsphi_HLMW.dev=0.;
	nuisparam->a0A1_Bsphi_HLMW.dev=0.;
	nuisparam->a1A1_Bsphi_HLMW.dev=0.;
	nuisparam->a0A12_Bsphi_HLMW.dev=0.;
	nuisparam->a1A12_Bsphi_HLMW.dev=0.;
	nuisparam->a0V_Bsphi_HLMW.dev=0.;
	nuisparam->a1V_Bsphi_HLMW.dev=0.;
	nuisparam->a0T1_Bsphi_HLMW.dev=0.;
	nuisparam->a1T1_Bsphi_HLMW.dev=0.;
	nuisparam->a0T2_Bsphi_HLMW.dev=0.;
	nuisparam->a1T2_Bsphi_HLMW.dev=0.;
	nuisparam->a0T23_Bsphi_HLMW.dev=0.;
	nuisparam->a1T23_Bsphi_HLMW.dev=0.;
	nuisparam->systErr_Bsphi_HLMW.dev=0.;

	/* Lambda_b -> Lambda l+l- */
	nuisparam->life_Lb.dev=0.;
    nuisparam->alphaL_LbLll.dev=0.;
    nuisparam->a0_HO_fplus_LbLll.dev=0.;
    nuisparam->a1_HO_fplus_LbLll.dev=0.;
    nuisparam->a2_HO_fplus_LbLll.dev=0.;
    nuisparam->a0_HO_fperp_LbLll.dev=0.;
    nuisparam->a1_HO_fperp_LbLll.dev=0.;
    nuisparam->a2_HO_fperp_LbLll.dev=0.;
    nuisparam->a0_HO_gpp_LbLll.dev=0.;
    nuisparam->a1_HO_gplus_LbLll.dev=0.;
    nuisparam->a2_HO_gplus_LbLll.dev=0.;
    nuisparam->a1_HO_gperp_LbLll.dev=0.;
    nuisparam->a2_HO_gperp_LbLll.dev=0.;
    nuisparam->a0_HO_hplus_LbLll.dev=0.;
    nuisparam->a1_HO_hplus_LbLll.dev=0.;
    nuisparam->a2_HO_hplus_LbLll.dev=0.;
    nuisparam->a0_HO_hperp_LbLll.dev=0.;
    nuisparam->a1_HO_hperp_LbLll.dev=0.;
    nuisparam->a2_HO_hperp_LbLll.dev=0.;
    nuisparam->a0_HO_htildepp_LbLll.dev=0.;
    nuisparam->a1_HO_htildeplus_LbLll.dev=0.;
    nuisparam->a2_HO_htildeplus_LbLll.dev=0.;
    nuisparam->a1_HO_htildeperp_LbLll.dev=0.;
    nuisparam->a2_HO_htildeperp_LbLll.dev=0.;

	/* Kaons */
	nuisparam->deltaPcu_Kppipnunu.dev=0.;
	nuisparam->err_Pc_Xlambda_Kppipnunu.dev=0.;
	nuisparam->BR_KLgammagamma_exp.dev=0.;
	nuisparam->Aterm_mu_KLmumu.dev=0.;
	nuisparam->chi_gg_Mrho.dev=0.;
	nuisparam->BR_KSgammagamma_exp.dev=0.;
	nuisparam->Iterm_mu_KSmumu.dev=0.;
	nuisparam->log_muK_1GeV.dev=0.;
	nuisparam->KLpill_Ce_Cdir.dev=0.;
    nuisparam->KLpill_Cmu_Cdir.dev=0.;
	nuisparam->KLpill_Ce_Cint.dev=0.;
	nuisparam->KLpill_Cmu_Cint.dev=0.;
	nuisparam->KLpill_Ce_Cmix.dev=0.;
	nuisparam->KLpill_Cmu_Cmix.dev=0.;
	nuisparam->KLpill_Cmu_CPC.dev=0.;
	nuisparam->KLpill_abs_aS.dev=0.;
	/* Scalar and Pseudoscalar contributions included - 08/10/2025 */
	nuisparam->KLpill_Ce_Sgg.dev=0.;
	nuisparam->KLpill_Cmu_Sgg.dev=0.;
	/* AFBmu_KLpi0mumu  - 09/10/2025 */
    nuisparam->KLpill_AFB_a1pi_gg_7V.dev=0.;
	nuisparam->KLpill_AFB_a1pi_gg_aS.dev=0.;

	return;
}

/*--------------------------------------------------------------------*/

void set_nuisance_value_from_param(struct nuisance* nuisparam, struct parameters* param)
/* Sets the values of the nuisance parameters from the param structure */
{
	#define setparam( myparam ) nuisparam->myparam.cent=param->myparam;
	
	setparam(alphas_MZ);

	setparam(mass_b);
	setparam(mass_c);
	setparam(mass_s);
	setparam(mass_top_pole);

	setparam(mass_h0);

	setparam(CKM_lambda);
	setparam(CKM_A);
	setparam(CKM_rhobar);
	setparam(CKM_etabar);

	setparam(log_mu_W_mass_W);
	setparam(log_mu_b_mass_b);

/* inclusive b -> s */
	setparam(BR_BXclnu_exp);
/* b -> s gamma */ 
	setparam(mu_G2_bsg);
	setparam(rho_D3_bsg);
	setparam(rho_LS3_bsg);
	setparam(bsgamma_rand);
	setparam(mu_c_bsg);
/* b -> s mu mu */
	setparam(BRBXsmumu_lowq2_rand);
	setparam(BRBXsmumu_highq2_rand);
	setparam(BRBXsmumu_full_rand);
/* b -> s e e */
	setparam(BRBXsee_lowq2_rand);
	setparam(BRBXsee_highq2_rand);
	setparam(BRBXsee_full_rand);
/* b -> s tau tau */
	setparam(BRBXstautau_lowq2_rand);
	setparam(BRBXstautau_highq2_rand);
	setparam(BRBXstautau_full_rand);

/* B */	
	setparam(f_B);
	setparam(lambda_Bp);
/* B -> K* */
	setparam(f_Kstar_par);
	setparam(f_Kstar_perp);
	setparam(a1perp);
	setparam(a2perp);
	setparam(a1par);
	setparam(a2par);
/* B -> K* gamma */
	setparam(T1_BKstar);
	setparam(log_mu_spec_lambda_h_mass_b);
	/* low */
	setparam(BtoKstarlow_ALperp_err_noq2);
	setparam(BtoKstarlow_ARperp_err_noq2);
	setparam(BtoKstarlow_ALpar_err_noq2);
	setparam(BtoKstarlow_ARpar_err_noq2);
	setparam(BtoKstarlow_AL0_err_noq2);
	setparam(BtoKstarlow_AR0_err_noq2);
	setparam(BtoKstarlow_At_err_noq2);
	setparam(BtoKstarlow_AS_err_noq2);

	setparam(BtoKstarlow_ALperp_err_q2);
	setparam(BtoKstarlow_ARperp_err_q2);
	setparam(BtoKstarlow_ALpar_err_q2);
	setparam(BtoKstarlow_ARpar_err_q2);
	setparam(BtoKstarlow_AL0_err_q2);
	setparam(BtoKstarlow_AR0_err_q2);
	setparam(BtoKstarlow_At_err_q2);
	setparam(BtoKstarlow_AS_err_q2);

	/* SN: Result from B->K* Hadronic Fit (to be used to check how it impacts predictions) - 13/04/2025 */
	setparam(hplus0r);
	setparam(hminus0r);
	setparam(hplus1r);
	setparam(hminus1r);
	setparam(hplus2r);
	setparam(hminus2r);
	setparam(hzero0r);
	setparam(hzero1r);
	setparam(hzero2r);
	setparam(hplus0i);
	setparam(hminus0i);
	setparam(hplus1i);
	setparam(hminus1i);
	setparam(hplus2i);
	setparam(hminus2i);
	setparam(hzero0i);
	setparam(hzero1i);
	setparam(hzero2i);

	setparam(real_alpha_perp0);
	setparam(real_alpha_perp1);
	setparam(real_alpha_perp2);
	setparam(real_alpha_par0);
	setparam(real_alpha_par1);
	setparam(real_alpha_par2);
	setparam(real_alpha_zero0);
	setparam(real_alpha_zero1);
	setparam(imag_alpha_perp0);
	setparam(imag_alpha_perp1);
	setparam(imag_alpha_perp2);
	setparam(imag_alpha_par0);
	setparam(imag_alpha_par1);
	setparam(imag_alpha_par2);
	setparam(imag_alpha_zero0);
	setparam(imag_alpha_zero1);
	
	setparam(DeltaC9_M1_q2bar);
	setparam(r1_M1);
	setparam(r2_M1);
	setparam(DeltaC9_M2_q2bar);
	setparam(r1_M2);
	setparam(r2_M2);
	setparam(DeltaC9_M3_q2bar);
	setparam(r1_M3);
	setparam(r2_M3);
	
	/* high */
	setparam(BtoKstarhigh_ALperp_err);
	setparam(BtoKstarhigh_ARperp_err);
	setparam(BtoKstarhigh_ALpar_err);
	setparam(BtoKstarhigh_ARpar_err);
	setparam(BtoKstarhigh_AL0_err);
	setparam(BtoKstarhigh_AR0_err);
	setparam(BtoKstarhigh_At_err);
	setparam(BtoKstarhigh_AS_err);

	/* Form factors B->K* ll */	
	setparam(a0A0_BKstar);
	setparam(a1A0_BKstar);
	setparam(a2A0_BKstar);
	setparam(a0A1_BKstar);
	setparam(a1A1_BKstar);
	setparam(a2A1_BKstar);
	setparam(a0A12_BKstar);
	setparam(a1A12_BKstar);
	setparam(a2A12_BKstar);
	setparam(a0V_BKstar);
	setparam(a1V_BKstar);
	setparam(a2V_BKstar);
	setparam(a0T1_BKstar);
	setparam(a1T1_BKstar);
	setparam(a2T1_BKstar);
	setparam(a0T2_BKstar);
	setparam(a1T2_BKstar);
	setparam(a2T2_BKstar);
	setparam(a0T23_BKstar);
	setparam(a1T23_BKstar);
	setparam(a2T23_BKstar);
	/* SN: B -> Kstar FF van Dyk et al. 2305.06301 (BKstar_FormFactor_choice:2)*/
	setparam(a0A0_BKstar_GRvDV_BSZ);
	setparam(a1A0_BKstar_GRvDV_BSZ);
	setparam(a2A0_BKstar_GRvDV_BSZ);
	setparam(a0A1_BKstar_GRvDV_BSZ);
	setparam(a1A1_BKstar_GRvDV_BSZ);
	setparam(a2A1_BKstar_GRvDV_BSZ);
	setparam(a1A12_BKstar_GRvDV_BSZ);
	setparam(a2A12_BKstar_GRvDV_BSZ);
	setparam(a0V_BKstar_GRvDV_BSZ);
	setparam(a1V_BKstar_GRvDV_BSZ);
	setparam(a2V_BKstar_GRvDV_BSZ);
	setparam(a0T1_BKstar_GRvDV_BSZ);
	setparam(a1T1_BKstar_GRvDV_BSZ);
	setparam(a2T1_BKstar_GRvDV_BSZ);
	setparam(a1T2_BKstar_GRvDV_BSZ);
	setparam(a2T2_BKstar_GRvDV_BSZ);
	setparam(a0T23_BKstar_GRvDV_BSZ);
	setparam(a1T23_BKstar_GRvDV_BSZ);
	setparam(a2T23_BKstar_GRvDV_BSZ);
	/* SN: B -> Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:3)*/
	setparam(a0A0_BKstar_GKvD_LCSR_Lattice);
	setparam(a1A0_BKstar_GKvD_LCSR_Lattice);
	setparam(a2A0_BKstar_GKvD_LCSR_Lattice);
	setparam(a0A1_BKstar_GKvD_LCSR_Lattice);
	setparam(a1A1_BKstar_GKvD_LCSR_Lattice);
	setparam(a2A1_BKstar_GKvD_LCSR_Lattice);
	setparam(a1A12_BKstar_GKvD_LCSR_Lattice);
	setparam(a2A12_BKstar_GKvD_LCSR_Lattice);
	setparam(a0V_BKstar_GKvD_LCSR_Lattice);
	setparam(a1V_BKstar_GKvD_LCSR_Lattice);
	setparam(a2V_BKstar_GKvD_LCSR_Lattice);
	setparam(a0T1_BKstar_GKvD_LCSR_Lattice);
	setparam(a1T1_BKstar_GKvD_LCSR_Lattice);
	setparam(a2T1_BKstar_GKvD_LCSR_Lattice);
	setparam(a1T2_BKstar_GKvD_LCSR_Lattice);
	setparam(a2T2_BKstar_GKvD_LCSR_Lattice);
	setparam(a0T23_BKstar_GKvD_LCSR_Lattice);
	setparam(a1T23_BKstar_GKvD_LCSR_Lattice);
	setparam(a2T23_BKstar_GKvD_LCSR_Lattice);
	/* SN: B -> Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:4)*/
	setparam(a0A0_BKstar_GKvD_LCSRonly);
	setparam(a1A0_BKstar_GKvD_LCSRonly);
	setparam(a2A0_BKstar_GKvD_LCSRonly);
	setparam(a0A1_BKstar_GKvD_LCSRonly);
	setparam(a1A1_BKstar_GKvD_LCSRonly);
	setparam(a2A1_BKstar_GKvD_LCSRonly);
	setparam(a1A12_BKstar_GKvD_LCSRonly);
	setparam(a2A12_BKstar_GKvD_LCSRonly);
	setparam(a0V_BKstar_GKvD_LCSRonly);
	setparam(a1V_BKstar_GKvD_LCSRonly);
	setparam(a2V_BKstar_GKvD_LCSRonly);
	setparam(a0T1_BKstar_GKvD_LCSRonly);
	setparam(a1T1_BKstar_GKvD_LCSRonly);
	setparam(a2T1_BKstar_GKvD_LCSRonly);
	setparam(a1T2_BKstar_GKvD_LCSRonly);
	setparam(a2T2_BKstar_GKvD_LCSRonly);
	setparam(a0T23_BKstar_GKvD_LCSRonly);
	setparam(a1T23_BKstar_GKvD_LCSRonly);
	setparam(a2T23_BKstar_GKvD_LCSRonly);
	/* SN: LCSR-only B->Kstar FF Zwicky et al. 1503.05534 (BKstar_FormFactor_choice:5)*/
	setparam(a0A0_BKstar_BSZ_LCSRonly);
	setparam(a1A0_BKstar_BSZ_LCSRonly);
	setparam(a2A0_BKstar_BSZ_LCSRonly);
	setparam(a0A1_BKstar_BSZ_LCSRonly);
	setparam(a1A1_BKstar_BSZ_LCSRonly);
	setparam(a2A1_BKstar_BSZ_LCSRonly);
	setparam(a0A12_BKstar_BSZ_LCSRonly);
	setparam(a1A12_BKstar_BSZ_LCSRonly);
	setparam(a2A12_BKstar_BSZ_LCSRonly);
	setparam(a0V_BKstar_BSZ_LCSRonly);
	setparam(a1V_BKstar_BSZ_LCSRonly);
	setparam(a2V_BKstar_BSZ_LCSRonly);
	setparam(a0T1_BKstar_BSZ_LCSRonly);
	setparam(a1T1_BKstar_BSZ_LCSRonly);
	setparam(a2T1_BKstar_BSZ_LCSRonly);
	setparam(a0T2_BKstar_BSZ_LCSRonly);
	setparam(a1T2_BKstar_BSZ_LCSRonly);
	setparam(a2T2_BKstar_BSZ_LCSRonly);
	setparam(a0T23_BKstar_BSZ_LCSRonly);
	setparam(a1T23_BKstar_BSZ_LCSRonly);
	setparam(a2T23_BKstar_BSZ_LCSRonly);
	/* SN: Lattice-only B->Kstar FF Wingate et al. 1310.3722 + 1501.00367 (BKstar_FormFactor_choice:6*/
	setparam(a0V_BKstar_HLMW);
	setparam(a1V_BKstar_HLMW);
	setparam(a0A0_BKstar_HLMW);
	setparam(a1A0_BKstar_HLMW);
	setparam(a0A1_BKstar_HLMW);
	setparam(a1A1_BKstar_HLMW);
	setparam(a0A12_BKstar_HLMW);
	setparam(a1A12_BKstar_HLMW);
	setparam(a0T1_BKstar_HLMW);
	setparam(a1T1_BKstar_HLMW);
	setparam(a0T2_BKstar_HLMW);
	setparam(a1T2_BKstar_HLMW);
	setparam(a0T23_BKstar_HLMW);
	setparam(a1T23_BKstar_HLMW);
	setparam(systErr_BKstar_HLMW);


/* B -> K */
	setparam(f_K);
	setparam(a1K);
	setparam(a2K);
	/* Form factors B->K ll */
	/* SN: Replacing with more explicit naming - 28/03/2025 */
	// setparam(a00_BK);
	// setparam(a10_BK);
	// setparam(a20_BK);
	// setparam(a0p_BK);
	// setparam(a1p_BK);
	// setparam(a2p_BK);
	// setparam(a0T_BK);
	// setparam(a1T_BK);
	// setparam(a2T_BK);
	setparam(a0f0_BK_AS_LCSR_Lattice);
	setparam(a1f0_BK_AS_LCSR_Lattice);
	setparam(a2f0_BK_AS_LCSR_Lattice);
	setparam(a3f0_BK_AS_LCSR_Lattice);
	setparam(a0fp_BK_AS_LCSR_Lattice);
	setparam(a1fp_BK_AS_LCSR_Lattice);
	setparam(a2fp_BK_AS_LCSR_Lattice);
	setparam(a0fT_BK_AS_LCSR_Lattice);
	setparam(a1fT_BK_AS_LCSR_Lattice);
	setparam(a2fT_BK_AS_LCSR_Lattice);
	/* SN: B->K FF van Dyk et al. 2305.06301 (BK_FormFactor_choice:2)*/
	setparam(a0fp_BK_GRvDV_BSZ);
	setparam(a1fp_BK_GRvDV_BSZ);
	setparam(a2fp_BK_GRvDV_BSZ);
	setparam(a1f0_BK_GRvDV_BSZ);
	setparam(a2f0_BK_GRvDV_BSZ);
	setparam(a0fT_BK_GRvDV_BSZ);
	setparam(a1fT_BK_GRvDV_BSZ);
	setparam(a2fT_BK_GRvDV_BSZ);
	/* SN: LCSR+Lattice B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:3)*/
	setparam(a0fp_BK_GKvD_LCSR_Lattice);
	setparam(a1fp_BK_GKvD_LCSR_Lattice);
	setparam(a2fp_BK_GKvD_LCSR_Lattice);
	setparam(a1f0_BK_GKvD_LCSR_Lattice);
	setparam(a2f0_BK_GKvD_LCSR_Lattice);
	setparam(a0fT_BK_GKvD_LCSR_Lattice);
	setparam(a1fT_BK_GKvD_LCSR_Lattice);
	setparam(a2fT_BK_GKvD_LCSR_Lattice);
	/* SN: LCSR-only B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:4)*/
	setparam(a0fp_BK_GKvD_LCSRonly);
	setparam(a1fp_BK_GKvD_LCSRonly);
	setparam(a2fp_BK_GKvD_LCSRonly);
	setparam(a1f0_BK_GKvD_LCSRonly);
	setparam(a2f0_BK_GKvD_LCSRonly);
	setparam(a0fT_BK_GKvD_LCSRonly);
	setparam(a1fT_BK_GKvD_LCSRonly);
	setparam(a2fT_BK_GKvD_LCSRonly);
	/* SN: Lattice FLAG24 average 2411.04268 of HPQCD'13 and FNAL/MILC'15 (BK_FormFactor_choice:6) - 26/04/2025*/
	setparam(a0fp_BK_FLAG24);
	setparam(a1fp_BK_FLAG24);
	setparam(a2fp_BK_FLAG24);
	setparam(a0f0_BK_FLAG24);
	setparam(a1f0_BK_FLAG24);
	setparam(a0fT_BK_FLAG24);
	setparam(a1fT_BK_FLAG24);
	setparam(a2fT_BK_FLAG24);
	/* SN: Lattice full q2-region HPQCD22 2207.12468 (BK_FormFactor_choice:7) - 26/04/2025*/
	setparam(a0f0_BK_HPQCD22);
	setparam(a1f0_BK_HPQCD22);
	setparam(a2f0_BK_HPQCD22);
	setparam(a1fp_BK_HPQCD22);
	setparam(a2fp_BK_HPQCD22);
	setparam(a0fT_BK_HPQCD22);
	setparam(a1fT_BK_HPQCD22);
	setparam(a2fT_BK_HPQCD22);
	setparam(L_BK_HPQCD22);
	/* SN: LCSR-only Khodjamirian & Rusov 1703.04765 (BK_FormFactor_choice:8) - 02/08/2025*/
	setparam(f0fp_BK_KR);
	setparam(b1fp_BK_KR);
	setparam(f0fT_BK_KR);
	setparam(b1fT_BK_KR);


	/* low */
	setparam(BtoKlow_FV_err_noq2);
	setparam(BtoKlow_FA_err_noq2);
	setparam(BtoKlow_FS_err_noq2);
	setparam(BtoKlow_FP_err_noq2);

	setparam(BtoKlow_FV_err_q2);
	setparam(BtoKlow_FA_err_q2);
	setparam(BtoKlow_FS_err_q2);
	setparam(BtoKlow_FP_err_q2);

	/* high */
	setparam(BtoKhigh_FV_err);
	setparam(BtoKhigh_FA_err);
	setparam(BtoKhigh_FS_err);
	setparam(BtoKhigh_FP_err);

/* Bs */
	setparam(life_Bs);
	setparam(f_Bs);
	setparam(lambda_Bsp);
	setparam(ys_Bs);
	
/* Bs -> phi */
	setparam(f_phi_par);
	setparam(f_phi_perp);
	setparam(a1phi_perp);
	setparam(a1phi_par);
	setparam(a2phi_perp);
	setparam(a2phi_par);

	/* low */
	setparam(Bstophilow_ALperp_err_noq2);
	setparam(Bstophilow_ARperp_err_noq2);
	setparam(Bstophilow_ALpar_err_noq2);
	setparam(Bstophilow_ARpar_err_noq2);
	setparam(Bstophilow_AL0_err_noq2);
	setparam(Bstophilow_AR0_err_noq2);
	setparam(Bstophilow_At_err_noq2);
	setparam(Bstophilow_AS_err_noq2);

	setparam(Bstophilow_ALperp_err_q2);
	setparam(Bstophilow_ARperp_err_q2);
	setparam(Bstophilow_ALpar_err_q2);
	setparam(Bstophilow_ARpar_err_q2);
	setparam(Bstophilow_AL0_err_q2);
	setparam(Bstophilow_AR0_err_q2);
	setparam(Bstophilow_At_err_q2);
	setparam(Bstophilow_AS_err_q2);
	
	/* high */
	setparam(Bstophihigh_ALperp_err);
	setparam(Bstophihigh_ARperp_err);
	setparam(Bstophihigh_ALpar_err);
	setparam(Bstophihigh_ARpar_err);
	setparam(Bstophihigh_AL0_err);
	setparam(Bstophihigh_AR0_err);
	setparam(Bstophihigh_At_err);
	setparam(Bstophihigh_AS_err);
			
	/* Form factors Bs->phi ll */	
    setparam(life_Lb);
    setparam(alphaL_LbLll);
	setparam(a0A0_Bsphi);
	setparam(a1A0_Bsphi);
	setparam(a2A0_Bsphi);
	setparam(a0A1_Bsphi);
	setparam(a1A1_Bsphi);
	setparam(a2A1_Bsphi);
	setparam(a0A12_Bsphi);
	setparam(a1A12_Bsphi);
	setparam(a2A12_Bsphi);
	setparam(a0V_Bsphi);
	setparam(a1V_Bsphi);
	setparam(a2V_Bsphi);
	setparam(a0T1_Bsphi);
	setparam(a1T1_Bsphi);
	setparam(a2T1_Bsphi);
	setparam(a0T2_Bsphi);
	setparam(a1T2_Bsphi);
	setparam(a2T2_Bsphi);
	setparam(a0T23_Bsphi);
	setparam(a1T23_Bsphi);
	setparam(a2T23_Bsphi);
	/* SN: Bs->phi FF van Dyk et al. 2305.06301 (Bsphi_FormFactor_choice:2)*/
	setparam(a0A0_Bsphi_GRvDV_BSZ);
	setparam(a1A0_Bsphi_GRvDV_BSZ);
	setparam(a2A0_Bsphi_GRvDV_BSZ);
	setparam(a0A1_Bsphi_GRvDV_BSZ);
	setparam(a1A1_Bsphi_GRvDV_BSZ);
	setparam(a2A1_Bsphi_GRvDV_BSZ);
	setparam(a1A12_Bsphi_GRvDV_BSZ);
	setparam(a2A12_Bsphi_GRvDV_BSZ);
	setparam(a0V_Bsphi_GRvDV_BSZ);
	setparam(a1V_Bsphi_GRvDV_BSZ);
	setparam(a2V_Bsphi_GRvDV_BSZ);
	setparam(a0T1_Bsphi_GRvDV_BSZ);
	setparam(a1T1_Bsphi_GRvDV_BSZ);
	setparam(a2T1_Bsphi_GRvDV_BSZ);
	setparam(a1T2_Bsphi_GRvDV_BSZ);
	setparam(a2T2_Bsphi_GRvDV_BSZ);
	setparam(a0T23_Bsphi_GRvDV_BSZ);
	setparam(a1T23_Bsphi_GRvDV_BSZ);
	setparam(a2T23_Bsphi_GRvDV_BSZ);
	/* SN: Bs->phi FF Zwicky et al. 1503.05534 (Bsphi_FormFactor_choice:5)*/
	setparam(a0A0_Bsphi_BSZ_LCSRonly);
	setparam(a1A0_Bsphi_BSZ_LCSRonly);
	setparam(a2A0_Bsphi_BSZ_LCSRonly);
	setparam(a0A1_Bsphi_BSZ_LCSRonly);
	setparam(a1A1_Bsphi_BSZ_LCSRonly);
	setparam(a2A1_Bsphi_BSZ_LCSRonly);
	setparam(a0A12_Bsphi_BSZ_LCSRonly);
	setparam(a1A12_Bsphi_BSZ_LCSRonly);
	setparam(a2A12_Bsphi_BSZ_LCSRonly);
	setparam(a0V_Bsphi_BSZ_LCSRonly);
	setparam(a1V_Bsphi_BSZ_LCSRonly);
	setparam(a2V_Bsphi_BSZ_LCSRonly);
	setparam(a0T1_Bsphi_BSZ_LCSRonly);
	setparam(a1T1_Bsphi_BSZ_LCSRonly);
	setparam(a2T1_Bsphi_BSZ_LCSRonly);
	setparam(a0T2_Bsphi_BSZ_LCSRonly);
	setparam(a1T2_Bsphi_BSZ_LCSRonly);
	setparam(a2T2_Bsphi_BSZ_LCSRonly);
	setparam(a0T23_Bsphi_BSZ_LCSRonly);
	setparam(a1T23_Bsphi_BSZ_LCSRonly);
	setparam(a2T23_Bsphi_BSZ_LCSRonly);
	/* SN: Bs->phi FF Wingate et al. 1310.3722 + 1501.00367 (Bsphi_FormFactor_choice:6)*/
	setparam(a0V_Bsphi_HLMW);
	setparam(a1V_Bsphi_HLMW);
	setparam(a0A0_Bsphi_HLMW);
	setparam(a1A0_Bsphi_HLMW);
	setparam(a0A1_Bsphi_HLMW);
	setparam(a1A1_Bsphi_HLMW);
	setparam(a0A12_Bsphi_HLMW);
	setparam(a1A12_Bsphi_HLMW);
	setparam(a0T1_Bsphi_HLMW);
	setparam(a1T1_Bsphi_HLMW);
	setparam(a0T2_Bsphi_HLMW);
	setparam(a1T2_Bsphi_HLMW);
	setparam(a0T23_Bsphi_HLMW);
	setparam(a1T23_Bsphi_HLMW);
	setparam(systErr_Bsphi_HLMW);

	/* Lambda_b -> Lambda l+l- */	
	setparam(life_Lb);
    setparam(alphaL_LbLll);
    setparam(a0_HO_fplus_LbLll);
    setparam(a1_HO_fplus_LbLll);
    setparam(a2_HO_fplus_LbLll);
    setparam(a0_HO_fperp_LbLll);
    setparam(a1_HO_fperp_LbLll);
    setparam(a2_HO_fperp_LbLll);
    setparam(a0_HO_gpp_LbLll);
    setparam(a1_HO_gplus_LbLll);
    setparam(a2_HO_gplus_LbLll);
    setparam(a1_HO_gperp_LbLll);
    setparam(a2_HO_gperp_LbLll);
    setparam(a0_HO_hplus_LbLll);
    setparam(a1_HO_hplus_LbLll);
    setparam(a2_HO_hplus_LbLll);
    setparam(a0_HO_hperp_LbLll);
    setparam(a1_HO_hperp_LbLll);
    setparam(a2_HO_hperp_LbLll);
    setparam(a0_HO_htildepp_LbLll);
    setparam(a1_HO_htildeplus_LbLll);
    setparam(a2_HO_htildeplus_LbLll);
    setparam(a1_HO_htildeperp_LbLll);
    setparam(a2_HO_htildeperp_LbLll);
    
	/* Kaons */
	setparam(deltaPcu_Kppipnunu);
	setparam(err_Pc_Xlambda_Kppipnunu);
	setparam(BR_KLgammagamma_exp);
	setparam(Aterm_mu_KLmumu);
	setparam(chi_gg_Mrho);
	setparam(BR_KSgammagamma_exp);
	setparam(Iterm_mu_KSmumu);
	setparam(log_muK_1GeV);
	setparam(KLpill_Ce_Cdir);
	setparam(KLpill_Cmu_Cdir);
	setparam(KLpill_Ce_Cint);
	setparam(KLpill_Cmu_Cint);
	setparam(KLpill_Ce_Cmix);
	setparam(KLpill_Cmu_Cmix);
	setparam(KLpill_Cmu_CPC);
	setparam(KLpill_abs_aS);
	/* Scalar and Pseudoscalar contributions included - 08/10/2025 */
	setparam(KLpill_Ce_Sgg);
    setparam(KLpill_Cmu_Sgg);
	/* AFBmu_KLpi0mumu  - 09/10/2025 */
    setparam(KLpill_AFB_a1pi_gg_7V);
	setparam(KLpill_AFB_a1pi_gg_aS);

	#undef setparam


	return;
}

/*--------------------------------------------------------------------*/

void write_nuisance(struct nuisance* nuisparam, char name[])
/* Write the nuisparam structure into file "name" */
{
	FILE *file;
	
	file=fopen(name,"w");

#define writeparam( myparam ) fprintf(file,"%-30s\t%.4e\t\t%.4e\t\t%d\n",nuisparam->myparam.name,nuisparam->myparam.cent,nuisparam->myparam.dev,nuisparam->myparam.type);
	
	writeparam(alphas_MZ);

	writeparam(mass_b);
	writeparam(mass_c);
	writeparam(mass_s);
	writeparam(mass_top_pole);

	writeparam(mass_h0);

	writeparam(CKM_lambda);
	writeparam(CKM_A);
	writeparam(CKM_rhobar);
	writeparam(CKM_etabar);

	writeparam(log_mu_W_mass_W);
	writeparam(log_mu_b_mass_b);

/* inclusive b -> s */
	writeparam(BR_BXclnu_exp);
/* b -> s gamma */ 
	writeparam(mu_G2_bsg);
	writeparam(rho_D3_bsg);
	writeparam(rho_LS3_bsg);
	writeparam(bsgamma_rand);
	writeparam(mu_c_bsg);
/* b -> s mu mu */
	writeparam(BRBXsmumu_lowq2_rand);
	writeparam(BRBXsmumu_highq2_rand);
	writeparam(BRBXsmumu_full_rand);
/* b -> s e e */
	writeparam(BRBXsee_lowq2_rand);
	writeparam(BRBXsee_highq2_rand);
	writeparam(BRBXsee_full_rand);
/* b -> s tau tau */
	writeparam(BRBXstautau_lowq2_rand);
	writeparam(BRBXstautau_highq2_rand);
	writeparam(BRBXstautau_full_rand);

/* B */	
	writeparam(f_B);
	writeparam(lambda_Bp);
/* B -> K* */
	writeparam(f_Kstar_par);
	writeparam(f_Kstar_perp);
	writeparam(a1perp);
	writeparam(a2perp);
	writeparam(a1par);
	writeparam(a2par);
/* B -> K* gamma */
	writeparam(T1_BKstar);
	writeparam(log_mu_spec_lambda_h_mass_b);
	/* low */
	writeparam(BtoKstarlow_ALperp_err_noq2);
	writeparam(BtoKstarlow_ARperp_err_noq2);
	writeparam(BtoKstarlow_ALpar_err_noq2);
	writeparam(BtoKstarlow_ARpar_err_noq2);
	writeparam(BtoKstarlow_AL0_err_noq2);
	writeparam(BtoKstarlow_AR0_err_noq2);
	writeparam(BtoKstarlow_At_err_noq2);
	writeparam(BtoKstarlow_AS_err_noq2);

	writeparam(BtoKstarlow_ALperp_err_q2);
	writeparam(BtoKstarlow_ARperp_err_q2);
	writeparam(BtoKstarlow_ALpar_err_q2);
	writeparam(BtoKstarlow_ARpar_err_q2);
	writeparam(BtoKstarlow_AL0_err_q2);
	writeparam(BtoKstarlow_AR0_err_q2);
	writeparam(BtoKstarlow_At_err_q2);
	writeparam(BtoKstarlow_AS_err_q2);

	/* SN: Result from B->K* Hadronic Fit (to be used to check how it impacts predictions) - 13/04/2025 */
	writeparam(hplus0r);
	writeparam(hminus0r);
	writeparam(hplus1r);
	writeparam(hminus1r);
	writeparam(hplus2r);
	writeparam(hminus2r);
	writeparam(hzero0r);
	writeparam(hzero1r);
	writeparam(hzero2r);
	writeparam(hplus0i);
	writeparam(hminus0i);
	writeparam(hplus1i);
	writeparam(hminus1i);
	writeparam(hplus2i);
	writeparam(hminus2i);
	writeparam(hzero0i);
	writeparam(hzero1i);
	writeparam(hzero2i);

	writeparam(real_alpha_perp0);
	writeparam(real_alpha_perp1);
	writeparam(real_alpha_perp2);
	writeparam(real_alpha_par0);
	writeparam(real_alpha_par1);
	writeparam(real_alpha_par2);
	writeparam(real_alpha_zero0);
	writeparam(real_alpha_zero1);
	writeparam(imag_alpha_perp0);
	writeparam(imag_alpha_perp1);
	writeparam(imag_alpha_perp2);
	writeparam(imag_alpha_par0);
	writeparam(imag_alpha_par1);
	writeparam(imag_alpha_par2);
	writeparam(imag_alpha_zero0);
	writeparam(imag_alpha_zero1);
	
	writeparam(DeltaC9_M1_q2bar);
	writeparam(r1_M1);
	writeparam(r2_M1);
	writeparam(DeltaC9_M2_q2bar);
	writeparam(r1_M2);
	writeparam(r2_M2);
	writeparam(DeltaC9_M3_q2bar);
	writeparam(r1_M3);
	writeparam(r2_M3);
	
	/* high */
	writeparam(BtoKstarhigh_ALperp_err);
	writeparam(BtoKstarhigh_ARperp_err);
	writeparam(BtoKstarhigh_ALpar_err);
	writeparam(BtoKstarhigh_ARpar_err);
	writeparam(BtoKstarhigh_AL0_err);
	writeparam(BtoKstarhigh_AR0_err);
	writeparam(BtoKstarhigh_At_err);
	writeparam(BtoKstarhigh_AS_err);

	/* Form factors B->K* ll */	
	writeparam(a0A0_BKstar);
	writeparam(a1A0_BKstar);
	writeparam(a2A0_BKstar);
	writeparam(a0A1_BKstar);
	writeparam(a1A1_BKstar);
	writeparam(a2A1_BKstar);
	writeparam(a0A12_BKstar);
	writeparam(a1A12_BKstar);
	writeparam(a2A12_BKstar);
	writeparam(a0V_BKstar);
	writeparam(a1V_BKstar);
	writeparam(a2V_BKstar);
	writeparam(a0T1_BKstar);
	writeparam(a1T1_BKstar);
	writeparam(a2T1_BKstar);
	writeparam(a0T2_BKstar);
	writeparam(a1T2_BKstar);
	writeparam(a2T2_BKstar);
	writeparam(a0T23_BKstar);
	writeparam(a1T23_BKstar);
	writeparam(a2T23_BKstar);

	/* SN: B->Kstar FF van Dyk et al. 2305.06301 (BKstar_FormFactor_choice:2)*/
	writeparam(a0A0_BKstar_GRvDV_BSZ);
	writeparam(a1A0_BKstar_GRvDV_BSZ);
	writeparam(a2A0_BKstar_GRvDV_BSZ);
	writeparam(a0A1_BKstar_GRvDV_BSZ);
	writeparam(a1A1_BKstar_GRvDV_BSZ);
	writeparam(a2A1_BKstar_GRvDV_BSZ);
	writeparam(a1A12_BKstar_GRvDV_BSZ);
	writeparam(a2A12_BKstar_GRvDV_BSZ);
	writeparam(a0V_BKstar_GRvDV_BSZ);
	writeparam(a1V_BKstar_GRvDV_BSZ);
	writeparam(a2V_BKstar_GRvDV_BSZ);
	writeparam(a0T1_BKstar_GRvDV_BSZ);
	writeparam(a1T1_BKstar_GRvDV_BSZ);
	writeparam(a2T1_BKstar_GRvDV_BSZ);
	writeparam(a1T2_BKstar_GRvDV_BSZ);
	writeparam(a2T2_BKstar_GRvDV_BSZ);
	writeparam(a0T23_BKstar_GRvDV_BSZ);
	writeparam(a1T23_BKstar_GRvDV_BSZ);
	writeparam(a2T23_BKstar_GRvDV_BSZ);

	/* SN: LCSR+Lattice B->Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:3)*/
	writeparam(a0A0_BKstar_GKvD_LCSR_Lattice);
	writeparam(a1A0_BKstar_GKvD_LCSR_Lattice);
	writeparam(a2A0_BKstar_GKvD_LCSR_Lattice);
	writeparam(a0A1_BKstar_GKvD_LCSR_Lattice);
	writeparam(a1A1_BKstar_GKvD_LCSR_Lattice);
	writeparam(a2A1_BKstar_GKvD_LCSR_Lattice);
	writeparam(a1A12_BKstar_GKvD_LCSR_Lattice);
	writeparam(a2A12_BKstar_GKvD_LCSR_Lattice);
	writeparam(a0V_BKstar_GKvD_LCSR_Lattice);
	writeparam(a1V_BKstar_GKvD_LCSR_Lattice);
	writeparam(a2V_BKstar_GKvD_LCSR_Lattice);
	writeparam(a0T1_BKstar_GKvD_LCSR_Lattice);
	writeparam(a1T1_BKstar_GKvD_LCSR_Lattice);
	writeparam(a2T1_BKstar_GKvD_LCSR_Lattice);
	// writeparam(a0T2_BKstar_GKvD_LCSR_Lattice);
	writeparam(a1T2_BKstar_GKvD_LCSR_Lattice);
	writeparam(a2T2_BKstar_GKvD_LCSR_Lattice);
	writeparam(a0T23_BKstar_GKvD_LCSR_Lattice);
	writeparam(a1T23_BKstar_GKvD_LCSR_Lattice);
	writeparam(a2T23_BKstar_GKvD_LCSR_Lattice);

	/* SN: LCSR-only B->Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:4)*/
	writeparam(a0A0_BKstar_GKvD_LCSRonly);
	writeparam(a1A0_BKstar_GKvD_LCSRonly);
	writeparam(a2A0_BKstar_GKvD_LCSRonly);
	writeparam(a0A1_BKstar_GKvD_LCSRonly);
	writeparam(a1A1_BKstar_GKvD_LCSRonly);
	writeparam(a2A1_BKstar_GKvD_LCSRonly);
	writeparam(a1A12_BKstar_GKvD_LCSRonly);
	writeparam(a2A12_BKstar_GKvD_LCSRonly);
	writeparam(a0V_BKstar_GKvD_LCSRonly);
	writeparam(a1V_BKstar_GKvD_LCSRonly);
	writeparam(a2V_BKstar_GKvD_LCSRonly);
	writeparam(a0T1_BKstar_GKvD_LCSRonly);
	writeparam(a1T1_BKstar_GKvD_LCSRonly);
	writeparam(a2T1_BKstar_GKvD_LCSRonly);
	// writeparam(a0T2_BKstar_GKvD_LCSRonly);
	writeparam(a1T2_BKstar_GKvD_LCSRonly);
	writeparam(a2T2_BKstar_GKvD_LCSRonly);
	writeparam(a0T23_BKstar_GKvD_LCSRonly);
	writeparam(a1T23_BKstar_GKvD_LCSRonly);
	writeparam(a2T23_BKstar_GKvD_LCSRonly);

	/* SN: LCSR-only B->Kstar FF Zwicky et al. 1503.05534 (BKstar_FormFactor_choice:5)*/
	writeparam(a0A0_BKstar_BSZ_LCSRonly);
	writeparam(a1A0_BKstar_BSZ_LCSRonly);
	writeparam(a2A0_BKstar_BSZ_LCSRonly);
	writeparam(a0A1_BKstar_BSZ_LCSRonly);
	writeparam(a1A1_BKstar_BSZ_LCSRonly);
	writeparam(a2A1_BKstar_BSZ_LCSRonly);
	writeparam(a0A12_BKstar_BSZ_LCSRonly);
	writeparam(a1A12_BKstar_BSZ_LCSRonly);
	writeparam(a2A12_BKstar_BSZ_LCSRonly);
	writeparam(a0V_BKstar_BSZ_LCSRonly);
	writeparam(a1V_BKstar_BSZ_LCSRonly);
	writeparam(a2V_BKstar_BSZ_LCSRonly);
	writeparam(a0T1_BKstar_BSZ_LCSRonly);
	writeparam(a1T1_BKstar_BSZ_LCSRonly);
	writeparam(a2T1_BKstar_BSZ_LCSRonly);
	writeparam(a0T2_BKstar_BSZ_LCSRonly);
	writeparam(a1T2_BKstar_BSZ_LCSRonly);
	writeparam(a2T2_BKstar_BSZ_LCSRonly);
	writeparam(a0T23_BKstar_BSZ_LCSRonly);
	writeparam(a1T23_BKstar_BSZ_LCSRonly);
	writeparam(a2T23_BKstar_BSZ_LCSRonly);

	/* SN: Lattice-only B->Kstar FF Wingate et al. 1310.3722 + 1501.00367 (BKstar_FormFactor_choice:6*/
	writeparam(a0V_BKstar_HLMW);
	writeparam(a1V_BKstar_HLMW);
	writeparam(a0A0_BKstar_HLMW);
	writeparam(a1A0_BKstar_HLMW);
	writeparam(a0A1_BKstar_HLMW);
	writeparam(a1A1_BKstar_HLMW);
	writeparam(a0A12_BKstar_HLMW);
	writeparam(a1A12_BKstar_HLMW);
	writeparam(a0T1_BKstar_HLMW);
	writeparam(a1T1_BKstar_HLMW);
	writeparam(a0T2_BKstar_HLMW);
	writeparam(a1T2_BKstar_HLMW);
	writeparam(a0T23_BKstar_HLMW);
	writeparam(a1T23_BKstar_HLMW);
	writeparam(systErr_BKstar_HLMW);

	/* B -> K */
	writeparam(f_K);
	writeparam(a1K);
	writeparam(a2K);
	/* Form factors B->K ll */
	/* SN: Replacing with more explicit naming and adding missing a3f0 - 28/03/2025 */
	// writeparam(a00_BK);
	// writeparam(a10_BK);
	// writeparam(a20_BK);
	// writeparam(a0p_BK);
	// writeparam(a1p_BK);
	// writeparam(a2p_BK);
	// writeparam(a0T_BK);
	// writeparam(a1T_BK);
	// writeparam(a2T_BK);
	/* SN: B->K FF Altmannshofer, Straub 1411.3161 (BK_FormFactor_choice:1)*/
	writeparam(a0f0_BK_AS_LCSR_Lattice);
	writeparam(a1f0_BK_AS_LCSR_Lattice);
	writeparam(a2f0_BK_AS_LCSR_Lattice);
	writeparam(a3f0_BK_AS_LCSR_Lattice);
	writeparam(a0fp_BK_AS_LCSR_Lattice);
	writeparam(a1fp_BK_AS_LCSR_Lattice);
	writeparam(a2fp_BK_AS_LCSR_Lattice);
	writeparam(a0fT_BK_AS_LCSR_Lattice);
	writeparam(a1fT_BK_AS_LCSR_Lattice);
	writeparam(a2fT_BK_AS_LCSR_Lattice);
	/* SN: B->K FF van Dyk et al. 2305.06301 (BK_FormFactor_choice:2)*/
	writeparam(a0fp_BK_GRvDV_BSZ);
	writeparam(a1fp_BK_GRvDV_BSZ);
	writeparam(a2fp_BK_GRvDV_BSZ);
	writeparam(a1f0_BK_GRvDV_BSZ);
	writeparam(a2f0_BK_GRvDV_BSZ);
	writeparam(a0fT_BK_GRvDV_BSZ);
	writeparam(a1fT_BK_GRvDV_BSZ);
	writeparam(a2fT_BK_GRvDV_BSZ);
	/* SN: LCSR+Lattice B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:3)*/
	writeparam(a0fp_BK_GKvD_LCSR_Lattice);
	writeparam(a1fp_BK_GKvD_LCSR_Lattice);
	writeparam(a2fp_BK_GKvD_LCSR_Lattice);
	writeparam(a1f0_BK_GKvD_LCSR_Lattice);
	writeparam(a2f0_BK_GKvD_LCSR_Lattice);
	writeparam(a0fT_BK_GKvD_LCSR_Lattice);
	writeparam(a1fT_BK_GKvD_LCSR_Lattice);
	writeparam(a2fT_BK_GKvD_LCSR_Lattice);
	/* SN: LCSR-only B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:4)*/
	writeparam(a0fp_BK_GKvD_LCSRonly);
	writeparam(a1fp_BK_GKvD_LCSRonly);
	writeparam(a2fp_BK_GKvD_LCSRonly);
	writeparam(a1f0_BK_GKvD_LCSRonly);
	writeparam(a2f0_BK_GKvD_LCSRonly);
	writeparam(a0fT_BK_GKvD_LCSRonly);
	writeparam(a1fT_BK_GKvD_LCSRonly);
	writeparam(a2fT_BK_GKvD_LCSRonly);
	/* SN: Lattice FLAG24 average 2411.04268 of HPQCD'13 and FNAL/MILC'15 (BK_FormFactor_choice:6) - 26/04/2025*/
	writeparam(a0fp_BK_FLAG24);
	writeparam(a1fp_BK_FLAG24);
	writeparam(a2fp_BK_FLAG24);
	writeparam(a0f0_BK_FLAG24);
	writeparam(a1f0_BK_FLAG24);
	writeparam(a0fT_BK_FLAG24);
	writeparam(a1fT_BK_FLAG24);
	writeparam(a2fT_BK_FLAG24);
	/* SN: Lattice full q2-region HPQCD22 2207.12468 (BK_FormFactor_choice:7) - 26/04/2025*/
	writeparam(a0f0_BK_HPQCD22);
	writeparam(a1f0_BK_HPQCD22);
	writeparam(a2f0_BK_HPQCD22);
	writeparam(a1fp_BK_HPQCD22);
	writeparam(a2fp_BK_HPQCD22);
	writeparam(a0fT_BK_HPQCD22);
	writeparam(a1fT_BK_HPQCD22);
	writeparam(a2fT_BK_HPQCD22);
	writeparam(L_BK_HPQCD22);
	/* SN: LCSR-only Khodjamirian & Rusov 1703.04765 (BK_FormFactor_choice:8) - 02/08/2025*/
	writeparam(f0fp_BK_KR);
	writeparam(b1fp_BK_KR);
	writeparam(f0fT_BK_KR);
	writeparam(b1fT_BK_KR);

	/* low */
	writeparam(BtoKlow_FV_err_noq2);
	writeparam(BtoKlow_FA_err_noq2);
	writeparam(BtoKlow_FS_err_noq2);
	writeparam(BtoKlow_FP_err_noq2);

	writeparam(BtoKlow_FV_err_q2);
	writeparam(BtoKlow_FA_err_q2);
	writeparam(BtoKlow_FS_err_q2);
	writeparam(BtoKlow_FP_err_q2);

	/* high */
	writeparam(BtoKhigh_FV_err);
	writeparam(BtoKhigh_FA_err);
	writeparam(BtoKhigh_FS_err);
	writeparam(BtoKhigh_FP_err);

/* Bs */
	writeparam(life_Bs);
	writeparam(f_Bs);
	writeparam(lambda_Bsp);
	writeparam(ys_Bs);
	
/* Bs -> phi */
	writeparam(f_phi_par);
	writeparam(f_phi_perp);
	writeparam(a1phi_perp);
	writeparam(a1phi_par);
	writeparam(a2phi_perp);
	writeparam(a2phi_par);

	/* low */
	writeparam(Bstophilow_ALperp_err_noq2);
	writeparam(Bstophilow_ARperp_err_noq2);
	writeparam(Bstophilow_ALpar_err_noq2);
	writeparam(Bstophilow_ARpar_err_noq2);
	writeparam(Bstophilow_AL0_err_noq2);
	writeparam(Bstophilow_AR0_err_noq2);
	writeparam(Bstophilow_At_err_noq2);
	writeparam(Bstophilow_AS_err_noq2);

	writeparam(Bstophilow_ALperp_err_q2);
	writeparam(Bstophilow_ARperp_err_q2);
	writeparam(Bstophilow_ALpar_err_q2);
	writeparam(Bstophilow_ARpar_err_q2);
	writeparam(Bstophilow_AL0_err_q2);
	writeparam(Bstophilow_AR0_err_q2);
	writeparam(Bstophilow_At_err_q2);
	writeparam(Bstophilow_AS_err_q2);
	
	/* high */
	writeparam(Bstophihigh_ALperp_err);
	writeparam(Bstophihigh_ARperp_err);
	writeparam(Bstophihigh_ALpar_err);
	writeparam(Bstophihigh_ARpar_err);
	writeparam(Bstophihigh_AL0_err);
	writeparam(Bstophihigh_AR0_err);
	writeparam(Bstophihigh_At_err);
	writeparam(Bstophihigh_AS_err);
			
	/* Form factors Bs->phi ll */	
	writeparam(a0A0_Bsphi);
	writeparam(a1A0_Bsphi);
	writeparam(a2A0_Bsphi);
	writeparam(a0A1_Bsphi);
	writeparam(a1A1_Bsphi);
	writeparam(a2A1_Bsphi);
	writeparam(a0A12_Bsphi);
	writeparam(a1A12_Bsphi);
	writeparam(a2A12_Bsphi);
	writeparam(a0V_Bsphi);
	writeparam(a1V_Bsphi);
	writeparam(a2V_Bsphi);
	writeparam(a0T1_Bsphi);
	writeparam(a1T1_Bsphi);
	writeparam(a2T1_Bsphi);
	writeparam(a0T2_Bsphi);
	writeparam(a1T2_Bsphi);
	writeparam(a2T2_Bsphi);
	writeparam(a0T23_Bsphi);
	writeparam(a1T23_Bsphi);
	writeparam(a2T23_Bsphi);

	/* SN: Bs->phi FF van Dyk et al. 2305.06301 (Bsphi_FormFactor_choice:2)*/
	writeparam(a0A0_Bsphi_GRvDV_BSZ);
	writeparam(a1A0_Bsphi_GRvDV_BSZ);
	writeparam(a2A0_Bsphi_GRvDV_BSZ);
	writeparam(a0A1_Bsphi_GRvDV_BSZ);
	writeparam(a1A1_Bsphi_GRvDV_BSZ);
	writeparam(a2A1_Bsphi_GRvDV_BSZ);
	writeparam(a1A12_Bsphi_GRvDV_BSZ);
	writeparam(a2A12_Bsphi_GRvDV_BSZ);
	writeparam(a0V_Bsphi_GRvDV_BSZ);
	writeparam(a1V_Bsphi_GRvDV_BSZ);
	writeparam(a2V_Bsphi_GRvDV_BSZ);
	writeparam(a0T1_Bsphi_GRvDV_BSZ);
	writeparam(a1T1_Bsphi_GRvDV_BSZ);
	writeparam(a2T1_Bsphi_GRvDV_BSZ);
	writeparam(a1T2_Bsphi_GRvDV_BSZ);
	writeparam(a2T2_Bsphi_GRvDV_BSZ);
	writeparam(a0T23_Bsphi_GRvDV_BSZ);
	writeparam(a1T23_Bsphi_GRvDV_BSZ);
	writeparam(a2T23_Bsphi_GRvDV_BSZ);

	/* SN: Bs->phi FF Zwicky et al. 1503.05534 (Bsphi_FormFactor_choice:5)*/
	writeparam(a0A0_Bsphi_BSZ_LCSRonly);
	writeparam(a1A0_Bsphi_BSZ_LCSRonly);
	writeparam(a2A0_Bsphi_BSZ_LCSRonly);
	writeparam(a0A1_Bsphi_BSZ_LCSRonly);
	writeparam(a1A1_Bsphi_BSZ_LCSRonly);
	writeparam(a2A1_Bsphi_BSZ_LCSRonly);
	writeparam(a0A12_Bsphi_BSZ_LCSRonly);
	writeparam(a1A12_Bsphi_BSZ_LCSRonly);
	writeparam(a2A12_Bsphi_BSZ_LCSRonly);
	writeparam(a0V_Bsphi_BSZ_LCSRonly);
	writeparam(a1V_Bsphi_BSZ_LCSRonly);
	writeparam(a2V_Bsphi_BSZ_LCSRonly);
	writeparam(a0T1_Bsphi_BSZ_LCSRonly);
	writeparam(a1T1_Bsphi_BSZ_LCSRonly);
	writeparam(a2T1_Bsphi_BSZ_LCSRonly);
	writeparam(a0T2_Bsphi_BSZ_LCSRonly);
	writeparam(a1T2_Bsphi_BSZ_LCSRonly);
	writeparam(a2T2_Bsphi_BSZ_LCSRonly);
	writeparam(a0T23_Bsphi_BSZ_LCSRonly);
	writeparam(a1T23_Bsphi_BSZ_LCSRonly);
	writeparam(a2T23_Bsphi_BSZ_LCSRonly);

	/* SN: Bs->phi FF Wingate et al. 1310.3722 + 1501.00367 (Bsphi_FormFactor_choice:6)*/
	writeparam(a0V_Bsphi_HLMW);
	writeparam(a1V_Bsphi_HLMW);
	writeparam(a0A0_Bsphi_HLMW);
	writeparam(a1A0_Bsphi_HLMW);
	writeparam(a0A1_Bsphi_HLMW);
	writeparam(a1A1_Bsphi_HLMW);
	writeparam(a0A12_Bsphi_HLMW);
	writeparam(a1A12_Bsphi_HLMW);
	writeparam(a0T1_Bsphi_HLMW);
	writeparam(a1T1_Bsphi_HLMW);
	writeparam(a0T2_Bsphi_HLMW);
	writeparam(a1T2_Bsphi_HLMW);
	writeparam(a0T23_Bsphi_HLMW);
	writeparam(a1T23_Bsphi_HLMW);
	writeparam(systErr_Bsphi_HLMW);


	/* Lambda_b -> Lambda l+l- */	
	writeparam(life_Lb);
    writeparam(alphaL_LbLll);
    writeparam(a0_HO_fplus_LbLll);
    writeparam(a1_HO_fplus_LbLll);
    writeparam(a2_HO_fplus_LbLll);
    writeparam(a0_HO_fperp_LbLll);
    writeparam(a1_HO_fperp_LbLll);
    writeparam(a2_HO_fperp_LbLll);
    writeparam(a0_HO_gpp_LbLll);
    writeparam(a1_HO_gplus_LbLll);
    writeparam(a2_HO_gplus_LbLll);
    writeparam(a1_HO_gperp_LbLll);
    writeparam(a2_HO_gperp_LbLll);
    writeparam(a0_HO_hplus_LbLll);
    writeparam(a1_HO_hplus_LbLll);
    writeparam(a2_HO_hplus_LbLll);
    writeparam(a0_HO_hperp_LbLll);
    writeparam(a1_HO_hperp_LbLll);
    writeparam(a2_HO_hperp_LbLll);
    writeparam(a0_HO_htildepp_LbLll);
    writeparam(a1_HO_htildeplus_LbLll);
    writeparam(a2_HO_htildeplus_LbLll);
    writeparam(a1_HO_htildeperp_LbLll);
    writeparam(a2_HO_htildeperp_LbLll);

	/* Kaons */
	writeparam(deltaPcu_Kppipnunu);
	writeparam(err_Pc_Xlambda_Kppipnunu);
	writeparam(BR_KLgammagamma_exp);
	writeparam(Aterm_mu_KLmumu);
	writeparam(chi_gg_Mrho);
	writeparam(BR_KSgammagamma_exp);
	writeparam(Iterm_mu_KSmumu);
	writeparam(log_muK_1GeV);
	writeparam(KLpill_Ce_Cdir);
	writeparam(KLpill_Cmu_Cdir);
	writeparam(KLpill_Ce_Cint);
	writeparam(KLpill_Cmu_Cint);
	writeparam(KLpill_Ce_Cmix);
	writeparam(KLpill_Cmu_Cmix);
	writeparam(KLpill_Cmu_CPC);
	writeparam(KLpill_abs_aS);
	/* Scalar and Pseudoscalar contributions included - 08/10/2025 */
	writeparam(KLpill_Ce_Sgg);
    writeparam(KLpill_Cmu_Sgg);
	/* AFBmu_KLpi0mumu  - 09/10/2025 */
    writeparam(KLpill_AFB_a1pi_gg_7V);
	writeparam(KLpill_AFB_a1pi_gg_aS);

	fclose(file);

#undef writeparam

	return;
}

/*--------------------------------------------------------------------*/

void read_nuisance(char name[], struct nuisance* nuisparam)
/* Reads file "name" and assign the values into the the nuisparam structure */
{
	if(!test_file(name)) {printf("%s does not exist!\n",name); exit(1);}
	
	set_nuisance(nuisparam);
	set_nuisance_deviation_to_zero(nuisparam);

	char dummy[50];
	
	FILE *file;
	
	file=fopen(name,"r");

#define readparam( myparam ) if(!strcmp(dummy,#myparam)){nuisparam->myparam.name=#myparam; fscanf(file,"%lf",&nuisparam->myparam.cent); fscanf(file,"%lf",&nuisparam->myparam.dev); fscanf(file,"%d",&nuisparam->myparam.type);}

	while(EOF != fscanf(file,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(file,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		
		readparam(alphas_MZ)
		else readparam(mass_b)
		else readparam(mass_c)
		else readparam(mass_s)
		else readparam(mass_top_pole)
		else readparam(mass_h0)
		else readparam(CKM_lambda)
		else readparam(CKM_A)
		else readparam(CKM_rhobar)
		else readparam(CKM_etabar)
		else readparam(log_mu_W_mass_W)
		else readparam(log_mu_b_mass_b)
		else readparam(BR_BXclnu_exp)
		else readparam(mu_G2_bsg)
		else readparam(rho_D3_bsg)
		else readparam(rho_LS3_bsg)
		else readparam(bsgamma_rand)
		else readparam(mu_c_bsg)
		else readparam(BRBXsmumu_lowq2_rand)
		else readparam(BRBXsmumu_highq2_rand)
		else readparam(BRBXsmumu_full_rand)
		else readparam(BRBXsee_lowq2_rand)
		else readparam(BRBXsee_highq2_rand)
		else readparam(BRBXsee_full_rand)
		else readparam(BRBXstautau_lowq2_rand)
		else readparam(BRBXstautau_highq2_rand)
		else readparam(BRBXstautau_full_rand)
		else readparam(f_B)
		else readparam(lambda_Bp)
		else readparam(f_Kstar_par)
		else readparam(f_Kstar_perp)
		else readparam(a1perp)
		else readparam(a2perp)
		else readparam(a1par)
		else readparam(a2par)
		else readparam(T1_BKstar)
		else readparam(log_mu_spec_lambda_h_mass_b)
		else readparam(BtoKstarlow_ALperp_err_noq2)
		else readparam(BtoKstarlow_ARperp_err_noq2)
		else readparam(BtoKstarlow_ALpar_err_noq2)
		else readparam(BtoKstarlow_ARpar_err_noq2)
		else readparam(BtoKstarlow_AL0_err_noq2)
		else readparam(BtoKstarlow_AR0_err_noq2)
		else readparam(BtoKstarlow_At_err_noq2)
		else readparam(BtoKstarlow_AS_err_noq2)
		else readparam(BtoKstarlow_ALperp_err_q2)
		else readparam(BtoKstarlow_ARperp_err_q2)
		else readparam(BtoKstarlow_ALpar_err_q2)
		else readparam(BtoKstarlow_ARpar_err_q2)
		else readparam(BtoKstarlow_AL0_err_q2)
		else readparam(BtoKstarlow_AR0_err_q2)
		else readparam(BtoKstarlow_At_err_q2)
		else readparam(BtoKstarlow_AS_err_q2)
		/* SN: Result from B->K* Hadronic Fit (to be used to check how it impacts predictions) - 13/04/2025 */
		else readparam(hplus0r)
		else readparam(hminus0r)
		else readparam(hplus1r)
		else readparam(hminus1r)
		else readparam(hplus2r)
		else readparam(hminus2r)
		else readparam(hzero0r)
		else readparam(hzero1r)
		else readparam(hzero2r)
		else readparam(hplus0i)
		else readparam(hminus0i)
		else readparam(hplus1i)
		else readparam(hminus1i)
		else readparam(hplus2i)
		else readparam(hminus2i)
		else readparam(hzero0i)
		else readparam(hzero1i)
		else readparam(hzero2i)
		/**********************/
		else readparam(real_alpha_perp0)
		else readparam(real_alpha_perp1)
		else readparam(real_alpha_perp2)
		else readparam(real_alpha_par0)
		else readparam(real_alpha_par1)
		else readparam(real_alpha_par2)
		else readparam(real_alpha_zero0)
		else readparam(real_alpha_zero1)
		else readparam(imag_alpha_perp0)
		else readparam(imag_alpha_perp1)
		else readparam(imag_alpha_perp2)
		else readparam(imag_alpha_par0)
		else readparam(imag_alpha_par1)
		else readparam(imag_alpha_par2)
		else readparam(imag_alpha_zero0)
		else readparam(imag_alpha_zero1)
		else readparam(DeltaC9_M1_q2bar)
		else readparam(r1_M1)
		else readparam(r2_M1)
		else readparam(DeltaC9_M2_q2bar)
		else readparam(r1_M2)
		else readparam(r2_M2)
		else readparam(DeltaC9_M3_q2bar)
		else readparam(r1_M3)
		else readparam(r2_M3)
		else readparam(BtoKstarhigh_ALperp_err)
		else readparam(BtoKstarhigh_ARperp_err)
		else readparam(BtoKstarhigh_ALpar_err)
		else readparam(BtoKstarhigh_ARpar_err)
		else readparam(BtoKstarhigh_AL0_err)
		else readparam(BtoKstarhigh_AR0_err)
		else readparam(BtoKstarhigh_At_err)
		else readparam(BtoKstarhigh_AS_err)
		else readparam(a0A0_BKstar)
		else readparam(a1A0_BKstar)
		else readparam(a2A0_BKstar)
		else readparam(a0A1_BKstar)
		else readparam(a1A1_BKstar)
		else readparam(a2A1_BKstar)
		else readparam(a0A12_BKstar)
		else readparam(a1A12_BKstar)
		else readparam(a2A12_BKstar)
		else readparam(a0V_BKstar)
		else readparam(a1V_BKstar)
		else readparam(a2V_BKstar)
		else readparam(a0T1_BKstar)
		else readparam(a1T1_BKstar)
		else readparam(a2T1_BKstar)
		else readparam(a0T2_BKstar)
		else readparam(a1T2_BKstar)
		else readparam(a2T2_BKstar)
		else readparam(a0T23_BKstar)
		else readparam(a1T23_BKstar)
		else readparam(a2T23_BKstar)
		/* SN: B -> Kstar FF van Dyk et al. 2305.06301 (BKstar_FormFactor_choice:2)*/
		else readparam(a0A0_BKstar_GRvDV_BSZ)
		else readparam(a1A0_BKstar_GRvDV_BSZ)
		else readparam(a2A0_BKstar_GRvDV_BSZ)
		else readparam(a0A1_BKstar_GRvDV_BSZ)
		else readparam(a1A1_BKstar_GRvDV_BSZ)
		else readparam(a2A1_BKstar_GRvDV_BSZ)
		else readparam(a1A12_BKstar_GRvDV_BSZ)
		else readparam(a2A12_BKstar_GRvDV_BSZ)
		else readparam(a0V_BKstar_GRvDV_BSZ)
		else readparam(a1V_BKstar_GRvDV_BSZ)
		else readparam(a2V_BKstar_GRvDV_BSZ)
		else readparam(a0T1_BKstar_GRvDV_BSZ)
		else readparam(a1T1_BKstar_GRvDV_BSZ)
		else readparam(a2T1_BKstar_GRvDV_BSZ)
		else readparam(a1T2_BKstar_GRvDV_BSZ)
		else readparam(a2T2_BKstar_GRvDV_BSZ)
		else readparam(a0T23_BKstar_GRvDV_BSZ)
		else readparam(a1T23_BKstar_GRvDV_BSZ)
		else readparam(a2T23_BKstar_GRvDV_BSZ)
		/* SN: LCSR+Lattice B -> Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:3)*/
		else readparam(a0A0_BKstar_GKvD_LCSR_Lattice)
		else readparam(a1A0_BKstar_GKvD_LCSR_Lattice)
		else readparam(a2A0_BKstar_GKvD_LCSR_Lattice)
		else readparam(a0A1_BKstar_GKvD_LCSR_Lattice)
		else readparam(a1A1_BKstar_GKvD_LCSR_Lattice)
		else readparam(a2A1_BKstar_GKvD_LCSR_Lattice)
		else readparam(a1A12_BKstar_GKvD_LCSR_Lattice)
		else readparam(a2A12_BKstar_GKvD_LCSR_Lattice)
		else readparam(a0V_BKstar_GKvD_LCSR_Lattice)
		else readparam(a1V_BKstar_GKvD_LCSR_Lattice)
		else readparam(a2V_BKstar_GKvD_LCSR_Lattice)
		else readparam(a0T1_BKstar_GKvD_LCSR_Lattice)
		else readparam(a1T1_BKstar_GKvD_LCSR_Lattice)
		else readparam(a2T1_BKstar_GKvD_LCSR_Lattice)
		// else readparam(a0T2_BKstar_GKvD_LCSR_Lattice)
		else readparam(a1T2_BKstar_GKvD_LCSR_Lattice)
		else readparam(a2T2_BKstar_GKvD_LCSR_Lattice)
		else readparam(a0T23_BKstar_GKvD_LCSR_Lattice)
		else readparam(a1T23_BKstar_GKvD_LCSR_Lattice)
		else readparam(a2T23_BKstar_GKvD_LCSR_Lattice)
		/* SN: LCSR-only B -> Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:4)*/
		else readparam(a0A0_BKstar_GKvD_LCSRonly)
		else readparam(a1A0_BKstar_GKvD_LCSRonly)
		else readparam(a2A0_BKstar_GKvD_LCSRonly)
		else readparam(a0A1_BKstar_GKvD_LCSRonly)
		else readparam(a1A1_BKstar_GKvD_LCSRonly)
		else readparam(a2A1_BKstar_GKvD_LCSRonly)
		else readparam(a1A12_BKstar_GKvD_LCSRonly)
		else readparam(a2A12_BKstar_GKvD_LCSRonly)
		else readparam(a0V_BKstar_GKvD_LCSRonly)
		else readparam(a1V_BKstar_GKvD_LCSRonly)
		else readparam(a2V_BKstar_GKvD_LCSRonly)
		else readparam(a0T1_BKstar_GKvD_LCSRonly)
		else readparam(a1T1_BKstar_GKvD_LCSRonly)
		else readparam(a2T1_BKstar_GKvD_LCSRonly)
		// else readparam(a0T2_BKstar_GKvD_LCSRonly)
		else readparam(a1T2_BKstar_GKvD_LCSRonly)
		else readparam(a2T2_BKstar_GKvD_LCSRonly)
		else readparam(a0T23_BKstar_GKvD_LCSRonly)
		else readparam(a1T23_BKstar_GKvD_LCSRonly)
		else readparam(a2T23_BKstar_GKvD_LCSRonly)
	/* SN: LCSR-only B->Kstar FF Zwicky et al. 1503.05534 (BKstar_FormFactor_choice:5)*/
		else readparam(a0A0_BKstar_BSZ_LCSRonly)
		else readparam(a1A0_BKstar_BSZ_LCSRonly)
		else readparam(a2A0_BKstar_BSZ_LCSRonly)
		else readparam(a0A1_BKstar_BSZ_LCSRonly)
		else readparam(a1A1_BKstar_BSZ_LCSRonly)
		else readparam(a2A1_BKstar_BSZ_LCSRonly)
		else readparam(a0A12_BKstar_BSZ_LCSRonly)
		else readparam(a1A12_BKstar_BSZ_LCSRonly)
		else readparam(a2A12_BKstar_BSZ_LCSRonly)
		else readparam(a0V_BKstar_BSZ_LCSRonly)
		else readparam(a1V_BKstar_BSZ_LCSRonly)
		else readparam(a2V_BKstar_BSZ_LCSRonly)
		else readparam(a0T1_BKstar_BSZ_LCSRonly)
		else readparam(a1T1_BKstar_BSZ_LCSRonly)
		else readparam(a2T1_BKstar_BSZ_LCSRonly)
		else readparam(a0T2_BKstar_BSZ_LCSRonly)
		else readparam(a1T2_BKstar_BSZ_LCSRonly)
		else readparam(a2T2_BKstar_BSZ_LCSRonly)
		else readparam(a0T23_BKstar_BSZ_LCSRonly)
		else readparam(a1T23_BKstar_BSZ_LCSRonly)
		else readparam(a2T23_BKstar_BSZ_LCSRonly)
	/* SN: Lattice-only B->Kstar FF Wingate et al. 1310.3722 + 1501.00367 (BKstar_FormFactor_choice:6*/
		else readparam(a0V_BKstar_HLMW)
		else readparam(a1V_BKstar_HLMW)
		else readparam(a0A0_BKstar_HLMW)
		else readparam(a1A0_BKstar_HLMW)
		else readparam(a0A1_BKstar_HLMW)
		else readparam(a1A1_BKstar_HLMW)
		else readparam(a0A12_BKstar_HLMW)
		else readparam(a1A12_BKstar_HLMW)
		else readparam(a0T1_BKstar_HLMW)
		else readparam(a1T1_BKstar_HLMW)
		else readparam(a0T2_BKstar_HLMW)
		else readparam(a1T2_BKstar_HLMW)
		else readparam(a0T23_BKstar_HLMW)
		else readparam(a1T23_BKstar_HLMW)
		else readparam(systErr_BKstar_HLMW)

		else readparam(f_K)
		else readparam(a1K)
		else readparam(a2K)
		/* SN: Replacing with more explicit naming and adding missing a3f0 - 28/03/2025 */
		// else readparam(a00_BK)
		// else readparam(a10_BK)
		// else readparam(a20_BK)
		// else readparam(a0p_BK)
		// else readparam(a1p_BK)
		// else readparam(a2p_BK)
		// else readparam(a0T_BK)
		// else readparam(a1T_BK)
		// else readparam(a2T_BK)
		else readparam(a0f0_BK_AS_LCSR_Lattice)
		else readparam(a1f0_BK_AS_LCSR_Lattice)
		else readparam(a2f0_BK_AS_LCSR_Lattice)
		else readparam(a3f0_BK_AS_LCSR_Lattice)
		else readparam(a0fp_BK_AS_LCSR_Lattice)
		else readparam(a1fp_BK_AS_LCSR_Lattice)
		else readparam(a2fp_BK_AS_LCSR_Lattice)
		else readparam(a0fT_BK_AS_LCSR_Lattice)
		else readparam(a1fT_BK_AS_LCSR_Lattice)
		else readparam(a2fT_BK_AS_LCSR_Lattice)
		/* SN: B -> K FF van Dyk et al. 2305.06301 (BK_FormFactor_choice:2)*/
		else readparam(a0fp_BK_GRvDV_BSZ)
		else readparam(a1fp_BK_GRvDV_BSZ)
		else readparam(a2fp_BK_GRvDV_BSZ)
		else readparam(a1f0_BK_GRvDV_BSZ)
		else readparam(a2f0_BK_GRvDV_BSZ)
		else readparam(a0fT_BK_GRvDV_BSZ)
		else readparam(a1fT_BK_GRvDV_BSZ)
		else readparam(a2fT_BK_GRvDV_BSZ)
		/* SN: LCSR+Lattice B -> K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:3)*/
		else readparam(a0fp_BK_GKvD_LCSR_Lattice)
		else readparam(a1fp_BK_GKvD_LCSR_Lattice)
		else readparam(a2fp_BK_GKvD_LCSR_Lattice)
		else readparam(a1f0_BK_GKvD_LCSR_Lattice)
		else readparam(a2f0_BK_GKvD_LCSR_Lattice)
		else readparam(a0fT_BK_GKvD_LCSR_Lattice)
		else readparam(a1fT_BK_GKvD_LCSR_Lattice)
		else readparam(a2fT_BK_GKvD_LCSR_Lattice)
		/* SN: LCSR-only B -> K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:4)*/
		else readparam(a0fp_BK_GKvD_LCSRonly)
		else readparam(a1fp_BK_GKvD_LCSRonly)
		else readparam(a2fp_BK_GKvD_LCSRonly)
		else readparam(a1f0_BK_GKvD_LCSRonly)
		else readparam(a2f0_BK_GKvD_LCSRonly)
		else readparam(a0fT_BK_GKvD_LCSRonly)
		else readparam(a1fT_BK_GKvD_LCSRonly)
		else readparam(a2fT_BK_GKvD_LCSRonly)
		/* SN: Lattice FLAG24 average 2411.04268 of HPQCD'13 and FNAL/MILC'15 (BK_FormFactor_choice:6) - 26/04/2025*/
		else readparam(a0fp_BK_FLAG24)
		else readparam(a1fp_BK_FLAG24)
		else readparam(a2fp_BK_FLAG24)
		else readparam(a0f0_BK_FLAG24)
		else readparam(a1f0_BK_FLAG24)
		else readparam(a0fT_BK_FLAG24)
		else readparam(a1fT_BK_FLAG24)
		else readparam(a2fT_BK_FLAG24)
		/* SN: Lattice full q2-region HPQCD22 2207.12468 (BK_FormFactor_choice:7) - 26/04/2025*/
		else readparam(a0f0_BK_HPQCD22)
		else readparam(a1f0_BK_HPQCD22)
		else readparam(a2f0_BK_HPQCD22)
		else readparam(a1fp_BK_HPQCD22)
		else readparam(a2fp_BK_HPQCD22)
		else readparam(a0fT_BK_HPQCD22)
		else readparam(a1fT_BK_HPQCD22)
		else readparam(a2fT_BK_HPQCD22)
		else readparam(L_BK_HPQCD22)
		/* SN: LCSR-only Khodjamirian & Rusov 1703.04765 (BK_FormFactor_choice:8) - 02/08/2025*/
		else readparam(f0fp_BK_KR)
		else readparam(b1fp_BK_KR)
		else readparam(f0fT_BK_KR)
		else readparam(b1fT_BK_KR)
		/* ***************************** */
		else readparam(BtoKlow_FV_err_noq2)
		else readparam(BtoKlow_FA_err_noq2)
		else readparam(BtoKlow_FS_err_noq2)
		else readparam(BtoKlow_FP_err_noq2)
		else readparam(BtoKlow_FV_err_q2)
		else readparam(BtoKlow_FA_err_q2)
		else readparam(BtoKlow_FS_err_q2)
		else readparam(BtoKlow_FP_err_q2)
		else readparam(BtoKhigh_FV_err)
		else readparam(BtoKhigh_FA_err)
		else readparam(BtoKhigh_FS_err)
		else readparam(BtoKhigh_FP_err)
		else readparam(life_Bs)
		else readparam(f_Bs)
		else readparam(lambda_Bsp)
		else readparam(ys_Bs)
		else readparam(f_phi_par)
		else readparam(f_phi_perp)
		else readparam(a1phi_perp)
		else readparam(a1phi_par)
		else readparam(a2phi_perp)
		else readparam(a2phi_par)
		else readparam(Bstophilow_ALperp_err_noq2)
		else readparam(Bstophilow_ARperp_err_noq2)
		else readparam(Bstophilow_ALpar_err_noq2)
		else readparam(Bstophilow_ARpar_err_noq2)
		else readparam(Bstophilow_AL0_err_noq2)
		else readparam(Bstophilow_AR0_err_noq2)
		else readparam(Bstophilow_At_err_noq2)
		else readparam(Bstophilow_AS_err_noq2)
		else readparam(Bstophilow_ALperp_err_q2)
		else readparam(Bstophilow_ARperp_err_q2)
		else readparam(Bstophilow_ALpar_err_q2)
		else readparam(Bstophilow_ARpar_err_q2)
		else readparam(Bstophilow_AL0_err_q2)
		else readparam(Bstophilow_AR0_err_q2)
		else readparam(Bstophilow_At_err_q2)
		else readparam(Bstophilow_AS_err_q2)
		else readparam(Bstophihigh_ALperp_err)
		else readparam(Bstophihigh_ARperp_err)
		else readparam(Bstophihigh_ALpar_err)
		else readparam(Bstophihigh_ARpar_err)
		else readparam(Bstophihigh_AL0_err)
		else readparam(Bstophihigh_AR0_err)
		else readparam(Bstophihigh_At_err)
		else readparam(Bstophihigh_AS_err)
		else readparam(a0A0_Bsphi)
		else readparam(a1A0_Bsphi)
		else readparam(a2A0_Bsphi)
		else readparam(a0A1_Bsphi)
		else readparam(a1A1_Bsphi)
		else readparam(a2A1_Bsphi)
		else readparam(a0A12_Bsphi)
		else readparam(a1A12_Bsphi)
		else readparam(a2A12_Bsphi)
		else readparam(a0V_Bsphi)
		else readparam(a1V_Bsphi)
		else readparam(a2V_Bsphi)
		else readparam(a0T1_Bsphi)
		else readparam(a1T1_Bsphi)
		else readparam(a2T1_Bsphi)
		else readparam(a0T2_Bsphi)
		else readparam(a1T2_Bsphi)
		else readparam(a2T2_Bsphi)
		else readparam(a0T23_Bsphi)
		else readparam(a1T23_Bsphi)
		else readparam(a2T23_Bsphi)
        /* SN: Bs -> phi FF van Dyk et al. 2305.06301 (Bsphi_FormFactor_choice:2)*/
		else readparam(a0A0_Bsphi_GRvDV_BSZ)
		else readparam(a1A0_Bsphi_GRvDV_BSZ)
		else readparam(a2A0_Bsphi_GRvDV_BSZ)
		else readparam(a0A1_Bsphi_GRvDV_BSZ)
		else readparam(a1A1_Bsphi_GRvDV_BSZ)
		else readparam(a2A1_Bsphi_GRvDV_BSZ)
		else readparam(a1A12_Bsphi_GRvDV_BSZ)
		else readparam(a2A12_Bsphi_GRvDV_BSZ)
		else readparam(a0V_Bsphi_GRvDV_BSZ)
		else readparam(a1V_Bsphi_GRvDV_BSZ)
		else readparam(a2V_Bsphi_GRvDV_BSZ)
		else readparam(a0T1_Bsphi_GRvDV_BSZ)
		else readparam(a1T1_Bsphi_GRvDV_BSZ)
		else readparam(a2T1_Bsphi_GRvDV_BSZ)
		else readparam(a1T2_Bsphi_GRvDV_BSZ)
		else readparam(a2T2_Bsphi_GRvDV_BSZ)
		else readparam(a0T23_Bsphi_GRvDV_BSZ)
		else readparam(a1T23_Bsphi_GRvDV_BSZ)
		else readparam(a2T23_Bsphi_GRvDV_BSZ)
		/* SN: Bs->phi FF Zwicky et al. 1503.05534 (Bsphi_FormFactor_choice:5)*/
		else readparam(a0A0_Bsphi_BSZ_LCSRonly)
		else readparam(a1A0_Bsphi_BSZ_LCSRonly)
		else readparam(a2A0_Bsphi_BSZ_LCSRonly)
		else readparam(a0A1_Bsphi_BSZ_LCSRonly)
		else readparam(a1A1_Bsphi_BSZ_LCSRonly)
		else readparam(a2A1_Bsphi_BSZ_LCSRonly)
		else readparam(a0A12_Bsphi_BSZ_LCSRonly)
		else readparam(a1A12_Bsphi_BSZ_LCSRonly)
		else readparam(a2A12_Bsphi_BSZ_LCSRonly)
		else readparam(a0V_Bsphi_BSZ_LCSRonly)
		else readparam(a1V_Bsphi_BSZ_LCSRonly)
		else readparam(a2V_Bsphi_BSZ_LCSRonly)
		else readparam(a0T1_Bsphi_BSZ_LCSRonly)
		else readparam(a1T1_Bsphi_BSZ_LCSRonly)
		else readparam(a2T1_Bsphi_BSZ_LCSRonly)
		else readparam(a0T2_Bsphi_BSZ_LCSRonly)
		else readparam(a1T2_Bsphi_BSZ_LCSRonly)
		else readparam(a2T2_Bsphi_BSZ_LCSRonly)
		else readparam(a0T23_Bsphi_BSZ_LCSRonly)
		else readparam(a1T23_Bsphi_BSZ_LCSRonly)
		else readparam(a2T23_Bsphi_BSZ_LCSRonly)
		/* SN: Bs->phi FF Wingate et al. 1310.3722 + 1501.00367 (Bsphi_FormFactor_choice:6)*/
		else readparam(a0V_Bsphi_HLMW)
		else readparam(a1V_Bsphi_HLMW)
		else readparam(a0A0_Bsphi_HLMW)
		else readparam(a1A0_Bsphi_HLMW)
		else readparam(a0A1_Bsphi_HLMW)
		else readparam(a1A1_Bsphi_HLMW)
		else readparam(a0A12_Bsphi_HLMW)
		else readparam(a1A12_Bsphi_HLMW)
		else readparam(a0T1_Bsphi_HLMW)
		else readparam(a1T1_Bsphi_HLMW)
		else readparam(a0T2_Bsphi_HLMW)
		else readparam(a1T2_Bsphi_HLMW)
		else readparam(a0T23_Bsphi_HLMW)
		else readparam(a1T23_Bsphi_HLMW)
		else readparam(systErr_Bsphi_HLMW)
		else readparam(life_Lb)
		else readparam(alphaL_LbLll)
		else readparam(a0_HO_fplus_LbLll)
		else readparam(a1_HO_fplus_LbLll)
		else readparam(a2_HO_fplus_LbLll)
		else readparam(a0_HO_fperp_LbLll)
		else readparam(a1_HO_fperp_LbLll)
		else readparam(a2_HO_fperp_LbLll)
		else readparam(a0_HO_gpp_LbLll)
		else readparam(a1_HO_gplus_LbLll)
		else readparam(a2_HO_gplus_LbLll)
		else readparam(a1_HO_gperp_LbLll)
		else readparam(a2_HO_gperp_LbLll)
		else readparam(a0_HO_hplus_LbLll)
		else readparam(a1_HO_hplus_LbLll)
		else readparam(a2_HO_hplus_LbLll)
		else readparam(a0_HO_hperp_LbLll)
		else readparam(a1_HO_hperp_LbLll)
		else readparam(a2_HO_hperp_LbLll)
		else readparam(a0_HO_htildepp_LbLll)
		else readparam(a1_HO_htildeplus_LbLll)
		else readparam(a2_HO_htildeplus_LbLll)
		else readparam(a1_HO_htildeperp_LbLll)
		else readparam(a2_HO_htildeperp_LbLll)
		else readparam(deltaPcu_Kppipnunu)
		else readparam(err_Pc_Xlambda_Kppipnunu)
		else readparam(BR_KLgammagamma_exp)
		else readparam(Aterm_mu_KLmumu)
		else readparam(chi_gg_Mrho)
		else readparam(BR_KSgammagamma_exp)
		else readparam(Iterm_mu_KSmumu)
		else readparam(log_muK_1GeV)
		else readparam(KLpill_Ce_Cdir)
		else readparam(KLpill_Cmu_Cdir)
		else readparam(KLpill_Ce_Cint)
		else readparam(KLpill_Cmu_Cint)
		else readparam(KLpill_Ce_Cmix)
		else readparam(KLpill_Cmu_Cmix)
		else readparam(KLpill_Cmu_CPC)
		else readparam(KLpill_abs_aS)
		/* Scalar and Pseudoscalar contributions included - 08/10/2025 */
		else readparam(KLpill_Ce_Sgg)
		else readparam(KLpill_Cmu_Sgg)
		/* AFBmu_KLpi0mumu  - 09/10/2025 */
		else readparam(KLpill_AFB_a1pi_gg_7V)
		else readparam(KLpill_AFB_a1pi_gg_aS)
	}
	
	fclose(file);

#undef readparam

	return;
}

/*---------------------------------------------------------------------*/

int read_nameobs(char name[], obsname* obs)
/* Reads a string "name" and converts it into observable name "obs" */
{	
	int length=strlen(name);    
	
	if(length==0) return 0;
	
	sprintf(obs->name,"%s",name);

    char *start, *end;        
    start=name;

    if ( (end = strstr(start,"_")) )
    {
		char *temp;
		temp = (char *)malloc(end-start+1);
		memcpy(temp, start, end - start );
		temp[end-start] = '\0';
		sprintf(obs->type,"%s",temp);
		free(temp);

        start=end+1;
    }
    else 
    {
		sprintf(obs->type,"%s",start);
		sprintf(obs->decay,"");
		obs->low=obs->high=0.;
		sprintf(obs->other,"");
		
		return 1;
    }

    if ( (end = strstr(start,"_")) )
    {

		char *temp;
		temp = (char *)malloc(end-start+1);
		memcpy(temp, start, end - start );
		temp[end - start] = '\0';
		sprintf(obs->decay,"%s",temp);
		free(temp);

        start=end+1;
    }
    else 
    {
		sprintf(obs->decay,"%s",start);
		obs->low=obs->high=0.;
		sprintf(obs->other,"");
	
		return 2;
    }

    if ( (end = strstr(start,"_")) )
    {
		char *temp;
		temp = (char *)malloc(end-start+1);
		memcpy(temp, start, end - start );
		temp[end - start] = '\0';
		if(atof(temp)!=0.)
		{
			obs->low=atof(temp);
		}
		else
		{
			obs->low=obs->high=0.;
			sprintf(obs->other,"%s",temp);
		}
		free(temp);

        start=end+1;
    }
    else 
    {
		if(atof(start)!=0.)
		{
			obs->low=atof(start);		
			obs->high=0.;
			sprintf(obs->other,"");
		}
		else
		{
			obs->low=obs->high=0.;
			sprintf(obs->other,"%s",start);
		}
		
		return 3;
    }

    if ( (end = strstr(start,"_")) )
    {
		char *temp;
		temp = (char *)malloc(end-start+1);
		memcpy(temp, start, end - start );
		temp[end - start] = '\0';
		if(atof(temp)!=0.)
		{
			obs->high=atof(temp);
		}
		else
		{
			sprintf(obs->other,"%s",temp);
		}
		free(temp);

        start=end+1;
    }
    else 
    {
		if(atof(start)!=0.)
		{
			obs->high=atof(start);		
			sprintf(obs->other,"");
		}
		else
		{
			obs->high=0.;
			sprintf(obs->other,"%s",start);
		}
		
		return 4;
    }

	sprintf(obs->other,"%s",start);

	return 5;
}

/*--------------------------------------------------------------------*/

void observables(int ke, obsname obs[], int nobs, double values[], double values_ref[], struct nuisance* nuisparam, char namenuisance[][50], struct parameters* param0)
/* Computes the "values" of the "nobs" observables "obs" for the ke-th iteration of the theory error calculation, using the nuisance parameters "nuisparam" (named as "namenuisance"). "values_ref" contains the values of the observables obtained at the 0th (reference) iteration */
{		
	int nbobs=0;
	
	int ie;
	
	int method; /* Remark: myrand(X,-1)=1, myrand(X,0)=0, myrand(X,1)=random number */
	if(ke==0) method=0;
	else if(ke>0) method=-1;
	else if(ke<0) method=1;

	if(ke>0) {for(ie=0;ie<nobs;ie++) values[ie]=values_ref[ie];}
	else {for(ie=0;ie<nobs;ie++) values[ie]=0.;}

	struct parameters param=*param0;
	
#define setparam( myparam , number ) if(ke==0){strcpy(namenuisance[number],nuisparam->myparam.name);} if(ke==number&&nuisparam->myparam.dev==0.){for(ie=0;ie<nobs;ie++) values[ie]=values_ref[ie]; return;} param.myparam=(nuisparam->myparam.cent+nuisparam->myparam.dev*myrand(nuisparam->myparam.type,method)*(method==1||ke==number))

#define setstdzero( myparam ) nuisparam->myparam.dev=0.;

/* nuisance for all observables */
	setparam(alphas_MZ,1);
	
	setparam(mass_b,2);
	setparam(mass_c,3);
	setparam(mass_s,4);
	setparam(mass_top_pole,5);

	setparam(mass_h0,6);

	setparam(CKM_lambda,7);
	setparam(CKM_A,8);
	setparam(CKM_rhobar,9);
	setparam(CKM_etabar,10);

	slha_adjust(&param);

	setparam(log_mu_W_mass_W,11);
	setparam(log_mu_b_mass_b,12);

/* inclusive b -> s */
	setparam(BR_BXclnu_exp,13);
/* b -> s gamma */ 
	setparam(mu_G2_bsg,14);
	setparam(rho_D3_bsg,15);
	setparam(rho_LS3_bsg,16);
	setparam(bsgamma_rand,17);
	setparam(mu_c_bsg,18);
/* b -> s mu mu */
	setparam(BRBXsmumu_lowq2_rand,19);
	setparam(BRBXsmumu_highq2_rand,20);
	setparam(BRBXsmumu_full_rand,21);
/* b -> s e e */
	setparam(BRBXsee_lowq2_rand,22);
	setparam(BRBXsee_highq2_rand,23);
	setparam(BRBXsee_full_rand,24);
/* b -> s tau tau */
	setparam(BRBXstautau_lowq2_rand,25);
	setparam(BRBXstautau_highq2_rand,26);
	setparam(BRBXstautau_full_rand,27);

/* B */	
	setparam(f_B,28);
	setparam(lambda_Bp,29);
/* B -> K* */
	setparam(f_Kstar_par,30);
	setparam(f_Kstar_perp,31);
	setparam(a1perp,32);
	setparam(a2perp,33);
	setparam(a1par,34);
	setparam(a2par,35);
/* B -> K* gamma */
	setparam(T1_BKstar,36);
	setparam(log_mu_spec_lambda_h_mass_b,37);
/* low */
	switch(param.BKstar_implementation)
	{
		default:
		{
			setparam(BtoKstarlow_ALperp_err_noq2,38)*cexp(I*pi/4.);
			setparam(BtoKstarlow_ARperp_err_noq2,39)*cexp(I*pi/4.);
			setparam(BtoKstarlow_ALpar_err_noq2,40)*cexp(I*pi/4.);
			setparam(BtoKstarlow_ARpar_err_noq2,41)*cexp(I*pi/4.);
			setparam(BtoKstarlow_AL0_err_noq2,42)*cexp(I*pi/4.);
			setparam(BtoKstarlow_AR0_err_noq2,43)*cexp(I*pi/4.);
			setparam(BtoKstarlow_At_err_noq2,44)*cexp(I*pi/4.);
			setparam(BtoKstarlow_AS_err_noq2,45)*cexp(I*pi/4.);

			setparam(BtoKstarlow_ALperp_err_q2,46)*cexp(I*pi/4.);
			setparam(BtoKstarlow_ARperp_err_q2,47)*cexp(I*pi/4.);
			setparam(BtoKstarlow_ALpar_err_q2,48)*cexp(I*pi/4.);
			setparam(BtoKstarlow_ARpar_err_q2,49)*cexp(I*pi/4.);
			setparam(BtoKstarlow_AL0_err_q2,50)*cexp(I*pi/4.);
			setparam(BtoKstarlow_AR0_err_q2,51)*cexp(I*pi/4.);
			setparam(BtoKstarlow_At_err_q2,52)*cexp(I*pi/4.);
			setparam(BtoKstarlow_AS_err_q2,53)*cexp(I*pi/4.);

			setstdzero(real_alpha_perp0);
			setstdzero(real_alpha_perp1);
			setstdzero(real_alpha_perp2);
			setstdzero(real_alpha_par0);
			setstdzero(real_alpha_par1);
			setstdzero(real_alpha_par2);
			setstdzero(real_alpha_zero0);
			setstdzero(real_alpha_zero1);
			setstdzero(imag_alpha_perp0);
			setstdzero(imag_alpha_perp1);
			setstdzero(imag_alpha_perp2);
			setstdzero(imag_alpha_par0);
			setstdzero(imag_alpha_par1);
			setstdzero(imag_alpha_par2);
			setstdzero(imag_alpha_zero0);
			setstdzero(imag_alpha_zero1);

			setstdzero(DeltaC9_M1_q2bar);
			setstdzero(r1_M1);
			setstdzero(r2_M1);
			setstdzero(DeltaC9_M2_q2bar);
			setstdzero(r1_M2);
			setstdzero(r2_M2);
			setstdzero(DeltaC9_M3_q2bar);
			setstdzero(r1_M3);
			setstdzero(r2_M3);
				
			break;
		}
		
		case 2:
		{
			setparam(real_alpha_perp0,38);
			setparam(real_alpha_perp1,39);
			setparam(real_alpha_perp2,40);
			setparam(real_alpha_par0,41);
			setparam(real_alpha_par1,42);
			setparam(real_alpha_par2,43);
			setparam(real_alpha_zero0,44);
			setparam(real_alpha_zero1,45);
			setparam(imag_alpha_perp0,46);
			setparam(imag_alpha_perp1,47);
			setparam(imag_alpha_perp2,48);
			setparam(imag_alpha_par0,49);
			setparam(imag_alpha_par1,50);
			setparam(imag_alpha_par2,51);
			setparam(imag_alpha_zero0,52);
			setparam(imag_alpha_zero1,53);

			setstdzero(BtoKstarlow_ALperp_err_noq2);
			setstdzero(BtoKstarlow_ARperp_err_noq2);
			setstdzero(BtoKstarlow_ALpar_err_noq2);
			setstdzero(BtoKstarlow_ARpar_err_noq2);
			setstdzero(BtoKstarlow_AL0_err_noq2);
			setstdzero(BtoKstarlow_AR0_err_noq2);
			setstdzero(BtoKstarlow_At_err_noq2);
			setstdzero(BtoKstarlow_AS_err_noq2);

			setstdzero(BtoKstarlow_ALperp_err_q2);
			setstdzero(BtoKstarlow_ARperp_err_q2);
			setstdzero(BtoKstarlow_ALpar_err_q2);
			setstdzero(BtoKstarlow_ARpar_err_q2);
			setstdzero(BtoKstarlow_AL0_err_q2);
			setstdzero(BtoKstarlow_AR0_err_q2);
			setstdzero(BtoKstarlow_At_err_q2);
			setstdzero(BtoKstarlow_AS_err_q2);

			setstdzero(DeltaC9_M1_q2bar);
			setstdzero(r1_M1);
			setstdzero(r2_M1);
			setstdzero(DeltaC9_M2_q2bar);
			setstdzero(r1_M2);
			setstdzero(r2_M2);
			setstdzero(DeltaC9_M3_q2bar);
			setstdzero(r1_M3);
			setstdzero(r2_M3);
			
			break;
		}
		
		case 3:
		{
			setparam(DeltaC9_M1_q2bar,38);
			setparam(r1_M1,39);
			setparam(r2_M1,40);
			setparam(DeltaC9_M2_q2bar,41);
			setparam(r1_M2,42);
			setparam(r2_M2,43);
			setparam(DeltaC9_M3_q2bar,44);
			setparam(r1_M3,45);
			setparam(r2_M3,46);

			setstdzero(BtoKstarlow_ALperp_err_noq2);
			setstdzero(BtoKstarlow_ARperp_err_noq2);
			setstdzero(BtoKstarlow_ALpar_err_noq2);
			setstdzero(BtoKstarlow_ARpar_err_noq2);
			setstdzero(BtoKstarlow_AL0_err_noq2);
			setstdzero(BtoKstarlow_AR0_err_noq2);
			setstdzero(BtoKstarlow_At_err_noq2);
			setstdzero(BtoKstarlow_AS_err_noq2);

			setstdzero(BtoKstarlow_ALperp_err_q2);
			setstdzero(BtoKstarlow_ARperp_err_q2);
			setstdzero(BtoKstarlow_ALpar_err_q2);
			setstdzero(BtoKstarlow_ARpar_err_q2);
			setstdzero(BtoKstarlow_AL0_err_q2);
			setstdzero(BtoKstarlow_AR0_err_q2);
			setstdzero(BtoKstarlow_At_err_q2);
			setstdzero(BtoKstarlow_AS_err_q2);
			setstdzero(real_alpha_perp0);
			setstdzero(real_alpha_perp1);
			setstdzero(real_alpha_perp2);
			setstdzero(real_alpha_par0);
			setstdzero(real_alpha_par1);
			setstdzero(real_alpha_par2);
			setstdzero(real_alpha_zero0);
			setstdzero(real_alpha_zero1);
			setstdzero(imag_alpha_perp0);
			setstdzero(imag_alpha_perp1);
			setstdzero(imag_alpha_perp2);
			setstdzero(imag_alpha_par0);
			setstdzero(imag_alpha_par1);
			setstdzero(imag_alpha_par2);
			setstdzero(imag_alpha_zero0);
			setstdzero(imag_alpha_zero1);
			
			break;
		}
		/* SN: Result from B->K* Hadronic Fit (to be used to check how it impacts predictions) - 13/04/2025 */
		case 4:
		{
			setparam(hplus0r,38);
			setparam(hminus0r,39);
			setparam(hplus1r,40);
			setparam(hminus1r,41);
			setparam(hplus2r,42);
			setparam(hminus2r,43);
			setparam(hzero0r,44);
			setparam(hzero1r,45);
			setparam(hzero2r,46);
			setparam(hplus0i,47);
			setparam(hminus0i,48);
			setparam(hplus1i,49);
			setparam(hminus1i,50);
			setparam(hplus2i,51);
			setparam(hminus2i,52);
			setparam(hzero0i,53);
			setparam(hzero1i,194);
			setparam(hzero2i,195);

			setstdzero(DeltaC9_M1_q2bar);
			setstdzero(r1_M1);
			setstdzero(r2_M1);
			setstdzero(DeltaC9_M2_q2bar);
			setstdzero(r1_M2);
			setstdzero(r2_M2);
			setstdzero(DeltaC9_M3_q2bar);
			setstdzero(r1_M3);
			setstdzero(r2_M3);

			setstdzero(BtoKstarlow_ALperp_err_noq2);
			setstdzero(BtoKstarlow_ARperp_err_noq2);
			setstdzero(BtoKstarlow_ALpar_err_noq2);
			setstdzero(BtoKstarlow_ARpar_err_noq2);
			setstdzero(BtoKstarlow_AL0_err_noq2);
			setstdzero(BtoKstarlow_AR0_err_noq2);
			setstdzero(BtoKstarlow_At_err_noq2);
			setstdzero(BtoKstarlow_AS_err_noq2);

			setstdzero(BtoKstarlow_ALperp_err_q2);
			setstdzero(BtoKstarlow_ARperp_err_q2);
			setstdzero(BtoKstarlow_ALpar_err_q2);
			setstdzero(BtoKstarlow_ARpar_err_q2);
			setstdzero(BtoKstarlow_AL0_err_q2);
			setstdzero(BtoKstarlow_AR0_err_q2);
			setstdzero(BtoKstarlow_At_err_q2);
			setstdzero(BtoKstarlow_AS_err_q2);
			setstdzero(real_alpha_perp0);
			setstdzero(real_alpha_perp1);
			setstdzero(real_alpha_perp2);
			setstdzero(real_alpha_par0);
			setstdzero(real_alpha_par1);
			setstdzero(real_alpha_par2);
			setstdzero(real_alpha_zero0);
			setstdzero(real_alpha_zero1);
			setstdzero(imag_alpha_perp0);
			setstdzero(imag_alpha_perp1);
			setstdzero(imag_alpha_perp2);
			setstdzero(imag_alpha_par0);
			setstdzero(imag_alpha_par1);
			setstdzero(imag_alpha_par2);
			setstdzero(imag_alpha_zero0);
			setstdzero(imag_alpha_zero1);

			break;
		}
		/************************************/
	}
	
/* high */
	setparam(BtoKstarhigh_ALperp_err,54)*cexp(I*pi/4.);
	setparam(BtoKstarhigh_ARperp_err,55)*cexp(I*pi/4.);
	setparam(BtoKstarhigh_ALpar_err,56)*cexp(I*pi/4.);
	setparam(BtoKstarhigh_ARpar_err,57)*cexp(I*pi/4.);
	setparam(BtoKstarhigh_AL0_err,58)*cexp(I*pi/4.);
	setparam(BtoKstarhigh_AR0_err,59)*cexp(I*pi/4.);
	setparam(BtoKstarhigh_At_err,60)*cexp(I*pi/4.);
	setparam(BtoKstarhigh_AS_err,61)*cexp(I*pi/4.);		

/* B -> K */
	setparam(f_K,62);
	setparam(a1K,63);
	setparam(a2K,64);
	/* Form factors B->K ll */
	/* Form factors B->K* ll */
	/* SN: similar to the switch for BKstar_implementation, I'm doing a BK_FormFactor_choice */
	switch(param.BK_FormFactor_choice)
	{
		default:
		{
			/* SN: Replacing with more explicit naming and adding missing a3f0 - 28/03/2025 */
			// setparam(a00_BK,65);
			// setparam(a10_BK,66);
			// setparam(a20_BK,67);
			// setparam(a0p_BK,68);
			// setparam(a1p_BK,69);
			// setparam(a2p_BK,70);
			// setparam(a0T_BK,71);
			// setparam(a1T_BK,72);
			// setparam(a2T_BK,73);
			/* SN: B->K FF Altmannshofer, Straub 1411.3161 (BK_FormFactor_choice:1)*/
			setparam(a0f0_BK_AS_LCSR_Lattice,192);
			setparam(a1f0_BK_AS_LCSR_Lattice,66);
			setparam(a2f0_BK_AS_LCSR_Lattice,67);
			setparam(a3f0_BK_AS_LCSR_Lattice,193);
			setparam(a0fp_BK_AS_LCSR_Lattice,68);
			setparam(a1fp_BK_AS_LCSR_Lattice,69);
			setparam(a2fp_BK_AS_LCSR_Lattice,70);
			setparam(a0fT_BK_AS_LCSR_Lattice,71);
			setparam(a1fT_BK_AS_LCSR_Lattice,72);
			setparam(a2fT_BK_AS_LCSR_Lattice,73);

			break;
		}

		case 2:
		/* SN: B->K FF van Dyk et al. 2305.06301 (BK_FormFactor_choice:2)*/
		{
			// setparam(a0f0_BK_GRvDV_BSZ,192);
			setparam(a1f0_BK_GRvDV_BSZ,66);
			setparam(a2f0_BK_GRvDV_BSZ,67);
			setparam(a0fp_BK_GRvDV_BSZ,68);
			setparam(a1fp_BK_GRvDV_BSZ,69);
			setparam(a2fp_BK_GRvDV_BSZ,70);
			setparam(a0fT_BK_GRvDV_BSZ,71);
			setparam(a1fT_BK_GRvDV_BSZ,72);
			setparam(a2fT_BK_GRvDV_BSZ,73);

			break;
		}

		case 3:
        /* SN: LCSR+Lattice B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:3)*/
		{
			// setparam(a0f0_BK_GKvD_LCSR_Lattice,192);
			setparam(a1f0_BK_GKvD_LCSR_Lattice,66);
			setparam(a2f0_BK_GKvD_LCSR_Lattice,67);
			setparam(a0fp_BK_GKvD_LCSR_Lattice,68);
			setparam(a1fp_BK_GKvD_LCSR_Lattice,69);
			setparam(a2fp_BK_GKvD_LCSR_Lattice,70);
			setparam(a0fT_BK_GKvD_LCSR_Lattice,71);
			setparam(a1fT_BK_GKvD_LCSR_Lattice,72);
			setparam(a2fT_BK_GKvD_LCSR_Lattice,73);

			break;
		}

		case 4:
        /* SN: LCSR-only B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:4)*/
		{
			// setparam(a0f0_BK_GKvD_LCSRonly,192);
			setparam(a1f0_BK_GKvD_LCSRonly,66);
			setparam(a2f0_BK_GKvD_LCSRonly,67);
			setparam(a0fp_BK_GKvD_LCSRonly,68);
			setparam(a1fp_BK_GKvD_LCSRonly,69);
			setparam(a2fp_BK_GKvD_LCSRonly,70);
			setparam(a0fT_BK_GKvD_LCSRonly,71);
			setparam(a1fT_BK_GKvD_LCSRonly,72);
			setparam(a2fT_BK_GKvD_LCSRonly,73);

			break;
		}

		case 6:
		/* SN: Lattice FLAG24 average 2411.04268 of HPQCD'13 and FNAL/MILC'15 (BK_FormFactor_choice:6) - 26/04/2025*/
		{
			setparam(a0f0_BK_FLAG24,192);
			setparam(a1f0_BK_FLAG24,66);
			// setparam(a2f0_BK_FLAG24,67);
			setparam(a0fp_BK_FLAG24,68);
			setparam(a1fp_BK_FLAG24,69);
			setparam(a2fp_BK_FLAG24,70);
			setparam(a0fT_BK_FLAG24,71);
			setparam(a1fT_BK_FLAG24,72);
			setparam(a2fT_BK_FLAG24,73);

			break;
		}

		case 7:
		/* SN: Lattice full q2-region HPQCD22 2207.12468 (BK_FormFactor_choice:7) - 26/04/2025*/
		{
			setparam(a0f0_BK_HPQCD22,192);
			setparam(a1f0_BK_HPQCD22,66);
			setparam(a2f0_BK_HPQCD22,67);
			// setparam(a0fp_BK_HPQCD22,68);
			setparam(a1fp_BK_HPQCD22,69);
			setparam(a2fp_BK_HPQCD22,70);
			setparam(a0fT_BK_HPQCD22,71);
			setparam(a1fT_BK_HPQCD22,72);
			setparam(a2fT_BK_HPQCD22,73);
			setparam(L_BK_HPQCD22,193);

			break;
		}

		case 8:
		/* SN: LCSR-only Khodjamirian & Rusov 1703.04765 (BK_FormFactor_choice:8) - 02/08/2025*/
		{
			setparam(f0fp_BK_KR,68);
			setparam(b1fp_BK_KR,69);
			setparam(f0fT_BK_KR,71);
			setparam(b1fT_BK_KR,72);

			break;
		}
	}



/* low */
	setparam(BtoKlow_FV_err_noq2,74)*cexp(I*pi/4.);
	setparam(BtoKlow_FA_err_noq2,75)*cexp(I*pi/4.);
	setparam(BtoKlow_FS_err_noq2,76)*cexp(I*pi/4.);
	setparam(BtoKlow_FP_err_noq2,77)*cexp(I*pi/4.);

	setparam(BtoKlow_FV_err_q2,78)*cexp(I*pi/4.);
	setparam(BtoKlow_FA_err_q2,79)*cexp(I*pi/4.);
	setparam(BtoKlow_FS_err_q2,80)*cexp(I*pi/4.);
	setparam(BtoKlow_FP_err_q2,81)*cexp(I*pi/4.);

/* high */
	setparam(BtoKhigh_FV_err,82)*cexp(I*pi/4.);
	setparam(BtoKhigh_FA_err,83)*cexp(I*pi/4.);
	setparam(BtoKhigh_FS_err,84)*cexp(I*pi/4.);
	setparam(BtoKhigh_FP_err,85)*cexp(I*pi/4.);

/* Bs */
	setparam(life_Bs,86);
	setparam(f_Bs,87);
	setparam(lambda_Bsp,88);
	setparam(ys_Bs,89);
	
/* Bs -> phi */
	setparam(f_phi_par,90);
	setparam(f_phi_perp,91);
	setparam(a1phi_perp,92);
	setparam(a1phi_par,93);
	setparam(a2phi_perp,94);
	setparam(a2phi_par,95);
	
/* low */
	setparam(Bstophilow_ALperp_err_noq2,96)*cexp(I*pi/4.);
	setparam(Bstophilow_ARperp_err_noq2,97)*cexp(I*pi/4.);
	setparam(Bstophilow_ALpar_err_noq2,98)*cexp(I*pi/4.);
	setparam(Bstophilow_ARpar_err_noq2,99)*cexp(I*pi/4.);
	setparam(Bstophilow_AL0_err_noq2,100)*cexp(I*pi/4.);
	setparam(Bstophilow_AR0_err_noq2,101)*cexp(I*pi/4.);
	setparam(Bstophilow_At_err_noq2,102)*cexp(I*pi/4.);
	setparam(Bstophilow_AS_err_noq2,103)*cexp(I*pi/4.);

	setparam(Bstophilow_ALperp_err_q2,104)*cexp(I*pi/4.);
	setparam(Bstophilow_ARperp_err_q2,105)*cexp(I*pi/4.);
	setparam(Bstophilow_ALpar_err_q2,106)*cexp(I*pi/4.);
	setparam(Bstophilow_ARpar_err_q2,107)*cexp(I*pi/4.);
	setparam(Bstophilow_AL0_err_q2,108)*cexp(I*pi/4.);
	setparam(Bstophilow_AR0_err_q2,109)*cexp(I*pi/4.);
	setparam(Bstophilow_At_err_q2,110)*cexp(I*pi/4.);
	setparam(Bstophilow_AS_err_q2,111)*cexp(I*pi/4.);	
/* high */
	setparam(Bstophihigh_ALperp_err,112)*cexp(I*pi/4.);
	setparam(Bstophihigh_ARperp_err,113)*cexp(I*pi/4.);
	setparam(Bstophihigh_ALpar_err,114)*cexp(I*pi/4.);
	setparam(Bstophihigh_ARpar_err,115)*cexp(I*pi/4.);
	setparam(Bstophihigh_AL0_err,116)*cexp(I*pi/4.);
	setparam(Bstophihigh_AR0_err,117)*cexp(I*pi/4.);
	setparam(Bstophihigh_At_err,118)*cexp(I*pi/4.);
	setparam(Bstophihigh_AS_err,119)*cexp(I*pi/4.);
			
	/* Form factors B->K* ll */	
	/* SN: similar to the switch for BKstar_implementation, I'm doing a BKstar_FormFactor_choice */
	switch(param.BKstar_FormFactor_choice)
	{
		default:
		{
			setparam(a0A0_BKstar,120);
			setparam(a1A0_BKstar,121);
			setparam(a2A0_BKstar,122);
			setparam(a0A1_BKstar,123);
			setparam(a1A1_BKstar,124);
			setparam(a2A1_BKstar,125);
			setparam(a0A12_BKstar,126);
			setparam(a1A12_BKstar,127);
			setparam(a2A12_BKstar,128);
			setparam(a0V_BKstar,129);
			setparam(a1V_BKstar,130);
			setparam(a2V_BKstar,131);
			setparam(a0T1_BKstar,132);
			setparam(a1T1_BKstar,133);
			setparam(a2T1_BKstar,134);
			setparam(a0T2_BKstar,135);
			setparam(a1T2_BKstar,136);
			setparam(a2T2_BKstar,137);
			setparam(a0T23_BKstar,138);
			setparam(a1T23_BKstar,139);
			setparam(a2T23_BKstar,140);

			break;
		}

		case 2:
		/* SN: B->Kstar FF van Dyk et al. 2305.06301 (BKstar_FormFactor_choice:2)*/
		{
			setparam(a0A0_BKstar_GRvDV_BSZ,120);
			setparam(a1A0_BKstar_GRvDV_BSZ,121);
			setparam(a2A0_BKstar_GRvDV_BSZ,122);
			setparam(a0A1_BKstar_GRvDV_BSZ,123);
			setparam(a1A1_BKstar_GRvDV_BSZ,124);
			setparam(a2A1_BKstar_GRvDV_BSZ,125);
			// setparam(a0A12_BKstar_GRvDV_BSZ,126);
			setparam(a1A12_BKstar_GRvDV_BSZ,127);
			setparam(a2A12_BKstar_GRvDV_BSZ,128);
			setparam(a0V_BKstar_GRvDV_BSZ,129);
			setparam(a1V_BKstar_GRvDV_BSZ,130);
			setparam(a2V_BKstar_GRvDV_BSZ,131);
			setparam(a0T1_BKstar_GRvDV_BSZ,132);
			setparam(a1T1_BKstar_GRvDV_BSZ,133);
			setparam(a2T1_BKstar_GRvDV_BSZ,134);
			// setparam(a0T2_BKstar_GRvDV_BSZ,135);
			setparam(a1T2_BKstar_GRvDV_BSZ,136);
			setparam(a2T2_BKstar_GRvDV_BSZ,137);
			setparam(a0T23_BKstar_GRvDV_BSZ,138);
			setparam(a1T23_BKstar_GRvDV_BSZ,139);
			setparam(a2T23_BKstar_GRvDV_BSZ,140);

			break;
		}

		case 3:
		/* SN: LCSR+Lattice B->Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:3)*/
		{
			setparam(a0A0_BKstar_GKvD_LCSR_Lattice,120);
			setparam(a1A0_BKstar_GKvD_LCSR_Lattice,121);
			setparam(a2A0_BKstar_GKvD_LCSR_Lattice,122);
			setparam(a0A1_BKstar_GKvD_LCSR_Lattice,123);
			setparam(a1A1_BKstar_GKvD_LCSR_Lattice,124);
			setparam(a2A1_BKstar_GKvD_LCSR_Lattice,125);
			// setparam(a0A12_BKstar_GKvD_LCSR_Lattice,126);
			setparam(a1A12_BKstar_GKvD_LCSR_Lattice,127);
			setparam(a2A12_BKstar_GKvD_LCSR_Lattice,128);
			setparam(a0V_BKstar_GKvD_LCSR_Lattice,129);
			setparam(a1V_BKstar_GKvD_LCSR_Lattice,130);
			setparam(a2V_BKstar_GKvD_LCSR_Lattice,131);
			setparam(a0T1_BKstar_GKvD_LCSR_Lattice,132);
			setparam(a1T1_BKstar_GKvD_LCSR_Lattice,133);
			setparam(a2T1_BKstar_GKvD_LCSR_Lattice,134);
			// setparam(a0T2_BKstar_GKvD_LCSR_Lattice,135);
			setparam(a1T2_BKstar_GKvD_LCSR_Lattice,136);
			setparam(a2T2_BKstar_GKvD_LCSR_Lattice,137);
			setparam(a0T23_BKstar_GKvD_LCSR_Lattice,138);
			setparam(a1T23_BKstar_GKvD_LCSR_Lattice,139);
			setparam(a2T23_BKstar_GKvD_LCSR_Lattice,140);

			break;
		}

		case 4:
		/* SN: LCSR-only B->Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:4)*/
		{
			setparam(a0A0_BKstar_GKvD_LCSRonly,120);
			setparam(a1A0_BKstar_GKvD_LCSRonly,121);
			setparam(a2A0_BKstar_GKvD_LCSRonly,122);
			setparam(a0A1_BKstar_GKvD_LCSRonly,123);
			setparam(a1A1_BKstar_GKvD_LCSRonly,124);
			setparam(a2A1_BKstar_GKvD_LCSRonly,125);
			// setparam(a0A12_BKstar_GKvD_LCSRonly,126);
			setparam(a1A12_BKstar_GKvD_LCSRonly,127);
			setparam(a2A12_BKstar_GKvD_LCSRonly,128);
			setparam(a0V_BKstar_GKvD_LCSRonly,129);
			setparam(a1V_BKstar_GKvD_LCSRonly,130);
			setparam(a2V_BKstar_GKvD_LCSRonly,131);
			setparam(a0T1_BKstar_GKvD_LCSRonly,132);
			setparam(a1T1_BKstar_GKvD_LCSRonly,133);
			setparam(a2T1_BKstar_GKvD_LCSRonly,134);
			// setparam(a0T2_BKstar_GKvD_LCSRonly,135);
			setparam(a1T2_BKstar_GKvD_LCSRonly,136);
			setparam(a2T2_BKstar_GKvD_LCSRonly,137);
			setparam(a0T23_BKstar_GKvD_LCSRonly,138);
			setparam(a1T23_BKstar_GKvD_LCSRonly,139);
			setparam(a2T23_BKstar_GKvD_LCSRonly,140);

			break;
		}

		case 5:
		/* SN: LCSR-only B->Kstar FF Zwicky et al. 1503.05534 (BKstar_FormFactor_choice:5)*/
		{
			setparam(a0A0_BKstar_BSZ_LCSRonly,120);
			setparam(a1A0_BKstar_BSZ_LCSRonly,121);
			setparam(a2A0_BKstar_BSZ_LCSRonly,122);
			setparam(a0A1_BKstar_BSZ_LCSRonly,123);
			setparam(a1A1_BKstar_BSZ_LCSRonly,124);
			setparam(a2A1_BKstar_BSZ_LCSRonly,125);
			setparam(a0A12_BKstar_BSZ_LCSRonly,126);
			setparam(a1A12_BKstar_BSZ_LCSRonly,127);
			setparam(a2A12_BKstar_BSZ_LCSRonly,128);
			setparam(a0V_BKstar_BSZ_LCSRonly,129);
			setparam(a1V_BKstar_BSZ_LCSRonly,130);
			setparam(a2V_BKstar_BSZ_LCSRonly,131);
			setparam(a0T1_BKstar_BSZ_LCSRonly,132);
			setparam(a1T1_BKstar_BSZ_LCSRonly,133);
			setparam(a2T1_BKstar_BSZ_LCSRonly,134);
			setparam(a0T2_BKstar_BSZ_LCSRonly,135);
			setparam(a1T2_BKstar_BSZ_LCSRonly,136);
			setparam(a2T2_BKstar_BSZ_LCSRonly,137);
			setparam(a0T23_BKstar_BSZ_LCSRonly,138);
			setparam(a1T23_BKstar_BSZ_LCSRonly,139);
			setparam(a2T23_BKstar_BSZ_LCSRonly,140);

			break;
		}

		case 6:
		/* SN: Lattice-only B->Kstar FF Wingate et al. 1310.3722 + 1501.00367 (BKstar_FormFactor_choice:6*/
		{
			setparam(a0A0_BKstar_HLMW,120);
			setparam(a1A0_BKstar_HLMW,121);
			setparam(a0A1_BKstar_HLMW,123);
			setparam(a1A1_BKstar_HLMW,124);
			setparam(a0A12_BKstar_HLMW,126);
			setparam(a1A12_BKstar_HLMW,127);
			setparam(a0V_BKstar_HLMW,129);
			setparam(a1V_BKstar_HLMW,130);
			setparam(a0T1_BKstar_HLMW,132);
			setparam(a1T1_BKstar_HLMW,133);
			setparam(a0T2_BKstar_HLMW,135);
			setparam(a1T2_BKstar_HLMW,136);
			setparam(a0T23_BKstar_HLMW,138);
			setparam(a1T23_BKstar_HLMW,139);
			setparam(systErr_BKstar_HLMW,122);

			break;
		}
	}


	/* Form factors Bs->phi ll */	
	/* SN: similar to the switch for BKstar_implementation, I'm doing a Bsphi_FormFactor_choice */
	switch(param.Bsphi_FormFactor_choice)
	{
		default:
		{
			setparam(a0A0_Bsphi,141);
			setparam(a1A0_Bsphi,142);
			setparam(a2A0_Bsphi,143);
			setparam(a0A1_Bsphi,144);
			setparam(a1A1_Bsphi,145);
			setparam(a2A1_Bsphi,146);
			setparam(a0A12_Bsphi,147);
			setparam(a1A12_Bsphi,148);
			setparam(a2A12_Bsphi,149);
			setparam(a0V_Bsphi,150);
			setparam(a1V_Bsphi,151);
			setparam(a2V_Bsphi,152);
			setparam(a0T1_Bsphi,153);
			setparam(a1T1_Bsphi,154);
			setparam(a2T1_Bsphi,155);
			setparam(a0T2_Bsphi,156);
			setparam(a1T2_Bsphi,157);
			setparam(a2T2_Bsphi,158);
			setparam(a0T23_Bsphi,159);
			setparam(a1T23_Bsphi,160);
			setparam(a2T23_Bsphi,161);

			break;
		}

		case 2:
		/* SN: Bs->phi FF van Dyk et al. 2305.06301 (Bsphi_FormFactor_choice:2)*/
		{
			setparam(a0A0_Bsphi_GRvDV_BSZ,141);
			setparam(a1A0_Bsphi_GRvDV_BSZ,142);
			setparam(a2A0_Bsphi_GRvDV_BSZ,143);
			setparam(a0A1_Bsphi_GRvDV_BSZ,144);
			setparam(a1A1_Bsphi_GRvDV_BSZ,145);
			setparam(a2A1_Bsphi_GRvDV_BSZ,146);
			// setparam(a0A12_Bsphi_GRvDV_BSZ,147);
			setparam(a1A12_Bsphi_GRvDV_BSZ,148);
			setparam(a2A12_Bsphi_GRvDV_BSZ,149);
			setparam(a0V_Bsphi_GRvDV_BSZ,150);
			setparam(a1V_Bsphi_GRvDV_BSZ,151);
			setparam(a2V_Bsphi_GRvDV_BSZ,152);
			setparam(a0T1_Bsphi_GRvDV_BSZ,153);
			setparam(a1T1_Bsphi_GRvDV_BSZ,154);
			setparam(a2T1_Bsphi_GRvDV_BSZ,155);
			// setparam(a0T2_Bsphi_GRvDV_BSZ,156);
			setparam(a1T2_Bsphi_GRvDV_BSZ,157);
			setparam(a2T2_Bsphi_GRvDV_BSZ,158);
			setparam(a0T23_Bsphi_GRvDV_BSZ,159);
			setparam(a1T23_Bsphi_GRvDV_BSZ,160);
			setparam(a2T23_Bsphi_GRvDV_BSZ,161);
			break;
		}

		case 5:
		/* SN: Bs->phi FF Zwicky et al. 1503.05534 (Bsphi_FormFactor_choice:5)*/
		{
			setparam(a0A0_Bsphi_BSZ_LCSRonly,141);
			setparam(a1A0_Bsphi_BSZ_LCSRonly,142);
			setparam(a2A0_Bsphi_BSZ_LCSRonly,143);
			setparam(a0A1_Bsphi_BSZ_LCSRonly,144);
			setparam(a1A1_Bsphi_BSZ_LCSRonly,145);
			setparam(a2A1_Bsphi_BSZ_LCSRonly,146);
			setparam(a0A12_Bsphi_BSZ_LCSRonly,147);
			setparam(a1A12_Bsphi_BSZ_LCSRonly,148);
			setparam(a2A12_Bsphi_BSZ_LCSRonly,149);
			setparam(a0V_Bsphi_BSZ_LCSRonly,150);
			setparam(a1V_Bsphi_BSZ_LCSRonly,151);
			setparam(a2V_Bsphi_BSZ_LCSRonly,152);
			setparam(a0T1_Bsphi_BSZ_LCSRonly,153);
			setparam(a1T1_Bsphi_BSZ_LCSRonly,154);
			setparam(a2T1_Bsphi_BSZ_LCSRonly,155);
			setparam(a0T2_Bsphi_BSZ_LCSRonly,156);
			setparam(a1T2_Bsphi_BSZ_LCSRonly,157);
			setparam(a2T2_Bsphi_BSZ_LCSRonly,158);
			setparam(a0T23_Bsphi_BSZ_LCSRonly,159);
			setparam(a1T23_Bsphi_BSZ_LCSRonly,160);
			setparam(a2T23_Bsphi_BSZ_LCSRonly,161);
			break;
		}

		case 6:
		/* SN: Lattice-only Bs->phi FF Wingate et al. 1310.3722 + 1501.00367 (Bsphi_FormFactor_choice:6*/
		{
			setparam(a0A0_Bsphi_HLMW,120);
			setparam(a1A0_Bsphi_HLMW,121);
			setparam(a0A1_Bsphi_HLMW,123);
			setparam(a1A1_Bsphi_HLMW,124);
			setparam(a0A12_Bsphi_HLMW,126);
			setparam(a1A12_Bsphi_HLMW,127);
			setparam(a0V_Bsphi_HLMW,129);
			setparam(a1V_Bsphi_HLMW,130);
			setparam(a0T1_Bsphi_HLMW,132);
			setparam(a1T1_Bsphi_HLMW,133);
			setparam(a0T2_Bsphi_HLMW,135);
			setparam(a1T2_Bsphi_HLMW,136);
			setparam(a0T23_Bsphi_HLMW,138);
			setparam(a1T23_Bsphi_HLMW,139);
			setparam(systErr_Bsphi_HLMW,122);

			break;
		}
	}

	/* Lambda_b -> Lambda l+l- */	
	setparam(life_Lb,162);
    setparam(alphaL_LbLll,163);
    setparam(a0_HO_fplus_LbLll,164);
    setparam(a1_HO_fplus_LbLll,165);
    setparam(a2_HO_fplus_LbLll,166);
    setparam(a0_HO_fperp_LbLll,167);
    setparam(a1_HO_fperp_LbLll,168);
    setparam(a2_HO_fperp_LbLll,169);
    setparam(a0_HO_gpp_LbLll,170);
    setparam(a1_HO_gplus_LbLll,171);
    setparam(a2_HO_gplus_LbLll,172);
    setparam(a1_HO_gperp_LbLll,173);
    setparam(a2_HO_gperp_LbLll,174);
    setparam(a0_HO_hplus_LbLll,175);
    setparam(a1_HO_hplus_LbLll,176);
    setparam(a2_HO_hplus_LbLll,177);
    setparam(a0_HO_hperp_LbLll,178);
    setparam(a1_HO_hperp_LbLll,179);
    setparam(a2_HO_hperp_LbLll,180);
    setparam(a0_HO_htildepp_LbLll,181);
    setparam(a1_HO_htildeplus_LbLll,182);
    setparam(a2_HO_htildeplus_LbLll,183);
    setparam(a1_HO_htildeperp_LbLll,184);
    setparam(a2_HO_htildeperp_LbLll,185);
	setparam(deltaPcu_Kppipnunu,186);
	setparam(err_Pc_Xlambda_Kppipnunu,187);
	setparam(BR_KLgammagamma_exp,188);
	setparam(Aterm_mu_KLmumu,189);
	setparam(chi_gg_Mrho,189);// Replaces Aterm_mu_KLmumu in the generalized calculation to also be used for electron mode
	setparam(BR_KSgammagamma_exp,190);
	setparam(Iterm_mu_KSmumu,191);
	setparam(log_muK_1GeV,196);
	setparam(KLpill_Ce_Cdir,197);
    setparam(KLpill_Cmu_Cdir,198);
	setparam(KLpill_Ce_Cint,199);
	setparam(KLpill_Cmu_Cint,200);
	setparam(KLpill_Ce_Cmix,201);
	setparam(KLpill_Cmu_Cmix,202);
	setparam(KLpill_Cmu_CPC,203);
	setparam(KLpill_abs_aS,204);
	/* Scalar and Pseudoscalar contributions included - 08/10/2025 */
	setparam(KLpill_Ce_Sgg,205);
	setparam(KLpill_Cmu_Sgg,206);
	/* AFBmu_KLpi0mumu  - 09/10/2025 */
	setparam(KLpill_AFB_a1pi_gg_7V,207);
	setparam(KLpill_AFB_a1pi_gg_aS,208);

	if(!param.fullFF)
	{
		/* Form factors B->K* ll */	
		setstdzero(a0A0_BKstar);
		setstdzero(a1A0_BKstar);
		setstdzero(a2A0_BKstar);
		setstdzero(a0T1_BKstar);
		setstdzero(a1T1_BKstar);
		setstdzero(a2T1_BKstar);
		setstdzero(a0T2_BKstar);
		setstdzero(a1T2_BKstar);
		setstdzero(a2T2_BKstar);
		setstdzero(a0T23_BKstar);
		setstdzero(a1T23_BKstar);
		setstdzero(a2T23_BKstar);
		
		/* Form factors Bs->phi ll */	
		setstdzero(a0A0_Bsphi);
		setstdzero(a1A0_Bsphi);
		setstdzero(a2A0_Bsphi);
		setstdzero(a0T1_Bsphi);
		setstdzero(a1T1_Bsphi);
		setstdzero(a2T1_Bsphi);
		setstdzero(a0T2_Bsphi);
		setstdzero(a1T2_Bsphi);
		setstdzero(a2T2_Bsphi);
		setstdzero(a0T23_Bsphi);
		setstdzero(a1T23_Bsphi);
		setstdzero(a2T23_Bsphi);
		
	}

#undef setparam
#undef setstdzero

	if(nobs==0) return;

#define mystrcmp( myparam ) (!strcmp(namenuisance[ke], #myparam))
	
	int test_general=(ke==0||mystrcmp(alphas_MZ)||mystrcmp(mass_b)||mystrcmp(mass_c)||mystrcmp(mass_s)||mystrcmp(mass_top_pole)||mystrcmp(mass_h0)||mystrcmp(CKM_lambda)||mystrcmp(CKM_A)||mystrcmp(CKM_rhobar)||mystrcmp(CKM_etabar)||mystrcmp(log_mu_W_mass_W)||mystrcmp(log_mu_b_mass_b));

	int test_delta0Bkstargamma=(mystrcmp(f_B)||mystrcmp(lambda_Bp)||mystrcmp(f_Kstar_par)||mystrcmp(f_Kstar_perp)||mystrcmp(a1perp)||mystrcmp(a2perp)||mystrcmp(a1par)||mystrcmp(a2par)||mystrcmp(T1_BKstar)||mystrcmp(log_mu_spec_lambda_h_mass_b));
	
	int test_BKstargamma=(mystrcmp(f_B)||mystrcmp(lambda_Bp)||mystrcmp(f_Kstar_par)||mystrcmp(f_Kstar_perp)||mystrcmp(a1perp)||mystrcmp(a2perp)||mystrcmp(a1par)||mystrcmp(a2par)||mystrcmp(T1_BKstar));
	
	int test_bsgamma=(mystrcmp(BR_BXclnu_exp)||mystrcmp(mu_G2_bsg)||mystrcmp(rho_D3_bsg)||mystrcmp(rho_LS3_bsg)||mystrcmp(bsgamma_rand)||mystrcmp(mu_c_bsg));
	
	int test_Bsmumu=(mystrcmp(life_Bs)||mystrcmp(f_Bs)||mystrcmp(lambda_Bsp)||mystrcmp(ys_Bs));
	
	int test_BXsll=(mystrcmp(BR_BXclnu_exp)||mystrcmp(BRBXsmumu_lowq2_rand)||mystrcmp(BRBXsmumu_highq2_rand)||mystrcmp(BRBXsmumu_full_rand)||mystrcmp(BRBXsee_lowq2_rand)||mystrcmp(BRBXsee_highq2_rand)||mystrcmp(BRBXsee_full_rand)||mystrcmp(BRBXstautau_lowq2_rand)||mystrcmp(BRBXstautau_highq2_rand)||mystrcmp(BRBXstautau_full_rand));

	/*------------------ B->Kstarll -----------------*/
	/* General uncertainties, considered in the full q2 region */
	int test_BKstar=(mystrcmp(f_B)||mystrcmp(lambda_Bp)||mystrcmp(f_Kstar_par)||mystrcmp(f_Kstar_perp)||mystrcmp(a1perp)||mystrcmp(a2perp)||mystrcmp(a1par)||mystrcmp(a2par)
	/* B->Kstar FF Zwicky et al. 1503.05534 (BKstar_FormFactor_choice:1)*/
	||mystrcmp(a0A0_BKstar)||mystrcmp(a1A0_BKstar)||mystrcmp(a2A0_BKstar)||mystrcmp(a0A1_BKstar)||mystrcmp(a1A1_BKstar)||mystrcmp(a2A1_BKstar)||mystrcmp(a0A12_BKstar)||mystrcmp(a1A12_BKstar)||mystrcmp(a2A12_BKstar)||mystrcmp(a0V_BKstar)||mystrcmp(a1V_BKstar)||mystrcmp(a2V_BKstar)||mystrcmp(a0T1_BKstar)||mystrcmp(a1T1_BKstar)||mystrcmp(a2T1_BKstar)||mystrcmp(a0T2_BKstar)||mystrcmp(a1T2_BKstar)||mystrcmp(a2T2_BKstar)||mystrcmp(a0T23_BKstar)||mystrcmp(a1T23_BKstar)||mystrcmp(a2T23_BKstar)
	/* SN: B->Kstar FF van Dyk et al. 2305.06301 (BKstar_FormFactor_choice:2)*/
	||mystrcmp(a0A0_BKstar_GRvDV_BSZ)||mystrcmp(a1A0_BKstar_GRvDV_BSZ)||mystrcmp(a2A0_BKstar_GRvDV_BSZ)||mystrcmp(a0A1_BKstar_GRvDV_BSZ)||mystrcmp(a1A1_BKstar_GRvDV_BSZ)||mystrcmp(a2A1_BKstar_GRvDV_BSZ)||mystrcmp(a1A12_BKstar_GRvDV_BSZ)||mystrcmp(a2A12_BKstar_GRvDV_BSZ)||mystrcmp(a0V_BKstar_GRvDV_BSZ)||mystrcmp(a1V_BKstar_GRvDV_BSZ)||mystrcmp(a2V_BKstar_GRvDV_BSZ)||mystrcmp(a0T1_BKstar_GRvDV_BSZ)||mystrcmp(a1T1_BKstar_GRvDV_BSZ)||mystrcmp(a2T1_BKstar_GRvDV_BSZ)||mystrcmp(a0T2_BKstar_GRvDV_BSZ)||mystrcmp(a1T2_BKstar_GRvDV_BSZ)||mystrcmp(a2T2_BKstar_GRvDV_BSZ)||mystrcmp(a0T23_BKstar_GRvDV_BSZ)||mystrcmp(a1T23_BKstar_GRvDV_BSZ)||mystrcmp(a2T23_BKstar_GRvDV_BSZ)
	/* SN: LCSR+Lattice B->Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:3)*/
	||mystrcmp(a0A0_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a1A0_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a2A0_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a0A1_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a1A1_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a2A1_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a1A12_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a2A12_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a0V_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a1V_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a2V_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a0T1_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a1T1_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a2T1_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a0T2_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a1T2_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a2T2_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a0T23_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a1T23_BKstar_GKvD_LCSR_Lattice)||mystrcmp(a2T23_BKstar_GKvD_LCSR_Lattice)
	/* SN: LCSR-only B->Kstar FF van Dyk et al. 1811.00983 (BKstar_FormFactor_choice:4)*/
	||mystrcmp(a0A0_BKstar_GKvD_LCSRonly)||mystrcmp(a1A0_BKstar_GKvD_LCSRonly)||mystrcmp(a2A0_BKstar_GKvD_LCSRonly)||mystrcmp(a0A1_BKstar_GKvD_LCSRonly)||mystrcmp(a1A1_BKstar_GKvD_LCSRonly)||mystrcmp(a2A1_BKstar_GKvD_LCSRonly)||mystrcmp(a1A12_BKstar_GKvD_LCSRonly)||mystrcmp(a2A12_BKstar_GKvD_LCSRonly)||mystrcmp(a0V_BKstar_GKvD_LCSRonly)||mystrcmp(a1V_BKstar_GKvD_LCSRonly)||mystrcmp(a2V_BKstar_GKvD_LCSRonly)||mystrcmp(a0T1_BKstar_GKvD_LCSRonly)||mystrcmp(a1T1_BKstar_GKvD_LCSRonly)||mystrcmp(a2T1_BKstar_GKvD_LCSRonly)||mystrcmp(a0T2_BKstar_GKvD_LCSRonly)||mystrcmp(a1T2_BKstar_GKvD_LCSRonly)||mystrcmp(a2T2_BKstar_GKvD_LCSRonly)||mystrcmp(a0T23_BKstar_GKvD_LCSRonly)||mystrcmp(a1T23_BKstar_GKvD_LCSRonly)||mystrcmp(a2T23_BKstar_GKvD_LCSRonly)
	/* SN: LCSR-only B->Kstar FF Zwicky et al. 1503.05534 (BKstar_FormFactor_choice:5)*/
	||mystrcmp(a0A0_BKstar_BSZ_LCSRonly)||mystrcmp(a1A0_BKstar_BSZ_LCSRonly)||mystrcmp(a2A0_BKstar_BSZ_LCSRonly)||mystrcmp(a0A1_BKstar_BSZ_LCSRonly)||mystrcmp(a1A1_BKstar_BSZ_LCSRonly)||mystrcmp(a2A1_BKstar_BSZ_LCSRonly)||mystrcmp(a0A12_BKstar_BSZ_LCSRonly)||mystrcmp(a1A12_BKstar_BSZ_LCSRonly)||mystrcmp(a2A12_BKstar_BSZ_LCSRonly)||mystrcmp(a0V_BKstar_BSZ_LCSRonly)||mystrcmp(a1V_BKstar_BSZ_LCSRonly)||mystrcmp(a2V_BKstar_BSZ_LCSRonly)||mystrcmp(a0T1_BKstar_BSZ_LCSRonly)||mystrcmp(a1T1_BKstar_BSZ_LCSRonly)||mystrcmp(a2T1_BKstar_BSZ_LCSRonly)||mystrcmp(a0T2_BKstar_BSZ_LCSRonly)||mystrcmp(a1T2_BKstar_BSZ_LCSRonly)||mystrcmp(a2T2_BKstar_BSZ_LCSRonly)||mystrcmp(a0T23_BKstar_BSZ_LCSRonly)||mystrcmp(a1T23_BKstar_BSZ_LCSRonly)||mystrcmp(a2T23_BKstar_BSZ_LCSRonly)
	/* SN: Lattice-only B->Kstar FF Wingate et al. 1310.3722 + 1501.00367 (BKstar_FormFactor_choice:6*/
	||mystrcmp(a0A0_BKstar_HLMW)||mystrcmp(a1A0_BKstar_HLMW)||mystrcmp(a0A1_BKstar_HLMW)||mystrcmp(a1A1_BKstar_HLMW)||mystrcmp(a0A12_BKstar_HLMW)||mystrcmp(a1A12_BKstar_HLMW)||mystrcmp(a0V_BKstar_HLMW)||mystrcmp(a1V_BKstar_HLMW)||mystrcmp(a0T1_BKstar_HLMW)||mystrcmp(a1T1_BKstar_HLMW)||mystrcmp(a0T2_BKstar_HLMW)||mystrcmp(a1T2_BKstar_HLMW)||mystrcmp(a0T23_BKstar_HLMW)||mystrcmp(a1T23_BKstar_HLMW)||mystrcmp(systErr_BKstar_HLMW)
	/* SN: Result from B->K* Hadronic Fit to be used to check how it impacts predictions using BKstar_implementation:4
	 * with modified nuiscance.c and nuisance_corr.c containting the fitted values  - 13/04/2025 */
	||mystrcmp(hplus0r)||mystrcmp(hminus0r)||mystrcmp(hplus1r)||mystrcmp(hminus1r)||mystrcmp(hplus2r)||mystrcmp(hminus2r)||mystrcmp(hzero0r)||mystrcmp(hzero1r)||mystrcmp(hzero2r)||mystrcmp(hplus0i)||mystrcmp(hminus0i)||mystrcmp(hplus1i)||mystrcmp(hminus1i)||mystrcmp(hplus2i)||mystrcmp(hminus2i)||mystrcmp(hzero0i)||mystrcmp(hzero1i)||mystrcmp(hzero2i)
	);

	/* Uncertainties considered in the low q2 region */
	int test_lowBKstar=(test_BKstar
	||mystrcmp(BtoKstarlow_ALperp_err_noq2)||mystrcmp(BtoKstarlow_ARperp_err_noq2)||mystrcmp(BtoKstarlow_ALpar_err_noq2)||mystrcmp(BtoKstarlow_ARpar_err_noq2)||mystrcmp(BtoKstarlow_AL0_err_noq2)||mystrcmp(BtoKstarlow_AR0_err_noq2)||mystrcmp(BtoKstarlow_At_err_noq2)||mystrcmp(BtoKstarlow_AS_err_noq2)||mystrcmp(BtoKstarlow_ALperp_err_q2)||mystrcmp(BtoKstarlow_ARperp_err_q2)||mystrcmp(BtoKstarlow_ALpar_err_q2)||mystrcmp(BtoKstarlow_ARpar_err_q2)||mystrcmp(BtoKstarlow_AL0_err_q2)||mystrcmp(BtoKstarlow_AR0_err_q2)||mystrcmp(BtoKstarlow_At_err_q2)||mystrcmp(BtoKstarlow_AS_err_q2)
	//
	||mystrcmp(real_alpha_perp0)||mystrcmp(real_alpha_perp1)||mystrcmp(real_alpha_perp2)||mystrcmp(real_alpha_par0)||mystrcmp(real_alpha_par1)||mystrcmp(real_alpha_par2)||mystrcmp(real_alpha_zero0)||mystrcmp(real_alpha_zero1)||mystrcmp(imag_alpha_perp0)||mystrcmp(imag_alpha_perp1)||mystrcmp(imag_alpha_perp2)||mystrcmp(imag_alpha_par0)||mystrcmp(imag_alpha_par1)||mystrcmp(imag_alpha_par2)||mystrcmp(imag_alpha_zero0)||mystrcmp(imag_alpha_zero1)
	//
	||mystrcmp(DeltaC9_M1_q2bar)||mystrcmp(r1_M1)||mystrcmp(r2_M1)||mystrcmp(DeltaC9_M2_q2bar)||mystrcmp(r1_M2)||mystrcmp(r2_M2)||mystrcmp(DeltaC9_M3_q2bar)||mystrcmp(r1_M3)||mystrcmp(r2_M3)
	);

	/* Uncertainties considered in the intermediate q2 region */
	int test_intermBKstar=(test_BKstar
	||mystrcmp(BtoKstarlow_ALperp_err_noq2)||mystrcmp(BtoKstarlow_ARperp_err_noq2)||mystrcmp(BtoKstarlow_ALpar_err_noq2)||mystrcmp(BtoKstarlow_ARpar_err_noq2)||mystrcmp(BtoKstarlow_AL0_err_noq2)||mystrcmp(BtoKstarlow_AR0_err_noq2)||mystrcmp(BtoKstarlow_At_err_noq2)||mystrcmp(BtoKstarlow_AS_err_noq2)||mystrcmp(BtoKstarlow_ALperp_err_q2)||mystrcmp(BtoKstarlow_ARperp_err_q2)||mystrcmp(BtoKstarlow_ALpar_err_q2)||mystrcmp(BtoKstarlow_ARpar_err_q2)||mystrcmp(BtoKstarlow_AL0_err_q2)||mystrcmp(BtoKstarlow_AR0_err_q2)||mystrcmp(BtoKstarlow_At_err_q2)||mystrcmp(BtoKstarlow_AS_err_q2)
	//
	||mystrcmp(real_alpha_perp0)||mystrcmp(real_alpha_perp1)||mystrcmp(real_alpha_perp2)||mystrcmp(real_alpha_par0)||mystrcmp(real_alpha_par1)||mystrcmp(real_alpha_par2)||mystrcmp(real_alpha_zero0)||mystrcmp(real_alpha_zero1)||mystrcmp(imag_alpha_perp0)||mystrcmp(imag_alpha_perp1)||mystrcmp(imag_alpha_perp2)||mystrcmp(imag_alpha_par0)||mystrcmp(imag_alpha_par1)||mystrcmp(imag_alpha_par2)||mystrcmp(imag_alpha_zero0)||mystrcmp(imag_alpha_zero1)
	//
	||mystrcmp(DeltaC9_M1_q2bar)||mystrcmp(r1_M1)||mystrcmp(r2_M1)||mystrcmp(DeltaC9_M2_q2bar)||mystrcmp(r1_M2)||mystrcmp(r2_M2)||mystrcmp(DeltaC9_M3_q2bar)||mystrcmp(r1_M3)||mystrcmp(r2_M3)
	//
	||mystrcmp(BtoKstarhigh_ALperp_err)||mystrcmp(BtoKstarhigh_ARperp_err)||mystrcmp(BtoKstarhigh_ALpar_err)||mystrcmp(BtoKstarhigh_ARpar_err)||mystrcmp(BtoKstarhigh_AL0_err)||mystrcmp(BtoKstarhigh_AR0_err)||mystrcmp(BtoKstarhigh_At_err)||mystrcmp(BtoKstarhigh_AS_err));

	/* Uncertainties considered in the high q2 region */
	int test_highBKstar=(test_BKstar
	||mystrcmp(BtoKstarhigh_ALperp_err)||mystrcmp(BtoKstarhigh_ARperp_err)||mystrcmp(BtoKstarhigh_ALpar_err)||mystrcmp(BtoKstarhigh_ARpar_err)||mystrcmp(BtoKstarhigh_AL0_err)||mystrcmp(BtoKstarhigh_AR0_err)||mystrcmp(BtoKstarhigh_At_err)||mystrcmp(BtoKstarhigh_AS_err));

	/*------------------ B->Kll --------------------*/
	/* General uncertainties, considered in the full q2 region */
	int test_BK=(mystrcmp(f_B)||mystrcmp(lambda_Bp)
	// For croschecking form factos, all parameters are for the moment in the overall bins -  should be changed
	/* SN: B->K FF Altmannshofer, Straub 1411.3161 (BK_FormFactor_choice:1)*/
	||mystrcmp(a0f0_BK_AS_LCSR_Lattice)||mystrcmp(a1f0_BK_AS_LCSR_Lattice)||mystrcmp(a2f0_BK_AS_LCSR_Lattice)||mystrcmp(a3fp_BK_AS_LCSR_Lattice)||mystrcmp(a0fp_BK_AS_LCSR_Lattice)||mystrcmp(a1fp_BK_AS_LCSR_Lattice)||mystrcmp(a2fp_BK_AS_LCSR_Lattice)||mystrcmp(a0fT_BK_AS_LCSR_Lattice)||mystrcmp(a1fT_BK_AS_LCSR_Lattice)||mystrcmp(a2fT_BK_AS_LCSR_Lattice)
	/* SN: B->K FF van Dyk et al. 2305.06301 (BK_FormFactor_choice:2)*/
	||mystrcmp(a1f0_BK_GRvDV_BSZ)||mystrcmp(a2f0_BK_GRvDV_BSZ)||mystrcmp(a0fp_BK_GRvDV_BSZ)||mystrcmp(a1fp_BK_GRvDV_BSZ)||mystrcmp(a2fp_BK_GRvDV_BSZ)||mystrcmp(a0fT_BK_GRvDV_BSZ)||mystrcmp(a1fT_BK_GRvDV_BSZ)||mystrcmp(a2fT_BK_GRvDV_BSZ)
	/* SN: LCSR+Lattice B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:3)*/
	||mystrcmp(a1f0_BK_GKvD_LCSR_Lattice)||mystrcmp(a2f0_BK_GKvD_LCSR_Lattice)||mystrcmp(a0fp_BK_GKvD_LCSR_Lattice)||mystrcmp(a1fp_BK_GKvD_LCSR_Lattice)||mystrcmp(a2fp_BK_GKvD_LCSR_Lattice)||mystrcmp(a0fT_BK_GKvD_LCSR_Lattice)||mystrcmp(a1fT_BK_GKvD_LCSR_Lattice)||mystrcmp(a2fT_BK_GKvD_LCSR_Lattice)
	/* SN: LCSR-only B->K FF van Dyk et al. 1811.00983 (BK_FormFactor_choice:4)*/
	||mystrcmp(a1f0_BK_GKvD_LCSRonly)||mystrcmp(a2f0_BK_GKvD_LCSRonly)||mystrcmp(a0fp_BK_GKvD_LCSRonly)||mystrcmp(a1fp_BK_GKvD_LCSRonly)||mystrcmp(a2fp_BK_GKvD_LCSRonly)||mystrcmp(a0fT_BK_GKvD_LCSRonly)||mystrcmp(a1fT_BK_GKvD_LCSRonly)||mystrcmp(a2fT_BK_GKvD_LCSRonly)
		/* SN: Lattice FLAG24 average 2411.04268 of HPQCD'13 and FNAL/MILC'15 (BK_FormFactor_choice:6) - 26/04/2025*/
	||mystrcmp(a0fp_BK_FLAG24)||mystrcmp(a1fp_BK_FLAG24)||mystrcmp(a2fp_BK_FLAG24)||mystrcmp(a0f0_BK_FLAG24)||mystrcmp(a1f0_BK_FLAG24)||mystrcmp(a0fT_BK_FLAG24)||mystrcmp(a1fT_BK_FLAG24)||mystrcmp(a2fT_BK_FLAG24)
	/* SN: Lattice full q2-region HPQCD22 2207.12468 (BK_FormFactor_choice:7) - 26/04/2025*/
	||mystrcmp(a0f0_BK_HPQCD22)||mystrcmp(a1f0_BK_HPQCD22)||mystrcmp(a2f0_BK_HPQCD22)||mystrcmp(a1fp_BK_HPQCD22)||mystrcmp(a2fp_BK_HPQCD22)||mystrcmp(a0fT_BK_HPQCD22)||mystrcmp(a1fT_BK_HPQCD22)||mystrcmp(a2fT_BK_HPQCD22)||mystrcmp(L_BK_HPQCD22)
	/* SN: LCSR-only Khodjamirian & Rusov 1703.04765 (BK_FormFactor_choice:8) - 02/08/2025*/
	||mystrcmp(f0fp_BK_KR)||mystrcmp(b1fp_BK_KR)||mystrcmp(f0fT_BK_KR)||mystrcmp(b1fT_BK_KR)
	);

	/* Uncertainties considered in the low q2 region */
	int test_lowBK=test_BK||(mystrcmp(f_K)||mystrcmp(a1K)||mystrcmp(a2K)
	//
	||mystrcmp(BtoKlow_FV_err_noq2)||mystrcmp(BtoKlow_FA_err_noq2)||mystrcmp(BtoKlow_FS_err_noq2)||mystrcmp(BtoKlow_FP_err_noq2)||mystrcmp(BtoKlow_FV_err_q2)||mystrcmp(BtoKlow_FA_err_q2)||mystrcmp(BtoKlow_FS_err_q2)||mystrcmp(BtoKlow_FP_err_q2));

	/* Uncertainties considered in the intermediate q2 region */
	int test_intermBK=test_BK||(mystrcmp(f_K)||mystrcmp(a1K)||mystrcmp(a2K)
	//
	||mystrcmp(BtoKlow_FV_err_noq2)||mystrcmp(BtoKlow_FA_err_noq2)||mystrcmp(BtoKlow_FS_err_noq2)||mystrcmp(BtoKlow_FP_err_noq2)||mystrcmp(BtoKlow_FV_err_q2)||mystrcmp(BtoKlow_FA_err_q2)||mystrcmp(BtoKlow_FS_err_q2)||mystrcmp(BtoKlow_FP_err_q2)
	//
	||mystrcmp(BtoKhigh_FV_err)||mystrcmp(BtoKhigh_FA_err)||mystrcmp(BtoKhigh_FS_err)||mystrcmp(BtoKhigh_FP_err));

	/* Uncertainties considered in the high q2 region */
	int test_highBK=test_BK||(mystrcmp(f_K)||mystrcmp(a1K)
	||mystrcmp(BtoKhigh_FV_err)||mystrcmp(BtoKhigh_FA_err)||mystrcmp(BtoKhigh_FS_err)||mystrcmp(BtoKhigh_FP_err));

	/*------------------ Bs->phill -----------------*/
	/* General uncertainties, considered in the full q2 region */
	int test_Bsphi=(mystrcmp(life_Bs)||mystrcmp(f_Bs)||mystrcmp(lambda_Bsp)||mystrcmp(ys_Bs)||mystrcmp(f_phi_par)||mystrcmp(f_phi_perp)||mystrcmp(a1phi_perp)||mystrcmp(a1phi_par)||mystrcmp(a2phi_perp)||mystrcmp(a2phi_par)
	/* Bs->phi FF BSZ 1503.05534 (Bsphi_FormFactor_choice:1)*/
	||mystrcmp(a0A0_Bsphi)||mystrcmp(a1A0_Bsphi)||mystrcmp(a2A0_Bsphi)||mystrcmp(a0A1_Bsphi)||mystrcmp(a1A1_Bsphi)||mystrcmp(a2A1_Bsphi)||mystrcmp(a0A12_Bsphi)||mystrcmp(a1A12_Bsphi)||mystrcmp(a2A12_Bsphi)||mystrcmp(a0V_Bsphi)||mystrcmp(a1V_Bsphi)||mystrcmp(a2V_Bsphi)||mystrcmp(a0T1_Bsphi)||mystrcmp(a1T1_Bsphi)||mystrcmp(a2T1_Bsphi)||mystrcmp(a0T2_Bsphi)||mystrcmp(a1T2_Bsphi)||mystrcmp(a2T2_Bsphi)||mystrcmp(a0T23_Bsphi)||mystrcmp(a1T23_Bsphi)||mystrcmp(a2T23_Bsphi)
	/* SN: Bs -> phi FF van Dyk et al. 2305.06301 (Bsphi_FormFactor_choice:2)*/
	||mystrcmp(a0A0_Bsphi_GRvDV_BSZ)||mystrcmp(a1A0_Bsphi_GRvDV_BSZ)||mystrcmp(a2A0_Bsphi_GRvDV_BSZ)||mystrcmp(a0A1_Bsphi_GRvDV_BSZ)||mystrcmp(a1A1_Bsphi_GRvDV_BSZ)||mystrcmp(a2A1_Bsphi_GRvDV_BSZ)||mystrcmp(a1A12_Bsphi_GRvDV_BSZ)||mystrcmp(a2A12_Bsphi_GRvDV_BSZ)||mystrcmp(a0V_Bsphi_GRvDV_BSZ)||mystrcmp(a1V_Bsphi_GRvDV_BSZ)||mystrcmp(a2V_Bsphi_GRvDV_BSZ)||mystrcmp(a0T1_Bsphi_GRvDV_BSZ)||mystrcmp(a1T1_Bsphi_GRvDV_BSZ)||mystrcmp(a2T1_Bsphi_GRvDV_BSZ)||mystrcmp(a0T2_Bsphi_GRvDV_BSZ)||mystrcmp(a1T2_Bsphi_GRvDV_BSZ)||mystrcmp(a2T2_Bsphi_GRvDV_BSZ)||mystrcmp(a0T23_Bsphi_GRvDV_BSZ)||mystrcmp(a1T23_Bsphi_GRvDV_BSZ)||mystrcmp(a2T23_Bsphi_GRvDV_BSZ)
	/* SN: Bs->phi FF Zwicky et al. 1503.05534 (Bsphi_FormFactor_choice:5)*/
	||mystrcmp(a0A0_Bsphi_BSZ_LCSRonly)||mystrcmp(a1A0_Bsphi_BSZ_LCSRonly)||mystrcmp(a2A0_Bsphi_BSZ_LCSRonly)||mystrcmp(a0A1_Bsphi_BSZ_LCSRonly)||mystrcmp(a1A1_Bsphi_BSZ_LCSRonly)||mystrcmp(a2A1_Bsphi_BSZ_LCSRonly)||mystrcmp(a0A12_Bsphi_BSZ_LCSRonly)||mystrcmp(a1A12_Bsphi_BSZ_LCSRonly)||mystrcmp(a2A12_Bsphi_BSZ_LCSRonly)||mystrcmp(a0V_Bsphi_BSZ_LCSRonly)||mystrcmp(a1V_Bsphi_BSZ_LCSRonly)||mystrcmp(a2V_Bsphi_BSZ_LCSRonly)||mystrcmp(a0T1_Bsphi_BSZ_LCSRonly)||mystrcmp(a1T1_Bsphi_BSZ_LCSRonly)||mystrcmp(a2T1_Bsphi_BSZ_LCSRonly)||mystrcmp(a0T2_Bsphi_BSZ_LCSRonly)||mystrcmp(a1T2_Bsphi_BSZ_LCSRonly)||mystrcmp(a2T2_Bsphi_BSZ_LCSRonly)||mystrcmp(a0T23_Bsphi_BSZ_LCSRonly)||mystrcmp(a1T23_Bsphi_BSZ_LCSRonly)||mystrcmp(a2T23_Bsphi_BSZ_LCSRonly)
	/* SN: Bs->phi FF Wingate et al. 1310.3722 + 1501.00367 (Bsphi_FormFactor_choice:6)*/
	||mystrcmp(a0A0_Bsphi_HLMW)||mystrcmp(a1A0_Bsphi_HLMW)||mystrcmp(a0A1_Bsphi_HLMW)||mystrcmp(a1A1_Bsphi_HLMW)||mystrcmp(a0A12_Bsphi_HLMW)||mystrcmp(a1A12_Bsphi_HLMW)||mystrcmp(a0V_Bsphi_HLMW)||mystrcmp(a1V_Bsphi_HLMW)||mystrcmp(a0T1_Bsphi_HLMW)||mystrcmp(a1T1_Bsphi_HLMW)||mystrcmp(a0T2_Bsphi_HLMW)||mystrcmp(a1T2_Bsphi_HLMW)||mystrcmp(a0T23_Bsphi_HLMW)||mystrcmp(a1T23_Bsphi_HLMW)||mystrcmp(systErr_Bsphi_HLMW)
	);

	/* Uncertainties considered in the low q2 region */
	int test_lowBsphi=(test_Bsphi||mystrcmp(Bstophilow_ALperp_err_noq2)||mystrcmp(Bstophilow_ARperp_err_noq2)||mystrcmp(Bstophilow_ALpar_err_noq2)||mystrcmp(Bstophilow_ARpar_err_noq2)||mystrcmp(Bstophilow_AL0_err_noq2)||mystrcmp(Bstophilow_AR0_err_noq2)||mystrcmp(Bstophilow_At_err_noq2)||mystrcmp(Bstophilow_AS_err_noq2)||mystrcmp(Bstophilow_ALperp_err_q2)||mystrcmp(Bstophilow_ARperp_err_q2)||mystrcmp(Bstophilow_ALpar_err_q2)||mystrcmp(Bstophilow_ARpar_err_q2)||mystrcmp(Bstophilow_AL0_err_q2)||mystrcmp(Bstophilow_AR0_err_q2)||mystrcmp(Bstophilow_At_err_q2)||mystrcmp(Bstophilow_AS_err_q2));

	/* Uncertainties considered in the intermediate q2 region */
	int test_intermBsphi=(test_Bsphi
	//
	||mystrcmp(Bstophilow_ALperp_err_noq2)||mystrcmp(Bstophilow_ARperp_err_noq2)||mystrcmp(Bstophilow_ALpar_err_noq2)||mystrcmp(Bstophilow_ARpar_err_noq2)||mystrcmp(Bstophilow_AL0_err_noq2)||mystrcmp(Bstophilow_AR0_err_noq2)||mystrcmp(Bstophilow_At_err_noq2)||mystrcmp(Bstophilow_AS_err_noq2)||mystrcmp(Bstophilow_ALperp_err_q2)||mystrcmp(Bstophilow_ARperp_err_q2)||mystrcmp(Bstophilow_ALpar_err_q2)||mystrcmp(Bstophilow_ARpar_err_q2)||mystrcmp(Bstophilow_AL0_err_q2)||mystrcmp(Bstophilow_AR0_err_q2)||mystrcmp(Bstophilow_At_err_q2)||mystrcmp(Bstophilow_AS_err_q2)
	//
	||mystrcmp(Bstophihigh_ALperp_err)||mystrcmp(Bstophihigh_ARperp_err)||mystrcmp(Bstophihigh_ALpar_err)||mystrcmp(Bstophihigh_ARpar_err)||mystrcmp(Bstophihigh_AL0_err)||mystrcmp(Bstophihigh_AR0_err)||mystrcmp(Bstophihigh_At_err)||mystrcmp(Bstophihigh_AS_err));

	/* Uncertainties considered in the high q2 region */
	int test_highBsphi=(test_Bsphi
	||mystrcmp(Bstophihigh_ALperp_err)||mystrcmp(Bstophihigh_ARperp_err)||mystrcmp(Bstophihigh_ALpar_err)||mystrcmp(Bstophihigh_ARpar_err)||mystrcmp(Bstophihigh_AL0_err)||mystrcmp(Bstophihigh_AR0_err)||mystrcmp(Bstophihigh_At_err)||mystrcmp(Bstophihigh_AS_err));


    int test_LbLll=(mystrcmp(life_Lb)||mystrcmp(alphaL_LbLll)||mystrcmp(a0_HO_fplus_LbLll)||mystrcmp(a1_HO_fplus_LbLll)||mystrcmp(a2_HO_fplus_LbLL)||mystrcmp(a0_HO_fperp_LbLll)||mystrcmp(a1_HO_fperp_LbLll)||mystrcmp(a2_HO_fperp_LbLll)||mystrcmp(a0_HO_gpp_LbLll)||mystrcmp(a1_HO_gplus_LbLll)||mystrcmp(a2_HO_gplus_LbLll)||mystrcmp(a1_HO_gperp_LbLll)||mystrcmp(a2_HO_gperp_LbLll)||mystrcmp(a0_HO_hplus_LbLll)||mystrcmp(a1_HO_hplus_LbLll)||mystrcmp(a2_HO_hplus_LbLll)||mystrcmp(a0_HO_hperp_LbLll)||mystrcmp(a1_HO_hperp_LbLll)||mystrcmp(a2_HO_hperp_LbLll)||mystrcmp(a0_HO_htildepp_LbLll)||mystrcmp(a1_HO_htildeplus_LbLll)||mystrcmp(a2_HO_htildeplus_LbLll)||mystrcmp(a1_HO_htildeperp_LbLll)||mystrcmp(a2_HO_htildeperp_LbLll));
    
    int test_Kaons=(mystrcmp(deltaPcu_Kppipnunu)||mystrcmp(err_Pc_Xlambda_Kppipnunu)||mystrcmp(f_K)
    ||mystrcmp(BR_KLgammagamma_exp)||mystrcmp(Aterm_mu_KLmumu)||mystrcmp(chi_gg_Mrho)||mystrcmp(BR_KSgammagamma_exp)||mystrcmp(Iterm_mu_KSmumu)
	||mystrcmp(log_muK_1GeV)||mystrcmp(KLpill_Ce_Cdir)||mystrcmp(KLpill_Cmu_Cdir)||mystrcmp(KLpill_Ce_Cint)||mystrcmp(KLpill_Cmu_Cint)||mystrcmp(KLpill_Ce_Cmix)||mystrcmp(KLpill_Cmu_Cmix)||mystrcmp(KLpill_Cmu_CPC)||mystrcmp(KLpill_abs_aS)||mystrcmp(KLpill_Ce_Sgg)||mystrcmp(KLpill_Cmu_Sgg)||mystrcmp(KLpill_AFB_a1pi_gg_7V)||mystrcmp(KLpill_AFB_a1pi_gg_aS)
    );

#undef mystrcmp

	for(ie=0;ie<nobs;ie++)
	{
		if(test_general) values[ie]=compute_nameobs(&obs[ie],ke,&param);
		
		else if(!strcmp(obs[ie].type,"AI")&&!strcmp(obs[ie].decay,"BKstargamma")&&test_delta0Bkstargamma) values[ie]=compute_nameobs(&obs[ie],ke,&param);
		
		else if(!strcmp(obs[ie].type,"BR")&&!strcmp(obs[ie].decay,"B0Kstar0gamma")&&test_BKstargamma) values[ie]=compute_nameobs(&obs[ie],ke,&param);

		else if(!strcmp(obs[ie].type,"BR")&&!strcmp(obs[ie].decay,"BXsgamma")&&test_bsgamma) values[ie]=compute_nameobs(&obs[ie],ke,&param);

		else if((!strcmp(obs[ie].decay,"Bsmumu")||!strcmp(obs[ie].decay,"Bdmumu")||!strcmp(obs[ie].decay,"Bsee")||!strcmp(obs[ie].decay,"Bdee")||!strcmp(obs[ie].decay,"Bstautau")||!strcmp(obs[ie].decay,"Bdtautau")||!strcmp(obs[ie].decay,"Bsll")||!strcmp(obs[ie].decay,"Bdll"))&&test_Bsmumu) values[ie]=compute_nameobs(&obs[ie],ke,&param);

		else if((!strcmp(obs[ie].decay,"BXsmumu")||!strcmp(obs[ie].decay,"BXsee")||!strcmp(obs[ie].decay,"BXstautau")||!strcmp(obs[ie].decay,"BXsll"))&&test_BXsll) values[ie]=compute_nameobs(&obs[ie],ke,&param);

		/* SN: Generalizing the low- and high-q2 limits - 28/04/205 */
		// else if((!strncmp(obs[ie].decay,"BKstar",6)||!strncmp(obs[ie].decay,"B0Kstar0",8))&&
		// ((obs[ie].high<=6.&&test_lowBKstar)||(obs[ie].low>=15.&&test_highBKstar)||((obs[ie].high>6.||obs[ie].low<15.)&&test_intermBKstar))) values[ie]=compute_nameobs(&obs[ie],ke,&param);
  //
		// else if((!strncmp(obs[ie].decay,"BK",2)||!strncmp(obs[ie].decay,"B0K0",4))&&
		// ((obs[ie].high<=6.&&test_lowBK)||(obs[ie].low>=15.&&test_highBK)||((obs[ie].high>6.||obs[ie].low<15.)&&test_intermBK))) values[ie]=compute_nameobs(&obs[ie],ke,&param);
  //
		// else if(!strncmp(obs[ie].decay,"Bsphi",5)&&
		// ((obs[ie].high<=6.&&test_lowBsphi)||(obs[ie].low>=15.&&test_highBsphi)||((obs[ie].high>6.||obs[ie].low<15.)&&test_intermBsphi))) values[ie]=compute_nameobs(&obs[ie],ke,&param);
		else if((!strncmp(obs[ie].decay,"BKstar",6)||!strncmp(obs[ie].decay,"B0Kstar0",8))&&
		((obs[ie].high<=LOW_q2&&test_lowBKstar)||(obs[ie].low>=HIGH_q2&&test_highBKstar)||((obs[ie].high>LOW_q2||obs[ie].low<HIGH_q2)&&test_intermBKstar))) values[ie]=compute_nameobs(&obs[ie],ke,&param);
		/* SN: for parenthesis with test_intermBKstar, shouldn't it be && instead of && and similarly for Bsphill and BKll */

		else if((!strncmp(obs[ie].decay,"BK",2)||!strncmp(obs[ie].decay,"B0K0",4))&&
		((obs[ie].high<=LOW_q2&&test_lowBK)||(obs[ie].low>=HIGH_q2&&test_highBK)||((obs[ie].high>LOW_q2||obs[ie].low<HIGH_q2)&&test_intermBK))) values[ie]=compute_nameobs(&obs[ie],ke,&param);

		else if(!strncmp(obs[ie].decay,"Bsphi",5)&&
		((obs[ie].high<=LOW_q2&&test_lowBsphi)||(obs[ie].low>=HIGH_q2&&test_highBsphi)||((obs[ie].high>LOW_q2||obs[ie].low<HIGH_q2)&&test_intermBsphi))) values[ie]=compute_nameobs(&obs[ie],ke,&param);


		else if(!(strncmp(obs[ie].decay,"LambdabLambdamumu",17)||!strncmp(obs[ie].decay,"LambdabLambdaee",15)||!strncmp(obs[ie].decay,"LambdabLambdatautau",19)||!strncmp(obs[ie].decay,"LambdabLambdall",15))&&test_LbLll) values[ie]=compute_nameobs(&obs[ie],ke,&param);
		
		else if((!strcmp(obs[ie].decay,"KLpi0nunu")||!strcmp(obs[ie].decay,"Kpinunu")||!strcmp(obs[ie].decay,"KSmumu")||!strcmp(obs[ie].decay,"KLmumu")
		||!strcmp(obs[ie].decay,"KLee")
		||!strcmp(obs[ie].decay,"KLpi0mumu")
		||!strcmp(obs[ie].decay,"KLpi0ee")
		)&&test_Kaons) values[ie]=compute_nameobs(&obs[ie],ke,&param);
		
		else values[ie]=values_ref[ie];
	}
	
	return;
}

/*---------------------------------------------------------------------*/

void write_correlation_nuisance(char name[], double **corr, char nameparam[][50], int n)
/* Writes the correlation matrix "corr" into file "name" (diagonal elements are not written) */ 
{
	FILE *output;
	
	output=fopen(name,"w");
	
	int ie,je;
			
	for(ie=1;ie<=n+1;ie++) for(je=ie+1;je<=n+1;je++) if(corr[ie-1][je-1]!=0.) fprintf(output,"%-30s\t%-30s\t%g\n",nameparam[ie],nameparam[je],corr[ie-1][je-1]);

	fclose(output);
	
	return;
}

/*---------------------------------------------------------------------*/

void convert_correlation(nuiscorr corrnuis[], int na, double **corr, char nameparam[][50], int n)
/* Converts correlation array "corrnuis" of size na into matrix "corr" of size n*n */ 
{
	int ie,je,ke;
	
	for(ie=0;ie<=n;ie++) for(je=0;je<=n;je++) corr[ie][je]=(ie==je);
	
	for(ke=0;ke<na;ke++)
	{
		for(ie=0;ie<=n;ie++) for(je=0;je<=n;je++) if(!strcmp(corrnuis[ke].obs1,nameparam[ie])&&!strcmp(corrnuis[ke].obs2,nameparam[je])){corr[ie][je]=corrnuis[ke].value;ie=n+1;je=n+1;}
	}
	
	for(ie=0;ie<=n;ie++) for(je=0;je<=n;je++)
	{
		if(ie==je) corr[ie][je]=1.;
		else if(corr[ie][je]!=0.&&corr[je][ie]==0.) corr[je][ie]=corr[ie][je];
	}
	
	return;
}

/*---------------------------------------------------------------------*/

void read_correlation(char name[], double **corr, char nameparam[][50], int n)
/* Reads file "name" into the correlation matrix "corr" of size n*n */ 
{
	FILE *input;

	int ie,je;
	
	if(!test_file(name)) {printf("%s does not exist!\n",name); exit(1);}
	
	
	struct nuiscorr corrnuis;
	
	for(ie=0;ie<=n;ie++) for(je=0;je<=n;je++) corr[ie][je]=(ie==je);
	
	input=fopen(name,"r");

	while(EOF!=fscanf(input,"%s",corrnuis.obs1))
	{
		if(!strncasecmp("#",corrnuis.obs1,1)) while ((EOF!=fscanf(input,"%c",corrnuis.obs1))&&(strncasecmp("\n",corrnuis.obs1,1)));
		else
		{
			fscanf(input,"%s",corrnuis.obs2);
			fscanf(input,"%lf",&corrnuis.value);

			for(ie=0;ie<=n;ie++) for(je=0;je<=n;je++) if(!strcmp(corrnuis.obs1,nameparam[ie])&&!strcmp(corrnuis.obs2,nameparam[je])){corr[ie][je]=corrnuis.value;ie=n+1;je=n+1;}
		}
	}

	fclose(input);
	
	for(ie=0;ie<=n;ie++) for(je=0;je<=n;je++)
	{
		if(ie==je) corr[ie][je]=1.;
		else if(corr[ie][je]!=0.&&corr[je][ie]==0.) corr[je][ie]=corr[ie][je];
	}
	
	return;
}

/*---------------------------------------------------------------------*/

int read_experimental_covariance(char name_val_exp[], char name_corr_exp[], char namesin[][50], int nbobsin, double* central_exp, double* errors_exp, double **correlations)
/* Reads the experimental values, errors and correlations from files name_val_exp and name_corr_exp, if they match the names given in the string array namesin, into central_exp, errors_exp and correlations */ 
{	
	int ie,je;
	
	int translate[nbobsin];
		
	char names[NBOBSMAX][50];
	double central[NBOBSMAX],errors[NBOBSMAX],stat,syst;
	
	char dummy[50];

	FILE *input;
	
	int nbobs=0;
	
	input=fopen(name_val_exp,"r");
	
	while(EOF!=fscanf(input,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(input,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		else
		{
			sprintf(names[nbobs],"%s",dummy);
			fscanf(input,"%lf",&central[nbobs]);
			fscanf(input,"%lf",&stat);
			fscanf(input,"%lf",&syst);
			errors[nbobs]=sqrt(stat*stat+syst*syst);
			nbobs++;
		}
	}
	
	fclose(input);

	if(nbobs>NBOBSMAX){printf("Memory problem: please increase NBOBSMAX to at least %d in src/include.h and recompile.",nbobs);exit(1);}

	for(ie=0;ie<nbobsin;ie++) translate[ie]=-1;

	for(ie=0;ie<nbobsin;ie++)
	{
		je=0;
		while(je<nbobs)
		{
			if(names[je]!=NULL)
			{
				if(!strcmp(namesin[ie],names[je]))
				{
					translate[ie]=je;
					je=nbobs;
				}
				else je++;
			} 
			else je++;
		}
	}

	for(ie=0;ie<nbobsin;ie++)
	{
		if(translate[ie]>-1&&errors[translate[ie]]!=0.)
		{
			central_exp[ie]=central[translate[ie]];
			errors_exp[ie]=errors[translate[ie]];
		}
		else 
		{
			printf("Warning! %s not present (or experimental error set to 0) in %s. Error set to 10^30...\n",name_val_exp,namesin[ie]);  
			central_exp[ie]=0.;
			errors_exp[ie]=1.e30;
		}
	}
	
	double **correlations_exp=(double **) malloc(NBOBSMAX*sizeof(double *));
	for(ie=0;ie<NBOBSMAX;ie++) correlations_exp[ie]=(double *) malloc(NBOBSMAX*sizeof(double));		
	
	read_correlation(name_corr_exp,correlations_exp,names,NBOBSMAX-1);
	
	for(ie=0;ie<nbobsin;ie++) for(je=0;je<nbobsin;je++) correlations[ie][je]=(ie==je);

	for(ie=0;ie<nbobsin;ie++)
	{
		for(je=0;je<nbobsin;je++)
		{
			if((translate[ie]==-1||translate[je]==-1))
			{
				correlations[ie][je]=(ie==je);
			}
			else
			{
				correlations[ie][je]=correlations_exp[translate[ie]][translate[je]];
			}
		}
	}
	
	for (ie = 0; ie < NBOBSMAX; ie++) free(correlations_exp[ie]);
	free(correlations_exp);
	
	return nbobs;
}

/*---------------------------------------------------------------------*/

void get_th_covariance_nuisance(double ***covariance_th, char names[][50], int *nbobs, struct parameters* param, struct nuisance* nuisparam, double **nuiscorr)
/* Computes the theory covariance matrix for the "nbobs" observables with names given in "names", using the nuisance parameters contained in the structure nuisparam and array nuiscorr */
{		
	double values_ref[*nbobs];
		
	int n=NNUIS; /* number of nuisance parameters */
	
	obsname obs[*nbobs];

	make_obslist(names,obs,nbobs);
	
	char namenuisance[NNUIS+1][50];
	
	observables(0,obs,*nbobs,values_ref,values_ref,nuisparam,namenuisance,param);
		
	int ie,je,ke,le;
	
	double **values_mod;
	values_mod=(double **) malloc((n+1)*sizeof(double *));
	for(ie=0;ie<=n;ie++) values_mod[ie]=(double *) malloc(*nbobs*sizeof(double));	

	for(ie=0;ie<*nbobs;ie++) values_mod[0][ie]=values_ref[ie];		
	
	for(ke=1;ke<=n;ke++)
	{			
		double values[*nbobs];
					
		observables(ke,obs,*nbobs,values,values_ref,nuisparam,namenuisance,param);
		
		for(ie=0;ie<*nbobs;ie++) values_mod[ke][ie]=values[ie];
	}	
		
	*covariance_th=NULL;
	*covariance_th=(double **) malloc(*nbobs*sizeof(double *));
	for(ie=0;ie<*nbobs;ie++) (*covariance_th)[ie]=(double *) malloc(*nbobs*sizeof(double));		

	for(ie=0;ie<*nbobs;ie++) for(je=0;je<*nbobs;je++)
	{
		(*covariance_th)[ie][je]=0.;
		for(ke=1;ke<=n;ke++) for(le=1;le<=n;le++) (*covariance_th)[ie][je]+=nuiscorr[ke][le]*(values_mod[ke][ie]-values_ref[ie])*(values_mod[le][je]-values_ref[je]);
	}

	for (ie = 0; ie < n+1; ie++) free(values_mod[ie]);
	free(values_mod);

	return;
}

/*---------------------------------------------------------------------*/

void get_exp_covariance(double ***covariance_exp, double **central_exp, char names[][50], int *nbobs, struct parameters* param)
{
	int ie,je;
	
	*covariance_exp=NULL;
	*covariance_exp=(double **) malloc(*nbobs*sizeof(double *));
	for(ie=0;ie<*nbobs;ie++) (*covariance_exp)[ie]=(double *) malloc(*nbobs*sizeof(double));		

	*central_exp=NULL;
	*central_exp=(double *) malloc(*nbobs*sizeof(double));		

	double *sigma_exp=NULL;	
	sigma_exp=(double *) malloc(*nbobs*sizeof(double));		
 
	read_experimental_covariance(param->exp_values,param->exp_corr,names,*nbobs,*central_exp,sigma_exp,*covariance_exp);
				
	for(ie=0;ie<*nbobs;ie++) for(je=0;je<*nbobs;je++) (*covariance_exp)[ie][je]*=sigma_exp[ie]*sigma_exp[je];

	return;
}

/*---------------------------------------------------------------------*/

void get_exp_values(double **central_exp, char names[][50], int *nbobs, struct parameters* param)
{
	int ie,je;
	
	*central_exp=NULL;
	*central_exp=(double *) malloc(*nbobs*sizeof(double));		

	read_experimental_values(param->exp_values,names,*nbobs,*central_exp);
				
	return;
}
/*---------------------------------------------------------------------*/

void get_covariance(double ***covariance_th, double ***covariance_exp, double **central_exp, char names[][50], int *nbobs, struct parameters* param)
/* Gives the experimental abd theoretical covariance matrices, computing first the theory covariance matrix for the "nbobs" observables with names given in "names", using the nuisance parameters contained in the file "param->nuisance_values", and reading the experimental values from the files param->exp_values_xxx and param->exp_corr_xxx */
{				
	get_th_covariance(covariance_th,names,nbobs,param);
	
	get_exp_covariance(covariance_exp,central_exp,names,nbobs,param);
	
	return;
}

/*---------------------------------------------------------------------*/

void get_predictions_nuisance(char names[][50], int *nbobs, double** predictions, struct parameters* param, struct nuisance* nuisparam)
/* Computes the predictions/reference values of observables in string array "names" */
{		
 	double values_ref[*nbobs];
	
	int ie;
	
	obsname obs[*nbobs];

	make_obslist(names,obs,nbobs);
	
	char namenuisance[NNUIS+1][50];

	observables(0,obs,*nbobs,values_ref,values_ref,nuisparam,namenuisance,param);
		
	*predictions=NULL;
	*predictions=(double *) malloc(*nbobs*sizeof(double));	
	
	for(ie=0;ie<*nbobs;ie++) (*predictions)[ie]=values_ref[ie];

	return;
}

/*---------------------------------------------------------------------*/

double get_chi2(double **inv_cov_tot, double *predictions, double *central_exp, int nbobs)
/* Computes chi2 using the inverse of the covariance matrix, experimental values and theoretical predictions */
{		
	int ie,je;
	double chi2=0.;
	
	for(ie=0;ie<nbobs;ie++) for(je=0;je<nbobs;je++) chi2+=(predictions[ie]-central_exp[ie])*inv_cov_tot[ie][je]*(predictions[je]-central_exp[je]);
		
	return chi2;
}

/*---------------------------------------------------------------------*/

double chi2(char names[][50], int nbobs, struct parameters* param)
/* Automatic chi2 calculation for the list of observables "names", using input from files */
{
	int ie;
	double *predictions;
	
	get_predictions(names,&nbobs,&predictions,param);
	
	double **covariance_th,**covariance_exp,*central_exp;
	
	get_covariance(&covariance_th,&covariance_exp,&central_exp,names,&nbobs,param);

#ifdef DEBUG
	for(ie=0;ie<nbobs;ie++) printf("\t%-25s:\tprediction %.3e\ttherr %.3e\t\t(central %.3e\texperr %.3e)\n",names[ie],predictions[ie],sqrt(covariance_th[ie][ie]),central_exp[ie],sqrt(covariance_exp[ie][ie]));

	int je;
	FILE *file;
	file=fopen("correlation.out","w");	
	for(ie=0;ie<nbobs;ie++) for(je=ie+1;je<nbobs;je++) fprintf(file,"%s-%s: %.4e\n",names[ie],names[je],covariance_th[ie][je]/sqrt(covariance_th[ie][ie])/sqrt(covariance_th[je][je]));
	fclose(file);
#endif

	double **covariance_tot;
	get_covtot(&covariance_th,&covariance_exp,&covariance_tot,nbobs);

	double **inv_cov_tot;	
	if(get_invcovtot(&covariance_tot,&inv_cov_tot,nbobs)<0) return -3.;
			
	double chi_2=get_chi2(inv_cov_tot,predictions,central_exp,nbobs);
	
	for (ie = 0; ie < nbobs; ie++) free(covariance_th[ie]);
	free(covariance_th);

	for (ie = 0; ie < nbobs; ie++) free(covariance_exp[ie]);
	free(covariance_exp);
	
	free(central_exp);

	for (ie = 0; ie < nbobs; ie++) free(covariance_tot[ie]);
	free(covariance_tot);
	
	for (ie = 0; ie < nbobs; ie++) free(inv_cov_tot[ie]);
	free(inv_cov_tot);
	
	return chi_2;	
}

/*---------------------------------------------------------------------*/

int get_invcovtot(double ***covariance_tot, double ***inv_cov_tot, int nbobs)
/* Computes the inverse covariance matrix inv_cov_tot from the covariance matrix covariance_tot */
{	
	int ie;
				
	*inv_cov_tot=NULL;
	*inv_cov_tot=(double **) malloc(nbobs*sizeof(double *));
	for(ie=0;ie<nbobs;ie++) (*inv_cov_tot)[ie]=(double *) malloc(nbobs*sizeof(double));		

	if(invert_matrix(nbobs,*covariance_tot,*inv_cov_tot)==0) return -1; /* inversion of the total covariance matrix */
	else return 1;
}

/*---------------------------------------------------------------------*/

void get_covtot(double ***covariance_th, double ***covariance_exp, double ***covariance_tot, int nbobs)
/* Sums the theoretical and experimental covariance matrices into covariance_tot */
{		
	int ie,je;
	
	*covariance_tot=NULL;
	*covariance_tot=(double **) malloc(nbobs*sizeof(double *));
	for(ie=0;ie<nbobs;ie++) (*covariance_tot)[ie]=(double *) malloc(nbobs*sizeof(double));		

	
	for(ie=0;ie<nbobs;ie++) for(je=0;je<nbobs;je++) (*covariance_tot)[ie][je]=(*covariance_exp)[ie][je]+(*covariance_th)[ie][je];	
		
	return;
}

/*---------------------------------------------------------------------*/

int reduce_covariance(double ***covariance_in, char* namesin[], int nbobsin, double ***covariance_out, char* namesout[], int nbobsout)
{		
	if(nbobsout>nbobsin) return 0;

	int ie,je;
	
	int translate[nbobsout];
	for(ie=0;ie<nbobsout;ie++) translate[ie]=-1;
	
	for(ie=0;ie<nbobsout;ie++)
	{
		je=0;
		while(je<nbobsin)
		{
			if(namesin[je]!=NULL)
			{
				if(!strcmp(namesout[ie],namesin[je]))
				{
					translate[ie]=je;
					je=nbobsin;
				}
				else je++;
			} 
			else je++;
		}
	}

	for(ie=0;ie<nbobsout;ie++) if(translate[ie]==-1) return 0;
	
	*covariance_out=NULL;
	*covariance_out=(double **) malloc(nbobsout*sizeof(double *));
	for(ie=0;ie<nbobsout;ie++) (*covariance_out)[ie]=(double *) malloc(nbobsout*sizeof(double));		

	for(ie=0;ie<nbobsout;ie++)  for(je=0;je<nbobsout;je++) (*covariance_out)[ie][je]=(*covariance_in)[translate[ie]][translate[je]];
	return 1;
}

/*---------------------------------------------------------------------*/

int reduce_values(double **values_in, char* namesin[], int nbobsin, double **values_out, char* namesout[], int nbobsout)
{		
	if(nbobsout>nbobsin) return 0;

	int ie,je;
	
	int translate[nbobsout];
	for(ie=0;ie<nbobsout;ie++) translate[ie]=-1;
	
	for(ie=0;ie<nbobsout;ie++)
	{
		je=0;
		while(je<nbobsin)
		{
			if(namesin[je]!=NULL)
			{
				if(!strcmp(namesout[ie],namesin[je]))
				{
					translate[ie]=je;
					je=nbobsin;
				}
				else je++;
			} 
			else je++;
		}
	}

	for(ie=0;ie<nbobsout;ie++) if(translate[ie]==-1) return 0;

	*values_out=NULL;
	*values_out=(double *) malloc(nbobsout*sizeof(double));		

	for(ie=0;ie<nbobsout;ie++)  (*values_out)[ie]=(*values_in)[translate[ie]];
	
	return 1;
}

/*---------------------------------------------------------------------*/

void read_covariance(char name[], double **cov, char nameparam[][50], int n)
/* Reads file "name" into the covariance matrix "cov" of size n*n */ 
{
	FILE *input;

	int ie,je;
	
	if(!test_file(name)) {printf("%s does not exist!\n",name); exit(1);}
			
	struct nuiscorr corrnuis;
		
	input=fopen(name,"r");
	
	while(EOF!=fscanf(input,"%s",corrnuis.obs1))
	{
		if(!strncasecmp("#",corrnuis.obs1,1)) while ((EOF!=fscanf(input,"%c",corrnuis.obs1))&&(strncasecmp("\n",corrnuis.obs1,1)));
		else
		{
			fscanf(input,"%s",corrnuis.obs2);
			fscanf(input,"%lf",&corrnuis.value);
			
			for(ie=0;ie<n;ie++) for(je=0;je<n;je++) if(!strcmp(corrnuis.obs1,nameparam[ie])&&!strcmp(corrnuis.obs2,nameparam[je])){cov[ie][je]=cov[je][ie]=corrnuis.value;ie=n;je=n;}			
		}
	}

	fclose(input);
	
	return;
}

/*---------------------------------------------------------------------*/

void write_covariance(char name[], double **cov, char nameparam[][50], int n)
/* Writes the covariance matrix "corr" into file "name"  (diagonal elements are written) */ 
{
	FILE *output;
	
	output=fopen(name,"w");
	
	int ie,je;
			
	for(ie=0;ie<n;ie++) for(je=ie;je<n;je++) if(cov[ie][je]!=0.) fprintf(output,"%-30s\t%-30s\t%g\n",nameparam[ie],nameparam[je],cov[ie][je]);

	fclose(output);
	
	return;
}

/*---------------------------------------------------------------------*/

int read_experimental_values(char name_val_exp[], char namesin[][50], int nbobsin, double* central_exp)
/* Reads the experimental values from files name_val_exp, if they match the names given in the string array namesin, into central_exp */ 
{	
	int ie,je;
	
	int translate[nbobsin];
		
	char names[NBOBSMAX][100];
	double central[NBOBSMAX],errors[NBOBSMAX],stat,syst;
	
	char dummy[100];

	FILE *input;
	
	int nbobs=0;
	
	input=fopen(name_val_exp,"r");
	
	while(EOF!=fscanf(input,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(input,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		else
		{
			sprintf(names[nbobs],"%s",dummy);
			fscanf(input,"%lf",&central[nbobs]);
			fscanf(input,"%lf",&stat);
			fscanf(input,"%lf",&syst);
			errors[nbobs]=sqrt(stat*stat+syst*syst);
			nbobs++;
		}
	}
	
	fclose(input);

	if(nbobs>NBOBSMAX){printf("Memory problem: please increase NBOBSMAX to at least %d in src/include.h and recompile.",nbobs);exit(1);}

	for(ie=0;ie<nbobsin;ie++) translate[ie]=-1;

	for(ie=0;ie<nbobsin;ie++)
	{
		je=0;
		while(je<nbobs)
		{
			if(names[je]!=NULL)
			{
				if(!strcmp(namesin[ie],names[je]))
				{
					translate[ie]=je;
					je=nbobs;
				}
				else je++;
			} 
			else je++;
		}
	}

	for(ie=0;ie<nbobsin;ie++)
	{
		if(translate[ie]>-1&&errors[translate[ie]]!=0.)
		{
			central_exp[ie]=central[translate[ie]];
		}
		else 
		{
			printf("Warning! %s not present (or experimental error set to 0) in %s. Value set to zero.\n",name_val_exp,namesin[ie]);  
			central_exp[ie]=0.;
		}
	}

	return nbobs;
}
