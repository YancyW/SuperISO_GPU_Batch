#include "include.h"

/**************************************************************/
/*                KL->pi0nunu  and K+->pi+nunu                */
/**************************************************************/

/*----------------------------------------------------------------------------*/

double X0(double xt)
{
  return xt/8. * ((xt + 2.)/(xt - 1.) + (3.*xt - 6.)/pow(xt - 1.,2.) * log(xt));
}

/*----------------------------------------------------------------------------*/

double X1(double xt, double mu, double mass_W)
{
	return -(29.*xt - xt*xt - 4.*pow(xt,3.))/(3.*pow(1.-xt,2.)) - (xt + 9.*xt*xt - pow(xt,3.) - pow(xt,4.))/(pow(1.-xt,3.)) * log(xt)
	+ (8.*xt + 4.*xt*xt + xt*xt*xt - pow(xt,4.))/(2*pow(1.-xt,3)) * log(xt)*log(xt) - (4.*xt - pow(xt,3.))/(pow(1.-xt,2.)) * Li2(1.-xt)
	+ 8.*xt*( (8. - 9.*xt + pow(xt,3.) +6.*log(xt))/(8.*pow(-1+xt,3.)) )*log(mu*mu/(mass_W*mass_W));
}

/*----------------------------------------------------------------------------*/

double Pc_Xlambda(double lambda, double mc, double alphas_Z)
// Returns the value of the charm sector contribution see equation (50) of 0805.4119v1 /*err: 1.15<mc<1.45, 0.1150<alphas<0.1230, 1.0<mu_c<3.0, 2.5<mu_b<10.0 */
{
  double kappa10 = 1.6624;
  double kappa01 = -2.3537;
  double kappa11 = -1.5862;
  double kappa20 = 1.5036;
  double kappa02 = -4.3477;
  double Lmc = log(mc/1.3);
  double Lalphas = log(alphas_Z/0.1176);

  return pow(0.2255/lambda,4.)*0.38049*pow(mc/1.3,0.5081)*pow(alphas_Z/0.1176,1.0192)*(1. + (Lalphas*kappa01+Lalphas*Lalphas*kappa02+Lmc*kappa10 + Lmc*Lalphas*kappa11 + Lmc*Lmc*kappa20));
}

double Pc_Xlambda_err_mc_alphas(double mc, double alphas_Z)
// Returns the value of the error of the charm sector contribution but does not include the error from lambda -- see equation (50) of 0805.4119v1  xxxx negligible error??
{
  double epsilon10 = -0.3537;
  double epsilon01 = 0.6003;
  double epsilon11 = -4.7652;
  double epsilon20 = 1.0253;
  double epsilon02 = 0.8866;
  double Lmc = log(mc/1.3);
  double Lalphas = log(alphas_Z/0.1176);

  return 0.008707*pow(mc/1.3,0.5276)*pow(alphas_Z/0.1176,1.8970)*(1. + (Lalphas*epsilon01+Lalphas*Lalphas*epsilon02+Lmc*epsilon10 + Lmc*Lalphas*epsilon11 + Lmc*Lmc*epsilon20));
}


/*----------------------------------------------------------------------------*/

double BRKLpinunu(double complex CLe, double complex CRe, double complex CLmu, double complex CRmu, double complex CLtau, double complex CRtau, struct parameters* param)
{
	if(CLe==0.) CLe+=param->KCLe;
	if(CRe==0.) CRe+=param->KCRe;
	if(CLmu==0.) CLmu+=param->KCLmu;
	if(CRmu==0.) CRmu+=param->KCRmu;
	if(CLtau==0.) CLtau+=param->KCLtau;
	if(CRtau==0.) CRtau+=param->KCRtau;

	double complex Vts,Vtd,lambdat;
	Vts = param->Vts;
	Vtd = param->Vtd;
	lambdat = Vtd*conj(Vts);
	double lambda = cabs(param->Vus);
	
	double mtmt = param->mtmt;

	double xt = pow(mtmt/param->mass_W,2.);
	double alphas_mut = alphas_running(mtmt,mtmt,param->mass_b_pole,param);

	double XtSM = X0(xt) + alphas_mut/(4.*pi) * X1(xt,mtmt,param->mass_W);
	
	double sw2=pow(sin(atan(param->gp/param->g2)),2.);


	double kappa_L;/*(Table 1 in 1009.0947 or Table 2 in 0705.2025 for details also see eq.12 in 2105.02868)*/
	kappa_L = 2.231*pow(10.,-10.)*pow(lambda/0.225,8.)*pow(127.9/param->inv_alpha_em,2.)*pow(0.231/sw2,2); /* err: +-0.013*pow(10.,-10.)*/


	double complex CNPe, CNPmu, CNPtau;
	CNPe = CLe + CRe;
	CNPmu = CLmu + CRmu;
	CNPtau = CLtau + CRtau;


	
// For SM we have XtSM = -sw2*CSM, and the NP contributions are are added as C = CSM + CNP  
	double complex XNPe, XNPmu, XNPtau;
	XNPe = -sw2*CNPe;
	XNPmu = -sw2*CNPmu;
	XNPtau = -sw2*CNPtau;
	
	return kappa_L/pow(lambda,10.)/3.*(pow(cimag(lambdat*(XtSM + XNPe)),2.) + pow(cimag(lambdat*(XtSM + XNPmu)),2.) + pow(cimag(lambdat*(XtSM + XNPtau)),2.));

}

/*----------------------------------------------------------------------------*/

double BRKppipnunu(double complex CLe, double complex CRe, double complex CLmu, double complex CRmu, double complex CLtau, double complex CRtau, struct parameters* param)
{
	if(CLe==0.) CLe+=param->KCLe;
	if(CRe==0.) CRe+=param->KCRe;
	if(CLmu==0.) CLmu+=param->KCLmu;
	if(CRmu==0.) CRmu+=param->KCRmu;
	if(CLtau==0.) CLtau+=param->KCLtau;
	if(CRtau==0.) CRtau+=param->KCRtau;
	
	double deltaPcu = param->deltaPcu_Kppipnunu; // 0.04+-0.02   see 0805.4119v1 xxxxxxxx
		
	double complex Vts,Vtd,Vcs,Vcd,lambdat,lambdac;
	Vts = param->Vts;
	Vtd = param->Vtd;
	Vcs = param->Vcs;
	Vcd = param->Vcd;
	lambdat = Vtd*conj(Vts);
	lambdac = Vcd*conj(Vcs);
	double lambda = cabs(param->Vus);
	
	double PcX = Pc_Xlambda(lambda, param->mass_c, param->alphas_MZ) + param->err_Pc_Xlambda_Kppipnunu*Pc_Xlambda_err_mc_alphas(param->mass_c, param->alphas_MZ);
	

	double mtmt = param->mtmt;

	double xt = pow(mtmt/param->mass_W,2.);
	double alphas_mut = alphas_running(mtmt,mtmt,param->mass_b_pole,param);

	double XtSM = X0(xt) + alphas_mut/(4.*pi) * X1(xt,mtmt,param->mass_W);
	double XcSM = pow(lambda,4.)*(PcX + deltaPcu);// (see e.g. eq.145 1609.04783)
	double XcSMemu = 10.05e-4;/* eq. 147 in 1609.04783 */
	double XcSMtau = 6.64e-4;/* eq. 147 in 1609.04783 */
	
	double sw2=pow(sin(atan(param->gp/param->g2)),2.);

	double DeltaEM = -0.003;
	double kappap = 0.5173*pow(10.,-10.)*pow(lambda/0.225,8.)*pow(127.9/param->inv_alpha_em,2.)*pow(0.231/sw2,2);/* err: 0.0025*pow(10.,-10.); (Table 1 in 1009.0947 or Table 2 in 0705.2025 for details also see eq.12 in 2105.02868)*/

	double complex CNPe, CNPmu, CNPtau;
	CNPe = CLe + CRe;
	CNPmu = CLmu + CRmu;
	CNPtau = CLtau + CRtau;

// For SM we have XtSM = -sw2*CSM, and the NP contributions are are added as C = CSM + CNP  
	double complex XNPe, XNPmu, XNPtau;
	XNPe = -sw2*CNPe;
	XNPmu = -sw2*CNPmu;
	XNPtau = -sw2*CNPtau;
	
/* NP contributions coded separately since we only have NLO results for neutrino flavour separated Xc; used for the NPxSN interference term.*/
	return kappap*(1. + DeltaEM)/pow(lambda,10.)*( pow(cimag(lambdat*XtSM),2.) + pow(creal(lambdac*XcSM + lambdat*XtSM),2.) )
	       + kappap*(1. + DeltaEM)/pow(lambda,10.)/3.*( 
		  2.*cimag(lambdat*XtSM)*cimag(lambdat*(XNPe + XNPmu + XNPtau))
		+ pow(cimag(lambdat*XNPe),2.) + pow(cimag(lambdat*XNPmu),2.) + pow(cimag(lambdat*XNPtau),2.) 
		+ 2.*creal(lambdat*XtSM)*creal(lambdat*(XNPe + XNPmu + XNPtau))
		+ pow(creal(lambdat*XNPe),2.) + pow(creal(lambdat*XNPmu),2.) + pow(creal(lambdat*XNPtau),2.)
		+ 2.*( creal(lambdac*XcSMemu)*creal(lambdat*XNPe) + creal(lambdac*XcSMemu)*creal(lambdat*XNPmu) + creal(lambdac*XcSMtau)*creal(lambdat*XNPtau) )   
		 );
}


/*----------------------------------------------------------------------------*/

/**************************************************************/
/*                      K_{L,S}->mumu                         */
/**************************************************************/

/* hep-ph/9901288 */
double Y0(double xt)
{
	return xt/8. * ((4.-xt)/(1.-xt) + (3.*xt)/pow(1.-xt,2.) * log(xt));
}

/**************************************************************/

/* hep-ph/9901288 */
double Y1(double xt, double mu, double mass_W)
{
	return (10.*xt + 10*xt*xt + 4.*pow(xt,3.))/(3.*pow(1.-xt,2.)) - (2.*xt - 8.*xt*xt - pow(xt,3.) - pow(xt,4.))/(pow(1.-xt,3.)) * log(xt)
	+ (2.*xt - 14.*xt*xt + xt*xt*xt - pow(xt,4.))/(2*pow(1.-xt,3)) * log(xt)*log(xt) + (2.*xt + pow(xt,3.))/(pow(1.-xt,2.)) * Li2(1.-xt)
	+ 8.*xt*log(mu*mu/(mass_W*mass_W))*( (-4. + 3.*xt + xt*xt*xt -6.*xt*log(xt))/(8.*pow(xt-1.,3.)) );
}

/**************************************************************/

/* hep-ph/0605203, eq.6 */ /*err: 1.15<mc<1.45, 0.1150<alphas<0.1230, 1.0<mu_c<3.0, 2.5<mu_b<10.0 */
double Pc_Yscales(double lambda, double mc, double alphas_Z, double muc, double mub)
{
	double kappa1000 = -0.5373;
	double kappa0100 = -6.0472;
	double kappa0010 = -0.0956;
	double kappa0001 = 0.0114;
//
	double kappa1100 = 3.9957;
	double kappa1010 = 0.3604;
	double kappa0110 = 0.0516;
	double kappa0101 = -0.0658;
//
	double kappa2000 = -0.1767;
	double kappa0200 = 16.4465;
	double kappa0020 = -0.1294;
	double kappa0030 = 0.0725;

	double Lmc = log(mc/1.3);
	double Las = log(alphas_Z/0.1187);
	double Lmuc = log(muc/1.5);
	double Lmub = log(mub/5.);
	
	return pow(0.225/lambda,4.)*(0.1198*pow(mc/1.3,2.3595)*pow(alphas_Z/0.1187,6.6055)*(1. 
	+ (kappa1000*Lmc + kappa0100*Las + kappa0010*Lmuc + kappa0001*Lmub 
	+ kappa1100*Lmc*Las + kappa1010*Lmc*Lmuc + kappa0110*Las*Lmuc + kappa0101*Las*Lmub
	+ kappa2000*Lmc*Lmc + kappa0200*Las*Las + kappa0020*Lmuc*Lmuc + kappa0030*Lmuc*Lmuc*Lmuc )));
}


/* hep-ph/0605203, eq.8 */
double Pc_Ylambda(double lambda)
{
	return 0.115*pow(0.225/lambda,4.); /* err: +-0.008_theory +-0.008_mc +-0.001_alphas */
}

/*----------------------------------------------------------------------*/
// One-loop amplitudes (Fkaon, Gkaon, Hkaon) for rare and radiative K decays; appendix of 1107.6001 (all not needed, but added for completeness)
double complex Fkaon(double z)
{
	if(z <= 4.) return 1. - 4./z*pow(asin(sqrt(z)/2.),2.);
	
	else return 1.+1./z*cpow(log( (1.-sqrt(1.-4./z))/(1.+sqrt(1.-4./z)) ) +I*pi,2.);
}
/*----------------------------------------------------------------------*/

double complex Gkaon(double z)
{
	if(z <= 4.) return sqrt(4./z-1.)*asin(sqrt(z)/2.);
	
	else return 1./2.*sqrt(1-4./z)*(log( (1.+sqrt(1.-4./z))/(1.-sqrt(1.-4./z)) ) -I*pi);
}
/*----------------------------------------------------------------------*/

double complex Hkaon(double z, double rpi)
{
	if(z == 0.) return -1./2*Fkaon(1./rpi/rpi);
		
	else return 1./(2.*pow(1.-z,2.))*( z*Fkaon(z/rpi/rpi) -Fkaon(1./rpi/rpi) -2.*z*(Gkaon(z/rpi/rpi) -Gkaon(1/rpi/rpi)) );
	
}

/*----------------------------------------------------------------------*/



/**************************************************************/

/*arXiv: 1711.11030*/
double BRKLll(int gen, int SL, double complex C10, double complex C10p, double complex CQ1, double complex CQ1p, double complex CQ2, double complex CQ2p, struct parameters* param)
{

	if(C10==0.) C10+=param->KC10mu;
	if(C10p==0.) C10p+=param->KC10pmu;
	if(CQ1==0.) CQ1+=param->KCQ1mu;
	if(CQ1p==0.) CQ1p+=param->KCQ1pmu;
	if(CQ2==0.) CQ2+=param->KCQ2mu;
	if(CQ2p==0.) CQ2p+=param->KCQ2pmu;
	
	double ml;
	if(gen==2)
	{
		ml = param->mass_mu;
		if(C10==0.) C10+=param->KC10mu;
		if(C10p==0.) C10p+=param->KC10pmu;
		if(CQ1==0.) CQ1+=param->KCQ1mu;
		if(CQ1p==0.) CQ1p+=param->KCQ1pmu;
		if(CQ2==0.) CQ2+=param->KCQ2mu;
		if(CQ2p==0.) CQ2p+=param->KCQ2pmu;
	}
	else if(gen==1)
	{
		ml = param->mass_e;
		if(C10==0.) C10+=param->KC10e;
		if(C10p==0.) C10p+=param->KC10pe;
		if(CQ1==0.) CQ1+=param->KCQ1e;
		if(CQ1p==0.) CQ1p+=param->KCQ1pe;
		if(CQ2==0.) CQ2+=param->KCQ2e;
		if(CQ2p==0.) CQ2p+=param->KCQ2pe;
	}

	int SDswitch, LDswitch;
	/* SL switch for having:
	 * 1) only short-distance (SD)
	 * 2) only long-distance (LD)
	 * 3) full amplitude (both SD and LD)
	 */
	if(SL==1) {SDswitch = 1; LDswitch = 0;} 
	else if(SL==2) {SDswitch = 0; LDswitch = 1;}
	else {SDswitch = 1; LDswitch = 1;}
	
	int sign;
	if (param->LDsign_KLmumu < 0) sign = -1;
	else sign = 1;

	double life_KL = 5.116e-8; /* err: +-0.0021e-8 (PDG)*/ /*To be added to leshouches.c*/

	double complex Vts,Vtd,Vcs,Vcd,lambdat,lambdac;
	Vts = param->Vts;
	Vtd = param->Vtd;
	Vcs = param->Vcs;
	Vcd = param->Vcd;
	lambdat = Vtd*conj(Vts);
	lambdac = Vcd*conj(Vcs);
	double lambda = cabs(param->Vus);

	double sw2=pow(sin(atan(param->gp/param->g2)),2.);

	double mtmt = param->mtmt;
	double alphas_mut = alphas_running(mtmt,mtmt,param->mass_b_pole,param);
	double xt = pow(mtmt/param->mass_W,2.);

	double YtSM = Y0(xt) + alphas_mut/(4.*pi) * Y1(xt,mtmt,param->mass_W);
	double C10SM = -YtSM/sw2;

// 	/* method 1 for PcY; needs integration; alternatively method 2 can be used below */
	double PcY, s, mucmin, mucmax;
	PcY = 0.;
	mucmin = 1.;
	mucmax = 3.; 
	
	double epsilon = 1.e-10;
	int ie, nmax;
	nmax = 20;
	for(ie=1;ie<=nmax;ie++) {
		s=mucmin+(mucmax-mucmin)*ie/nmax;
		PcY+=Pc_Yscales(lambda, param->mass_c, param->alphas_MZ, s, 5.+epsilon );
// 		PcY+=Pc_Yscales(0.225+epsilon, 1.3+epsilon, 0.1187+epsilon, s, 5.+epsilon);/* For  crosscheck; 0605203 inputs */
	}
	PcY*=1./nmax;
// 	/* method 2 for PcY */
// 	double PcY = Pc_Ylambda(lambda);
// 	printf("PcY\t\t\t%.4f\n",PcY);

	double YcSM = pow(lambda,4.)*PcY;
	
	
	double complex CASM = 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)/sw2*(lambdat*C10SM*sw2 - lambdac*YcSM);/*eq.2.28 in 1711.11030, but with the prefactor defined such that there is full agreement with eq.2.31 in 1712.01295 */

	double complex CA, CAp;
	//CA, CAp expressed in terms of C10, C10p
	CA = CASM + 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*C10;
	CAp =  1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*C10p;
	
	double complex CS, CSp, CP, CPp;
	CS = 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*(CQ1/param->mass_s);
	CSp = 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*(CQ1p/param->mass_s);
	CP = 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*(CQ2/param->mass_s);
	CPp = 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*(CQ2p/param->mass_s);
	
	double complex ASmu, BSmu, NLDSmu;
	
	double complex ALmu, BLmu, NLDLmu;

	/*Long-distance effect in expanded form*/
	double BRKLggEXP = param->BR_KLgammagamma_exp; /*5.47e-4; PDG2020; err: +-0.04e-4 */
	double beta_muK = sqrt(1.-4.*pow(param->mass_mu,2.)/pow(param->m_K0,2.));
	double chi_abs = pi/2./beta_muK *log((1.-beta_muK)/(1.+beta_muK));





	double Mrho = 0.770;
	double chi_gg_Mrho = param->chi_gg_Mrho;
	double beta_l = sqrt(1.-4.*pow(ml,2.)/pow(param->m_K0,2.));
	double complex C_gg = 1./beta_l * (Li2( (beta_l-1.)/(beta_l+1.) ) + pi*pi/3. + 1./4.*cpow(clog((beta_l-1.)/(beta_l+1.)),2.));
	double chi_abs_l = cimag(C_gg);
	double chi_disp_l = creal(C_gg)- 5./2. + 3./2.*log(ml*ml/Mrho/Mrho)+ chi_gg_Mrho;
	double complex Aterm_l = chi_disp_l + I*chi_abs_l;
	double alpha_0 = 1./137.04;
	double complex NLDL_l = sign* (4.*alpha_0*ml)/(pi*param->f_K*param->m_K0*param->m_K0)*sqrt(2.*pi/param->m_K0*BRKLggEXP*hbar/life_KL)*Aterm_l;
	double complex AL_l = LDswitch*NLDL_l + SDswitch*(- param->mass_s*param->m_K0/(param->mass_s+param->mass_d)*creal(CP - CPp) - 2.*ml/param->m_K0*creal(CA - CAp));
	double complex BL_l = SDswitch*(param->mass_s*param->m_K0/(param->mass_s+param->mass_d)*cimag(CS - CSp));

	return life_KL/hbar*pow(param->f_K,2.)*pow(param->m_K0,3.)*beta_l/16./pi*( pow(cabs(AL_l),2.) + pow(beta_l,2.)*pow(cabs(BL_l),2.) );
}

/**************************************************************/

/*arXiv: 1711.11030*/
double BRKSmumu(int SL, double complex C10, double complex C10p, double complex CQ1, double complex CQ1p, double complex CQ2, double complex CQ2p, struct parameters* param)
{
	if(C10==0.) C10+=param->KC10mu;
	if(C10p==0.) C10p+=param->KC10pmu;
	if(CQ1==0.) CQ1+=param->KCQ1mu;
	if(CQ1p==0.) CQ1p+=param->KCQ1pmu;
	if(CQ2==0.) CQ2+=param->KCQ2mu;
	if(CQ2p==0.) CQ2p+=param->KCQ2pmu;
	
	int SDswitch, LDswitch;
	if(SL==1) {SDswitch = 1; LDswitch = 0;} 
	else if(SL==2) {SDswitch = 0; LDswitch = 1;}
	else {SDswitch = 1; LDswitch = 1;}

	double life_KS = 0.8954e-10; /* err: +-0.0004e-10 (PDG)*/ /*To be added to leshouches.c*/

	double complex Vts,Vtd,Vcs,Vcd,lambdat,lambdac;
	Vts = param->Vts;
	Vtd = param->Vtd;
	Vcs = param->Vcs;
	Vcd = param->Vcd;
	lambdat = Vtd*conj(Vts);
	lambdac = Vcd*conj(Vcs);
	double lambda = cabs(param->Vus);

	double sw2=pow(sin(atan(param->gp/param->g2)),2.);
	double mtmt = param->mtmt;
		
	double alphas_mut = alphas_running(mtmt,mtmt,param->mass_b_pole,param);
	
	double xt = pow(mtmt/param->mass_W,2.);
	double YtSM = Y0(xt) + alphas_mut/(4.*pi) * Y1(xt,mtmt,param->mass_W);
	double C10SM = -YtSM/sw2;

// 	/* method 1 for PcY; needs integration; alternatively method 2 can be used below */
	double PcY, s, mucmin, mucmax;
	PcY = 0.;
	mucmin = 1.;
	mucmax = 3.;

	double epsilon = 1.e-10;
	int ie, nmax;
	nmax = 20;
	double step = (mucmax-mucmin)/nmax;
	for(ie=1;ie<=nmax;ie++) {
		s=mucmin+step*ie;
		PcY+=Pc_Yscales(lambda, param->mass_c, param->alphas_MZ, s, 5.+epsilon );
// 		PcY+=Pc_Yscales(0.225+epsilon, 1.3+epsilon, 0.1187+epsilon, s, 5.+epsilon);/* For  crosscheck; 0605203 inputs */
	}
	PcY*=step/(mucmax-mucmin);
// 	/* method 2 for PcY */
// 	double PcY = Pc_Ylambda(lambda);
// 	printf("PcY\t\t\t%.4f\n",PcY);

	double YcSM = pow(lambda,4.)*PcY;
	
	double complex CASM = 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)/sw2*(lambdat*C10SM*sw2 - lambdac*YcSM);/*eq.2.28 in 1711.11030, but with the prefactor defined such that there is full agreement with eq.2.31 in 1712.01295 */

	double complex CA, CAp;
	//CA, CAp expressed in terms of C10, C10p
	CA = CASM + 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*C10;
	CAp =  1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*C10p;
	
	
	double complex CS, CSp, CP, CPp;
	CS = 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*(CQ1/param->mass_s);
	CSp = 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*(CQ1p/param->mass_s);
	CP = 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*(CQ2/param->mass_s);
	CPp = 1./param->inv_alpha_em/pi*param->Gfermi/sqrt(2.)*lambdat*(CQ2p/param->mass_s);
	
	double complex ASmu, BSmu, NLDSmu;
	
// 	NLDSmu = -2.65e-11 + I*1.14e-11;
	/*Long-distance effect in expanded form*/
	double BRKSggEXP = param->BR_KSgammagamma_exp; /* 2.63e-6; PDG2020; err: +-0.17e-6 + 30% UNCERTAINTY ON THE CENTRAL VALUE (see texts after eq. in hep-ph/0311084 and eq.13 in 1707.06999)*/
	double complex Iterm_mu = param->Iterm_mu_KSmumu /* -2.821 */ + I*1.216; /* err: +- 0.001; see eq. 34 in 0404127 (in general a complicated  3d integral)*/
	double r_pi= param->m_pi/param->m_K0;
	double alpha_0 = 1./137.04; /*in principle this should be equal to 1./param->inv_alpha_em but gives different result (scale unclear to me)*/
	NLDSmu = (2.*alpha_0*param->mass_mu)/(pi*param->f_K*param->m_K0*param->m_K0*cabs(Hkaon(0.,r_pi)))*sqrt(2.*pi/param->m_K0*BRKSggEXP*hbar/life_KS)*Iterm_mu;
// 	printf("NLDSmu = %.2e + i*(%.2e)\n", creal(NLDSmu), cimag(NLDSmu));

	ASmu = SDswitch*(param->mass_s*param->m_K0/(param->mass_s+param->mass_d)*cimag(CP - CPp) +2.*param->mass_mu/param->m_K0*cimag(CA - CAp));
	BSmu = LDswitch*NLDSmu +SDswitch*(- param->mass_s*param->m_K0/(param->mass_s+param->mass_d)*creal(CS - CSp));

	double beta_mu=sqrt(1.-4.*pow(param->mass_mu,2.)/pow(param->m_K0,2.));

	return life_KS/hbar*pow(param->f_K,2.)*pow(param->m_K0,3.)*beta_mu/16./pi*( pow(cabs(ASmu),2.) + pow(beta_mu,2.)*pow(cabs(BSmu),2.) );
}

/*----------------------------------------------------------------------*/


/**************************************************************/
/***                   KL->pi0ll                            ***/
/**************************************************************/

double BRKLpill(int gen, double complex C9, double complex C9p, double complex C10, double complex C10p, double complex CQ1, double complex CQ1p, double complex CQ2, double complex CQ2p, struct parameters* param)
{
	double ml=0.;
	if(gen==2)
	{
		ml = param->mass_mu;// For (pseudo)scalar contributions - 06/09/2025
		if(C9==0.) C9+=param->KC9mu;
		if(C9p==0.) C9p+=param->KC9pmu;
		if(C10==0.) C10+=param->KC10mu;
		if(C10p==0.) C10p+=param->KC10pmu;
		/* Scalar and Pseudoscalar contributions - 06/09/2025 */
		if(CQ1==0.) CQ1+=param->KCQ1mu;
		if(CQ1p==0.) CQ1p+=param->KCQ1pmu;
		if(CQ2==0.) CQ2+=param->KCQ2mu;
		if(CQ2p==0.) CQ2p+=param->KCQ2pmu;

	}
	else if(gen==1)
	{
		ml = param->mass_e;// For (pseudo)scalar contributions - 06/09/2025
		if(C9==0.) C9+=param->KC9e;
		if(C9p==0.) C9p+=param->KC9pe;
		if(C10==0.) C10+=param->KC10e;
		if(C10p==0.) C10p+=param->KC10pe;
		/* Scalar and Pseudoscalar contributions - 06/09/2025 */
		if(CQ1==0.) CQ1+=param->KCQ1e;
		if(CQ1p==0.) CQ1p+=param->KCQ1pe;
		if(CQ2==0.) CQ2+=param->KCQ2e;
		if(CQ2p==0.) CQ2p+=param->KCQ2pe;
	}
	
// 	sign = + or -;	The +/- sign correponds to whether consntructive or destructive interference between the direct and indirect
// CP-violating contributions is assumed. Constructive interference is theoretically preferred (see 0308008).
	double complex Vts,Vtd,lambdat;
	Vts = param->Vts;
	Vtd = param->Vtd;
	lambdat = Vtd*conj(Vts);

//	double C9SM_muK = 4.587; /* In principle when the full RGE running is implemented, the value can be taken from Wilson.c */
//	double C10SM_muK = -4.291; /* C10 has no running and only the EW corrections would be added, so already can be taken from Wilson.c */

	double complex C0K[11],C1K[11],C2K[11],C0w[11],C1w[11],C2w[11];
	
	double mu_W=param->mass_W;
	double mu_K=pow(2.,param->log_muK_1GeV);
				
	CW_calculator(gen,C0w,C1w,C2w,mu_W,param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0K,C1K,C2K,mu_K,param);
	
	double alphas_muK=alphas_running(mu_K,param->mass_top_pole,param->mass_b_pole,param);
	
	double complex C9_muK=C0K[9]+alphas_muK/4./pi*C1K[9]+alphas_muK*alphas_muK/16./pi/pi*C2K[9];
	double complex C10_muK=C0K[10]+alphas_muK/4./pi*C1K[10]+alphas_muK*alphas_muK/16./pi/pi*C2K[10];


	double Ce_CPC = 0.;
	int interf = param->sign_KLpill;
	double abs_aS = param->KLpill_abs_aS;


	double omega7V = 1./2./pi * cimag((C9_muK + C9 + C9p) * lambdat / 1.407e-4);
	double omega7A = 1./2./pi * cimag((C10_muK + C10 + C10p) * lambdat / 1.407e-4);
	double Cdir,Cint,Cmix,CPC;


	double sw2=pow(sin(atan(param->gp/param->g2)),2.);
	double complex yS = - pow(param->mass_W,2.)*sw2/ml/param->mass_s * lambdat * (CQ1 + CQ1p)/2. ;
	double complex yP = - pow(param->mass_W,2.)*sw2/ml/param->mass_s * lambdat * (CQ2 + CQ2p)/2. ;


	double Cmix_Sgg, Cdir_S2, Cdir_PA, Cdir_P2;


	if(gen==1)
	{
		Cdir = param->KLpill_Ce_Cdir*(pow(omega7V,2.) + pow(omega7A,2.));
		Cint = param->KLpill_Ce_Cint*omega7V;
		Cmix = param->KLpill_Ce_Cmix;
		CPC = Ce_CPC;

		Cdir_PA = param->KLpill_Ce_PA;
		Cdir_P2 = param->KLpill_Ce_PP;
		Cmix_Sgg = param->KLpill_Ce_Sgg;
		Cdir_S2 = param->KLpill_Ce_SS;

	}
	else
	{
		Cdir = param->KLpill_Cmu_Cdir*(pow(omega7V,2.) + 2.32*pow(omega7A,2.));
		Cint = param->KLpill_Cmu_Cint*omega7V;
		Cmix = param->KLpill_Cmu_Cmix;
		CPC = param->KLpill_Cmu_CPC;

		Cdir_PA = param->KLpill_Cmu_PA;
		Cdir_P2 = param->KLpill_Cmu_PP;
		Cmix_Sgg = param->KLpill_Cmu_Sgg;
		Cdir_S2 = param->KLpill_Cmu_SS;
	}


	double BR_VA = (Cdir + interf*Cint*abs_aS + Cmix*pow(abs_aS,2.) + CPC)*1.e-12;
	double BR_S = (Cmix_Sgg * creal(yS) + Cdir_S2 * pow(creal(yS),2.))*1.e-12;
	double BR_P = (Cdir_PA * cimag(yP) + Cdir_P2 * pow(cimag(yP),2.))*1.e-12;


	return BR_VA + BR_S + BR_P;
}

/*----------------------------------------------------------------------*/


double AFBmu_KLpill(double complex C9, double complex C9p, double complex C10, double complex C10p, double complex CQ1, double complex CQ1p, double complex CQ2, double complex CQ2p, struct parameters* param)
{
	double ml = param->mass_mu;
	if(C9==0.) C9+=param->KC9mu;
	if(C9p==0.) C9p+=param->KC9pmu;
	if(C10==0.) C10+=param->KC10mu;
	if(C10p==0.) C10p+=param->KC10pmu;
	/* Scalar and Pseudoscalar contributions - 06/09/2025 */
	if(CQ1==0.) CQ1+=param->KCQ1mu;
	if(CQ1p==0.) CQ1p+=param->KCQ1pmu;
	if(CQ2==0.) CQ2+=param->KCQ2mu;
	if(CQ2p==0.) CQ2p+=param->KCQ2pmu;


	double complex Vts,Vtd,lambdat;
	Vts = param->Vts;
	Vtd = param->Vtd;
	lambdat = Vtd*conj(Vts);


	double complex C0K[11],C1K[11],C2K[11],C0w[11],C1w[11],C2w[11];

	double mu_W=param->mass_W;
	double mu_K=pow(2.,param->log_muK_1GeV);

	int gen = 2;
	CW_calculator(gen,C0w,C1w,C2w,mu_W,param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0K,C1K,C2K,mu_K,param);

	double alphas_muK=alphas_running(mu_K,param->mass_top_pole,param->mass_b_pole,param);

	double complex C9_muK=C0K[9]+alphas_muK/4./pi*C1K[9]+alphas_muK*alphas_muK/16./pi/pi*C2K[9];
	double complex C10_muK=C0K[10]+alphas_muK/4./pi*C1K[10]+alphas_muK*alphas_muK/16./pi/pi*C2K[10];


	int interf = param->sign_KLpill;
	double abs_aS = param->KLpill_abs_aS;


	double omega7V = 1./2./pi * cimag((C9_muK + C9 + C9p) * lambdat / 1.407e-4);
	double omega7A = 1./2./pi * cimag((C10_muK + C10 + C10p) * lambdat / 1.407e-4);
	double Cdir,Cint,Cmix,CPC;


	double sw2=pow(sin(atan(param->gp/param->g2)),2.);
	double complex yS = - pow(param->mass_W,2.)*sw2/ml/param->mass_s * lambdat * (CQ1 + CQ1p)/2. ;
	double complex yP = - pow(param->mass_W,2.)*sw2/ml/param->mass_s * lambdat * (CQ2 + CQ2p)/2. ;


	double Cmix_Sgg, Cdir_S2, Cdir_PA, Cdir_P2;


	Cdir = param->KLpill_Cmu_Cdir*(pow(omega7V,2.) + 2.32*pow(omega7A,2.));
	Cint = param->KLpill_Cmu_Cint*omega7V;
	Cmix = param->KLpill_Cmu_Cmix;
	CPC = param->KLpill_Cmu_CPC;

	Cdir_PA = param->KLpill_Cmu_PA;
	Cdir_P2 = param->KLpill_Cmu_PP;
	Cmix_Sgg = param->KLpill_Cmu_Sgg;
	Cdir_S2 = param->KLpill_Cmu_SS;


	double BR_VA = (Cdir + interf*Cint*abs_aS + Cmix*pow(abs_aS,2.) + CPC)*1.e-12;
	double BR_S = (Cmix_Sgg * creal(yS) + Cdir_S2 * pow(creal(yS),2.))*1.e-12;
	double BR_P = (Cdir_PA * cimag(yP) + Cdir_P2 * pow(cimag(yP),2.))*1.e-12;

	double Cmu_gg_7V = param->KLpill_AFB_a1pi_gg_7V;
	double Cmu_gg_aS = param->KLpill_AFB_a1pi_gg_aS;

	return (Cmu_gg_7V*omega7V + interf*Cmu_gg_aS*abs_aS - interf*abs_aS*yS)*1.e-12 / (BR_VA + BR_S + BR_P);
}



