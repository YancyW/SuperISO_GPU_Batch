#include "include.h"

/*----------------------------------------------------------------------*/
/*------------------------------ SOFT ----------------------------------*/
/*----------------------------------------------------------------------*/

double dGamma_Bsphill_dq2_soft(int gen, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{	
	double ml;
	if(gen==1) ml=param->mass_e;
	else if(gen==3) ml=param->mass_tau;
	else ml=param->mass_mu;
	
	double mBs=param->m_Bs;
	double mphi=param->m_phi;
	double eq=-1./3.;

	double mc=mc_pole_1loop(param);
	double mbpole=mb_pole_1loop(param);
	double mb_mub=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);

	int ie,je;
	
	double beta_l=sqrt(1.-4.*ml*ml/q2);

	double alpha_em=1./133.;

	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	
	double mu_f=sqrt(mu_b*0.5);
	
	double alphas_muf=alphas_running(mu_f,param->mass_top_pole,param->mass_b_pole,param);
	double eta=alphas_muf/alphas_running(1.,param->mass_top_pole,param->mass_b_pole,param);

	double Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_phi=(mBs*mBs+mphi*mphi-q2)/2./mBs;

	int nf=5;
	double f_phi_perp=param->f_phi_perp;
	f_phi_perp*=pow(eta,4./3./(11.-2./3.*nf));

	double f_phi_par=param->f_phi_par;
	
	double complex ALperp=0.;
	double complex ARperp=0.;
	double complex ALpar=0.;
	double complex ARpar=0.;
	double complex AL0=0.;
	double complex AR0=0.;
	double complex At=0.;
	double complex AS=0.;

	double complex ALperp_bar=0.;
	double complex ARperp_bar=0.;
	double complex ALpar_bar=0.;
	double complex ARpar_bar=0.;
	double complex AL0_bar=0.;
	double complex AR0_bar=0.;
	double complex At_bar=0.;
	double complex AS_bar=0.;

	double complex ALperp_high=0.;
	double complex ARperp_high=0.;
	double complex ALpar_high=0.;
	double complex ARpar_high=0.;
	double complex AL0_high=0.;
	double complex AR0_high=0.;
	double complex At_high=0.;
	double complex AS_high=0.;

	double complex ALperp_bar_high=0.;
	double complex ARperp_bar_high=0.;
	double complex ALpar_bar_high=0.;
	double complex ARpar_bar_high=0.;
	double complex AL0_bar_high=0.;
	double complex AR0_bar_high=0.;
	double complex At_bar_high=0.;
	double complex AS_bar_high=0.;
		

	double V = 0.0;
	double A1 = 0.0;
	double A2 = 0.0;
	double A0 = 0.0;
	double T1 = 0.0;
	double T2 = 0.0;
	double T23 = 0.0;
	double T3 = 0.0;
	double A12 = 0.0;

	double MB, MV;
	double a0V, a1V, a2V, a0A0, a1A0, a2A0, a0A1, a1A1, a2A1, a0A12, a1A12, a2A12, a0T1, a1T1, a2T1, a0T2, a1T2, a2T2, a0T23, a1T23, a2T23;
	MB=MV=a0V=a1V=a2V=a0A0=a1A0=a2A0=a0A1=a1A1=a2A1=a0A12=a1A12=a2A12=a0T1=a1T1=a2T1=a0T2=a1T2=a2T2=a0T23=a1T23=a2T23 = 0.;

	/********LCSR+Lattice fit from Barucha,Straub,Zwicky 1503.05534***************/
	if(param->Bsphi_FormFactor_choice ==1)
	{
		double P_V=1./(1.-q2/pow(param->MV_Bsphi,2.));
		double P_A1=1./(1.-q2/pow(param->MA1_Bsphi,2.));
		double P_A12=1./(1.-q2/pow(param->MA12_Bsphi,2.));

		double tau_plus=pow(mBs+mphi,2.);
		double tau_minus=pow(mBs-mphi,2.);
		double tau_0=tau_plus*(1.-sqrt(1.-tau_minus/tau_plus));
		double z_q2=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));
		double z_0=(sqrt(tau_plus)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus)+sqrt(tau_plus-tau_0));

		V=P_V*(param->a0V_Bsphi+param->a1V_Bsphi*(z_q2-z_0)+param->a2V_Bsphi*(pow(z_q2-z_0,2.)));
		A1=P_A1*(param->a0A1_Bsphi+param->a1A1_Bsphi*(z_q2-z_0)+param->a2A1_Bsphi*(pow(z_q2-z_0,2.)));
		A12=P_A12*(param->a0A12_Bsphi+param->a1A12_Bsphi*(z_q2-z_0)+param->a2A12_Bsphi*(pow(z_q2-z_0,2.)));
		A2=(pow(mBs+mphi,2.)*(mBs*mBs-mphi*mphi-q2)*A1-16.*mBs*mphi*mphi*(mBs+mphi)*A12)/((pow(mBs+mphi,2.)-q2)*(pow(mBs-mphi,2.)-q2));

		double P_A0=1./(1.-q2/pow(param->MA0_Bsphi,2.));
		double P_T1=1./(1.-q2/pow(param->MT1_Bsphi,2.));
		double P_T2=1./(1.-q2/pow(param->MT2_Bsphi,2.));
		double P_T23=1./(1.-q2/pow(param->MT23_Bsphi,2.));

		A0=P_A0*(param->a0A0_Bsphi+param->a1A0_Bsphi*(z_q2-z_0)+param->a2A0_Bsphi*(pow(z_q2-z_0,2.)));
		T1=P_T1*(param->a0T1_Bsphi+param->a1T1_Bsphi*(z_q2-z_0)+param->a2T1_Bsphi*(pow(z_q2-z_0,2.)));
		T2=P_T2*(param->a0T2_Bsphi+param->a1T2_Bsphi*(z_q2-z_0)+param->a2T2_Bsphi*(pow(z_q2-z_0,2.)));
		T23=P_T23*(param->a0T23_Bsphi+param->a1T23_Bsphi*(z_q2-z_0)+param->a2T23_Bsphi*(pow(z_q2-z_0,2.)));

		// printf("[%.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f],\n",q2, V, A0, A1, A12, T1, T2, T23, A2);

	}

	/********Fit parameters from Gubernari, Reboud, van Dyk, Virto, 2305.06301***************/
	if(param->Bsphi_FormFactor_choice ==2)
	{
		MB = mBs;
		MV = mphi;
		double MVT1_Bsphi = 5.415;
		double MA0_Bsphi = 5.367;
		double MAT_Bsphi = 5.829;
		a0V = param->a0V_Bsphi_GRvDV_BSZ;
		a1V = param->a1V_Bsphi_GRvDV_BSZ;
		a2V = param->a2V_Bsphi_GRvDV_BSZ;
		a0A0 = param->a0A0_Bsphi_GRvDV_BSZ;
		a1A0 = param->a1A0_Bsphi_GRvDV_BSZ;
		a2A0 = param->a2A0_Bsphi_GRvDV_BSZ;
		a0A1 = param->a0A1_Bsphi_GRvDV_BSZ;
		a1A1 = param->a1A1_Bsphi_GRvDV_BSZ;
		a2A1 = param->a2A1_Bsphi_GRvDV_BSZ;
		a0A12 = (MB*MB - MV*MV)/(8.*MB*MV)*param->a0A0_Bsphi_GRvDV_BSZ;// Eq.17 in BSZ [1503.05534]
		a1A12 = param->a1A12_Bsphi_GRvDV_BSZ;
		a2A12 = param->a2A12_Bsphi_GRvDV_BSZ;
		a0T1 = param->a0T1_Bsphi_GRvDV_BSZ;
		a1T1 = param->a1T1_Bsphi_GRvDV_BSZ;
		a2T1 = param->a2T1_Bsphi_GRvDV_BSZ;
		a0T2 = param->a0T1_Bsphi_GRvDV_BSZ;// Eq.17 in BSZ [1503.05534]
		a1T2 = param->a1T2_Bsphi_GRvDV_BSZ;
		a2T2 = param->a2T2_Bsphi_GRvDV_BSZ;
		a0T23 = param->a0T23_Bsphi_GRvDV_BSZ;
		a1T23 = param->a1T23_Bsphi_GRvDV_BSZ;
		a2T23 = param->a2T23_Bsphi_GRvDV_BSZ;

		double P_V=1./(1.-q2/pow(MVT1_Bsphi,2.));
		double P_A0=1./(1.-q2/pow(MA0_Bsphi,2.));
		double P_A1=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_A12=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_T1=1./(1.-q2/pow(MVT1_Bsphi,2.));
		double P_T2=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_T23=1./(1.-q2/pow(MAT_Bsphi,2.));

		double tau_plus=pow(MB+MV,2.);
		double tau_minus=pow(MB-MV,2.);
		double tau_0=tau_plus-sqrt((tau_plus-tau_minus)*tau_plus);
		double z_q2tau0=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));
		double z_0tau0=(sqrt(tau_plus)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus)+sqrt(tau_plus-tau_0));
		double z_diff = z_q2tau0 - z_0tau0;


		V = P_V * ( a0V + a1V*z_diff + a2V* (z_diff*z_diff) );
		A0 = P_A0 * ( a0A0 + a1A0*z_diff + a2A0*(z_diff*z_diff) );
		A1 = P_A1 * ( a0A1 + a1A1*z_diff + a2A1*(z_diff*z_diff) );
		A12 = P_A12 * ( a0A12 + a1A12*z_diff + a2A12*(z_diff*z_diff) );
		T1 = P_T1 * ( a0T1 + a1T1*z_diff + a2T1*(z_diff*z_diff) );
		T2 = P_T2 * ( a0T2 + a1T2*z_diff + a2T2*(z_diff*z_diff) );
		T23 = P_T23 * ( a0T23 + a1T23*z_diff + a2T23*(z_diff*z_diff) );

		A2=(pow(MB+MV,2.)*(MB*MB-MV*MV-q2)*A1-16.*MB*MV*MV*(MB+MV)*A12)/((pow(MB+MV,2.)-q2)*(pow(MB-MV,2.)-q2));

		// printf("[%.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f],\n",q2, V, A0, A1, A12, T1, T2, T23, A2);
	}

	/********LCSR-only from Barucha,Straub,Zwicky 1503.05534***************/
	if(param->Bsphi_FormFactor_choice ==5)
	{
		MB = mBs;
		MV = mphi;
		a0V = param->a0V_Bsphi_BSZ_LCSRonly;
		a1V = param->a1V_Bsphi_BSZ_LCSRonly;
		a2V = param->a2V_Bsphi_BSZ_LCSRonly;
		a0A0 = param->a0A0_Bsphi_BSZ_LCSRonly;
		a1A0 = param->a1A0_Bsphi_BSZ_LCSRonly;
		a2A0 = param->a2A0_Bsphi_BSZ_LCSRonly;
		a0A1 = param->a0A1_Bsphi_BSZ_LCSRonly;
		a1A1 = param->a1A1_Bsphi_BSZ_LCSRonly;
		a2A1 = param->a2A1_Bsphi_BSZ_LCSRonly;
		a0A12 = param->a0A12_Bsphi_BSZ_LCSRonly;
		a1A12 = param->a1A12_Bsphi_BSZ_LCSRonly;
		a2A12 = param->a2A12_Bsphi_BSZ_LCSRonly;
		a0T1 = param->a0T1_Bsphi_BSZ_LCSRonly;
		a1T1 = param->a1T1_Bsphi_BSZ_LCSRonly;
		a2T1 = param->a2T1_Bsphi_BSZ_LCSRonly;
		a0T2 = param->a0T2_Bsphi_BSZ_LCSRonly;
		a1T2 = param->a1T2_Bsphi_BSZ_LCSRonly;
		a2T2 = param->a2T2_Bsphi_BSZ_LCSRonly;
		a0T23 = param->a0T23_Bsphi_BSZ_LCSRonly;
		a1T23 = param->a1T23_Bsphi_BSZ_LCSRonly;
		a2T23 = param->a2T23_Bsphi_BSZ_LCSRonly;

		double P_V=1./(1.-q2/pow(param->MV_Bsphi,2.));
		double P_A0=1./(1.-q2/pow(param->MA0_Bsphi,2.));
		double P_A1=1./(1.-q2/pow(param->MA1_Bsphi,2.));
		double P_A12=1./(1.-q2/pow(param->MA12_Bsphi,2.));
		double P_T1=1./(1.-q2/pow(param->MT1_Bsphi,2.));
		double P_T2=1./(1.-q2/pow(param->MT2_Bsphi,2.));
		double P_T23=1./(1.-q2/pow(param->MT23_Bsphi,2.));

		double tau_plus=pow(MB+MV,2.);
		double tau_minus=pow(MB-MV,2.);
		double tau_0=tau_plus-sqrt((tau_plus-tau_minus)*tau_plus);
		double z_q2tau0=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));
		double z_0tau0=(sqrt(tau_plus)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus)+sqrt(tau_plus-tau_0));
		double z_diff = z_q2tau0 - z_0tau0;


		V = P_V * ( a0V + a1V*z_diff + a2V* (z_diff*z_diff) );
		A0 = P_A0 * ( a0A0 + a1A0*z_diff + a2A0* (z_diff*z_diff) );
		A1 = P_A1 * ( a0A1 + a1A1*z_diff + a2A1* (z_diff*z_diff) );
		A12 = P_A12 * ( a0A12 + a1A12*z_diff + a2A12* (z_diff*z_diff) );
		T1 = P_T1 * ( a0T1 + a1T1*z_diff + a2T1* (z_diff*z_diff) );
		T2 = P_T2 * ( a0T2 + a1T2*z_diff + a2T2* (z_diff*z_diff) );
		T23 = P_T23 * ( a0T23 + a1T23*z_diff + a2T23* (z_diff*z_diff) );

		A2=(pow(MB+MV,2.)*(MB*MB-MV*MV-q2)*A1-16.*MB*MV*MV*(MB+MV)*A12)/((pow(MB+MV,2.)-q2)*(pow(MB-MV,2.)-q2));

		// printf("[%.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f,  %.3f],\n",q2, V, A0, A1, A12, T1, T2, T23, A2);

	}

	/********Fit parameters from Horgan, Liu, Meinel, Wingate 1310.3722 + 1501.00367 ***************/
	if(param->Bsphi_FormFactor_choice == 6)
	{
		MB = mBs;
		MV = mphi;
		double MVT1_Bsphi = 5.411;//Mres=mBs+deltam = 5.366 + 0.045
		double MA0_Bsphi = 5.366;//Mres=mBs+deltam = 5.366 + 0.000
		double MAT_Bsphi = 5.806;//Mres=mBs+deltam = 5.366 + 0.440
		a0V   = param->a0V_Bsphi_HLMW;
		a1V   = param->a1V_Bsphi_HLMW;
		a0A0  = param->a0A0_Bsphi_HLMW;
		a1A0  = param->a1A0_Bsphi_HLMW;
		a0A1  = param->a0A1_Bsphi_HLMW;
		a1A1  = param->a1A1_Bsphi_HLMW;
		a0A12 = param->a0A12_Bsphi_HLMW;
		a1A12 = param->a1A12_Bsphi_HLMW;
		a0T1  = param->a0T1_Bsphi_HLMW;
		a1T1  = param->a1T1_Bsphi_HLMW;
		a0T2  = param->a0T2_Bsphi_HLMW;
		a1T2  = param->a1T2_Bsphi_HLMW;
		a0T23 = param->a0T23_Bsphi_HLMW;
		a1T23 = param->a1T23_Bsphi_HLMW;
		double sysErr = param->systErr_Bsphi_HLMW;

		double P_V=1./(1.-q2/pow(MVT1_Bsphi,2.));
		double P_A0=1./(1.-q2/pow(MA0_Bsphi,2.));
		double P_A1=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_A12=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_T1=1./(1.-q2/pow(MVT1_Bsphi,2.));
		double P_T2=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_T23=1./(1.-q2/pow(MAT_Bsphi,2.));

		double tau_plus=pow(MB+MV,2.);
		double tau_minus=pow(MB-MV,2.);
		double tau_0=12.;// Page 9 in [1310.3722]
		double z_q2tau0=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));

		/* The (1. + sysErr) is there to take into account 5% systematic error (e.g. see p.18 in 1310.3722)*/
		V = P_V * ( a0V + a1V*z_q2tau0  ) * (1. + sysErr);
		A0 = P_A0 * ( a0A0 + a1A0*z_q2tau0  ) * (1. + sysErr);
		A1 = P_A1 * ( a0A1 + a1A1*z_q2tau0  ) * (1. + sysErr);
		A12 = P_A12 * ( a0A12 + a1A12*z_q2tau0  ) * (1. + sysErr);
		T1 = P_T1 * ( a0T1 + a1T1*z_q2tau0  ) * (1. + sysErr);
		T2 = P_T2 * ( a0T2 + a1T2*z_q2tau0  ) * (1. + sysErr);
		T23 = P_T23 * ( a0T23 + a1T23*z_q2tau0 ) * (1. + sysErr);

		A2=(pow(MB+MV,2.)*(MB*MB-MV*MV-q2)*A1-16.*MB*MV*MV*(MB+MV)*A12)/((pow(MB+MV,2.)-q2)*(pow(MB-MV,2.)-q2));

	}



/****************************************************************************/

	double xi_par=(mBs+mphi)/2./E_phi*A1-(mBs-mphi)/mBs*A2;
	double xi_perp=mBs/(mBs+mphi)*V;
		

	if(q2<HIGH_q2)
	{
	 	double complex C7eff=Cmub[7];
		double complex C8eff=Cmub[8];
		double complex C9=Cmub[9];
		double complex C10=Cmub[10];
	
		double complex C7effp=Cpb[7];
		double complex C9p=Cpb[9];
		double complex C10p=Cpb[10];
	
		double complex C1bar=Cmub[1]/2.;
		double complex C2bar=Cmub[2]-Cmub[1]/6.;
		double complex C3bar=Cmub[3]-Cmub[4]/6.+16.*Cmub[5]-8./3.*Cmub[6];
		double complex C4bar=Cmub[4]/2.+8.*Cmub[6];
		double complex C5bar=Cmub[3]-Cmub[4]/6.+4.*Cmub[5]-2./3.*Cmub[6];
		double complex C6bar=Cmub[4]/2.+2.*Cmub[6];
	
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];

		double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
		double mb=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */
				
					
		double complex h_mc=h_bkll(q2,mc,mu_b);
		double complex h_mb=h_bkll(q2,mbpole,mu_b);
		double complex h_0=h_bkll(q2,0.,mu_b);

		double complex Y=4./3.*Cmub[3]+64./9.*Cmub[5]+64./27.*Cmub[6]
		+h_mc*(4./3.*Cmub[1]+Cmub[2]+6.*Cmub[3]+60.*Cmub[5])
		+h_mb*(-7./2.*Cmub[3]-2./3.*Cmub[4]-38.*Cmub[5]-32./3.*Cmub[6])
		+h_0*(-1./2.*Cmub[3]-2./3.*Cmub[4]-8.*Cmub[5]-32./3.*Cmub[6]);

		double Yu=(h_mc-h_0)*(4./3.*Cmub[1]+Cmub[2]);

		double complex Cperpp0=C7eff+C7effp+q2/2./mb/mBs*Y;
		double complex Cperpm0=C7eff-C7effp+q2/2./mb/mBs*Y;
		//double complex Cparp0=-(C7eff+C7effp)-mBs/2./mb*Y;
		double complex Cparm0=-(C7eff-C7effp)-mBs/2./mb*Y;
		double complex Cperpp0_bar=conj(C7eff+C7effp)+q2/2./mb/mBs*Y;
		double complex Cperpm0_bar=conj(C7eff-C7effp)+q2/2./mb/mBs*Y;
		//double complex Cparp0_bar=-conj(C7eff+C7effp)-mBs/2./mb*Y;
		double complex Cparm0_bar=-conj(C7eff-C7effp)-mBs/2./mb*Y;


// 		double complex Cperpp0u=q2/2./mb/mBs*Yu;
// 		double complex Cperpm0u=q2/2./mb/mBs*Yu;
// 		//double complex Cparp0u=-mBs/2./mb*Yu;
// 		double complex Cparm0u=-mBs/2./mb*Yu;
	        /*******Used in redefinition of Tauperpu**********/
		double complex Cperp0u=q2/2./mb/mBs*Yu;
		double complex Cpar0u=-mBs/2./mb*Yu;

		double logb=log(mu_b/mb);
		double DeltaM=-6.*logb-4.*(1.-mu_f/mb);
		double L=-(mb*mb-q2)/q2*log(1.-q2/mb/mb);
	
		double shat=q2/mb/mb; 

		double mchat=mc/mb; 
		double z=mchat*mchat;	

		double complex Cperppf=(C7eff+C7effp)*(-2.*logb-L+DeltaM);
		double complex Cperpmf=(C7eff-C7effp)*(-2.*logb-L+DeltaM);
		//double complex Cparpf=-(C7eff+C7effp)*(-2.*logb+2.*L+DeltaM);
		double complex Cparmf=-(C7eff-C7effp)*(-2.*logb+2.*L+DeltaM);
		double complex Cperppf_bar=conj(C7eff+C7effp)*(-2.*logb-L+DeltaM);
		double complex Cperpmf_bar=conj(C7eff-C7effp)*(-2.*logb-L+DeltaM);
		//double complex Cparpf_bar=-conj(C7eff+C7effp)*(-2.*logb+2.*L+DeltaM);
		double complex Cparmf_bar=-conj(C7eff-C7effp)*(-2.*logb+2.*L+DeltaM);

		double complex F27=F27_bsll(shat,z,logb);
		double complex F87=F87_bsll(shat,logb);
		double complex F29=F29_bkll(shat,z,logb);
		double complex F19=F19_bkll(shat,z,logb);
		double complex F89=F89_bsll(shat);
		double complex F27_u=F27u(shat,logb);
		double complex F29_u=F29u(shat,logb);
		double complex F19_u=F19u(shat,logb);

		double complex Cperpnf=(-C2bar*F27-C8eff*F87
		-q2/2./mb/mBs*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;
		
		double complex Cparnf=(C2bar*F27+C8eff*F87
		+mBs/2./mb*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;

		double complex Cperpnfu=(-C2bar*(F27+F27_u)
		-q2/2./mb/mBs*(C2bar*(F29+F29_u)+2.*C1bar*(F19+F19_u+1./6.*(F29+F29_u))))/4.*3.;
		
		double complex Cparnfu=(C2bar*(F27+F27_u)
		+mBs/2./mb*(C2bar*(F29+F29_u)+2.*C1bar*(F19+F19_u+1./6.*(F29+F29_u))))/4.*3.;
		
		double complex Cperpp1=Cperppf+Cperpnf;
		double complex Cperpm1=Cperpmf+Cperpnf;
		//double complex Cparp1=Cparpf+Cparnf;
		double complex Cparm1=Cparmf+Cparnf;
		double complex Cperpp1_bar=Cperppf_bar+Cperpnf;
		double complex Cperpm1_bar=Cperpmf_bar+Cperpnf;
		//double complex Cparp1_bar=Cparpf_bar+Cparnf;
		double complex Cparm1_bar=Cparmf_bar+Cparnf;

		double complex Cperpp=Cperpp0+alphas_mub*4./3./4./pi*Cperpp1;
		double complex Cperpm=Cperpm0+alphas_mub*4./3./4./pi*Cperpm1;
		//double complex Cparp=Cparp0+alphas_mub*4./3./4./pi*Cparp1;
		double complex Cparm=Cparm0+alphas_mub*4./3./4./pi*Cparm1;
		double complex Cperpp_bar=Cperpp0_bar+alphas_mub*4./3./4./pi*Cperpp1_bar;
		double complex Cperpm_bar=Cperpm0_bar+alphas_mub*4./3./4./pi*Cperpm1_bar;
		//double complex Cparp=Cparp0_bar+alphas_mub*4./3./4./pi*Cparp1_bar;
		double complex Cparm_bar=Cparm0_bar+alphas_mub*4./3./4./pi*Cparm1_bar;

// 		double complex Cperppu=Cperpp0u+alphas_mub*4./3./4./pi*Cperpnfu; 
// 		double complex Cperpmu=Cperpm0u+alphas_mub*4./3./4./pi*Cperpnfu; 
// 		//double complex Cparpu=Cparp0u+alphas_mub*4./3./4./pi*Cparnfu; 
// 		double complex Cparmu=Cparm0u+alphas_mub*4./3./4./pi*Cparnfu; 
	        /*******Used in redefinition of Tauperpu**********/
		double complex Cperpu=Cperp0u+alphas_mub*4./3./4./pi*Cperpnfu; 
		double complex Cparu=Cpar0u+alphas_mub*4./3./4./pi*Cparnfu; 
		
		
		double Xi_perp=1.;
		double Xi_par=mphi/E_phi;
	
		double eu=2./3.;
		double ed=-1./3.;
	
		double a1phi_perp=param->a1phi_perp;
		double a2phi_perp=param->a2phi_perp;
		double a1phi_par=param->a1phi_par;
		double a2phi_par=param->a2phi_par;

		a1phi_perp*=pow(eta,4./3.*(4.*1./2.)/(11.-2./3.*nf));
		a2phi_perp*=pow(eta,4./3.*(4.*(1./2.+1./3.))/(11.-2./3.*nf));

		a1phi_par*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
		a2phi_par*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

		double u;
	
		double complex int_perppp,int_perppm,int_perpmp,int_perpmm;
		double complex /*int_parpp,*/int_parpm,int_parmp,int_parmm;//int_parpp commented; not used
		double complex int_perppu,int_parpu,int_parmu;
		double complex Tperppp0,Tperpppf,Tperpppnf,Tperppp;
		double complex Tperppm0,Tperppmf,Tperppmnf,Tperppm;
		double complex Tperpmp0,Tperpmpf,Tperpmpnf,Tperpmp;
		double complex Tperpmm0,Tperpmmf,Tperpmmnf,Tperpmm;
		double complex Tparpp0/*,Tparppf*/,Tparppnf/*,Tparpp*/;// Tparppf and Tparpp commented; not used
		double complex Tparpm0,Tparpmf,Tparpmnf,Tparpm;
		double complex Tparmp0,Tparmpf,Tparmpnf,Tparmp;
		double complex Tparmm0,Tparmmf,Tparmmnf,Tparmm;
		double complex Tperppnfu,Tperppu;
		double complex Tparp0u,Tparpnfu,Tparpu;
		double complex Tparm0u,Tparmnfu,Tparmu;
		/* perp_bar */
		double complex int_perppp_bar;
		double complex int_perpmp_bar;
		double complex Tperppp_bar,Tperpmp_bar;
		double complex Tperpppf_bar,Tperpmpf_bar;
		/* par_bar */
		// double complex int_parpp_bar;// commented; not used
		double complex int_parmp_bar;
		double complex /*Tparpp_bar,*/Tparmp_bar;//  Tparpp_bar commented; not used
		double complex /*Tparppf_bar,*/Tparmpf_bar;// Tparppf_bar commented; not used

		int_perppp=int_perppm=int_perpmp=int_perpmm=0.;
		int_perppu=0.;
		/*int_parpp=*/int_parpm=int_parmp=int_parmm=0.;// int_parpp commented; not used
		int_parpu=int_parmu=0.;
		int_perppp_bar=int_perpmp_bar=0.;
		/*int_parpp_bar=*/int_parmp_bar=0.;// int_parpp_bar commented; not

		double lambda_Bp=param->lambda_Bsp; 
		lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_f,2.))*(1.-2.*1.4);

		double omega0=2.*(mBs-mb)/3.; 
		double complex lambda_Bm=1./(exp(-q2/mBs/omega0)/omega0*(-Ei(q2/mBs/omega0)+I*pi));

		double phiphi_perp,phiphi_par;
		double complex tperp_mb,tperp_mc,tperp_0;
		double complex tpar_mb,tpar_mc,tpar_0;

		double complex integ3=0;

		double complex Fperp=0.;
		double complex Xperp=0.;
		double x;
		double complex integ4=0.;
		double complex FV;
		double complex integ4u=0.;
		double complex FVu;
		double complex integ5=0.;
		double complex integ5u=0.;

		double zeta3A=param->zeta3A;
		double zeta3V=param->zeta3V;
		double wA10=param->wA10;
		double deltatp=param->deltatp_phi;
		double deltatm=param->deltatm_phi;
		
		int n1=10;
		int n1sav=n1;
		for(ie=0;ie<=n1;ie++)
		{
			u=(double)ie/n1;
			if(ie==0) n1*=2;
			if(ie==n1){u=0.99;n1*=2;}

		/* Tperp */		
			Tperppp0=Tperpmp0=0.;
			Tperpppf=(C7eff+C7effp)*2.*mBs/E_phi/(1.-u);
			Tperpmpf=(C7eff-C7effp)*2.*mBs/E_phi/(1.-u);
			Tperpppf_bar=conj(C7eff+C7effp)*2.*mBs/E_phi/(1.-u);
			Tperpmpf_bar=conj(C7eff-C7effp)*2.*mBs/E_phi/(1.-u);

			Tperppm0=Tperpmm0=0.;
			Tperppmf=Tperpmmf=0.;
			
			phiphi_perp=phi_Kstar(u,a1phi_perp,a2phi_perp);
			tperp_mc=tperp_bkll(u,mc,q2,E_phi,param);
			tperp_mb=tperp_bkll(u,mb,q2,E_phi,param);
			tperp_0=tperp_bkll(u,0.,q2,E_phi,param);

			Tperpppnf=Tperpmpnf=-4.*ed*C8eff/(u+(1.-u)*q2/mBs/mBs)
			+mBs/2./mb*(eu*tperp_mc*(C2bar+C4bar-C6bar)
			+ed*tperp_mb*(C3bar+C4bar-C6bar-4.*mb/mBs*C5bar)
			+ed*tperp_0*C3bar);
				
			Tperppmnf=Tperpmmnf=0.;

			Tperppnfu=mBs/2./mb*eu*(tperp_mc-tperp_0)*(Cmub[2]-Cmub[1]/6.);
		
			Tperppp=Tperppp0+alphas_muf*4./3./4./pi*(Tperpppf+Tperpppnf); 
			Tperppm=Tperppm0+alphas_muf*4./3./4./pi*(Tperppmf+Tperppmnf);
			Tperpmp=Tperpmp0+alphas_muf*4./3./4./pi*(Tperpmpf+Tperpmpnf);
			Tperpmm=Tperpmm0+alphas_muf*4./3./4./pi*(Tperpmmf+Tperpmmnf);
			Tperppp_bar=Tperppp0+alphas_muf*4./3./4./pi*(Tperpppf_bar+Tperpppnf);
			Tperpmp_bar=Tperpmp0+alphas_muf*4./3./4./pi*(Tperpmpf_bar+Tperpmpnf);

			Tperppu=alphas_muf*4./3./4./pi*Tperppnfu;
	
			int_perppp+=phiphi_perp*Tperppp/n1/lambda_Bp; 
			int_perppm+=phiphi_perp*Tperppm/n1/lambda_Bm;
			int_perpmp+=phiphi_perp*Tperpmp/n1/lambda_Bp;
			int_perpmm+=phiphi_perp*Tperpmm/n1/lambda_Bm;
			int_perppp_bar+=phiphi_perp*Tperppp_bar/n1/lambda_Bp;
			int_perpmp_bar+=phiphi_perp*Tperpmp_bar/n1/lambda_Bp;

			int_perppu+=phiphi_perp*Tperppu/n1/lambda_Bp;


		/* Tpar */		

			phiphi_par=phi_Kstar(u,a1phi_par,a2phi_par);
			tpar_mc=tpar_bkll(u,mc,q2,E_phi,param);
			tpar_mb=tpar_bkll(u,mb,q2,E_phi,param);
			tpar_0=tpar_bkll(u,0.,q2,E_phi,param);
		
			Tparpp0=Tparmp0=0.;
		
			// Tparppf=(C7eff+C7effp)*4.*mBs/E_phi/(1.-u);// commented; not used
			Tparmpf=(C7eff-C7effp)*4.*mBs/E_phi/(1.-u);
			// Tparppf_bar=conj(C7eff+C7effp)*4.*mBs/E_phi/(1.-u);// commented; not used
			Tparmpf_bar=conj(C7eff-C7effp)*4.*mBs/E_phi/(1.-u);

			Tparppnf=Tparmpnf=mBs/mb*(eu*tpar_mc*(C2bar+C4bar-C6bar)
			+ed*tpar_mb*(C3bar+C4bar-C6bar)
			+ed*tpar_0*C3bar);
				
			Tparpnfu=mBs/mb*eu*(tpar_mc-tpar_0)*(Cmub[2]-Cmub[1]/6.);
	
			Tparpu=alphas_muf*4./3./4./pi*Tparpnfu;
	
			Tparpm0=Tparmm0=-eq*4.*mBs/mb*((C3bar+3.*C4bar)+12.*(Cmub[3]+10.*Cmub[5]));
			
			Tparp0u=0.;

			Tparm0u=-eq*4.*mBs/mb*(-(4./3.*Cmub[1]+Cmub[2]));
		
			Tparpmf=Tparmmf=0.;

			h_mc=h_bkll((1.-u)*mBs*mBs+u*q2,mc,mu_b);
			h_mb=h_bkll((1.-u)*mBs*mBs+u*q2,mbpole,mu_b);
			h_0=h_bkll((1.-u)*mBs*mBs+u*q2,0.,mu_b);

			Tparpmnf=Tparmmnf=eq*(8.*C8eff/((1.-u)+u*q2/mBs/mBs)
			+6.*mBs/mb*(h_mc*(C2bar+C4bar+C6bar)
			+h_mb*(C3bar+C4bar+C6bar)
			+h_0*(C3bar+3.*C4bar+3.*C6bar)
			-8./27.*(C3bar-C5bar-15.*C6bar)));
	
			// Tparpp=Tparpp0+alphas_muf*4./3./4./pi*(Tparppf+Tparppnf);// commented; not used
			Tparpm=Tparpm0+alphas_muf*4./3./4./pi*(Tparpmf+Tparpmnf);
			Tparmp=Tparmp0+alphas_muf*4./3./4./pi*(Tparmpf+Tparmpnf);
			Tparmm=Tparmm0+alphas_muf*4./3./4./pi*(Tparmmf+Tparmmnf);
			// Tparpp_bar=Tparpp0+alphas_muf*4./3./4./pi*(Tparppf_bar+Tparppnf);// commented; not used
			Tparmp_bar=Tparmp0+alphas_muf*4./3./4./pi*(Tparmpf_bar+Tparmpnf);

			Tparmnfu=eq*(6.*mBs/mb*(h_mc-h_0)*(Cmub[2]-Cmub[1]/6.));

			Tparpu=Tparp0u+alphas_muf*4./3./4./pi*Tparpnfu;
	
			Tparmu=Tparm0u+alphas_muf*4./3./4./pi*Tparmnfu;


			// int_parpp+=(phiphi_par*Tparpp/lambda_Bp)/n1;// commented out; not used
			// int_parpp_bar+=(phiphi_par*Tparpp_bar/lambda_Bp)/n1;// commented out; not used
			int_parpm+=(phiphi_par*Tparpm/lambda_Bm)/n1;
			
			int_parmp+=(phiphi_par*Tparmp/lambda_Bp)/n1;
			int_parmp_bar+=(phiphi_par*Tparmp_bar/lambda_Bp)/n1;
			int_parmm+=(phiphi_par*Tparmm/lambda_Bm)/n1;
			
			int_parpu+=(phiphi_par*Tparpu/lambda_Bp)/n1;
			
			int_parmu+=(phiphi_par*Tparmu/lambda_Bm)/n1;


			integ3+=phiphi_perp/((1.-u)+u*q2/mBs/mBs)/n1;

			x=(1.-u)*mBs*mBs+u*q2;
			
			h_mc=h_bkll(x,mc,mu_b);
			h_mb=h_bkll(x,mbpole,mu_b);
			h_0=h_bkll(x,0.,mu_b);
			
			FV=3./4.*(h_mc*(C2bar+C4bar+C6bar)+h_mb*(C3bar+C4bar+C6bar)+h_0*(C3bar+3.*C4bar+3.*C6bar)-8./27.*(C3bar-C5bar-15.*C6bar));

			FVu=3./4.*(h_mc-h_0)*(Cmub[2]-Cmub[1]/6.);
		

			integ4+=phiphi_perp/((1.-u)+u*q2/mBs/mBs)*FV/n1;
			integ4u+=phiphi_perp/((1.-u)+u*q2/mBs/mBs)*FVu/n1;
		
			Fperp+=phiphi_perp/((1.-u)+u*q2/mBs/mBs)/3./n1;
			Xperp+=(u<=1.-0.5/mBs)*phiphi_perp/pow((1.-u)+u*q2/mBs/mBs,2.)/3./n1;
		
			integ5+=((3./4.*(1.+pow(2.*u-1.,2.))+a1phi_par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2phi_par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2phi_par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1phi_par*(2.*u-1.)+(a2phi_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1phi_par*u+(a2phi_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FV/n1;

			integ5u+=((3./4.*(1.+pow(2.*u-1.,2.))+a1phi_par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2phi_par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2phi_par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1phi_par*(2.*u-1.)+(a2phi_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1phi_par*u+(a2phi_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FVu/n1;

			if(ie==0||ie==n1sav) n1=n1sav;
		}
		
		/* Tau_perp */		

		double complex Tauperpp=xi_perp*Cperpp+pi*pi/3.*param->f_Bs*f_phi_perp/mBs*Xi_perp*(int_perppp+int_perppm); 
		double complex Tauperpp_bar=xi_perp*Cperpp_bar+pi*pi/3.*param->f_Bs*f_phi_perp/mBs*Xi_perp*(int_perppp_bar+int_perppm);

		double complex Tauperpm=xi_perp*Cperpm+pi*pi/3.*param->f_Bs*f_phi_perp/mBs*Xi_perp*(int_perpmp+int_perpmm); 
		double complex Tauperpm_bar=xi_perp*Cperpm_bar+pi*pi/3.*param->f_Bs*f_phi_perp/mBs*Xi_perp*(int_perpmp_bar+int_perpmm);

		double complex Tauperpu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_perp*Cperpu+pi*pi/3.*param->f_Bs*f_phi_perp/mBs*Xi_perp*int_perppu);
		double complex Tauperpu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_perp*Cperpu+pi*pi/3.*param->f_Bs*f_phi_perp/mBs*Xi_perp*int_perppu);

		
		/* Tau_par */		

		double complex Tauparm=xi_par*Cparm+pi*pi/3.*param->f_Bs*f_phi_par/mBs*Xi_par*(int_parmp+int_parmm);
		double complex Tauparm_bar=xi_par*Cparm_bar+pi*pi/3.*param->f_Bs*f_phi_par/mBs*Xi_par*(int_parmp_bar+int_parmm);

		double complex Tauparu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_par*Cparu+pi*pi/3.*param->f_Bs*f_phi_par/mBs*Xi_par*(int_parpu+int_parmu));
		double complex Tauparu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_par*Cparu+pi*pi/3.*param->f_Bs*f_phi_par/mBs*Xi_par*(int_parpu+int_parmu));
				
				
		double complex DeltaTauperpWA=-eq*4.*pi*pi/3.*param->f_Bs*f_phi_perp/mb/mBs*(Cmub[3]+4./3.*Cmub[4]+4.*Cmub[5]+16./3.*Cmub[6])*integ3
		+eq*2.*pi*pi/3.*param->f_Bs*f_phi_par/mb/mBs*mphi/(1.-q2/mBs/mBs)/lambda_Bp*(Cmub[3]+4./3.*Cmub[4]+16.*Cmub[5]+64./3.*Cmub[6]+12.*(Cmub[3]+10.*Cmub[5]));

		double complex DeltaTauperpuWA=-eq*2.*pi*pi/3.*param->f_Bs*f_phi_par/mb/mBs*mphi/(1.-q2/mBs/mBs)/lambda_Bp*(4./3.*Cmub[1]+Cmub[2]);


		double rho=0.;
		double phi=0.;
		Xperp=Fperp+(1.+rho*(cos(phi)+I*sin(phi)))*Xperp;


		double complex DeltaTauperpHSA=eq*4./3.*alphas_muf/4./pi*pi*pi*param->f_Bs/3./mb/mBs*(12.*C8eff*mb/mBs*f_phi_perp*Xperp
		+8.*f_phi_perp*integ4-4.*mphi*f_phi_par/(1.-q2/mBs/mBs)/lambda_Bp*integ5);
		
		double complex DeltaTauperpuHSA=eq*4./3.*alphas_muf/4./pi*pi*pi*param->f_Bs/3./mb/mBs*
		(8.*f_phi_perp*integ4u-4.*mphi*f_phi_par/(1.-q2/mBs/mBs)/lambda_Bp*integ5u);

		Tauperpp+=DeltaTauperpWA+DeltaTauperpHSA;
		Tauperpm+=DeltaTauperpWA+DeltaTauperpHSA;
		Tauperpp_bar+=DeltaTauperpWA+DeltaTauperpHSA;
		Tauperpm_bar+=DeltaTauperpWA+DeltaTauperpHSA;


		/*******Using Tauperpu**********/
		Tauperpu+=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(DeltaTauperpuWA+DeltaTauperpuHSA);
		Tauperpu_bar+=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(DeltaTauperpuWA+DeltaTauperpuHSA);

		double lambda=pow(mBs,4.)+pow(mphi,4.)+q2*q2-2.*(mBs*mBs*mphi*mphi+mphi*mphi*q2+mBs*mBs*q2);
	
		double complex N=param->Vtb*conj(param->Vts)*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)/pow(mBs,3.)*q2*sqrt(lambda)*beta_l);
		
		double complex Nbar=conj(param->Vtb)*param->Vts*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)/pow(mBs,3.)*q2*sqrt(lambda)*beta_l);
		
		int n9=10;
		double integ9=0.;
		for(ie=1;ie<=n9-1;ie++)
		{
			u=(double)ie/n9;	
			integ9+=phi_Kstar(u,a1phi_par,a2phi_par)/(1.-u)/n9;
		}
		double Delta_par=1.+alphas_mub*4./3./4./pi*(-2.+2.*L)-alphas_mub*4./3./4./pi*2.*q2/E_phi/E_phi*pi*pi*param->f_Bs*f_phi_par/lambda_Bp/3./param->m_Bs/(E_phi/param->m_phi)/xi_par*integ9;
	 
		
		/* A and A_bar choosen such that we have agreement with the convention of 1502.05509 & 22210.11995 */
		ALperp_bar=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)-(C10+C10p))*V/(mBs+mphi)+2.*mb/q2*(Tauperpp+Tauperpu));
		ARperp_bar=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)+(C10+C10p))*V/(mBs+mphi)+2.*mb/q2*(Tauperpp+Tauperpu));

		ALpar_bar=-N*sqrt(2.)*(mBs*mBs-mphi*mphi)*(((C9-C9p)-(C10-C10p))*(2.*E_phi/(mBs+mphi)*xi_perp)/(mBs-mphi)+4.*mb/mBs*E_phi/q2*(Tauperpm+Tauperpu));
		ARpar_bar=-N*sqrt(2.)*(mBs*mBs-mphi*mphi)*(((C9-C9p)+(C10-C10p))*(2.*E_phi/(mBs+mphi)*xi_perp)/(mBs-mphi)+4.*mb/mBs*E_phi/q2*(Tauperpm+Tauperpu));

		AL0_bar=-N/2./mphi/sqrt(q2)*(((C9-C9p)-(C10-C10p))*((mBs*mBs-mphi*mphi-q2)*(mBs+mphi)*(2.*E_phi/(mBs+mphi)*xi_perp)-lambda*(mBs/(mBs-mphi)*(xi_perp-xi_par))/(mBs+mphi))+2.*mb*(2.*E_phi/mBs*(mBs*mBs+3.*mphi*mphi-q2)*(Tauperpm+Tauperpu)-lambda/(mBs*mBs-mphi*mphi)*(Tauperpm+Tauparm+Tauperpu+Tauparu)));
		AR0_bar=-N/2./mphi/sqrt(q2)*(((C9-C9p)+(C10-C10p))*((mBs*mBs-mphi*mphi-q2)*(mBs+mphi)*(2.*E_phi/(mBs+mphi)*xi_perp)-lambda*(mBs/(mBs-mphi)*(xi_perp-xi_par))/(mBs+mphi))+2.*mb*(2.*E_phi/mBs*(mBs*mBs+3.*mphi*mphi-q2)*(Tauperpm+Tauperpu)-lambda/(mBs*mBs-mphi*mphi)*(Tauperpm+Tauparm+Tauperpu+Tauparu)));

		At_bar=N/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/ml*(CQ2-CQ2p)/(mb_mub+param->mass_s))*(E_phi/mphi*xi_par/Delta_par);

		AS_bar=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(mb_mub+param->mass_s)*(E_phi/mphi*xi_par/Delta_par);



		ALperp=Nbar*sqrt(2.)*sqrt(lambda)*(conj((C9+C9p)-(C10+C10p))*V/(mBs+mphi)+2.*mb/q2*(Tauperpp_bar+Tauperpu_bar));
		ARperp=Nbar*sqrt(2.)*sqrt(lambda)*(conj((C9+C9p)+(C10+C10p))*V/(mBs+mphi)+2.*mb/q2*(Tauperpp_bar+Tauperpu_bar));

		ALpar=-Nbar*sqrt(2.)*(mBs*mBs-mphi*mphi)*(conj((C9-C9p)-(C10-C10p))*(2.*E_phi/(mBs+mphi)*xi_perp)/(mBs-mphi)+4.*mb/mBs*E_phi/q2*(Tauperpm_bar+Tauperpu_bar));
		ARpar=-Nbar*sqrt(2.)*(mBs*mBs-mphi*mphi)*(conj((C9-C9p)+(C10-C10p))*(2.*E_phi/(mBs+mphi)*xi_perp)/(mBs-mphi)+4.*mb/mBs*E_phi/q2*(Tauperpm_bar+Tauperpu_bar));

		AL0=-Nbar/2./mphi/sqrt(q2)*(conj((C9-C9p)-(C10-C10p))*((mBs*mBs-mphi*mphi-q2)*(mBs+mphi)*(2.*E_phi/(mBs+mphi)*xi_perp)-lambda*(mBs/(mBs-mphi)*(xi_perp-xi_par))/(mBs+mphi))+2.*mb*(2.*E_phi/mBs*(mBs*mBs+3.*mphi*mphi-q2)*(Tauperpm_bar+Tauperpu_bar)-lambda/(mBs*mBs-mphi*mphi)*(Tauperpm_bar+Tauparm_bar+Tauperpu_bar+Tauparu_bar)));
		AR0=-Nbar/2./mphi/sqrt(q2)*(conj((C9-C9p)+(C10-C10p))*((mBs*mBs-mphi*mphi-q2)*(mBs+mphi)*(2.*E_phi/(mBs+mphi)*xi_perp)-lambda*(mBs/(mBs-mphi)*(xi_perp-xi_par))/(mBs+mphi))+2.*mb*(2.*E_phi/mBs*(mBs*mBs+3.*mphi*mphi-q2)*(Tauperpm_bar+Tauperpu_bar)-lambda/(mBs*mBs-mphi*mphi)*(Tauperpm_bar+Tauparm_bar+Tauperpu_bar+Tauparu_bar)));

		At=Nbar/sqrt(q2)*sqrt(lambda)*(2.*conj(C10-C10p)+q2/ml*conj(CQ2-CQ2p)/(mb_mub+param->mass_s))*(E_phi/mphi*xi_par/Delta_par);
		AS=-2.*Nbar*sqrt(lambda)*conj(CQ1-CQ1p)/(mb_mub+param->mass_s)*(E_phi/mphi*xi_par/Delta_par);



		/* hadronic uncertainties */
		ALperp*=1.+param->Bstophilow_ALperp_err_noq2+q2/6.*param->Bstophilow_ALperp_err_q2;
		ARperp*=1.+param->Bstophilow_ARperp_err_noq2+q2/6.*param->Bstophilow_ARperp_err_q2;
		ALpar*=1.+param->Bstophilow_ALpar_err_noq2+q2/6.*param->Bstophilow_ALpar_err_q2;
		ARpar*=1.+param->Bstophilow_ARpar_err_noq2+q2/6.*param->Bstophilow_ARpar_err_q2;
		AL0*=1.+param->Bstophilow_AL0_err_noq2+q2/6.*param->Bstophilow_AL0_err_q2;
		AR0*=1.+param->Bstophilow_AR0_err_noq2+q2/6.*param->Bstophilow_AR0_err_q2;

		ALperp_bar*=1.+param->Bstophilow_ALperp_err_noq2+q2/6.*param->Bstophilow_ALperp_err_q2;
		ARperp_bar*=1.+param->Bstophilow_ARperp_err_noq2+q2/6.*param->Bstophilow_ARperp_err_q2;
		ALpar_bar*=1.+param->Bstophilow_ALpar_err_noq2+q2/6.*param->Bstophilow_ALpar_err_q2;
		ARpar_bar*=1.+param->Bstophilow_ARpar_err_noq2+q2/6.*param->Bstophilow_ARpar_err_q2;
		AL0_bar*=1.+param->Bstophilow_AL0_err_noq2+q2/6.*param->Bstophilow_AL0_err_q2;
		AR0_bar*=1.+param->Bstophilow_AR0_err_noq2+q2/6.*param->Bstophilow_AR0_err_q2;
	}
	
	if(q2>LOW_q2)
	{
		// double mb=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);
		/* Loop function are calculated with Pole mass [S04] which should be replaced by PS mass [BF01] */
		double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
		double mb_PS=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */ //mb=>4.64
		double sh=q2/mb_PS/mb_PS;

		double shat=q2/mBs/mBs;

		double lambda_hat=1.+shat*shat+pow(mphi/mBs,4.)-2.*(shat+shat*mphi*mphi/mBs/mBs+mphi*mphi/mBs/mBs);


		double complex A = A_Seidel(q2, mb_PS, mu_b);

		double complex B = B_Seidel(q2, mb_PS, mu_b);

		double complex C = C_Seidel(q2, mu_b);

		double kappa=1.-2.*alphas_mub/3./pi*log(mu_b/param->mass_b);// MSbar mass should be used for mb

		double complex C9eff=Cmub[9]
		+h_bkll(q2,0.,mu_b)*(4./3.*Cmub[1]+Cmub[2]+11./2.*Cmub[3]-2./3.*Cmub[4]+52.*Cmub[5]-32./3.*Cmub[6])
		-1./2.*h_bkll(q2,mb_PS,mu_b)*(7.*Cmub[3]+4./3.*Cmub[4]+76.*Cmub[5]+64./3.*Cmub[6])
		+4./3.*(Cmub[3]+16./3.*Cmub[5]+16./9.*Cmub[6])
		+alphas_mub/4./pi*(Cmub[1]*(B+4.*C)-3.*Cmub[2]*(2.*B-C)-Cmub[8]*F89_bsll(sh))
		+8.*param->mass_c*param->mass_c/q2*((4./9.*Cmub[1]+1./3.*Cmub[2])*(1.+param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts))+2.*Cmub[3]+20.*Cmub[5]);

		double complex C9eff_bar=Cmub[9]
		+h_bkll(q2,0.,mu_b)*(4./3.*Cmub[1]+Cmub[2]+11./2.*Cmub[3]-2./3.*Cmub[4]+52.*Cmub[5]-32./3.*Cmub[6])
		-1./2.*h_bkll(q2,mb_PS,mu_b)*(7.*Cmub[3]+4./3.*Cmub[4]+76.*Cmub[5]+64./3.*Cmub[6])
		+4./3.*(Cmub[3]+16./3.*Cmub[5]+16./9.*Cmub[6])
		+alphas_mub/4./pi*(Cmub[1]*(B+4.*C)-3.*Cmub[2]*(2.*B-C)-Cmub[8]*F89_bsll(sh))
		+8.*param->mass_c*param->mass_c/q2*((4./9.*Cmub[1]+1./3.*Cmub[2])*(1.+param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb))+2.*Cmub[3]+20.*Cmub[5]);


		double complex C7eff=Cmub[7]
		+alphas_mub/4./pi*((Cmub[1]-6.*Cmub[2])*A-Cmub[8]*F87_bsll(sh,log(mu_b/mb_PS)));

		double complex C7eff_bar=conj(Cmub[7])
		+alphas_mub/4./pi*((Cmub[1]-6.*Cmub[2])*A-Cmub[8]*F87_bsll(sh,log(mu_b/mb_PS)));

	    double complex C7effp=Cpb[7];
        double complex C9p=Cpb[9];
        double complex C10p=Cpb[10];

		double complex N=param->Vtb*conj(param->Vts)*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em*beta_l/3./1024./pow(pi,5.)*mBs*shat*sqrt(lambda_hat));// Missing beta_l corrected

		double complex Nbar=conj(param->Vtb)*param->Vts*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em*beta_l/3./1024./pow(pi,5.)*mBs*shat*sqrt(lambda_hat));// Missing beta_l corrected

		double complex f_perp=N*mBs*sqrt(2.*lambda_hat)/(1.+mphi/mBs)*V;

		double complex f_par=N*mBs*sqrt(2.)*(1.+mphi/mBs)*A1;

		double complex f_0=N*mBs*((1.-shat-mphi*mphi/mBs/mBs)*pow(1.+mphi/mBs,2.)*A1-lambda_hat*A2)/(2.*mphi/mBs*(1.+mphi/mBs)*sqrt(shat));

		double complex f_perp_bar=Nbar*mBs*sqrt(2.*lambda_hat)/(1.+mphi/mBs)*V;

		double complex f_par_bar=Nbar*mBs*sqrt(2.)*(1.+mphi/mBs)*A1;

		double complex f_0_bar=Nbar*mBs*((1.-shat-mphi*mphi/mBs/mBs)*pow(1.+mphi/mBs,2.)*A1-lambda_hat*A2)/(2.*mphi/mBs*(1.+mphi/mBs)*sqrt(shat));

		/* A and A_bar choosen such that we have agreement with the convention of 1502.05509 & 22210.11995 */
		ALperp_bar_high=(((C9eff+C9p)-(Cmub[10]+C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff+C7effp))*f_perp;
		ARperp_bar_high=(((C9eff+C9p)+(Cmub[10]+C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff+C7effp))*f_perp;

		ALpar_bar_high=-(((C9eff-C9p)-(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff-C7effp))*f_par;
		ARpar_bar_high=-(((C9eff-C9p)+(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff-C7effp))*f_par;

		AL0_bar_high=-(((C9eff-C9p)-(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff-C7effp))*f_0;
		AR0_bar_high=-(((C9eff-C9p)+(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff-C7effp))*f_0;


		/************************************************************************************************/
		ALperp_high=(((C9eff_bar+conj(C9p))-conj(Cmub[10]+C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar+conj(C7effp)))*f_perp_bar;
		ARperp_high=(((C9eff_bar+conj(C9p))+conj(Cmub[10]+C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar+conj(C7effp)))*f_perp_bar;

		ALpar_high=-(((C9eff_bar-conj(C9p))-conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar-conj(C7effp)))*f_par_bar;
		ARpar_high=-(((C9eff_bar-conj(C9p))+conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar-conj(C7effp)))*f_par_bar;

		AL0_high=-(((C9eff_bar-conj(C9p))-conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar-conj(C7effp)))*f_0_bar;
		AR0_high=-(((C9eff_bar-conj(C9p))+conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar-conj(C7effp)))*f_0_bar;


		double complex C10=Cmub[10];
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];

		double lambda=pow(mBs,4.)+pow(mphi,4.)+q2*q2-2.*(mBs*mBs*mphi*mphi+mphi*mphi*q2+mBs*mBs*q2);

		At_bar_high=N/sqrt(q2)*sqrt(lambda)*( 2.*conj(C10-C10p)+q2/ml*conj(CQ2-CQ2p)/(param->mass_b+param->mass_s))*A0;
		AS_bar_high=-2.*N*sqrt(lambda)*conj(CQ1-CQ1p)/(param->mass_b+param->mass_s)*A0;

		At_high=Nbar/sqrt(q2)*sqrt(lambda)*( 2.*(C10-C10p)+q2/ml*(CQ2-CQ2p)/(param->mass_b+param->mass_s))*A0;
		AS_high=-2.*Nbar*sqrt(lambda)*(CQ1-CQ1p)/(param->mass_b+param->mass_s)*A0;

		/* hadronic uncertainties */
		ALperp_high*=1.+param->Bstophihigh_ALperp_err;
		ARperp_high*=1.+param->Bstophihigh_ARperp_err;
		ALpar_high*=1.+param->Bstophihigh_ALpar_err;
		ARpar_high*=1.+param->Bstophihigh_ARpar_err;
		AL0_high*=1.+param->Bstophihigh_AL0_err;
		AR0_high*=1.+param->Bstophihigh_AR0_err;
		At_high*=1.+param->Bstophihigh_At_err;
		AS_high*=1.+param->Bstophihigh_AS_err;
		ALperp_bar_high*=1.+param->Bstophihigh_ALperp_err;
		ARperp_bar_high*=1.+param->Bstophihigh_ARperp_err;
		ALpar_bar_high*=1.+param->Bstophihigh_ALpar_err;
		ARpar_bar_high*=1.+param->Bstophihigh_ARpar_err;
		AL0_bar_high*=1.+param->Bstophihigh_AL0_err;
		AR0_bar_high*=1.+param->Bstophihigh_AR0_err;
		At_bar_high*=1.+param->Bstophihigh_At_err;
		AS_bar_high*=1.+param->Bstophihigh_AS_err;


		if(q2>HIGH_q2)
		{
			ALperp=ALperp_high;
			ARperp=ARperp_high;
			ALpar=ALpar_high;
			ARpar=ARpar_high;
			AL0=AL0_high;
			AR0=AR0_high;
			At=At_high;
			AS=AS_high;

			ALperp_bar=ALperp_bar_high;
			ARperp_bar=ARperp_bar_high;
			ALpar_bar=ALpar_bar_high;
			ARpar_bar=ARpar_bar_high;
			AL0_bar=AL0_bar_high;
			AR0_bar=AR0_bar_high;
			At_bar=At_bar_high;
			AS_bar=AS_bar_high;
		}
		else
		{
			ALperp=ALperp*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ALperp_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ARperp=ARperp*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ARperp_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ALpar=ALpar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ALpar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ARpar=ARpar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ARpar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AL0=AL0*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AL0_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AR0=AR0*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AR0_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			At=At*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+At_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AS=AS*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AS_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);

			ALperp_bar=ALperp_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ALperp_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ARperp_bar=ARperp_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ARperp_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ALpar_bar=ALpar_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ALpar_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ARpar_bar=ARpar_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ARpar_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AL0_bar=AL0_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AL0_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AR0_bar=AR0_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AR0_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			At_bar=At_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+At_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AS_bar=AS_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AS_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
		}
	}

	double A02=AL0*conj(AL0)+AR0*conj(AR0);
	double Apar2=ALpar*conj(ALpar)+ARpar*conj(ARpar);
	double Aperp2=ALperp*conj(ALperp)+ARperp*conj(ARperp);
	
	double A02_bar=AL0_bar*conj(AL0_bar)+AR0_bar*conj(AR0_bar);
	double Apar2_bar=ALpar_bar*conj(ALpar_bar)+ARpar_bar*conj(ARpar_bar);
	double Aperp2_bar=ALperp_bar*conj(ALperp_bar)+ARperp_bar*conj(ARperp_bar);
	
	
	double J1s=0.25*(2.+beta_l*beta_l)*(Aperp2 + Apar2) + 4.*ml*ml/q2*creal(ALperp*conj(ARperp)+ALpar*conj(ARpar));

	double J1c=A02 + 4.*ml*ml/q2*(At*conj(At)+2.*creal(AL0*conj(AR0)))+beta_l*beta_l*AS*conj(AS);

	double J2s=0.25*beta_l*beta_l*(Aperp2+Apar2);

	double J2c=-beta_l*beta_l*A02;
	
	double J3=0.5*beta_l*beta_l*(Aperp2-Apar2);
	
	double J4=1./sqrt(2.)*beta_l*beta_l*creal(AL0*conj(ALpar)+AR0*conj(ARpar));
	
	double J5=sqrt(2.)*beta_l*(creal(AL0*conj(ALperp)-AR0*conj(ARperp))-ml/sqrt(q2)*(creal(ALpar*conj(AS)+ARpar*conj(AS))));
	
	double J6s=2.*beta_l*creal(ALpar*conj(ALperp)-ARpar*conj(ARperp));

	double J6c=4.*beta_l*ml/sqrt(q2)*creal(AL0*conj(AS)+AR0*conj(AS));

	double J7=sqrt(2.)*beta_l*(cimag(AL0*conj(ALpar)-AR0*conj(ARpar))+ml/sqrt(q2)*(cimag(ALperp*conj(AS)+ARperp*conj(AS))));

	double J8=1./sqrt(2.)*beta_l*beta_l*cimag(AL0*conj(ALperp)+AR0*conj(ARperp));
	
	double J9=beta_l*beta_l*cimag(conj(ALpar)*ALperp+conj(ARpar)*ARperp);

	double J1s_bar=0.25*(2.+beta_l*beta_l)*(Aperp2_bar + Apar2_bar) + 4.*ml*ml/q2*creal(ALperp_bar*conj(ARperp_bar)+ALpar_bar*conj(ARpar_bar));

	double J1c_bar=A02_bar + 4.*ml*ml/q2*(At_bar*conj(At_bar)+2.*creal(AL0_bar*conj(AR0_bar)))+beta_l*beta_l*AS_bar*conj(AS_bar);

	double J2s_bar=0.25*beta_l*beta_l*(Aperp2_bar+Apar2_bar);

	double J2c_bar=-beta_l*beta_l*A02_bar;
	
	double J3_bar=0.5*beta_l*beta_l*(Aperp2_bar-Apar2_bar);
	
	double J4_bar=1./sqrt(2.)*beta_l*beta_l*creal(AL0_bar*conj(ALpar_bar)+AR0_bar*conj(ARpar_bar));
	
	double J5_bar=sqrt(2.)*beta_l*(creal(AL0_bar*conj(ALperp_bar)-AR0_bar*conj(ARperp_bar))-ml/sqrt(q2)*(creal(ALpar_bar*conj(AS_bar)+ARpar_bar*conj(AS_bar))));
	
	double J6s_bar=2.*beta_l*creal(ALpar_bar*conj(ALperp_bar)-ARpar_bar*conj(ARperp_bar));

	double J6c_bar=4.*beta_l*ml/sqrt(q2)*creal(AL0_bar*conj(AS_bar)+AR0_bar*conj(AS_bar));

	double J7_bar=sqrt(2.)*beta_l*(cimag(AL0_bar*conj(ALpar_bar)-AR0_bar*conj(ARpar_bar))+ml/sqrt(q2)*(cimag(ALperp_bar*conj(AS_bar)+ARperp_bar*conj(AS_bar))));

	double J8_bar=1./sqrt(2.)*beta_l*beta_l*cimag(AL0_bar*conj(ALperp_bar)+AR0_bar*conj(ARperp_bar));
	
	double J9_bar=beta_l*beta_l*cimag(conj(ALpar_bar)*ALperp_bar+conj(ARpar_bar)*ARperp_bar);
	
	
	double complex exppIphi=(1.+0.04*I);
	double complex expmIphi=(1.-0.04*I);
	
	double ys=param->ys_Bs;
	
	double AL0_tilde = AL0_bar;
	double ALpar_tilde = ALpar_bar;
	double AR0_tilde = AR0_bar;
	double ARpar_tilde = ARpar_bar;
	double At_tilde = At_bar;


	double ALperp_tilde = -ALperp_bar;
	double ARperp_tilde = -ARperp_bar;
	double AS_tilde = -AS_bar;

	double complex sh1s_A = exppIphi*(ALperp_tilde*conj(ALperp)+ALpar_tilde*conj(ALpar)+ARperp_tilde*conj(ARperp)+ARpar_tilde*conj(ARpar));
	double complex sh1s_B = exppIphi*(ALperp_tilde*conj(ARperp)+ALpar_tilde*conj(ARpar));
	double complex sh1s_C = expmIphi*(ALperp*conj(ARperp_tilde)+ALpar*conj(ARpar_tilde));
	double h1s = (2.+beta_l*beta_l)/2.*creal(sh1s_A) + 4.*ml*ml/q2*creal(sh1s_B + sh1s_C);
	double s1s = (2.+beta_l*beta_l)/2.*cimag(sh1s_A) + 4.*ml*ml/q2*cimag(sh1s_B - sh1s_C);

	double complex sh1c_A = exppIphi*(AL0_tilde*conj(AL0)+AR0_tilde*conj(AR0));
	double complex sh1c_B = exppIphi*At_tilde*conj(At);
	double complex sh1c_C = exppIphi*AL0_tilde*conj(AR0);
	double complex sh1c_D = expmIphi*AL0*conj(AR0_tilde);
	double complex sh1c_E = exppIphi*AS_tilde*conj(AS);
	double h1c = 2.*creal(sh1c_A)+8.*ml*ml/q2*(creal(sh1c_B)+creal(sh1c_C + sh1c_D))+2.*beta_l*beta_l*creal(sh1c_E);
	double s1c = 2.*cimag(sh1c_A)+8.*ml*ml/q2*(cimag(sh1c_B)+cimag(sh1c_C - sh1c_D))+2.*beta_l*beta_l*cimag(sh1c_E);

	double h2s = beta_l*beta_l/2.*creal(sh1s_A);
	double s2s = beta_l*beta_l/2.*cimag(sh1s_A);

	double h2c = -2.*beta_l*beta_l*creal(sh1c_A);
	double s2c = -2.*beta_l*beta_l*cimag(sh1c_A);

	double complex sh3_A = exppIphi*(ALperp_tilde*conj(ALperp)-ALpar_tilde*conj(ALpar)+ARperp_tilde*conj(ARperp)-ARpar_tilde*conj(ARpar));
	double h3 = beta_l*beta_l*creal(sh3_A);
	double s3 = beta_l*beta_l*cimag(sh3_A);


	double complex sh4_A = exppIphi*(AL0_tilde*conj(ALpar)+AR0_tilde*conj(ARpar));
	double complex sh4_B = expmIphi*(AL0*conj(ALpar_tilde)+AR0*conj(ARpar_tilde));
	double h4 = 1./(sqrt(2.))*beta_l*beta_l*creal(sh4_A+sh4_B);
	double s4 = 1./(sqrt(2.))*beta_l*beta_l*cimag(sh4_A-sh4_B);

	double complex sh5_A = exppIphi*(AL0_tilde*conj(ALperp)-AR0_tilde*conj(ARperp));
	double complex sh5_B = expmIphi*(AL0*conj(ALperp_tilde)-AR0*conj(ARpar_tilde));
	double complex sh5_C = exppIphi*(ALpar_tilde*conj(AS)+ARpar_tilde*conj(AS));
	double complex sh5_D = expmIphi*(ALpar*conj(AS_tilde)+ARpar*conj(AS_tilde));
	double h5 = sqrt(2.)*beta_l*( creal(sh5_A + sh5_B) -ml/sqrt(q2)*creal(sh5_C + sh5_D) );
	double s5 = sqrt(2.)*beta_l*( cimag(sh5_A - sh5_B) -ml/sqrt(q2)*creal(sh5_C - sh5_D) );

	double complex sh6s_A = exppIphi*(ALpar_tilde*conj(ALperp)-ARpar_tilde*conj(ARperp));
	double complex sh6s_B = expmIphi*(ALpar*conj(ALperp_tilde)-ARpar*conj(ARperp_tilde));
	double h6s = 2.*beta_l*creal(sh6s_A + sh6s_B);
	double s6s = 2.*beta_l*cimag(sh6s_A - sh6s_B);

	double complex sh6c_A = exppIphi*(AL0_tilde*conj(AS)+AR0_tilde*conj(AS));
	double complex sh6c_B = expmIphi*(AL0*conj(AS_tilde)+AR0*conj(AS_tilde));
	double h6c = 4.*beta_l*ml/sqrt(q2)*creal(sh6c_A + sh6c_B);
	double s6c = 4.*beta_l*ml/sqrt(q2)*cimag(sh6c_A - sh6c_B);


	double complex sh7_A = exppIphi*(AL0_tilde*conj(ALpar)-AR0_tilde*conj(ARpar));
	double complex sh7_B = expmIphi*(AL0*conj(ALpar_tilde)-AR0*conj(ARpar_tilde));
	double complex sh7_C = exppIphi*(ALperp_tilde*conj(AS)+ARperp_tilde*conj(AS));
	double complex sh7_D = expmIphi*(ALperp*conj(AS_tilde)+ARperp*conj(AS_tilde));
	double h7 = sqrt(2.)*beta_l*(cimag(sh7_A + sh7_B)+ml/(sqrt(q2))*cimag(sh7_C + sh7_D));
	double s7 = -sqrt(2.)*beta_l*(creal(sh7_A - sh7_B)+ml/(sqrt(q2))*creal(sh7_C - sh7_D));

	double complex sh8_A = exppIphi*(AL0_tilde*conj(ALperp)+AR0_tilde*conj(ARperp));
	double complex sh8_B = expmIphi*(AL0*conj(ALperp_tilde)+AR0*conj(ARperp_tilde));
	double h8 = 1./sqrt(q2)*beta_l*beta_l*cimag(sh8_A + sh8_B);
	double s8 = -1./sqrt(q2)*beta_l*beta_l*creal(sh8_A - sh8_B);

	double complex sh9_A = exppIphi*(ALpar_tilde*conj(ALperp)+ARpar_tilde*conj(ARperp));
	double complex sh9_B = expmIphi*(ALpar*conj(ALperp_tilde)+ARpar*conj(ARperp_tilde));
	double h9 = -beta_l*beta_l*cimag(sh9_A + sh9_B);
	double s9 = beta_l*beta_l*creal(sh9_A - sh9_B);


	double dGamma_bsphill_dq2=3./4.*(2.*(J1s-ys/2.*h1s)+(J1c-ys/2.*h1c)-(2.*(J2s-ys/2.*h2s)+(J2c-ys/2.*h2c))/3.);

	double dGamma_bar_bsphill_dq2=3./4.*(2.*(J1s_bar-ys/2.*h1s)+(J1c_bar-ys/2.*h1c)-(2.*(J2s_bar-ys/2.*h2s)+(J2c_bar-ys/2.*h2c))/3.);
	

	double FL[3],AT2[3],S3[3],S4[3],P4prime[3],S7[3];
	
	FL[0]=-((J2c+J2c_bar)-ys*h2c)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	FL[1]=-(J2c+J2c_bar-ys*h2c);
	FL[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	AT2[0]=(J3+J3_bar-ys*h3)/2./(J2s+J2s_bar-ys*h2s);
	AT2[1]=J3+J3_bar-ys*h3;
	AT2[2]=2.*(J2s+J2s_bar-ys*h2s);
	/*********************************************/

	S3[0]=(J3+J3_bar-ys*h3)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	S3[1]=J3+J3_bar-ys*h3;
	S3[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	S4[0]=((J4+J4_bar)-ys*h4)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	S4[1]=(J4+J4_bar-ys*h4);
	S4[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	P4prime[0]=(J4+J4_bar-ys*h4)/sqrt(fabs(-(J2c+J2c_bar-ys*h2c)*(J2s+J2s_bar-ys*h2s)));
	P4prime[1]=J4+J4_bar-ys*h4;
	P4prime[2]=-(J2c+J2c_bar-ys*h2c);

	S7[0]=((J7+J7_bar)-ys*h7)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	S7[1]=(J7+J7_bar-ys*h7);
	S7[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	double P6prime[3];
	P6prime[0]=-(J7+J7_bar-ys*h7)/2./sqrt(fabs(-(J2c+J2c_bar-ys*h2c)*(J2s+J2s_bar-ys*h2s)));
	P6prime[1]=-(J7+J7_bar-ys*h7)/2.;
	P6prime[2]=-(J2c+J2c_bar-ys*h2c);

	double ATReCP[3], ATImCP[3], A5[3], AFBCP[3], A8[3], A9[3], P2CP[3], P3CP[3], P5primeCP[3], P8primeCP[3], A6c[3];
	/* From 1106.3283 which includes the beta_l function; ATRe = beta_l*2*P2 */
	ATReCP[0]=beta_l*(J6s-J6s_bar-ys*h6s)/(4.*(J2s+J2s_bar-ys*h2s));
	ATReCP[1]=beta_l*(J6s-J6s_bar-ys*h6s);
	ATReCP[2]=4.*(J2s+J2s_bar-ys*h2s);

	ATImCP[0]=(J9-J9_bar-ys*h9)/(2.*(J2s+J2s_bar-ys*h2s));
	ATImCP[1]=(J9-J9_bar-ys*h9);
	ATImCP[2]=2.*(J2s+J2s_bar-ys*h2s);

	A5[0]=(J5-J5_bar-ys*h5)/dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;
	A5[1]=(J5-J5_bar-ys*h5);
	A5[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	int scalar_switch = 0;// 0 if we want to put scalar to zero in AFB, else 1
	AFBCP[0]=-3./8.*(2.*(J6s-J6s_bar-ys*h6s)+scalar_switch*(J6c-J6c_bar-ys*h6c))/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	AFBCP[1]=-3./8.*(2.*(J6s-J6s_bar-ys*h6s)+scalar_switch*(J6c-J6c_bar-ys*h6c));
	AFBCP[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	A8[0]=(J8-J8_bar-ys*h8)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	A8[1]=(J8-J8_bar-ys*h8);
	A8[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	A9[0]=(J9-J9_bar-ys*h9)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	A9[1]=(J9-J9_bar-ys*h9);
	A9[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	P2CP[0]=(J6s-J6s_bar-ys*h6s)/8./(J2s+J2s_bar-ys*h2s);// Definition from 1303.5794 without beta_l and not from 1202.4266
	P2CP[1]=(J6s-J6s_bar-ys*h6s);
	P2CP[2]=8.*(J2s+J2s_bar-ys*h2s);

	P3CP[0]=-(J9-J9_bar-ys*h9)/(4.*(J2s+J2s_bar-ys*h2s));
	P3CP[1]=-(J9-J9_bar-ys*h9);
	P3CP[2]=4.*(J2s+J2s_bar-ys*h2s);

	P5primeCP[0]=(J5-J5_bar-ys*h5)/2./sqrt(fabs(-(J2c+J2c_bar-ys*h2c)*(J2s+J2s_bar-ys*h2s)));
	P5primeCP[1]=(J5-J5_bar-ys*h5)/2.;
	P5primeCP[2]=-(J2c+J2c_bar-ys*h2c);

	P8primeCP[0]=-(J8-J8_bar-ys*h8)/sqrt(fabs(-(J2c+J2c_bar-ys*h2c)*(J2s+J2s_bar-ys*h2s)));
	P8primeCP[1]=-(J8-J8_bar-ys*h8);
	P8primeCP[2]=-(J2c+J2c_bar-ys*h2c);

	A6c[0]=-(J6c-J6c_bar-ys*h6c)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	A6c[1]=-(J6c-J6c_bar-ys*h6c);
	A6c[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	/* testing */
	double Q8m[3], Q8p[3], Q9[3], S2s[3];
	Q8m[0]=s8/sqrt(fabs( -2.*(J2c+J2c_bar-ys*h2c)*( 2.*(J2s+J2s_bar-ys*h2s) - (J3+J3_bar-ys*h3) ) ));
	Q8m[1]=s8;
	Q8m[2]=sqrt(fabs( -2.*(J2c+J2c_bar-ys*h2c)*( 2.*(J2s+J2s_bar-ys*h2s) - (J3+J3_bar-ys*h3) ) ));

	Q8p[0]=s8/sqrt(fabs( -2.*(J2c+J2c_bar-ys*h2c)*( 2.*(J2s+J2s_bar-ys*h2s) + (J3+J3_bar-ys*h3) ) ));
	Q8p[1]=s8;
	Q8p[2]=sqrt(fabs( -2.*(J2c+J2c_bar-ys*h2c)*( 2.*(J2s+J2s_bar-ys*h2s) + (J3+J3_bar-ys*h3) ) ));

	Q9[0]=s9/(2.*(J2s+J2s_bar-ys*h2s));
	Q9[1]=s9;
	Q9[2]=2.*(J2s+J2s_bar-ys*h2s);

	S2s[0]=(J2s+J2s_bar-ys*h2s)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	S2s[1]=(J2s+J2s_bar-ys*h2s);
	S2s[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;
	/*********************************************/

	for(je=0;je<=Nobs_Bsphill;je++) for(ie=0;ie<=2;ie++) obs[je][ie]=0.;
	
	for(ie=0;ie<=2;ie++)
	{
		obs[1][ie]=FL[ie];
		obs[2][ie]=AT2[ie];// = P1
		obs[3][ie]=S3[ie];
		obs[4][ie]=S4[ie];
		obs[5][ie]=P4prime[ie];
		obs[6][ie]=S7[ie];
		obs[7][ie]=P6prime[ie];
		obs[8][ie]=ATReCP[ie];
		obs[9][ie]=ATImCP[ie];
		obs[10][ie]=A5[ie];
		obs[11][ie]=AFBCP[ie];
		obs[12][ie]=A8[ie];
		obs[13][ie]=A9[ie];
		obs[14][ie]=P2CP[ie];
		obs[15][ie]=P3CP[ie];
		obs[16][ie]=P5primeCP[ie];
		obs[17][ie]=P8primeCP[ie];
		obs[18][ie]=A6c[ie];
		obs[19][ie]=Q8m[ie];
		obs[20][ie]=Q8p[ie];
		obs[21][ie]=Q9[ie];
		obs[22][ie]=S2s[ie];
		/*********************************************/
	}

	return 1./(1.-ys*ys)*(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2)/2.;
}

/*----------------------------------------------------------------------*/

double dGamma_Bsphill_dq2_soft_calculator(int gen, double q2, double obs[][3], char name[])
/* "container" function scanning the SLHA file "name" and calculating dGamma/dq2(Bs->phi mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(gen,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(gen,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(gen,Cpb,CQpb,mu_W,mu_b,&param);

	return dGamma_Bsphill_dq2_soft(gen,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBsphill_soft(int gen, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	int ie,je;
	int nmax=10;
	if((smin<0.099)||(smax-smin>10.)) nmax=100;
	double Gamma=0.;
	double s;
	
	double obs_num[Nobs_Bsphill+1],obs_den[Nobs_Bsphill+1];
	for(je=0;je<=Nobs_Bsphill;je++) obs_num[je]=obs_den[je]=0.;

	obs[0]=0.; /* none */
	obs[1]=0.; /* integrated FL */
	obs[2]=0.; /* integrated AT2 */
	obs[3]=0.; /* integrated S3 */
	obs[4]=0.; /* integrated S4 */
	obs[5]=0.; /* integrated P4prime */
	obs[6]=0.; /* integrated S7 */
	obs[7]=0.; /* integrated P6prime */
	obs[8]=0.; /* integrated ATReCP */
	obs[9]=0.; /* integrated ATReIm */
	obs[10]=0.; /* integrated A5 */
	obs[11]=0.; /* integrated AFBCP */
	obs[12]=0.; /* integrated A8 */
	obs[13]=0.; /* integrated A9 */
	obs[14]=0.; /* integrated P2CP */
	obs[15]=0.; /* integrated P3CP */
	obs[16]=0.; /* integrated P5primeCP */
	obs[17]=0.; /* integrated P8primeCP */
	obs[18]=0.; /* integrated A6c */
	obs[19]=0.; /* integrated Q8m */
	obs[20]=0.; /* integrated Q8p */
	obs[21]=0.; /* integrated Q9 */
	obs[22]=0.; /* integrated S2s */
	/*********************************************/

	double dobs[Nobs_Bsphill+1][3];
	double s0m,s0p,s0;
		
	s0=s0m=0.;
	s0p=1.;
	
	double h=(smax-smin)/nmax;	
	s=smin;
	Gamma=dGamma_Bsphill_dq2_soft(gen,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_Bsphill;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}
	
	for(ie=1;ie<nmax;ie++)	
	{
		s+=h;

		Gamma+=4.*dGamma_Bsphill_dq2_soft(gen,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_Bsphill;je++) 
		{
			obs_num[je]+=4.*dobs[je][1];
			obs_den[je]+=4.*dobs[je][2];
		}

		Gamma+=2.*dGamma_Bsphill_dq2_soft(gen,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_Bsphill;je++) 
		{
			obs_num[je]+=2.*dobs[je][1];
			obs_den[je]+=2.*dobs[je][2];
		}
	}
	
	s=smax;
	Gamma+=4.*dGamma_Bsphill_dq2_soft(gen,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_Bsphill;je++) 
	{
		obs_num[je]+=4.*dobs[je][1];
		obs_den[je]+=4.*dobs[je][2];
	}	
	Gamma+=dGamma_Bsphill_dq2_soft(gen,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_Bsphill;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}	
	
	Gamma*=h/6.;
	for(je=1;je<=Nobs_Bsphill;je++) 
	{
		obs_num[je]*=h/6.;
		obs_den[je]*=h/6.;
	}
	
	obs[0]=s0;
	for(je=1;je<=Nobs_Bsphill;je++)
	if(je==5 || je==7 || je==16 || je==17)
	{
		obs[je]=obs_num[je]/sqrt(fabs(obs_den[je]*obs_den[2]/2.));
	}
	else obs[je]=obs_num[je]/obs_den[je];
	
	for(je=1;je<=Nobs_Bsphill;je++) if(fabs(obs[je])<1.e-15) obs[je]=0.;

	return param->life_Bs/hbar*Gamma;
}

/*----------------------------------------------------------------------*/

double BRBsphill_lowq2_soft(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBsphill_soft(gen,1.,6.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBsphill_highq2_soft(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBsphill_soft(gen,14.18,16.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_Bsphill_lowq2_soft_calculator(char name[], double obs[])
/* "container" function scanning the SLHA file "name" and calculating BR(Bs->phi mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBsphill_lowq2_soft(2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_Bsphill_highq2_soft_calculator(char name[], double obs[])
/* "container" function scanning the SLHA file "name" and calculating BR(Bs->phi mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBsphill_highq2_soft(2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/
/*------------------------------ FULL ----------------------------------*/
/*----------------------------------------------------------------------*/

double dGamma_Bsphill_dq2_full(int gen, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{	
	double ml;
	if(gen==1) ml=param->mass_e;
	else if(gen==3) ml=param->mass_tau;
	else ml=param->mass_mu;
	
	double mBs=param->m_Bs;
	double mphi=param->m_phi;
	double eq=-1./3.;

	double mc=mc_pole_1loop(param);
	double mbpole=mb_pole_1loop(param);
	double mb_mub=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);

	int ie,je;
	
	double beta_l=sqrt(1.-4.*ml*ml/q2);

	double alpha_em=1./133.;

	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	
	double mu_f=sqrt(mu_b*0.5);
	
	double alphas_muf=alphas_running(mu_f,param->mass_top_pole,param->mass_b_pole,param);
	double eta=alphas_muf/alphas_running(1.,param->mass_top_pole,param->mass_b_pole,param);

	double Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_phi=(mBs*mBs+mphi*mphi-q2)/2./mBs;

	int nf=5;
	double f_phi_perp=param->f_phi_perp;
	f_phi_perp*=pow(eta,4./3./(11.-2./3.*nf));

	double f_phi_par=param->f_phi_par;
	
	double complex ALperp=0.;
	double complex ARperp=0.;
	double complex ALpar=0.;
	double complex ARpar=0.;
	double complex AL0=0.;
	double complex AR0=0.;
	double complex At=0.;
	double complex AS=0.;

	double complex ALperp_bar=0.;
	double complex ARperp_bar=0.;
	double complex ALpar_bar=0.;
	double complex ARpar_bar=0.;
	double complex AL0_bar=0.;
	double complex AR0_bar=0.;
	double complex At_bar=0.;
	double complex AS_bar=0.;

	double complex ALperp_high=0.;
	double complex ARperp_high=0.;
	double complex ALpar_high=0.;
	double complex ARpar_high=0.;
	double complex AL0_high=0.;
	double complex AR0_high=0.;
	double complex At_high=0.;
	double complex AS_high=0.;

	double complex ALperp_bar_high=0.;
	double complex ARperp_bar_high=0.;
	double complex ALpar_bar_high=0.;
	double complex ARpar_bar_high=0.;
	double complex AL0_bar_high=0.;
	double complex AR0_bar_high=0.;
	double complex At_bar_high=0.;
	double complex AS_bar_high=0.;
		

	double V = 0.0;
	double A1 = 0.0;
	double A2 = 0.0;
	double A0 = 0.0;
	double T1 = 0.0;
	double T2 = 0.0;
	double T23 = 0.0;
	double T3 = 0.0;
	double A12 = 0.0;

	double MB, MV;
	double a0V, a1V, a2V, a0A0, a1A0, a2A0, a0A1, a1A1, a2A1, a0A12, a1A12, a2A12, a0T1, a1T1, a2T1, a0T2, a1T2, a2T2, a0T23, a1T23, a2T23;
	MB=MV=a0V=a1V=a2V=a0A0=a1A0=a2A0=a0A1=a1A1=a2A1=a0A12=a1A12=a2A12=a0T1=a1T1=a2T1=a0T2=a1T2=a2T2=a0T23=a1T23=a2T23 = 0.;

	/********LCSR+Lattice fit from Barucha,Straub,Zwicky 1503.05534***************/
	if(param->Bsphi_FormFactor_choice ==1)
	{
		double P_V=1./(1.-q2/pow(param->MV_Bsphi,2.));
		double P_A1=1./(1.-q2/pow(param->MA1_Bsphi,2.));
		double P_A12=1./(1.-q2/pow(param->MA12_Bsphi,2.));

		double tau_plus=pow(mBs+mphi,2.);
		double tau_minus=pow(mBs-mphi,2.);
		double tau_0=tau_plus*(1.-sqrt(1.-tau_minus/tau_plus));
		double z_q2=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));
		double z_0=(sqrt(tau_plus)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus)+sqrt(tau_plus-tau_0));

		V=P_V*(param->a0V_Bsphi+param->a1V_Bsphi*(z_q2-z_0)+param->a2V_Bsphi*(pow(z_q2-z_0,2.)));
		A1=P_A1*(param->a0A1_Bsphi+param->a1A1_Bsphi*(z_q2-z_0)+param->a2A1_Bsphi*(pow(z_q2-z_0,2.)));
		A12=P_A12*(param->a0A12_Bsphi+param->a1A12_Bsphi*(z_q2-z_0)+param->a2A12_Bsphi*(pow(z_q2-z_0,2.)));
		A2=(pow(mBs+mphi,2.)*(mBs*mBs-mphi*mphi-q2)*A1-16.*mBs*mphi*mphi*(mBs+mphi)*A12)/((pow(mBs+mphi,2.)-q2)*(pow(mBs-mphi,2.)-q2));

		double P_A0=1./(1.-q2/pow(param->MA0_Bsphi,2.));
		double P_T1=1./(1.-q2/pow(param->MT1_Bsphi,2.));
		double P_T2=1./(1.-q2/pow(param->MT2_Bsphi,2.));
		double P_T23=1./(1.-q2/pow(param->MT23_Bsphi,2.));

		A0=P_A0*(param->a0A0_Bsphi+param->a1A0_Bsphi*(z_q2-z_0)+param->a2A0_Bsphi*(pow(z_q2-z_0,2.)));
		T1=P_T1*(param->a0T1_Bsphi+param->a1T1_Bsphi*(z_q2-z_0)+param->a2T1_Bsphi*(pow(z_q2-z_0,2.)));
		T2=P_T2*(param->a0T2_Bsphi+param->a1T2_Bsphi*(z_q2-z_0)+param->a2T2_Bsphi*(pow(z_q2-z_0,2.)));
		T23=P_T23*(param->a0T23_Bsphi+param->a1T23_Bsphi*(z_q2-z_0)+param->a2T23_Bsphi*(pow(z_q2-z_0,2.)));


	}

	/********Fit parameters from Gubernari, Reboud, van Dyk, Virto, 2305.06301***************/
	if(param->Bsphi_FormFactor_choice ==2)
	{
		MB = mBs;
		MV = mphi;
		double MVT1_Bsphi = 5.415;
		double MA0_Bsphi = 5.367;
		double MAT_Bsphi = 5.829;
		a0V = param->a0V_Bsphi_GRvDV_BSZ;
		a1V = param->a1V_Bsphi_GRvDV_BSZ;
		a2V = param->a2V_Bsphi_GRvDV_BSZ;
		a0A0 = param->a0A0_Bsphi_GRvDV_BSZ;
		a1A0 = param->a1A0_Bsphi_GRvDV_BSZ;
		a2A0 = param->a2A0_Bsphi_GRvDV_BSZ;
		a0A1 = param->a0A1_Bsphi_GRvDV_BSZ;
		a1A1 = param->a1A1_Bsphi_GRvDV_BSZ;
		a2A1 = param->a2A1_Bsphi_GRvDV_BSZ;
		a0A12 = (MB*MB - MV*MV)/(8.*MB*MV)*param->a0A0_Bsphi_GRvDV_BSZ;// Eq.17 in BSZ [1503.05534]
		a1A12 = param->a1A12_Bsphi_GRvDV_BSZ;
		a2A12 = param->a2A12_Bsphi_GRvDV_BSZ;
		a0T1 = param->a0T1_Bsphi_GRvDV_BSZ;
		a1T1 = param->a1T1_Bsphi_GRvDV_BSZ;
		a2T1 = param->a2T1_Bsphi_GRvDV_BSZ;
		a0T2 = param->a0T1_Bsphi_GRvDV_BSZ;// Eq.17 in BSZ [1503.05534]
		a1T2 = param->a1T2_Bsphi_GRvDV_BSZ;
		a2T2 = param->a2T2_Bsphi_GRvDV_BSZ;
		a0T23 = param->a0T23_Bsphi_GRvDV_BSZ;
		a1T23 = param->a1T23_Bsphi_GRvDV_BSZ;
		a2T23 = param->a2T23_Bsphi_GRvDV_BSZ;

		double P_V=1./(1.-q2/pow(MVT1_Bsphi,2.));
		double P_A0=1./(1.-q2/pow(MA0_Bsphi,2.));
		double P_A1=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_A12=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_T1=1./(1.-q2/pow(MVT1_Bsphi,2.));
		double P_T2=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_T23=1./(1.-q2/pow(MAT_Bsphi,2.));

		double tau_plus=pow(MB+MV,2.);
		double tau_minus=pow(MB-MV,2.);
		double tau_0=tau_plus-sqrt((tau_plus-tau_minus)*tau_plus);
		double z_q2tau0=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));
		double z_0tau0=(sqrt(tau_plus)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus)+sqrt(tau_plus-tau_0));
		double z_diff = z_q2tau0 - z_0tau0;


		V = P_V * ( a0V + a1V*z_diff + a2V* (z_diff*z_diff) );
		A0 = P_A0 * ( a0A0 + a1A0*z_diff + a2A0*(z_diff*z_diff) );
		A1 = P_A1 * ( a0A1 + a1A1*z_diff + a2A1*(z_diff*z_diff) );
		A12 = P_A12 * ( a0A12 + a1A12*z_diff + a2A12*(z_diff*z_diff) );
		T1 = P_T1 * ( a0T1 + a1T1*z_diff + a2T1*(z_diff*z_diff) );
		T2 = P_T2 * ( a0T2 + a1T2*z_diff + a2T2*(z_diff*z_diff) );
		T23 = P_T23 * ( a0T23 + a1T23*z_diff + a2T23*(z_diff*z_diff) );

		A2=(pow(MB+MV,2.)*(MB*MB-MV*MV-q2)*A1-16.*MB*MV*MV*(MB+MV)*A12)/((pow(MB+MV,2.)-q2)*(pow(MB-MV,2.)-q2));

	}

	/********LCSR-only from Barucha,Straub,Zwicky 1503.05534***************/
	if(param->Bsphi_FormFactor_choice ==5)
	{
		MB = mBs;
		MV = mphi;
		a0V = param->a0V_Bsphi_BSZ_LCSRonly;
		a1V = param->a1V_Bsphi_BSZ_LCSRonly;
		a2V = param->a2V_Bsphi_BSZ_LCSRonly;
		a0A0 = param->a0A0_Bsphi_BSZ_LCSRonly;
		a1A0 = param->a1A0_Bsphi_BSZ_LCSRonly;
		a2A0 = param->a2A0_Bsphi_BSZ_LCSRonly;
		a0A1 = param->a0A1_Bsphi_BSZ_LCSRonly;
		a1A1 = param->a1A1_Bsphi_BSZ_LCSRonly;
		a2A1 = param->a2A1_Bsphi_BSZ_LCSRonly;
		a0A12 = param->a0A12_Bsphi_BSZ_LCSRonly;
		a1A12 = param->a1A12_Bsphi_BSZ_LCSRonly;
		a2A12 = param->a2A12_Bsphi_BSZ_LCSRonly;
		a0T1 = param->a0T1_Bsphi_BSZ_LCSRonly;
		a1T1 = param->a1T1_Bsphi_BSZ_LCSRonly;
		a2T1 = param->a2T1_Bsphi_BSZ_LCSRonly;
		a0T2 = param->a0T2_Bsphi_BSZ_LCSRonly;
		a1T2 = param->a1T2_Bsphi_BSZ_LCSRonly;
		a2T2 = param->a2T2_Bsphi_BSZ_LCSRonly;
		a0T23 = param->a0T23_Bsphi_BSZ_LCSRonly;
		a1T23 = param->a1T23_Bsphi_BSZ_LCSRonly;
		a2T23 = param->a2T23_Bsphi_BSZ_LCSRonly;

		double P_V=1./(1.-q2/pow(param->MV_Bsphi,2.));
		double P_A0=1./(1.-q2/pow(param->MA0_Bsphi,2.));
		double P_A1=1./(1.-q2/pow(param->MA1_Bsphi,2.));
		double P_A12=1./(1.-q2/pow(param->MA12_Bsphi,2.));
		double P_T1=1./(1.-q2/pow(param->MT1_Bsphi,2.));
		double P_T2=1./(1.-q2/pow(param->MT2_Bsphi,2.));
		double P_T23=1./(1.-q2/pow(param->MT23_Bsphi,2.));

		double tau_plus=pow(MB+MV,2.);
		double tau_minus=pow(MB-MV,2.);
		double tau_0=tau_plus-sqrt((tau_plus-tau_minus)*tau_plus);
		double z_q2tau0=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));
		double z_0tau0=(sqrt(tau_plus)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus)+sqrt(tau_plus-tau_0));
		double z_diff = z_q2tau0 - z_0tau0;


		V = P_V * ( a0V + a1V*z_diff + a2V* (z_diff*z_diff) );
		A0 = P_A0 * ( a0A0 + a1A0*z_diff + a2A0* (z_diff*z_diff) );
		A1 = P_A1 * ( a0A1 + a1A1*z_diff + a2A1* (z_diff*z_diff) );
		A12 = P_A12 * ( a0A12 + a1A12*z_diff + a2A12* (z_diff*z_diff) );
		T1 = P_T1 * ( a0T1 + a1T1*z_diff + a2T1* (z_diff*z_diff) );
		T2 = P_T2 * ( a0T2 + a1T2*z_diff + a2T2* (z_diff*z_diff) );
		T23 = P_T23 * ( a0T23 + a1T23*z_diff + a2T23* (z_diff*z_diff) );

		A2=(pow(MB+MV,2.)*(MB*MB-MV*MV-q2)*A1-16.*MB*MV*MV*(MB+MV)*A12)/((pow(MB+MV,2.)-q2)*(pow(MB-MV,2.)-q2));


	}

	/********Fit parameters from Horgan, Liu, Meinel, Wingate 1310.3722 + 1501.00367 ***************/
	if(param->Bsphi_FormFactor_choice == 6)
	{
		MB = mBs;
		MV = mphi;
		double MVT1_Bsphi = 5.411;//Mres=mBs+deltam = 5.366 + 0.045
		double MA0_Bsphi = 5.366;//Mres=mBs+deltam = 5.366 + 0.000
		double MAT_Bsphi = 5.806;//Mres=mBs+deltam = 5.366 + 0.440
		a0V   = param->a0V_Bsphi_HLMW;
		a1V   = param->a1V_Bsphi_HLMW;
		a0A0  = param->a0A0_Bsphi_HLMW;
		a1A0  = param->a1A0_Bsphi_HLMW;
		a0A1  = param->a0A1_Bsphi_HLMW;
		a1A1  = param->a1A1_Bsphi_HLMW;
		a0A12 = param->a0A12_Bsphi_HLMW;
		a1A12 = param->a1A12_Bsphi_HLMW;
		a0T1  = param->a0T1_Bsphi_HLMW;
		a1T1  = param->a1T1_Bsphi_HLMW;
		a0T2  = param->a0T2_Bsphi_HLMW;
		a1T2  = param->a1T2_Bsphi_HLMW;
		a0T23 = param->a0T23_Bsphi_HLMW;
		a1T23 = param->a1T23_Bsphi_HLMW;
		double sysErr = param->systErr_Bsphi_HLMW;

		double P_V=1./(1.-q2/pow(MVT1_Bsphi,2.));
		double P_A0=1./(1.-q2/pow(MA0_Bsphi,2.));
		double P_A1=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_A12=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_T1=1./(1.-q2/pow(MVT1_Bsphi,2.));
		double P_T2=1./(1.-q2/pow(MAT_Bsphi,2.));
		double P_T23=1./(1.-q2/pow(MAT_Bsphi,2.));

		double tau_plus=pow(MB+MV,2.);
		double tau_minus=pow(MB-MV,2.);
		double tau_0=12.;// Page 9 in [1310.3722]
		double z_q2tau0=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));

		/* The (1. + sysErr) is there to take into account 5% systematic error (e.g. see p.18 in 1310.3722)*/
		V = P_V * ( a0V + a1V*z_q2tau0  ) * (1. + sysErr);
		A0 = P_A0 * ( a0A0 + a1A0*z_q2tau0  ) * (1. + sysErr);
		A1 = P_A1 * ( a0A1 + a1A1*z_q2tau0  ) * (1. + sysErr);
		A12 = P_A12 * ( a0A12 + a1A12*z_q2tau0  ) * (1. + sysErr);
		T1 = P_T1 * ( a0T1 + a1T1*z_q2tau0  ) * (1. + sysErr);
		T2 = P_T2 * ( a0T2 + a1T2*z_q2tau0  ) * (1. + sysErr);
		T23 = P_T23 * ( a0T23 + a1T23*z_q2tau0 ) * (1. + sysErr);

		A2=(pow(MB+MV,2.)*(MB*MB-MV*MV-q2)*A1-16.*MB*MV*MV*(MB+MV)*A12)/((pow(MB+MV,2.)-q2)*(pow(MB-MV,2.)-q2));

	}



/****************************************************************************/

	double xi_par=(mBs+mphi)/2./E_phi*A1-(mBs-mphi)/mBs*A2;
	double xi_perp=mBs/(mBs+mphi)*V;
		

	if(q2<HIGH_q2)
	{
	 	double complex C7eff=Cmub[7];
		double complex C8eff=Cmub[8];
		double complex C9=Cmub[9];
		double complex C10=Cmub[10];
	
		double complex C7effp=Cpb[7];
		double complex C9p=Cpb[9];
		double complex C10p=Cpb[10];
	
		double complex C1bar=Cmub[1]/2.;
		double complex C2bar=Cmub[2]-Cmub[1]/6.;
		double complex C3bar=Cmub[3]-Cmub[4]/6.+16.*Cmub[5]-8./3.*Cmub[6];
		double complex C4bar=Cmub[4]/2.+8.*Cmub[6];
		double complex C5bar=Cmub[3]-Cmub[4]/6.+4.*Cmub[5]-2./3.*Cmub[6];
		double complex C6bar=Cmub[4]/2.+2.*Cmub[6];
	
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];

		double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
		double mb=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */
				
					
		double complex h_mc=h_bkll(q2,mc,mu_b);
		double complex h_mb=h_bkll(q2,mbpole,mu_b);
		double complex h_0=h_bkll(q2,0.,mu_b);

		double complex Y=4./3.*Cmub[3]+64./9.*Cmub[5]+64./27.*Cmub[6]
		+h_mc*(4./3.*Cmub[1]+Cmub[2]+6.*Cmub[3]+60.*Cmub[5])
		+h_mb*(-7./2.*Cmub[3]-2./3.*Cmub[4]-38.*Cmub[5]-32./3.*Cmub[6])
		+h_0*(-1./2.*Cmub[3]-2./3.*Cmub[4]-8.*Cmub[5]-32./3.*Cmub[6]);

		double Yu=(h_mc-h_0)*(4./3.*Cmub[1]+Cmub[2]);

		double complex Cperpp0=0.;
		double complex Cperpm0=0.;
		//double complex Cparp0=0.;
		double complex Cparm0=0.;

// 		double complex Cperpp0u=q2/2./mb/mBs*Yu;
// 		double complex Cperpm0u=q2/2./mb/mBs*Yu;
// 		//double complex Cparp0u=-mBs/2./mb*Yu;
// 		double complex Cparm0u=-mBs/2./mb*Yu;
		/*******Used in redefinition of Tauperpu**********/
		double complex Cperp0u=q2/2./mb/mBs*Yu;
		double complex Cpar0u=-mBs/2./mb*Yu;

		double logb=log(mu_b/mb);
		double DeltaM=-6.*logb-4.*(1.-mu_f/mb);
		//double L=-(mb*mb-q2)/q2*log(1.-q2/mb/mb);
	
		double shat=q2/mb/mb; 

		double mchat=mc/mb; 
		double z=mchat*mchat;	

		double complex Cperppf=0.;
		double complex Cperpmf=0.;
		//double complex Cparpf=0.;
		double complex Cparmf=0.;

		double complex F27=F27_bsll(shat,z,logb);
		double complex F87=F87_bsll(shat,logb);
		double complex F29=F29_bkll(shat,z,logb);
		double complex F19=F19_bkll(shat,z,logb);
		double complex F89=F89_bsll(shat);
		double complex F27_u=F27u(shat,logb);
		double complex F29_u=F29u(shat,logb);
		double complex F19_u=F19u(shat,logb);

		double complex Cperpnf=(-C2bar*F27-C8eff*F87
		-q2/2./mb/mBs*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;
		
		double complex Cparnf=(C2bar*F27+C8eff*F87
		+mBs/2./mb*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;

		double complex Cperpnfu=(-C2bar*(F27+F27_u)
		-q2/2./mb/mBs*(C2bar*(F29+F29_u)+2.*C1bar*(F19+F19_u+1./6.*(F29+F29_u))))/4.*3.;
		
		double complex Cparnfu=(C2bar*(F27+F27_u)
		+mBs/2./mb*(C2bar*(F29+F29_u)+2.*C1bar*(F19+F19_u+1./6.*(F29+F29_u))))/4.*3.;
		
		double complex Cperpp1=Cperppf+Cperpnf;
		double complex Cperpm1=Cperpmf+Cperpnf;
		//double complex Cparp1=Cparpf+Cparnf;
		double complex Cparm1=Cparmf+Cparnf;
	
		double complex Cperpp=Cperpp0+alphas_mub*4./3./4./pi*Cperpp1;
		double complex Cperpm=Cperpm0+alphas_mub*4./3./4./pi*Cperpm1;
		//double complex Cparp=Cparp0+alphas_mub*4./3./4./pi*Cparp1;
		double complex Cparm=Cparm0+alphas_mub*4./3./4./pi*Cparm1;
				
// 		double complex Cperppu=Cperpp0u+alphas_mub*4./3./4./pi*Cperpnfu; 
// 		double complex Cperpmu=Cperpm0u+alphas_mub*4./3./4./pi*Cperpnfu; 
// 		//double complex Cparpu=Cparp0u+alphas_mub*4./3./4./pi*Cparnfu; 
// 		double complex Cparmu=Cparm0u+alphas_mub*4./3./4./pi*Cparnfu; 
		/*******Used in redefinition of Tauperpu**********/
		double complex Cperpu=Cperp0u+alphas_mub*4./3./4./pi*Cperpnfu; 
		double complex Cparu=Cpar0u+alphas_mub*4./3./4./pi*Cparnfu; 
		
		
		double Xi_perp=1.;
		double Xi_par=mphi/E_phi;
	
		double eu=2./3.;
		double ed=-1./3.;
	
		double a1phi_perp=param->a1phi_perp;
		double a2phi_perp=param->a2phi_perp;
		double a1phi_par=param->a1phi_par;
		double a2phi_par=param->a2phi_par;

		a1phi_perp*=pow(eta,4./3.*(4.*1./2.)/(11.-2./3.*nf));
		a2phi_perp*=pow(eta,4./3.*(4.*(1./2.+1./3.))/(11.-2./3.*nf));

		a1phi_par*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
		a2phi_par*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

		double u;
	
		double complex int_perppp,int_perppm,int_perpmp,int_perpmm;
		double complex int_parpp,int_parpm,int_parmp,int_parmm;
		double complex int_perppu,int_parpu,int_parmu;
		double complex Tperppp0,Tperpppf,Tperpppnf,Tperppp;
		double complex Tperppm0,Tperppmf,Tperppmnf,Tperppm;
		double complex Tperpmp0,Tperpmpf,Tperpmpnf,Tperpmp;
		double complex Tperpmm0,Tperpmmf,Tperpmmnf,Tperpmm;
		double complex Tparpp0,Tparppf,Tparppnf,Tparpp;
		double complex Tparpm0,Tparpmf,Tparpmnf,Tparpm;
		double complex Tparmp0,Tparmpf,Tparmpnf,Tparmp;
		double complex Tparmm0,Tparmmf,Tparmmnf,Tparmm;
		double complex Tperppnfu,Tperppu;
		double complex Tparp0u,Tparpnfu,Tparpu;
		double complex Tparm0u,Tparmnfu,Tparmu;
		
		int_perppp=int_perppm=int_perpmp=int_perpmm=0.;
		int_perppu=0.;
		int_parpp=int_parpm=int_parmp=int_parmm=0.;
		int_parpu=int_parmu=0.;


		double lambda_Bp=param->lambda_Bsp; 
		lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_f,2.))*(1.-2.*1.4);

		double omega0=2.*(mBs-mb)/3.; 
		double complex lambda_Bm=1./(exp(-q2/mBs/omega0)/omega0*(-Ei(q2/mBs/omega0)+I*pi));

		double phiphi_perp,phiphi_par;
		double complex tperp_mb,tperp_mc,tperp_0;
		double complex tpar_mb,tpar_mc,tpar_0;

		double complex integ3=0;

		double complex Fperp=0.;
		double complex Xperp=0.;
		double x;
		double complex integ4=0.;
		double complex FV;
		double complex integ4u=0.;
		double complex FVu;
		double complex integ5=0.;
		double complex integ5u=0.;

		double zeta3A=param->zeta3A;
		double zeta3V=param->zeta3V;
		double wA10=param->wA10;
		double deltatp=param->deltatp_phi;
		double deltatm=param->deltatm_phi;
		
		int n1=10;
		int n1sav=n1;
		for(ie=0;ie<=n1;ie++)
		{
			u=(double)ie/n1;
			if(ie==0) n1*=2;
			if(ie==n1){u=0.99;n1*=2;}

		/* Tperp */		
			Tperppp0=Tperpmp0=0.;
			Tperpppf=0.;
			Tperpmpf=0.;

			Tperppm0=Tperpmm0=0.;
			Tperppmf=Tperpmmf=0.;
			
			phiphi_perp=phi_Kstar(u,a1phi_perp,a2phi_perp);
			tperp_mc=tperp_bkll(u,mc,q2,E_phi,param);
			tperp_mb=tperp_bkll(u,mb,q2,E_phi,param);
			tperp_0=tperp_bkll(u,0.,q2,E_phi,param);

			Tperpppnf=Tperpmpnf=-4.*ed*C8eff/(u+(1.-u)*q2/mBs/mBs)
			+mBs/2./mb*(eu*tperp_mc*(C2bar+C4bar-C6bar)
			+ed*tperp_mb*(C3bar+C4bar-C6bar-4.*mb/mBs*C5bar)
			+ed*tperp_0*C3bar);
				
			Tperppmnf=Tperpmmnf=0.;

			Tperppnfu=mBs/2./mb*eu*(tperp_mc-tperp_0)*(Cmub[2]-Cmub[1]/6.);
		
			Tperppp=Tperppp0+alphas_muf*4./3./4./pi*(Tperpppf+Tperpppnf); 
			Tperppm=Tperppm0+alphas_muf*4./3./4./pi*(Tperppmf+Tperppmnf);
			Tperpmp=Tperpmp0+alphas_muf*4./3./4./pi*(Tperpmpf+Tperpmpnf);
			Tperpmm=Tperpmm0+alphas_muf*4./3./4./pi*(Tperpmmf+Tperpmmnf);
	
			Tperppu=alphas_muf*4./3./4./pi*Tperppnfu;
	
			int_perppp+=phiphi_perp*Tperppp/n1/lambda_Bp; 
			int_perppm+=phiphi_perp*Tperppm/n1/lambda_Bm;
			int_perpmp+=phiphi_perp*Tperpmp/n1/lambda_Bp;
			int_perpmm+=phiphi_perp*Tperpmm/n1/lambda_Bm;

			int_perppu+=phiphi_perp*Tperppu/n1/lambda_Bp;


		/* Tpar */		

			phiphi_par=phi_Kstar(u,a1phi_par,a2phi_par);
			tpar_mc=tpar_bkll(u,mc,q2,E_phi,param);
			tpar_mb=tpar_bkll(u,mb,q2,E_phi,param);
			tpar_0=tpar_bkll(u,0.,q2,E_phi,param);
		
			Tparpp0=Tparmp0=0.;
		
			Tparppf=0.;
			Tparmpf=0.;

			Tparppnf=Tparmpnf=mBs/mb*(eu*tpar_mc*(C2bar+C4bar-C6bar)
			+ed*tpar_mb*(C3bar+C4bar-C6bar)
			+ed*tpar_0*C3bar);
				
			Tparpnfu=mBs/mb*eu*(tpar_mc-tpar_0)*(Cmub[2]-Cmub[1]/6.);
	
			Tparpu=alphas_muf*4./3./4./pi*Tparpnfu;
	
			Tparpm0=Tparmm0=-eq*4.*mBs/mb*((C3bar+3.*C4bar)+12.*(Cmub[3]+10.*Cmub[5]));
			
			Tparp0u=0.;

			Tparm0u=-eq*4.*mBs/mb*(-(4./3.*Cmub[1]+Cmub[2]));
		
			Tparpmf=Tparmmf=0.;

			h_mc=h_bkll((1.-u)*mBs*mBs+u*q2,mc,mu_b);
			h_mb=h_bkll((1.-u)*mBs*mBs+u*q2,mbpole,mu_b);
			h_0=h_bkll((1.-u)*mBs*mBs+u*q2,0.,mu_b);

			Tparpmnf=Tparmmnf=eq*(8.*C8eff/((1.-u)+u*q2/mBs/mBs)
			+6.*mBs/mb*(h_mc*(C2bar+C4bar+C6bar)
			+h_mb*(C3bar+C4bar+C6bar)
			+h_0*(C3bar+3.*C4bar+3.*C6bar)
			-8./27.*(C3bar-C5bar-15.*C6bar)));
	
			Tparpp=Tparpp0+alphas_muf*4./3./4./pi*(Tparppf+Tparppnf);
			Tparpm=Tparpm0+alphas_muf*4./3./4./pi*(Tparpmf+Tparpmnf);
			Tparmp=Tparmp0+alphas_muf*4./3./4./pi*(Tparmpf+Tparmpnf);
			Tparmm=Tparmm0+alphas_muf*4./3./4./pi*(Tparmmf+Tparmmnf);
			
			Tparmnfu=eq*(6.*mBs/mb*(h_mc-h_0)*(Cmub[2]-Cmub[1]/6.));

			Tparpu=Tparp0u+alphas_muf*4./3./4./pi*Tparpnfu;
	
			Tparmu=Tparm0u+alphas_muf*4./3./4./pi*Tparmnfu;


			int_parpp+=(phiphi_par*Tparpp/lambda_Bp)/n1;
			int_parpm+=(phiphi_par*Tparpm/lambda_Bm)/n1;
			
			int_parmp+=(phiphi_par*Tparmp/lambda_Bp)/n1;
			int_parmm+=(phiphi_par*Tparmm/lambda_Bm)/n1;
			
			int_parpu+=(phiphi_par*Tparpu/lambda_Bp)/n1;
			
			int_parmu+=(phiphi_par*Tparmu/lambda_Bm)/n1;


			integ3+=phiphi_perp/((1.-u)+u*q2/mBs/mBs)/n1;

			x=(1.-u)*mBs*mBs+u*q2;
			
			h_mc=h_bkll(x,mc,mu_b);
			h_mb=h_bkll(x,mbpole,mu_b);
			h_0=h_bkll(x,0.,mu_b);
			
			FV=3./4.*(h_mc*(C2bar+C4bar+C6bar)+h_mb*(C3bar+C4bar+C6bar)+h_0*(C3bar+3.*C4bar+3.*C6bar)-8./27.*(C3bar-C5bar-15.*C6bar));

			FVu=3./4.*(h_mc-h_0)*(Cmub[2]-Cmub[1]/6.);
		

			integ4+=phiphi_perp/((1.-u)+u*q2/mBs/mBs)*FV/n1;
			integ4u+=phiphi_perp/((1.-u)+u*q2/mBs/mBs)*FVu/n1;
		
			Fperp+=phiphi_perp/((1.-u)+u*q2/mBs/mBs)/3./n1;
			Xperp+=(u<=1.-0.5/mBs)*phiphi_perp/pow((1.-u)+u*q2/mBs/mBs,2.)/3./n1;
		
			integ5+=((3./4.*(1.+pow(2.*u-1.,2.))+a1phi_par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2phi_par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2phi_par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1phi_par*(2.*u-1.)+(a2phi_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1phi_par*u+(a2phi_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FV/n1;

			integ5u+=((3./4.*(1.+pow(2.*u-1.,2.))+a1phi_par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2phi_par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2phi_par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1phi_par*(2.*u-1.)+(a2phi_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1phi_par*u+(a2phi_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FVu/n1;

			if(ie==0||ie==n1sav) n1=n1sav;
		}
		
		/* Tau_perp */		

		double complex Tauperpp=xi_perp*Cperpp+pi*pi/3.*param->f_Bs*f_phi_perp/mBs*Xi_perp*(int_perppp+int_perppm); 
		
		double complex Tauperpm=xi_perp*Cperpm+pi*pi/3.*param->f_Bs*f_phi_perp/mBs*Xi_perp*(int_perpmp+int_perpmm); 
		
		double complex Tauperpu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_perp*Cperpu+pi*pi/3.*param->f_Bs*f_phi_perp/mBs*Xi_perp*int_perppu);
		double complex Tauperpu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_perp*Cperpu+pi*pi/3.*param->f_Bs*f_phi_perp/mBs*Xi_perp*int_perppu);

		
		/* Tau_par */		

		//double complex Tauparp=xi_par*Cparp+pi*pi/3.*param->f_Bs*f_phi_par/mBs*Xi_par*(int_parpp+int_parpm);
		
		double complex Tauparm=xi_par*Cparm+pi*pi/3.*param->f_Bs*f_phi_par/mBs*Xi_par*(int_parmp+int_parmm);
		
		double complex Tauparu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_par*Cparu+pi*pi/3.*param->f_Bs*f_phi_par/mBs*Xi_par*(int_parpu+int_parmu));
		double complex Tauparu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_par*Cparu+pi*pi/3.*param->f_Bs*f_phi_par/mBs*Xi_par*(int_parpu+int_parmu));
				
				
		double complex DeltaTauperpWA=-eq*4.*pi*pi/3.*param->f_Bs*f_phi_perp/mb/mBs*(Cmub[3]+4./3.*Cmub[4]+4.*Cmub[5]+16./3.*Cmub[6])*integ3
		+eq*2.*pi*pi/3.*param->f_Bs*f_phi_par/mb/mBs*mphi/(1.-q2/mBs/mBs)/lambda_Bp*(Cmub[3]+4./3.*Cmub[4]+16.*Cmub[5]+64./3.*Cmub[6]+12.*(Cmub[3]+10.*Cmub[5]));

		double complex DeltaTauperpuWA=-eq*2.*pi*pi/3.*param->f_Bs*f_phi_par/mb/mBs*mphi/(1.-q2/mBs/mBs)/lambda_Bp*(4./3.*Cmub[1]+Cmub[2]);


		double rho=0.;
		double phi=0.;
		Xperp=Fperp+(1.+rho*(cos(phi)+I*sin(phi)))*Xperp;


		double complex DeltaTauperpHSA=eq*4./3.*alphas_muf/4./pi*pi*pi*param->f_Bs/3./mb/mBs*(12.*C8eff*mb/mBs*f_phi_perp*Xperp
		+8.*f_phi_perp*integ4-4.*mphi*f_phi_par/(1.-q2/mBs/mBs)/lambda_Bp*integ5);
		
		double complex DeltaTauperpuHSA=eq*4./3.*alphas_muf/4./pi*pi*pi*param->f_Bs/3./mb/mBs*
		(8.*f_phi_perp*integ4u-4.*mphi*f_phi_par/(1.-q2/mBs/mBs)/lambda_Bp*integ5u);

		Tauperpp+=DeltaTauperpWA+DeltaTauperpHSA;
		Tauperpm+=DeltaTauperpWA+DeltaTauperpHSA;

		/*******Using Tauperpu**********/
		Tauperpu+=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(DeltaTauperpuWA+DeltaTauperpuHSA);
		Tauperpu_bar+=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(DeltaTauperpuWA+DeltaTauperpuHSA);

		double lambda=pow(mBs,4.)+pow(mphi,4.)+q2*q2-2.*(mBs*mBs*mphi*mphi+mphi*mphi*q2+mBs*mBs*q2);
	
		double complex N=param->Vtb*conj(param->Vts)*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)/pow(mBs,3.)*q2*sqrt(lambda)*beta_l);
		
		double complex Nbar=conj(param->Vtb)*param->Vts*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)/pow(mBs,3.)*q2*sqrt(lambda)*beta_l);
		
		/* hadronic uncertainties */
		double ALperp_err=1.+param->Bstophilow_ALperp_err_noq2+q2/6.*param->Bstophilow_ALperp_err_q2;
		double ARperp_err=1.+param->Bstophilow_ARperp_err_noq2+q2/6.*param->Bstophilow_ARperp_err_q2;
		double ALpar_err=1.+param->Bstophilow_ALpar_err_noq2+q2/6.*param->Bstophilow_ALpar_err_q2;
		double ARpar_err=1.+param->Bstophilow_ARpar_err_noq2+q2/6.*param->Bstophilow_ARpar_err_q2;
		double AL0_err=1.+param->Bstophilow_AL0_err_noq2+q2/6.*param->Bstophilow_AL0_err_q2;
		double AR0_err=1.+param->Bstophilow_AR0_err_noq2+q2/6.*param->Bstophilow_AR0_err_q2;
		double At_err=1.+param->Bstophilow_At_err_noq2+q2/6.*param->Bstophilow_At_err_q2;
		double AS_err=1.+param->Bstophilow_AS_err_noq2+q2/6.*param->Bstophilow_AS_err_q2;
		double ALperp_bar_err=1.+param->Bstophilow_ALperp_err_noq2+q2/6.*param->Bstophilow_ALperp_err_q2;
		double ARperp_bar_err=1.+param->Bstophilow_ARperp_err_noq2+q2/6.*param->Bstophilow_ARperp_err_q2;
		double ALpar_bar_err=1.+param->Bstophilow_ALpar_err_noq2+q2/6.*param->Bstophilow_ALpar_err_q2;
		double ARpar_bar_err=1.+param->Bstophilow_ARpar_err_noq2+q2/6.*param->Bstophilow_ARpar_err_q2;
		double AL0_bar_err=1.+param->Bstophilow_AL0_err_noq2+q2/6.*param->Bstophilow_AL0_err_q2;
		double AR0_bar_err=1.+param->Bstophilow_AR0_err_noq2+q2/6.*param->Bstophilow_AR0_err_q2;
		double At_bar_err=1.+param->Bstophilow_At_err_noq2+q2/6.*param->Bstophilow_At_err_q2;
		double AS_bar_err=1.+param->Bstophilow_AS_err_noq2+q2/6.*param->Bstophilow_AS_err_q2;
	
		ALperp_bar=N*sqrt(2.)*sqrt(lambda)*(((C9+Y+C9p)-(C10+C10p))*V/(mBs+mphi)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff+C7effp)*T1)+sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu)*(ALperp_err);
		
		ARperp_bar=N*sqrt(2.)*sqrt(lambda)*(((C9+Y+C9p)+(C10+C10p))*V/(mBs+mphi)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff+C7effp)*T1)+sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu)*(ARperp_err);
		
		ALpar_bar=-N*sqrt(2.)*(mBs*mBs-mphi*mphi)*(((C9+Y-C9p)-(C10-C10p))*A1/(mBs-mphi)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff-C7effp)*T2)-sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mBs*mBs-mphi*mphi)*2.*E_phi/mBs*(Tauperpm+Tauperpu)*(ALpar_err);
	
		ARpar_bar=-N*sqrt(2.)*(mBs*mBs-mphi*mphi)*(((C9+Y-C9p)+(C10-C10p))*A1/(mBs-mphi)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff-C7effp)*T2)-sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mBs*mBs-mphi*mphi)*2.*E_phi/mBs*(Tauperpm+Tauperpu)*(ARpar_err);
			
		AL0_bar=-N/2./mphi/sqrt(q2)*(((C9+Y-C9p)-(C10-C10p))*(16.*mBs*mphi*mphi*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*(C7eff-C7effp)*(8.*mBs*mphi*mphi/(mBs+mphi)*T23))-N*(mb+alphas_mub/3./pi*DeltaM)/mphi/sqrt(q2)*((mBs*mBs+3.*mphi*mphi-q2)*2.*E_phi/mBs-lambda/(mBs*mBs-mphi*mphi))*(Tauperpm+Tauperpu)*(AL0_err)+N*(mb+alphas_mub/3./pi*DeltaM)/mphi/sqrt(q2)*lambda/(mBs*mBs-mphi*mphi)*(Tauparm+Tauparu)*(AL0_err);
		
		AR0_bar=-N/2./mphi/sqrt(q2)*(((C9+Y-C9p)+(C10-C10p))*(16.*mBs*mphi*mphi*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*(C7eff-C7effp)*(8.*mBs*mphi*mphi/(mBs+mphi)*T23))-N*(mb+alphas_mub/3./pi*DeltaM)/mphi/sqrt(q2)*((mBs*mBs+3.*mphi*mphi-q2)*2.*E_phi/mBs-lambda/(mBs*mBs-mphi*mphi))*(Tauperpm+Tauperpu)*(AR0_err)+N*(mb+alphas_mub/3./pi*DeltaM)/mphi/sqrt(q2)*lambda/(mBs*mBs-mphi*mphi)*(Tauparm+Tauparu)*(AR0_err);
		
		At_bar=N/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/ml*(CQ2-CQ2p)/(mb_mub+param->mass_s))*A0;
	
		AS_bar=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(mb_mub+param->mass_s)*A0;

		
		/************************************************************************************************/
		ALperp=Nbar*sqrt(2.)*sqrt(lambda)*(((conj(C9)+Y+conj(C9p))-conj(C10+C10p))*V/(mBs+mphi)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff+C7effp)*T1)+sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu_bar)*(ALperp_bar_err);

		ARperp=Nbar*sqrt(2.)*sqrt(lambda)*(((conj(C9)+Y+conj(C9p))+conj(C10+C10p))*V/(mBs+mphi)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff+C7effp)*T1)+sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu_bar)*(ARperp_bar_err);

		ALpar=-Nbar*sqrt(2.)*(mBs*mBs-mphi*mphi)*(((conj(C9)+Y-conj(C9p))-conj(C10-C10p))*A1/(mBs-mphi)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff-C7effp)*T2)-sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mBs*mBs-mphi*mphi)*2.*E_phi/mBs*(Tauperpm+Tauperpu_bar)*(ALpar_bar_err);

		ARpar=-Nbar*sqrt(2.)*(mBs*mBs-mphi*mphi)*(((conj(C9)+Y-conj(C9p))+conj(C10-C10p))*A1/(mBs-mphi)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff-C7effp)*T2)-sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mBs*mBs-mphi*mphi)*2.*E_phi/mBs*(Tauperpm+Tauperpu_bar)*(ARpar_bar_err);

		AL0=-Nbar/2./mphi/sqrt(q2)*(((conj(C9)+Y-conj(C9p))-conj(C10-C10p))*(16.*mBs*mphi*mphi*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*conj(C7eff-C7effp)*(8.*mBs*mphi*mphi/(mBs+mphi)*T23))-Nbar*(mb+alphas_mub/3./pi*DeltaM)/mphi/sqrt(q2)*((mBs*mBs+3.*mphi*mphi-q2)*2.*E_phi/mBs-lambda/(mBs*mBs-mphi*mphi))*(Tauperpm+Tauperpu_bar)*(AL0_bar_err)+Nbar*(mb+alphas_mub/3./pi*DeltaM)/mphi/sqrt(q2)*lambda/(mBs*mBs-mphi*mphi)*(Tauparm+Tauparu_bar)*(AL0_bar_err);

		AR0=-Nbar/2./mphi/sqrt(q2)*(((conj(C9)+Y-conj(C9p))+conj(C10-C10p))*(16.*mBs*mphi*mphi*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*conj(C7eff-C7effp)*(8.*mBs*mphi*mphi/(mBs+mphi)*T23))-Nbar*(mb+alphas_mub/3./pi*DeltaM)/mphi/sqrt(q2)*((mBs*mBs+3.*mphi*mphi-q2)*2.*E_phi/mBs-lambda/(mBs*mBs-mphi*mphi))*(Tauperpm+Tauperpu_bar)*(AR0_bar_err)+Nbar*(mb+alphas_mub/3./pi*DeltaM)/mphi/sqrt(q2)*lambda/(mBs*mBs-mphi*mphi)*(Tauparm+Tauparu_bar)*(AR0_bar_err);

		At=Nbar/sqrt(q2)*sqrt(lambda)*(2.*conj(C10-C10p)+q2/ml*conj(CQ2-CQ2p)/(mb_mub+param->mass_s))*A0;

		AS=-2.*Nbar*sqrt(lambda)*conj(CQ1-CQ1p)/(mb_mub+param->mass_s)*A0;

	}
	
	if(q2>LOW_q2)
	{
		// double mb=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);
		/* Loop function are calculated with Pole mass [S04] which should be replaced by PS mass [BF01] */
		double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
		double mb_PS=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */ //mb=>4.64
		double sh=q2/mb_PS/mb_PS;

		double shat=q2/mBs/mBs;

		double lambda_hat=1.+shat*shat+pow(mphi/mBs,4.)-2.*(shat+shat*mphi*mphi/mBs/mBs+mphi*mphi/mBs/mBs);

		double complex A = A_Seidel(q2, mb_PS, mu_b);

		double complex B = B_Seidel(q2, mb_PS, mu_b);

		double complex C = C_Seidel(q2, mu_b);

		double kappa=1.-2.*alphas_mub/3./pi*log(mu_b/param->mass_b);// MSbar mass should be used for mb

		double complex C9eff=Cmub[9]
		+h_bkll(q2,0.,mu_b)*(4./3.*Cmub[1]+Cmub[2]+11./2.*Cmub[3]-2./3.*Cmub[4]+52.*Cmub[5]-32./3.*Cmub[6])
		-1./2.*h_bkll(q2,mb_PS,mu_b)*(7.*Cmub[3]+4./3.*Cmub[4]+76.*Cmub[5]+64./3.*Cmub[6])
		+4./3.*(Cmub[3]+16./3.*Cmub[5]+16./9.*Cmub[6])
		+alphas_mub/4./pi*(Cmub[1]*(B+4.*C)-3.*Cmub[2]*(2.*B-C)-Cmub[8]*F89_bsll(sh))
		+8.*param->mass_c*param->mass_c/q2*((4./9.*Cmub[1]+1./3.*Cmub[2])*(1.+param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts))+2.*Cmub[3]+20.*Cmub[5]);

		double complex C9eff_bar=conj(Cmub[9])
		+h_bkll(q2,0.,mu_b)*(4./3.*Cmub[1]+Cmub[2]+11./2.*Cmub[3]-2./3.*Cmub[4]+52.*Cmub[5]-32./3.*Cmub[6])
		-1./2.*h_bkll(q2,mb_PS,mu_b)*(7.*Cmub[3]+4./3.*Cmub[4]+76.*Cmub[5]+64./3.*Cmub[6])
		+4./3.*(Cmub[3]+16./3.*Cmub[5]+16./9.*Cmub[6])
		+alphas_mub/4./pi*(Cmub[1]*(B+4.*C)-3.*Cmub[2]*(2.*B-C)-Cmub[8]*F89_bsll(sh))
		+8.*param->mass_c*param->mass_c/q2*((4./9.*Cmub[1]+1./3.*Cmub[2])*(1.+param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb))+2.*Cmub[3]+20.*Cmub[5]);


		double complex C7eff=Cmub[7]
		+alphas_mub/4./pi*((Cmub[1]-6.*Cmub[2])*A-Cmub[8]*F87_bsll(sh,log(mu_b/mb_PS)));

		double complex C7eff_bar=conj(Cmub[7])
		+alphas_mub/4./pi*((Cmub[1]-6.*Cmub[2])*A-Cmub[8]*F87_bsll(sh,log(mu_b/mb_PS)));

	    double complex C7effp=Cpb[7];
        double complex C9p=Cpb[9];
        double complex C10p=Cpb[10];

		double complex N=param->Vtb*conj(param->Vts)*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em*beta_l/3./1024./pow(pi,5.)*mBs*shat*sqrt(lambda_hat));

		double complex Nbar=conj(param->Vtb)*param->Vts*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em*beta_l/3./1024./pow(pi,5.)*mBs*shat*sqrt(lambda_hat));

		double complex f_perp=N*mBs*sqrt(2.*lambda_hat)/(1.+mphi/mBs)*V;

		double complex f_par=N*mBs*sqrt(2.)*(1.+mphi/mBs)*A1;

		double complex f_0=N*mBs*((1.-shat-mphi*mphi/mBs/mBs)*pow(1.+mphi/mBs,2.)*A1-lambda_hat*A2)/(2.*mphi/mBs*(1.+mphi/mBs)*sqrt(shat));

		double complex f_perp_bar=Nbar*mBs*sqrt(2.*lambda_hat)/(1.+mphi/mBs)*V;

		double complex f_par_bar=Nbar*mBs*sqrt(2.)*(1.+mphi/mBs)*A1;

		double complex f_0_bar=Nbar*mBs*((1.-shat-mphi*mphi/mBs/mBs)*pow(1.+mphi/mBs,2.)*A1-lambda_hat*A2)/(2.*mphi/mBs*(1.+mphi/mBs)*sqrt(shat));

		/* A and A_bar choosen such that we have agreement with the convention of 1502.05509 & 22210.11995 */
		ALperp_bar_high=(((C9eff+C9p)-(Cmub[10]+C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff+C7effp))*f_perp;
		ARperp_high=(((C9eff+C9p)+(Cmub[10]+C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff+C7effp))*f_perp;

		ALpar_bar_high=-(((C9eff-C9p)-(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff-C7effp))*f_par;
		ARpar_bar_high=-(((C9eff-C9p)+(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff-C7effp))*f_par;

		AL0_bar_high=-(((C9eff-C9p)-(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff-C7effp))*f_0;
		AR0_bar_high=-(((C9eff-C9p)+(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff-C7effp))*f_0;


		/************************************************************************************************/
		ALperp_high=(((C9eff_bar+conj(C9p))-conj(Cmub[10]+C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar+conj(C7effp)))*f_perp_bar;
		ARperp_high=(((C9eff_bar+conj(C9p))+conj(Cmub[10]+C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar+conj(C7effp)))*f_perp_bar;

		ALpar_high=-(((C9eff_bar-conj(C9p))-conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar-conj(C7effp)))*f_par_bar;
		ARpar_high=-(((C9eff_bar-conj(C9p))+conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar-conj(C7effp)))*f_par_bar;

		AL0_high=-(((C9eff_bar-conj(C9p))-conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar-conj(C7effp)))*f_0_bar;
		AR0_high=-(((C9eff_bar-conj(C9p))+conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mBs/shat*(C7eff_bar-conj(C7effp)))*f_0_bar;

		double complex C10=Cmub[10];
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];

		double lambda=pow(mBs,4.)+pow(mphi,4.)+q2*q2-2.*(mBs*mBs*mphi*mphi+mphi*mphi*q2+mBs*mBs*q2);

		/* A and A_bar choosen such that we have agreement with the convention of 1502.05509 & 22210.11995 */
		At_bar_high=N/sqrt(q2)*sqrt(lambda)*( 2.*(C10-C10p)+q2/ml*(CQ2-CQ2p)/(param->mass_b+param->mass_s))*A0;
		AS_bar_high=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(param->mass_b+param->mass_s)*A0;

		/************************************************************************************************/
		At_high=Nbar/sqrt(q2)*sqrt(lambda)*( 2.*conj(C10-C10p)+q2/ml*conj(CQ2-CQ2p)/(param->mass_b+param->mass_s))*A0;
		AS_high=-2.*Nbar*sqrt(lambda)*conj(CQ1-CQ1p)/(param->mass_b+param->mass_s)*A0;

		/* hadronic uncertainties */
		ALperp_high*=1.+param->Bstophihigh_ALperp_err;
		ARperp_high*=1.+param->Bstophihigh_ARperp_err;
		ALpar_high*=1.+param->Bstophihigh_ALpar_err;
		ARpar_high*=1.+param->Bstophihigh_ARpar_err;
		AL0_high*=1.+param->Bstophihigh_AL0_err;
		AR0_high*=1.+param->Bstophihigh_AR0_err;
		At_high*=1.+param->Bstophihigh_At_err;
		AS_high*=1.+param->Bstophihigh_AS_err;
		ALperp_bar_high*=1.+param->Bstophihigh_ALperp_err;
		ARperp_bar_high*=1.+param->Bstophihigh_ARperp_err;
		ALpar_bar_high*=1.+param->Bstophihigh_ALpar_err;
		ARpar_bar_high*=1.+param->Bstophihigh_ARpar_err;
		AL0_bar_high*=1.+param->Bstophihigh_AL0_err;
		AR0_bar_high*=1.+param->Bstophihigh_AR0_err;
		At_bar_high*=1.+param->Bstophihigh_At_err;
		AS_bar_high*=1.+param->Bstophihigh_AS_err;


		if(q2>HIGH_q2)
		{
			ALperp=ALperp_high;
			ARperp=ARperp_high;
			ALpar=ALpar_high;
			ARpar=ARpar_high;
			AL0=AL0_high;
			AR0=AR0_high;
			At=At_high;
			AS=AS_high;

			ALperp_bar=ALperp_bar_high;
			ARperp_bar=ARperp_bar_high;
			ALpar_bar=ALpar_bar_high;
			ARpar_bar=ARpar_bar_high;
			AL0_bar=AL0_bar_high;
			AR0_bar=AR0_bar_high;
			At_bar=At_bar_high;
			AS_bar=AS_bar_high;
		}
		else
		{
			ALperp=ALperp*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ALperp_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ARperp=ARperp*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ARperp_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ALpar=ALpar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ALpar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ARpar=ARpar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ARpar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AL0=AL0*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AL0_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AR0=AR0*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AR0_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			At=At*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+At_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AS=AS*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AS_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);

			ALperp_bar=ALperp_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ALperp_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ARperp_bar=ARperp_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ARperp_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ALpar_bar=ALpar_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ALpar_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			ARpar_bar=ARpar_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+ARpar_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AL0_bar=AL0_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AL0_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AR0_bar=AR0_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AR0_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			At_bar=At_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+At_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
			AS_bar=AS_bar*(HIGH_q2-q2)/(HIGH_q2-LOW_q2)+AS_bar_high*(q2-LOW_q2)/(HIGH_q2-LOW_q2);
		}
	}

	double A02=AL0*conj(AL0)+AR0*conj(AR0);
	double Apar2=ALpar*conj(ALpar)+ARpar*conj(ARpar);
	double Aperp2=ALperp*conj(ALperp)+ARperp*conj(ARperp);
	
	double A02_bar=AL0_bar*conj(AL0_bar)+AR0_bar*conj(AR0_bar);
	double Apar2_bar=ALpar_bar*conj(ALpar_bar)+ARpar_bar*conj(ARpar_bar);
	double Aperp2_bar=ALperp_bar*conj(ALperp_bar)+ARperp_bar*conj(ARperp_bar);
	
	
	double J1s=0.25*(2.+beta_l*beta_l)*(Aperp2 + Apar2) + 4.*ml*ml/q2*creal(ALperp*conj(ARperp)+ALpar*conj(ARpar));

	double J1c=A02 + 4.*ml*ml/q2*(At*conj(At)+2.*creal(AL0*conj(AR0)))+beta_l*beta_l*AS*conj(AS);

	double J2s=0.25*beta_l*beta_l*(Aperp2+Apar2);

	double J2c=-beta_l*beta_l*A02;
	
	double J3=0.5*beta_l*beta_l*(Aperp2-Apar2);
	
	double J4=1./sqrt(2.)*beta_l*beta_l*creal(AL0*conj(ALpar)+AR0*conj(ARpar));
	
	double J5=sqrt(2.)*beta_l*(creal(AL0*conj(ALperp)-AR0*conj(ARperp))-ml/sqrt(q2)*(creal(ALpar*conj(AS)+ARpar*conj(AS))));
	
	double J6s=2.*beta_l*creal(ALpar*conj(ALperp)-ARpar*conj(ARperp));

	double J6c=4.*beta_l*ml/sqrt(q2)*creal(AL0*conj(AS)+AR0*conj(AS));

	double J7=sqrt(2.)*beta_l*(cimag(AL0*conj(ALpar)-AR0*conj(ARpar))+ml/sqrt(q2)*(cimag(ALperp*conj(AS)+ARperp*conj(AS))));

	double J8=1./sqrt(2.)*beta_l*beta_l*cimag(AL0*conj(ALperp)+AR0*conj(ARperp));
	
	double J9=beta_l*beta_l*cimag(conj(ALpar)*ALperp+conj(ARpar)*ARperp);

	double J1s_bar=0.25*(2.+beta_l*beta_l)*(Aperp2_bar + Apar2_bar) + 4.*ml*ml/q2*creal(ALperp_bar*conj(ARperp_bar)+ALpar_bar*conj(ARpar_bar));

	double J1c_bar=A02_bar + 4.*ml*ml/q2*(At_bar*conj(At_bar)+2.*creal(AL0_bar*conj(AR0_bar)))+beta_l*beta_l*AS_bar*conj(AS_bar);

	double J2s_bar=0.25*beta_l*beta_l*(Aperp2_bar+Apar2_bar);

	double J2c_bar=-beta_l*beta_l*A02_bar;
	
	double J3_bar=0.5*beta_l*beta_l*(Aperp2_bar-Apar2_bar);
	
	double J4_bar=1./sqrt(2.)*beta_l*beta_l*creal(AL0_bar*conj(ALpar_bar)+AR0_bar*conj(ARpar_bar));
	
	double J5_bar=sqrt(2.)*beta_l*(creal(AL0_bar*conj(ALperp_bar)-AR0_bar*conj(ARperp_bar))-ml/sqrt(q2)*(creal(ALpar_bar*conj(AS_bar)+ARpar_bar*conj(AS_bar))));
	
	double J6s_bar=2.*beta_l*creal(ALpar_bar*conj(ALperp_bar)-ARpar_bar*conj(ARperp_bar));

	double J6c_bar=4.*beta_l*ml/sqrt(q2)*creal(AL0_bar*conj(AS_bar)+AR0_bar*conj(AS_bar));

	double J7_bar=sqrt(2.)*beta_l*(cimag(AL0_bar*conj(ALpar_bar)-AR0_bar*conj(ARpar_bar))+ml/sqrt(q2)*(cimag(ALperp_bar*conj(AS_bar)+ARperp_bar*conj(AS_bar))));

	double J8_bar=1./sqrt(2.)*beta_l*beta_l*cimag(AL0_bar*conj(ALperp_bar)+AR0_bar*conj(ARperp_bar));
	
	double J9_bar=beta_l*beta_l*cimag(conj(ALpar_bar)*ALperp_bar+conj(ARpar_bar)*ARperp_bar);


	double complex exppIphi=(1.+0.04*I);
	double complex expmIphi=(1.-0.04*I);
	
	double ys=param->ys_Bs;
	

	double AL0_tilde = AL0_bar;
	double ALpar_tilde = ALpar_bar;
	double AR0_tilde = AR0_bar;
	double ARpar_tilde = ARpar_bar;
	double At_tilde = At_bar;


	double ALperp_tilde = -ALperp_bar;
	double ARperp_tilde = -ARperp_bar;
	double AS_tilde = -AS_bar;

	double complex sh1s_A = exppIphi*(ALperp_tilde*conj(ALperp)+ALpar_tilde*conj(ALpar)+ARperp_tilde*conj(ARperp)+ARpar_tilde*conj(ARpar));
	double complex sh1s_B = exppIphi*(ALperp_tilde*conj(ARperp)+ALpar_tilde*conj(ARpar));
	double complex sh1s_C = expmIphi*(ALperp*conj(ARperp_tilde)+ALpar*conj(ARpar_tilde));
	double h1s = (2.+beta_l*beta_l)/2.*creal(sh1s_A) + 4.*ml*ml/q2*creal(sh1s_B + sh1s_C);
	double s1s = (2.+beta_l*beta_l)/2.*cimag(sh1s_A) + 4.*ml*ml/q2*cimag(sh1s_B - sh1s_C);

	double complex sh1c_A = exppIphi*(AL0_tilde*conj(AL0)+AR0_tilde*conj(AR0));
	double complex sh1c_B = exppIphi*At_tilde*conj(At);
	double complex sh1c_C = exppIphi*AL0_tilde*conj(AR0);
	double complex sh1c_D = expmIphi*AL0*conj(AR0_tilde);
	double complex sh1c_E = exppIphi*AS_tilde*conj(AS);
	double h1c = 2.*creal(sh1c_A)+8.*ml*ml/q2*(creal(sh1c_B)+creal(sh1c_C + sh1c_D))+2.*beta_l*beta_l*creal(sh1c_E);
	double s1c = 2.*cimag(sh1c_A)+8.*ml*ml/q2*(cimag(sh1c_B)+cimag(sh1c_C - sh1c_D))+2.*beta_l*beta_l*cimag(sh1c_E);

	double h2s = beta_l*beta_l/2.*creal(sh1s_A);
	double s2s = beta_l*beta_l/2.*cimag(sh1s_A);

	double h2c = -2.*beta_l*beta_l*creal(sh1c_A);
	double s2c = -2.*beta_l*beta_l*cimag(sh1c_A);

	double complex sh3_A = exppIphi*(ALperp_tilde*conj(ALperp)-ALpar_tilde*conj(ALpar)+ARperp_tilde*conj(ARperp)-ARpar_tilde*conj(ARpar));
	double h3 = beta_l*beta_l*creal(sh3_A);
	double s3 = beta_l*beta_l*cimag(sh3_A);


	double complex sh4_A = exppIphi*(AL0_tilde*conj(ALpar)+AR0_tilde*conj(ARpar));
	double complex sh4_B = expmIphi*(AL0*conj(ALpar_tilde)+AR0*conj(ARpar_tilde));
	double h4 = 1./(sqrt(2.))*beta_l*beta_l*creal(sh4_A+sh4_B);
	double s4 = 1./(sqrt(2.))*beta_l*beta_l*cimag(sh4_A-sh4_B);

	double complex sh5_A = exppIphi*(AL0_tilde*conj(ALperp)-AR0_tilde*conj(ARperp));
	double complex sh5_B = expmIphi*(AL0*conj(ALperp_tilde)-AR0*conj(ARpar_tilde));
	double complex sh5_C = exppIphi*(ALpar_tilde*conj(AS)+ARpar_tilde*conj(AS));
	double complex sh5_D = expmIphi*(ALpar*conj(AS_tilde)+ARpar*conj(AS_tilde));
	double h5 = sqrt(2.)*beta_l*( creal(sh5_A + sh5_B) -ml/sqrt(q2)*creal(sh5_C + sh5_D) );
	double s5 = sqrt(2.)*beta_l*( cimag(sh5_A - sh5_B) -ml/sqrt(q2)*creal(sh5_C - sh5_D) );

	double complex sh6s_A = exppIphi*(ALpar_tilde*conj(ALperp)-ARpar_tilde*conj(ARperp));
	double complex sh6s_B = expmIphi*(ALpar*conj(ALperp_tilde)-ARpar*conj(ARperp_tilde));
	double h6s = 2.*beta_l*creal(sh6s_A + sh6s_B);
	double s6s = 2.*beta_l*cimag(sh6s_A - sh6s_B);

	double complex sh6c_A = exppIphi*(AL0_tilde*conj(AS)+AR0_tilde*conj(AS));
	double complex sh6c_B = expmIphi*(AL0*conj(AS_tilde)+AR0*conj(AS_tilde));
	double h6c = 4.*beta_l*ml/sqrt(q2)*creal(sh6c_A + sh6c_B);
	double s6c = 4.*beta_l*ml/sqrt(q2)*cimag(sh6c_A - sh6c_B);

	double complex sh7_A = exppIphi*(AL0_tilde*conj(ALpar)-AR0_tilde*conj(ARpar));
	double complex sh7_B = expmIphi*(AL0*conj(ALpar_tilde)-AR0*conj(ARpar_tilde));
	double complex sh7_C = exppIphi*(ALperp_tilde*conj(AS)+ARperp_tilde*conj(AS));
	double complex sh7_D = expmIphi*(ALperp*conj(AS_tilde)+ARperp*conj(AS_tilde));
	double h7 = sqrt(2.)*beta_l*(cimag(sh7_A + sh7_B)+ml/(sqrt(q2))*cimag(sh7_C + sh7_D));
	double s7 = -sqrt(2.)*beta_l*(creal(sh7_A - sh7_B)+ml/(sqrt(q2))*creal(sh7_C - sh7_D));

	double complex sh8_A = exppIphi*(AL0_tilde*conj(ALperp)+AR0_tilde*conj(ARperp));
	double complex sh8_B = expmIphi*(AL0*conj(ALperp_tilde)+AR0*conj(ARperp_tilde));
	double h8 = 1./sqrt(q2)*beta_l*beta_l*cimag(sh8_A + sh8_B);
	double s8 = -1./sqrt(q2)*beta_l*beta_l*creal(sh8_A - sh8_B);

	double complex sh9_A = exppIphi*(ALpar_tilde*conj(ALperp)+ARpar_tilde*conj(ARperp));
	double complex sh9_B = expmIphi*(ALpar*conj(ALperp_tilde)+ARpar*conj(ARperp_tilde));
	double h9 = -beta_l*beta_l*cimag(sh9_A + sh9_B);
	double s9 = beta_l*beta_l*creal(sh9_A - sh9_B);




	double dGamma_bsphill_dq2=3./4.*(2.*(J1s-ys/2.*h1s)+(J1c-ys/2.*h1c)-(2.*(J2s-ys/2.*h2s)+(J2c-ys/2.*h2c))/3.);

	double dGamma_bar_bsphill_dq2=3./4.*(2.*(J1s_bar-ys/2.*h1s)+(J1c_bar-ys/2.*h1c)-(2.*(J2s_bar-ys/2.*h2s)+(J2c_bar-ys/2.*h2c))/3.);
	

	double FL[3],AT2[3],S3[3],S4[3],P4prime[3],S7[3];
	

	FL[0]=-((J2c+J2c_bar)-ys*h2c)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	FL[1]=-(J2c+J2c_bar-ys*h2c);
	FL[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	AT2[0]=(J3+J3_bar-ys*h3)/2./(J2s+J2s_bar-ys*h2s);
	AT2[1]=J3+J3_bar-ys*h3;
	AT2[2]=2.*(J2s+J2s_bar-ys*h2s);

	S3[0]=(J3+J3_bar-ys*h3)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	S3[1]=J3+J3_bar-ys*h3;
	S3[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	S4[0]=((J4+J4_bar)-ys*h4)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	S4[1]=(J4+J4_bar-ys*h4);
	S4[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	P4prime[0]=(J4+J4_bar-ys*h4)/sqrt(fabs(-(J2c+J2c_bar-ys*h2c)*(J2s+J2s_bar-ys*h2s)));
	P4prime[1]=J4+J4_bar-ys*h4;
	P4prime[2]=-(J2c+J2c_bar-ys*h2c);

    S7[0]=((J7+J7_bar)-ys*h7)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
    S7[1]=(J7+J7_bar-ys*h7);
    S7[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	double P6prime[3];
	P6prime[0]=-(J7+J7_bar-ys*h7)/2./sqrt(fabs(-(J2c+J2c_bar-ys*h2c)*(J2s+J2s_bar-ys*h2s)));
	P6prime[1]=-(J7+J7_bar-ys*h7)/2.;
	P6prime[2]=-(J2c+J2c_bar-ys*h2c);

	double ATReCP[3], ATImCP[3], A5[3], AFBCP[3], A8[3], A9[3], P2CP[3], P3CP[3], P5primeCP[3], P8primeCP[3], A6c[3];
	/* From 1106.3283 which includes the beta_l function; ATRe = beta_l*2*P2 */
	ATReCP[0]=beta_l*(J6s-J6s_bar-ys*h6s)/(4.*(J2s+J2s_bar-ys*h2s));
	ATReCP[1]=beta_l*(J6s-J6s_bar-ys*h6s);
	ATReCP[2]=4.*(J2s+J2s_bar-ys*h2s);

	ATImCP[0]=(J9-J9_bar-ys*h9)/(2.*(J2s+J2s_bar-ys*h2s));
	ATImCP[1]=(J9-J9_bar-ys*h9);
	ATImCP[2]=2.*(J2s+J2s_bar-ys*h2s);

	A5[0]=(J5-J5_bar-ys*h5)/dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;
	A5[1]=(J5-J5_bar-ys*h5);
	A5[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	int scalar_switch = 0;// 0 if we want to put scalar to zero in AFB, else 1
	AFBCP[0]=-3./8.*(2.*(J6s-J6s_bar-ys*h6s)+scalar_switch*(J6c-J6c_bar-ys*h6c))/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	AFBCP[1]=-3./8.*(2.*(J6s-J6s_bar-ys*h6s)+scalar_switch*(J6c-J6c_bar-ys*h6c));
	AFBCP[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	A8[0]=(J8-J8_bar-ys*h8)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	A8[1]=(J8-J8_bar-ys*h8);
	A8[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	A9[0]=(J9-J9_bar-ys*h9)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	A9[1]=(J9-J9_bar-ys*h9);
	A9[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	P2CP[0]=(J6s-J6s_bar-ys*h6s)/8./(J2s+J2s_bar-ys*h2s);// Definition from 1303.5794 without beta_l and not from 1202.4266
	P2CP[1]=(J6s-J6s_bar-ys*h6s);
	P2CP[2]=8.*(J2s+J2s_bar-ys*h2s);

	P3CP[0]=-(J9-J9_bar-ys*h9)/(4.*(J2s+J2s_bar-ys*h2s));
	P3CP[1]=-(J9-J9_bar-ys*h9);
	P3CP[2]=4.*(J2s+J2s_bar-ys*h2s);

	P5primeCP[0]=(J5-J5_bar-ys*h5)/2./sqrt(fabs(-(J2c+J2c_bar-ys*h2c)*(J2s+J2s_bar-ys*h2s)));
	P5primeCP[1]=(J5-J5_bar-ys*h5)/2.;
	P5primeCP[2]=-(J2c+J2c_bar-ys*h2c);

	P8primeCP[0]=-(J8-J8_bar-ys*h8)/sqrt(fabs(-(J2c+J2c_bar-ys*h2c)*(J2s+J2s_bar-ys*h2s)));
	P8primeCP[1]=-(J8-J8_bar-ys*h8);
	P8primeCP[2]=-(J2c+J2c_bar-ys*h2c);

	A6c[0]=-(J6c-J6c_bar-ys*h6c)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	A6c[1]=-(J6c-J6c_bar-ys*h6c);
	A6c[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;

	double Q8m[3], Q8p[3], Q9[3], S2s[3];
	Q8m[0]=s8/sqrt(fabs( -2.*(J2c+J2c_bar-ys*h2c)*( 2.*(J2s+J2s_bar-ys*h2s) - (J3+J3_bar-ys*h3) ) ));
	Q8m[1]=s8;
	Q8m[2]=sqrt(fabs( -2.*(J2c+J2c_bar-ys*h2c)*( 2.*(J2s+J2s_bar-ys*h2s) - (J3+J3_bar-ys*h3) ) ));

	Q8p[0]=s8/sqrt(fabs( -2.*(J2c+J2c_bar-ys*h2c)*( 2.*(J2s+J2s_bar-ys*h2s) + (J3+J3_bar-ys*h3) ) ));
	Q8p[1]=s8;
	Q8p[2]=sqrt(fabs( -2.*(J2c+J2c_bar-ys*h2c)*( 2.*(J2s+J2s_bar-ys*h2s) + (J3+J3_bar-ys*h3) ) ));

	Q9[0]=s9/(2.*(J2s+J2s_bar-ys*h2s));
	Q9[1]=s9;
	Q9[2]=2.*(J2s+J2s_bar-ys*h2s);

	S2s[0]=(J2s+J2s_bar-ys*h2s)/(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2);
	S2s[1]=(J2s+J2s_bar-ys*h2s);
	S2s[2]=dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2;
	/*********************************************/

	for(je=0;je<=Nobs_Bsphill;je++) for(ie=0;ie<=2;ie++) obs[je][ie]=0.;
	
	for(ie=0;ie<=2;ie++)
	{
		obs[1][ie]=FL[ie];
		obs[2][ie]=AT2[ie];
		obs[3][ie]=S3[ie];
		obs[4][ie]=S4[ie];
		obs[5][ie]=P4prime[ie];
		obs[6][ie]=S7[ie];
		obs[7][ie]=P6prime[ie];
		obs[8][ie]=ATReCP[ie];
		obs[9][ie]=ATImCP[ie];
		obs[10][ie]=A5[ie];
		obs[11][ie]=AFBCP[ie];
		obs[12][ie]=A8[ie];
		obs[13][ie]=A9[ie];
		obs[14][ie]=P2CP[ie];
		obs[15][ie]=P3CP[ie];
		obs[16][ie]=P5primeCP[ie];
		obs[17][ie]=P8primeCP[ie];
		obs[18][ie]=A6c[ie];
		obs[19][ie]=Q8m[ie];
		obs[20][ie]=Q8p[ie];
		obs[21][ie]=Q9[ie];
		obs[22][ie]=S2s[ie];
		// If you intend to add physical observables, add them above following 22 and change Nobs_Bsphill_phys in `include.c` accordingly
	}
	/****** Outputting non-physical quantities - 20/04/2025 ******/
	// If you intend to add non-physical observables, add them below after last line:
	obs[Nobs_Bsphill_phys+1][0] = V;
	obs[Nobs_Bsphill_phys+2][0] = A0;
	obs[Nobs_Bsphill_phys+3][0] = A1;
	obs[Nobs_Bsphill_phys+4][0] = A12;
	obs[Nobs_Bsphill_phys+5][0] = T1;
	obs[Nobs_Bsphill_phys+6][0] = T2;
	obs[Nobs_Bsphill_phys+7][0] = T23;
	obs[Nobs_Bsphill_phys+8][0] = A2;

	return 1./(1.-ys*ys)*(dGamma_bsphill_dq2+dGamma_bar_bsphill_dq2)/2.;
}

/*----------------------------------------------------------------------*/

double dGamma_Bsphill_dq2_full_calculator(int gen, double q2, double obs[][3], char name[])
/* "container" function scanning the SLHA file "name" and calculating dGamma/dq2(Bs->phi mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(gen,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(gen,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(gen,Cpb,CQpb,mu_W,mu_b,&param);

	return dGamma_Bsphill_dq2_full(gen,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBsphill_full(int gen, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	int ie,je;
	int nmax=10;
	if((smin<0.099)||(smax-smin>10.)) nmax=100;
	double Gamma=0.;
	double s;
	
	double obs_num[Nobs_Bsphill+1],obs_den[Nobs_Bsphill+1];
	for(je=0;je<=Nobs_Bsphill;je++) obs_num[je]=obs_den[je]=0.;

	/****** Outputting non-physical quantities - 20/04/2025 ******/
	double obs_min[Nobs_Bsphill+1],obs_max[Nobs_Bsphill+1];
	for(je=0;je<=Nobs_Bsphill;je++) obs_max[je]=obs_min[je]=0.;
	/*************************************************************/

	obs[0]=0.; /* none */
	obs[1]=0.; /* integrated FL */
	obs[2]=0.; /* integrated AT2 */
	obs[3]=0.; /* integrated S3 */
	obs[4]=0.; /* integrated S4 */
	obs[5]=0.; /* integrated P4prime */
	obs[6]=0.; /* integrated S7 */
	obs[7]=0.; /* integrated P6prime */
	obs[8]=0.; /* integrated ATReCP */
	obs[9]=0.; /* integrated ATReIm */
	obs[10]=0.; /* integrated A5 */
	obs[11]=0.; /* integrated AFBCP */
	obs[12]=0.; /* integrated A8 */
	obs[13]=0.; /* integrated A9 */
	obs[14]=0.; /* integrated P2CP */
	obs[15]=0.; /* integrated P3CP */
	obs[16]=0.; /* integrated P5primeCP */
	obs[17]=0.; /* integrated P8primeCP */
	obs[18]=0.; /* integrated A6c */
	obs[19]=0.; /* integrated Q8m */
	obs[20]=0.; /* integrated Q8p */
	obs[21]=0.; /* integrated Q9 */
	obs[22]=0.; /* integrated S2s */
	/****** Outputting non-physical quantities - 20/04/2025 ******/
	obs[Nobs_Bsphill_phys+1]=0.; /* non-integrated (averaged) V */
	obs[Nobs_Bsphill_phys+2]=0.; /* non-integrated (averaged) A0 */
	obs[Nobs_Bsphill_phys+3]=0.; /* non-integrated (averaged) A1 */
	obs[Nobs_Bsphill_phys+4]=0.; /* non-integrated (averaged) A12 */
	obs[Nobs_Bsphill_phys+5]=0.; /* non-integrated (averaged) T1 */
	obs[Nobs_Bsphill_phys+6]=0.; /* non-integrated (averaged) T2 */
	obs[Nobs_Bsphill_phys+7]=0.; /* non-integrated (averaged) T23 */
	obs[Nobs_Bsphill_phys+8]=0.; /* non-integrated (averaged) A2 */
	/************************************************************/

	double dobs[Nobs_Bsphill+1][3];
	double s0m,s0p,s0;
		
	s0=s0m=0.;
	s0p=1.;
	
	double h=(smax-smin)/nmax;	
	s=smin;
	Gamma=dGamma_Bsphill_dq2_full(gen,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_Bsphill;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}
	/****** Outputting non-physical quantities ******/
	for(je=1;je<=Nobs_Bsphill;je++){
		obs_min[je] = dobs[je][0];
	}
	/*************************************************************/

	for(ie=1;ie<nmax;ie++)	
	{
		s+=h;

		Gamma+=4.*dGamma_Bsphill_dq2_full(gen,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_Bsphill;je++) 
		{
			obs_num[je]+=4.*dobs[je][1];
			obs_den[je]+=4.*dobs[je][2];
		}

		Gamma+=2.*dGamma_Bsphill_dq2_full(gen,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_Bsphill;je++) 
		{
			obs_num[je]+=2.*dobs[je][1];
			obs_den[je]+=2.*dobs[je][2];
		}
	}
	
	s=smax;
	Gamma+=4.*dGamma_Bsphill_dq2_full(gen,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_Bsphill;je++) 
	{
		obs_num[je]+=4.*dobs[je][1];
		obs_den[je]+=4.*dobs[je][2];
	}	
	Gamma+=dGamma_Bsphill_dq2_full(gen,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_Bsphill;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}	
	
	Gamma*=h/6.;
	for(je=1;je<=Nobs_Bsphill;je++) 
	{
		obs_num[je]*=h/6.;
		obs_den[je]*=h/6.;
	}
	/****** Outputting non-physical quantities ******/
	for(je=1;je<=Nobs_Bsphill;je++){
		obs_max[je] = dobs[je][0];
	}
	/*************************************************************/

	obs[0]=s0;
	for(je=1;je<=Nobs_Bsphill;je++) 
	if(je==5 || je==7 || je==16 || je==17)
	{
		obs[je]=obs_num[je]/sqrt(fabs(obs_den[je]*obs_den[2]/2.));
	}
	/****** Outputting non-physical quantities ******/
	else if (je>Nobs_Bsphill_phys)// Catching the non-physical observables (last lines appearing in `dGamma_Bsphill_dq2_full` function)
	{
		obs[je] = (obs_min[je] + obs_max[je]) / 2.;
	}
	/*************************************************************/
	else obs[je]=obs_num[je]/obs_den[je];
	
	for(je=1;je<=Nobs_Bsphill;je++) if(fabs(obs[je])<1.e-15) obs[je]=0.;
	
	return param->life_Bs/hbar*Gamma;
}

/*----------------------------------------------------------------------*/

double BRBsphill_lowq2_full(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBsphill_full(gen,1.,6.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBsphill_highq2_full(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBsphill_full(gen,14.18,16.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_Bsphill_lowq2_full_calculator(char name[], double obs[])
/* "container" function scanning the SLHA file "name" and calculating BR(Bs->phi mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
			
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBsphill_lowq2_full(2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_Bsphill_highq2_full_calculator(char name[], double obs[])
/* "container" function scanning the SLHA file "name" and calculating BR(Bs->phi mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBsphill_highq2_full(2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/
/*---------------------------- WRAPPER ---------------------------------*/
/*----------------------------------------------------------------------*/

double dGamma_Bsphill_dq2(int gen, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return dGamma_Bsphill_dq2_full(gen,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return dGamma_Bsphill_dq2_soft(gen,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double dGamma_Bsphill_dq2_calculator(int gen, double q2, double obs[][3], char name[])
{
	return dGamma_Bsphill_dq2_full_calculator(gen,q2,obs,name);
}

/*----------------------------------------------------------------------*/

double BRBsphill(int gen, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return BRBsphill_full(gen,smin,smax,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return BRBsphill_soft(gen,smin,smax,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBsphill_lowq2(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return BRBsphill_lowq2_full(gen,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return BRBsphill_lowq2_soft(gen,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBsphill_highq2(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return BRBsphill_highq2_full(gen,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return BRBsphill_highq2_soft(gen,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_Bsphill_lowq2_calculator(char name[],double obs[])
{
	return BRobs_Bsphill_lowq2_full_calculator(name,obs);
}

/*----------------------------------------------------------------------*/

double BRobs_Bsphill_highq2_calculator(char name[],double obs[])
{
	return BRobs_Bsphill_highq2_full_calculator(name,obs);
}
