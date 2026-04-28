#include "include.h"

double dGamma_LbLll_dq2(int gen, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex Cpb[],  struct parameters* param, double mu_b)
{
    double ml;
    if(gen==1) ml=param->mass_e;
    else if(gen==3) ml=param->mass_tau;
    else ml=param->mass_mu;

    double mLb,mL,mB,mK;

    mLb=param->m_Lambdab;
    mL=param->m_Lambda;
    mB=param->m_B;
    mK=param->m_K;

    double mc=mc_pole_1loop(param);
    double mbpole = mb_pole_1loop(param);

    int ie;

	double beta_l=sqrt(1.-4.*ml*ml/q2);

    double alpha_em=1./133.;

    double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);

    double mu_f=sqrt(mu_b*0.5);

    double complex Cmub[11];
    for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];

    double tau_plus=pow(mB+mK,2.);

    double tau_0=pow(mLb-mL,2.);

    double s_plus = pow(mLb+mL,2.)-q2;
    double s_minus = pow(mLb-mL,2.)-q2;

	double lambda = pow(mLb,4.)+pow(mL,4.)+pow(q2,2.)-2.*(mLb*mLb*mL*mL+mLb*mLb*q2+mL*mL*q2);

	double complex N=param->Gfermi*param->Vtb*conj(param->Vts)*alpha_em*sqrt(q2*sqrt(lambda)*beta_l/(3.*2048.*pow(mLb,3.)*pow(pi,5.)));
	double complex N_bar=param->Gfermi*conj(param->Vtb)*param->Vts*alpha_em*sqrt(q2*sqrt(lambda)*beta_l/(3.*2048.*pow(mLb,3.)*pow(pi,5.)));

	double z= (sqrt(tau_plus -q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));

	double fperp = (param->a0_HO_fperp_LbLll+param->a1_HO_fperp_LbLll*z+param->a2_HO_fperp_LbLll*z*z)/(1.-q2/pow(param->mpole_fh_LbLll,2));
	double hperp = (param->a0_HO_hperp_LbLll+param->a1_HO_hperp_LbLll*z+param->a2_HO_hperp_LbLll*z*z)/(1.-q2/pow(param->mpole_fh_LbLll,2.));
	double gperp = (param->a0_HO_gpp_LbLll+param->a1_HO_gperp_LbLll*z+param->a2_HO_gperp_LbLll*z*z)/(1.-q2/pow(param->mpole_g_LbLll,2.));
	double fplus = (param->a0_HO_fplus_LbLll+param->a1_HO_fplus_LbLll*z+param->a2_HO_fplus_LbLll*z*z)/(1.-q2/pow(param->mpole_fh_LbLll,2.));
	double hplus = (param->a0_HO_hplus_LbLll+param->a1_HO_hplus_LbLll*z+param->a2_HO_hplus_LbLll*z*z)/(1.-q2/pow(param->mpole_fh_LbLll,2.));
	double gplus = (param->a0_HO_gpp_LbLll+param->a1_HO_gplus_LbLll*z+param->a2_HO_gplus_LbLll*z*z)/(1.-q2/pow(param->mpole_g_LbLll,2.));
	double htildeperp = (param->a0_HO_htildepp_LbLll+param->a1_HO_htildeperp_LbLll*z+param->a2_HO_htildeperp_LbLll*z*z)/(1.-q2/pow(param->mpole_g_LbLll,2.));
	double htildeplus = (param->a0_HO_htildepp_LbLll+param->a1_HO_htildeplus_LbLll*z+param->a2_HO_htildeplus_LbLll*z*z)/(1.-q2/pow(param->mpole_g_LbLll,2.));

	double HVplus = -fperp*sqrt(2.*s_minus);
	double HTplus = hperp*(mLb+mL)*sqrt(2.*s_minus);
	double HAplus = -gperp*sqrt(2.*s_plus);
	double HT5plus = -htildeperp*(mLb-mL)*sqrt(2.*s_plus);
	double HV0 = fplus*(mLb+mL)*sqrt(s_minus/q2);
	double HT0 = -hplus*sqrt(q2*s_minus);
	double HA0 = gplus*(mLb-mL)*sqrt(s_plus/q2);
	double HT50 = htildeplus*sqrt(q2*s_plus);

	double mb_mub=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);            
	
	double complex ALperp1 = sqrt(2.)*N*(((Cmub[9]-Cmub[10])+(Cpb[9]-Cpb[10]))*HVplus-2.*mb_mub*(Cmub[7]+Cpb[7])/q2*HTplus);
	double complex ARperp1 = sqrt(2.)*N*(((Cmub[9]+Cmub[10])+(Cpb[9]+Cpb[10]))*HVplus-2.*mb_mub*(Cmub[7]+Cpb[7])/q2*HTplus);
	double complex ALpar1 = -sqrt(2.)*N*(((Cmub[9]-Cmub[10])-(Cpb[9]-Cpb[10]))*HAplus+2.*mb_mub*(Cmub[7]-Cpb[7])/q2*HT5plus);
	double complex ARpar1 = -sqrt(2.)*N*(((Cmub[9]+Cmub[10])-(Cpb[9]+Cpb[10]))*HAplus+2.*mb_mub*(Cmub[7]-Cpb[7])/q2*HT5plus);
	double complex ALperp0 = sqrt(2.)*N*(((Cmub[9]-Cmub[10])+(Cpb[9]-Cpb[10]))*HV0-2.*mb_mub*(Cmub[7]+Cpb[7])/q2*HT0);
	double complex ARperp0 = sqrt(2.)*N*(((Cmub[9]+Cmub[10])+(Cpb[9]+Cpb[10]))*HV0-2.*mb_mub*(Cmub[7]+Cpb[7])/q2*HT0);
	double complex ALpar0 = -sqrt(2.)*N*(((Cmub[9]-Cmub[10])-(Cpb[9]-Cpb[10]))*HA0+2.*mb_mub*(Cmub[7]-Cpb[7])/q2*HT50);
	double complex ARpar0 = -sqrt(2.)*N*(((Cmub[9]+Cmub[10])-(Cpb[9]+Cpb[10]))*HA0+2*mb_mub*(Cmub[7]-Cpb[7])/q2*HT50);

	double complex ALperp1_bar = sqrt(2.)*N_bar*(((Cmub[9]-Cmub[10])+(Cpb[9]-Cpb[10]))*HVplus-2.*mb_mub*(Cmub[7]+Cpb[7])/q2*HTplus);
	double complex ARperp1_bar = sqrt(2.)*N_bar*(((Cmub[9]+Cmub[10])+(Cpb[9]+Cpb[10]))*HVplus-2.*mb_mub*(Cmub[7]+Cpb[7])/q2*HTplus);
	double complex ALpar1_bar = -sqrt(2.)*N_bar*(((Cmub[9]-Cmub[10])-(Cpb[9]-Cpb[10]))*HAplus+2.*mb_mub*(Cmub[7]-Cpb[7])/q2*HT5plus);
	double complex ARpar1_bar = -sqrt(2.)*N_bar*(((Cmub[9]+Cmub[10])-(Cpb[9]+Cpb[10]))*HAplus+2.*mb_mub*(Cmub[7]-Cpb[7])/q2*HT5plus);
	double complex ALperp0_bar = sqrt(2.)*N_bar*(((Cmub[9]-Cmub[10])+(Cpb[9]-Cpb[10]))*HV0-2.*mb_mub*(Cmub[7]+Cpb[7])/q2*HT0);
	double complex ARperp0_bar = sqrt(2.)*N_bar*(((Cmub[9]+Cmub[10])+(Cpb[9]+Cpb[10]))*HV0-2.*mb_mub*(Cmub[7]+Cpb[7])/q2*HT0);
	double complex ALpar0_bar = -sqrt(2.)*N_bar*(((Cmub[9]-Cmub[10])-(Cpb[9]-Cpb[10]))*HA0+2.*mb_mub*(Cmub[7]-Cpb[7])/q2*HT50);
	double complex ARpar0_bar = -sqrt(2.)*N_bar*(((Cmub[9]+Cmub[10])-(Cpb[9]+Cpb[10]))*HA0+2*mb_mub*(Cmub[7]-Cpb[7])/q2*HT50);
	
	/* see 1710.00746 for lepton mass dependance: K1ss = K1, K1cc = K2, K1c = K3, K2ss = K4, K2cc = K5, K2c = K6 -- A_t kept at 0 */
	double K1ss = 	0.25*(ALpar1*conj(ALpar1)+ALperp1*conj(ALperp1)+ARpar1*conj(ARpar1)+ARperp1*conj(ARperp1))
					+0.25*(1.+beta_l*beta_l)*(ALpar0*conj(ALpar0)+ALperp0*conj(ALperp0)+ARpar0*conj(ARpar0)+ARperp0*conj(ARperp0))
					+0.5*(1.-beta_l*beta_l)*creal(ARpar1*conj(ALpar1)+ARperp1*conj(ALperp1)+ARpar0*conj(ALpar0)+ARperp0*conj(ALperp0));
					
	double K1cc = 	0.25*(1.+beta_l*beta_l)*(ARpar1*conj(ARpar1)+ARperp1*conj(ARperp1)+ALpar1*conj(ALpar1)+ALperp1*conj(ALperp1))
					+0.25*(1.-beta_l*beta_l)*(ARpar0*conj(ARpar0)+ARperp0*conj(ARperp0)+ALpar0*conj(ALpar0)+ALperp0*conj(ALperp0))
					+0.5*(1.-beta_l*beta_l)*creal(ARpar1*conj(ALpar1)+ARperp1*conj(ALperp1)+ARpar0*conj(ALpar0)+ARperp0*conj(ALperp0));
					
	double K1c = -beta_l*creal(ARperp1*conj(ARpar1)-ALperp1*conj(ALpar1));
	
	double K2ss = 	0.5*param->alphaL_LbLll*creal(ARperp1*conj(ARpar1)+ALperp1*conj(ALpar1))
					+0.5*param->alphaL_LbLll*(1.+beta_l*beta_l)*creal(ARperp0*conj(ARpar0)+ALperp0*conj(ALpar0))
					+0.5*param->alphaL_LbLll*(1.-beta_l*beta_l)*creal(ARperp1*conj(ALpar1)+ARpar1*conj(ALperp1)+ARperp0*conj(ALpar0)+ARpar0*conj(ALperp0));
	
	double K2cc = 	0.5*param->alphaL_LbLll*(1.+beta_l*beta_l)*creal(ARperp1*conj(ARpar1)+ALperp1*conj(ALpar1))
					+0.5*param->alphaL_LbLll*(1.-beta_l*beta_l)*creal(ARpar0*conj(ARperp0)+ALpar0*conj(ALperp0))
					+0.5*param->alphaL_LbLll*(1.-beta_l*beta_l)*creal(ARperp1*conj(ALpar1)+ARpar1*conj(ALperp1)+ARperp0*conj(ALpar0)+ARpar0*conj(ALperp0));
	
	double K2c = -0.5*param->alphaL_LbLll*beta_l*(ARpar1*conj(ARpar1)+ARperp1*conj(ARperp1)-ALpar1*conj(ALpar1)-ALperp1*conj(ALperp1));


	double K1ss_bar = 	0.25*(ALpar1_bar*conj(ALpar1_bar)+ALperp1_bar*conj(ALperp1_bar)+ARpar1_bar*conj(ARpar1_bar)+ARperp1_bar*conj(ARperp1_bar))
						+0.25*(1.+beta_l*beta_l)*(ALpar0_bar*conj(ALpar0_bar)+ALperp0_bar*conj(ALperp0_bar)+ARpar0_bar*conj(ARpar0_bar)+ARperp0_bar*conj(ARperp0_bar))
						+0.5*(1.-beta_l*beta_l)*creal(ARpar1_bar*conj(ALpar1_bar)+ARperp1_bar*conj(ALperp1_bar)+ARpar0_bar*conj(ALpar0_bar)+ARperp0_bar*conj(ALperp0_bar));
					
	double K1cc_bar = 	0.25*(1.+beta_l*beta_l)*(ARpar1_bar*conj(ARpar1_bar)+ARperp1_bar*conj(ARperp1_bar)+ALpar1_bar*conj(ALpar1_bar)+ALperp1_bar*conj(ALperp1_bar))
						+0.25*(1.-beta_l*beta_l)*(ARpar0_bar*conj(ARpar0_bar)+ARperp0_bar*conj(ARperp0_bar)+ALpar0_bar*conj(ALpar0_bar)+ALperp0_bar*conj(ALperp0_bar))
						+0.5*(1.-beta_l*beta_l)*creal(ARpar1_bar*conj(ALpar1_bar)+ARperp1_bar*conj(ALperp1_bar)+ARpar0_bar*conj(ALpar0_bar)+ARperp0_bar*conj(ALperp0_bar));
					
	double K1c_bar = -beta_l*creal(ARperp1_bar*conj(ARpar1_bar)-ALperp1_bar*conj(ALpar1_bar));
	
	double K2ss_bar = 	0.5*param->alphaL_LbLll*creal(ARperp1_bar*conj(ARpar1_bar)+ALperp1_bar*conj(ALpar1_bar))
						+0.5*param->alphaL_LbLll*(1.+beta_l*beta_l)*creal(ARperp0_bar*conj(ARpar0_bar)+ALperp0_bar*conj(ALpar0_bar))
						+0.5*param->alphaL_LbLll*(1.-beta_l*beta_l)*creal(ARperp1_bar*conj(ALpar1_bar)+ARpar1_bar*conj(ALperp1_bar)+ARperp0_bar*conj(ALpar0_bar)+ARpar0_bar*conj(ALperp0_bar));
	
	double K2cc_bar = 	0.5*param->alphaL_LbLll*(1.+beta_l*beta_l)*creal(ARperp1_bar*conj(ARpar1_bar)+ALperp1_bar*conj(ALpar1_bar))
						+0.5*param->alphaL_LbLll*(1.-beta_l*beta_l)*creal(ARpar0_bar*conj(ARperp0_bar)+ALpar0_bar*conj(ALperp0_bar))
						+0.5*param->alphaL_LbLll*(1.-beta_l*beta_l)*creal(ARperp1_bar*conj(ALpar1_bar)+ARpar1_bar*conj(ALperp1_bar)+ARperp0_bar*conj(ALpar0_bar)+ARpar0_bar*conj(ALperp0_bar));
	
	double K2c_bar = -0.5*param->alphaL_LbLll*beta_l*(ARpar1_bar*conj(ARpar1_bar)+ARperp1_bar*conj(ARperp1_bar)-ALpar1_bar*conj(ALpar1_bar)-ALperp1_bar*conj(ALperp1_bar));


	double dGamma_dq2=0;

	for(int je=0;je<=Nobs_LbLll;je++) for(int ie=0;ie<=2;ie++) obs[je][ie]=0.;

	dGamma_dq2 = (2.*K1ss+2.*K1ss_bar+K1cc+K1cc_bar)/2.;

	obs[1][0] = 3./4.*(K1c+K1c_bar)/dGamma_dq2; // AlFB
	obs[1][1] = 3./4.*(K1c+K1c_bar);
	obs[1][2] = dGamma_dq2;

	obs[2][0] = 1./4.*(2.*K2ss+2.*K2ss_bar+K2cc+K2cc_bar)/dGamma_dq2; //AhFB
	obs[2][1] = 1./4.*(2.*K2ss+2.*K2ss_bar+K2cc+K2cc_bar);
	obs[2][2] = dGamma_dq2;

	obs[3][0] = 3./8.*(K2c+K2c_bar)/dGamma_dq2; //AlhFB
	obs[3][1] = 3./8.*(K2c+K2c_bar);
	obs[3][2] = dGamma_dq2;

	obs[4][0] = (2.*K1ss+2.*K1ss_bar-K1cc-K1cc_bar)/2./dGamma_dq2; // FL
	obs[4][1] = (2.*K1ss+2.*K1ss_bar-K1cc-K1cc_bar)/2.;
	obs[4][2] = dGamma_dq2;

	obs[5][0] = (K1cc+K1cc_bar)/dGamma_dq2; //FT
	obs[5][1] = (K1cc+K1cc_bar);
	obs[5][2] = dGamma_dq2;

    return dGamma_dq2;
}

double dGamma_LbLll_dq2_calculator(int gen, double q2, double obs[][3],char name[])
/* "container" function scanning the SLHA file "name" and calculating dGamma/dq2(Lb->L mu+ mu-) */
{
    double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
    double complex CQ0b[3],CQ1b[3],CQpb[3];

    struct parameters param;

    Init_param(&param);

    if(!Les_Houches_Reader(name,&param)) return 0.;

    double mu_W=param.mass_W;
    double mu_b=param.mass_b_pole;

    CW_calculator(gen,C0w,C1w,C2w,mu_W,&param);
    C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
    CQ_calculator(gen,CQ0b,CQ1b,mu_W,mu_b,&param);
    Cprime_calculator(gen,Cpb,CQpb,mu_W,mu_b,&param);

    return dGamma_LbLll_dq2(gen,q2,obs,C0b,C1b,C2b,Cpb,&param,mu_b);
}

double BRLbLll(int gen, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex Cpb[], struct parameters* param, double mu_b)
{
    int ie, je;
    int nmax=10;
    if((smin<0.099)||(smax-smin>10.)) nmax=100;
    double Gamma=0.;
    double s;

    double obs_num[Nobs_LbLll+1],obs_den[Nobs_LbLll+1];
    for(je=0;je<=Nobs_LbLll;je++) obs_num[je]=obs_den[je]=0.;

	obs[1]=0.; /* integrated AlFB */
	obs[2]=0.; /* integrated AhFB */
	obs[3]=0.; /* integrated AlhFB */
	obs[4]=0.; /* integrated FL */
	obs[5]=0.; /* integrated FT */

    double dobs[Nobs_LbLll+1][3];

    double h=(smax-smin)/nmax;
    s=smin;
    Gamma=dGamma_LbLll_dq2(gen,s,dobs,C0b,C1b,C2b,Cpb,param,mu_b);
    for(je=0;je<=Nobs_LbLll;je++)
    {
        obs_num[je]+=dobs[je][1];
        obs_den[je]+=dobs[je][2];
    }

    for(ie=1;ie<nmax;ie++)
    {
        s+=h;

        Gamma+=4.*dGamma_LbLll_dq2(gen,s-h/2.,dobs,C0b,C1b,C2b,Cpb,param,mu_b);
        for(je=0;je<=Nobs_LbLll;je++)
        {
            obs_num[je]+=4.*dobs[je][1];
            obs_den[je]+=4.*dobs[je][2];
        }

        Gamma+=2.*dGamma_LbLll_dq2(gen,s,dobs,C0b,C1b,C2b,Cpb,param,mu_b);
        for(je=0;je<=Nobs_LbLll;je++)
        {
            obs_num[je]+=2.*dobs[je][1];
            obs_den[je]+=2.*dobs[je][2];
        }
    }

    s=smax;
    Gamma+=4.*dGamma_LbLll_dq2(gen,s-h/2.,dobs,C0b,C1b,C2b,Cpb,param,mu_b);
    for(je=0;je<=Nobs_LbLll;je++)
    {
        obs_num[je]+=4.*dobs[je][1];
        obs_den[je]+=4.*dobs[je][2];
    }
    Gamma+=dGamma_LbLll_dq2(gen,s,dobs,C0b,C1b,C2b,Cpb,param,mu_b);
    for(je=0;je<=Nobs_LbLll;je++)
    {
        obs_num[je]+=dobs[je][1];
        obs_den[je]+=dobs[je][2];
    }

    Gamma*=h/6.;
    for(je=0;je<=Nobs_LbLll;je++)
    {
        obs_num[je]*=h/6.;
        obs_den[je]*=h/6.;
    }

    for(je=1;je<=Nobs_LbLll;je++)
    {
        obs[je]=obs_num[je]/obs_den[je];
    }
    for(je=1;je<=Nobs_LbLll;je++) if(fabs(obs[je])<1.e-15) obs[je]=0.;

    return param->life_Lb/hbar*Gamma; // BR
}
