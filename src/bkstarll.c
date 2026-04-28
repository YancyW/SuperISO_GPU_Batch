#include "include.h"



/*----------------------------------------------------------------------*/

double complex h_bkll(double q2, double mq, double mu)
{
	if(mq==0.) return 4./9.*(2./3.+I*pi+log(mu*mu/q2));
	
	double z=4.*mq*mq/q2;
	
	if(z>1.) return -4./9.*(log(mq*mq/mu/mu)-2./3.-z)
	-4./9.*(2.+z)*sqrt(z-1.)*atan(1./sqrt(z-1.));
	
	else return -4./9.*(log(mq*mq/mu/mu)-2./3.-z)
	-4./9.*(2.+z)*sqrt(1.-z)*(log((1.+sqrt(1.-z))/sqrt(z))-I*pi/2.);
}

/*----------------------------------------------------------------------*/

double phi_Kstar(double u, double a1, double a2)
{
	double x=2.*u-1.;
	double C1=3.*x;
	double C2=-1.5+15./2.*x*x;

	return 6.*u*(1.-u)*(1.+a1*C1+a2*C2);
}

/*----------------------------------------------------------------------*/

double complex B0_bkll(double s, double mq)
{
	double epsilon=1.e-10;
	return -2.*csqrt(4.*(mq*mq-I*epsilon)/s-1.)*catan(1./csqrt(4.*(mq*mq-I*epsilon)/s-1.));
}

/*----------------------------------------------------------------------*/

double complex L1_bkll(double complex x)
{
	return clog((x-1.)/x)*clog(1.-x)-pi*pi/6.+CLi2(x/(x-1.));
}

/*----------------------------------------------------------------------*/

double complex I1_bkll(double u, double mq, double q2, struct parameters* param)
{
	if(mq==0.) return 1.;
	
	double epsilon=1.e-10;

	double complex xp=0.5+csqrt(0.25-(mq*mq-I*epsilon)/((1.-u)*param->m_Bd*param->m_Bd+u*q2));
	double complex xm=0.5-csqrt(0.25-(mq*mq-I*epsilon)/((1.-u)*param->m_Bd*param->m_Bd+u*q2));
	double complex yp=0.5+csqrt(0.25-(mq*mq-I*epsilon)/q2);
	double complex ym=0.5-csqrt(0.25-(mq*mq-I*epsilon)/q2);

	return 1.+2.*mq*mq/(1.-u)/(param->m_Bd*param->m_Bd-q2)*(L1_bkll(xp)+L1_bkll(xm)-L1_bkll(yp)-L1_bkll(ym));
}

/*----------------------------------------------------------------------*/

double complex tperp_bkll(double u, double mq, double q2, double E_Kstar, struct parameters* param)
{
	if(q2==0.)
	{
		if (mq==0.) return 4./(1.-u);
		double epsilon=1.e-10;
		double complex xp=0.5+csqrt(0.25-(mq*mq-I*epsilon)/((1.-u)*param->m_Bd*param->m_Bd));
		double complex xm=0.5-csqrt(0.25-(mq*mq-I*epsilon)/((1.-u)*param->m_Bd*param->m_Bd));
		return 4./(1.-u)*(1.+2.*mq*mq/(1.-u)/param->m_Bd/param->m_Bd*(L1_bkll(xp)+L1_bkll(xm)));
	}
	else return 2.*param->m_Bd/(1.-u)/E_Kstar*I1_bkll(u,mq,q2,param)+q2/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(B0_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,mq)-B0_bkll(q2,mq));
}

/*----------------------------------------------------------------------*/

double complex tpar_bkll(double u, double mq, double q2, double E_Kstar, struct parameters* param)
{
	return 2.*param->m_Bd/(1.-u)/E_Kstar*I1_bkll(u,mq,q2,param)+((1.-u)*param->m_Bd*param->m_Bd+u*q2)/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(B0_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,mq)-B0_bkll(q2,mq));
}

/*----------------------------------------------------------------------*/

double complex F27u(double shat, double l)
{
	double z=4./shat;

	double complex A=
	208./243.*l+4.*shat/27./(1.-shat)*(Li2(shat)+log(shat)*log(1.-shat))
	+1./729./(1.-shat)/(1.-shat)*(6.*shat*(29.-47.*shat)*log(shat)+785.-1600.*shat+833.*shat*shat+6.*pi*I*(20.-49.*shat+47.*shat*shat))
	-2./243./pow(1.-shat,3.)*(2.*csqrt(z-1.)*(-4.+9.*shat-15.*shat*shat+4.*shat*shat*shat)*(pi/2.-catan(csqrt(z-1.)))+9.*shat*shat*shat*log(shat)*log(shat)+18.*pi*I*shat*(1.-2.*shat)*log(shat))
	+2.*shat/243./pow(1.-shat,4.)*(36.*cpow(pi/2.-catan(csqrt(z-1.)),2.)+pi*pi*(-4.+9.*shat-9.*shat*shat+3.*shat*shat*shat));
	
	return -6.*A;
}

/*----------------------------------------------------------------------*/

double complex F19u(double shat, double l)
{
	double z=4./shat;
	
	double complex x1=0.5+0.5*I*csqrt(z-1.);
	double complex x2=0.5-0.5*I*csqrt(z-1.);
	double complex x3=0.5+0.5*I/csqrt(z-1.);
	double complex x4=0.5-0.5*I/csqrt(z-1.);

	double complex B=
	8./243./shat*(-2.*(4.-34.*shat-17.*pi*I*shat)*l+8.*shat*pow(2.*l,2.)-17.*shat*log(shat)*2.*l)
	+(2.+shat)*csqrt(z-1.)/729./shat*(48.*2.*l*(pi/2.-catan(csqrt(z-1.)))-18.*pi*clog(z-1.)+3.*I*clog(z-1.)*clog(z-1.)
	-24.*I*CLi2(-x2/x1)-5.*pi*pi*I+6.*I*(-9.*clog(x1)*clog(x1)+clog(x2)*clog(x2)-2.*clog(x4)*clog(x4)+6.*clog(x1)*clog(x2)-4.*clog(x1)*clog(x3)+8.*clog(x1)*clog(x4))
	-12.*pi*(2.*clog(x1)+clog(x3)+clog(x4)))
	-2./243./shat/(1.-shat)*(4.*shat*(-8.+17.*shat)*(Li2(shat)+log(shat)*log(1.-shat))
	+3.*(2.+shat)*(3.-shat)*clog(x2/x1)*clog(x2/x1)+12.*pi*(-6.-shat+shat*shat)*(pi/2.-catan(csqrt(z-1.))))
	+2./2187./shat/(1.-shat)/(1.-shat)*(-18.*shat*(120.-211.*shat+73.*shat*shat)*log(shat)-288.-8.*shat+934.*shat*shat-692.*shat*shat*shat+18.*pi*I*shat*(82.-173.*shat+73.*shat*shat))
	-4./243./shat/pow(1.-shat,3.)*(-2.*csqrt(z-1.)*(4.-3.*shat-18.*shat*shat+16.*shat*shat*shat-5.*pow(shat,4.))*(pi/2.-catan(csqrt(z-1.)))-9.*shat*shat*shat*log(shat)*log(shat)+2.*pi*I*shat*(8.-33.*shat+51.*shat*shat-17.*shat*shat*shat)*log(shat))
	+2./729./shat/pow(1.-shat,4.)*(72.*(3.-8.*shat+2.*shat*shat)*cpow(pi/2.-catan(csqrt(z-1.)),2.)-pi*pi*(54.-53.*shat-286.*shat*shat+612.*pow(shat,3.)-446.*pow(shat,4.)+113.*pow(shat,5.)));
		
	double complex C=-16./81.*(log(shat)-2.*l)+428./243.-64./27.*zeta3+16./81.*pi*I;
		
	return B+4.*C;
}

/*----------------------------------------------------------------------*/

double complex F29u(double shat, double l)
{
	double z=4./shat;

	double complex x1=0.5+0.5*I*csqrt(z-1.);
	double complex x2=0.5-0.5*I*csqrt(z-1.);
	double complex x3=0.5+0.5*I/csqrt(z-1.);
	double complex x4=0.5-0.5*I/csqrt(z-1.);

	double complex B=
	8./243./shat*(-2.*(4.-34.*shat-17.*pi*I*shat)*l+8.*shat*pow(2.*l,2.)-17.*shat*log(shat)*2.*l)
	+(2.+shat)*csqrt(z-1.)/729./shat*(48.*2.*l*(pi/2.-catan(csqrt(z-1.)))-18.*pi*clog(z-1.)+3.*I*clog(z-1.)*clog(z-1.)
	-24.*I*CLi2(-x2/x1)-5.*pi*pi*I+6.*I*(-9.*clog(x1)*clog(x1)+clog(x2)*clog(x2)-2.*clog(x4)*clog(x4)+6.*clog(x1)*clog(x2)-4.*clog(x1)*clog(x3)+8.*clog(x1)*clog(x4))
	-12.*pi*(2.*clog(x1)+clog(x3)+clog(x4)))
	-2./243./shat/(1.-shat)*(4.*shat*(-8.+17.*shat)*(Li2(shat)+log(shat)*log(1.-shat))
	+3.*(2.+shat)*(3.-shat)*clog(x2/x1)*clog(x2/x1)+12.*pi*(-6.-shat+shat*shat)*(pi/2.-catan(csqrt(z-1.))))
	+2./2187./shat/(1.-shat)/(1.-shat)*(-18.*shat*(120.-211.*shat+73.*shat*shat)*log(shat)-288.-8.*shat+934.*shat*shat-692.*shat*shat*shat+18.*pi*I*shat*(82.-173.*shat+73.*shat*shat))
	-4./243./shat/pow(1.-shat,3.)*(-2.*csqrt(z-1.)*(4.-3.*shat-18.*shat*shat+16.*shat*shat*shat-5.*pow(shat,4.))*(pi/2.-catan(csqrt(z-1.)))-9.*shat*shat*shat*log(shat)*log(shat)+2.*pi*I*shat*(8.-33.*shat+51.*shat*shat-17.*shat*shat*shat)*log(shat))
	+2./729./shat/pow(1.-shat,4.)*(72.*(3.-8.*shat+2.*shat*shat)*cpow(pi/2.-catan(csqrt(z-1.)),2.)-pi*pi*(54.-53.*shat-286.*shat*shat+612.*pow(shat,3.)-446.*pow(shat,4.)+113.*pow(shat,5.)));
		
	double complex C=-16./81.*(log(shat)-2.*l)+428./243.-64./27.*zeta3+16./81.*pi*I;	
	return -6.*B+3.*C;
}

/*----------------------------------------------------------------------*/
/* A, B, C functions for loop contributions by Seidel */
/*----------------------------------------------------------------------*/
double complex A_Seidel(double s, double mb, double mu) {

    double shat = s/mb/mb;
    double z = 4./shat;

    // In the limit s -> 0 limit
    if (fabs(s) <= 1e-6) {
    return (1.1426611796982167 + 0.517134593183505 * I) +
    shat * ((-0.3221817635475234 - 0.23271056693257725 * I) +
    shat * ((-0.16999092571433885 + 0.23271056693257727 * I) +
    ((-0.13262865040070346 + 0.6981317007977318 * I) -
    (0.11612865336448869 - 1.1635528346628865 * I) * shat) * shat)) -
    0.8559670781893004 * log(mb/mu) -
    shat * ((0.23868312757201646 - 0.4654211338651545 * I) +
    shat * ((-0.05761316872427984 - 0.46542113386515455 * I) +
    (-0.27983539094650206 - (0.4773662551440329 -
    0.9308422677303091 * I) * shat) * shat) +
    (-0.07407407407407407 - 0.22222222222222224 * shat) * (-shat * shat / 2));
    //     return -104.0 / 243.0 * 2.0 * log(mb / mu);
    }


    /* Taylor approximation near shat = 1 */
    double complex mu_term = -104. / 243. * 2. * log(mb/mu);
    if (fabs(shat - 1.0) < 1e-2) {
        const double complex c0 = (997. + 18. * sqrt(3.) * pi) / 1458. + I * (64./243. * pi);
        const double complex c1 = (215. + 9. * sqrt(3.) * pi) / 1215. + I * (-1./27. * pi);
        const double complex c2 = (95. + 12. * sqrt(3.) * pi) / 2430. + I * (-7./405. * pi);

        return mu_term + c0 + c1 * (1.-shat) + c2 * pow(1.-shat, 2);
    }

    /* In the A expression, for (CLi2(shat)+log(shat)*clog(1.-shat)), the real part is calculated correctly
     * but the imaginary part above the branch cut is not (in principle it should get cancelled).
     * I fixed the CLi2 to use the same branch cut as clog, so it can be used directly */
    double complex Li2log_term = CLi2(shat)+clog(shat)*clog(1.-shat);

    return -104./243.*log(mb*mb/mu/mu)+4.*shat/27./(1.-shat)*Li2log_term
		+1./729./(1.-shat)/(1.-shat)*(6.*shat*(29.-47.*shat)*clog(shat)+785.-1600.*shat+833.*shat*shat+6.*pi*I*(20.-49.*shat+47.*shat*shat))
		-2./243./pow(1.-shat,3.)*(2.*csqrt(z-1.)*(-4.+9.*shat-15.*shat*shat+4.*shat*shat*shat)*(pi/2.-catan(csqrt(z-1.)))+9.*shat*shat*shat*clog(shat)*clog(shat)+18.*pi*I*shat*(1.-2.*shat)*clog(shat))
		+2.*shat/243./pow(1.-shat,4.)*(36.*cpow(pi/2.-catan(csqrt(z-1.)),2.)+pi*pi*(-4.+9.*shat-9.*shat*shat+3.*shat*shat*shat));

}

complex double B_Seidel(double s, double mb, double mu) {
	double shat = s/mb/mb;
	double z = 4./shat;


	double complex x1=0.5+0.5*I*csqrt(z-1.);
	double complex x2=0.5-0.5*I*csqrt(z-1.);
	double complex x3=0.5+0.5*I/csqrt(z-1.);
	double complex x4=0.5-0.5*I/csqrt(z-1.);

	if (fabs(shat - 1.0) < 1e-2) {
		double complex mu_terms =
		8./243./shat * (4.-34.*shat-17.*pi*I*shat) * log(mb*mb/mu/mu)
		+ 17.*shat * 8./243./shat * log(shat) * log(mb*mb/mu/mu)
		+ (2.+shat) * csqrt(z-1.)/729./shat * (-48.) * (pi/2. - catan(csqrt(z-1.))) * log(mb*mb/mu/mu);

		double complex nonmu_terms = -1.2534705628994441 + 3.1545210184193809 * I +
		(-1.1399966466176837 - 1.3704066719362884 * I) * (shat - 1.0) +
		(0.77575942579740349 + 0.59987612809286587 * I) * pow(shat - 1.0, 2);

		double complex result = mu_terms + nonmu_terms;

		return result;

	}


	/* In the A expression, for (CLi2(shat)+log(shat)*clog(1.-shat)), the real part is calculated correctly
	 * but the imaginary part above the branch cut is not (in principle it should get cancelled).
	 * I fixed the CLi2 to use the same branch cut as clog, so it can be used directly */
	double complex Li2log_term = CLi2(shat)+clog(shat)*clog(1.-shat);

	return 8./243./shat*((4.-34.*shat-17.*pi*I*shat)*log(mb*mb/mu/mu)+8.*shat*pow(log(mb*mb/mu/mu),2.)+17.*shat*clog(shat)*log(mb*mb/mu/mu))
	+(2.+shat)*csqrt(z-1.)/729./shat*(-48.*log(mb*mb/mu/mu)*(pi/2.-catan(csqrt(z-1.)))-18.*pi*clog(z-1.)+3.*I*clog(z-1.)*clog(z-1.)
	-24.*I*CLi2(-x2/x1)-5.*pi*pi*I+6.*I*(-9.*clog(x1)*clog(x1)+clog(x2)*clog(x2)-2.*clog(x4)*clog(x4)+6.*clog(x1)*clog(x2)-4.*clog(x1)*clog(x3)+8.*clog(x1)*clog(x4))
	-12.*pi*(2.*clog(x1)+clog(x3)+clog(x4)))
	-2./243./shat/(1.-shat)*(4.*shat*(-8.+17.*shat)*Li2log_term
	+3.*(2.+shat)*(3.-shat)*clog(x2/x1)*clog(x2/x1)+12.*pi*(-6.-shat+shat*shat)*(pi/2.-catan(csqrt(z-1.))))
	+2./2187./shat/(1.-shat)/(1.-shat)*(-18.*shat*(120.-211.*shat+73.*shat*shat)*clog(shat)-288.-8.*shat+934.*shat*shat-692.*shat*shat*shat+18.*pi*I*shat*(82.-173.*shat+73.*shat*shat))
	-4./243./shat/pow(1.-shat,3.)*(-2.*csqrt(z-1.)*(4.-3.*shat-18.*shat*shat+16.*shat*shat*shat-5.*pow(shat,4.))*(pi/2.-catan(csqrt(z-1.)))-9.*shat*shat*shat*clog(shat)*clog(shat)+2.*pi*I*shat*(8.-33.*shat+51.*shat*shat-17.*shat*shat*shat)*clog(shat))
	+2./729./shat/pow(1.-shat,4.)*(72.*(3.-8.*shat+2.*shat*shat)*cpow(pi/2.-catan(csqrt(z-1.)),2.)-pi*pi*(54.-53.*shat-286.*shat*shat+612.*pow(shat,3.)-446.*pow(shat,4.)+113.*pow(shat,5.)));
}

complex double C_Seidel(double s, double mu) {

    return -16./81.*clog(s/mu/mu)+428./243.-64./27.*zeta3+16./81.*pi*I;
}


/*----------------------------------------------------------------------*/
/*------------------------ ISOSPIN ASYMMETRY ---------------------------*/
/*----------------------------------------------------------------------*/

double dAI_BKstarmumu_dq2(double q2, double complex C0b[], double complex C1b[], double complex C2b[], struct parameters* param, double mu_b)
{
	double mc=mc_pole(param);

	int ie;
	double shat=q2/param->m_B/param->m_B;
	
	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);

	double mu_f=sqrt(mu_b*0.5);
	double alphas_muf=alphas_running(mu_f,param->mass_top_pole,param->mass_b_pole,param);
	double eta=alphas_muf/alphas_running(1.,param->mass_top_pole,param->mass_b_pole,param);

	double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
	double mb=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */

	double complex Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_Kstar=(param->m_B*param->m_B+param->m_Kstar0*param->m_Kstar0-q2)/2./param->m_B;

	int nf=5;
	double f_Kstar_perp=param->f_Kstar_perp;
	f_Kstar_perp*=pow(eta,4./3./(11.-2./3.*nf));

	double f_Kstar_par=param->f_Kstar_par;

	/********LCSR+Lattice fit parameters from Barucha,Straub,Zwicky 1503.05534***************/
	double tau_plus=pow(param->m_B+param->m_Kstar0,2.);
	double tau_minus=pow(param->m_B-param->m_Kstar0,2.);
	double tau_0=tau_plus-sqrt((tau_plus-tau_minus)*tau_plus);
	double z_q2tau0=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));
	double z_0tau0=(sqrt(tau_plus)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus)+sqrt(tau_plus-tau_0));
	
	double P_V=1./(1.-q2/pow(param->MV_BKstar,2.));
	double P_A1=1./(1.-q2/pow(param->MA1_BKstar,2.));
	double P_A12=1./(1.-q2/pow(param->MA12_BKstar,2.));

	double V=P_V*(param->a0V_BKstar+param->a1V_BKstar*(z_q2tau0-z_0tau0)+param->a2V_BKstar*(pow(z_q2tau0-z_0tau0,2.)));
	double A1=P_A1*(param->a0A1_BKstar+param->a1A1_BKstar*(z_q2tau0-z_0tau0)+param->a2A1_BKstar*(pow(z_q2tau0-z_0tau0,2.)));
	double A12=P_A12*(param->a0A12_BKstar+param->a1A12_BKstar*(z_q2tau0-z_0tau0)+param->a2A12_BKstar*(pow(z_q2tau0-z_0tau0,2.)));

	double A2=(pow(param->m_B+param->m_Kstar0,2.)*(param->m_B*param->m_B-param->m_Kstar0*param->m_Kstar0-q2)*A1-16.*param->m_B*param->m_Kstar0*param->m_Kstar0*(param->m_B+param->m_Kstar0)*A12)/((pow(param->m_B+param->m_Kstar0,2.)-q2)*(pow(param->m_B-param->m_Kstar0,2.)-q2));
		
	double xi_par=(param->m_B+param->m_Kstar0)/2./E_Kstar*A1-(param->m_B-param->m_Kstar0)/param->m_B*A2;
	double xi_perp=param->m_B/(param->m_B+param->m_Kstar0)*V;

	double complex C7eff=Cmub[7];
	double complex C8eff=Cmub[8];
	double complex C9=Cmub[9];
	double complex C10=Cmub[10];
		
	double complex C1bar=Cmub[1]/2.;
	double complex C2bar=Cmub[2]-Cmub[1]/6.;
	double complex C3bar=Cmub[3]-Cmub[4]/6.+16.*Cmub[5]-8./3.*Cmub[6];
	double complex C4bar=Cmub[4]/2.+8.*Cmub[6];
	double complex C5bar=Cmub[3]-Cmub[4]/6.+4.*Cmub[5]-2./3.*Cmub[6];
	double complex C6bar=Cmub[4]/2.+2.*Cmub[6];
	
	double complex Y=4./3.*Cmub[3]+64./9.*Cmub[5]+64./27.*Cmub[6]
	+h_bkll(q2,mc,mu_b)*(4./3.*Cmub[1]+Cmub[2]+6.*Cmub[3]+60.*Cmub[5])
	+h_bkll(q2,param->mass_b_pole,mu_b)*(-7./2.*Cmub[3]-2./3.*Cmub[4]-38.*Cmub[5]-32./3.*Cmub[6])
	+h_bkll(q2,0.,mu_b)*(-1./2.*Cmub[3]-2./3.*Cmub[4]-8.*Cmub[5]-32./3.*Cmub[6]);

	double complex C90perp=C9+Y+2.*mb*param->m_B/q2*C7eff;
	double complex C90par=C9+Y+2.*mb/param->m_B*C7eff;

	double a1perp=param->a1perp;
	double a2perp=param->a2perp;
	double a1par=param->a1par;
	double a2par=param->a2par;

	a1perp*=pow(eta,4./3.*(4.*1./2.)/(11.-2./3.*nf));
	a2perp*=pow(eta,4./3.*(4.*(1./2.+1./3.))/(11.-2./3.*nf));

	a1par*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
	a2par*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

	double u;

	double lambda_Bp=param->lambda_Bp;
	lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_f,2.))*(1.-2.*1.4);

	double omega0=2.*(param->m_Bd-mb)/3.;
	double complex lambda_Bm=1./(exp(-q2/param->m_Bd/omega0)/omega0*(-Ei(q2/param->m_Bd/omega0)+I*pi));

	double complex integ1,integ2,integ3;
	integ1=integ2=integ3=0.;
	double complex Fperp=0.;
	double complex Xperp=0.;
	double complex Fpar=0.;
	double complex FV;
	double x;
	
	double zeta3A=param->zeta3A;
	double zeta3V=param->zeta3V;
	double wA10=param->wA10;
	double deltatp=param->deltatp;
	double deltatm=param->deltatm;
	
	double phiKstar_perp,phiKstar_par;

	int n1=10;
	int n1sav=n1;
	for(ie=0;ie<=n1;ie++)
	{
		u=(double)ie/n1;
			if(ie==0) n1*=2;
			if(ie==n1){u=0.99;n1*=2;}
		
		phiKstar_perp=phi_Kstar(u,a1perp,a2perp);
		phiKstar_par=phi_Kstar(u,a1par,a2par);
		
		x=(1.-u)*param->m_B*param->m_B+u*q2;
		FV=3./4.*(h_bkll(x,mc,mu_b)*(C2bar+C4bar+C6bar)+h_bkll(x,param->mass_b_pole,mu_b)*(C3bar+C4bar+C6bar)+h_bkll(x,0.,mu_b)*(C3bar+3.*C4bar+3.*C6bar)-8./27.*(C3bar-C5bar-15.*C6bar));
		
		integ1+=phiKstar_perp/((1.-u)+u*shat)*FV/n1;
		
		integ2+=((3./4.*(1.+pow(2.*u-1.,2.))+a1par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1par*(2.*u-1.)+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1par*u+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FV/n1;
	
		integ3+=phiKstar_par*FV/n1;
		
		Fperp+=phiKstar_perp/((1.-u)+u*shat)/3./n1;
		Xperp+=(u<=1.-0.5/param->m_B)*phiKstar_perp/pow((1.-u)+u*shat,2.)/3./n1;
		Fpar+=2.*phiKstar_par/((1.-u)+u*shat)/n1;

		if(ie==0||ie==n1sav) n1=n1sav;
	}
	double rho=0.;
	double phi=0.;
	Xperp=Fperp+(1.+rho*(cos(phi)+I*sin(phi)))*Xperp;
	
	double complex lambda_u=creal(conj(param->Vus)*param->Vub); /* only the real part is contributing */
	double complex lambda_t=creal(conj(param->Vts)*param->Vtb);

	double complex K1upar_a=-lambda_u/lambda_t*(C1bar/3.+C2bar)+(C4bar+C3bar/3.);
	double complex K1dpar_a=C4bar+C3bar/3.;
	
	double complex K1perp_a=-(C6bar+C5bar/3.)*Fperp;
	
	double complex K2uperp_a=-lambda_u/lambda_t*(C1bar/3.+C2bar)+(C4bar+C3bar/3.);
	double complex K2dperp_a=C4bar+C3bar/3.;
	
	double complex K1perp_b=C8eff*mb/param->m_B*4./9.*alphas_mub/4./pi*Xperp;
	//double complex K2perp_b=0.;
	double complex K1par_b=-C8eff*mb/param->m_B*4./9.*alphas_mub/4./pi*Fpar;
			
	double complex K1perp_c=4./9.*alphas_mub/4./pi*2./3.*integ1;
	double complex K2perp_c=-4./9.*alphas_mub/4./pi*1./2.*integ2;
	
	double complex K1par_c=-4./9.*alphas_mub/4./pi*2.*integ3;
	
	double complex K1perp=K1perp_a+K1perp_b+K1perp_c;
	double complex K2uperp=K2uperp_a+K2perp_c;
	double complex K2dperp=K2dperp_a+K2perp_c;
	double complex K1upar=K1upar_a+K1par_b+K1par_c;
	double complex K1dpar=K1dpar_a+K1par_b+K1par_c;
	
	double eu=2./3.;
	double ed=-1./3.;
	
	double complex bd_perp=24.*pi*pi*param->m_B*param->f_B*ed/q2/xi_perp/C90perp*(f_Kstar_perp/param->m_B*K1perp+f_Kstar_par*param->m_Kstar0/6./lambda_Bp/param->m_B*K2dperp/(1.-q2/param->m_B/param->m_B));
	
	double complex bu_perp=24.*pi*pi*param->m_B*param->f_B*eu/q2/xi_perp/C90perp*(f_Kstar_perp/param->m_B*K1perp+f_Kstar_par*param->m_Kstar0/6./lambda_Bp/param->m_B*K2uperp/(1.-q2/param->m_B/param->m_B));
	
	double complex bd_par=24.*pi*pi*param->f_B*ed*param->m_Kstar0/param->m_B/E_Kstar/xi_par/C90par*(f_Kstar_par/3./lambda_Bm*K1dpar);
	
	double complex bu_par=24.*pi*pi*param->f_B*eu*param->m_Kstar0/param->m_B/E_Kstar/xi_par/C90par*(f_Kstar_par/3./lambda_Bm*K1upar);	
	
	double dAI_dq2=creal(bd_perp-bu_perp)*pow(cabs(C90perp),2.)/(pow(cabs(C90perp),2.)+C10*C10)*(1.+0.25/q2*pow(E_Kstar*param->m_B/param->m_Kstar0*xi_par/xi_perp*cabs(C90par/C90perp),2.)*creal(bd_par-bu_par)/creal(bd_perp-bu_perp))/(1.+0.25/q2*pow(E_Kstar*param->m_B/param->m_Kstar0*xi_par/xi_perp,2.)*(pow(cabs(C90par),2.)+C10*C10)/(pow(cabs(C90perp),2.)+C10*C10));
		
	return dAI_dq2;
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu(double smin, double smax, double complex C0b[], double complex C1b[], double complex C2b[], struct parameters* param, double mu_b)
{
	int ie;
	int nmax=10;
	if((smin<0.099)||(smax-smin>10.)) nmax=100;
	double AI=0.;
	double s;
	
	/*AI=dAI_BKstarmumu_dq2(smin,C0b,C1b,C2b,param,mu_b)/2.;
	for(ie=1;ie<nmax;ie++)
	{
		s=smin+(smax-smin)*ie/nmax;
		AI+=dAI_BKstarmumu_dq2(s,C0b,C1b,C2b,param,mu_b);
	}
	AI+=dAI_BKstarmumu_dq2(smax,C0b,C1b,C2b,param,mu_b)/2.;
	AI*=(smax-smin)/nmax;*/
	
	double h=(smax-smin)/nmax;	
	s=smin;
	AI=dAI_BKstarmumu_dq2(smin,C0b,C1b,C2b,param,mu_b);
	for(ie=1;ie<nmax;ie++)
	{
		s+=h;
		AI+=4.*dAI_BKstarmumu_dq2(s-h/2.,C0b,C1b,C2b,param,mu_b);
		AI+=2.*dAI_BKstarmumu_dq2(s,C0b,C1b,C2b,param,mu_b);
	}
	s=smax;
	AI+=4.*dAI_BKstarmumu_dq2(s-h/2.,C0b,C1b,C2b,param,mu_b);
	AI+=dAI_BKstarmumu_dq2(s,C0b,C1b,C2b,param,mu_b);
	
	AI*=h/6.;

	return AI;
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_lowq2(double complex C0b[], double complex C1b[], double complex C2b[], struct parameters* param, double mu_b)
{
	return AI_BKstarmumu(1.,6.,C0b,C1b,C2b,param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_highq2(double complex C0b[], double complex C1b[], double complex C2b[], struct parameters* param, double mu_b)
{
	return AI_BKstarmumu(14.18,16.,C0b,C1b,C2b,param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_zero(double complex C0b[], double complex C1b[], double complex C2b[], struct parameters* param, double mu_b)
{
	double smin=pow(2.*param->mass_mu,2.);
	double smax=pow(param->m_Bd-param->m_Kstar0,2.)*0.999; 

	double stemp=0.;

	while(fabs(1.-smin/smax)>1.e-4)
	{
		stemp=(smax+smin)/2.;
		
		if(dAI_BKstarmumu_dq2(stemp,C0b,C1b,C2b,param,mu_b)>0.) smin=stemp; else smax=stemp;
	}
	return stemp;
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_lowq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating AI(B->Kstar mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);

	return AI_BKstarmumu_lowq2(C0b,C1b,C2b,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_highq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating AI(B->Kstar mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);

	return AI_BKstarmumu_highq2(C0b,C1b,C2b,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_zero_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the zero of AI(B->Kstar mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);

	return AI_BKstarmumu_zero(C0b,C1b,C2b,&param,mu_b);
}

/*----------------------------------------------------------------------*/
/*------------------------------ SOFT ----------------------------------*/
/*----------------------------------------------------------------------*/

double dGamma_BKstarll_dq2_soft(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{	
	double ml;
	if(gen==1) ml=param->mass_e;
	else if(gen==3) ml=param->mass_tau;
	else ml=param->mass_mu;
	
	double mB,mKstar,eq;
	if(charge==0) 
	{
		mB=param->m_Bd;
		mKstar=param->m_Kstar0;
		eq=-1./3.;
	}
	else
	{
		mB=param->m_B;
		mKstar=param->m_Kstar;
		eq=2./3.;
	}

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

	double complex Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_Kstar=(mB*mB+mKstar*mKstar-q2)/2./mB;

	int nf=5;
	double f_Kstar_perp=param->f_Kstar_perp;
	f_Kstar_perp*=pow(eta,4./3./(11.-2./3.*nf));

	double f_Kstar_par=param->f_Kstar_par;
	
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

	/********LCSR+Lattice fit parameters from Barucha,Straub,Zwicky 1503.05534***************/
	if(param->BKstar_FormFactor_choice ==1)
	{
		MB = mB;
		MV = mKstar;
		a0V = param->a0V_BKstar;
		a1V = param->a1V_BKstar;
		a2V = param->a2V_BKstar;
		a0A0 = param->a0A0_BKstar;
		a1A0 = param->a1A0_BKstar;
		a2A0 = param->a2A0_BKstar;
		a0A1 = param->a0A1_BKstar;
		a1A1 = param->a1A1_BKstar;
		a2A1 = param->a2A1_BKstar;
		a0A12 = param->a0A12_BKstar;
		a1A12 = param->a1A12_BKstar;
		a2A12 = param->a2A12_BKstar;
		a0T1 = param->a0T1_BKstar;
		a1T1 = param->a1T1_BKstar;
		a2T1 = param->a2T1_BKstar;
		a0T2 = param->a0T2_BKstar;
		a1T2 = param->a1T2_BKstar;
		a2T2 = param->a2T2_BKstar;
		a0T23 = param->a0T23_BKstar;
		a1T23 = param->a1T23_BKstar;
		a2T23 = param->a2T23_BKstar;

		double P_V=1./(1.-q2/pow(param->MV_BKstar,2.));
		double P_A0=1./(1.-q2/pow(param->MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(param->MA1_BKstar,2.));
		double P_A12=1./(1.-q2/pow(param->MA12_BKstar,2.));
		double P_T1=1./(1.-q2/pow(param->MT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(param->MT2_BKstar,2.));
		double P_T23=1./(1.-q2/pow(param->MT23_BKstar,2.));

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


	/********Fit parameters from Gubernari, Reboud, van Dyk, Virto, 2305.06301***************/
	if(param->BKstar_FormFactor_choice ==2)
	{
		MB = mB;
		MV = mKstar;
		double MVT1_BKstar = 5.415;
		double MA0_BKstar = 5.367;
		double MAT_BKstar = 5.829;
		a0V = param->a0V_BKstar_GRvDV_BSZ;
		a1V = param->a1V_BKstar_GRvDV_BSZ;
		a2V = param->a2V_BKstar_GRvDV_BSZ;
		a0A0 = param->a0A0_BKstar_GRvDV_BSZ;
		a1A0 = param->a1A0_BKstar_GRvDV_BSZ;
		a2A0 = param->a2A0_BKstar_GRvDV_BSZ;
		a0A1 = param->a0A1_BKstar_GRvDV_BSZ;
		a1A1 = param->a1A1_BKstar_GRvDV_BSZ;
		a2A1 = param->a2A1_BKstar_GRvDV_BSZ;
		a0A12 = (MB*MB - MV*MV)/(8.*MB*MV)*param->a0A0_BKstar_GRvDV_BSZ;// Eq.17 in BSZ [1503.05534]
		a1A12 = param->a1A12_BKstar_GRvDV_BSZ;
		a2A12 = param->a2A12_BKstar_GRvDV_BSZ;
		a0T1 = param->a0T1_BKstar_GRvDV_BSZ;
		a1T1 = param->a1T1_BKstar_GRvDV_BSZ;
		a2T1 = param->a2T1_BKstar_GRvDV_BSZ;
		a0T2 = param->a0T1_BKstar_GRvDV_BSZ;// Eq.17 in BSZ [1503.05534]
		a1T2 = param->a1T2_BKstar_GRvDV_BSZ;
		a2T2 = param->a2T2_BKstar_GRvDV_BSZ;
		a0T23 = param->a0T23_BKstar_GRvDV_BSZ;
		a1T23 = param->a1T23_BKstar_GRvDV_BSZ;
		a2T23 = param->a2T23_BKstar_GRvDV_BSZ;

		double P_V=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_A0=1./(1.-q2/pow(MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_A12=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T1=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T23=1./(1.-q2/pow(MAT_BKstar,2.));

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


	/********LCSR+Lattice fit parameters from Gubernari van Dyk Kokulu 1811.00983***************/
	if(param->BKstar_FormFactor_choice == 3)
	{
		MB = mB;
		MV = mKstar;
		double MVT1_BKstar = 5.412;
		double MA0_BKstar = 5.336;
		double MAT_BKstar = 5.829;
		a0V = param->a0V_BKstar_GKvD_LCSR_Lattice;
		a1V = param->a1V_BKstar_GKvD_LCSR_Lattice;
		a2V = param->a2V_BKstar_GKvD_LCSR_Lattice;
		a0A0 = param->a0A0_BKstar_GKvD_LCSR_Lattice;
		a1A0 = param->a1A0_BKstar_GKvD_LCSR_Lattice;
		a2A0 = param->a2A0_BKstar_GKvD_LCSR_Lattice;
		a0A1 = param->a0A1_BKstar_GKvD_LCSR_Lattice;
		a1A1 = param->a1A1_BKstar_GKvD_LCSR_Lattice;
		a2A1 = param->a2A1_BKstar_GKvD_LCSR_Lattice;
		a0A12 = (MB*MB - MV*MV)/(8.*MB*MV)*param->a0A0_BKstar_GKvD_LCSR_Lattice;// Eq.17 in BSZ [1503.05534]
		a1A12 = param->a1A12_BKstar_GKvD_LCSR_Lattice;
		a2A12 = param->a2A12_BKstar_GKvD_LCSR_Lattice;
		a0T1 = param->a0T1_BKstar_GKvD_LCSR_Lattice;
		a1T1 = param->a1T1_BKstar_GKvD_LCSR_Lattice;
		a2T1 = param->a2T1_BKstar_GKvD_LCSR_Lattice;
		a0T2 = param->a0T1_BKstar_GKvD_LCSR_Lattice;// Eq.17 in BSZ [1503.05534]
		a1T2 = param->a1T2_BKstar_GKvD_LCSR_Lattice;
		a2T2 = param->a2T2_BKstar_GKvD_LCSR_Lattice;
		a0T23 = param->a0T23_BKstar_GKvD_LCSR_Lattice;
		a1T23 = param->a1T23_BKstar_GKvD_LCSR_Lattice;
		a2T23 = param->a2T23_BKstar_GKvD_LCSR_Lattice;

		double P_V=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_A0=1./(1.-q2/pow(MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_A12=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T1=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T23=1./(1.-q2/pow(MAT_BKstar,2.));

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


	/********LCSR-only fit parameters from Gubernari van Dyk Kokulu 1811.00983***************/
	if(param->BKstar_FormFactor_choice == 4)
	{
		MB = mB;
		MV = mKstar;
		double MVT1_BKstar = 5.412;
		double MA0_BKstar = 5.336;
		double MAT_BKstar = 5.829;
		a0V = param->a0V_BKstar_GKvD_LCSRonly;
		a1V = param->a1V_BKstar_GKvD_LCSRonly;
		a2V = param->a2V_BKstar_GKvD_LCSRonly;
		a0A0 = param->a0A0_BKstar_GKvD_LCSRonly;
		a1A0 = param->a1A0_BKstar_GKvD_LCSRonly;
		a2A0 = param->a2A0_BKstar_GKvD_LCSRonly;
		a0A1 = param->a0A1_BKstar_GKvD_LCSRonly;
		a1A1 = param->a1A1_BKstar_GKvD_LCSRonly;
		a2A1 = param->a2A1_BKstar_GKvD_LCSRonly;
		a0A12 = (MB*MB - MV*MV)/(8.*MB*MV)*param->a0A0_BKstar_GKvD_LCSRonly;// Eq.17 in BSZ [1503.05534]
		a1A12 = param->a1A12_BKstar_GKvD_LCSRonly;
		a2A12 = param->a2A12_BKstar_GKvD_LCSRonly;
		a0T1 = param->a0T1_BKstar_GKvD_LCSRonly;
		a1T1 = param->a1T1_BKstar_GKvD_LCSRonly;
		a2T1 = param->a2T1_BKstar_GKvD_LCSRonly;
		a0T2 = param->a0T1_BKstar_GKvD_LCSRonly;// Eq.17 in BSZ [1503.05534]
		a1T2 = param->a1T2_BKstar_GKvD_LCSRonly;
		a2T2 = param->a2T2_BKstar_GKvD_LCSRonly;
		a0T23 = param->a0T23_BKstar_GKvD_LCSRonly;
		a1T23 = param->a1T23_BKstar_GKvD_LCSRonly;
		a2T23 = param->a2T23_BKstar_GKvD_LCSRonly;

		double P_V=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_A0=1./(1.-q2/pow(MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_A12=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T1=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T23=1./(1.-q2/pow(MAT_BKstar,2.));

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


	/********LCSR-only fit parameters from Zwicky et al. 1503.05534***************/
	if(param->BKstar_FormFactor_choice ==5)
	{
		MB = mB;
		MV = mKstar;
		a0V = param->a0V_BKstar_BSZ_LCSRonly;
		a1V = param->a1V_BKstar_BSZ_LCSRonly;
		a2V = param->a2V_BKstar_BSZ_LCSRonly;
		a0A0 = param->a0A0_BKstar_BSZ_LCSRonly;
		a1A0 = param->a1A0_BKstar_BSZ_LCSRonly;
		a2A0 = param->a2A0_BKstar_BSZ_LCSRonly;
		a0A1 = param->a0A1_BKstar_BSZ_LCSRonly;
		a1A1 = param->a1A1_BKstar_BSZ_LCSRonly;
		a2A1 = param->a2A1_BKstar_BSZ_LCSRonly;
		a0A12 = param->a0A12_BKstar_BSZ_LCSRonly;
		a1A12 = param->a1A12_BKstar_BSZ_LCSRonly;
		a2A12 = param->a2A12_BKstar_BSZ_LCSRonly;
		a0T1 = param->a0T1_BKstar_BSZ_LCSRonly;
		a1T1 = param->a1T1_BKstar_BSZ_LCSRonly;
		a2T1 = param->a2T1_BKstar_BSZ_LCSRonly;
		a0T2 = param->a0T2_BKstar_BSZ_LCSRonly;
		a1T2 = param->a1T2_BKstar_BSZ_LCSRonly;
		a2T2 = param->a2T2_BKstar_BSZ_LCSRonly;
		a0T23 = param->a0T23_BKstar_BSZ_LCSRonly;
		a1T23 = param->a1T23_BKstar_BSZ_LCSRonly;
		a2T23 = param->a2T23_BKstar_BSZ_LCSRonly;

		double P_V=1./(1.-q2/pow(param->MV_BKstar,2.));
		double P_A0=1./(1.-q2/pow(param->MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(param->MA1_BKstar,2.));
		double P_A12=1./(1.-q2/pow(param->MA12_BKstar,2.));
		double P_T1=1./(1.-q2/pow(param->MT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(param->MT2_BKstar,2.));
		double P_T23=1./(1.-q2/pow(param->MT23_BKstar,2.));

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
	if(param->BKstar_FormFactor_choice == 6)
	{
		MB = mB;
		MV = mKstar;
		double MVT1_BKstar = 5.414;//Mres=mB+deltam = 5.279 + 0.135
		double MA0_BKstar = 5.366;//Mres=mB+deltam = 5.279 + 0.087
		double MAT_BKstar = 5.829;//Mres=mB+deltam = 5.279 + 0.550
		a0V   = param->a0V_BKstar_HLMW;
		a1V   = param->a1V_BKstar_HLMW;
		a0A0  = param->a0A0_BKstar_HLMW;
		a1A0  = param->a1A0_BKstar_HLMW;
		a0A1  = param->a0A1_BKstar_HLMW;
		a1A1  = param->a1A1_BKstar_HLMW;
		a0A12 = param->a0A12_BKstar_HLMW;
		a1A12 = param->a1A12_BKstar_HLMW;
		a0T1  = param->a0T1_BKstar_HLMW;
		a1T1  = param->a1T1_BKstar_HLMW;
		a0T2  = param->a0T2_BKstar_HLMW;
		a1T2  = param->a1T2_BKstar_HLMW;
		a0T23 = param->a0T23_BKstar_HLMW;
		a1T23 = param->a1T23_BKstar_HLMW;
		double sysErr = param->systErr_BKstar_HLMW;

		double P_V=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_A0=1./(1.-q2/pow(MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_A12=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T1=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T23=1./(1.-q2/pow(MAT_BKstar,2.));

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


	double xi_par=(mB+mKstar)/2./E_Kstar*A1-(mB-mKstar)/mB*A2;
	double xi_perp=mB/(mB+mKstar)*V;


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

		double complex Cperpp0=C7eff+C7effp+q2/2./mb/mB*Y;
		double complex Cperpm0=C7eff-C7effp+q2/2./mb/mB*Y;
		//double complex Cparp0=-(C7eff+C7effp)-mB/2./mb*Y;
		double complex Cparm0=-(C7eff-C7effp)-mB/2./mb*Y;
		double complex Cperpp0_bar=conj(C7eff+C7effp)+q2/2./mb/mB*Y;
		double complex Cperpm0_bar=conj(C7eff-C7effp)+q2/2./mb/mB*Y;
		//double complex Cparp0_bar=-conj(C7eff+C7effp)-mB/2./mb*Y;
		double complex Cparm0_bar=-conj(C7eff-C7effp)-mB/2./mb*Y;

// 		double complex Cperpp0u=q2/2./mb/mB*Yu;
// 		double complex Cperpm0u=q2/2./mb/mB*Yu;
// 		//double complex Cparp0u=-mB/2./mb*Yu;
// 		double complex Cparm0u=-mB/2./mb*Yu;
	        /*******Used in redefinition of Tauperpu**********/
		double complex Cperp0u=q2/2./mb/mB*Yu;
		double complex Cpar0u=-mB/2./mb*Yu;

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
		-q2/2./mb/mB*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;
		
		double complex Cparnf=(C2bar*F27+C8eff*F87
		+mB/2./mb*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;

		double complex Cperpnfu=(-C2bar*(F27+F27_u)
		-q2/2./mb/mB*(C2bar*(F29+F29_u)+2.*C1bar*(F19+F19_u+1./6.*(F29+F29_u))))/4.*3.;
		
		double complex Cparnfu=(C2bar*(F27+F27_u)
		+mB/2./mb*(C2bar*(F29+F29_u)+2.*C1bar*(F19+F19_u+1./6.*(F29+F29_u))))/4.*3.;
		
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
		//double complex Cparp_bar=Cparp0_bar+alphas_mub*4./3./4./pi*Cparp1_bar;
		double complex Cparm_bar=Cparm0_bar+alphas_mub*4./3./4./pi*Cparm1_bar;

// 		double complex Cperppu=Cperpp0u+alphas_mub*4./3./4./pi*Cperpnfu; 
// 		double complex Cperpmu=Cperpm0u+alphas_mub*4./3./4./pi*Cperpnfu; 
// 		//double complex Cparpu=Cparp0u+alphas_mub*4./3./4./pi*Cparnfu; 
// 		double complex Cparmu=Cparm0u+alphas_mub*4./3./4./pi*Cparnfu; 
		/*******Used in redefinition of Tauperpu**********/
		double complex Cperpu=Cperp0u+alphas_mub*4./3./4./pi*Cperpnfu; 
		double complex Cparu=Cpar0u+alphas_mub*4./3./4./pi*Cparnfu; 
		
		
		double Xi_perp=1.;
		double Xi_par=mKstar/E_Kstar;
	
		double eu=2./3.;
		double ed=-1./3.;
	
		double a1perp=param->a1perp;
		double a2perp=param->a2perp;
		double a1par=param->a1par;
		double a2par=param->a2par;

		a1perp*=pow(eta,4./3.*(4.*1./2.)/(11.-2./3.*nf));
		a2perp*=pow(eta,4./3.*(4.*(1./2.+1./3.))/(11.-2./3.*nf));

		a1par*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
		a2par*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

		double u;
	
		double complex int_perppp,int_perppm,int_perpmp,int_perpmm;
		double complex /*int_parpp,*/int_parpm,int_parmp,int_parmm;// int_parpp commented; not used
		double complex int_perppu,int_parpu,int_parmu;
		double complex Tperppp0,Tperpppf,Tperpppnf,Tperppp;
		double complex Tperppm0,Tperppmf,Tperppmnf,Tperppm;
		double complex Tperpmp0,Tperpmpf,Tperpmpnf,Tperpmp;
		double complex Tperpmm0,Tperpmmf,Tperpmmnf,Tperpmm;
		double complex Tparpp0/*,Tparppf*/,Tparppnf/*,Tparpp*/;/*,Tparpp*/;// Tparppf and Tparpp commented; not used
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
		double complex /*Tparpp_bar,*/Tparmp_bar;// Tparpp_bar commented; not used
		double complex /*Tparppf_bar,*/Tparmpf_bar;// Tparppf_bar commented; not used
		
		int_perppp=int_perppm=int_perpmp=int_perpmm=0.;
		int_perppu=0.;
		/*int_parpp=*/int_parpm=int_parmp=int_parmm=0.;// int_parpp commented; not used
		int_parpu=int_parmu=0.;
		int_perppp_bar=int_perpmp_bar=0.;
		/*int_parpp_bar=*/int_parmp_bar=0.;// int_parpp_bar commented; not used


		double lambda_Bp=param->lambda_Bp;
		lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_f,2.))*(1.-2.*1.4);

		double omega0=2.*(mB-mb)/3.;
		double complex lambda_Bm=1./(exp(-q2/mB/omega0)/omega0*(-Ei(q2/mB/omega0)+I*pi));

		double phiKstar_perp,phiKstar_par;
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
		double deltatp=param->deltatp;
		double deltatm=param->deltatm;
		
		int n1=10;
		int n1sav=n1;
		for(ie=0;ie<=n1;ie++)
		{
			u=(double)ie/n1;
			if(ie==0) n1*=2;
			if(ie==n1){u=0.99;n1*=2;}

		/* Tperp */		
			Tperppp0=Tperpmp0=0.;
			Tperpppf=(C7eff+C7effp)*2.*mB/E_Kstar/(1.-u);
			Tperpmpf=(C7eff-C7effp)*2.*mB/E_Kstar/(1.-u);
			Tperpppf_bar=conj(C7eff+C7effp)*2.*mB/E_Kstar/(1.-u);
			Tperpmpf_bar=conj(C7eff-C7effp)*2.*mB/E_Kstar/(1.-u);


			Tperppm0=Tperpmm0=0.;
			Tperppmf=Tperpmmf=0.;
			
			phiKstar_perp=phi_Kstar(u,a1perp,a2perp);
			tperp_mc=tperp_bkll(u,mc,q2,E_Kstar,param);
			tperp_mb=tperp_bkll(u,mb,q2,E_Kstar,param);
			tperp_0=tperp_bkll(u,0.,q2,E_Kstar,param);

			Tperpppnf=Tperpmpnf=-4.*ed*C8eff/(u+(1.-u)*q2/mB/mB)
			+mB/2./mb*(eu*tperp_mc*(C2bar+C4bar-C6bar)
			+ed*tperp_mb*(C3bar+C4bar-C6bar-4.*mb/mB*C5bar)
			+ed*tperp_0*C3bar);
				
			Tperppmnf=Tperpmmnf=0.;

			Tperppnfu=mB/2./mb*eu*(tperp_mc-tperp_0)*(Cmub[2]-Cmub[1]/6.);
		
			Tperppp=Tperppp0+alphas_muf*4./3./4./pi*(Tperpppf+Tperpppnf); 
			Tperppm=Tperppm0+alphas_muf*4./3./4./pi*(Tperppmf+Tperppmnf);
			Tperpmp=Tperpmp0+alphas_muf*4./3./4./pi*(Tperpmpf+Tperpmpnf);
			Tperpmm=Tperpmm0+alphas_muf*4./3./4./pi*(Tperpmmf+Tperpmmnf);
			Tperppp_bar=Tperppp0+alphas_muf*4./3./4./pi*(Tperpppf_bar+Tperpppnf);
			Tperpmp_bar=Tperpmp0+alphas_muf*4./3./4./pi*(Tperpmpf_bar+Tperpmpnf);

			Tperppu=alphas_muf*4./3./4./pi*Tperppnfu;
	
			int_perppp+=phiKstar_perp*Tperppp/n1/lambda_Bp; 
			int_perppm+=phiKstar_perp*Tperppm/n1/lambda_Bm;
			int_perpmp+=phiKstar_perp*Tperpmp/n1/lambda_Bp;
			int_perpmm+=phiKstar_perp*Tperpmm/n1/lambda_Bm;
			int_perppp_bar+=phiKstar_perp*Tperppp_bar/n1/lambda_Bp;
			int_perpmp_bar+=phiKstar_perp*Tperpmp_bar/n1/lambda_Bp;

			int_perppu+=phiKstar_perp*Tperppu/n1/lambda_Bp;


		/* Tpar */		

			phiKstar_par=phi_Kstar(u,a1par,a2par);
			tpar_mc=tpar_bkll(u,mc,q2,E_Kstar,param);
			tpar_mb=tpar_bkll(u,mb,q2,E_Kstar,param);
			tpar_0=tpar_bkll(u,0.,q2,E_Kstar,param);
		
			Tparpp0=Tparmp0=0.;
		
			// Tparppf=(C7eff+C7effp)*4.*mB/E_Kstar/(1.-u);// commented; not used
			Tparmpf=(C7eff-C7effp)*4.*mB/E_Kstar/(1.-u);
			// Tparppf_bar=conj(C7eff+C7effp)*4.*mB/E_Kstar/(1.-u);// commented; not used
			Tparmpf_bar=conj(C7eff-C7effp)*4.*mB/E_Kstar/(1.-u);

			Tparppnf=Tparmpnf=mB/mb*(eu*tpar_mc*(C2bar+C4bar-C6bar)
			+ed*tpar_mb*(C3bar+C4bar-C6bar)
			+ed*tpar_0*C3bar);
				
			Tparpnfu=mB/mb*eu*(tpar_mc-tpar_0)*(Cmub[2]-Cmub[1]/6.);
	
			Tparpu=alphas_muf*4./3./4./pi*Tparpnfu;
	
	
			Tparpm0=Tparmm0=-eq*4.*mB/mb*(C3bar+3.*C4bar);
			
			Tparp0u=0.;

			if(charge==0) Tparm0u=0.; else Tparm0u=eq*12.*mB/mb*Cmub[2];
		
			Tparpmf=Tparmmf=0.;

			h_mc=h_bkll((1.-u)*mB*mB+u*q2,mc,mu_b);
			h_mb=h_bkll((1.-u)*mB*mB+u*q2,mbpole,mu_b);
			h_0=h_bkll((1.-u)*mB*mB+u*q2,0.,mu_b);

			Tparpmnf=Tparmmnf=eq*(8.*C8eff/((1.-u)+u*q2/mB/mB)
			+6.*mB/mb*(h_mc*(C2bar+C4bar+C6bar)
			+h_mb*(C3bar+C4bar+C6bar)
			+h_0*(C3bar+3.*C4bar+3.*C6bar)
			-8./27.*(C3bar-C5bar-15.*C6bar)));
	
			// Tparpp=Tparpp0+alphas_muf*4./3./4./pi*(Tparppf+Tparppnf);// commented; not used
			Tparpm=Tparpm0+alphas_muf*4./3./4./pi*(Tparpmf+Tparpmnf);
			Tparmp=Tparmp0+alphas_muf*4./3./4./pi*(Tparmpf+Tparmpnf);
			Tparmm=Tparmm0+alphas_muf*4./3./4./pi*(Tparmmf+Tparmmnf);
			// Tparpp_bar=Tparpp0+alphas_muf*4./3./4./pi*(Tparppf_bar+Tparppnf);// commented; not used
			Tparmp_bar=Tparmp0+alphas_muf*4./3./4./pi*(Tparmpf_bar+Tparmpnf);

			Tparmnfu=eq*(6.*mB/mb*(h_mc-h_0)*(Cmub[2]-Cmub[1]/6.));

			Tparpu=Tparp0u+alphas_muf*4./3./4./pi*Tparpnfu;
	
			Tparmu=Tparm0u+alphas_muf*4./3./4./pi*Tparmnfu;


			// int_parpp+=(phiKstar_par*Tparpp/lambda_Bp)/n1;// commented out; not used
			// int_parpp_bar+=(phiKstar_par*Tparpp_bar/lambda_Bp)/n1;// commented out; not used
			int_parpm+=(phiKstar_par*Tparpm/lambda_Bm)/n1;

			int_parmp+=(phiKstar_par*Tparmp/lambda_Bp)/n1;
			int_parmp_bar+=(phiKstar_par*Tparmp_bar/lambda_Bp)/n1;
			int_parmm+=(phiKstar_par*Tparmm/lambda_Bm)/n1;
			
			int_parpu+=(phiKstar_par*Tparpu/lambda_Bp)/n1;
			
			int_parmu+=(phiKstar_par*Tparmu/lambda_Bm)/n1;


			integ3+=phiKstar_perp/((1.-u)+u*q2/mB/mB)/n1;

			x=(1.-u)*mB*mB+u*q2;
			
			h_mc=h_bkll(x,mc,mu_b);
			h_mb=h_bkll(x,mbpole,mu_b);
			h_0=h_bkll(x,0.,mu_b);
			
			FV=3./4.*(h_mc*(C2bar+C4bar+C6bar)+h_mb*(C3bar+C4bar+C6bar)+h_0*(C3bar+3.*C4bar+3.*C6bar)-8./27.*(C3bar-C5bar-15.*C6bar));

			FVu=3./4.*(h_mc-h_0)*(Cmub[2]-Cmub[1]/6.);
		

			integ4+=phiKstar_perp/((1.-u)+u*q2/mB/mB)*FV/n1;
			integ4u+=phiKstar_perp/((1.-u)+u*q2/mB/mB)*FVu/n1;
		
			Fperp+=phiKstar_perp/((1.-u)+u*q2/mB/mB)/3./n1;
			Xperp+=(u<=1.-0.5/mB)*phiKstar_perp/pow((1.-u)+u*q2/mB/mB,2.)/3./n1;
		
			integ5+=((3./4.*(1.+pow(2.*u-1.,2.))+a1par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./112.*a2par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1par*(2.*u-1.)+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1par*u+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FV/n1;

			integ5u+=((3./4.*(1.+pow(2.*u-1.,2.))+a1par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./112.*a2par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1par*(2.*u-1.)+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1par*u+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FVu/n1;

			if(ie==0||ie==n1sav) n1=n1sav;
		}
		
		/* Tau_perp */		

		double complex Tauperpp=xi_perp*Cperpp+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*(int_perppp+int_perppm); 
		double complex Tauperpp_bar=xi_perp*Cperpp_bar+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*(int_perppp_bar+int_perppm);

		double complex Tauperpm=xi_perp*Cperpm+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*(int_perpmp+int_perpmm); 
		double complex Tauperpm_bar=xi_perp*Cperpm_bar+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*(int_perpmp_bar+int_perpmm);

		double complex Tauperpu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_perp*Cperpu+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*int_perppu);
		double complex Tauperpu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_perp*Cperpu+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*int_perppu);


		/* Tau_par */		


		double complex Tauparm=xi_par*Cparm+pi*pi/3.*param->f_B*f_Kstar_par/mB*Xi_par*(int_parmp+int_parmm);
		double complex Tauparm_bar=xi_par*Cparm_bar+pi*pi/3.*param->f_B*f_Kstar_par/mB*Xi_par*(int_parmp_bar+int_parmm);

		double complex Tauparu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_par*Cparu+pi*pi/3.*param->f_B*f_Kstar_par/mB*Xi_par*(int_parpu+int_parmu));
		double complex Tauparu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_par*Cparu+pi*pi/3.*param->f_B*f_Kstar_par/mB*Xi_par*(int_parpu+int_parmu));
				
				
		double complex DeltaTauperpWA=-eq*4.*pi*pi/3.*param->f_B*f_Kstar_perp/mb/mB*(Cmub[3]+4./3.*Cmub[4]+4.*Cmub[5]+16./3.*Cmub[6])*integ3
		+eq*2.*pi*pi/3.*param->f_B*f_Kstar_par/mb/mB*mKstar/(1.-q2/mB/mB)/lambda_Bp*(Cmub[3]+4./3.*Cmub[4]+16.*Cmub[5]+64./3.*Cmub[6]);

		double complex DeltaTauperpuWA=0.;
		if(charge!=0) DeltaTauperpuWA=-eq*2.*pi*pi*param->f_B*f_Kstar_par/mb/mB*mKstar/(1.-q2/mB/mB)/lambda_Bp*Cmub[2];


		double rho=0.;
		double phi=0.;
		Xperp=Fperp+(1.+rho*(cos(phi)+I*sin(phi)))*Xperp;


		double complex DeltaTauperpHSA=eq*4./3.*alphas_muf/4./pi*pi*pi*param->f_B/3./mb/mB*(12.*C8eff*mb/mB*f_Kstar_perp*Xperp
		+8.*f_Kstar_perp*integ4-4.*mKstar*f_Kstar_par/(1.-q2/mB/mB)/lambda_Bp*integ5);
		
		double complex DeltaTauperpuHSA=eq*4./3.*alphas_muf/4./pi*pi*pi*param->f_B/3./mb/mB*
		(8.*f_Kstar_perp*integ4u-4.*mKstar*f_Kstar_par/(1.-q2/mB/mB)/lambda_Bp*integ5u);

		Tauperpp+=DeltaTauperpWA+DeltaTauperpHSA;
		Tauperpm+=DeltaTauperpWA+DeltaTauperpHSA;
		Tauperpp_bar+=DeltaTauperpWA+DeltaTauperpHSA;
		Tauperpm_bar+=DeltaTauperpWA+DeltaTauperpHSA;

// 		Tauperppu+=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(DeltaTauperpuWA+DeltaTauperpuHSA);
// 		Tauperppu_bar+=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(DeltaTauperpuWA+DeltaTauperpuHSA);	
// 		
// 		Tauperpmu+=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(DeltaTauperpuWA+DeltaTauperpuHSA);
// 		Tauperpmu_bar+=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(DeltaTauperpuWA+DeltaTauperpuHSA);		
	        /*******Using Tauperpu**********/
	        Tauperpu+=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(DeltaTauperpuWA+DeltaTauperpuHSA);
	        Tauperpu_bar+=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(DeltaTauperpuWA+DeltaTauperpuHSA);    

		double lambda=pow(mB,4.)+pow(mKstar,4.)+q2*q2-2.*(mB*mB*mKstar*mKstar+mKstar*mKstar*q2+mB*mB*q2);
	
		double complex N=param->Vtb*conj(param->Vts)*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)/pow(mB,3.)*q2*sqrt(lambda)*beta_l);
		
		double complex Nbar=conj(param->Vtb)*param->Vts*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)/pow(mB,3.)*q2*sqrt(lambda)*beta_l);
		

		int n9=10;
		double integ9=0.;
		for(ie=1;ie<=n9-1;ie++)
		{
			u=(double)ie/n9;	
			integ9+=phi_Kstar(u,a1par,a2par)/(1.-u)/n9;
		}
		double Delta_par=1.+alphas_mub*4./3./4./pi*(-2.+2.*L)-alphas_mub*4./3./4./pi*2.*q2/E_Kstar/E_Kstar*pi*pi*param->f_B*f_Kstar_par/lambda_Bp/3./param->m_Bd/(E_Kstar/param->m_Kstar)/xi_par*integ9;

		 	
		ALperp=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)-(C10+C10p))*V/(mB+mKstar)+2.*mb/q2*(Tauperpp+Tauperpu));
		ARperp=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)+(C10+C10p))*V/(mB+mKstar)+2.*mb/q2*(Tauperpp+Tauperpu));
		
		ALpar=-N*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((C9-C9p)-(C10-C10p))*((2.*E_Kstar/(mB+mKstar)*xi_perp))/(mB-mKstar)+4.*mb/mB*E_Kstar/q2*(Tauperpm+Tauperpu));
		ARpar=-N*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((C9-C9p)+(C10-C10p))*((2.*E_Kstar/(mB+mKstar)*xi_perp))/(mB-mKstar)+4.*mb/mB*E_Kstar/q2*(Tauperpm+Tauperpu));
	
		AL0=-N/2./mKstar/sqrt(q2)*(((C9-C9p)-(C10-C10p))*((mB*mB-mKstar*mKstar-q2)*(mB+mKstar)*((2.*E_Kstar/(mB+mKstar)*xi_perp))-lambda*((mB/(mB-mKstar)*(xi_perp-xi_par)))/(mB+mKstar))+2.*mb*(2.*E_Kstar/mB*(mB*mB+3.*mKstar*mKstar-q2)*(Tauperpm+Tauperpu)-lambda/(mB*mB-mKstar*mKstar)*(Tauperpm+Tauparm+Tauperpu+Tauparu)));
		AR0=-N/2./mKstar/sqrt(q2)*(((C9-C9p)+(C10-C10p))*((mB*mB-mKstar*mKstar-q2)*(mB+mKstar)*((2.*E_Kstar/(mB+mKstar)*xi_perp))-lambda*((mB/(mB-mKstar)*(xi_perp-xi_par)))/(mB+mKstar))+2.*mb*(2.*E_Kstar/mB*(mB*mB+3.*mKstar*mKstar-q2)*(Tauperpm+Tauperpu)-lambda/(mB*mB-mKstar*mKstar)*(Tauperpm+Tauparm+Tauperpu+Tauparu)));

		At=N/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/ml*(CQ2-CQ2p)/(mb_mub+param->mass_s))*(E_Kstar/mKstar*xi_par/Delta_par);
	
		AS=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(mb_mub+param->mass_s)*(E_Kstar/mKstar*xi_par/Delta_par);
			
				
		ALperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(conj((C9+C9p)-(C10+C10p))*V/(mB+mKstar)+2.*mb/q2*(Tauperpp_bar+Tauperpu_bar));
		ARperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(conj((C9+C9p)+(C10+C10p))*V/(mB+mKstar)+2.*mb/q2*(Tauperpp_bar+Tauperpu_bar));
	
		ALpar_bar=-Nbar*sqrt(2.)*(mB*mB-mKstar*mKstar)*(conj((C9-C9p)-(C10-C10p))*((2.*E_Kstar/(mB+mKstar)*xi_perp))/(mB-mKstar)+4.*mb/mB*E_Kstar/q2*(Tauperpm_bar+Tauperpu_bar));
		ARpar_bar=-Nbar*sqrt(2.)*(mB*mB-mKstar*mKstar)*(conj((C9-C9p)+(C10-C10p))*((2.*E_Kstar/(mB+mKstar)*xi_perp))/(mB-mKstar)+4.*mb/mB*E_Kstar/q2*(Tauperpm_bar+Tauperpu_bar));
	
		AL0_bar=-Nbar/2./mKstar/sqrt(q2)*(conj((C9-C9p)-(C10-C10p))*((mB*mB-mKstar*mKstar-q2)*(mB+mKstar)*((2.*E_Kstar/(mB+mKstar)*xi_perp))-lambda*((mB/(mB-mKstar)*(xi_perp-xi_par)))/(mB+mKstar))+2.*mb*(2.*E_Kstar/mB*(mB*mB+3.*mKstar*mKstar-q2)*(Tauperpm_bar+Tauperpu_bar)-lambda/(mB*mB-mKstar*mKstar)*(Tauperpm_bar+Tauparm+Tauperpu_bar+Tauparu_bar)));
		AR0_bar=-Nbar/2./mKstar/sqrt(q2)*(conj((C9-C9p)+(C10-C10p))*((mB*mB-mKstar*mKstar-q2)*(mB+mKstar)*((2.*E_Kstar/(mB+mKstar)*xi_perp))-lambda*((mB/(mB-mKstar)*(xi_perp-xi_par)))/(mB+mKstar))+2.*mb*(2.*E_Kstar/mB*(mB*mB+3.*mKstar*mKstar-q2)*(Tauperpm_bar+Tauperpu_bar)-lambda/(mB*mB-mKstar*mKstar)*(Tauperpm_bar+Tauparm_bar+Tauperpu_bar+Tauparu_bar)));

		At_bar=Nbar/sqrt(q2)*sqrt(lambda)*(2.*conj(C10-C10p)+q2/ml*conj(CQ2-CQ2p)/(mb_mub+param->mass_s))*(E_Kstar/mKstar*xi_par/Delta_par);
	
		AS_bar=-2.*Nbar*sqrt(lambda)*conj(CQ1-CQ1p)/(mb_mub+param->mass_s)*(E_Kstar/mKstar*xi_par/Delta_par);
			
		/* hadronic uncertainties */
		ALperp*=1.+param->BtoKstarlow_ALperp_err_noq2+q2/6.*param->BtoKstarlow_ALperp_err_q2;
		ARperp*=1.+param->BtoKstarlow_ARperp_err_noq2+q2/6.*param->BtoKstarlow_ARperp_err_q2;
		ALpar*=1.+param->BtoKstarlow_ALpar_err_noq2+q2/6.*param->BtoKstarlow_ALpar_err_q2;
		ARpar*=1.+param->BtoKstarlow_ARpar_err_noq2+q2/6.*param->BtoKstarlow_ARpar_err_q2;
		AL0*=1.+param->BtoKstarlow_AL0_err_noq2+q2/6.*param->BtoKstarlow_AL0_err_q2;
		AR0*=1.+param->BtoKstarlow_AR0_err_noq2+q2/6.*param->BtoKstarlow_AR0_err_q2;

		ALperp_bar*=1.+param->BtoKstarlow_ALperp_err_noq2+q2/6.*param->BtoKstarlow_ALperp_err_q2;
		ARperp_bar*=1.+param->BtoKstarlow_ARperp_err_noq2+q2/6.*param->BtoKstarlow_ARperp_err_q2;
		ALpar_bar*=1.+param->BtoKstarlow_ALpar_err_noq2+q2/6.*param->BtoKstarlow_ALpar_err_q2;
		ARpar_bar*=1.+param->BtoKstarlow_ARpar_err_noq2+q2/6.*param->BtoKstarlow_ARpar_err_q2;
		AL0_bar*=1.+param->BtoKstarlow_AL0_err_noq2+q2/6.*param->BtoKstarlow_AL0_err_q2;
		AR0_bar*=1.+param->BtoKstarlow_AR0_err_noq2+q2/6.*param->BtoKstarlow_AR0_err_q2;
	}
	
	if(q2>LOW_q2)
	{
		// double mb=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);// => mb=4.08
		double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
		double mb_PS=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */ //mb=>4.64
		double sh=q2/mb_PS/mb_PS;

		double shat=q2/mB/mB;

		double lambda_hat=1.+shat*shat+pow(mKstar/mB,4.)-2.*(shat+shat*mKstar*mKstar/mB/mB+mKstar*mKstar/mB/mB);

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

		double complex N=param->Vtb*conj(param->Vts)*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em*beta_l/3./1024./pow(pi,5.)*mB*shat*sqrt(lambda_hat));// Missing beta_l corrected

		double complex Nbar=conj(param->Vtb)*param->Vts*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em*beta_l/3./1024./pow(pi,5.)*mB*shat*sqrt(lambda_hat));// Missing beta_l corrected

		double complex f_perp=N*mB*sqrt(2.*lambda_hat)/(1.+mKstar/mB)*V;

		double complex f_par=N*mB*sqrt(2.)*(1.+mKstar/mB)*A1;

		double complex f_0=N*mB*((1.-shat-mKstar*mKstar/mB/mB)*pow(1.+mKstar/mB,2.)*A1-lambda_hat*A2)/(2.*mKstar/mB*(1.+mKstar/mB)*sqrt(shat));

		double complex f_perp_bar=Nbar*mB*sqrt(2.*lambda_hat)/(1.+mKstar/mB)*V;

		double complex f_par_bar=Nbar*mB*sqrt(2.)*(1.+mKstar/mB)*A1;

		double complex f_0_bar=Nbar*mB*((1.-shat-mKstar*mKstar/mB/mB)*pow(1.+mKstar/mB,2.)*A1-lambda_hat*A2)/(2.*mKstar/mB*(1.+mKstar/mB)*sqrt(shat));

		ALperp_high=(((C9eff+C9p)-(Cmub[10]+C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff+C7effp))*f_perp;
		ARperp_high=(((C9eff+C9p)+(Cmub[10]+C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff+C7effp))*f_perp;

		ALpar_high=-(((C9eff-C9p)-(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff-C7effp))*f_par;
		ARpar_high=-(((C9eff-C9p)+(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff-C7effp))*f_par;

		AL0_high=-(((C9eff-C9p)-(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff-C7effp))*f_0;
		AR0_high=-(((C9eff-C9p)+(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff-C7effp))*f_0;


		ALperp_bar_high=(((C9eff_bar+conj(C9p))-conj(Cmub[10]+C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar+conj(C7effp)))*f_perp_bar;
		ARperp_bar_high=(((C9eff_bar+conj(C9p))+conj(Cmub[10]+C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar+conj(C7effp)))*f_perp_bar;

		ALpar_bar_high=-(((C9eff_bar-conj(C9p))-conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar-conj(C7effp)))*f_par_bar;
		ARpar_bar_high=-(((C9eff_bar-conj(C9p))+conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar-conj(C7effp)))*f_par_bar;

		AL0_bar_high=-(((C9eff_bar-conj(C9p))-conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar-conj(C7effp)))*f_0_bar;
		AR0_bar_high=-(((C9eff_bar-conj(C9p))+conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar-conj(C7effp)))*f_0_bar;

		double complex C10=Cmub[10];
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];

		double lambda=pow(mB,4.)+pow(mKstar,4.)+q2*q2-2.*(mB*mB*mKstar*mKstar+mKstar*mKstar*q2+mB*mB*q2);

		At_high=N/sqrt(q2)*sqrt(lambda)*( 2.*(C10-C10p)+q2/ml*(CQ2-CQ2p)/(param->mass_b+param->mass_s))*A0;
		AS_high=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(param->mass_b+param->mass_s)*A0;

		At_bar_high=Nbar/sqrt(q2)*sqrt(lambda)*( 2.*conj(C10-C10p)+q2/ml*conj(CQ2-CQ2p)/(param->mass_b+param->mass_s))*A0;
		AS_bar_high=-2.*Nbar*sqrt(lambda)*conj(CQ1-CQ1p)/(param->mass_b+param->mass_s)*A0;


		/* hadronic uncertainties */
		ALperp_high*=1.+param->BtoKstarhigh_ALperp_err;
		ARperp_high*=1.+param->BtoKstarhigh_ARperp_err;
		ALpar_high*=1.+param->BtoKstarhigh_ALpar_err;
		ARpar_high*=1.+param->BtoKstarhigh_ARpar_err;
		AL0_high*=1.+param->BtoKstarhigh_AL0_err;
		AR0_high*=1.+param->BtoKstarhigh_AR0_err;
		At_high*=1.+param->BtoKstarhigh_At_err;
		AS_high*=1.+param->BtoKstarhigh_AS_err;
		ALperp_bar_high*=1.+param->BtoKstarhigh_ALperp_err;
		ARperp_bar_high*=1.+param->BtoKstarhigh_ARperp_err;
		ALpar_bar_high*=1.+param->BtoKstarhigh_ALpar_err;
		ARpar_bar_high*=1.+param->BtoKstarhigh_ARpar_err;
		AL0_bar_high*=1.+param->BtoKstarhigh_AL0_err;
		AR0_bar_high*=1.+param->BtoKstarhigh_AR0_err;
		At_bar_high*=1.+param->BtoKstarhigh_At_err;
		AS_bar_high*=1.+param->BtoKstarhigh_AS_err;

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
	
	double dGamma_BKstarll_dq2=3./4.*(2.*J1s+J1c-(2.*J2s+J2c)/3.);

	double dGamma_bar_BKstarll_dq2=3./4.*(2.*J1s_bar+J1c_bar-(2.*J2s_bar+J2c_bar)/3.);

	double AFB[3],FL[3],FT[3],AT1[3],AT2[3],AT3[3],AT4[3],AT5[3],HT1[3],HT2[3],HT3[3],alpha_Kstar[3],AIm[3],P2[3],P3[3],P6[3],P8[3],P4prime[3],P5prime[3],P6prime[3],P8prime[3],A3[3],A4[3],A5[3],A7[3],A8[3],A9[3],S3[3],S4[3],S5[3],S7[3],S8[3],S9[3];


	int scalar_switch = 1;// 0 if we want to put scalar to zero in AFB, else 1
	AFB[0]=-3./8.*(2.*(J6s+J6s_bar)+scalar_switch*(J6c+J6c_bar))/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	AFB[1]=-3./8.*(2.*(J6s+J6s_bar)+scalar_switch*(J6c+J6c_bar));
	AFB[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	FL[0]=1./4.*(3.*(J1c+J1c_bar)-(J2c+J2c_bar))/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	FL[1]=1./4.*(3.*(J1c+J1c_bar)-(J2c+J2c_bar));
	FL[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	FT[0]=4.*(J2s+J2s_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	FT[1]=4.*(J2s+J2s_bar);
	FT[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	AT1[0]=-2.*creal(ALpar*conj(ALperp)+ARpar*conj(ARperp)+ALpar_bar*conj(ALperp_bar)+ARpar_bar*conj(ARperp_bar))/(Apar2+Aperp2+Apar2_bar+Aperp2_bar);
	AT1[1]=-2.*creal(ALpar*conj(ALperp)+ARpar*conj(ARperp)+ALpar_bar*conj(ALperp_bar)+ARpar_bar*conj(ARperp_bar));
	AT1[2]=Apar2+Aperp2+Apar2_bar+Aperp2_bar;

	AT2[0]=(J3+J3_bar)/2./(J2s+J2s_bar);
	AT2[1]=J3+J3_bar;
	AT2[2]=2.*(J2s+J2s_bar);

	AT3[0]=sqrt(fabs((4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar))/(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar))));
	AT3[1]=sqrt(fabs(4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar)));
	AT3[2]=sqrt(fabs(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	
	AT4[0]=sqrt((beta_l*beta_l*(J5+J5_bar)*(J5+J5_bar)+4.*(J8+J8_bar)*(J8+J8_bar))/(4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar)));
	AT4[1]=sqrt(beta_l*beta_l*(J5+J5_bar)*(J5+J5_bar)+4.*(J8+J8_bar)*(J8+J8_bar));
	AT4[2]=sqrt(4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar));
	
	AT5[0]=cabs(ALperp*conj(ARpar)+ALpar*conj(ARperp)+ALperp_bar*conj(ARpar_bar)+ALpar_bar*conj(ARperp_bar))/(Apar2+Aperp2+Apar2_bar+Aperp2_bar);
	AT5[1]=cabs(ALperp*conj(ARpar)+ALpar*conj(ARperp)+ALperp_bar*conj(ARpar_bar)+ALpar_bar*conj(ARperp_bar));
	AT5[2]=Apar2+Aperp2+Apar2_bar+Aperp2_bar;
	
	HT1[0]=sqrt(2.)*(J4+J4_bar)/sqrt(fabs(-(J2c+J2c_bar)*(2.*J2s-J3+2.*J2s_bar-J3_bar)));
	HT1[1]=sqrt(2.)*(J4+J4_bar);
	HT1[2]=sqrt(fabs(-(J2c+J2c_bar)*(2.*J2s-J3+2.*J2s_bar-J3_bar)));
	
	HT2[0]=beta_l*(J5+J5_bar)/sqrt(fabs(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	HT2[1]=beta_l*(J5+J5_bar);
	HT2[2]=sqrt(fabs(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	
	HT3[0]=(J6s+J6s_bar)/2./sqrt(fabs(4.*(J2s+J2s_bar)*(J2s+J2s_bar)-(J3+J3_bar)*(J3+J3_bar)));
	HT3[1]=(J6s+J6s_bar)/2.;
	HT3[2]=J3+J3_bar;

	alpha_Kstar[0]=-(2.*J2s+J2c+2.*J2s_bar+J2c_bar)/2./(J2s+J2s_bar);
	alpha_Kstar[1]=-(2.*J2s+J2c+2.*J2s_bar+J2c_bar);
	alpha_Kstar[2]=2.*(J2s+J2s_bar);

	AIm[0]=(J9+J9_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	AIm[1]=J9+J9_bar;
	AIm[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	P2[0]=(J6s+J6s_bar)/8./(J2s+J2s_bar);// Definition is from 1303.5794 without beta_l and not from 1202.4266
	P2[1]=J6s+J6s_bar;
	P2[2]=8.*(J2s+J2s_bar);

	P3[0]=-(J9+J9_bar)/4./(J2s+J2s_bar);
	P3[1]=-(J9+J9_bar);
	P3[2]=4.*(J2s+J2s_bar);
	
	P6[0]=-beta_l*(J7+J7_bar)/sqrt(fabs(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	P6[1]=-beta_l*(J7+J7_bar);
	P6[2]=sqrt(fabs(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	
	P4prime[0]=(J4+J4_bar)/sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P4prime[1]=J4+J4_bar;
	P4prime[2]=-(J2c+J2c_bar);
	
	P5prime[0]=(J5+J5_bar)/2./sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P5prime[1]=(J5+J5_bar)/2.;
	P5prime[2]=-(J2c+J2c_bar);
	
	P6prime[0]=-(J7+J7_bar)/2./sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P6prime[1]=-(J7+J7_bar)/2.;
	P6prime[2]=-(J2c+J2c_bar);

	P8[0]=-sqrt(2.)*(J8+J8_bar)/sqrt(fabs(-(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	P8[1]=-sqrt(2.)*(J8+J8_bar);
	P8[2]=2.*J2s+J3+2.*J2s_bar+J3_bar;
	
	P8prime[0]=-(J8+J8_bar)/sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P8prime[1]=-(J8+J8_bar);
	P8prime[2]=-(J2c+J2c_bar);
	
	A3[0]=(J3-J3_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A3[1]=J3-J3_bar;
	A3[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	A4[0]=(J4-J4_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A4[1]=J4-J4_bar;
	A4[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	
	A5[0]=(J5-J5_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A5[1]=J5-J5_bar;
	A5[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	A7[0]=(J7-J7_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A7[1]=J7-J7_bar;
	A7[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	
	A8[0]=(J8-J8_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A8[1]=J8-J8_bar;
	A8[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	
	A9[0]=(J9-J9_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A9[1]=J9-J9_bar;
	A9[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	
	S3[0]=(J3+J3_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S3[1]=J3+J3_bar;
	S3[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S4[0]=(J4+J4_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S4[1]=J4+J4_bar;
	S4[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	
	S5[0]=(J5+J5_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S5[1]=J5+J5_bar;
	S5[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S7[0]=(J7+J7_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S7[1]=J7+J7_bar;
	S7[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S8[0]=(J8+J8_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S8[1]=J8+J8_bar;
	S8[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S9[0]=(J9+J9_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S9[1]=J9+J9_bar;
	S9[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	double ATRe[3], ATReCP[3], P1CP[3], P2CP[3], P3CP[3], P4primeCP[3], P5primeCP[3], P6primeCP[3], P8primeCP[3], S6c[3], A6s[3], ACP[3];

	ATRe[0]=beta_l*(J6s+J6s_bar)/(4.*(J2s+J2s_bar));//based on 1106.3283 which includes the beta_l function; ATRe = beta_l*2*P2
	ATRe[1]=beta_l*(J6s+J6s_bar);
	ATRe[2]=4.*(J2s+J2s_bar);

	ATReCP[0]=beta_l*(J6s-J6s_bar)/(4.*(J2s+J2s_bar));//based on 1106.3283 which includes the beta_l function; ATRe = beta_l*2*P2CP
	ATReCP[1]=beta_l*(J6s-J6s_bar);
	ATReCP[2]=4.*(J2s+J2s_bar);

	P1CP[0]=(J3-J3_bar)/(2.*(J2s+J2s_bar));
	P1CP[1]=(J3-J3_bar);
	P1CP[2]=2.*(J2s+J2s_bar);

	P2CP[0]=(J6s-J6s_bar)/8./(J2s+J2s_bar);// Definition from 1303.5794 without beta_l and not from 1202.4266
	P2CP[1]=J6s-J6s_bar;
	P2CP[2]=8.*(J2s+J2s_bar);

	P3CP[0]=-(J9-J9_bar)/4./(J2s+J2s_bar);
	P3CP[1]=-(J9-J9_bar);
	P3CP[2]=4.*(J2s+J2s_bar);

	P4primeCP[0]=(J4-J4_bar)/sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P4primeCP[1]=J4-J4_bar;
	P4primeCP[2]=-(J2c+J2c_bar);

	P5primeCP[0]=(J5-J5_bar)/2./sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P5primeCP[1]=(J5-J5_bar)/2.;
	P5primeCP[2]=-(J2c+J2c_bar);

	P6primeCP[0]=-(J7-J7_bar)/2./sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P6primeCP[1]=-(J7-J7_bar)/2.;
	P6primeCP[2]=-(J2c+J2c_bar);

	P8primeCP[0]=-(J8-J8_bar)/sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P8primeCP[1]=-(J8-J8_bar);
	P8primeCP[2]=-(J2c+J2c_bar);

	/* We use the same convention as AFB above regarding the sign */
	S6c[0]=-(J6c+J6c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S6c[1]=-(J6c+J6c_bar);
	S6c[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	/* We use the same convention as AFB above regarding the sign (no need for change of sign compared to LHCb) */
	A6s[0]=-(J6s-J6s_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A6s[1]=-(J6s-J6s_bar);
	A6s[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	ACP[0]=(dGamma_BKstarll_dq2-dGamma_bar_BKstarll_dq2)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	ACP[1]=(dGamma_BKstarll_dq2-dGamma_bar_BKstarll_dq2);
	ACP[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;


	double S1c[3], S2c[3], S2s[3], AFL[3], A1c[3], A2s[3], A6c[3], AFBCP[3];

	S1c[0]=(J1c+J1c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S1c[1]=(J1c+J1c_bar);
	S1c[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S2c[0]=(J2c+J2c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S2c[1]=(J2c+J2c_bar);
	S2c[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S2s[0]=(J2s+J2s_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S2s[1]=(J2s+J2s_bar);
	S2s[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	AFL[0]=-(J2c-J2c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	AFL[1]=-(J2c-J2c_bar);
	AFL[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	A1c[0]=(J1c-J1c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A1c[1]=(J1c-J1c_bar);
	A1c[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	A2s[0]=(J2s-J2s_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A2s[1]=(J2s-J2s_bar);
	A2s[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	/* See explanation above AFB for `scalar_switch` */
	AFBCP[0]=-3./8.*(2.*(J6s-J6s_bar)+scalar_switch*(J6c-J6c_bar))/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	AFBCP[1]=-3./8.*(2.*(J6s-J6s_bar)+scalar_switch*(J6c-J6c_bar));
	AFBCP[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	A6c[0]=-(J6c-J6c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A6c[1]=-(J6c-J6c_bar);
	A6c[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	/*********************************************/

	for(je=0;je<=Nobs_BKsll;je++) for(ie=0;ie<=2;ie++) obs[je][ie]=0.;
	
	for(ie=0;ie<=2;ie++)
	{
		obs[1][ie]=AFB[ie];
		obs[2][ie]=FL[ie];
		obs[3][ie]=FT[ie];
		obs[4][ie]=AT1[ie];
		obs[5][ie]=AT2[ie];
		obs[6][ie]=AT3[ie];
		obs[7][ie]=AT4[ie];
		obs[8][ie]=AT5[ie];
		obs[9][ie]=HT1[ie];
		obs[10][ie]=HT2[ie];
		obs[11][ie]=HT3[ie];
		obs[12][ie]=alpha_Kstar[ie];
		obs[13][ie]=AIm[ie];
		obs[14][ie]=P2[ie];
		obs[15][ie]=P3[ie];
		obs[16][ie]=P6[ie];
		obs[17][ie]=P4prime[ie];
		obs[18][ie]=P5prime[ie];
		obs[19][ie]=P6prime[ie];
		obs[20][ie]=P8[ie];
		obs[21][ie]=P8prime[ie];
		obs[22][ie]=A7[ie];
		obs[23][ie]=A8[ie];
		obs[24][ie]=A9[ie];
		obs[25][ie]=S3[ie];
		obs[26][ie]=S4[ie];
		obs[27][ie]=S5[ie];
		obs[28][ie]=S7[ie];
		obs[29][ie]=S8[ie];
		obs[30][ie]=S9[ie];
		obs[31][ie]=A3[ie];
		obs[32][ie]=A4[ie];
		obs[33][ie]=A5[ie];
		obs[34][ie]=ATRe[ie];
		obs[35][ie]=ATReCP[ie];
		obs[36][ie]=P1CP[ie];
		obs[37][ie]=P2CP[ie];
		obs[38][ie]=P3CP[ie];
		obs[39][ie]=P4primeCP[ie];
		obs[40][ie]=P5primeCP[ie];
		obs[41][ie]=P6primeCP[ie];
		obs[42][ie]=P8primeCP[ie];
		obs[43][ie]=S6c[ie];
		obs[44][ie]=A6s[ie];
		obs[45][ie]=ACP[ie];
		obs[46][ie]=S1c[ie];
		obs[47][ie]=S2c[ie];
		obs[48][ie]=S2s[ie];
		obs[49][ie]=AFL[ie];
		obs[50][ie]=A1c[ie];
		obs[51][ie]=A2s[ie];
		obs[52][ie]=AFBCP[ie];
		obs[53][ie]=A6c[ie];
	}

	return (dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2)/2.;
}

/*----------------------------------------------------------------------*/

double dGamma_BKstarll_dq2_soft_calculator(int gen, int charge, double q2, double obs[][3], char name[])
/* "container" function scanning the SLHA file "name" and calculating dGamma/dq2(B->Kstar mu+ mu-) */
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

	return dGamma_BKstarll_dq2_soft(gen,charge,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarll_soft(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	int ie,je;
	int nmax=10;
	if((smin<0.099)||(smax-smin>10.)) nmax=100;
	double Gamma=0.;
	double s;
	
	double obs_num[Nobs_BKsll+1],obs_den[Nobs_BKsll+1];
	for(je=0;je<=Nobs_BKsll;je++) obs_num[je]=obs_den[je]=0.;

	obs[0]=0.; /* zero AFB */
	obs[1]=0.; /* integrated AFB */
	obs[2]=0.; /* integrated FL */
	obs[3]=0.; /* integrated FT */
	obs[4]=0.; /* integrated AT1 */
	obs[5]=0.; /* integrated AT2 */
	obs[6]=0.; /* integrated AT3 */
	obs[7]=0.; /* integrated AT4 */
	obs[8]=0.; /* integrated AT5 */
	obs[9]=0.; /* integrated HT1 */
	obs[10]=0.; /* integrated HT2 */
	obs[11]=0.; /* integrated HT3 */
	obs[12]=0.; /* integrated alpha */
	obs[13]=0.; /* integrated AIm */
	obs[14]=0.; /* integrated P2 */
	obs[15]=0.; /* integrated P3 */
	obs[16]=0.; /* integrated P6 */
	obs[17]=0.; /* integrated P4prime */
	obs[18]=0.; /* integrated P5prime */
	obs[19]=0.; /* integrated P6prime */
	obs[20]=0.; /* integrated P8 */
	obs[21]=0.; /* integrated P8prime */
	obs[22]=0.; /* integrated A7 */
	obs[23]=0.; /* integrated A8 */
	obs[24]=0.; /* integrated A9 */
	obs[25]=0.; /* integrated S3 */
	obs[26]=0.; /* integrated S4 */
	obs[27]=0.; /* integrated S5 */
	obs[28]=0.; /* integrated S7 */
	obs[29]=0.; /* integrated S8 */
	obs[30]=0.; /* integrated S9 */
	obs[31]=0.; /* integrated A3 */
	obs[32]=0.; /* integrated A4 */
	obs[33]=0.; /* integrated A5 */
	obs[34]=0.; /* integrated ATRe */
	obs[35]=0.; /* integrated ATReCP */
	obs[36]=0.; /* integrated P1CP */
	obs[37]=0.; /* integrated P2CP */
	obs[38]=0.; /* integrated P3CP */
	obs[39]=0.; /* integrated P4primeCP */
	obs[40]=0.; /* integrated P5primeCP */
	obs[41]=0.; /* integrated P6primeCP */
	obs[42]=0.; /* integrated P8primeCP */
	obs[43]=0.; /* integrated S6c */
	obs[44]=0.; /* integrated A6s */
	obs[45]=0.; /* integrated ACP */
	obs[46]=0.; /* integrated S1c */
	obs[47]=0.; /* integrated S2c */
	obs[48]=0.; /* integrated S2s */
	obs[49]=0.; /* integrated AFL */
	obs[50]=0.; /* integrated A1c */
	obs[51]=0.; /* integrated A2s */
	obs[52]=0.; /* integrated AFBCP */
	obs[53]=0.; /* integrated A6c */

	double dobs[Nobs_BKsll+1][3],dAFBtmp;
	double s0m,s0p,s0;
		
	dAFBtmp=0.;
	s0=s0m=0.;
	s0p=1.;
	
	double h=(smax-smin)/nmax;	
	s=smin;
	Gamma=dGamma_BKstarll_dq2_soft(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKsll;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}
	
	for(ie=1;ie<nmax;ie++)	
	{
		s+=h;

		Gamma+=4.*dGamma_BKstarll_dq2_soft(gen,charge,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_BKsll;je++) 
		{
			obs_num[je]+=4.*dobs[je][1];
			obs_den[je]+=4.*dobs[je][2];
		}

		if(ie>1){if(dAFBtmp/dobs[1][0]<0.) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);}

		Gamma+=2.*dGamma_BKstarll_dq2_soft(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_BKsll;je++) 
		{
			obs_num[je]+=2.*dobs[je][1];
			obs_den[je]+=2.*dobs[je][2];
		}

		if(dAFBtmp/dobs[1][0]<0.) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);
	}
	
	s=smax;
	Gamma+=4.*dGamma_BKstarll_dq2_soft(gen,charge,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKsll;je++) 
	{
		obs_num[je]+=4.*dobs[je][1];
		obs_den[je]+=4.*dobs[je][2];
	}	
	Gamma+=dGamma_BKstarll_dq2_soft(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKsll;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}	
	
	Gamma*=h/6.;
	for(je=1;je<=Nobs_BKsll;je++) 
	{
		obs_num[je]*=h/6.;
		obs_den[je]*=h/6.;
	}
	
	obs[0]=s0;
	for(je=1;je<=Nobs_BKsll;je++) 
	if((je==17)||(je==18)||(je==19)||(je==21)||(je==39)||(je==40)||(je==41)||(je==42))
	{
		obs[je]=obs_num[je]/sqrt(fabs(obs_den[je]*obs_den[15]/4.));
	}
	else if(je==11)
	{
		obs[je]=obs_num[je]/sqrt(fabs(obs_den[5]*obs_den[5]-obs_den[je]*obs_den[je]));
	}
	else if(je==20)
	{
		obs[je]=obs_num[je]/sqrt(fabs(obs_den[je]*obs_den[19]));
	}
	else obs[je]=obs_num[je]/obs_den[je];
	
	for(je=1;je<=Nobs_BKsll;je++) if(fabs(obs[je])<1.e-15) obs[je]=0.;
	
	if(charge==0) return param->life_Bd/hbar*Gamma;
	else return param->life_B/hbar*Gamma;
}

/*----------------------------------------------------------------------*/

double BRBKstarll_lowq2_soft(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKstarll_soft(gen,charge,1.,6.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarll_highq2_soft(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKstarll_soft(gen,charge,14.18,16.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_BKstarmumu_lowq2_soft_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) and all the other observables */

	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarll_lowq2_soft(2,0,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_BKstarmumu_highq2_soft_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) and all the other observables */
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarll_highq2_soft(2,0,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/
/*------------------------------ FULL ----------------------------------*/
/*----------------------------------------------------------------------*/

double dGamma_BKstarll_dq2_full(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{	
	double ml;
	if(gen==1) ml=param->mass_e;
	else if(gen==3) ml=param->mass_tau;
	else ml=param->mass_mu;
	
	double mB,mKstar,eq;
	if(charge==0) 
	{
		mB=param->m_Bd;
		mKstar=param->m_Kstar0;
		eq=-1./3.;
	}
	else
	{
		mB=param->m_B;
		mKstar=param->m_Kstar;
		eq=2./3.;
	}

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

	double complex Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_Kstar=(mB*mB+mKstar*mKstar-q2)/2./mB;

	int nf=5;
	double f_Kstar_perp=param->f_Kstar_perp;
	f_Kstar_perp*=pow(eta,4./3./(11.-2./3.*nf));

	double f_Kstar_par=param->f_Kstar_par;
	
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

	/********Defining deltaA***************/
	double complex deltaAperp=0.;
	double complex deltaApar=0.;
	double complex deltaA0=0.;
	double complex deltaAperp_bar=0.;
	double complex deltaApar_bar=0.;
	double complex deltaA0_bar=0.;
	/***********************************/
	
	
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

	/********LCSR+Lattice fit parameters from Barucha,Straub,Zwicky 1503.05534***************/
	if(param->BKstar_FormFactor_choice ==1)
	{
		MB = mB;
		MV = mKstar;
		a0V = param->a0V_BKstar;
		a1V = param->a1V_BKstar;
		a2V = param->a2V_BKstar;
		a0A0 = param->a0A0_BKstar;
		a1A0 = param->a1A0_BKstar;
		a2A0 = param->a2A0_BKstar;
		a0A1 = param->a0A1_BKstar;
		a1A1 = param->a1A1_BKstar;
		a2A1 = param->a2A1_BKstar;
		a0A12 = param->a0A12_BKstar;
		a1A12 = param->a1A12_BKstar;
		a2A12 = param->a2A12_BKstar;
		a0T1 = param->a0T1_BKstar;
		a1T1 = param->a1T1_BKstar;
		a2T1 = param->a2T1_BKstar;
		a0T2 = param->a0T2_BKstar;
		a1T2 = param->a1T2_BKstar;
		a2T2 = param->a2T2_BKstar;
		a0T23 = param->a0T23_BKstar;
		a1T23 = param->a1T23_BKstar;
		a2T23 = param->a2T23_BKstar;

		double P_V=1./(1.-q2/pow(param->MV_BKstar,2.));
		double P_A0=1./(1.-q2/pow(param->MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(param->MA1_BKstar,2.));
		double P_A12=1./(1.-q2/pow(param->MA12_BKstar,2.));
		double P_T1=1./(1.-q2/pow(param->MT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(param->MT2_BKstar,2.));
		double P_T23=1./(1.-q2/pow(param->MT23_BKstar,2.));

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


	/********Fit parameters from Gubernari, Reboud, van Dyk, Virto, 2305.06301***************/
	if(param->BKstar_FormFactor_choice ==2)
	{
		MB = mB;
		MV = mKstar;
		double MVT1_BKstar = 5.415;
		double MA0_BKstar = 5.367;
		double MAT_BKstar = 5.829;
		a0V = param->a0V_BKstar_GRvDV_BSZ;
		a1V = param->a1V_BKstar_GRvDV_BSZ;
		a2V = param->a2V_BKstar_GRvDV_BSZ;
		a0A0 = param->a0A0_BKstar_GRvDV_BSZ;
		a1A0 = param->a1A0_BKstar_GRvDV_BSZ;
		a2A0 = param->a2A0_BKstar_GRvDV_BSZ;
		a0A1 = param->a0A1_BKstar_GRvDV_BSZ;
		a1A1 = param->a1A1_BKstar_GRvDV_BSZ;
		a2A1 = param->a2A1_BKstar_GRvDV_BSZ;
		a0A12 = (MB*MB - MV*MV)/(8.*MB*MV)*param->a0A0_BKstar_GRvDV_BSZ;// Eq.17 in BSZ [1503.05534]
		a1A12 = param->a1A12_BKstar_GRvDV_BSZ;
		a2A12 = param->a2A12_BKstar_GRvDV_BSZ;
		a0T1 = param->a0T1_BKstar_GRvDV_BSZ;
		a1T1 = param->a1T1_BKstar_GRvDV_BSZ;
		a2T1 = param->a2T1_BKstar_GRvDV_BSZ;
		a0T2 = param->a0T1_BKstar_GRvDV_BSZ;// Eq.17 in BSZ [1503.05534]
		a1T2 = param->a1T2_BKstar_GRvDV_BSZ;
		a2T2 = param->a2T2_BKstar_GRvDV_BSZ;
		a0T23 = param->a0T23_BKstar_GRvDV_BSZ;
		a1T23 = param->a1T23_BKstar_GRvDV_BSZ;
		a2T23 = param->a2T23_BKstar_GRvDV_BSZ;

		double P_V=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_A0=1./(1.-q2/pow(MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_A12=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T1=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T23=1./(1.-q2/pow(MAT_BKstar,2.));

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


	/********LCSR+Lattice fit parameters from Gubernari van Dyk Kokulu 1811.00983***************/
	if(param->BKstar_FormFactor_choice == 3)
	{
		MB = mB;
		MV = mKstar;
		double MVT1_BKstar = 5.412;
		double MA0_BKstar = 5.336;
		double MAT_BKstar = 5.829;
		a0V = param->a0V_BKstar_GKvD_LCSR_Lattice;
		a1V = param->a1V_BKstar_GKvD_LCSR_Lattice;
		a2V = param->a2V_BKstar_GKvD_LCSR_Lattice;
		a0A0 = param->a0A0_BKstar_GKvD_LCSR_Lattice;
		a1A0 = param->a1A0_BKstar_GKvD_LCSR_Lattice;
		a2A0 = param->a2A0_BKstar_GKvD_LCSR_Lattice;
		a0A1 = param->a0A1_BKstar_GKvD_LCSR_Lattice;
		a1A1 = param->a1A1_BKstar_GKvD_LCSR_Lattice;
		a2A1 = param->a2A1_BKstar_GKvD_LCSR_Lattice;
		a0A12 = (MB*MB - MV*MV)/(8.*MB*MV)*param->a0A0_BKstar_GKvD_LCSR_Lattice;// Eq.17 in BSZ [1503.05534]
		a1A12 = param->a1A12_BKstar_GKvD_LCSR_Lattice;
		a2A12 = param->a2A12_BKstar_GKvD_LCSR_Lattice;
		a0T1 = param->a0T1_BKstar_GKvD_LCSR_Lattice;
		a1T1 = param->a1T1_BKstar_GKvD_LCSR_Lattice;
		a2T1 = param->a2T1_BKstar_GKvD_LCSR_Lattice;
		a0T2 = param->a0T1_BKstar_GKvD_LCSR_Lattice;// Eq.17 in BSZ [1503.05534]
		a1T2 = param->a1T2_BKstar_GKvD_LCSR_Lattice;
		a2T2 = param->a2T2_BKstar_GKvD_LCSR_Lattice;
		a0T23 = param->a0T23_BKstar_GKvD_LCSR_Lattice;
		a1T23 = param->a1T23_BKstar_GKvD_LCSR_Lattice;
		a2T23 = param->a2T23_BKstar_GKvD_LCSR_Lattice;

		double P_V=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_A0=1./(1.-q2/pow(MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_A12=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T1=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T23=1./(1.-q2/pow(MAT_BKstar,2.));

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


	/********LCSR-only fit parameters from Gubernari van Dyk Kokulu 1811.00983***************/
	if(param->BKstar_FormFactor_choice == 4)
	{
		MB = mB;
		MV = mKstar;
		double MVT1_BKstar = 5.412;
		double MA0_BKstar = 5.336;
		double MAT_BKstar = 5.829;
		a0V = param->a0V_BKstar_GKvD_LCSRonly;
		a1V = param->a1V_BKstar_GKvD_LCSRonly;
		a2V = param->a2V_BKstar_GKvD_LCSRonly;
		a0A0 = param->a0A0_BKstar_GKvD_LCSRonly;
		a1A0 = param->a1A0_BKstar_GKvD_LCSRonly;
		a2A0 = param->a2A0_BKstar_GKvD_LCSRonly;
		a0A1 = param->a0A1_BKstar_GKvD_LCSRonly;
		a1A1 = param->a1A1_BKstar_GKvD_LCSRonly;
		a2A1 = param->a2A1_BKstar_GKvD_LCSRonly;
		a0A12 = (MB*MB - MV*MV)/(8.*MB*MV)*param->a0A0_BKstar_GKvD_LCSRonly;// Eq.17 in BSZ [1503.05534]
		a1A12 = param->a1A12_BKstar_GKvD_LCSRonly;
		a2A12 = param->a2A12_BKstar_GKvD_LCSRonly;
		a0T1 = param->a0T1_BKstar_GKvD_LCSRonly;
		a1T1 = param->a1T1_BKstar_GKvD_LCSRonly;
		a2T1 = param->a2T1_BKstar_GKvD_LCSRonly;
		a0T2 = param->a0T1_BKstar_GKvD_LCSRonly;// Eq.17 in BSZ [1503.05534]
		a1T2 = param->a1T2_BKstar_GKvD_LCSRonly;
		a2T2 = param->a2T2_BKstar_GKvD_LCSRonly;
		a0T23 = param->a0T23_BKstar_GKvD_LCSRonly;
		a1T23 = param->a1T23_BKstar_GKvD_LCSRonly;
		a2T23 = param->a2T23_BKstar_GKvD_LCSRonly;

		double P_V=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_A0=1./(1.-q2/pow(MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_A12=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T1=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T23=1./(1.-q2/pow(MAT_BKstar,2.));

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


	/********LCSR-only fit parameters from Zwicky et al. 1503.05534***************/
	if(param->BKstar_FormFactor_choice ==5)
	{
		MB = mB;
		MV = mKstar;
		a0V = param->a0V_BKstar_BSZ_LCSRonly;
		a1V = param->a1V_BKstar_BSZ_LCSRonly;
		a2V = param->a2V_BKstar_BSZ_LCSRonly;
		a0A0 = param->a0A0_BKstar_BSZ_LCSRonly;
		a1A0 = param->a1A0_BKstar_BSZ_LCSRonly;
		a2A0 = param->a2A0_BKstar_BSZ_LCSRonly;
		a0A1 = param->a0A1_BKstar_BSZ_LCSRonly;
		a1A1 = param->a1A1_BKstar_BSZ_LCSRonly;
		a2A1 = param->a2A1_BKstar_BSZ_LCSRonly;
		a0A12 = param->a0A12_BKstar_BSZ_LCSRonly;
		a1A12 = param->a1A12_BKstar_BSZ_LCSRonly;
		a2A12 = param->a2A12_BKstar_BSZ_LCSRonly;
		a0T1 = param->a0T1_BKstar_BSZ_LCSRonly;
		a1T1 = param->a1T1_BKstar_BSZ_LCSRonly;
		a2T1 = param->a2T1_BKstar_BSZ_LCSRonly;
		a0T2 = param->a0T2_BKstar_BSZ_LCSRonly;
		a1T2 = param->a1T2_BKstar_BSZ_LCSRonly;
		a2T2 = param->a2T2_BKstar_BSZ_LCSRonly;
		a0T23 = param->a0T23_BKstar_BSZ_LCSRonly;
		a1T23 = param->a1T23_BKstar_BSZ_LCSRonly;
		a2T23 = param->a2T23_BKstar_BSZ_LCSRonly;

		double P_V=1./(1.-q2/pow(param->MV_BKstar,2.));
		double P_A0=1./(1.-q2/pow(param->MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(param->MA1_BKstar,2.));
		double P_A12=1./(1.-q2/pow(param->MA12_BKstar,2.));
		double P_T1=1./(1.-q2/pow(param->MT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(param->MT2_BKstar,2.));
		double P_T23=1./(1.-q2/pow(param->MT23_BKstar,2.));

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
	if(param->BKstar_FormFactor_choice == 6)
	{
		MB = mB;
		MV = mKstar;
		double MVT1_BKstar = 5.414;//Mres=mB+deltam = 5.279 + 0.135
		double MA0_BKstar = 5.366;//Mres=mB+deltam = 5.279 + 0.087
		double MAT_BKstar = 5.829;//Mres=mB+deltam = 5.279 + 0.550
		a0V   = param->a0V_BKstar_HLMW;
		a1V   = param->a1V_BKstar_HLMW;
		a0A0  = param->a0A0_BKstar_HLMW;
		a1A0  = param->a1A0_BKstar_HLMW;
		a0A1  = param->a0A1_BKstar_HLMW;
		a1A1  = param->a1A1_BKstar_HLMW;
		a0A12 = param->a0A12_BKstar_HLMW;
		a1A12 = param->a1A12_BKstar_HLMW;
		a0T1  = param->a0T1_BKstar_HLMW;
		a1T1  = param->a1T1_BKstar_HLMW;
		a0T2  = param->a0T2_BKstar_HLMW;
		a1T2  = param->a1T2_BKstar_HLMW;
		a0T23 = param->a0T23_BKstar_HLMW;
		a1T23 = param->a1T23_BKstar_HLMW;
		double sysErr = param->systErr_BKstar_HLMW;

		double P_V=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_A0=1./(1.-q2/pow(MA0_BKstar,2.));
		double P_A1=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_A12=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T1=1./(1.-q2/pow(MVT1_BKstar,2.));
		double P_T2=1./(1.-q2/pow(MAT_BKstar,2.));
		double P_T23=1./(1.-q2/pow(MAT_BKstar,2.));

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

	double xi_par=(mB+mKstar)/2./E_Kstar*A1-(mB-mKstar)/mB*A2;
	double xi_perp=mB/(mB+mKstar)*V;
//
	/***param->BKstar_implementation=1 -> QCDf; implementation=2 -> van Dyk et al.; implementation=3 -> Khodjamirian; implementation=4 -> Hadronic Fit ****/
	
	double m_Jpsi=3.096907; /*PDG 2016 page 1371*/
	double m_psi2S=3.686097; /*PDG 2016 page 1419*/

	/* Needed in the case of BKstar_implementation:4 to check the uncertainty of the result of hadronic fit */
	double complex hplus_fit = 0.;
	double complex hminus_fit = 0.;
	double complex hzero_fit = 0.;



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

// 		double complex Cperpp0u=q2/2./mb/mB*Yu;
// 		double complex Cperpm0u=q2/2./mb/mB*Yu;
// 		//double complex Cparp0u=-mB/2./mb*Yu;
// 		double complex Cparm0u=-mB/2./mb*Yu;
		/*******Used in redefinition of Tauperpu**********/
		double complex Cperp0u=q2/2./mb/mB*Yu;
		double complex Cpar0u=-mB/2./mb*Yu;
		
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
		-q2/2./mb/mB*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;
		
		double complex Cparnf=(C2bar*F27+C8eff*F87
		+mB/2./mb*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;

		double complex Cperpnfu=(-C2bar*(F27+F27_u)
		-q2/2./mb/mB*(C2bar*(F29+F29_u)+2.*C1bar*(F19+F19_u+1./6.*(F29+F29_u))))/4.*3.;
		
		double complex Cparnfu=(C2bar*(F27+F27_u)
		+mB/2./mb*(C2bar*(F29+F29_u)+2.*C1bar*(F19+F19_u+1./6.*(F29+F29_u))))/4.*3.;
		
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
		double Xi_par=mKstar/E_Kstar;
	
		double eu=2./3.;
		double ed=-1./3.;
	
		double a1perp=param->a1perp;
		double a2perp=param->a2perp;
		double a1par=param->a1par;
		double a2par=param->a2par;

		a1perp*=pow(eta,4./3.*(4.*1./2.)/(11.-2./3.*nf));
		a2perp*=pow(eta,4./3.*(4.*(1./2.+1./3.))/(11.-2./3.*nf));

		a1par*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
		a2par*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

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


		double lambda_Bp=param->lambda_Bp;
		lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_f,2.))*(1.-2.*1.4);

		double omega0=2.*(mB-mb)/3.;
		double complex lambda_Bm=1./(exp(-q2/mB/omega0)/omega0*(-Ei(q2/mB/omega0)+I*pi));

		double phiKstar_perp,phiKstar_par;
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
		double deltatp=param->deltatp;
		double deltatm=param->deltatm;
		
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
			
			phiKstar_perp=phi_Kstar(u,a1perp,a2perp);
			tperp_mc=tperp_bkll(u,mc,q2,E_Kstar,param);
			tperp_mb=tperp_bkll(u,mb,q2,E_Kstar,param);
			tperp_0=tperp_bkll(u,0.,q2,E_Kstar,param);

			Tperpppnf=Tperpmpnf=-4.*ed*C8eff/(u+(1.-u)*q2/mB/mB)
			+mB/2./mb*(eu*tperp_mc*(C2bar+C4bar-C6bar)
			+ed*tperp_mb*(C3bar+C4bar-C6bar-4.*mb/mB*C5bar)
			+ed*tperp_0*C3bar);
				
			Tperppmnf=Tperpmmnf=0.;

			Tperppnfu=mB/2./mb*eu*(tperp_mc-tperp_0)*(Cmub[2]-Cmub[1]/6.);
		
			Tperppp=Tperppp0+alphas_muf*4./3./4./pi*(Tperpppf+Tperpppnf); 
			Tperppm=Tperppm0+alphas_muf*4./3./4./pi*(Tperppmf+Tperppmnf);
			Tperpmp=Tperpmp0+alphas_muf*4./3./4./pi*(Tperpmpf+Tperpmpnf);
			Tperpmm=Tperpmm0+alphas_muf*4./3./4./pi*(Tperpmmf+Tperpmmnf);
	
			Tperppu=alphas_muf*4./3./4./pi*Tperppnfu;
	
			int_perppp+=phiKstar_perp*Tperppp/n1/lambda_Bp; 
			int_perppm+=phiKstar_perp*Tperppm/n1/lambda_Bm;
			int_perpmp+=phiKstar_perp*Tperpmp/n1/lambda_Bp;
			int_perpmm+=phiKstar_perp*Tperpmm/n1/lambda_Bm;

			int_perppu+=phiKstar_perp*Tperppu/n1/lambda_Bp;


		/* Tpar */		

			phiKstar_par=phi_Kstar(u,a1par,a2par);
			tpar_mc=tpar_bkll(u,mc,q2,E_Kstar,param);
			tpar_mb=tpar_bkll(u,mb,q2,E_Kstar,param);
			tpar_0=tpar_bkll(u,0.,q2,E_Kstar,param);
		
			Tparpp0=Tparmp0=0.;
		
			Tparppf=0.;
			Tparmpf=0.;

			Tparppnf=Tparmpnf=mB/mb*(eu*tpar_mc*(C2bar+C4bar-C6bar)
			+ed*tpar_mb*(C3bar+C4bar-C6bar)
			+ed*tpar_0*C3bar);
				
			Tparpnfu=mB/mb*eu*(tpar_mc-tpar_0)*(Cmub[2]-Cmub[1]/6.);
	
			Tparp0u=0.;
			Tparpu=alphas_muf*4./3./4./pi*Tparpnfu;
	
	
			Tparpm0=Tparmm0=-eq*4.*mB/mb*(C3bar+3.*C4bar);
			
			Tparp0u=0.;
			
			if(charge==0) Tparm0u=0.; else Tparm0u=eq*12.*mB/mb*Cmub[2];
		
			Tparpmf=Tparmmf=0.;

			h_mc=h_bkll((1.-u)*mB*mB+u*q2,mc,mu_b);
			h_mb=h_bkll((1.-u)*mB*mB+u*q2,mbpole,mu_b);
			h_0=h_bkll((1.-u)*mB*mB+u*q2,0.,mu_b);

			Tparpmnf=Tparmmnf=eq*(8.*C8eff/((1.-u)+u*q2/mB/mB)
			+6.*mB/mb*(h_mc*(C2bar+C4bar+C6bar)
			+h_mb*(C3bar+C4bar+C6bar)
			+h_0*(C3bar+3.*C4bar+3.*C6bar)
			-8./27.*(C3bar-C5bar-15.*C6bar)));
	
			Tparpp=Tparpp0+alphas_muf*4./3./4./pi*(Tparppf+Tparppnf);
			Tparpm=Tparpm0+alphas_muf*4./3./4./pi*(Tparpmf+Tparpmnf);
			Tparmp=Tparmp0+alphas_muf*4./3./4./pi*(Tparmpf+Tparmpnf);
			Tparmm=Tparmm0+alphas_muf*4./3./4./pi*(Tparmmf+Tparmmnf);
			
			Tparmnfu=eq*(6.*mB/mb*(h_mc-h_0)*(Cmub[2]-Cmub[1]/6.));

			Tparpu=Tparp0u+alphas_muf*4./3./4./pi*Tparpnfu;
	
			Tparmu=Tparm0u+alphas_muf*4./3./4./pi*Tparmnfu;


			int_parpp+=(phiKstar_par*Tparpp/lambda_Bp)/n1;
			int_parpm+=(phiKstar_par*Tparpm/lambda_Bm)/n1;
			
			int_parmp+=(phiKstar_par*Tparmp/lambda_Bp)/n1;
			int_parmm+=(phiKstar_par*Tparmm/lambda_Bm)/n1;
			
			int_parpu+=(phiKstar_par*Tparpu/lambda_Bp)/n1;
			
			int_parmu+=(phiKstar_par*Tparmu/lambda_Bm)/n1;


			integ3+=phiKstar_perp/((1.-u)+u*q2/mB/mB)/n1;

			x=(1.-u)*mB*mB+u*q2;
			
			h_mc=h_bkll(x,mc,mu_b);
			h_mb=h_bkll(x,mbpole,mu_b);
			h_0=h_bkll(x,0.,mu_b);
			
			FV=3./4.*(h_mc*(C2bar+C4bar+C6bar)+h_mb*(C3bar+C4bar+C6bar)+h_0*(C3bar+3.*C4bar+3.*C6bar)-8./27.*(C3bar-C5bar-15.*C6bar));

			FVu=3./4.*(h_mc-h_0)*(Cmub[2]-Cmub[1]/6.);
		

			integ4+=phiKstar_perp/((1.-u)+u*q2/mB/mB)*FV/n1;
			integ4u+=phiKstar_perp/((1.-u)+u*q2/mB/mB)*FVu/n1;
		
			Fperp+=phiKstar_perp/((1.-u)+u*q2/mB/mB)/3./n1;
			Xperp+=(u<=1.-0.5/mB)*phiKstar_perp/pow((1.-u)+u*q2/mB/mB,2.)/3./n1;
		
			integ5+=((3./4.*(1.+pow(2.*u-1.,2.))+a1par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1par*(2.*u-1.)+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1par*u+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FV/n1;

			integ5u+=((3./4.*(1.+pow(2.*u-1.,2.))+a1par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1par*(2.*u-1.)+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1par*u+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FVu/n1;

			if(ie==0||ie==n1sav) n1=n1sav;
		}
		
		/* Tau_perp */		

		double complex Tauperpp=xi_perp*Cperpp+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*(int_perppp+int_perppm); 
		
		double complex Tauperpm=xi_perp*Cperpm+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*(int_perpmp+int_perpmm); 
		
// 		double complex Tauperppu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_perp*Cperppu+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*int_perppu);
// 		double complex Tauperppu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_perp*Cperppu+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*int_perppu);
// 		
// 		double complex Tauperpmu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*xi_perp*Cperpmu;
// 		double complex Tauperpmu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*xi_perp*Cperpmu;

		double complex Tauperpu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_perp*Cperpu+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*int_perppu);
		double complex Tauperpu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_perp*Cperpu+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*int_perppu);

		
		/* Tau_par */		

		//double complex Tauparp=xi_par*Cparp+pi*pi/3.*param->f_B*f_Kstar_par/mB*Xi_par*(int_parpp+int_parpm);
		
		double complex Tauparm=xi_par*Cparm+pi*pi/3.*param->f_B*f_Kstar_par/mB*Xi_par*(int_parmp+int_parmm);
		
		double complex Tauparu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_par*Cparu+pi*pi/3.*param->f_B*f_Kstar_par/mB*Xi_par*(int_parpu+int_parmu));
		double complex Tauparu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_par*Cparu+pi*pi/3.*param->f_B*f_Kstar_par/mB*Xi_par*(int_parpu+int_parmu));
				
				
		double complex DeltaTauperpWA=-eq*4.*pi*pi/3.*param->f_B*f_Kstar_perp/mb/mB*(Cmub[3]+4./3.*Cmub[4]+4.*Cmub[5]+16./3.*Cmub[6])*integ3
		+eq*2.*pi*pi/3.*param->f_B*f_Kstar_par/mb/mB*mKstar/(1.-q2/mB/mB)/lambda_Bp*(Cmub[3]+4./3.*Cmub[4]+16.*Cmub[5]+64./3.*Cmub[6]);

		double complex DeltaTauperpuWA=0.;
		if(charge!=0) DeltaTauperpuWA=-eq*2.*pi*pi*param->f_B*f_Kstar_par/mb/mB*mKstar/(1.-q2/mB/mB)/lambda_Bp*Cmub[2];


		double rho=0.;
		double phi=0.;
		Xperp=Fperp+(1.+rho*(cos(phi)+I*sin(phi)))*Xperp;


		double complex DeltaTauperpHSA=eq*4./3.*alphas_muf/4./pi*pi*pi*param->f_B/3./mb/mB*(12.*C8eff*mb/mB*f_Kstar_perp*Xperp
		+8.*f_Kstar_perp*integ4-4.*mKstar*f_Kstar_par/(1.-q2/mB/mB)/lambda_Bp*integ5);
		
		double complex DeltaTauperpuHSA=eq*4./3.*alphas_muf/4./pi*pi*pi*param->f_B/3./mb/mB*
		(8.*f_Kstar_perp*integ4u-4.*mKstar*f_Kstar_par/(1.-q2/mB/mB)/lambda_Bp*integ5u);

		Tauperpp+=DeltaTauperpWA+DeltaTauperpHSA;
		Tauperpm+=DeltaTauperpWA+DeltaTauperpHSA;

// 		Tauperppu+=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(DeltaTauperpuWA+DeltaTauperpuHSA);
// 		Tauperppu_bar+=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(DeltaTauperpuWA+DeltaTauperpuHSA);	
// 		
// 		Tauperpmu+=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(DeltaTauperpuWA+DeltaTauperpuHSA);
// 		Tauperpmu_bar+=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(DeltaTauperpuWA+DeltaTauperpuHSA);		

		/*******Using Tauperpu**********/
		Tauperpu+=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(DeltaTauperpuWA+DeltaTauperpuHSA);
		Tauperpu_bar+=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(DeltaTauperpuWA+DeltaTauperpuHSA);    

		double lambda=pow(mB,4.)+pow(mKstar,4.)+q2*q2-2.*(mB*mB*mKstar*mKstar+mKstar*mKstar*q2+mB*mB*q2);
	
		double complex N=param->Vtb*conj(param->Vts)*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)/pow(mB,3.)*q2*sqrt(lambda)*beta_l);
		
		double complex Nbar=conj(param->Vtb)*param->Vts*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)/pow(mB,3.)*q2*sqrt(lambda)*beta_l);

		/* hadronic uncertainties */
		double ALperp_err=1.+param->BtoKstarlow_ALperp_err_noq2+q2/6.*param->BtoKstarlow_ALperp_err_q2;
		double ARperp_err=1.+param->BtoKstarlow_ARperp_err_noq2+q2/6.*param->BtoKstarlow_ARperp_err_q2;
		double ALpar_err=1.+param->BtoKstarlow_ALpar_err_noq2+q2/6.*param->BtoKstarlow_ALpar_err_q2;
		double ARpar_err=1.+param->BtoKstarlow_ARpar_err_noq2+q2/6.*param->BtoKstarlow_ARpar_err_q2;
		double AL0_err=1.+param->BtoKstarlow_AL0_err_noq2+q2/6.*param->BtoKstarlow_AL0_err_q2;
		double AR0_err=1.+param->BtoKstarlow_AR0_err_noq2+q2/6.*param->BtoKstarlow_AR0_err_q2;
		double At_err=1.+param->BtoKstarlow_At_err_noq2+q2/6.*param->BtoKstarlow_At_err_q2;
		double AS_err=1.+param->BtoKstarlow_AS_err_noq2+q2/6.*param->BtoKstarlow_AS_err_q2;
		double ALperp_bar_err=1.+param->BtoKstarlow_ALperp_err_noq2+q2/6.*param->BtoKstarlow_ALperp_err_q2;
		double ARperp_bar_err=1.+param->BtoKstarlow_ARperp_err_noq2+q2/6.*param->BtoKstarlow_ARperp_err_q2;
		double ALpar_bar_err=1.+param->BtoKstarlow_ALpar_err_noq2+q2/6.*param->BtoKstarlow_ALpar_err_q2;
		double ARpar_bar_err=1.+param->BtoKstarlow_ARpar_err_noq2+q2/6.*param->BtoKstarlow_ARpar_err_q2;
		double AL0_bar_err=1.+param->BtoKstarlow_AL0_err_noq2+q2/6.*param->BtoKstarlow_AL0_err_q2;
		double AR0_bar_err=1.+param->BtoKstarlow_AR0_err_noq2+q2/6.*param->BtoKstarlow_AR0_err_q2;
		double At_bar_err=1.+param->BtoKstarlow_At_err_noq2+q2/6.*param->BtoKstarlow_At_err_q2;
		double AS_bar_err=1.+param->BtoKstarlow_AS_err_noq2+q2/6.*param->BtoKstarlow_AS_err_q2;

		/**** Expressions for Hadronic fit*****/
		double complex deltaAperpHad,deltaAparHad,deltaA0Had,deltaAperp_barHad,deltaApar_barHad,deltaA0_barHad;
		deltaAperpHad=deltaAparHad=deltaA0Had=deltaAperp_barHad=deltaApar_barHad=deltaA0_barHad=0.;




		double Fcoeff = 0.;
		double complex Nprimecoeff = 0.;
		double complex Nbarprimecoeff = 0.;
		double Vtilde_plus=0.;
		double Vtilde_minus=0.;
		double Vtilde_zero=0.;
		if(param->BKstar_implementation==1)
		{
			Fcoeff=sqrt(lambda)*beta_l*q2/(3.*pow(2.,5.)*pow(pi,3.)*pow(mB,3.));
			Nprimecoeff=(-4.*param->Gfermi)/sqrt(2.)*mB*alpha_em/4./pi*param->Vtb*conj(param->Vts);
			Nbarprimecoeff=(-4.*param->Gfermi)/sqrt(2.)*mB*alpha_em/4./pi*conj(param->Vtb)*param->Vts;

			Vtilde_plus=0.5*((1.+mKstar/mB)*A1-sqrt(lambda)/mB/(mB+mKstar)*V);
			Vtilde_minus=0.5*((1.+mKstar/mB)*A1+sqrt(lambda)/mB/(mB+mKstar)*V);
			Vtilde_zero=4.*mKstar/sqrt(q2)*A12;

			deltaAperpHad=I*sqrt(Fcoeff)/(2.*sqrt(2.))*(-I*Nprimecoeff*mB*mB/q2)*(-16.*pi*pi)*((param->hplus0-param->hminus0)+(param->hplus1-param->hminus1)*q2+(param->hplus2-param->hminus2)*q2*q2 + 
			(param->DeltaC9plus*Vtilde_plus-param->DeltaC9minus*Vtilde_minus)*q2/(-16.*pi*pi*mB*mB));
			
			deltaAparHad=I*sqrt(Fcoeff)/(2.*sqrt(2.))*(-I*Nprimecoeff*mB*mB/q2)*(-16.*pi*pi)*((param->hplus0+param->hminus0)+(param->hplus1+param->hminus1)*q2+(param->hplus2+param->hminus2)*q2*q2 +
			(param->DeltaC9plus*Vtilde_plus+param->DeltaC9minus*Vtilde_minus)*q2/(-16.*pi*pi*mB*mB));
			deltaA0Had=I*sqrt(Fcoeff)/2.*(-I*Nprimecoeff*mB*mB/q2)*(-16.*pi*pi)*sqrt(q2)*(param->hzero0+param->hzero1*q2+param->hzero2*q2*q2 +
			(param->DeltaC9zero*Vtilde_zero)*sqrt(q2)/(-16.*pi*pi*mB*mB));
			deltaAperp_barHad=I*sqrt(Fcoeff)/(2.*sqrt(2.))*(-I*Nbarprimecoeff*mB*mB/q2)*(-16.*pi*pi)*((param->hplus0-param->hminus0)+(param->hplus1-param->hminus1)*q2+(param->hplus2-param->hminus2)*q2*q2 + 
			(param->DeltaC9plus*Vtilde_plus-param->DeltaC9minus*Vtilde_minus)*q2/(-16.*pi*pi*mB*mB));
			
			deltaApar_barHad=I*sqrt(Fcoeff)/(2.*sqrt(2.))*(-I*Nbarprimecoeff*mB*mB/q2)*(-16.*pi*pi)*((param->hplus0+param->hminus0)+(param->hplus1+param->hminus1)*q2+(param->hplus2+param->hminus2)*q2*q2 +
			(param->DeltaC9plus*Vtilde_plus+param->DeltaC9minus*Vtilde_minus)*q2/(-16.*pi*pi*mB*mB));
			deltaA0_barHad=I*sqrt(Fcoeff)/2.*(-I*Nbarprimecoeff*mB*mB/q2)*(-16.*pi*pi)*sqrt(q2)*(param->hzero0+param->hzero1*q2+param->hzero2*q2*q2 +
			(param->DeltaC9zero*Vtilde_zero)*sqrt(q2)/(-16.*pi*pi*mB*mB));
		}
		
		/**************************************************************************************************************************************/
 		ALperp=N*sqrt(2.)*sqrt(lambda)*(((C9+Y+C9p)-(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff+C7effp)*T1)+sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu)*(ALperp_err)+deltaAperpHad;
		
		ARperp=N*sqrt(2.)*sqrt(lambda)*(((C9+Y+C9p)+(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff+C7effp)*T1)+sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu)*(ARperp_err)+deltaAperpHad;
		
		ALpar=-N*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((C9+Y-C9p)-(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff-C7effp)*T2)-sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mB*mB-mKstar*mKstar)*2.*E_Kstar/mB*(Tauperpm+Tauperpu)*(ALpar_err)+deltaAparHad;
	
		ARpar=-N*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((C9+Y-C9p)+(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff-C7effp)*T2)-sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mB*mB-mKstar*mKstar)*2.*E_Kstar/mB*(Tauperpm+Tauperpu)*(ARpar_err)+deltaAparHad;
			
		AL0=-N/2./mKstar/sqrt(q2)*(((C9+Y-C9p)-(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))-N*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*((mB*mB+3.*mKstar*mKstar-q2)*2.*E_Kstar/mB-lambda/(mB*mB-mKstar*mKstar))*(Tauperpm+Tauperpu)*(AL0_err)+N*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*lambda/(mB*mB-mKstar*mKstar)*(Tauparm+Tauparu)*(AL0_err)+deltaA0Had;
		
		AR0=-N/2./mKstar/sqrt(q2)*(((C9+Y-C9p)+(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))-N*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*((mB*mB+3.*mKstar*mKstar-q2)*2.*E_Kstar/mB-lambda/(mB*mB-mKstar*mKstar))*(Tauperpm+Tauperpu)*(AR0_err)+N*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*lambda/(mB*mB-mKstar*mKstar)*(Tauparm+Tauparu)*(AR0_err)+deltaA0Had;
	
		At=N/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/ml*(CQ2-CQ2p)/(mb_mub+param->mass_s))*A0;
	
		AS=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(mb_mub+param->mass_s)*A0;

		/***** The barred amplitudes are defined with conjugated Wilson coefficietns ************/
		ALperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(((conj(C9)+Y+conj(C9p))-conj(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff+C7effp)*T1)+sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu_bar)*(ALperp_bar_err)+deltaAperp_barHad;

		ARperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(((conj(C9)+Y+conj(C9p))+conj(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff+C7effp)*T1)+sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu_bar)*(ARperp_bar_err)+deltaAperp_barHad;

		ALpar_bar=-Nbar*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((conj(C9)+Y-conj(C9p))-conj(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff-C7effp)*T2)-sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mB*mB-mKstar*mKstar)*2.*E_Kstar/mB*(Tauperpm+Tauperpu_bar)*(ALpar_bar_err)+deltaApar_barHad;

		ARpar_bar=-Nbar*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((conj(C9)+Y-conj(C9p))+conj(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff-C7effp)*T2)-sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mB*mB-mKstar*mKstar)*2.*E_Kstar/mB*(Tauperpm+Tauperpu_bar)*(ARpar_bar_err)+deltaApar_barHad;

		AL0_bar=-Nbar/2./mKstar/sqrt(q2)*(((conj(C9)+Y-conj(C9p))-conj(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*conj(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))-Nbar*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*((mB*mB+3.*mKstar*mKstar-q2)*2.*E_Kstar/mB-lambda/(mB*mB-mKstar*mKstar))*(Tauperpm+Tauperpu_bar)*(AL0_bar_err)+Nbar*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*lambda/(mB*mB-mKstar*mKstar)*(Tauparm+Tauparu_bar)*(AL0_bar_err)+deltaA0_barHad;

		AR0_bar=-Nbar/2./mKstar/sqrt(q2)*(((conj(C9)+Y-C9p)+conj(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*conj(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))-Nbar*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*((mB*mB+3.*mKstar*mKstar-q2)*2.*E_Kstar/mB-lambda/(mB*mB-mKstar*mKstar))*(Tauperpm+Tauperpu_bar)*(AR0_bar_err)+Nbar*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*lambda/(mB*mB-mKstar*mKstar)*(Tauparm+Tauparu_bar)*(AR0_bar_err)+deltaA0_barHad;

		At_bar=Nbar/sqrt(q2)*sqrt(lambda)*(2.*conj(C10-C10p)+q2/ml*conj(CQ2-CQ2p)/(mb_mub+param->mass_s))*A0;

		AS_bar=-2.*Nbar*sqrt(lambda)*conj(CQ1-CQ1p)/(mb_mub+param->mass_s)*A0;







		/********{\cal F}_{\lambda} form factors as used in van Dyk et al. 1707.07305***************/
		if(param->BKstar_implementation==2)
		{
			double T3=((mB*mB-mKstar*mKstar)*(mB*mB+3.*mKstar*mKstar-q2)*T2-8.*mB*mKstar*mKstar*(mB-mKstar)*T23)/((pow(mB+mKstar,2.)-q2)*(pow(mB-mKstar,2.)-q2));
		
			double calF_perp   = sqrt(2.*lambda)/mB/(mB+mKstar)*V;
			double calF_perp_T = sqrt(2.*lambda)/mB/mB*T1;
			double calF_par    = sqrt(2.)*(mB+mKstar)/mB*A1;
			double calF_par_T  = sqrt(2.)*(mB*mB-mKstar*mKstar)/mB/mB*T2;
			double calF_0      = mKstar/sqrt(q2)*((mB*mB-mKstar*mKstar-q2)*(mB+mKstar)*(mB+mKstar)*A1-lambda*A2)/(2.*mKstar*mKstar*(mB+mKstar)*(mB+mKstar));
			double calF_zero_T = q2/(2.*mKstar*sqrt(q2)*mB*(mB+mKstar))*((mB*mB+3.*mKstar*mKstar-q2)*T2 - lambda/(mB*mB-mKstar*mKstar)*T3);
			
			/***Table 1 of 1707.07305***/
			double real_alpha_perp0=param->real_alpha_perp0;
			double real_alpha_perp1=param->real_alpha_perp1;
			double real_alpha_perp2=param->real_alpha_perp2;
		
			double real_alpha_par0=param->real_alpha_par0;
			double real_alpha_par1=param->real_alpha_par1;
			double real_alpha_par2=param->real_alpha_par2;
			
			double real_alpha_zero0=param->real_alpha_zero0;
			double real_alpha_zero1=param->real_alpha_zero1;
			
			double imag_alpha_perp0=param->imag_alpha_perp0;
			double imag_alpha_perp1=param->imag_alpha_perp1;
			double imag_alpha_perp2=param->imag_alpha_perp2;
		
			double imag_alpha_par0=param->imag_alpha_par0;
			double imag_alpha_par1=param->imag_alpha_par1;
			double imag_alpha_par2=param->imag_alpha_par2;

			double imag_alpha_zero0=param->imag_alpha_zero0;
			double imag_alpha_zero1=param->imag_alpha_zero1;


			double t_plus=4.*pow(param->m_D0,2.);
			double t_zero=t_plus-sqrt(t_plus*(t_plus-pow(m_psi2S,2.)));
			double complex z_q2=(sqrt(t_plus-q2)-sqrt(t_plus-t_zero))/(sqrt(t_plus-q2)+sqrt(t_plus-t_zero));
			double complex z_q2_0=(sqrt(t_plus)-sqrt(t_plus-t_zero))/(sqrt(t_plus)+sqrt(t_plus-t_zero));
			double complex z_Jpsi=(sqrt(t_plus-pow(m_Jpsi,2.))-sqrt(t_plus-t_zero))/(sqrt(t_plus-pow(m_Jpsi,2.))+sqrt(t_plus-t_zero));
			double complex z_psi2S=(sqrt(t_plus-pow(m_psi2S,2.))-sqrt(t_plus-t_zero))/(sqrt(t_plus-pow(m_psi2S,2.))+sqrt(t_plus-t_zero));
			
			double complex P_Hperp= (real_alpha_perp0+real_alpha_perp1*z_q2+real_alpha_perp2*pow(z_q2,2.))
			+I*(imag_alpha_perp0+imag_alpha_perp1*z_q2+imag_alpha_perp2*pow(z_q2,2.));

			double complex P_Hpar= (real_alpha_par0+real_alpha_par1*z_q2+real_alpha_par2*pow(z_q2,2.))
			+I*(imag_alpha_par0+imag_alpha_par1*z_q2+imag_alpha_par2*pow(z_q2,2.));

			double complex P_Hzero= (real_alpha_zero0*(z_q2-z_q2_0)+real_alpha_zero1*z_q2*(z_q2-z_q2_0))
			+I*(imag_alpha_zero0*(z_q2-z_q2_0)+imag_alpha_zero1*z_q2*(z_q2-z_q2_0));

			/********{\cal H}_{\lambda} of van Dyk et al. 1707.07305 describing the charm contributions***************/	
			double complex calHperp =(1.-z_q2*conj(z_Jpsi))/(z_q2-z_Jpsi)*(1.-z_q2*conj(z_psi2S))/(z_q2-z_psi2S)*P_Hperp*calF_perp;
			double complex calHpar  =(1.-z_q2*conj(z_Jpsi))/(z_q2-z_Jpsi)*(1.-z_q2*conj(z_psi2S))/(z_q2-z_psi2S)*P_Hpar*calF_perp;
			double complex calHzero =(1.-z_q2*conj(z_Jpsi))/(z_q2-z_Jpsi)*(1.-z_q2*conj(z_psi2S))/(z_q2-z_psi2S)*P_Hzero*calF_perp;
			
			/*****Contribution of van Dyk et al. to the transversity amplitudes***************/
			deltaAperp = -16.*pi*pi*N*mB*2.*mB*mB/q2*calHperp;
			deltaApar  =  16.*pi*pi*N*mB*2.*mB*mB/q2*calHpar;
			deltaA0    =  16.*pi*pi*N*(mB+mKstar)*2.*mB*mB/q2*calHzero;
			
			deltaAperp_bar = -16.*pi*pi*Nbar*mB*2.*mB*mB/q2*calHperp;
			deltaApar_bar  =  16.*pi*pi*Nbar*mB*2.*mB*mB/q2*calHpar;
			deltaA0_bar    =  16.*pi*pi*Nbar*(mB+mKstar)*2.*mB*mB/q2*calHzero;

			
			ALperp=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)-(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff+C7effp)*T1)+deltaAperp;
			
			ARperp=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)+(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff+C7effp)*T1)+deltaAperp;
			
			ALpar=-N*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((C9-C9p)-(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff-C7effp)*T2)+deltaApar;
		
			ARpar=-N*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((C9-C9p)+(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff-C7effp)*T2)+deltaApar;
				
			AL0=-N/2./mKstar/sqrt(q2)*(((C9-C9p)-(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))+deltaA0;
			
			AR0=-N/2./mKstar/sqrt(q2)*(((C9-C9p)+(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))+deltaA0;
		
			At=N/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/ml*(CQ2-CQ2p)/(mb_mub+param->mass_s))*A0;
		
			AS=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(mb_mub+param->mass_s)*A0;

			/************************************************************************************************/
			ALperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(conj((C9+C9p)-(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff+C7effp)*T1)+deltaAperp_bar;
		
			ARperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(conj((C9+C9p)+(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff+C7effp)*T1)+deltaAperp_bar;
		
			ALpar_bar=-Nbar*sqrt(2.)*(mB*mB-mKstar*mKstar)*(conj((C9-C9p)-(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff-C7effp)*T2)+deltaApar_bar;
		
			ARpar_bar=-Nbar*sqrt(2.)*(mB*mB-mKstar*mKstar)*(conj((C9-C9p)+(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff-C7effp)*T2)+deltaApar_bar;
		
			AL0_bar=-Nbar/2./mKstar/sqrt(q2)*(((C9-C9p)-(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*conj(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))+deltaA0_bar;
			
			AR0_bar=-Nbar/2./mKstar/sqrt(q2)*(conj((C9-C9p)+(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*conj(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))+deltaA0_bar;
			
			At_bar=Nbar/sqrt(q2)*sqrt(lambda)*(2.*conj(C10-C10p)+q2/ml*conj(CQ2-CQ2p)/(mb_mub+param->mass_s))*A0;
		
			AS_bar=-2.*Nbar*sqrt(lambda)*conj(CQ1-CQ1p)/(mb_mub+param->mass_s)*A0;
		}
		

		/* Khodjamirian et al. (the charm part of Y + soft gluon) phenomenological formula 7.14 of 1006.4945
		 * does not contain QCDfactorisation improved calculations of Beneke et al. */
		if(param->BKstar_implementation==3)
		{
			/*****Table 2 of 1006.4945****/
			double q2bar=1.;
		
			double DeltaC9_M1_q2bar = param->DeltaC9_M1_q2bar;
			double r1_M1            = param->r1_M1;
			double r2_M1            = param->r2_M1;
			
			double DeltaC9_M2_q2bar = param->DeltaC9_M2_q2bar;
			double r1_M2            = param->r1_M2;
			double r2_M2            = param->r2_M2;
			
			double DeltaC9_M3_q2bar = param->DeltaC9_M3_q2bar;
			double r1_M3            = param->r1_M3;
			double r2_M3            = param->r2_M3;
		
			double DeltaC9_M1=(r1_M1*(1.-q2bar/q2)+DeltaC9_M1_q2bar*q2bar/q2)/(1.+r2_M1*(q2bar-q2)/m_Jpsi/m_Jpsi);
			double DeltaC9_M2=(r1_M2*(1.-q2bar/q2)+DeltaC9_M2_q2bar*q2bar/q2)/(1.+r2_M2*(q2bar-q2)/m_Jpsi/m_Jpsi);
			double DeltaC9_M3=(r1_M3*(1.-q2bar/q2)+DeltaC9_M3_q2bar*q2bar/q2)/(1.+r2_M3*(q2bar-q2)/m_Jpsi/m_Jpsi);

			/*****Contribution of Khodjamirian et al. to the transversity amplitudes***************/
			double complex AperpKMPW =N*sqrt(2.)*sqrt(lambda)*DeltaC9_M1*V/(mB+mKstar);
			double complex AparKMPW  =-N*sqrt(2.)*(mB*mB-mKstar*mKstar)*DeltaC9_M2*A1/(mB-mKstar);
			double complex A0KMPW    =-N/2./mKstar/sqrt(q2)*(((mB*mB-mKstar*mKstar-q2)*(mB+mKstar)*A1*DeltaC9_M2-lambda*A2*DeltaC9_M3/(mB+mKstar)));

			double complex AperpKMPW_bar=Nbar*sqrt(2.)*sqrt(lambda)*DeltaC9_M1*V/(mB+mKstar);
			double complex AparKMPW_bar=-Nbar*sqrt(2.)*(mB*mB-mKstar*mKstar)*DeltaC9_M2*A1/(mB-mKstar);
			double complex A0KMPW_bar=-Nbar/2./mKstar/sqrt(q2)*(((mB*mB-mKstar*mKstar-q2)*(mB+mKstar)*A1*DeltaC9_M2-lambda*A2*DeltaC9_M3/(mB+mKstar)));
			
			deltaAperp = AperpKMPW;
			deltaApar  = AparKMPW;
			deltaA0    = A0KMPW;

			deltaAperp_bar = AperpKMPW_bar; 
			deltaApar_bar  = AparKMPW_bar;  
			deltaA0_bar    = A0KMPW_bar;

			
			ALperp=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)-(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff+C7effp)*T1)+deltaAperp;
			
			ARperp=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)+(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff+C7effp)*T1)+deltaAperp;
			
			ALpar=-N*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((C9-C9p)-(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff-C7effp)*T2)+deltaApar;
		
			ARpar=-N*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((C9-C9p)+(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff-C7effp)*T2)+deltaApar;
				
			AL0=-N/2./mKstar/sqrt(q2)*(((C9-C9p)-(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))+deltaA0;
			
			AR0=-N/2./mKstar/sqrt(q2)*(((C9-C9p)+(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))+deltaA0;
		
			At=N/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/ml*(CQ2-CQ2p)/(mb_mub+param->mass_s))*A0;
		
			AS=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(mb_mub+param->mass_s)*A0;

			
			/************************************************************************************************/
			ALperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(conj((C9+C9p)-(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff+C7effp)*T1)+deltaAperp_bar;
		
			ARperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(conj((C9+C9p)+(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff+C7effp)*T1)+deltaAperp_bar;
		
			ALpar_bar=-Nbar*sqrt(2.)*(mB*mB-mKstar*mKstar)*(conj((C9-C9p)-(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff-C7effp)*T2)+deltaApar_bar;
		
			ARpar_bar=-Nbar*sqrt(2.)*(mB*mB-mKstar*mKstar)*(conj((C9-C9p)+(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff-C7effp)*T2)+deltaApar_bar;
		
			AL0_bar=-Nbar/2./mKstar/sqrt(q2)*(conj((C9-C9p)-(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*conj(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))+deltaA0_bar;
			
			AR0_bar=-Nbar/2./mKstar/sqrt(q2)*(conj((C9-C9p)+(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*conj(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))+deltaA0_bar;
			
			At_bar=Nbar/sqrt(q2)*sqrt(lambda)*(2.*conj(C10-C10p)+q2/ml*conj(CQ2-CQ2p)/(mb_mub+param->mass_s))*A0;
		
			AS_bar=-2.*Nbar*sqrt(lambda)*conj(CQ1-CQ1p)/(mb_mub+param->mass_s)*A0;
		}
		

        /* Implementation to use result from B->K* Hadronic Fit */
		if(param->BKstar_implementation==4)
		{
			double Fcoeff=sqrt(lambda)*beta_l*q2/(3.*pow(2.,5.)*pow(pi,3.)*pow(mB,3.));
			double complex Nprimecoeff=(-4.*param->Gfermi)/sqrt(2.)*mB*alpha_em/4./pi*param->Vtb*conj(param->Vts);
			double complex Nbarprimecoeff=(-4.*param->Gfermi)/sqrt(2.)*mB*alpha_em/4./pi*conj(param->Vtb)*param->Vts;

			double Vtilde_plus=0.5*((1.+mKstar/mB)*A1-sqrt(lambda)/mB/(mB+mKstar)*V);
			double Vtilde_minus=0.5*((1.+mKstar/mB)*A1+sqrt(lambda)/mB/(mB+mKstar)*V);
			double Vtilde_zero=4.*mKstar/sqrt(q2)*A12;

            hplus_fit = (param->hplus0r+I*param->hplus0i) + (param->hplus1r+I*param->hplus1i)*q2 + (param->hplus2r+I*param->hplus2i)*q2*q2;
            hminus_fit =(param->hminus0r+I*param->hminus0i) + (param->hminus1r+I*param->hminus1i)*q2 + (param->hminus2r+I*param->hminus2i)*q2*q2;
            hzero_fit = (param->hzero0r+I*param->hzero0i)+(param->hzero1r+I*param->hzero1i)*q2+(param->hzero2r+I*param->hzero2i)*q2*q2;


			deltaAperpHad=I*sqrt(Fcoeff)/(2.*sqrt(2.))*(-I*Nprimecoeff*mB*mB/q2)*(-16.*pi*pi)*(hplus_fit-hminus_fit);
			deltaAparHad=I*sqrt(Fcoeff)/(2.*sqrt(2.))*(-I*Nprimecoeff*mB*mB/q2)*(-16.*pi*pi)*(hplus_fit+hminus_fit);
			deltaA0Had=I*sqrt(Fcoeff)/2.*(-I*Nprimecoeff*mB*mB/q2)*(-16.*pi*pi)*sqrt(q2)*(hzero_fit);

			deltaAperp_barHad=I*sqrt(Fcoeff)/(2.*sqrt(2.))*(-I*Nbarprimecoeff*mB*mB/q2)*(-16.*pi*pi)*(hplus_fit-hminus_fit);
			deltaApar_barHad=I*sqrt(Fcoeff)/(2.*sqrt(2.))*(-I*Nbarprimecoeff*mB*mB/q2)*(-16.*pi*pi)*(hplus_fit+hminus_fit);
            deltaA0_barHad=I*sqrt(Fcoeff)/2.*(-I*Nbarprimecoeff*mB*mB/q2)*(-16.*pi*pi)*sqrt(q2)*(hzero_fit);


            /**************************************************************************************************************************************/
			ALperp=N*sqrt(2.)*sqrt(lambda)*(((C9+Y+C9p)-(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff+C7effp)*T1)+sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu)+deltaAperpHad;

			ARperp=N*sqrt(2.)*sqrt(lambda)*(((C9+Y+C9p)+(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff+C7effp)*T1)+sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu)+deltaAperpHad;

			ALpar=-N*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((C9+Y-C9p)-(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff-C7effp)*T2)-sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mB*mB-mKstar*mKstar)*2.*E_Kstar/mB*(Tauperpm+Tauperpu)+deltaAparHad;

			ARpar=-N*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((C9+Y-C9p)+(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(C7eff-C7effp)*T2)-sqrt(2.)*N*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mB*mB-mKstar*mKstar)*2.*E_Kstar/mB*(Tauperpm+Tauperpu)+deltaAparHad;

			AL0=-N/2./mKstar/sqrt(q2)*(((C9+Y-C9p)-(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))-N*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*((mB*mB+3.*mKstar*mKstar-q2)*2.*E_Kstar/mB-lambda/(mB*mB-mKstar*mKstar))*(Tauperpm+Tauperpu)+N*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*lambda/(mB*mB-mKstar*mKstar)*(Tauparm+Tauparu)+deltaA0Had;

			AR0=-N/2./mKstar/sqrt(q2)*(((C9+Y-C9p)+(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))-N*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*((mB*mB+3.*mKstar*mKstar-q2)*2.*E_Kstar/mB-lambda/(mB*mB-mKstar*mKstar))*(Tauperpm+Tauperpu)+N*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*lambda/(mB*mB-mKstar*mKstar)*(Tauparm+Tauparu)+deltaA0Had;

			At=N/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/ml*(CQ2-CQ2p)/(mb_mub+param->mass_s))*A0;

			AS=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(mb_mub+param->mass_s)*A0;





			/***** The barred amplitudes are defined with conjugated Wilson coefficietns ************/
			ALperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(((conj(C9)+Y+conj(C9p))-conj(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff+C7effp)*T1)+sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu_bar)+deltaAperp_barHad;

			ARperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(((conj(C9)+Y+conj(C9p))+conj(C10+C10p))*V/(mB+mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff+C7effp)*T1)+sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*sqrt(lambda)*(Tauperpp+Tauperpu_bar)+deltaAperp_barHad;

			ALpar_bar=-Nbar*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((conj(C9)+Y-conj(C9p))-conj(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff-C7effp)*T2)-sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mB*mB-mKstar*mKstar)*2.*E_Kstar/mB*(Tauperpm+Tauperpu_bar)+deltaApar_barHad;

			ARpar_bar=-Nbar*sqrt(2.)*(mB*mB-mKstar*mKstar)*(((conj(C9)+Y-conj(C9p))+conj(C10-C10p))*A1/(mB-mKstar)+2.*(mb+alphas_mub/3./pi*DeltaM)/q2*conj(C7eff-C7effp)*T2)-sqrt(2.)*Nbar*2.*(mb+alphas_mub/3./pi*DeltaM)/q2*(mB*mB-mKstar*mKstar)*2.*E_Kstar/mB*(Tauperpm+Tauperpu_bar)+deltaApar_barHad;

			AL0_bar=-Nbar/2./mKstar/sqrt(q2)*(((conj(C9)+Y-conj(C9p))-conj(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*conj(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))-Nbar*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*((mB*mB+3.*mKstar*mKstar-q2)*2.*E_Kstar/mB-lambda/(mB*mB-mKstar*mKstar))*(Tauperpm+Tauperpu_bar)+Nbar*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*lambda/(mB*mB-mKstar*mKstar)*(Tauparm+Tauparu_bar)+deltaA0_barHad;

			AR0_bar=-Nbar/2./mKstar/sqrt(q2)*(((conj(C9)+Y-C9p)+conj(C10-C10p))*(16.*mB*mKstar*mKstar*A12)+2.*(mb+alphas_mub/3./pi*DeltaM)*conj(C7eff-C7effp)*(8.*mB*mKstar*mKstar/(mB+mKstar)*T23))-Nbar*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*((mB*mB+3.*mKstar*mKstar-q2)*2.*E_Kstar/mB-lambda/(mB*mB-mKstar*mKstar))*(Tauperpm+Tauperpu_bar)+Nbar*(mb+alphas_mub/3./pi*DeltaM)/mKstar/sqrt(q2)*lambda/(mB*mB-mKstar*mKstar)*(Tauparm+Tauparu_bar)+deltaA0_barHad;

			At_bar=Nbar/sqrt(q2)*sqrt(lambda)*(2.*conj(C10-C10p)+q2/ml*conj(CQ2-CQ2p)/(mb_mub+param->mass_s))*A0;

			AS_bar=-2.*Nbar*sqrt(lambda)*conj(CQ1-CQ1p)/(mb_mub+param->mass_s)*A0;
        }

		
	}
	
	if(q2>LOW_q2)
	{
		double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
		double mb_PS=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */ //mb=>4.64
		double sh=q2/mb_PS/mb_PS;

		double shat=q2/mB/mB;

		double lambda_hat=1.+shat*shat+pow(mKstar/mB,4.)-2.*(shat+shat*mKstar*mKstar/mB/mB+mKstar*mKstar/mB/mB);


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

		double complex N=param->Vtb*conj(param->Vts)*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em*beta_l/3./1024./pow(pi,5.)*mB*shat*sqrt(lambda_hat));// Missing beta_l corrected

		double complex Nbar=conj(param->Vtb)*param->Vts*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em*beta_l/3./1024./pow(pi,5.)*mB*shat*sqrt(lambda_hat));// Missing beta_l corrected

		double complex f_perp=N*mB*sqrt(2.*lambda_hat)/(1.+mKstar/mB)*V;

		double complex f_par=N*mB*sqrt(2.)*(1.+mKstar/mB)*A1;

		double complex f_0=N*mB*((1.-shat-mKstar*mKstar/mB/mB)*pow(1.+mKstar/mB,2.)*A1-lambda_hat*A2)/(2.*mKstar/mB*(1.+mKstar/mB)*sqrt(shat));

		double complex f_perp_bar=Nbar*mB*sqrt(2.*lambda_hat)/(1.+mKstar/mB)*V;

		double complex f_par_bar=Nbar*mB*sqrt(2.)*(1.+mKstar/mB)*A1;

		double complex f_0_bar=Nbar*mB*((1.-shat-mKstar*mKstar/mB/mB)*pow(1.+mKstar/mB,2.)*A1-lambda_hat*A2)/(2.*mKstar/mB*(1.+mKstar/mB)*sqrt(shat));

		ALperp_high=(((C9eff+C9p)-(Cmub[10]+C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff+C7effp))*f_perp;
		ARperp_high=(((C9eff+C9p)+(Cmub[10]+C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff+C7effp))*f_perp;

		ALpar_high=-(((C9eff-C9p)-(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff-C7effp))*f_par;
		ARpar_high=-(((C9eff-C9p)+(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff-C7effp))*f_par;

		AL0_high=-(((C9eff-C9p)-(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff-C7effp))*f_0;
		AR0_high=-(((C9eff-C9p)+(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff-C7effp))*f_0;

		/************************************************************************************************/
		ALperp_bar_high=(((C9eff_bar+conj(C9p))-conj(Cmub[10]+C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar+conj(C7effp)))*f_perp_bar;
		ARperp_bar_high=(((C9eff_bar+conj(C9p))+conj(Cmub[10]+C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar+conj(C7effp)))*f_perp_bar;

		ALpar_bar_high=-(((C9eff_bar-conj(C9p))-conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar-conj(C7effp)))*f_par_bar;
		ARpar_bar_high=-(((C9eff_bar-conj(C9p))+conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar-conj(C7effp)))*f_par_bar;

		AL0_bar_high=-(((C9eff_bar-conj(C9p))-conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar-conj(C7effp)))*f_0_bar;
		AR0_bar_high=-(((C9eff_bar-conj(C9p))+conj(Cmub[10]-C10p))+kappa*2.*param->mass_b/mB/shat*(C7eff_bar-conj(C7effp)))*f_0_bar;

		double complex C10=Cmub[10];
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];

		double lambda=pow(mB,4.)+pow(mKstar,4.)+q2*q2-2.*(mB*mB*mKstar*mKstar+mKstar*mKstar*q2+mB*mB*q2);

		At_high=N/sqrt(q2)*sqrt(lambda)*( 2.*(C10-C10p)+q2/ml*(CQ2-CQ2p)/(param->mass_b+param->mass_s))*A0;
		AS_high=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(param->mass_b+param->mass_s)*A0;

		At_bar_high=Nbar/sqrt(q2)*sqrt(lambda)*( 2.*conj(C10-C10p)+q2/ml*conj(CQ2-CQ2p)/(param->mass_b+param->mass_s))*A0;
		AS_bar_high=-2.*Nbar*sqrt(lambda)*conj(CQ1-CQ1p)/(param->mass_b+param->mass_s)*A0;


		/* hadronic uncertainties */
		ALperp_high*=1.+param->BtoKstarhigh_ALperp_err;
		ARperp_high*=1.+param->BtoKstarhigh_ARperp_err;
		ALpar_high*=1.+param->BtoKstarhigh_ALpar_err;
		ARpar_high*=1.+param->BtoKstarhigh_ARpar_err;
		AL0_high*=1.+param->BtoKstarhigh_AL0_err;
		AR0_high*=1.+param->BtoKstarhigh_AR0_err;
		At_high*=1.+param->BtoKstarhigh_At_err;
		AS_high*=1.+param->BtoKstarhigh_AS_err;
		ALperp_bar_high*=1.+param->BtoKstarhigh_ALperp_err;
		ARperp_bar_high*=1.+param->BtoKstarhigh_ARperp_err;
		ALpar_bar_high*=1.+param->BtoKstarhigh_ALpar_err;
		ARpar_bar_high*=1.+param->BtoKstarhigh_ARpar_err;
		AL0_bar_high*=1.+param->BtoKstarhigh_AL0_err;
		AR0_bar_high*=1.+param->BtoKstarhigh_AR0_err;
		At_bar_high*=1.+param->BtoKstarhigh_At_err;
		AS_bar_high*=1.+param->BtoKstarhigh_AS_err;


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
		
	double dGamma_BKstarll_dq2=3./4.*(2.*J1s+J1c-(2.*J2s+J2c)/3.);

	double dGamma_bar_BKstarll_dq2=3./4.*(2.*J1s_bar+J1c_bar-(2.*J2s_bar+J2c_bar)/3.);

	double AFB[3],FL[3],FT[3],AT1[3],AT2[3],AT3[3],AT4[3],AT5[3],HT1[3],HT2[3],HT3[3],alpha_Kstar[3],AIm[3],P2[3],P3[3],P6[3],P8[3],P4prime[3],P5prime[3],P6prime[3],P8prime[3],A3[3],A4[3],A5[3],A7[3],A8[3],A9[3],S3[3],S4[3],S5[3],S7[3],S8[3],S9[3];

	int scalar_switch = 1;// 0 if we want to put scalar to zero in AFB, else 1
	AFB[0]=-3./8.*(2.*(J6s+J6s_bar)+scalar_switch*(J6c+J6c_bar))/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	AFB[1]=-3./8.*(2.*(J6s+J6s_bar)+scalar_switch*(J6c+J6c_bar));
	AFB[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	FL[0]=1./4.*(3.*(J1c+J1c_bar)-(J2c+J2c_bar))/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	FL[1]=1./4.*(3.*(J1c+J1c_bar)-(J2c+J2c_bar));
	FL[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	FT[0]=4.*(J2s+J2s_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	FT[1]=4.*(J2s+J2s_bar);
	FT[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	AT1[0]=-2.*creal(ALpar*conj(ALperp)+ARpar*conj(ARperp)+ALpar_bar*conj(ALperp_bar)+ARpar_bar*conj(ARperp_bar))/(Apar2+Aperp2+Apar2_bar+Aperp2_bar);
	AT1[1]=-2.*creal(ALpar*conj(ALperp)+ARpar*conj(ARperp)+ALpar_bar*conj(ALperp_bar)+ARpar_bar*conj(ARperp_bar));
	AT1[2]=Apar2+Aperp2+Apar2_bar+Aperp2_bar;

	AT2[0]=(J3+J3_bar)/2./(J2s+J2s_bar);
	AT2[1]=J3+J3_bar;
	AT2[2]=2.*(J2s+J2s_bar);

	AT3[0]=sqrt(fabs((4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar))/(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar))));
	AT3[1]=sqrt(fabs(4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar)));
	AT3[2]=sqrt(fabs(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	
	AT4[0]=sqrt((beta_l*beta_l*(J5+J5_bar)*(J5+J5_bar)+4.*(J8+J8_bar)*(J8+J8_bar))/(4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar)));
	AT4[1]=sqrt(beta_l*beta_l*(J5+J5_bar)*(J5+J5_bar)+4.*(J8+J8_bar)*(J8+J8_bar));
	AT4[2]=sqrt(4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar));
	
	AT5[0]=cabs(ALperp*conj(ARpar)+ALpar*conj(ARperp)+ALperp_bar*conj(ARpar_bar)+ALpar_bar*conj(ARperp_bar))/(Apar2+Aperp2+Apar2_bar+Aperp2_bar);
	AT5[1]=cabs(ALperp*conj(ARpar)+ALpar*conj(ARperp)+ALperp_bar*conj(ARpar_bar)+ALpar_bar*conj(ARperp_bar));
	AT5[2]=Apar2+Aperp2+Apar2_bar+Aperp2_bar;
	
	HT1[0]=sqrt(2.)*(J4+J4_bar)/sqrt(fabs(-(J2c+J2c_bar)*(2.*J2s-J3+2.*J2s_bar-J3_bar)));
	HT1[1]=sqrt(2.)*(J4+J4_bar);
	HT1[2]=sqrt(fabs(-(J2c+J2c_bar)*(2.*J2s-J3+2.*J2s_bar-J3_bar)));
	
	HT2[0]=beta_l*(J5+J5_bar)/sqrt(fabs(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	HT2[1]=beta_l*(J5+J5_bar);
	HT2[2]=sqrt(fabs(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	
	HT3[0]=(J6s+J6s_bar)/2./sqrt(fabs(4.*(J2s+J2s_bar)*(J2s+J2s_bar)-(J3+J3_bar)*(J3+J3_bar)));
	HT3[1]=(J6s+J6s_bar)/2.;
	HT3[2]=J3+J3_bar;

	alpha_Kstar[0]=-(2.*J2s+J2c+2.*J2s_bar+J2c_bar)/2./(J2s+J2s_bar);
	alpha_Kstar[1]=-(2.*J2s+J2c+2.*J2s_bar+J2c_bar);
	alpha_Kstar[2]=2.*(J2s+J2s_bar);

	AIm[0]=(J9+J9_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	AIm[1]=J9+J9_bar;
	AIm[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	P2[0]=(J6s+J6s_bar)/8./(J2s+J2s_bar);// Definition is from 1303.5794 without beta_l and not from 1202.4266
	P2[1]=J6s+J6s_bar;
	P2[2]=8.*(J2s+J2s_bar);

	P3[0]=-(J9+J9_bar)/4./(J2s+J2s_bar);
	P3[1]=-(J9+J9_bar);
	P3[2]=4.*(J2s+J2s_bar);
	
	P6[0]=-beta_l*(J7+J7_bar)/sqrt(fabs(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	P6[1]=-beta_l*(J7+J7_bar);
	P6[2]=sqrt(fabs(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	
	P4prime[0]=(J4+J4_bar)/sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P4prime[1]=J4+J4_bar;
	P4prime[2]=-(J2c+J2c_bar);
	
	P5prime[0]=(J5+J5_bar)/2./sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P5prime[1]=(J5+J5_bar)/2.;
	P5prime[2]=-(J2c+J2c_bar);
	
	P6prime[0]=-(J7+J7_bar)/2./sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P6prime[1]=-(J7+J7_bar)/2.;
	P6prime[2]=-(J2c+J2c_bar);

	P8[0]=-sqrt(2.)*(J8+J8_bar)/sqrt(fabs(-(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	P8[1]=-sqrt(2.)*(J8+J8_bar);
	P8[2]=2.*J2s+J3+2.*J2s_bar+J3_bar;
	
	P8prime[0]=-(J8+J8_bar)/sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P8prime[1]=-(J8+J8_bar);
	P8prime[2]=-(J2c+J2c_bar);

	A3[0]=(J3-J3_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A3[1]=J3-J3_bar;
	A3[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	A4[0]=(J4-J4_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A4[1]=J4-J4_bar;
	A4[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	
	A5[0]=(J5-J5_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A5[1]=J5-J5_bar;
	A5[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	
	A7[0]=(J7-J7_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A7[1]=J7-J7_bar;
	A7[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	
	A8[0]=(J8-J8_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A8[1]=J8-J8_bar;
	A8[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	
	A9[0]=(J9-J9_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A9[1]=J9-J9_bar;
	A9[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	
	S3[0]=(J3+J3_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S3[1]=J3+J3_bar;
	S3[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S4[0]=(J4+J4_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S4[1]=J4+J4_bar;
	S4[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	
	S5[0]=(J5+J5_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S5[1]=J5+J5_bar;
	S5[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S7[0]=(J7+J7_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S7[1]=J7+J7_bar;
	S7[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S8[0]=(J8+J8_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S8[1]=J8+J8_bar;
	S8[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S9[0]=(J9+J9_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S9[1]=J9+J9_bar;
	S9[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	double ATRe[3], ATReCP[3], P1CP[3], P2CP[3], P3CP[3], P4primeCP[3], P5primeCP[3], P6primeCP[3], P8primeCP[3], S6c[3], A6s[3], ACP[3];

	ATRe[0]=beta_l*(J6s+J6s_bar)/(4.*(J2s+J2s_bar));//based on 1106.3283 which includes the beta_l function; ATRe = beta_l*2*P2
	ATRe[1]=beta_l*(J6s+J6s_bar);
	ATRe[2]=4.*(J2s+J2s_bar);

	ATReCP[0]=beta_l*(J6s-J6s_bar)/(4.*(J2s+J2s_bar));//based on 1106.3283 which includes the beta_l function; ATRe = beta_l*2*P2CP
	ATReCP[1]=beta_l*(J6s-J6s_bar);
	ATReCP[2]=4.*(J2s+J2s_bar);

	P1CP[0]=(J3-J3_bar)/(2.*(J2s+J2s_bar));
	P1CP[1]=(J3-J3_bar);
	P1CP[2]=2.*(J2s+J2s_bar);

	P2CP[0]=(J6s-J6s_bar)/8./(J2s+J2s_bar);// Definition from 1303.5794 without beta_l and not from 1202.4266
	P2CP[1]=J6s-J6s_bar;
	P2CP[2]=8.*(J2s+J2s_bar);

	P3CP[0]=-(J9-J9_bar)/4./(J2s+J2s_bar);
	P3CP[1]=-(J9-J9_bar);
	P3CP[2]=4.*(J2s+J2s_bar);

	P4primeCP[0]=(J4-J4_bar)/sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P4primeCP[1]=J4-J4_bar;
	P4primeCP[2]=-(J2c+J2c_bar);

	P5primeCP[0]=(J5-J5_bar)/2./sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P5primeCP[1]=(J5-J5_bar)/2.;
	P5primeCP[2]=-(J2c+J2c_bar);

	P6primeCP[0]=-(J7-J7_bar)/2./sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P6primeCP[1]=-(J7-J7_bar)/2.;
	P6primeCP[2]=-(J2c+J2c_bar);

	P8primeCP[0]=-(J8-J8_bar)/sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P8primeCP[1]=-(J8-J8_bar);
	P8primeCP[2]=-(J2c+J2c_bar);

	/* We use the same convention as AFB above regarding the sign */
	S6c[0]=-(J6c+J6c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S6c[1]=-(J6c+J6c_bar);
	S6c[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	/* We use the same convention as AFB above regarding the sign (no need for change of sign compared to LHCb) */
	A6s[0]=-(J6s-J6s_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A6s[1]=-(J6s-J6s_bar);
	A6s[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	ACP[0]=(dGamma_BKstarll_dq2-dGamma_bar_BKstarll_dq2)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	ACP[1]=(dGamma_BKstarll_dq2-dGamma_bar_BKstarll_dq2);
	ACP[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	double S1c[3], S2c[3], S2s[3], AFL[3], A1c[3], A2s[3], A6c[3], AFBCP[3];

	S1c[0]=(J1c+J1c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S1c[1]=(J1c+J1c_bar);
	S1c[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S2c[0]=(J2c+J2c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S2c[1]=(J2c+J2c_bar);
	S2c[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	S2s[0]=(J2s+J2s_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	S2s[1]=(J2s+J2s_bar);
	S2s[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	AFL[0]=-(J2c-J2c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	AFL[1]=-(J2c-J2c_bar);
	AFL[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	A1c[0]=(J1c-J1c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A1c[1]=(J1c-J1c_bar);
	A1c[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	A2s[0]=(J2s-J2s_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A2s[1]=(J2s-J2s_bar);
	A2s[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	/* See explanation above AFB for `scalar_switch` */
	AFBCP[0]=-3./8.*(2.*(J6s-J6s_bar)+scalar_switch*(J6c-J6c_bar))/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	AFBCP[1]=-3./8.*(2.*(J6s-J6s_bar)+scalar_switch*(J6c-J6c_bar));
	AFBCP[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;

	A6c[0]=-(J6c-J6c_bar)/(dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2);
	A6c[1]=-(J6c-J6c_bar);
	A6c[2]=dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2;
	/*********************************************/


	for(je=0;je<=Nobs_BKsll;je++) for(ie=0;ie<=2;ie++) obs[je][ie]=0.;
	
	for(ie=0;ie<=2;ie++)
	{
		obs[1][ie]=AFB[ie];
		obs[2][ie]=FL[ie];
		obs[3][ie]=FT[ie];
		obs[4][ie]=AT1[ie];
		obs[5][ie]=AT2[ie];
		obs[6][ie]=AT3[ie];
		obs[7][ie]=AT4[ie];
		obs[8][ie]=AT5[ie];
		obs[9][ie]=HT1[ie];
		obs[10][ie]=HT2[ie];
		obs[11][ie]=HT3[ie];
		obs[12][ie]=alpha_Kstar[ie];
		obs[13][ie]=AIm[ie];
		obs[14][ie]=P2[ie];
		obs[15][ie]=P3[ie];
		obs[16][ie]=P6[ie];
		obs[17][ie]=P4prime[ie];
		obs[18][ie]=P5prime[ie];
		obs[19][ie]=P6prime[ie];
		obs[20][ie]=P8[ie];
		obs[21][ie]=P8prime[ie];
		obs[22][ie]=A7[ie];
		obs[23][ie]=A8[ie];
		obs[24][ie]=A9[ie];
		obs[25][ie]=S3[ie];
		obs[26][ie]=S4[ie];
		obs[27][ie]=S5[ie];
		obs[28][ie]=S7[ie];
		obs[29][ie]=S8[ie];
		obs[30][ie]=S9[ie];
		obs[31][ie]=A3[ie];
		obs[32][ie]=A4[ie];
		obs[33][ie]=A5[ie];
		obs[34][ie]=ATRe[ie];
		obs[35][ie]=ATReCP[ie];
		obs[36][ie]=P1CP[ie];
		obs[37][ie]=P2CP[ie];
		obs[38][ie]=P3CP[ie];
		obs[39][ie]=P4primeCP[ie];
		obs[40][ie]=P5primeCP[ie];
		obs[41][ie]=P6primeCP[ie];
		obs[42][ie]=P8primeCP[ie];
		obs[43][ie]=S6c[ie];
		obs[44][ie]=A6s[ie];
		obs[45][ie]=ACP[ie];
		obs[46][ie]=S1c[ie];
		obs[47][ie]=S2c[ie];
		obs[48][ie]=S2s[ie];
		obs[49][ie]=AFL[ie];
		obs[50][ie]=A1c[ie];
		obs[51][ie]=A2s[ie];
		obs[52][ie]=AFBCP[ie];
		obs[53][ie]=A6c[ie];
		// If you intend to add physical observables, add them above and change Nobs_BKsll_phys in `include.c` accordingly
	}

	// If you intend to add non-physical observables, add them below after last line and change Nobs_BKsll_nonphys in `include.c` accordingly
	/* Form Factors */
	obs[Nobs_BKsll_phys+1][0] = V;
	obs[Nobs_BKsll_phys+2][0] = A0;
	obs[Nobs_BKsll_phys+3][0] = A1;
	obs[Nobs_BKsll_phys+4][0] = A12;
	obs[Nobs_BKsll_phys+5][0] = T1;
	obs[Nobs_BKsll_phys+6][0] = T2;
	obs[Nobs_BKsll_phys+7][0] = T23;
	obs[Nobs_BKsll_phys+8][0] = A2;
	/******/
	/* hadronic fit results */
	obs[Nobs_BKsll_phys+9][0] = creal(hplus_fit);
	obs[Nobs_BKsll_phys+10][0] = cimag(hplus_fit);
	obs[Nobs_BKsll_phys+11][0] = creal(hminus_fit);
	obs[Nobs_BKsll_phys+12][0] = cimag(hminus_fit);
	obs[Nobs_BKsll_phys+13][0] = creal(hzero_fit);
	obs[Nobs_BKsll_phys+14][0] = cimag(hzero_fit);

	return (dGamma_BKstarll_dq2+dGamma_bar_BKstarll_dq2)/2.;
}

/*----------------------------------------------------------------------*/

double dGamma_BKstarll_dq2_full_calculator(int gen, int charge, double q2, double obs[][3], char name[])
/* "container" function scanning the SLHA file "name" and calculating dGamma/dq2(B->Kstar mu+ mu-) */
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

	return dGamma_BKstarll_dq2_full(gen,charge,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}


/*----------------------------------------------------------------------*/

double BRBKstarll_full(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	int ie,je;
	int nmax=10;
	if((smin<0.099)||(smax-smin>10.)) nmax=100;
	double Gamma=0.;
	double s;
	
	double obs_num[Nobs_BKsll+1],obs_den[Nobs_BKsll+1];
	for(je=0;je<=Nobs_BKsll;je++) obs_num[je]=obs_den[je]=0.;

	/****** Outputting non-physical quantities - 20/04/2025 ******/
	double obs_min[Nobs_BKsll+1],obs_max[Nobs_BKsll+1];
	for(je=0;je<=Nobs_BKsll;je++) obs_max[je]=obs_min[je]=0.;
	/*************************************************************/

	obs[0]=0.; /* zero AFB */
	obs[1]=0.; /* integrated AFB */
	obs[2]=0.; /* integrated FL */
	obs[3]=0.; /* integrated FT */
	obs[4]=0.; /* integrated AT1 */
	obs[5]=0.; /* integrated AT2 */
	obs[6]=0.; /* integrated AT3 */
	obs[7]=0.; /* integrated AT4 */
	obs[8]=0.; /* integrated AT5 */
	obs[9]=0.; /* integrated HT1 */
	obs[10]=0.; /* integrated HT2 */
	obs[11]=0.; /* integrated HT3 */
	obs[12]=0.; /* integrated alpha */
	obs[13]=0.; /* integrated AIm */
	obs[14]=0.; /* integrated P2 */
	obs[15]=0.; /* integrated P3 */
	obs[16]=0.; /* integrated P6 */
	obs[17]=0.; /* integrated P4prime */
	obs[18]=0.; /* integrated P5prime */
	obs[19]=0.; /* integrated P6prime */
	obs[20]=0.; /* integrated P8 */
	obs[21]=0.; /* integrated P8prime */
	obs[22]=0.; /* integrated A7 */
	obs[23]=0.; /* integrated A8 */
	obs[24]=0.; /* integrated A9 */
	obs[25]=0.; /* integrated S3 */
	obs[26]=0.; /* integrated S4 */
	obs[27]=0.; /* integrated S5 */
	obs[28]=0.; /* integrated S7 */
	obs[29]=0.; /* integrated S8 */
	obs[30]=0.; /* integrated S9 */
	obs[31]=0.; /* integrated A3 */
	obs[32]=0.; /* integrated A4 */
	obs[33]=0.; /* integrated A5 */
	obs[34]=0.; /* integrated ATRe */
	obs[35]=0.; /* integrated ATReCP */
	obs[36]=0.; /* integrated P1CP */
	obs[37]=0.; /* integrated P2CP */
	obs[38]=0.; /* integrated P3CP */
	obs[39]=0.; /* integrated P4primeCP */
	obs[40]=0.; /* integrated P5primeCP */
	obs[41]=0.; /* integrated P6primeCP */
	obs[42]=0.; /* integrated P8primeCP */
	obs[43]=0.; /* integrated S6c */
	obs[44]=0.; /* integrated A6s */
	obs[45]=0.; /* integrated ACP */
	obs[46]=0.; /* integrated S1c */
	obs[47]=0.; /* integrated S2c */
	obs[48]=0.; /* integrated S2s */
	obs[49]=0.; /* integrated AFL */
	obs[50]=0.; /* integrated A1c */
	obs[51]=0.; /* integrated A2s */
	obs[52]=0.; /* integrated AFBCP */
	obs[53]=0.; /* integrated A6c */
	/****** Outputting non-physical quantities - 20/04/2025 ******/
	/* Form Factors */
	obs[Nobs_BKsll_phys+1]=0.; /* non-integrated (averaged) V */
	obs[Nobs_BKsll_phys+2]=0.; /* non-integrated (averaged) A0 */
	obs[Nobs_BKsll_phys+3]=0.; /* non-integrated (averaged) A1 */
	obs[Nobs_BKsll_phys+4]=0.; /* non-integrated (averaged) A12 */
	obs[Nobs_BKsll_phys+5]=0.; /* non-integrated (averaged) T1 */
	obs[Nobs_BKsll_phys+6]=0.; /* non-integrated (averaged) T2 */
	obs[Nobs_BKsll_phys+7]=0.; /* non-integrated (averaged) T23 */
	obs[Nobs_BKsll_phys+8]=0.; /* non-integrated (averaged) A2 */
	/* Result from Hadronic Fit */
	obs[Nobs_BKsll_phys+9]=0.; /* non-integrated (averaged) creal(hplus_fit) */
	obs[Nobs_BKsll_phys+10]=0.; /* non-integrated (averaged) cimag(hplus_fit) */
	obs[Nobs_BKsll_phys+11]=0.; /* non-integrated (averaged) creal(hminus_fit) */
	obs[Nobs_BKsll_phys+12]=0.; /* non-integrated (averaged) cimag(hminus_fit) */
	obs[Nobs_BKsll_phys+13]=0.; /* non-integrated (averaged) creal(hzero_fit) */
	obs[Nobs_BKsll_phys+14]=0.; /* non-integrated (averaged) cimag(hzero_fit) */
	/************************************************************/


	double dobs[Nobs_BKsll+1][3],dAFBtmp;
	double s0m,s0p,s0;
		
	dAFBtmp=0.;
	s0=s0m=0.;
	s0p=1.;
	
	double h=(smax-smin)/nmax;	
	s=smin;
	Gamma=dGamma_BKstarll_dq2_full(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKsll;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}
	/****** Outputting non-physical quantities  ******/
	/* Non-phyisal observables for any specific point in q2 is given as the
	 * average of the two points; define obs_min and obs_max to later average */
	for(je=1;je<=Nobs_BKsll;je++){
		obs_min[je] = dobs[je][0];
	}
	/*************************************************************/


	for(ie=1;ie<nmax;ie++)	
	{
		dAFBtmp=dobs[1][0];
		s0m=s0p;
		s+=h;
		s0p=s;
		
		Gamma+=4.*dGamma_BKstarll_dq2_full(gen,charge,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_BKsll;je++) 
		{
			obs_num[je]+=4.*dobs[je][1];
			obs_den[je]+=4.*dobs[je][2];
		}

		Gamma+=2.*dGamma_BKstarll_dq2_full(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_BKsll;je++) 
		{
			obs_num[je]+=2.*dobs[je][1];
			obs_den[je]+=2.*dobs[je][2];
		}

		if(dAFBtmp/dobs[1][0]<0.) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);
	}
	
	s=smax;
	Gamma+=4.*dGamma_BKstarll_dq2_full(gen,charge,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKsll;je++)
	{
		obs_num[je]+=4.*dobs[je][1];
		obs_den[je]+=4.*dobs[je][2];
	}
	Gamma+=dGamma_BKstarll_dq2_full(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKsll;je++)
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}

	Gamma*=h/6.;
	for(je=1;je<=Nobs_BKsll;je++)
	{
		obs_num[je]*=h/6.;
		obs_den[je]*=h/6.;
	}
	/****** Outputting non-physical quantities ******/
	for(je=1;je<=Nobs_BKsll;je++){
		obs_max[je] = dobs[je][0];
	}
	/*************************************************************/

	obs[0]=s0;
	for(je=1;je<=Nobs_BKsll;je++)
	if((je==17)||(je==18)||(je==19)||(je==21)||(je==39)||(je==40)||(je==41)||(je==42))
	{
		obs[je]=obs_num[je]/sqrt(fabs(obs_den[je]*obs_den[15]/4.));
	}
	else if(je==11)
	{
		obs[je]=obs_num[je]/sqrt(fabs(obs_den[5]*obs_den[5]-obs_den[je]*obs_den[je]));
	}
	else if(je==20)
	{
		obs[je]=obs_num[je]/sqrt(fabs(obs_den[je]*obs_den[19]));
	}
	/****** Outputting non-physical quantities ******/
	else if (je>Nobs_BKsll_phys)// Catching the non-physical observables (last lines appearing in `dGamma_BKstarll_dq2_full` function)
	{
		obs[je] = (obs_min[je] + obs_max[je]) / 2.;
	}
	/*************************************************************/
	else obs[je]=obs_num[je]/obs_den[je];

	for(je=1;je<=Nobs_BKsll;je++) if(fabs(obs[je])<1.e-15) obs[je]=0.;

	if(charge==0) return param->life_Bd/hbar*Gamma;
	else return param->life_B/hbar*Gamma;
}

/*----------------------------------------------------------------------*/

double BRBKstarll_lowq2_full(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKstarll_full(gen,charge,1.,6.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarll_highq2_full(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKstarll_full(gen,charge,14.18,16.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/

double BRobs_BKstarmumu_lowq2_full_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) and all the other observables */

	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarll_lowq2_full(2,0,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_BKstarmumu_highq2_full_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) and all the other observables */
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarll_highq2_full(2,0,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/
/*---------------------------- WRAPPER ---------------------------------*/
/*----------------------------------------------------------------------*/

double dGamma_BKstarll_dq2(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return dGamma_BKstarll_dq2_full(gen,charge,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return dGamma_BKstarll_dq2_soft(gen,charge,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double dGamma_BKstarll_dq2_calculator(int gen, int charge, double q2, double obs[][3], char name[])
{
	return dGamma_BKstarll_dq2_full_calculator(gen,charge,q2,obs,name);
}

/*----------------------------------------------------------------------*/

double BRBKstarll(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return BRBKstarll_full(gen,charge,smin,smax,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return BRBKstarll_soft(gen,charge,smin,smax,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarll_lowq2(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return BRBKstarll_lowq2_full(gen,charge,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return BRBKstarll_lowq2_soft(gen,charge,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarll_highq2(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return BRBKstarll_highq2_full(gen,charge,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return BRBKstarll_highq2_soft(gen,charge,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_BKstarmumu_lowq2_calculator(char name[], double obs[])
{
	return BRobs_BKstarmumu_lowq2_full_calculator(name,obs);
}
	
/*----------------------------------------------------------------------*/

double BRobs_BKstarmumu_highq2_calculator(char name[], double obs[])
{
	return BRobs_BKstarmumu_highq2_full_calculator(name,obs);
}

