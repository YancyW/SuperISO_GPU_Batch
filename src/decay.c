/* This code was written by Mathieu Boudaud and Pierre Salati and then adapted to fit SuperIso Relic */

#include "include.h"
#include "include_dm.h"

#include "propagation.h"
#include "IDetection.h"



#define NAMS02 57

/********************************************************************************************/
/**************************Primary_pbar******************************************************/
/********************************************************************************************/

void calculation_BESSEL_PBAR_PRIMARY_Epbar_i_decay(long n_vert, long n_rad, struct Structure_Bessel* pt_Bessel,
                                             struct Structure_Pbar* pt_Pbar,
                                             struct Structure_Propagation* pt_Propagation,
                                             struct Structure_Primary_Source_Term* pt_Primary_Source_Term, struct propagation_parameters* pparam)
/* This routine calculates Bessel coefficients pt_Pbar->BESSEL_PBAR_PRI_Epbar_i
   for antiproton kinetic energy T_pbar et and Bessel coefficient i. */
{
  long i_pbar,i;
  double T_pbar,E_pbar,impulsion_pbar,v_pbar,K_pbar;
  double Si,Abar_i;

  long i_vert,i_rad;
  double x_vert,z_vert,dx_vert,weight_SIMPSON_vert;
  double x_rad, r_rad, dx_rad, weight_SIMPSON_rad;
  static double q_pbar_primary_i_z[2001];
  
  /* n_vert <= 1000 but never greater! */
  if (n_vert > 1000)
  {
    printf(" WARNING! n_vert > 1000 !");
    exit (0);
  }

  /* pt_Pbar->BESSEL_PBAR_PRI_Epbar_i is set to zero. */
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    for (i=0;i<=N_BESSEL;i++)
    {
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] = 0.0;
    }
  }
  
  /* pt_Pbar->BESSEL_PBAR_PRI_Epbar_i is filled. */
  for (i=1;i<=N_BESSEL;i++)
  {
/* The radial integral q_pbar_primary_i_z is first calculated.
   \beq
   q_pbar_primary_i_z = \frac{2}{J_{1}^{2}(\alpha_{i})} \times
   \int_{0}^{1} \, x \, dx \, J_{0} \left( \alpha_{i} x \right) \,
   \left\{ \frac{\rho_{\chi}}{\rho_0} \right\}^{2} \;\; .
   \eeq

   This integral depends on i and z. Results corresponding to the ith Bessel order are saved 
   in array q_pbar_primary_i_z which depends only on z.

   Loop over z followed by a loop over r */
  
    dx_vert = 1. / (double) (2*n_vert);
    x_vert  = 0.;
    for (i_vert=0;i_vert<=(2*n_vert);i_vert++)
    {
      
      z_vert = pt_Propagation->E_DIFFUS * x_vert; /* z_vert is the vertical coordinate expressed in [kpc]. */
      q_pbar_primary_i_z[i_vert] = 0.0;
      dx_rad = 1. / (double) (2*n_rad);
      x_rad  = 0.0;
       
      for (i_rad=0;i_rad<=(2*n_rad);i_rad++)
      {
       
        r_rad = R_GAL * x_rad; /* r_rad is the radial coordinate expressed in [kpc]. */
        if (i_rad==0 || i_rad==(2*n_rad)) {weight_SIMPSON_rad = 1./3.;}
        else {weight_SIMPSON_rad = (1. + (double)(i_rad % 2)) * 2. / 3.;}
/*		  printf("llol2 %e \n",einasto_CU10(1,1));*/
	    q_pbar_primary_i_z[i_vert] += (x_rad * dx_rad * weight_SIMPSON_rad) * besselj0(pt_Bessel->alpha_i[i]*x_rad) * pparam->Halo.profile(r_rad,z_vert); //[NO UNIT].
/*printf("llol\n");*/
		x_rad += dx_rad;
		
      }
      q_pbar_primary_i_z[i_vert] *= 2. / pow(besselj1(pt_Bessel->alpha_i[i]),2.0); /* [NO UNIT]. */

      x_vert += dx_vert;
    }
    
 
/* Loop over antiproton energy containing a loop over the vertical coordinate,
   allowing us to calculate calculer pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i]. */
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      T_pbar = T_PBAR_MIN *
      pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
      E_pbar = T_pbar + MASSE_PROTON;
      impulsion_pbar = sqrt(pow(E_pbar,2) - pow(MASSE_PROTON,2));
      v_pbar         = CELERITY_LIGHT * impulsion_pbar / E_pbar;
      K_pbar         = K_space_diffusion(E_pbar,MASSE_PROTON,1.0,pt_Propagation);
	  

      Si = sqrt(pow(2.0*pt_Bessel->alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K_pbar,2)); /* Si is expressed in [kpc^{-1}]. */

      Abar_i  = pt_Propagation->VENT_GALACTIQUE; /* Abar_i is expressed in [cm s^{-1}]. */
	  Abar_i += 2.0*E_DISC*CM_PAR_KPC *
      ((sigma_inelastic_pbarH_TAN_and_NG(E_pbar) - sigma_inelastic_NOANN_pbarH_TAN_and_NG(E_pbar)) * v_pbar * (DENSITE_H_DISC + pow(4.,(2./3.))*1.0*DENSITE_HE_DISC));
      Abar_i += K_pbar * Si / CM_PAR_KPC / tanh(Si*pt_Propagation->E_DIFFUS/2.);

	  pt_Pbar->TABLE_Abar_i[i_pbar][i] = Abar_i;


	  /* Integration over the vertical variable x_vert and over the radial variable x_rad. */
      dx_vert = 1. / (double) (2*n_vert);
      x_vert  = 0.;
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] = 0.0;
      for (i_vert=0;i_vert<=(2*n_vert);i_vert++)
      {
        z_vert = pt_Propagation->E_DIFFUS * x_vert;

        if (i_vert==0 || i_vert==(2*n_vert)) {weight_SIMPSON_vert = 1./3.;}
        else {weight_SIMPSON_vert = (1. + (double)(i_vert % 2)) * 2. / 3.;}
      
        pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] += (dx_vert * weight_SIMPSON_vert) *
        q_pbar_primary_i_z[i_vert] *
        pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar] *
        exp(- pt_Propagation->VENT_GALACTIQUE*z_vert*CM_PAR_KPC / (2.*K_pbar)) *
        sinh((Si/2.)*(pt_Propagation->E_DIFFUS-z_vert)) / sinh((Si/2.)*pt_Propagation->E_DIFFUS); //[antiprotons cm^{-3} s^{-1} GeV^{-1}].
		  
		x_vert += dx_vert;
      }
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] *= 2. * pt_Propagation->E_DIFFUS*CM_PAR_KPC / Abar_i; /* expressed in [antiprotons cm^{-3} GeV^{-1}]. */
/*      printf("grou %e\n",  pt_Propagation->VENT_GALACTIQUE);*/
    }
  }
  return;
}

/********************************************************************************************/





/********************************************************************************************/
/********************************************************************************************/

double dNpbar_on_dEpbar_primary_calculation_Br_decay(char* filename, struct Structure_Primary_Source_Term* pt_Primary_Source_Term)
/* Calculates the primary pbar spectrum for a DM particle whose mass and annihilation channels ared defined in spec and save it in pt_Primary_Source_Term */
{

struct array spectrum;
read_spectrumfile(filename, &spectrum, 7, 7);


	int    i_pbar, i_scan_pbar;
	double xi_x_pbar, T_pbar, x_pbar;
	double dNpbar_on_dEpbar_xi_m, dNpbar_on_dEpbar_xi_plus_un_m;
	
	int index;

	int found;

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar = T_PBAR_MIN * pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		
		found=0;
		index=0;
		while(!found&&index<spectrum.n1){
		if(spectrum.a[index][0]>T_pbar+mp) found=1;
		else index++;
		}
		
		if(index==0)index++;
		if(index==spectrum.n1) index--;
/*		printf("warning antiproton spectrum data missing at kinectic enegry %e\n", T_pbar);*/
/*		else{*/
		pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]  = spectrum.a[index-1][1]+
	  (spectrum.a[index][1]-spectrum.a[index-1][1])/(spectrum.a[index][0]-spectrum.a[index-1][0])*(T_pbar+mp-spectrum.a[index-1][0]);
/*	  }*/
/*	 printf(" test %e %e\n",T_pbar, pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]);	*/
	}

	return (1.0);
}
/********************************************************************************************/
/********************************************************************************************/

/* Warning! This function must be preceeded by 'DNPBAR_ON_DTPBAR_gaelle_read_file' and 'dNpbar_on_dEpbar_primary_calculation' */

void primary_source_calculation_decay (double mass_chi,double fraction, struct Structure_Primary_Source_Term* pt_Primary_Source_Term)
/* modulates the primary pbar spectrum at production with the total annihilation cross section and the DM number density */
{
	int i_pbar;

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
		pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar] = 0.0;
	}
	
	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{			
		if(pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]>1.0e-100)pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar] = fraction* (RHO_CHI_0/mass_chi)*pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]; // [#pbar cm{-3} s{-1} GeV{-1}]
/*		printf("hola %e\n", pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar]);*/
	}
		
	return;
}
	
/********************************************************************************************/



/********************************************************************************************/

void primary_spectra_BCGS_2014_Br_decay(char* filename,double mass_chi, double fraction,struct Structure_Pbar* pt_Pbar, struct Structure_Cross_Section* pt_Cross_Section, struct Structure_Propagation* pt_Propagation, struct Structure_Primary_Source_Term* pt_Primary_Source_Term, struct Structure_Bessel* pt_Bessel,struct propagation_parameters*pparam, struct array* PBAR_SPECTRUM)
/* Calculates the primary pbar spectrum at Earth position for a DM particle defined in structure spec and save the result in array PBAR_SPECTRUM */
{
/*	char file_name[1024];*/

	
	double perc=0;
	long i_data,i_iteration,i_pbar,i;
	
	double T_pbar_IS ,E_pbar_IS ,flux_antiproton_IS ,flux_proton_IS;
	
	double FLUX_PBAR;

	PBAR_SPECTRUM->n1=DIM_TAB_PBAR+1;
	PBAR_SPECTRUM->n2=2;
	
	DM_preliminary(pt_Primary_Source_Term);

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar] = 0.0;
	}
	dNpbar_on_dEpbar_primary_calculation_Br_decay(filename, pt_Primary_Source_Term);
	
	
	primary_source_calculation_decay(mass_chi, fraction, pt_Primary_Source_Term);
			
	/* We define the parameters determined by FIORENZA, DAVID et RICHARD -- hereafter called FDR. */
	
	/* pt_Pbar->BESSEL_PBAR_SEC_Epbar_i and pt_Pbar->BESSEL_PBAR_TER_Epbar_i are set to zero. */
	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		for (i=0;i<=N_BESSEL;i++)
		{
			pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] = 0.0;
			pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] = 0.0;
			pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i] = 0.0;

			pt_Pbar->BESSEL_PBAR_SEC_SNR_Epbar_i[i_pbar][i] = 0.0;
		}
	}
		
	/* CALCULATION OF THE PRIMARY CONTRIBUTION FROM NEUTRALINO ANNIHILATION */
	
	calculation_BESSEL_PBAR_PRIMARY_Epbar_i_decay(100,500, pt_Bessel, pt_Pbar, pt_Propagation, pt_Primary_Source_Term, pparam);
	
	/* CALCULATION OF SECONDARY CONTRIBUTION FROM THE INTERSTELLAR GAS SPALLATION BY COSMIC-RAY PROTONS AND HELIONS */
	calculation_BESSEL_PBAR_SUM_123_Epbar_i(pt_Pbar);

	/* CALCULATION OF FINAL ANTIPROTON SPECTRUM */
	for (i_iteration=1;i_iteration<=5;i_iteration++)
	{
		calculation_BESSEL_PBAR_TERTIARY_Epbar_i(pt_Bessel, pt_Pbar);
		calculation_BESSEL_PBAR_SUM_123_Epbar_i(pt_Pbar);
	}

	for (i_iteration=1;i_iteration<=5;i_iteration++)
	{
		calculation_BESSEL_PBAR_TERTIARY_Epbar_i(pt_Bessel, pt_Pbar);
		calculation_BESSEL_PBAR_TOT_direct_inversion_A(pt_Pbar, pt_Propagation);
	}

	/* We compute now the antiproton spectrum */
	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar_IS = T_PBAR_MIN *
		pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		PBAR_SPECTRUM->a[i_pbar][0]=log10(T_pbar_IS) ;
		E_pbar_IS = T_pbar_IS + MASSE_PROTON;
		for (i=1;i<=N_BESSEL;i++)
		{
			pt_Pbar->BESSEL_PBARi[i] = pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i];
		}
		flux_antiproton_IS = GENERIC_FLUX_04(pparam->Halo.r_earth,0.,E_pbar_IS,MASSE_PROTON,1.,pt_Bessel,pt_Pbar->BESSEL_PBARi, pt_Propagation);
		PBAR_SPECTRUM->a[i_pbar][1] = flux_antiproton_IS*1.0e+4/3.e-26;		
		
	}
				
	return;
}


void spectrum_pbar_after_propagation_decay(char*filename, double mass,double fraction, struct array* result, struct propagation_parameters* pparam)
/* Calculates pbar spectrum after propagation "result", from propagation parameters "pparam" and DM annihilation spectrum "spec" */
{
	
		struct Structure_Nuclei					Proton;
		struct Structure_Nuclei					Helium;
		struct Structure_Primary_Source_Term	Primary_Source_Term;
		struct Structure_Bessel					Bessel;

		int j;

		search_for_zero(&Bessel);
		production(&Bessel, &Proton);
		production(&Bessel, &Helium);

		struct Structure_Pbar					Pbar;
		init_pbar(&Pbar);
		struct Structure_Cross_Section			Cross_Section;
		init_cross_section(&Cross_Section);
			
		primary_spectra_BCGS_2014_Br_decay(filename, mass, fraction, &Pbar, &Cross_Section, &(pparam->Propagation), &Primary_Source_Term, &Bessel, pparam,  result);

		free_pbar(&Pbar);
		free_cross_section(&Cross_Section);
	
	
	int i;
	
	for(i=0;i<DIM_TAB_PBAR;i++)
	{
		if(result->a[i][1]<1.e-100)result->a[i][1]=0;
	}
	return;
}


double deltachi2_decay(char*filename, double mass,double fraction, struct propagation_parameters* pparam)
/* computes AMS-02 delta-chi2 for a DM particle with annihilation spectrum spec and halo and propagation parameters pparam */
{
	
	
	int i;
	double buff;
	
	double chi;
	double chiDM;
	double phi_f, phi_min, chi_min, phi_DMmin, chiDM_min,A, Amin;
	
	init_background(pparam);
	
	if(pparam->chi2_noDM<-1.e29) return -1.e30;

	struct array spect;
	 init_existing_model(pparam,"Einasto_CU10", "MAX");
/*	 printf("grou %e\n",  pparam->Propagation.VENT_GALACTIQUE);*/
	spectrum_pbar_after_propagation_decay(filename, mass, fraction, &spect, pparam);
/*	for(i=0;i<spect.n1;i++) printf("spect %e %e\n", spect.a[i][0], spect.a[i][1]);*/
	chiDM_min=0;
	double A_DM;
	double logE[NAMS02+1];
	for(i=0;i<NAMS02;i++) logE[i]=log10(table_AMS02(i,0));
	logE[NAMS02]=2.6;
   
	double xlimmin[2]={0,0.1};
	double xlimmax[2]={1,1.1};
	double xmin[2];

	struct fermi fe;
	IDpowell(2, dchi_tot, logE, spect, pparam, &fe, xlimmin, xlimmax, &chiDM_min, xmin, 1.0e-3, "min");
	
	return chiDM_min-pparam->chi2_noDM ; 
}



double deltalikelihood_AMS_decay(char*filename, double mass, double fraction,char* option)
/* Computes AMS-02 deltachi2 for option "conservative", "standard" or "stringent" */
{
	
	int ioption=readoption(option);
	struct propagation_parameters pparam;
	if(ioption==1){
	init_halo(&pparam, "Einasto_CU10");
	init_propagation(&pparam, "MAX");
	}
	else if(ioption==0){
	init_halo(&pparam, "Einasto_CU10");
	init_propagation(&pparam, "MED");
	}
	else if(ioption==-1){
	init_halo(&pparam, "Burkert_NS13");
	init_propagation(&pparam, "MED");
	}


	return -deltachi2_decay(filename, mass, fraction, &pparam);
}



int read_spectrumfile(char* filename, struct array* res, int ncolumn, int ncolumns){
//read spectrum from file filename containing ncolumns columns. ncolumn in the column numer we want to extract (ncolumn=2 for gammas)
	FILE* ff=fopen(filename, "r");
	if(ff==NULL)
	{
		printf("could not open file read_file %s\n", filename);
		exit(EXIT_FAILURE);
	}
	else{
	char c;
	int i, nline;
	double buff[ncolumns];
	
	nline=0;
	
	do{
	c=getc(ff);
	if(c=='#')while(c!='\n') c=getc(ff);
	else{
	ungetc(c,ff);
	for(i=0; i<ncolumns;i++) fscanf(ff,"%lf", &(buff[i]));
	res->a[nline][0]=buff[0];
	res->a[nline][1]=buff[ncolumn-1];
	nline++;
	}
	}while(!feof(ff));
	 nline--;
res->n1=nline;
	res->n2=2;
	
/*	for(i=0;i<res->n1;i++) printf("%e %e\n", res->a[i][0], res->a[i][1]);*/
	}


}

/***************************************************/
int test_fermi(char* filename, double mass, double fraction){
//tests a gamma ray spectrum from filename for blackholes of mass mass representing a fraction of DM
	int i;
	struct array res;
	read_spectrumfile(filename,&res,2,7);
	
	struct fermi fe;
	init_fermi(&fe);
	
	for (i=0; i<fe.Nenergybins; i++) fe.eflux[i]=Spect_integration( res, fe.Emin[i],  fe.Emax[i])/(4.*pi)/(mass)*fraction;
	
	double lJobs, deltaj;
	double tmp=0;
	for (i=0;i<fe.NdSphs;i++)
	{	
		lJobs=fe.dSphs[i].D;
		deltaj=fe.dSphs[i].deltaD;
		

		if(lJobs>-1) tmp+=maxlikelihood(&fe, i,lJobs, deltaj)-fe.dSphs[i].likelihood_noDM;
	}
	
	
	
/*	FILE* ff=fopen("test2.out", "w");*/
/*	for (i=0; i<fe.Nenergybins; i++) fprintf(ff, "%e\t %e\n", fe.Emin[i], fe.eflux[i]);*/
/*	fclose(ff);*/
/*	*/
	
	int dof=1;
	
	int sig=find_isigma("2sigma");	
		
	printf("\n deltaloglikelihood %e\n", 2*tmp);
	printf("Point excluded at 2 sigma if deltaloglikelihood < %f with %d degrees of freedom\n", -pvalue_fromchi2(sig,dof), dof);
	return (2*tmp<-pvalue_fromchi2(sig,dof));
		

}


int indirect_ams02_calculator_decay(char name[],double mass, double fraction)
/* "container" function scanning the SLHA file "name" and checking if the point is excluded by AMS-02 */
{

		double dchi2_standard=deltalikelihood_AMS_decay(name, mass,fraction, "standard");
		
		double dchi2lim;
		
		
	int dof=1;
	int sig=find_isigma("2sigma");	
		
		printf("\n deltaloglikelihood %e\n", dchi2_standard);
	printf("Point excluded at 2 sigma if deltaloglikelihood < %f with %d degrees of freedom\n", -pvalue_fromchi2(sig,dof), dof);
		return (dchi2_standard<-pvalue_fromchi2(sig,dof));
	
}

