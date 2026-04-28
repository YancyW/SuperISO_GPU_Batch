#include "src/include.h"

/*------------------------------------------------------*/
/* Calculation of chi2 for given NP Wilson coefficients */
/*------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50],totcovname[500];

	double dC9e,dC9mu,dC10e,dC10mu;
	int extcov=0;

  	if(argc<5) 
  	{ 
    		printf(" This program needs 4 parameters:\n"
           	"   dC9e, dC9mu, dC10e, dC10mu\n"
				"Optional parameter: total covariance matrix file name\n");
      		exit(1); 
  	} 
	else 
  	{
  		sscanf(argv[1],"%lf",&dC9e);
  		sscanf(argv[2],"%lf",&dC9mu);
  		sscanf(argv[3],"%lf",&dC10e);
  		sscanf(argv[4],"%lf",&dC10mu);
  	}
  	
  	if(argc>5)
  	{
		sscanf(argv[5],"%s",totcovname);
		extcov=1;
	}
  	  	
	printf("\n");
	
	printf("SuperIso " VERSION " - %d\n\n", YEAR);
	printf("Model-independent analysis\n\n");

	char obsfile[500],expfile[500],corrfile[500],nuisfile[500],nuiscorrfile[500];
	sprintf(obsfile,"chi2_input/myobs.in");
    sprintf(expfile,"chi2_input/exp_values.in");
    sprintf(corrfile,"chi2_input/exp_corr.in");
    sprintf(nuisfile,"chi2_input/nuisance.in");
    sprintf(nuiscorrfile,"chi2_input/nuisance_corr.in");
    
	struct parameters param;
	int ie=0;

	char names[NBOBSMAX][50];
	
	int nobs=read_obs_list(obsfile,names);
	printf("Reading observable list from %s\n",obsfile);
	
	Init_param(&param);
	
	slha_adjust(&param);
	param.SM=1;
	
	param.deltaC[9]=dC9mu;
	param.deltaC[10]=dC10mu;
	param.deltaC[10+9]=dC9e;
	param.deltaC[10+10]=dC10e;


	sprintf(param.exp_values,"%s",expfile);
	sprintf(param.exp_corr,"%s",corrfile);
	sprintf(param.nuisance_values,"%s",nuisfile);
	sprintf(param.nuisance_corr,"%s",nuiscorrfile);

	printf("\n");

	//param.LDsign_KLmumu=-1;

	//param.fullFF=0; /* uncomment to use soft FF approach for the inclusive b->s ll observables - full FF approach by default */
		
	double *predictions;
	double **covariance_th;
	get_predictions(names,&nobs,&predictions,&param); // Theory predictions

	double *central_exp; // Exp values
	double **covariance_exp;

	double **covariance_tot; // Total covariance matrix
		
	if(!extcov)
	{
		get_th_covariance(&covariance_th,names,&nobs,&param); // Theory covariance matrix calculation

		get_exp_covariance(&covariance_exp,&central_exp,names,&nobs,&param);

		for(ie=0;ie<nobs;ie++) printf("%-45s =\t% 6.5e +/- %.5e\t(exp: %4.5e +/- %.5e)\n",names[ie],predictions[ie],sqrt(covariance_th[ie][ie]),central_exp[ie],sqrt(covariance_exp[ie][ie]));

		get_covtot(&covariance_th,&covariance_exp,&covariance_tot,nobs);
	}
	else
	{
		get_exp_values(&central_exp,names,&nobs,&param);

		covariance_tot=malloc(nobs*sizeof(double *));
		for(ie=0;ie<nobs;ie++) covariance_tot[ie]=(double *) malloc(nobs*sizeof(double));
		read_covariance(totcovname,covariance_tot,names,nobs);

		for(ie=0;ie<nobs;ie++) printf("\t%-45s =\t%.5e +/- %.5e (th+exp)\t(exp: %.5e)\n",names[ie],predictions[ie],sqrt(covariance_tot[ie][ie]),central_exp[ie]);
	}

	double **inv_cov_tot; // Inverse of total covariance matrix
	if(get_invcovtot(&covariance_tot,&inv_cov_tot,nobs)<0) return -3.;
			
	double chi2=get_chi2(inv_cov_tot,predictions,central_exp,nobs);  // chi2 calculation
	printf("\nchi2=%.2f\n",chi2);
	printf("n_obs=%d\n",nobs);
	printf("with C9e=%.2f, C9mu=%.2f, C10e=%.2f, C10mu=%.2f\n\n",creal(dC9e),creal(dC9mu),creal(dC10e),creal(dC10mu));

	if(!extcov)
	{
		write_covariance("modelindep_covariance.out",covariance_tot,names,nobs);
		printf("Total covariance matrix written in modelindep_covariance.out\n");
	
		// Free variables
		for (ie = 0; ie < nobs; ie++) free(covariance_th[ie]);

		for (ie = 0; ie < nobs; ie++) free(covariance_exp[ie]);
	}
	
	free(covariance_th);
	free(covariance_exp);
	
	free(central_exp);

	for (ie = 0; ie < nobs; ie++) free(covariance_tot[ie]);
	free(covariance_tot);
	
	for (ie = 0; ie < nobs; ie++) free(inv_cov_tot[ie]);
	free(inv_cov_tot);
	
	return 1;
}

