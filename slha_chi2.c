#include "src/include.h"

/*--------------------------------------------------------*/
/* Calculation of the observables using a given SLHA file */
/*--------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	int test;

  	if(argc<2) 
  	{ 
    		printf(" This program needs 1 parameter:\n"
           	"   name    name of the SLHA file\n");
      		exit(1); 
  	} 
	else 
  	{
  		sscanf(argv[1],"%s",name);
  	}

	printf("\n");
	
	printf("SuperIso " VERSION " - %d\n\n", YEAR);
	printf("SLHA input file\n\n");

	test=test_slha(name);
	
	if(test>0)
	{
		if(test==2) printf("WARNING: only tested in the MFV scenario!\n\n");
	}
	else 
	{
		if(test==-1) printf("Invalid point\n\n");
		else if(test==-2) printf("Model not yet implemented\n\n");
		else if(test==-3) printf("Invalid SLHA file\n\n");
		else if(test==-4) printf("SLHA file absent\n\n");
		exit(1);
	}

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
	
	Les_Houches_Reader(name,&param);
	
	sprintf(param.exp_values,"%s",expfile);
	sprintf(param.exp_corr,"%s",corrfile);
	sprintf(param.nuisance_values,"%s",nuisfile);
	sprintf(param.nuisance_corr,"%s",nuiscorrfile);

	printf("\n");

	//param.LDsign_KLmumu=-1;

	//param.fullFF=0; /* uncomment to use soft FF approach for the inclusive b->s ll observables - full FF approach by default */
		
	double *predictions;
	get_predictions(names,&nobs,&predictions,&param); // Theory predictions

	double **covariance_th;
	get_th_covariance(&covariance_th,names,&nobs,&param); // Theory covariance matrix calculation

	double *central_exp,**covariance_exp; // Exp values
	get_exp_covariance(&covariance_exp,&central_exp,names,&nobs,&param);
	
	for(ie=0;ie<nobs;ie++) printf("\t%-25s =\t%.5e +/- %.5e\t(exp: %.5e +/- %.5e)\n",names[ie],predictions[ie],sqrt(covariance_th[ie][ie]),central_exp[ie],sqrt(covariance_exp[ie][ie]));
	
	double **covariance_tot; // Total covariance matrix
	get_covtot(&covariance_th,&covariance_exp,&covariance_tot,nobs);

	double **inv_cov_tot; // Inverse of total covariance matrix
	if(get_invcovtot(&covariance_tot,&inv_cov_tot,nobs)<0) return -3.;
			
	double chi2=get_chi2(inv_cov_tot,predictions,central_exp,nobs);  // chi2 calculation
	printf("\nchi2 = %.2f\n",chi2);
	printf("n_obs=%d\n\n",nobs);

	write_covariance("slha_covariance.out",covariance_tot,names,nobs);
	printf("Total covariance matrix written in slha_covariance.out\n");
	
	// Free variables
	for (ie = 0; ie < nobs; ie++) free(covariance_th[ie]);
	free(covariance_th);

	for (ie = 0; ie < nobs; ie++) free(covariance_exp[ie]);
	free(covariance_exp);
	
	free(central_exp);

	for (ie = 0; ie < nobs; ie++) free(covariance_tot[ie]);
	free(covariance_tot);
	
	for (ie = 0; ie < nobs; ie++) free(inv_cov_tot[ie]);
	free(inv_cov_tot);
	
	return 1;
}

