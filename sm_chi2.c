#include "src/include.h"
#include <stdbool.h>

/*--------------------------------------------------------*/
/* Calculation of the observables in the Standard Model   */
/*--------------------------------------------------------*/

int get_exponent(double number);
void format_with_matched_exponent(double main_value, double uncertainty, char *formatted_str, int precision);
void pretty_print_latex(char *str, double a, double b, double c, double d, int precision, int padding_size, bool latex_output);

#define PLAINDISPLAY 0.0050// smaller than this value is shown in scientific notation
#define PRECISION 2// precision of output
#define STRINGPADDING 35// padding of the names
#define LATEXOUPUT false// either `false` or `true` for latex ouput




int main(int argc,char** argv)
{
	printf("SuperIso " VERSION " - %d\n\n", YEAR);

	char obsfile[500],expfile[500],corrfile[500],nuisfile[500],nuiscorrfile[500];

	if(argc>=2)
    {
        sscanf(argv[1],"%s",obsfile);
    }
    else
    {
		sprintf(obsfile,"chi2_input/myobs.in");
    }

    if(argc>=3)
    {
        sscanf(argv[2],"%s",expfile);
    }
    else
    {
		sprintf(expfile,"chi2_input/exp_values.in");
    }

    if(argc>=4)
    {
        sscanf(argv[3],"%s",corrfile);
    }
    else
    {
		sprintf(corrfile,"chi2_input/exp_corr.in");
    }

    if(argc>=5)
    {
        sscanf(argv[4],"%s",nuisfile);
    }
    else
    {
		sprintf(nuisfile,"chi2_input/nuisance.in");
    }

    if(argc>=6)
    {
        sscanf(argv[5],"%s",nuiscorrfile);
    }
    else
    {
		sprintf(nuiscorrfile,"chi2_input/nuisance_corr.in");
    }


	struct parameters param;
	int ie=0;

	char names[NBOBSMAX][50];

	int nobs=read_obs_list(obsfile,names);
	printf("Reading observable list from %s\n",obsfile);

	Init_param(&param);

	slha_adjust(&param);
	param.SM=1;

	sprintf(param.exp_values,"%s",expfile);
	sprintf(param.exp_corr,"%s",corrfile);
	sprintf(param.nuisance_values,"%s",nuisfile);
	sprintf(param.nuisance_corr,"%s",nuiscorrfile);

	if(argc>=3) printf("Reading experimental values from %s and %s\n",expfile,corrfile);

	if(argc>=5) printf("Reading nuisance parameters from %s and %s\n",nuisfile,nuiscorrfile);

	printf("\n");

	//param.LDsign_KLmumu=-1;

	//param.fullFF=0; /* uncomment to use soft FF approach for the inclusive b->s ll observables - full FF approach by default */

	double *predictions;
	get_predictions(names,&nobs,&predictions,&param); // Theory predictions

	double **covariance_th;
	get_th_covariance(&covariance_th,names,&nobs,&param); // Theory covariance matrix calculation

	double *central_exp,**covariance_exp; // Exp values
	get_exp_covariance(&covariance_exp,&central_exp,names,&nobs,&param);

	for(ie=0;ie<nobs;ie++) printf("%-45s =\t% 6.5e +/- %.5e\t(exp: %4.5e +/- %.5e)\n",names[ie],predictions[ie],sqrt(covariance_th[ie][ie]),central_exp[ie],sqrt(covariance_exp[ie][ie]));

	double **covariance_tot; // Total covariance matrix
	get_covtot(&covariance_th,&covariance_exp,&covariance_tot,nobs);

	double **inv_cov_tot; // Inverse of total covariance matrix
	if(get_invcovtot(&covariance_tot,&inv_cov_tot,nobs)<0) return -3.;

	double chi2=get_chi2(inv_cov_tot,predictions,central_exp,nobs);  // chi2 calculation
	printf("\nchi2=%.3f\n",chi2);
	printf("n_obs=%d\n\n",nobs);

	write_covariance("sm_covariance.out",covariance_tot,names,nobs);
	printf("Total covariance matrix written in sm_covariance.out\n");

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


// Helper function to get the exponent of a number in scientific notation
int get_exponent(double number) {
    if (number == 0.0) return 0;
    return (int)floor(log10(fabs(number)));
}

// Helper function to format uncertainty with matched exponent and adjustable precision
void format_with_matched_exponent(double main_value, double uncertainty, char *formatted_str, int precision) {
    int main_exponent = get_exponent(main_value);
    double scaled_uncertainty = uncertainty / pow(10, main_exponent);

    // Create the format string with the specified precision
    char format[10];
    snprintf(format, sizeof(format), "%%.%df", precision);  // e.g., "%.3f" if precision is 3
    snprintf(formatted_str, 20, format, scaled_uncertainty);
}

// Main pretty-printing function with adjustable precision and LaTeX output option
void pretty_print_latex(char *str, double a, double b, double c, double d, int precision, int padding_size, bool latex_output) {
    if (latex_output) {
        // Determine if we should use plain or scientific notation based on `a`
        bool use_plain_notation = (fabs(a) >= PLAINDISPLAY);

        // Get the exponent for scaling (only used if scientific notation is needed)
        int exp_a = get_exponent(a);
        double scale_factor = use_plain_notation ? 1.0 : pow(10, exp_a);

        // Scale all values based on the notation decision
        double scaled_a = a / scale_factor;
        double scaled_b = b / scale_factor;
        double scaled_c = c / scale_factor;
        double scaled_d = d / scale_factor;

        // Determine if we need to add \phantom{-} for `a` and `c`
        const char *phantom_a = (scaled_a >= 0) ? "\\phantom{-}" : "";
        const char *phantom_c = (scaled_c >= 0) ? "\\phantom{-}" : "";

        // Print the first pair
        printf("%s = $ (%s%.*f \\pm %.*f)", str, phantom_a, precision, scaled_a, precision, scaled_b);

        // Add the scaling factor in scientific notation if needed
        if (!use_plain_notation) {
            printf(" \\times 10^{%d} $", exp_a);
        } else {
            printf(" $");
        }

        // Print the separator
        printf(" & ");

        // Print the second pair
        printf("$(%s%.*f \\pm %.*f)", phantom_c, precision, scaled_c, precision, scaled_d);

        // Add the scaling factor in scientific notation if needed
        if (!use_plain_notation) {
            printf(" \\times 10^{%d} $", exp_a);
        } else {
            printf(" $");
        }

        printf("\n");
    } else {
        // Regular format output (unchanged)
        char format[50];
        char uncertainty_str[20];

        if (padding_size > 0) {
            snprintf(format, sizeof(format), "%%-%ds = ", padding_size);
        } else {
            strcpy(format, "%s = ");
        }

        printf(format, str);

        char main_format[10];
        snprintf(main_format, sizeof(main_format), "%% .%df", precision);

        char scientific_format[10];
        snprintf(scientific_format, sizeof(scientific_format), "%% .%de", precision);

        if (fabs(a) >= PLAINDISPLAY) {
            printf(main_format, a);
        } else {
            printf(scientific_format, a);
        }

        printf(" ± ");

        if (fabs(a) >= PLAINDISPLAY) {
            printf("\b");
            printf(main_format, b);
        } else {
            format_with_matched_exponent(a, b, uncertainty_str, precision);
            printf("%se%d", uncertainty_str, get_exponent(a));
        }

        printf("\t(exp: %.4e ± %.4e)\n", c, d);
    }
}
